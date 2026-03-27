/*------------------------------------------------------------------------------
* postpos.c : post-processing positioning
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/05/08  1.0  new
*           2008/06/16  1.1  support binary inputs
*           2009/01/02  1.2  support new rtk positioning api
*           2009/09/03  1.3  fix bug on combined mode of moving-baseline
*           2009/12/04  1.4  fix bug on obs data buffer overflow
*           2010/07/26  1.5  support ppp-kinematic and ppp-static
*                            support multiple sessions
*                            support sbas positioning
*                            changed api:
*                                postpos()
*                            deleted api:
*                                postposopt()
*           2010/08/16  1.6  fix bug sbas message synchronization (2.4.0_p4)
*           2010/12/09  1.7  support qzss lex and ssr corrections
*           2011/02/07  1.8  fix bug on sbas navigation data conflict
*           2011/03/22  1.9  add function reading g_tec file
*           2011/08/20  1.10 fix bug on freez if solstatic=single and combined
*           2011/09/15  1.11 add function reading stec file
*           2012/02/01  1.12 support keyword expansion of rtcm ssr corrections
*           2013/03/11  1.13 add function reading otl and erp data
*           2014/06/29  1.14 fix problem on overflow of # of satellites
*           2015/03/23  1.15 fix bug on ant type replacement by rinex header
*                            fix bug on combined filter for moving-base mode
*           2015/04/29  1.16 fix bug on reading rtcm ssr corrections
*                            add function to read satellite fcb
*                            add function to read stec and troposphere file
*                            add keyword replacement in dcb, erp and ionos file
*           2015/11/13  1.17 add support of L5 antenna phase center parameters
*                            add *.stec and *.trp file for ppp correction
*           2015/11/26  1.18 support opt->freqopt(disable L2)
*           2016/01/12  1.19 add carrier-phase bias correction by ssr
*           2016/07/31  1.20 fix error message problem in rnx2rtkp
*           2016/08/29  1.21 suppress warnings
*           2016/10/10  1.22 fix bug on identification of file fopt->blq
*           2017/06/13  1.23 add smoother of velocity solution
*           2020/11/30  1.24 use API sat2freq() to get carrier frequency
*                            fix bug on select best solution in static mode
*                            delete function to use L2 instead of L5 PCV
*                            writing solution file in binary mode
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define MIN(x,y)    ((x)<(y)?(x):(y))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

#define MAXPRCDAYS  100          /* max days of continuous processing */
#define MAXINFILE   1000         /* max number of input files */
#define MAXINVALIDTM 100         /* max number of invalid time marks */

/* constants/global variables ------------------------------------------------*/
static spcvs_t pcvss={0};        /* satellite antenna parameters */
static rpcvs_t pcvsr={0};        /* receiver antenna parameters */
static imu_t imus={0};          /* imu data */
static pos_t poss={0};          /* pos data */
static obs_t obss={0};          /* observation data */
static nav_t navs={0};          /* navigation data */
static sbs_t sbss={0};          /* sbas messages */
static sta_t stas[MAXRCV];      /* station information */
static vrs_t vrs={0};           /* vrs data */
static int nepoch=0;            /* number of observation epochs */
static int nitm  =0;            /* number of invalid time marks */
static int iobsu =0;            /* current rover observation data index */
static int iobsr =0;            /* current reference observation data index */
static int iimu  =0;            /* current imu data index */
static int ipos  =0;            /* current pos data index */
static int isbs  =0;            /* current sbas message index */
static int aborts=0;            /* abort status */
static sol_t *solf;             /* forward solutions */
static sol_t *solb;             /* backward solutions */
static double *rbf;             /* forward base positions */
static double *rbb;             /* backward base positions */
static int isolf=0;             /* current forward solutions index */
static int isolb=0;             /* current backward solutions index */
static int reverse_flag=0;      /* reverse flag for GNSS/INS forward and backward smoothing */
static char proc_rov [64]="";   /* rover for current processing */
static char proc_base[64]="";   /* base station for current processing */
static char rtcm_file[1024]=""; /* rtcm data file */
static char rtcm_path[1024]=""; /* rtcm data path */
static rtcm_t rtcm;             /* rtcm control struct */
static FILE *fp_rtcm=NULL;      /* rtcm data file pointer */
static int checkbrk_counter=0;
static char last_buff[256]="";

/* show message and check break ----------------------------------------------*/
static int checkbrk(const char *format, ...)
{
    va_list arg;
    char buff[256],*p=buff;

    /* only check every 5th call */
    if (++checkbrk_counter%5!=0) return 0;

    if (!*format) return showmsg("");
    va_start(arg,format);
    p+=vsprintf(p,format,arg);
    va_end(arg);

    if (*proc_rov&&*proc_base) sprintf(p," (%s-%s)",proc_rov,proc_base);
    else if (*proc_rov ) sprintf(p," (%s)",proc_rov );
    else if (*proc_base) sprintf(p," (%s)",proc_base);

    /* only show when the message content changes */
    if (strcmp(buff,last_buff)!=0) {
        strcpy(last_buff, buff);
        return showmsg(buff);
    }

    return 0;
}
/* Solution option to field separator ----------------------------------------*/
/* Repeated from solution.c */
static const char *opt2sep(const solopt_t *opt)
{
    if (!*opt->sep) return " ";
    else if (!strcmp(opt->sep,"\\t")) return "\t";
    return opt->sep;
}
/* output reference position -------------------------------------------------*/
static void outrpos(FILE *fp, const double *r, const solopt_t *opt)
{
    double pos[3],dms1[3],dms2[3];

    trace(3,"outrpos :\n");

    const char *sep = opt2sep(opt);
    if (opt->posf==SOLF_LLH||opt->posf==SOLF_ENU) {
        ecef2pos(r,pos);
        if (opt->degf) {
            deg2dms(pos[0]*R2D,dms1,5);
            deg2dms(pos[1]*R2D,dms2,5);
            fprintf(fp,"%3.0f%s%02.0f%s%08.5f%s%4.0f%s%02.0f%s%08.5f%s%10.4f",
                    dms1[0],sep,dms1[1],sep,dms1[2],sep,dms2[0],sep,dms2[1],
                    sep,dms2[2],sep,pos[2]);
        }
        else {
            fprintf(fp,"%13.9f%s%14.9f%s%10.4f",pos[0]*R2D,sep,pos[1]*R2D,
                    sep,pos[2]);
        }
    }
    else if (opt->posf==SOLF_XYZ) {
        fprintf(fp,"%14.4f%s%14.4f%s%14.4f",r[0],sep,r[1],sep,r[2]);
    }
}
/* output header -------------------------------------------------------------*/
static void outheader(FILE *fp, const prcopt_t *popt, const solopt_t *sopt, const filopt_t *fopt)
{
    const char *s1[]={"GPST","UTC","JST"};
    gtime_t ts,te;
    double t1,t2;
    int w1,w2;
    char s2[32],s3[32];

    if (sopt->posf==SOLF_NMEA||sopt->posf==SOLF_STAT) {
        return;
    }
    if (sopt->outhead) {
        if (!*sopt->prog) {
            fprintf(fp,"%s program   : GINSLIB ver.%s %s\n",COMMENTH,VER_RTKLIB,PATCH_LEVEL);
        }
        else {
            fprintf(fp,"%s program   : %s\n",COMMENTH,sopt->prog);
        }

        if (PMODE_LC_POS!=popt->mode) {
            if (*fopt->obs_u) fprintf(fp,"%s inp file  : %s\n",COMMENTH,fopt->obs_u);
            if (PMODE_DGPS<=popt->mode&&PMODE_FIXED>=popt->mode&&*fopt->obs_b) fprintf(fp,"%s inp file  : %s\n",COMMENTH,fopt->obs_b);
            if (*fopt->nav) fprintf(fp,"%s inp file  : %s\n",COMMENTH,fopt->nav);            
        }
        else  {
            if (*fopt->pos) fprintf(fp,"%s inp file  : %s\n",COMMENTH,fopt->pos);
        }

        if (GINS_OFF!=popt->GI_mode) fprintf(fp,"%s inp file  : %s\n",COMMENTH,fopt->imu);
        if (EPHOPT_PREC==popt->sateph&&PMODE_LC_POS!=popt->mode) {
            if (*fopt->sp3) fprintf(fp,"%s inp file  : %s\n",COMMENTH,fopt->sp3); 
            if (*fopt->clk) fprintf(fp,"%s inp file  : %s\n",COMMENTH,fopt->clk);
        }
        if (*fopt->mgex_dcb&&PMODE_LC_POS!=popt->mode) fprintf(fp,"%s inp file  : %s\n",COMMENTH,fopt->mgex_dcb);

        ts=popt->ts;
        te=popt->te;
        t1=time2gpst(ts,&w1);
        t2=time2gpst(te,&w2);
        if (sopt->times>=1) {
            ts=gpst2utc(ts);
            te=gpst2utc(te);
        }
        if (sopt->times==2) {
            ts=timeadd(ts,9*3600.0);
            te=timeadd(te,9*3600.0);
        }
        time2str(ts,s2,1);
        time2str(te,s3,1);
        fprintf(fp,"%s obs start : %s  (%s %04d %8.1fs)\n",COMMENTH,s2,s1[sopt->times],w1,t1);
        fprintf(fp,"%s obs end   : %s  (%s %04d %8.1fs)\n",COMMENTH,s3,s1[sopt->times],w2,t2);
    }
    if (sopt->outopt) {
        outprcopt(fp,popt);
    }
    if (PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED&&popt->mode!=PMODE_MOVEB) {
        fprintf(fp,"%s ref pos   :",COMMENTH);
        outrpos(fp,popt->rb,sopt);
        fprintf(fp,"\n");
    }
    if (sopt->outhead||sopt->outopt) fprintf(fp,"%s\n",COMMENTH);

    outsolhead(fp,popt,sopt);
}
/* search next observation data index ----------------------------------------*/
static int nextobsf(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;

    for (;*i<obs->n;(*i)++) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i+n<obs->n;n++) {
        tt=timediff(obs->data[*i+n].time,obs->data[*i].time);
        if (obs->data[*i+n].rcv!=rcv||tt>DTTOL) break;
    }
    return n;
}
/* search previous observation data index ----------------------------------------*/
static int nextobsb(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;

    for (;*i>=0;(*i)--) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i-n>=0;n++) {
        tt=timediff(obs->data[*i-n].time,obs->data[*i].time);
        if (obs->data[*i-n].rcv!=rcv||tt<-DTTOL) break;
    }
    return n;
}
/* update rtcm ssr correction ------------------------------------------------*/
static void update_rtcm_ssr(gtime_t time)
{
    char path[1024];
    int i;

    /* open or swap rtcm file */
    reppath(rtcm_file,path,time,"","");

    if (strcmp(path,rtcm_path)) {
        strcpy(rtcm_path,path);

        if (fp_rtcm) fclose(fp_rtcm);
        fp_rtcm=fopen(path,"rb");
        if (fp_rtcm) {
            rtcm.time=time;
            input_rtcm3f(&rtcm,fp_rtcm);
            trace(2,"rtcm file open: %s\n",path);
        }
    }
    if (!fp_rtcm) return;

    /* read rtcm file until current time */
    while (timediff(rtcm.time,time)<1E-3) {
        if (input_rtcm3f(&rtcm,fp_rtcm)<-1) break;

        /* update ssr corrections */
        for (i=0;i<MAXSAT;i++) {
            if (!rtcm.ssr[i].update||
                rtcm.ssr[i].iod[0]!=rtcm.ssr[i].iod[1]||
                timediff(time,rtcm.ssr[i].t0[0])<-1E-3) continue;
            navs.ssr[i]=rtcm.ssr[i];
            rtcm.ssr[i].update=0;
        }
    }
}
/* input obs data, navigation messages and sbas correction -------------------*/
static int inputobs(rtk_t *rtk, obsd_t *obs, imud_t *imu, int stat, const prcopt_t *popt)
{
    gtime_t time={0},gnss_time={0};
    ins_t *ins=&rtk->ins;
    int i,nu=0,nr=0,npos=(PMODE_LC_POS==popt->mode)?1:0,n=0,nn=ins->nn;
    double dt,dt_next,GI_dt,sec,ndt;

    trace(3,"infunc  : dir=%d iobsu=%d iobsr=%d isbs=%d\n",popt->reverse,iobsu,iobsr,isbs);

    stat=(GINS_LC==popt->GI_mode||GINS_STC==popt->GI_mode||PMODE_LC_POS==popt->mode)?rtk->lcgins.sol.stat:rtk->sol.stat;

    if ((0<=iobsu&&iobsu<obss.n)||PMODE_LC_POS==popt->mode) {
        time=(GINS_OFF==popt->GI_mode)?obss.data[iobsu].time:imus.data[iimu].time;
        if (checkbrk("processing : %s Q=%d",time_str(time,0),stat)) {
            aborts=1; showmsg("aborted"); return -1;
        }            
    }
    /* input forward data */
    if (SOLTYPE_FORWARD==popt->reverse) 
    { 
        if (PMODE_LC_POS!=popt->mode) 
        {
            if ((nu=nextobsf(&obss,&iobsu,1))<=0) return -1;
            if (popt->intpref) {
                /* for interpolation, find first base timestamp after rover timestamp */
                for (;(nr=nextobsf(&obss,&iobsr,2))>0;iobsr+=nr)
                    if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)>-DTTOL) break;
            }
            else {
                /* if not interpolating, find closest timestamp */
                dt=timediff(obss.data[iobsr].time,obss.data[iobsu].time);
                for (i=iobsr;(nr=nextobsf(&obss,&i,2))>0;iobsr=i,i+=nr) {
                    dt_next=timediff(obss.data[i].time,obss.data[iobsu].time);
                    if (fabs(dt_next)>fabs(dt)) break;
                    dt=dt_next;
                }
            }
            nr=nextobsf(&obss,&iobsr,2);
            if (nr<=0) {
                nr=nextobsf(&obss,&iobsr,2);
            }
            /* NOTE: store the observations of rover and base in the obs structure in order */
            for (i=0;i<nu&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsu+i];
            for (i=0;i<nr&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsr+i];

            gnss_time=obss.data[iobsu].time;
        }
        else {
            pos2sol(poss,&rtk->sol,ipos);
            gnss_time=rtk->sol.time;
        }

        /* NOTE: GNSS/INS time synchronization */
        if (GINS_LC==popt->GI_mode||GINS_TC==popt->GI_mode||GINS_STC==popt->GI_mode)
        {
            if (iimu>=imus.n) return -1;

            /* calculate the difference between the current time and the nominal measurement update time */
            sec=imus.data[iimu+(nn-1)].time.sec;
            rtk->nominal_upte=NO;
            ndt=fabs(sec-round((sec+ins->interval/2.0)/rtk->interval)*rtk->interval);
            /* INS navigation time */
            ins->time=imus.data[iimu+(nn-1)].time;

            /* the second condition is used to process the IMU time on both sides being 0.5 IMU sampling interval away from GNSS time */
            if ((fabs(ndt)-ins->dttol)<=(ins->nn*ins->interval)/2.0
                &&(fabs(timediff(ins->time,rtk->upte_time))>=(rtk->interval-ins->interval))) {
                rtk->nominal_upte=YES;
                /* record the synchronization time */
                rtk->upte_time=ins->time;                
            }  

            /* calculate the difference between the current IMU and GNSS observation time */
            GI_dt=timediff(ins->time,gnss_time);
            rtk->upte=SYNC_NO;   
            imucpy(popt,imu,imus,iimu,ins->nn); 

            /* GNSS/INS matching and synchronization */
            if ((fabs(GI_dt)-ins->dttol)<=(ins->nn*ins->interval)/2.0){
                if (NO==rtk->match) rtk->match=YES;
                rtk->upte=SYNC_YES;
                /* record the synchronization time */
                rtk->upte_time=ins->time; 
                /* if GNSS is available, set the nominal IMU update flag to 0 */
                rtk->nominal_upte=NO;               
                iimu+=nn; iobsu+=nu; ipos+=npos;
            }
            else if (GI_dt<0){
                if (NO==rtk->match) {iimu+=nn; return 0;}
                else iimu+=nn;
            }
            else if (GI_dt>0){
                if (NO==rtk->match) {iobsu+=nu; ipos+=npos; return 0;}   
                else { iobsu+=nu; ipos+=npos;}    
            }            
        } 
        else {
            iobsu+=nu; ipos+=npos;
        }     

        /* update sbas corrections */
        while (isbs<sbss.n) {
            time=gpst2time(sbss.msgs[isbs].week,sbss.msgs[isbs].tow);

            if (getbitu(sbss.msgs[isbs].msg,8,6)!=9) { /* except for geo nav */
                sbsupdatecorr(sbss.msgs+isbs,&navs);
            }
            if (timediff(time,obs[0].time)>-1.0-DTTOL) break;
            isbs++;
        }
        /* update rtcm ssr corrections */
        if (*rtcm_file) {
            update_rtcm_ssr(obs[0].time);
        }
    }
    /* input backward data */
    else if (SOLTYPE_BACKWARD==popt->reverse) 
    { 
        if (PMODE_LC_POS!=popt->mode) 
        {
            if ((nu=nextobsb(&obss,&iobsu,1))<=0) return -1;
            if (popt->intpref) {
                /* for interpolation, find first base timestamp before rover timestamp */
                for (;(nr=nextobsb(&obss,&iobsr,2))>0;iobsr-=nr)
                    if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)<DTTOL) break;
            }
            else {
                /* if not interpolating, find closest timestamp */
                dt=iobsr>=0?timediff(obss.data[iobsr].time,obss.data[iobsu].time):0;
                for (i=iobsr;(nr=nextobsb(&obss,&i,2))>0;iobsr=i,i-=nr) {
                    dt_next=timediff(obss.data[i].time,obss.data[iobsu].time);
                    if (fabs(dt_next)>fabs(dt)) break;
                    dt=dt_next;
                }
            }
            nr=nextobsb(&obss,&iobsr,2);

            /* NOTE: store the observations of rover and base in the obs structure in order */
            for (i=0;i<nu&&n<MAXOBS*2;i++) obs[n++]=obss.data[(iobsu-nu+1)+i];
            for (i=0;i<nr&&n<MAXOBS*2;i++) obs[n++]=obss.data[(iobsr-nr+1)+i];

            gnss_time=obss.data[iobsu].time;
        }
        else {
            pos2sol(poss,&rtk->sol,ipos);
            gnss_time=rtk->sol.time;            
        }

        /* NOTE: GNSS/INS time synchronization */
        if (GINS_LC==popt->GI_mode||GINS_TC==popt->GI_mode||GINS_STC==popt->GI_mode)
        {
            if (iimu<0) return -1;

            /* calculate the difference between the current time and the nominal measurement update time */
            sec=imus.data[iimu-nn].time.sec;
            rtk->nominal_upte=NO;
            ndt=fabs(sec-round((sec+ins->interval/2.0)/rtk->interval)*rtk->interval);
            /* INS navigation time */
            ins->time=imus.data[iimu-nn].time;            

            /* the second condition is used to process the IMU time on both sides being 0.5 IMU sampling interval away from GNSS time */
            if ((fabs(ndt)-ins->dttol)<=(ins->nn*ins->interval)/2.0
                &&(fabs(timediff(ins->time,rtk->upte_time))>=(rtk->interval-ins->interval))) {
                rtk->nominal_upte=YES;
                /* record the synchronization time */
                rtk->upte_time=ins->time;                
            }  

            /* calculate the difference between the current IMU and GNSS observation time */
            GI_dt=timediff(ins->time,gnss_time);
            rtk->upte=SYNC_NO;   
            /* for backward processing mode, the sign of the INS velocity and gyroscope bias is inverted */
            if (!reverse_flag) pos_reverse(popt,ins,&reverse_flag);
            imucpy(popt,imu,imus,iimu,ins->nn); 

            /* GNSS/INS matching and synchronization */
            if ((fabs(GI_dt)-ins->dttol)<=(ins->nn*ins->interval)/2.0){
                if (NO==rtk->match) rtk->match=YES;
                rtk->upte=SYNC_YES;
                /* record the synchronization time */
                rtk->upte_time=ins->time; 
                /* if GNSS is available, set the nominal IMU update flag to 0 */
                rtk->nominal_upte=NO;               
                iimu-=nn; iobsu-=nu; ipos-=npos;
            }
            else if (GI_dt>0){
                if (NO==rtk->match) {iimu-=nn; return 0;}
                else iimu-=nn;
            }
            else if (GI_dt<0){
                if (NO==rtk->match) {iobsu-=nu; ipos-=npos; return 0;}   
                else { iobsu-=nu; ipos-=npos; }    
            }            
        } 
        else {
            iobsu-=nu; ipos-=npos;
        } 

        /* update sbas corrections */
        while (isbs>=0) {
            time=gpst2time(sbss.msgs[isbs].week,sbss.msgs[isbs].tow);

            if (getbitu(sbss.msgs[isbs].msg,8,6)!=9) { /* except for geo nav */
                sbsupdatecorr(sbss.msgs+isbs,&navs);
            }
            if (timediff(time,obs[0].time)<1.0+DTTOL) break;
            isbs--;
        }
    }
    return (PMODE_LC_POS!=popt->mode)?n:npos;
}
/* output to file message of invalid time mark -------------------------------*/
static void outinvalidtm(FILE *fptm, const solopt_t *opt, const gtime_t tm)
{
    gtime_t time = tm;
    double gpst;
    const double secondsInAWeek = 604800;
    int week,timeu;
    char s[100];

    timeu=opt->timeu<0?0:(opt->timeu>20?20:opt->timeu);

    if (opt->times>=TIMES_UTC) time=gpst2utc(time);
    if (opt->times==TIMES_JST) time=timeadd(time,9*3600.0);

    if (opt->timef) time2str(time,s,timeu);
    else {
        gpst=time2gpst(time,&week);
        if (secondsInAWeek-gpst < 0.5/pow(10.0,timeu)) {
            week++;
            gpst=0.0;
        }
        sprintf(s,"%4d   %*.*f",week,6+(timeu<=0?0:timeu+1),timeu,gpst);
    }
    strcat(s, "   Q=0, Time mark is not valid\n");

    fwrite(s,strlen(s),1,fptm);
}
/* fill structure sol_t for time mark ----------------------------------------*/
static sol_t fillsoltm(const sol_t solold, const sol_t solnew, const gtime_t tm)
{
    gtime_t t1={0},t2={0};
    sol_t sol=solold;
    int i=0;

    if (solold.stat == 0 || solnew.stat == 0) {
        sol.stat = 0;
    } else {
        sol.stat = (solold.stat > solnew.stat) ? solold.stat : solnew.stat;
    }
    sol.ns = (solold.ns < solnew.ns) ? solold.ns : solnew.ns;
    sol.ratio = (solold.ratio < solnew.ratio) ? solold.ratio : solnew.ratio;

    /* interpolation position and speed of time mark */
    t1 = solold.time;
    t2 = solnew.time;
    sol.time = tm;

    for (i=0;i<6;i++)
    {
        sol.rr[i] = solold.rr[i] + timediff(tm,t1) / timediff(t2,t1) * (solnew.rr[i] - solold.rr[i]);
    }

    return sol;
}

/* carrier-phase bias correction by ssr --------------------------------------*/
static void corr_phase_bias_ssr(obsd_t *obs, int n, const nav_t *nav)
{
    double freq;
    uint8_t code;
    int i,j;

    for (i=0;i<n;i++) for (j=0;j<NFREQ;j++) {
        code=obs[i].code[j];

        if ((freq=sat2freq(obs[i].sat,code,nav))==0.0) continue;

        /* correct phase bias (cyc) */
        obs[i].L[j]-=nav->ssr[obs[i].sat-1].pbias[code-1]*freq/CLIGHT;
    }
}

/* determine the position of the current reference station (vrs mode) */
static int vrs_pos(prcopt_t *popt, const obsd_t *obs, vrs_t *vrs)
{
    int i;
    gtime_t obs_time=obs[0].time;

    for (i=vrs->idx;i<vrs->nbase;i++) {
        if (timediff(vrs->time[i],obs_time)<=0&&timediff(vrs->time[i+1],obs_time)>0) break;
        /*the observation times of the base station and the rover station are not yet aligned */
        else if (timediff(vrs->time[i],obs_time)>0) return 0; 
    }
    vrs->idx=i;
    matcpy(popt->rb,vrs->pos[vrs->idx],3,1);

    return 1;
}

/* GNSS-assisted detection INS status */
static int gnss_aid_ins(rtk_t *rtk, const int stat, const double *rr)
{
    prcopt_t *opt=&rtk->opt;
    double rr_[3],ins_pos[3],dpos[3],thres_ins_ouj=50.0;
    int i;

    /* GNSS-assisted detection INS status count, used for INS reinitialization */
    if (stat&&rtk->align&&(GINS_LC==opt->GI_mode||GINS_TC==opt->GI_mode||GINS_STC==opt->GI_mode))
    {
        if (rtk->outage<MAX_OUTIME&&!outsim.valid_flag) {
            ins2gnss(opt,&rtk->ins,rr_,3);
            pos2ecef(rr_,ins_pos);
            for (i=0;i<3;i++) dpos[i]=ins_pos[i]-rr[i];
            if (norm(dpos,3)>thres_ins_ouj) {
                rtk->gnss_aid_age++;
            }
            else { rtk->gnss_aid_age=0; return 0; }       
        }
        else {
            rtk->gnss_aid_age=0;
            return 0;
        }
    }

    /* check GNSS-assisted INS status */
    if (rtk->gnss_aid_age>MAX_GNSS_AID_AGE&&!outsim.valid_flag) {
        rtk->outage+=(MAX_OUTIME+1); /* trigger INS reinitialization */
        rtk->gnss_aid_age=0;         /* reset GNSS-assisted INS status count */
        trace(7,"warning: The GNSS and INS positions differ too much!\n");
        return 0;
    }

    if (!stat) {
        rtk->gnss_aid_age=0;
        return 0;
    }
}

/* process positioning -------------------------------------------------------*/
static void procpos(FILE *fp, prcopt_t *popt, const solopt_t *sopt, rtk_t *rtk, int mode)
{
    sol_t sol={{0}},oldsol={{0}},newsol={{0}};
    rtk_t *rtk_tdcp=(rtk_t *)malloc(sizeof(rtk_t));   /* for tdcp module */
    obsd_t *obs=(obsd_t *)malloc(sizeof(obsd_t)*MAXOBS*2);     /* observations at the current epoch for rover and base */
    obsd_t *obs_old=(obsd_t *)malloc(sizeof(obsd_t)*MAXOBS*2); /* observations at the previous epoch for rover and base */
    imud_t *imu=(imud_t *)malloc(sizeof(imud_t)*MAXINS);
    double rb[3]={0};
    int i,nobs,num=0,align,vel_flag=0,stat;
    int n=0,n_old=0,nr=0,nr_old=0;

    trace(3,"procpos : mode=%d\n",mode); /* 0=forward or backward, 1=forward and backward smoothing */
    
    rtcm_path[0]='\0';
    vrs.idx=(PMODE_FIXED==popt->mode)?1:0; /* init vrs index */

    /* initialize rtk_tdcp */
    if (GINS_OFF!=popt->GI_mode) rtkinit(rtk_tdcp,popt,NULL);

    /* initialize GNSS sampling interval */
    if (!rtk->interval) gnss_intervel(rtk,&obss,&poss);

    /* epoch-by-epoch processing */
    while ((nobs=inputobs(rtk,obs,imu,stat,popt))>=0) {    
        /* DebugGlo initialization */
        if (GINS_OFF!=popt->GI_mode) Debug_Glo.tNow=rtk->ins.time; 
        else Debug_Glo.tNow=obs[0].time;           
        Debug_Glo=DebugGlo_init(Debug_Glo);     
        DebugTime(rtk,Debug_Glo.tNow,118221,2362); 

        /* determine the position of the current reference station (vrs mode) */
        if (PMODE_DGPS<=popt->mode&&PMODE_FIXED>=popt->mode&&STA_VRS==popt->statype) vrs_pos(&rtk->opt,obs,&vrs);

        /* vehicle zero speed detection for ZUPT and ZIHR */
        if (popt->constraint[1]||popt->constraint[2]) zerovel_detect(rtk,imu);

        if (PMODE_LC_POS==popt->mode) n=nobs;
        else { /* exclude satellites */
            for (i=n=0;i<nobs;i++) {
                if ((satsys(obs[i].sat,NULL)&popt->navsys)&&popt->exsats[obs[i].sat-1]!=1) obs[n++]=obs[i];
            }            
        }
       
        /* if no satellites are available in GNSS mode or initial alignment fails in GNSS/INS mode, exit current epoch processing */
        if (GINS_OFF==popt->GI_mode&&n<=0) continue;
        else if (GINS_OFF!=popt->GI_mode&&n<=0&&!rtk->align) continue;

        /* the sign of Doppler observations is determined based on pseudorange variation between adjacent epochs */
        if ((GINS_OFF==popt->GI_mode||SYNC_YES==rtk->upte)&&PMODE_LC_POS!=popt->mode&&!rtk->dopsgn) dopple_sgn(rtk,obs,obs_old,n,n_old);

        /* ins initial alignment */
        if (GINS_LC==popt->GI_mode||GINS_TC==popt->GI_mode||GINS_STC==popt->GI_mode){
            if (NO==rtk->match) continue; 
            /* TDCP estimated velocity */
            if (SYNC_YES==rtk->upte) {
                vel_flag=0; /* reset vel flag */
                /* determine the number of satellites of rover in the current epoch and the previous epoch */
                for (i=nr=0;i<n;i++)         if (obs[i].rcv==1) nr++;
                for (i=nr_old=0;i<n_old;i++) if (obs_old[i].rcv==1) nr_old++;
                /* initialize rtk_tdcp parameters */
                rtk_tdcp->interval=rtk->interval; rtk_tdcp->dopsgn=rtk->dopsgn;
                /* multi-strategy velocity estimation (TDCP/dopple/position difference) */
                if (nr_old&&nr) vel_flag=tdcp_vel(rtk_tdcp,rtk->align,obs,obs_old,nr,nr_old,&navs,popt);
                if (norm(rtk_tdcp->sol.rr+3,3)>0.0) matcpy(rtk->sol.rr+3,rtk_tdcp->sol.rr+3,3,1); /* copy TDCP estimated velocity to rtk struct */
                /* GNSS-assisted detection INS status */  
                gnss_aid_ins(rtk,rtk_tdcp->sol.stat,rtk_tdcp->sol.rr);     
                /* save the GNSS observations of the previous epoch */              
                n_old=n; 
                for (i=0;i<n;i++) obs_old[i]=obs[i];                 
            }  
            /* velocity vector assisted alignment */  
            if (!rtk->align||rtk->outage>MAX_OUTIME) {
                rtk->align=ins_align(rtk,obs,n,&navs,popt,vel_flag);
            }    
            if (!rtk->align) continue;  
        }

        /* INS mechanization and GNSS/INS time update */
        if (GINS_LC==popt->GI_mode||GINS_TC==popt->GI_mode||GINS_STC==popt->GI_mode){
            ins_mech(&rtk->ins,imu,popt);            
            ins_update(rtk);           
            if (SYNC_NO==rtk->upte&&NO==rtk->nominal_upte) continue;
        }

        /* for GNSS/INS integration navigation, when GNSS is not available, use motion constraints to assist */
        if (n<=0&&GINS_OFF!=popt->GI_mode&&(popt->constraint[0]||popt->constraint[1]||popt->constraint[2])) {
            motion_constraints(rtk,popt);
        }

        /* GNSS outage simulation */
        if ((outsim.valid_flag=isoutage(rtk,Debug_Glo.tNow,outsim))||YES==rtk->nominal_upte||0==n) {
            rtk->outage++;
            if (GINS_LC==popt->GI_mode||GINS_STC==popt->GI_mode) {
                rtk->lcgins.sol.stat=SOLQ_INS;
                update_instat(&rtk->opt,&rtk->ins,rtk->lcgins.P,&rtk->lcgins.sol,rtk->ins.nx);
                outsol(fp,&rtk->lcgins.sol,rtk->lcgins.sol.rr,popt,sopt);                
            }
            else if (GINS_TC==popt->GI_mode){
                rtk->sol.stat=SOLQ_INS;
                update_instat(&rtk->opt,&rtk->ins,rtk->P,&rtk->sol,rtk->nx); 
                outsol(fp,&rtk->sol,rtk->rb,popt,sopt);                 
            }
            continue;
        }

        /* GNSS/INS tightly coupled integration */
        if (PMODE_LC_POS!=popt->mode) {
            /* carrier-phase bias correction */
            if (!strstr(popt->pppopt,"-ENA_FCB")) {
                corr_phase_bias_ssr(obs,n,&navs);
            }

            /* multipath correction for BDS2 */
            if (popt->navsys&SYS_CMP) {
                BDmulCorr(rtk,obs,n); 
            }

            /* navigation processing */
            if (!rtkpos(rtk,obs,n,&navs)&&n>0) {
                /* whether to output a blank line when GNSS is outage  */
                /* if (GINS_OFF==popt->GI_mode) continue; */
            }            
        }

        /* GNSS/INS loosely coupled/semi-tight coupled integration */
        if ((GINS_LC==popt->GI_mode||GINS_STC==popt->GI_mode||PMODE_LC_POS==popt->mode)&&n>0){
            lc_gins(rtk);            
        }

        /* forward or backward mode */
        if (mode==SOLMODE_SINGLE_DIR) {    
            /* save the GNSS observations of the previous epoch */
            n_old=n; 
            for (i=0;i<n;i++) obs_old[i]=obs[i]; 
            /* output the GNSS/INS solution */
            if (GINS_LC==popt->GI_mode||GINS_STC==popt->GI_mode) outsol(fp,&rtk->lcgins.sol,rtk->lcgins.sol.rr,popt,sopt);
            else outsol(fp,&rtk->sol,rtk->rb,popt,sopt);
            rtk->sol.iFlag=SOLF_GNSS; /* reset solution flag */
            oldsol=rtk->sol;
        }
        /* combined-forward */
        else if (SOLTYPE_FORWARD==popt->reverse) { 
            if (isolf>=nepoch) {
                free(obs); return;
            }
            solf[isolf]=rtk->sol;
            for (i=0;i<3;i++) rbf[i+isolf*3]=rtk->rb[i];
            isolf++;
        }
        /* combined-backward */
        else if (SOLTYPE_BACKWARD==popt->reverse) { 
            if (isolb>=nepoch) {
                free(obs); return;
            }
            solb[isolb]=rtk->sol;
            for (i=0;i<3;i++) rbb[i+isolb*3]=rtk->rb[i];
            isolb++;
        }
    }

    /* obs and obs_old point to the same memory and only need to free once */
    rtkfree(rtk_tdcp);
    free(obs); free(imu); 
    obs=obs_old=NULL; /* free obs_old to avoid memory leak */
}
/* validation of combined solutions ------------------------------------------*/
static int valcomb(const sol_t *solf, const sol_t *solb, double *rbf, double *rbb, const prcopt_t *popt)
{
    double dr[3],var[3];
    int i;
    char tstr[32];

    trace(8,"valcomb :\n");

    /* compare forward and backward solution */
    for (i=0;i<3;i++) {
        dr[i]=solf->rr[i]-solb->rr[i];
        if (popt->mode==PMODE_MOVEB) dr[i]-=(rbf[i]-rbb[i]);
        var[i]=(double)solf->qr[i]+(double)solb->qr[i];
    }
    for (i=0;i<3;i++) {
        if (dr[i]*dr[i]<=16.0*var[i]) continue; /* ok if in 4-sigma */

        time2str(solf->time,tstr,2);
        trace(7,"Forward and backward smoothing (FBS) quality check failed: %s dr=%.3f %.3f %.3f std=%.3f %.3f %.3f\n",
              tstr+11,dr[0],dr[1],dr[2],SQRT(var[0]),SQRT(var[1]),SQRT(var[2]));
        return 0;
    }
    return 1;
}
/* combine forward/backward solutions and save results ---------------------*/
static void combres(FILE *fp, rtk_t *rtk, const prcopt_t *popt, const solopt_t *sopt)
{
    sol_t sols={{0}},sol={{0}},oldsol={{0}},newsol={{0}};
    double tt,Qf[9],Qf2[9],Qb[9],Qb2[9],Qs[9],Qs2[9],rbs[3]={0},rb[3]={0},rr_f[3],rr_b[3],rr_s[3];
    int i,j,k,num=0;
    int pri[]={7,1,2,3,4,5,1,6}; /* no:0,fix:1,float:2,sbas:3,dgps:4,single:5,ppp:6,ins:7,cons:8 */

    trace(3,"Forward and backward smoothing (FBS) : isolf=%d isolb=%d\n",isolf,isolb);


    /* set reference station position */
    for (i=0,j=isolb-1;i<isolf&&j>=0;i++,j--) {
        /* time debug */
        Debug_Glo.tNow=solf[i].time;           
        Debug_Glo=DebugGlo_init(Debug_Glo);
        DebugTime(rtk,solf[i].time,437459,2188);

        if ((tt=timediff(solf[i].time,solb[j].time))<-DTTOL) {
            sols=solf[i];
            for (k=0;k<3;k++) rbs[k]=rbf[k+i*3];
            j++;
        }
        else if (tt>DTTOL) {
            sols=solb[j];
            for (k=0;k<3;k++) rbs[k]=rbb[k+j*3];
            i--;
        }
        /* prioritize using reference station position with high solution quality */
        else if (pri[solf[i].stat]<pri[solb[j].stat]) {
            sols=solf[i];
            for (k=0;k<3;k++) rbs[k]=rbf[k+i*3];
        }
        else if (pri[solf[i].stat]>pri[solb[j].stat]) {
            sols=solb[j];
            for (k=0;k<3;k++) rbs[k]=rbb[k+j*3];
        }
        else {
            sols=solf[i];
            sols.time=timeadd(sols.time,-tt/2.0);

            /* forward and backward smoothing (FBS) quality check (only PPK mode)*/
            if ((popt->mode==PMODE_KINEMA||popt->mode==PMODE_MOVEB)&&sols.stat==SOLQ_FIX) {
                /* degrade fix to float if validation failed */
                if (!valcomb(solf+i,solb+j,rbf+i*3,rbb+j*3,popt)) sols.stat=SOLQ_FLOAT;
            }

            /* solution to position covariance matrix */
            soltocov(solf+i,Qf);
            soltocov(solb+j,Qb);

            /* NOTE: smoother for position solution */
            if (popt->mode==PMODE_MOVEB) { /* for the moving baseline mode, the baseline length is smoothed */
                for (k=0;k<3;k++) rr_f[k]=solf[i].rr[k]-rbf[k+i*3];
                for (k=0;k<3;k++) rr_b[k]=solb[j].rr[k]-rbb[k+j*3];
                if (smoother(rr_f,Qf,rr_b,Qb,3,rr_s,Qs)) continue;
                for (k=0;k<3;k++) sols.rr[k]=rbs[k]+rr_s[k];
            }
            else {
                if (smoother(solf[i].rr,Qf,solb[j].rr,Qb,3,sols.rr,Qs)) continue;
            }
            /* position covariance matrix to solution */
            covtosol(Qs,&sols);

            /* NOTE: smoother for velocity solution */
            if (popt->dynamics||GINS_OFF!=popt->GI_mode) {
                /* solution to velocity covariance matrix */
                soltocov_vel(solf+i,Qf);
                soltocov_vel(solb+j,Qb);
                if (GINS_LC==popt->GI_mode||GINS_STC==popt->GI_mode||GINS_TC==popt->GI_mode) {
                    if (smoother(solf[i].vel,Qf,solb[j].vel,Qb,3,sols.vel,Qs)) continue;
                }
                else {
                    if (smoother(solf[i].rr+3,Qf,solb[j].rr+3,Qb,3,sols.rr+3,Qs)) continue;                    
                }

                /* velocity covariance matrix to solution */
                covtosol_vel(Qs,&sols);
            }

            /* NOTE: smoother for attitude and bias solution */
            if (GINS_LC==popt->GI_mode||GINS_STC==popt->GI_mode||GINS_TC==popt->GI_mode) {
                /* solution to attitude covariance matrix */
                soltocov_att(solf+i,Qf);
                soltocov_att(solb+j,Qb);
                
                if (smoother_att(solf[i].qnb,Qf,solb[j].qnb,Qb,3,sols.att,Qs)) continue;

                /* velocity covariance matrix to solution */
                covtosol_att(Qs,&sols);

                /* solution to bias covariance matrix */
                soltocov_bga(solf+i,Qf,Qf2);
                soltocov_bga(solb+j,Qb,Qb2);

                /* convert bias unit from deg/h and ug to rad/s and g */
                for (k=0;k<3;k++) {
                    solf[i].bg[k]/=R2D*3600;
                    solb[j].bg[k]/=R2D*3600;
                    solf[i].ba[k]/=1E5;
                    solb[j].ba[k]/=1E5;
                }

                /* smoother for gyroscope and accelerometer bias solution */
                if (smoother(solf[i].bg,Qf ,solb[j].bg,Qb ,3,sols.bg,Qs)) continue;
                if (smoother(solf[i].ba,Qf2,solb[j].ba,Qb2,3,sols.ba,Qs2)) continue;

                /* convert bias unit from rad/s and g to deg/h and ug */
                for (k=0;k<3;k++) {
                    sols.bg[k]*=R2D*3600;
                    sols.ba[k]*=1E5;
                }

                /* bias covariance matrix to solution */
                covtosol_bga(Qs,Qs2,&sols);
            }
        }
        /* output the GNSS/INS solution */
        outsol(fp,&sols,rbs,popt,sopt);
        oldsol=sols;
    }
}
/* read prec ephemeris, sbas data, tec grid and open rtcm --------------------*/
static int readpreceph(const filopt_t *fopt, const prcopt_t *prcopt, nav_t *nav, sbs_t *sbs)
{
    seph_t seph0={0};
    int i;
    const char *ext;

    nav->ne=nav->nemax=0;
    nav->nc=nav->ncmax=0;
    sbs->n =sbs->nmax =0;

    /* read precise ephemeris files */
    if (EPHOPT_PREC==prcopt->sateph) {
        if (*fopt->sp3&&(ext=strrchr(fopt->sp3,'.'))) {   
            if (!strstr(ext,".sp3")&&!strstr(ext,".SP3")&&
                !strstr(ext,".eph")&&!strstr(ext,".EPH")) return 0;
            readsp3(fopt->sp3,nav,0);
        }
        else {
            trace(1,"missing precise ephemeris file!\n");
            return 0;
        }
        /* read precise clock files */
        if (*fopt->clk&&(ext=strrchr(fopt->clk,'.'))) {
            if (!strstr(ext,".clk")&&!strstr(ext,".CLK")) return 0;
            readrnxc(fopt->clk,nav);
        }  
        else {
            trace(1,"missing precise clock file!\n");
            return 0;            
        }      
    }

    /* read sbas message files */
    if (*fopt->sbs&&(ext=strrchr(fopt->sbs,'.'))) {     
        if (strstr(ext,".sbs")&&strstr(ext,".SBS")&&
            strstr(ext,".ems")&&strstr(ext,".EMS")) return 0;
        sbsreadmsg(fopt->sbs,prcopt->sbassatsel,sbs);
    }

    /* allocate sbas ephemeris */
    nav->ns=nav->nsmax=NSATSBS*2;
    if (!(nav->seph=(seph_t *)malloc(sizeof(seph_t)*nav->ns))) {
         showmsg("error : sbas ephem memory allocation");
         trace(1,"error : sbas ephem memory allocation");
         return 0;
    }
    for (i=0;i<nav->ns;i++) nav->seph[i]=seph0;

    /* set rtcm file and initialize rtcm struct */
    /* rtcm_file[0]=rtcm_path[0]='\0'; fp_rtcm=NULL;

    for (i=0;i<n;i++) {
        if ((ext=strrchr(infile[i],'.'))&&
            (!strcmp(ext,".rtcm3")||!strcmp(ext,".RTCM3"))) {
            strcpy(rtcm_file,infile[i]);
            init_rtcm(&rtcm);
            break;
        }
    } */
   return 1;
}
/* free prec ephemeris and sbas data -----------------------------------------*/
static void freepreceph(nav_t *nav, sbs_t *sbs)
{
    int i;

    trace(3,"freepreceph:\n");

    free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
    free(nav->pclk); nav->pclk=NULL; nav->nc=nav->ncmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
    free(sbs->msgs); sbs->msgs=NULL; sbs->n =sbs->nmax =0;
    for (i=0;i<nav->nt;i++) {
        free(nav->tec[i].data);
        free(nav->tec[i].rms );
    }
    free(nav->tec ); nav->tec =NULL; nav->nt=nav->ntmax=0;

    if (fp_rtcm) fclose(fp_rtcm);
    free_rtcm(&rtcm);
}
/* read obs and nav data -----------------------------------------------------*/
static int readobsnav(gtime_t ts, gtime_t te, double ti, const filopt_t *fopt, const prcopt_t *prcopt,
                      obs_t *obs, nav_t *nav, sta_t *sta)
{
    int i,j,ind=0,nobs=0,flag=1; /* flag: only precise ephemeris and precise clock*/
    const char *ext;

    trace(3,"readobsnav: ts=%s\n",time_str(ts,0));

    obs->data=NULL; obs->n =obs->nmax =0;
    nav->eph =NULL; nav->n =nav->nmax =0;
    nav->geph=NULL; nav->ng=nav->ngmax=0;
    /* free(nav->seph); */ /* is this needed to avoid memory leak??? */
    nav->seph=NULL; nav->ns=nav->nsmax=0;
    nepoch=0;

    /* read user rinex obs file */
    if (*fopt->obs_u&&(ext=strrchr(fopt->obs_u,'.'))) {
        if (ext[1]=='o'||ext[1]=='O'||ext[3]=='o'||ext[3]=='O')
        if (readrnxt(fopt->obs_u,1,ts,te,ti,prcopt,obs,nav,sta)<0) {
            trace(1,"error : insufficient memory of user observation file!\n");
            return 0;
        }
    }
    else {
        trace(1,"missing user observation file!\n");
        return 0;
    }

    /* read base station rinex obs file */
    if (PMODE_DGPS<=prcopt->mode&&PMODE_MOVEB>=prcopt->mode) {
        if (*fopt->obs_b&&(ext=strrchr(fopt->obs_b,'.'))) {
            if (ext[1]=='o'||ext[1]=='O'||ext[3]=='o'||ext[3]=='O')
            if (readrnxt(fopt->obs_b,2,ts,te,ti,prcopt,obs,nav,sta+1)<0) {
                trace(1,"error : insufficient memory of base observation file!\n");
                return 0;
            }
        }  
        else {
            trace(1,"differential mode but lacks base observation file!\n");
            return 0;
        }      
    }
   
    /* read nav file */
    if (*fopt->nav&&(ext=strrchr(fopt->nav,'.'))) {
        if (ext[1]=='p'||ext[1]=='P'||ext[1]=='n'||ext[3]=='p'||ext[3]=='P'||ext[3]=='n')
        if (readrnxt(fopt->nav,3,ts,te,ti,prcopt,obs,nav,NULL)<0) {
            trace(1,"error : insufficient memory of navigation file!\n");
            return 0;
        }
    } 
    else {
        trace(1,"missing navigation file!\n");
        return 0;
    }

    if (obs->n<=0) {
        trace(1,"error : no observation data!\n");
        return 0;
    }
    /* if have precise ephemeris and precise clock, no broadcast ephemeris is OK*/
    if (nav->n<=0&&nav->ng<=0&&nav->ns<=0&&nav->ne<=0&&nav->nc<=0) {
        trace(1,"error : no navigation data!\n");
        return 0;
    }
    else if (nav->n<=0&&nav->ng<=0&&nav->ns<=0) {
        flag=0;
        trace(1,"warning : missing navigation data but with precise orbit and clock data!\n");
    }
    /* sort observation data */
    nepoch=sortobs(obs);

    /* delete duplicated ephemeris */
    if (flag) {
       uniqnav(nav);
    }   

    /* set time span for progress display */
    if (ts.time==0||te.time==0) {
        for (i=0;   i<obs->n;i++) if (obs->data[i].rcv==1) break;
        for (j=obs->n-1;j>=0;j--) if (obs->data[j].rcv==1) break;
        if (i<j) {
            if (ts.time==0) ts=obs->data[i].time;
            if (te.time==0) te=obs->data[j].time;
            settspan(ts,te);
        }
    }
    return 1;
}
/* free obs and nav data -----------------------------------------------------*/
static void freeobsnav(obs_t *obs, nav_t *nav)
{
    trace(3,"freeobsnav:\n");

    free(obs->data); obs->data=NULL; obs->n =obs->nmax =0;
    free(nav->eph ); nav->eph =NULL; nav->n =nav->nmax =0;
    free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
}
/* average of single position ------------------------------------------------*/
static int avepos(double *ra, int rcv, const obs_t *obs, const nav_t *nav,
                  const prcopt_t *opt)
{
    obsd_t data[MAXOBS];
    gtime_t ts={0};
    sol_t sol={{0}};
    int i,j,n=0,m,iobs;

    trace(3,"avepos: rcv=%d obs.n=%d\n",rcv,obs->n);

    for (i=0;i<3;i++) ra[i]=0.0;

    for (iobs=0;(m=nextobsf(obs,&iobs,rcv))>0;iobs+=m) {

        for (i=j=0;i<m&&i<MAXOBS;i++) {
            data[j]=obs->data[iobs+i];
            if ((satsys(data[j].sat,NULL)&opt->navsys)&&
                opt->exsats[data[j].sat-1]!=1) j++;
        }
        if (j<=0||!screent(data[0].time,ts,ts,1.0)) continue; /* only 1 hz */

        if (!pntpos(NULL,data,j,nav,opt,&sol,NULL,NULL)) continue;

        for (i=0;i<3;i++) ra[i]+=sol.rr[i];
        n++;
    }
    if (n<=0) {
        trace(1,"no average of base station position\n");
        return 0;
    }
    for (i=0;i<3;i++) ra[i]/=n;
    return 1;
}
/* station position from file ------------------------------------------------*/
static int getstapos(const char *file, const char *name, vrs_t *vrs, int idx_sta)
{
    FILE *fp;
    char buff[256],sname[256],*p;
    const char *q;
    double pos[3],llh[3],data[7];
    int flag=-1;

    trace(3,"getstapos: file=%s name=%s\n",file,name);

    if (!(fp=fopen(file,"r"))) {
        trace(1,"station position file open error: %s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        if ((p=strchr(buff,'%'))) *p='\0';
        
        /* geodetic position (lat(dd mm ss.ss),lon(dd mm ss.ss),height(m)) */
        if (sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %255s",data,data+1,data+2,data+3,data+4,data+5,data+6,sname)==8) flag=3; 
        /* ecef/llh(lat(dd.dd),lon(dd.dd),height(m)) position */
        else if (sscanf(buff,"%lf %lf %lf %255s",pos,pos+1,pos+2,sname)==4) flag=0;  
        else continue;
        
        for (p=sname,q=name;*p&&*q;p++,q++) {
            if (toupper((int)*p)!=toupper((int)*q)) break;
        }
        if (!*p) {
            if (0==flag) {
                if (fabs(pos[0])<90&&fabs(pos[1])<180) flag=2; /* llh */
                else flag=1; /* ecef */ 
            } 

            if (1==flag) {
                vrs->pos[idx_sta][0]=pos[0];
                vrs->pos[idx_sta][1]=pos[1];
                vrs->pos[idx_sta][2]=pos[2];
                vrs->nbase++;                
            }
            else if (2==flag) {
                llh[0]=pos[0]*D2R;
                llh[1]=pos[1]*D2R;
                llh[2]=pos[2];
                pos2ecef(llh,vrs->pos[idx_sta]);
                vrs->nbase++;
            }
            else if (3==flag) {
                llh[0]=data[0]*D2R+data[1]*D2R/60.0+data[2]*D2R/3600.0;
                llh[1]=data[3]*D2R+data[4]*D2R/60.0+data[5]*D2R/3600.0;
                llh[2]=data[6];
                pos2ecef(llh,vrs->pos[idx_sta]);
                vrs->nbase++;
            }

            fclose(fp);
            return 1;
        }
    }
    fclose(fp);
    trace(1,"no station position: %s %s\n",name,file);
    return 0;
}
/* antenna phase center position ---------------------------------------------*/
static int antpos(prcopt_t *opt, int rcvno, const obs_t *obs, const nav_t *nav,
                  const sta_t *sta, const char *posfile)
{
    double *rr=rcvno==1?opt->ru:opt->rb,del[3],pos[3],dr[3]={0};
    int i=0,j,k=(PMODE_FIXED==opt->mode&&rcvno==2)?1:0,postype=rcvno==1?opt->rovpos:opt->refpos;
    char *name;

    trace(3,"antpos  : rcvno=%d\n",rcvno);

    if (postype==POSOPT_SINGLE) { /* average of single position ,don't support for vrs mode */
        if (!avepos(rr,rcvno,obs,nav,opt)) {
            showmsg("error : station pos computation");
            return 0;
        }
    }
    else if (postype==POSOPT_FILE) { /* read from position file */
        while (strlen(stas[rcvno==1?0:(i+1)].name)>0) {
            name=stas[rcvno==1?0:(i+1)].name;
            if (!getstapos(posfile,name,&vrs,k)) {
                showmsg("error : no position of %s in %s",name,posfile);
                return 0;
            }
            vrs.time[k]=stas[i+1].time;
            i++;k++;                 
        }     
        for (j=0;j<3;j++) rr[j]=vrs.pos[0][j]; /* init pos for first station */
    }
    else if (postype==POSOPT_RINEX) { /* get from rinex header */
        while (strlen(stas[rcvno==1?0:(i+1)].name)>0) {
            if (norm(stas[rcvno==1?0:(i+1)].pos,3)<=0.0) {
                showmsg("error : no position in rinex header");
                trace(1,"no position in rinex header\n");
                return 0;
            }
            /* add antenna delta unless already done in antpcv() */
            if (!strcmp(opt->anttype[rcvno],"*")) {
                if (stas[rcvno==1?0:(i+1)].deltype==0) { /* enu */
                    for (j=0;j<3;j++) del[j]=stas[rcvno==1?0:(i+1)].del[j];
                    del[2]+=stas[rcvno==1?0:(i+1)].hgt;
                    ecef2pos(stas[rcvno==1?0:(i+1)].pos,pos);
                    enu2ecef(pos,del,dr);
                }  else { /* xyz */
                    for (j=0;j<3;j++) dr[j]=stas[rcvno==1?0:(i+1)].del[j];
                }
            }
            for (j=0;j<3;j++) {
                vrs.pos[k][j]=stas[rcvno==1?0:(i+1)].pos[j]+dr[j]; 
                if (i==0) rr[j]=stas[rcvno==1?0:(i+1)].pos[j]+dr[j]; /* init pos for first station */ 
            } 
            vrs.time[k]=stas[i+1].time;
            vrs.nbase++;i++;k++; 
        }
    }
    return 1;
}
/* open processing session ----------------------------------------------------*/
static int openses(const prcopt_t *popt, const solopt_t *sopt,
                   const filopt_t *fopt, nav_t *nav, spcvs_t *pcvs, rpcvs_t *pcvr)
{
    trace(3,"openses :\n");

    /* read satellite and receiver antenna parameters */
    if (*fopt->antp&&!(readpcv(fopt->antp,pcvs,pcvr))) {
        showmsg("error : no sat ant pcv in %s",fopt->antp);
        trace(1,"sat antenna pcv read error: %s\n",fopt->antp);
        return 0;
    }
    /* open geoid data */
    if (sopt->geoid>0&&*fopt->geoid) {
        if (!opengeoid(sopt->geoid,fopt->geoid)) {
            showmsg("error : no geoid data %s",fopt->geoid);
            trace(2,"no geoid data %s\n",fopt->geoid);
        }
    }
    return 1;
}
/* close processing session ---------------------------------------------------*/
static void closeses(nav_t *nav, spcvs_t *pcvs, rpcvs_t *pcvr)
{
    trace(3,"closeses:\n");

    /* free antenna parameters */
    free(pcvs->pcv); pcvs->pcv=NULL; pcvs->n=pcvs->nmax=0;
    free(pcvr->pcv); pcvr->pcv=NULL; pcvr->n=pcvr->nmax=0;

    /* close geoid data */
    closegeoid();

    /* free erp data */
    free(nav->erp.data); nav->erp.data=NULL; nav->erp.n=nav->erp.nmax=0;

    /* close solution statistics and debug trace */
    rtkcloseoutfile();
    traceclose();
}
/* set antenna parameters ----------------------------------------------------*/
static void setpcv(gtime_t time, prcopt_t *popt, nav_t *nav, const spcvs_t *pcvs,
                   const rpcvs_t *pcvr, const sta_t *sta)
{
    spcv_t *spcv,spcv0={0};
    rpcv_t *rpcv,rpcv0={0};
    double pos[3],del[3];
    int i,j,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;
    char id[64];

    /* set satellite antenna parameters */
    for (i=0;i<MAXSAT;i++) {
        nav->spcvs[i]=spcv0;
        if (!(satsys(i+1,NULL)&popt->navsys)) continue;
        if (!(spcv=searchspcv(i+1,time,pcvs))) {
            satno2id(i+1,id);
            trace(7,"no satellite antenna pcv: %s\n",id);
            continue;
        }
        nav->spcvs[i]=*spcv;
    }
    /* set receiver antenna parameters */
    for (i=0;i<(mode?2:1);i++) {
        popt->pcvr[i]=rpcv0;
        if (!strcmp(popt->anttype[i],"*")) { /* set by station parameters */
            strcpy(popt->anttype[i],sta[i].antdes);
            if (sta[i].deltype==1) { /* xyz */
                if (norm(sta[i].pos,3)>0.0) {
                    ecef2pos(sta[i].pos,pos);
                    ecef2enu(pos,sta[i].del,del);
                    for (j=0;j<3;j++) popt->antdel[i][j]=del[j];
                }
            }
            else { /* enu */
                for (j=0;j<3;j++) popt->antdel[i][j]=stas[i].del[j];
            }
        }
        if (!(rpcv=searchrpcv(popt->anttype[i],time,pcvr))) {
            trace(7,"no receiver antenna pcv: %s,%s\n",sta[i].rectype,sta[i].antdes);
            *popt->anttype[i]='\0';
            continue;
        }
        strcpy(popt->anttype[i],rpcv->type);
        popt->pcvr[i]=*rpcv;
    }
}
/* read ocean tide loading parameters ----------------------------------------*/
static void readotl(prcopt_t *popt, const char *file, const sta_t *sta)
{
    int i,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;

    for (i=0;i<(mode?2:1);i++) {
        readblq(file,sta[i].name,popt->odisp[i]);
    }
}
/* open debug*/
static int opentrace(const prcopt_t *popt, const solopt_t *sopt, filopt_t *fopt)
{
    char tracefile[1024],current_dir[1024],path[1024];
    const char *p,*q,*sep=NULL;
    const char *s1[]={
        "SPP","PPD","PPK","PPS","GNSS-TC","Static-Start","Moving-Base","Fixed",
        "PPP_Kine","PPP_Static","PPP_Fixed","GNSS","",""
    };
    const char *s2[]={
        "F","B","FB"
    };

    /* set path separator */
    #ifdef WIN32
        sep = "\\";
    #else
        sep = "/";
    #endif

    /* get current directory */
    getcwd(current_dir, sizeof(current_dir));
    p=strrchr(current_dir, sep[0]);
    if (p) p++; /* move to the first character after the delimiter */

    q=strchr(p,'_');

    /* create the solution file(.pos) path */
    sprintf(path,"%s%s%s%s%s%s_%s_%s%s",current_dir,sep,"result",sep,(GINS_OFF==popt->GI_mode?"GNSS":fopt->ins_type),q,s1[popt->mode],s2[popt->soltype],
    (GINS_LC==popt->GI_mode?"_LC.pos":(GINS_TC==popt->GI_mode?"_TC.pos":(GINS_STC==popt->GI_mode?"_STC.pos":".pos")))); 
    
    strncpy(fopt->sol,path,1024);

     /* open debug trace */
    if (sopt->trace>0) {
        if (*fopt->sol) {
            strcpy(tracefile,fopt->sol);
            strcat(tracefile,".trace");
        }
        traceclose();
        traceopen(tracefile);
        tracelevel(sopt->trace);
    }

}
/* write header to output file -----------------------------------------------*/
static int outhead(const prcopt_t *popt, const solopt_t *sopt, filopt_t *fopt)
{
    FILE *fp=stdout;

    if (*fopt->sol) {
        createdir(fopt->sol);

        if (!(fp=fopen(fopt->sol,"wb"))) {
            showerr("error : open output file %s!",fopt->sol);
            return 0;
        }
    }
    /* output header */
    outheader(fp,popt,sopt,fopt);

    if (*fopt->sol) fclose(fp);

    return 1;
}
/* open output file for append -----------------------------------------------*/
static FILE *openfile(const char *outfile)
{
    trace(3,"openfile: outfile=%s\n",outfile);

    return !*outfile?stdout:fopen(outfile,"ab");
}
/* Name time marks file ------------------------------------------------------*/
static void namefiletm(char *outfiletm, const char *outfile)
{
    int i;

    for (i=(int)strlen(outfile);i>0;i--) {
        if (outfile[i] == '.') {
            break;
        }
    }
    /* if no file extension, then name time marks file as name of outfile + _events.pos */
    if (i == 0) {
        i = (int)strlen(outfile);
    }
    strncpy(outfiletm, outfile, i);
    strcat(outfiletm, "_events.pos");
}
/* post-processing positioning -------------------------------------------------
* post-processing positioning
* args   : gtime_t ts       I   processing start time (ts.time==0: no limit)
*        : gtime_t te       I   processing end time   (te.time==0: no limit)
*          double ti        I   processing interval  (s) (0:all)
*          double tu        I   processing unit time (s) (0:all)
*          prcopt_t *popt   I   processing options
*          solopt_t *sopt   I   solution options
*          filopt_t *fopt   I   file options
* return : status (0:ok,0>:error,1:aborted)
* notes  : input files should contain observation data, navigation data, precise
*          ephemeris/clock (optional), sbas log file (optional), ssr message
*          log file (optional) and tec grid file (optional). only the first
*          observation data file in the input files is recognized as the rover
*          data.
*
*          the type of an input file is recognized by the file extension as ]
*          follows:
*              .sp3,.SP3,.eph*,.EPH*: precise ephemeris (sp3c)
*              .sbs,.SBS,.ems,.EMS  : sbas message log files (rtklib or ems)
*              .rtcm3,.RTCM3        : ssr message log files (rtcm3)
*              .*i,.*I              : tec grid files (ionex)
*              others               : rinex obs, nav, gnav, hnav, qnav or clock
*
*          inputs files can include wild-cards (*). if an file includes
*          wild-cards, the wild-card expanded multiple files are used.
*
*          inputs files can include keywords. if an file includes keywords,
*          the keywords are replaced by date, time, rover id and base station
*          id and multiple session analyses run. refer reppath() for the
*          keywords.
*
*          the output file can also include keywords. if the output file does
*          not include keywords. the results of all multiple session analyses
*          are output to a single output file.
*
*          ssr corrections are valid only for forward estimation.
*-----------------------------------------------------------------------------*/
extern int execses(gtime_t ts, gtime_t te, double ti, prcopt_t *popt, const solopt_t *sopt, filopt_t *fopt)
{
    FILE *fp;
    /* moved from stack to heap to avoid stack overflow warning */
    rtk_t *rtk = (rtk_t *)malloc(sizeof(rtk_t)); 
    char tracefile[1024],statfile[1024],iposfile[1024],azelfile[1024],satdopfile[1024],path[1024],outfiletm[1024]={0};
    const char *ext;
    int i,j,k,week=0;  

    /* open debug trace */
    opentrace(popt,sopt,fopt);

    /* read obs and nav data */
    if (PMODE_LC_POS!=popt->mode) {
        if (!readobsnav(ts,te,ti,fopt,popt,&obss,&navs,stas)) {
            /* free obs and nav data */
            freeobsnav(&obss,&navs);
            free(rtk);
            return 0;
        } 

        /* read prec ephemeris and sbas data */
        if (!readpreceph(fopt,popt,&navs,&sbss)) {
            /* free prec ephemeris and clock data */
            freepreceph(&navs,&sbss); 
            return 0;       
        }

        /* read satellite and receiver antenna parameters */
        if (*fopt->antp&&!(readpcv(fopt->antp,&pcvss,&pcvsr))) {
            /* free antenna parameters */
            freeant(&pcvss,&pcvsr);
            return 0;
        }

        /* read ionosphere data file */
        if (*fopt->iono&&(ext=strrchr(fopt->iono,'.'))) {
            if (strlen(ext)==4&&(ext[3]=='i'||ext[3]=='I'||strcmp(ext,".INX")==0||strcmp(ext,".inx")==0)) {
                reppath(fopt->iono,path,ts,"","");
                readtec(path,&navs,1);
            }
        }

        /* read erp data */
        if (*fopt->eop) {
            free(navs.erp.data); navs.erp.data=NULL; navs.erp.n=navs.erp.nmax=0;
            reppath(fopt->eop,path,ts,"","");
            readerp(path,&navs.erp);
        }   

        /* read ocean tide loading parameters */
        if (popt->mode>PMODE_SINGLE&&*fopt->blq) {
            readotl(popt,fopt->blq,stas);
        }

        /* open geoid data */
        if (sopt->geoid>0&&*fopt->geoid) {
            if (!opengeoid(sopt->geoid,fopt->geoid)) {
                /* close geoid data */
                closegeoid();            
            }
        }

        /* read dcb parameters from DCB, BIA, BSX files */
        if (*fopt->mgex_dcb) {
            readdcb(popt,fopt->mgex_dcb,&navs,stas);  
        }

        /* unify the correction format of TGD and DCB/OBS */
        tgdarrge(popt,&navs);
    }

    if (PMODE_LC_POS==popt->mode) {
        week=popt->week;
        /* read pos data */
        if (*fopt->pos&&!(readpos(fopt->pos,popt,&poss,week))){
            /* free pos parameters */
            freepos(&poss);
        }
        if (!(*fopt->pos)) showerr("Error : pos file open failed %s",fopt->pos);
    }
 
    if (GINS_OFF!=popt->GI_mode) {
        if (PMODE_LC_POS==popt->mode) week=popt->week;
        else time2gpst(obss.data[0].time,&week);
        /* read imu data */
        if (*fopt->imu&&!(readimu(ts,te,fopt->imu,popt,&imus,week))) {
            /* free imu parameters */
            freeimu(&imus);
        }    
        if (!(*fopt->imu)) showerr("Error : imu file open failed %s",fopt->imu);    
    }

    /* set antenna parameters */
    if (PMODE_LC_POS!=popt->mode&&popt->mode!=PMODE_SINGLE) {
        setpcv(obss.n>0?obss.data[0].time:timeget(),popt,&navs,&pcvss,&pcvsr,
               stas);
    }

    /* rover/reference fixed position */
     if (popt->mode==PMODE_FIXED) {
        if (!antpos(popt,1,&obss,&navs,stas,fopt->stapos)) {
            freeobsnav(&obss,&navs);
            free(rtk);
            return 0;
        }
        if (!antpos(popt,2,&obss,&navs,stas,fopt->stapos)) {
            freeobsnav(&obss,&navs);
            free(rtk);
            return 0;
        }
    }
    else if (PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_STATIC_START) {
        if (!antpos(popt,2,&obss,&navs,stas,fopt->stapos)) {
            freeobsnav(&obss,&navs);
            free(rtk);
            return 0;
        }
    }

    /* write header to output file */
    if (!outhead(popt,sopt,fopt)) {
        freeobsnav(&obss,&navs);
        free(rtk);
        return 0;
    }

    if (PMODE_LC_POS!=popt->mode) {
        /* open solution statistics */
        if (sopt->sstat>0) {
            strcpy(statfile,fopt->sol);
            strcat(statfile,".stat");
            rtkcloseoutfile();
            rtkopenstat(statfile,sopt->sstat);
        }

        /* output ipos file */
        if (sopt->ipos) {
            strcpy(iposfile,fopt->sol);
            strcat(iposfile,".ipos");
            rtkcloseoutfile();
            rtkopenipos(popt,iposfile);
        }

        /* output azel file */
        if (sopt->azel) {
            strcpy(azelfile,fopt->sol);
            strcat(azelfile,".azel");
            rtkcloseoutfile();
            rtkopenazel(popt,azelfile);
        }

        /* output satdop file */
        if (sopt->satdop) {
            strcpy(satdopfile,fopt->sol);
            strcat(satdopfile,".satdop");
            rtkcloseoutfile();
            rtkopensatdop(popt,satdopfile);
        }        
    }

    /* index of obs/imu data */
    iobsu=iobsr=isbs=iimu=aborts=reverse_flag=0;

    /* forward solutions */
    if ((GINS_OFF==popt->GI_mode&&popt->mode==PMODE_SINGLE)||popt->soltype==SOLTYPE_FORWARD) {
        /* solution file (.pos) pointer */
        if (fp=openfile(fopt->sol)) {
            rtkinit(rtk,popt,sopt);
            procpos(fp,popt,sopt,rtk,SOLMODE_SINGLE_DIR);
            rtkfree(rtk); fclose(fp);
        }
    }
    /*  backward solutions */
    else if (popt->soltype==SOLTYPE_BACKWARD) {
        if (fp=openfile(fopt->sol)) {
            rtk->opt.reverse=popt->reverse=1; iobsu=iobsr=obss.n-1; isbs=sbss.n-1;
            if (GINS_OFF!=popt->GI_mode) iimu=imus.n-1;
            rtkinit(rtk,popt,sopt);
            procpos(fp,popt,sopt,rtk,SOLMODE_SINGLE_DIR);
            rtkfree(rtk); fclose(fp);
        }
    }
    /* forward and backward smoothing solutions (FBS) */
    else { 
        solf=(sol_t *)malloc(sizeof(sol_t)*nepoch);
        solb=(sol_t *)malloc(sizeof(sol_t)*nepoch);
        rbf=(double *)malloc(sizeof(double)*nepoch*3);
        rbb=(double *)malloc(sizeof(double)*nepoch*3);

        if (solf&&solb) {
            isolf=isolb=0;
            rtkinit(rtk,popt,sopt);
            procpos(NULL,popt,sopt,rtk,SOLMODE_COMBINED); /* forward */

            rtk->opt.reverse=popt->reverse=1; iobsu=iobsr=obss.n-1; isbs=sbss.n-1;
            if (GINS_OFF!=popt->GI_mode) iimu=imus.n-1;
            /* FBS with no state reset */
            if (popt->soltype!=SOLTYPE_COMBINED_NORESET) {
                /* reset */
                rtkfree(rtk);
                rtkinit(rtk,popt,sopt);
            }
            procpos(NULL,popt,sopt,rtk,SOLMODE_COMBINED); /* backward */
            rtkfree(rtk);

            if (!aborts&&(fp=openfile(fopt->sol))) {
                /* forward and backward smoothing */
                combres(fp,rtk,popt,sopt);
                fclose(fp);
            }
        }
        else showmsg("error : memory allocation");
        free(solf); free(solb); free(rbf); free(rbb);
    }

    /* free rtk, obs/nav , ant and imu data */
    free(rtk);          
    if (obss.n)            freeobsnav(&obss,&navs);
    if (imus.n)            freeimu(&imus);
    if (poss.n)            freepos(&poss);
    if (pcvss.n&&pcvsr.n)  freeant(&pcvss,&pcvsr);
    rtkcloseoutfile();

    return aborts?1:0;
}
