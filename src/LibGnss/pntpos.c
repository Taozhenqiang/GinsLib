/*------------------------------------------------------------------------------
* pntpos.c : standard positioning
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* version : $Revision:$ $Date:$
* history : 2010/07/28 1.0  moved from rtkcmn.c
*                           changed api:
*                               pntpos()
*                           deleted api:
*                               pntvel()
*           2011/01/12 1.1  add option to include unhealthy satellite
*                           reject duplicated observation data
*                           changed api: ionocorr()
*           2011/11/08 1.2  enable snr mask for single-mode (rtklib_2.4.1_p3)
*           2012/12/25 1.3  add variable snr mask
*           2014/05/26 1.4  support galileo and beidou
*           2015/03/19 1.5  fix bug on ionosphere correction for GLO and BDS
*           2018/10/10 1.6  support api change of satexclude()
*           2020/11/30 1.7  support NavIC/IRNSS in pntpos()
*                           no support IONOOPT_LEX option in ioncorr()
*                           improve handling of TGD correction for each system
*                           use E1-E5b for Galileo dual-freq iono-correction
*                           use API sat2freq() to get carrier frequency
*                           add output of velocity estimation error in estvel()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants/macros ----------------------------------------------------------*/
#define NX          (4+5+4)       /* # of estimated parameters (ecef pos, GPS rec clk, ISBs(GLO/Gal/BDS/IRNSS/QZSS), ECEF vel, GPS clk drift) */
#define ERR_ION     5.0         /* ionospheric delay std (m) */
#define ERR_TROP    3.0         /* tropspheric delay std (m) */
#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.2         /* broadcast ionosphere model error factor 0.5 */
#define ERR_CBIAS   0.3         /* code bias error std (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */
#define MAX_GDOP    30           /* max gdop for valid solution  */

/* add by tzq */
#define VAR_POS     SQR(30.0) /* initial variance of receiver pos (m^2) */
#define VAR_VEL     SQR(10.0) /* initial variance of receiver vel ((m/s)^2) */
#define VAR_CLK     SQR(60.0) /* init variance receiver clock (m^2) */
#define VAR_CLKD    SQR(30.0) /* init variance receiver clock drift (m^2/s^2) */

/* number of parameters (pos,vel) */
#define NP(opt)     ((opt)->GI_mode==GINS_TC?GINS_NX:6)
#define NC(opt)     (((opt)->GI_mode==GINS_TC&&(opt)->mode==PMODE_SINGLE)?7:0) /* clock (0-5,GPS,GLONASS,Galileo,BDS,IRNSS,QZSS); clock drift (6:GPS)*/

/* state variable index */
#define IC(s,opt)   (NP(opt)+(s))

/* check if the satellite is a BDS GEO satellite */
extern int isGEOsat(int sat)
{
    int sys,prn;

    sys=satsys(sat,&prn);

    if (sys!=SYS_CMP) return 0;

    if (prn<=5||prn>=59) return 1;
    else return 0;
}

/* determine the number of observations */
extern int obsNum(const rtk_t *rtk, obsd_t *obs, int nobs) 
{
    const prcopt_t *popt=&rtk->opt;
    int i,n;

    /* if GNSS and INS are not synchronized, return 0 observations */
    if (isGINS(popt)&&SYNC_NO==rtk->upte) return 0;

    if (PMODE_LC_POS==popt->mode) n=nobs;
    else { /* exclude satellites */
        for (i=n=0;i<nobs;i++) {
            if ((satsys(obs[i].sat,NULL)&popt->navsys)&&popt->exsats[obs[i].sat-1]!=1) obs[n++]=obs[i];
        }            
    }

    return n;
}

/* save old observation data */
extern void save_old_obs(const obsd_t *obs, obsd_t *obs_old, const int ns, int *nu_old_)
{
    int i,n_old=0,nu_old=0;

    if (ns<=0) return;

    /* save the GNSS observations of the previous epoch */              
    n_old=ns; for (i=0;i<ns;i++) obs_old[i]=obs[i];   

    /* determine the number of satellites of rover in the current epoch and the previous epoch */
    for (i=nu_old=0;i<n_old;i++) if (obs_old[i].rcv==1) nu_old++;

    *nu_old_=nu_old;
}

/* calculate GNSS sampling interval -------------------------------------------
*args  :  rtk_t    *rtk   IO   rtk structure
*         obs_t    *obss  I   observation data
*         pos_t    *poss  I   position data (LC mode)
*return:none
*-----------------------------------------------------------------------------*/
extern int gnss_intervel(rtk_t *rtk, const obs_t *obss, const pos_t *poss)
{
    int i,j,k;
    double t0=0.0,t[2]={0.0},dttol=1e-3;

    if ((!obss&&!poss)||!rtk) {
        return 0;  /* invalid input */
    }

    /* calculate interval from observation file */
    if (PMODE_LC_POS!=rtk->opt.mode) 
    {
        for (i=j=k=0;i<obss->n;i++) {
            if (obss->data[i].rcv!=1) continue;  
            
            for (j=i+1;j<obss->n;j++) {
                if (obss->data[j].rcv!= 1) continue;  /* skip no rover station data */
                
                t0=fabs(timediff(obss->data[i].time,obss->data[j].time));
                
                if (t0>dttol) {  
                    t[k++]=t0;
                    i=j;
                    /* check if intervals are equal */
                    if (k==2) {  
                        if (fabs(t[0]-t[1])<dttol) {  
                            rtk->interval=t[0];
                            return 1;
                        } else {
                            k=0; break;
                        }
                    }
                }
            }
        }        
    }
    else { /* calculate interval from pos file */
        for (i=j=k=0;i<poss->n;i++) {  

            for (j=i+1;j<poss->n;j++) {  

                t0=fabs(timediff(poss->data[i].time,poss->data[j].time));  

                if (t0>dttol) {  
                    t[k++]=t0;
                    i=j;
                    /* check if intervals are equal */
                    if (k==2) {  
                        if (fabs(t[0]-t[1])<dttol) {  
                            rtk->interval=t[0];
                            return 1;
                        } else {
                            k=0; break;
                        }
                    }
                }
            }
        }
    }

    return 1;
}

/* observation pre-check*/
extern int obsScan(const prcopt_t *opt, obsd_t *obs, const int n, int *nu_, int *nr_)
{
	int i,nu,nr,ns=0,sat,sys,fr2[2],prn;
    double threshold=100;
    char id[4];

    /* count rover/base station observations */
    for (nu=0;nu   <n&&obs[nu   ].rcv==1;nu++) ;
    for (nr=0;nu+nr<n&&obs[nu+nr].rcv==2;nr++) ;

    /* init the number of base and rover station observations */
    *nu_=nu;
    *nr_=nr;

	for (i=ns=0;i<n&&i<MAXOBS;i++) {
		sat=obs[i].sat;
        sys=satsys(sat,&prn);
        fr2[0]=sys2freid(sys,0,opt);
        fr2[1]=sys2freid(sys,1,opt);

        /* carrier integrity check for ppp */
        if (opt->mode>=PMODE_PPP_KINEMA) {
            if ((fabs(obs[i].L[fr2[0]])==0.0)&&(fabs(obs[i].L[fr2[1]])==0.0)) {
                (*nu_)--;
                continue;
            }
        }

        /* pseudorange outlier detection */
        if ((obs[i].P[fr2[0]]!=0.0&&fabs(obs[i].P[fr2[0]])<19e6)||(obs[i].P[fr2[1]]!=0.0&&fabs(obs[i].P[fr2[1]])<19e6)) {
            if (i<nu) (*nu_)--; else (*nr_)--;
            satno2id(sat,id);
            trace(6,"obsScan: abnormal pseudorange observations, less than 19000 km, sat=%s\n",id);
            continue;            
        }
        if (obs[i].P[fr2[0]]!=0.0&&obs[i].P[fr2[1]]!=0.0&&fabs(obs[i].P[fr2[0]]-obs[i].P[fr2[1]])>=threshold) {
            if (i<nu) (*nu_)--; else (*nr_)--;
            satno2id(sat,id);
            trace(6,"obsScan: dual-frequency pseudorange difference exceeded the limit, sat=%s\n",id);
            continue;
        }

        obs[ns]=obs[i];
        ns++;
	}

    if (!ns) return 0;

    return ns;
}

/* the sign of doppler observations is determined based on pseudorange variation between adjacent epochs */
extern int dopple_sgn(rtk_t *rtk, const obsd_t *obs, const obsd_t *obs_old, int nu, int n_old)
{   
    int i,j,fr,sys,nu_old;
    double dr;

    /* if no previous epoch data is available, exit current epoch processing */
    if (n_old<=0) return 0;
    
    /* determine the number of satellites of rover in the current epoch and the previous epoch */
    for (i=nu_old=0;i<n_old;i++) if (obs_old[i].rcv==1) nu_old++;

    for (i=0;i<nu;i++) {
        sys=satsys(obs[i].sat,NULL);
        fr=sys2freid(sys,0,&rtk->opt);
        for (j=0;j<nu_old;j++) {
            if (obs[i].sat==obs_old[j].sat&&obs[i].P[fr]&&obs_old[j].P[fr]&&obs[i].D[fr]) break;
        }
        if (j>=nu_old) continue;
        else break;
    }

    if (i>=nu) return 0;

    /* pseudorange variation between adjacent epochs */
    dr=obs[i].P[fr]-obs_old[j].P[fr];

    /* when the satellite is close to the receiver, the doppler sign is positive; otherwise, it is negative */
    if (SGN(dr)==-SGN(obs[i].D[fr])) {
        rtk->dopsgn=-1.0;
    }
    else {
        rtk->dopsgn=1.0;
    }

    /* in both forward and backward processing modes, the sign of doppler observations remains unchanged ? */
    if (SOLTYPE_BACKWARD==rtk->opt.reverse) {
        rtk->dopsgn=-rtk->dopsgn;
    }

    return 1;
}

/* carrier-phase bias correction by ssr --------------------------------------*/
extern void corr_phase_bias_ssr(obsd_t *obs, int n, const nav_t *nav)
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

/* mutipath correct-------------------------------------------------------------
* BeiDou satellite-induced code pseudorange variations correct
* args  :rtk_t *rtk       IO  rtk control/result struct
           obsd_t *obs      IO  observation data
           int    n         I   number of observation data
           nav_t  *nav      I   navigation messages
* note   :
* -----------------------------------------------------------------------------*/
extern void BDmulCorr(rtk_t *rtk, obsd_t *obs, int n) 
{
    int i,j,sat,prn,b,*f=(int *)&rtk->opt.fre[4],ix[3]={-1,-1,-1};
    double dp[3],elev,a;

    const static double IGSOCOEF[3][10]={
        /* m */
        {-0.55,-0.40,-0.34,-0.23,-0.15,-0.04,0.09,0.19,0.27,0.35},// B1
        {-0.71,-0.36,-0.33,-0.19,-0.14,-0.03,0.08,0.17,0.24,0.33},// B2
        {-0.27,-0.23,-0.21,-0.15,-0.11,-0.04,0.05,0.14,0.19,0.32},// B3
    };
    const static double MEOCOEF[3][10]={
        /* m */
        {-0.47,-0.38,-0.32,-0.23,-0.11,0.06,0.34,0.69,0.97,1.05},// B1
        {-0.40,-0.31,-0.26,-0.18,-0.06,0.09,0.28,0.48,0.64,0.69},// B2
        {-0.22,-0.15,-0.13,-0.10,-0.04,0.05,0.14,0.27,0.36,0.47},// B3
    };

    for (i=0;i<n&&i<MAXOBS;i++) {
        sat=obs[i].sat;

        if (satsys(sat,&prn)!=SYS_CMP) continue;
        if (prn<=5) continue;

        elev=rtk->ssat[sat-1].azel[1]*R2D;

        if (elev<=0.0) continue;

        for (j=0;j<3;j++) dp[j]=0.0;

        a=elev*0.1;
        b=(int) a;

        if (prn>=6&&prn<11) { // IGSO(C06,C07,C08,C09,C10)
            if (b<0) {
                for (j=0;j<3;j++)
                    dp[j]=IGSOCOEF[j][0];
            } else if (b >=9) {
                for (j=0;j<3;j++)
                    dp[j]=IGSOCOEF[j][9];
            } else {
                for (j=0;j<3;j++)
                    dp[j]=IGSOCOEF[j][b]*(1.0-a+b)+IGSOCOEF[j][b+1]*(a-b);
            }
        } else if (prn>=11&&prn<=14) { // MEO(C11,C12,C13,C14)
            if (b<0) {
                for (j=0;j<3;j++)
                    dp[j]=MEOCOEF[j][0];
            } else if (b >=9) {
                for (j=0;j<3;j++)
                    dp[j]=MEOCOEF[j][9];
            } else {
                for (j=0;j<3;j++)
                    dp[j]=MEOCOEF[j][b]*(1.0-a+b)+MEOCOEF[j][b+1]*(a-b);
            }
        } else
            continue;

        /* find idx of B1I,B2I,B3I */
        for (j=0;j<MAXFREQ;j++) {
            if (0==f[j]) ix[0]=j; /* B1I */
            if (1==f[j]) ix[1]=j; /* B2I */
            if (2==f[j]) ix[2]=j; /* B3I */
        }
        for (j=0;j<3;j++) {
            if (obs[i].P[ix[j]]>0.0&&ix[j]>0) obs[i].P[ix[j]]+=dp[j]; 
        }
    }
}

/* observation preprocessing */
extern int obsPreprocess(rtk_t *rtk, obsd_t *obs, obsd_t *obs_old, const nav_t *nav, const int n, const int n_old, int *nu_, int *nr_)
{
    const prcopt_t *popt=&rtk->opt;
    int ns;

    /* check observation data */
    if (n<=0||(isGINS(popt)&&SYNC_NO==rtk->upte)) return 0;

    /* for LC_POS mode, no observation preprocessing is required */
    if (PMODE_LC_POS==popt->mode) {
        ns=*nu_=n;
        return ns;
     }

    /* observation pre-check*/
    ns=obsScan(popt,obs,n,nu_,nr_);

    /* the sign of Doppler observations is determined based on pseudorange variation between adjacent epochs */
    if ((isGNSS(popt)||SYNC_YES==rtk->upte)&&PMODE_LC_POS!=popt->mode&&!rtk->dopsgn) dopple_sgn(rtk,obs,obs_old,*nu_,n_old);
    
    /* carrier-phase bias correction */
    if (PMODE_DGPS<popt->mode&&!strstr(popt->pppopt,"-ENA_FCB")) {
        corr_phase_bias_ssr(obs,ns,nav);
    }

    /* multipath correction for BDS2 */
    if (popt->navsys&SYS_CMP) {
        BDmulCorr(rtk,obs,ns); 
    }

    return ns;
}

/* pseudorange measurement error variance ------------------------------------*/
extern double varerr_spp(const prcopt_t *opt, const ssat_t *ssat, const obsd_t *obs, double el, int sys)
{
    int fr2[2]={0};
    double a,b,c,fact=1.0,varr,snr_rover,snr_thres=5.0,snr_slid=0.0,dsnr=0.0;

    fr2[0]=sys2freid(sys,0,opt);
    fr2[1]=sys2freid(sys,1,opt);
    if (ssat) {
        if (IONOOPT_IFLC==opt->ionoopt) {
            snr_slid=(ssat->snr_slidr[fr2[0]]+ssat->snr_slidr[fr2[1]])/2.0;    
            snr_rover=SNR_UNIT*(ssat->snr_rover[fr2[0]]+ssat->snr_rover[fr2[1]])/2.0;
        } 
        else {
            snr_slid=ssat->snr_slidr[fr2[0]];
            snr_rover=SNR_UNIT*ssat->snr_rover[fr2[0]];
        }   
        
        /* if window size is less than SNR_WINDOW, set snr_slid unabled */
        if (ssat->window_size[fr2[0]]<=SNR_WINDOW) snr_slid=0.0;
    }

    /* pseudoranges error (m) */
    fact=3;

    switch (sys) {
        case SYS_GPS: fact*=EFACT_GPS; break;
        case SYS_GLO: fact*=EFACT_GLO; break;
        case SYS_SBS: fact*=EFACT_SBS; break;
        case SYS_CMP: fact*=EFACT_CMP; break;
        case SYS_QZS: fact*=EFACT_QZS; break;
        case SYS_IRN: fact*=EFACT_IRN; break;
        default:      fact*=EFACT_GPS; break;
    }
    if (el<MIN_EL) el=MIN_EL;

    /* experience error factor */
    a=b=fact;  
    /* elevation term */
    /* var = (a^2 + b^2/sin(el) + c^2*(10^(0.1*(snr_max-snr_rover))) + (d*rcv_std)^2) */
    varr=a*a+b*b/sin(el);

    /* snr term */
    dsnr=(ssat)?(snr_slid-snr_rover):0.0;
    if (opt->err[6]>0.0&&dsnr>=snr_thres) {  /* if snr term not zero */
        c=3.0;
        varr+=c*c*pow(10,0.1*dsnr); 

        trace(7,"var_el=%10.4f, var_snr=%5.1f\n",a*a+b*b/sin(el),c*c*pow(10,0.1*dsnr));
    }

    /* iono-free scale factor */
    if (opt->ionoopt==IONOOPT_IFLC) varr*=SQR(3.0); 

    return varr;
}
/* get group delay parameter (m) ---------------------------------------------*/
extern double gettgd(int sat, const nav_t *nav, int type)
{
    int i,sys=satsys(sat,NULL);
    
    if (sys==SYS_GLO) {
        for (i=0;i<nav->ng;i++) {
            if (nav->geph[i].sat==sat) break;
        }
        return (i>=nav->ng)?0.0:-nav->geph[i].dtaun*CLIGHT;
    }
    else {
        for (i=0;i<nav->n;i++) {
            if (nav->eph[i].sat==sat) break;
        }
        return (i>=nav->n)?0.0:nav->eph[i].tgd[type]*CLIGHT;
    }
}
/* test SNR mask -------------------------------------------------------------*/
static int snrmask(const int iter, const obsd_t *obs, const double *azel, const prcopt_t *opt)
{
    int sys,fr;
    char id[4];
    fr=sys2freid(sys,0,opt);

    if (testsnr(0,0,azel[1],obs->SNR[fr]*SNR_UNIT,&opt->snrmask)) {
        if (iter<=1&&PMODE_SINGLE==opt->mode) {
            satno2id(obs->sat,id);
            trace(6,"SNR check failed(spp): %4s, el=%4.1f, SNR=%5.1f\n",id,azel[1]*R2D,obs->SNR[fr]*SNR_UNIT);            
        }
        return 0;
    }
    if (opt->ionoopt==IONOOPT_IFLC) {
        sys=satsys(obs->sat,NULL);
        fr=sys2freid(sys,1,opt);
        if (testsnr(0,1,azel[1],obs->SNR[fr]*SNR_UNIT,&opt->snrmask)) return 0;
    }
    return 1;
}
/* iono-free or "pseudo iono-free" pseudorange with code bias correction -----*/
static double prange(const int fr, const obsd_t *obs, const nav_t *nav, const prcopt_t *opt,
                     double *var, double *dcb)
{
    double P1,P2,gamma,b1=0.0,freq1=0.0,freq2=0.0;
    int sat,sys,id,f2,bias_ix[2],flag=0,fr2[2];

    sat=obs->sat;
    sys=satsys(sat,&id);
    if (IONOOPT_IFLC==opt->ionoopt) {
        fr2[0]=sys2freid(sys,0,opt);
        fr2[1]=sys2freid(sys,1,opt);        
    }
    else {
        fr2[0]=sys2freid(sys,fr,opt);
        fr2[1]=sys2freid(sys,1,opt);        
    }

    P1=obs->P[fr2[0]];
    P2=obs->P[fr2[1]];
    *var=0.0;*dcb=0.0;
    
    if (P1==0.0||(opt->ionoopt==IONOOPT_IFLC&&P2==0.0)) return 0.0;
    bias_ix[0]=code2bias_ix(sys,obs->code[fr2[0]]);  /* L1 code bias */
    bias_ix[1]=code2bias_ix(sys,obs->code[fr2[1]]);

    /* obias_flag 1:DCB product, 2:OSB product */
    /* BDS broadcast ephemeris DCB corrections */
    if (opt->sateph==EPHOPT_BRDC&&sys==SYS_CMP) { 
        if (bias_ix[0]>=0) {
            P1-=nav->bds_tgd[id-1][bias_ix[0]];
            *dcb=nav->bds_tgd[id-1][bias_ix[0]];
        }
        if (bias_ix[1]>=0) P2-=nav->bds_tgd[id-1][bias_ix[1]];
    }
    else if (nav->obias_flag>0) {
        if (bias_ix[0]>=0) {
            P1-=nav->obias[sat-1][bias_ix[0]];
            *dcb=nav->obias[sat-1][bias_ix[0]];
        }
        if (bias_ix[1]>=0) P2-=nav->obias[sat-1][bias_ix[1]];
    }

    /* P1-C1,P2-C2 DCB correction */
    if (sys==SYS_GPS||sys==SYS_GLO) {
        if (obs->code[fr2[0]]==CODE_L1C) P1+=nav->cbias[sat-1][0]; /* C1->P1 */
        if (obs->code[fr2[1]]==CODE_L2C) P2+=nav->cbias[sat-1][1]; /* C2->P2 */
    }

    if (opt->ionoopt==IONOOPT_IFLC) { /* dual-frequency */     
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2 or L1-L5 */
            freq1=code2freq(SYS_GPS,obs->code[fr2[0]],0);
            freq2=code2freq(SYS_GPS,obs->code[fr2[1]],0);
            gamma=SQR(freq1/freq2);
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_GLO) { /* G1-G2 or G1-G3 */
            gamma=f2==1?SQR(FREQ1_GLO/FREQ2_GLO):SQR(FREQ1_GLO/FREQ3_GLO);
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_GAL) { /* E1-E5b, E1-E5a */
            freq1=code2freq(SYS_GPS,obs->code[fr2[0]],0);
            freq2=code2freq(SYS_GPS,obs->code[fr2[1]],0);
            gamma=SQR(freq1/freq2);
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_CMP) { /* B1-B2 */
            freq1=code2freq(SYS_CMP,obs->code[fr2[0]],0);
            freq2=code2freq(SYS_CMP,obs->code[fr2[1]],0);
            gamma=SQR(freq1/freq2);
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_IRN) { /* L5-S */
            gamma=SQR(FREQL5/FREQs);
            return (P2-gamma*P1)/(1.0-gamma);
        }
    }
    else { /* single-freq (L1/E1/B1) */
        *var=SQR(ERR_CBIAS);
        
        if (sys==SYS_GPS||sys==SYS_QZS||SYS_GAL) { /* L1 */
            return P1;
        }
        else if (sys==SYS_GLO) { /* G1 */
            gamma=SQR(FREQ1_GLO/FREQ2_GLO);
            b1=gettgd(sat,nav,0); /* -dtaun (m) */
            return P1-b1/(gamma-1.0);
        }
        else if (sys==SYS_CMP) {
            return P1;
        }
        else if (sys==SYS_IRN) { /* L5 */
            gamma=SQR(FREQs/FREQL5);
            b1=gettgd(sat,nav,0); /* TGD (m) */
            return P1-gamma*b1;
        }
    }
}
/* ionospheric correction for spp ------------------------------------------------------
* compute ionospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          int    sat       I   satellite number
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    ionoopt   I   ionospheric correction option (IONOOPT_???)
*          double *ion      O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric delay (L1) variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int ionocorr(const obsd_t *obs, gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, const prcopt_t *opt, double *ion, double *var)
{
    int i,ionoopt=opt->ionoopt,sys,fr2[2],flag=0;
    double freq1=0.0,freq2=0.0;

    trace(9,"ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
          time_str(time,3),ionoopt,sat,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* SBAS ionosphere model */
    if (ionoopt==IONOOPT_SBAS) {
        if (sbsioncorr(time,nav,pos,azel,ion,var)) return 1;
    }
    /* IONEX TEC model */
    if (ionoopt==IONOOPT_TEC) {
        if (iontec(time,nav,pos,azel,1,ion,var)) return 1;
    }
    /* QZSS broadcast ionosphere model */
    if (ionoopt==IONOOPT_QZS&&norm(nav->ion_qzs,8)>0.0) {
        *ion=ionmodel(time,nav->ion_qzs,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* Double-frequency ionosphere model */
    if (ionoopt==IONOOPT_DF) {
        sys=satsys(sat,NULL);
        fr2[0]=sys2freid(sys,0,opt);
        fr2[1]=sys2freid(sys,1,opt);
        freq1=sat2freq(sat,obs->code[fr2[0]],nav);
        freq2=sat2freq(sat,obs->code[fr2[1]],nav);
        if (obs->P[fr2[0]]==0.0||obs->P[fr2[1]]==0.0||freq1==0.0||freq2==0.0) {
            flag=0;
        }
        else {
            /* slant ionospheric delay based GPS L1 frequency */
            *ion=(obs->P[fr2[0]]-obs->P[fr2[1]])/(SQR(FREQL1/freq1)-SQR(FREQL1/freq2));
            *var=SQR(*ion*ERR_BRDCI);
            return 1;            
        }
    }    
    /* GPS broadcast ionosphere model */
    if (ionoopt==IONOOPT_BRDC||ionoopt==IONOOPT_EST||!flag) {
        *ion=ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }

    *ion=0.0;
    *var=ionoopt==IONOOPT_OFF?SQR(ERR_ION):0.0;
    return 1;
}
/* tropospheric correction for spp -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    tropopt   I   tropospheric correction option (TROPOPT_???)
*          double *trp      O   tropospheric delay (m)
*          double *var      O   tropospheric delay variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var)
{
    trace(9,"tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f\n",
          time_str(time,3),tropopt,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* Saastamoinen model */
    if (tropopt==TROPOPT_SAAS||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=tropmodel(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        return 1;
    }
    /* SBAS (MOPS) troposphere model */
    if (tropopt==TROPOPT_SBAS) {
        *trp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    /* no correction */
    *trp=0.0;
    *var=tropopt==TROPOPT_OFF?SQR(ERR_TROP):0.0;
    return 1;
}
/* calculating expectation (exclude 0 elements) -----------------------------------------------------
* args   : double *data     I   data(nx1)
*          int    *size     I   size of data
* return : mean
*-----------------------------------------------------------------------------*/
extern double calexp(const double *data, int size) 
{
    int i=0,size_=size;
    double sum=0.0;

    for (i=0;i<size;i++) {
        if (fabs(data[i])<1E-4) {
            size_--;
            continue;
        }
        sum+=data[i];
    }
    return sum/size_;
}
/* calculate standard deviation (exclude 0 elements) -----------------------------------------------------
* args   : double *data     I   data(nx1)
*          int    *size     I   size of data
* return : std
*-----------------------------------------------------------------------------*/
extern double calstd(const double *data, int size) 
{
    int i=0,size_=size;
    double mean=calexp(data,size);
    double sum=0.0;

    for (i=0;i<size;i++) {
        if (fabs(data[i])<1E-4) {
            size_--;
            continue;
        }
        sum+=(data[i]-mean)*(data[i]-mean);
    }
    return sqrt(sum/size_);
}
/* pseudorange residuals -----------------------------------------------------*/
static int rescode(int iter, const obsd_t *obs, int n, int nx_code, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   ssat_t *ssat, double *v, double *H, double *var, int nx, int *mask, int *clock_idx,
                   double *azel, int *vsat, double *resp, int *ns, int *sati, int *vi)
{
    gtime_t time;
    double r,freq,dion=0.0,dtrp=0.0,vmeas,vion=0.0,vtrp=0.0,vobs=0.0,rr[3],pos[3],e[3],P,dtr=0.0,dcb=0.0;
    int i,j,k,nf=opt->mfspp?(IONOOPT_IFLC==opt->ionoopt?1:opt->nf):1,nv=0,sat,sys,fr=0,idx=0;
    char id[4];

    /* update receiver position and GPS receiver clock bias based on LS estimation */
    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[3];
    
    ecef2pos(rr,pos);
    trace(8,"rescode: rr=%.3f %.3f %.3f\n",rr[0], rr[1], rr[2]);

    /* frequency loop */
    for (k=0;k<nf;k++) {
        for (i=0;i<NX-3;i++) mask[i]=0; /* reset mask flag */
        /* observation loop */
        for (i=ns[k]=0;i<n&&i<MAXOBS;i++) {
            time=obs[i].time; sat=obs[i].sat;
            fr=sys2freid(sys,k,opt);
            idx=i+n*k;

            /* init azimuth/elevation angle and vsat flag */
            if (!iter&&!k) azel[i*2]=azel[1+i*2]=0.0;
            if (!iter)     vsat[idx]=0;
            resp[idx]=0.0;

            /* reset spp valid satellite flags (only for ipos result output) */
            if (ssat) ssat[sat-1].vs=0;
            
            /* exclude satellites with large residuals */
            if (iter&&!vsat[idx]) continue;
            if (!(sys=satsys(sat,NULL))) continue;
            
            /* reject duplicated observation data */
            if (i<n-1&&i<MAXOBS-1&&sat==obs[i+1].sat) {
                satno2id(sat,id); trace(7,"duplicated obs data sat=%s\n",id);
                i++; continue;
            }
            
            /* excluded satellite */
            if (satexclude(sat,vare[i],svh[i],opt)) continue;
            
            /* geometric distance and elevation mask*/
            if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
            
            if (iter>0) {      
                /* test elevation mask */
                if (satazel(pos,e,azel+i*2)<opt->elmin) continue;

                /* test SNR mask */
                if (!snrmask(iter,obs+i,azel+i*2,opt)) continue;

                /* ionospheric correction */
                if (!ionocorr(obs+i,time,nav,sat,pos,azel+i*2,opt,&dion,&vion)) {
                    continue;
                }
                if ((freq=sat2freq(sat,obs[i].code[fr],nav))==0.0) continue;
                /* convert from FREQL1 to freq */
                dion*=SQR(FREQL1/freq);
                vion*=SQR(SQR(FREQL1/freq));
            
                /* tropospheric correction */
                if (!tropcorr(time,nav,pos,azel+i*2,opt->tropopt,&dtrp,&vtrp)) {
                    continue;
                }
            }
            /* pseudorange with code bias correction */
            if ((P=prange(k,obs+i,nav,opt,&vmeas,&dcb))==0.0) continue;
            
            /* pseudorange residual */
            v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
            if (ssat) trace(8,"sat=%d: v=%.3f P=%.3f r=%.3f dts=%.6f dion=%.3f dtrp=%.3f\n",ssat[sat-1].id,v[nv],P,r,dts[i*2],dion,dtrp);
            
            /* design matrix */
            for (j=0;j<nx;j++) {
                H[j+nv*nx]=j<3?-e[j]:(j==3?1.0:0.0);
            }                

            /* time system offset and receiver bias correction */
            if      (sys==SYS_GLO) {v[nv]-=x[clock_idx[1]];  H[clock_idx[1]+nv*nx]=1.0; mask[1]=1;}
            else if (sys==SYS_GAL) {v[nv]-=x[clock_idx[2]];  H[clock_idx[2]+nv*nx]=1.0; mask[2]=1;}
            else if (sys==SYS_CMP) {v[nv]-=x[clock_idx[3]];  H[clock_idx[3]+nv*nx]=1.0; mask[3]=1;}
            else if (sys==SYS_IRN) {v[nv]-=x[clock_idx[4]];  H[clock_idx[4]+nv*nx]=1.0; mask[4]=1;}
            else if (sys==SYS_QZS) {v[nv]-=x[clock_idx[5]];  H[clock_idx[5]+nv*nx]=1.0; mask[5]=1;}
            else mask[0]=1;
            
            vsat[idx]=1; resp[idx]=v[nv]; ns[k]++; sati[nv]=sat; vi[nv]=idx;

            /* update solution status */
            if (ssat) update_ssat(ssat+(sat-1),opt,1,sat,fr,rs+i*6,rr,r,azel+i*2,resp[idx],dtr,dts[i*2],dtrp,dion,0.0,NULL,NULL,dcb);

            /* variance of pseudorange error */
            var[nv]=vare[i]+vmeas+vion+vtrp;
            if (ssat) {
                vobs=varerr_spp(opt,ssat+(sat-1),&obs[i],azel[1+i*2],sys); 
            }           
            else {
                vobs=varerr_spp(opt,NULL,&obs[i],azel[1+i*2],sys);
            }  
            var[nv++]+=vobs;
#if 0                 
            if (iter==3) trace(7,"sat=%3s el=%3.1f var=%5.1f vare=%5.1f vmeas=%5.1f vion=%5.1f vtrp=%5.1f vobs=%5.1f\n",ssat[sat-1].id,azel[1+i*2]*R2D,var[nv-1],vare[i],vmeas,vion,vtrp,vobs);
#endif
        }  
        /* constraint to avoid rank-deficient */
        /* multi system constraint */
        for (i=0;i<nx_code-4;i++) {
            if ((clock_idx[1]==i+4)&&mask[1]) continue; /* GLONASS */
            else if ((clock_idx[2]==i+4)&&mask[2]) continue; /* GALILEO */
            else if ((clock_idx[3]==i+4)&&mask[3]) continue; /* BDS */
            else if ((clock_idx[4]==i+4)&&mask[4]) continue; /* IRN */
            else if ((clock_idx[5]==i+4)&&mask[5]) continue; /* QZS */

            v[nv]=0.0;
            for (j=0;j<nx;j++) {
                if ((clock_idx[1]==i+4)&&!mask[1]) H[j+nv*nx]=j==clock_idx[1]?1.0:0.0;   
                else if ((clock_idx[2]==i+4)&&!mask[2]) H[j+nv*nx]=j==clock_idx[2]?1.0:0.0;   
                else if ((clock_idx[3]==i+4)&&!mask[3]) H[j+nv*nx]=j==clock_idx[3]?1.0:0.0;   
                else if ((clock_idx[4]==i+4)&&!mask[4]) H[j+nv*nx]=j==clock_idx[4]?1.0:0.0;   
                else if ((clock_idx[5]==i+4)&&!mask[5]) H[j+nv*nx]=j==clock_idx[5]?1.0:0.0;   
            }
            var[nv++]=0.01;
        } 
        /* single system constraint (exclude GPS) */ 
        if (opt->navsys==SYS_GLO) {
            for (j=0;j<nx;j++)  H[j+nv*nx]=j==clock_idx[1]?1.0:0.0;
            v[nv]=0.0; var[nv++]=0.01;
        }    
        else if (opt->navsys==SYS_GAL) {
            for (j=0;j<nx;j++)  H[j+nv*nx]=j==clock_idx[2]?1.0:0.0;
            v[nv]=0.0; var[nv++]=0.01;
        }
        else if (opt->navsys==SYS_CMP) {
            for (j=0;j<nx;j++)  H[j+nv*nx]=j==clock_idx[3]?1.0:0.0;
            v[nv]=0.0; var[nv++]=0.01;
        }
        else if (opt->navsys==SYS_IRN) {
            for (j=0;j<nx;j++)  H[j+nv*nx]=j==clock_idx[4]?1.0:0.0;
            v[nv]=0.0; var[nv++]=0.01;
        }
        else if (opt->navsys==SYS_QZS) {
            for (j=0;j<nx;j++)  H[j+nv*nx]=j==clock_idx[5]?1.0:0.0;
            v[nv]=0.0; var[nv++]=0.01;
        }
    } 

    return nv;
}

/* pseudorange residuals -----------------------------------------------------*/
static int rescode_filter(rtk_t *rtk, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   ssat_t *ssat, double *v, double *H, double *var,
                   double *azel, int *vsat, double *resp, int *ns, int *sati)
{
    gtime_t time;
    ins_t *ins=&rtk->ins;
    double r,freq,dion=0.0,dtrp=0.0,vmeas,vion=0.0,vtrp=0.0,rr[3],pos[3],e[3],P,dtr=0.0,dcb;
    double Hpa[3],Hpp[3],Cne[9],Cen[9],Cnb[9],lever_n[3]; /* spp/ins */
    int i,j,k,nx=(GINS_TC==opt->GI_mode)?rtk->nx:NX,nv=0,sat,sys,fr=0,ic,nf=opt->mfspp?(IONOOPT_IFLC==opt->ionoopt?1:opt->nf):1,idx;
    int spp_tc_flag=(GINS_TC==opt->GI_mode);
    char id[4];

    /* receiver position and clock error (m) */
    for (i=0;i<3;i++) rr[i]=x[i];
    ic=IC(0,opt); dtr=rtk->x[ic];

    if (spp_tc_flag) {
        xyz2enu(ins->pos,Cne); DCMT(Cne,Cen);
        Mat3mulv(1.0,ins->Cnb,ins->lever,lever_n);        
    }
    
    ecef2pos(rr,pos);
    trace(8,"rescode: rr=%.3f %.3f %.3f\n",rr[0], rr[1], rr[2]);
    
    for (k=0;k<nf;k++) {
        for (i=*ns=0;i<n&&i<MAXOBS;i++) { 
            time=obs[i].time; sat=obs[i].sat;
            fr=sys2freid(sys,k,opt);
            idx=i+n*k;

            if (!k) azel[i*2]=azel[1+i*2]=0.0;
            resp[idx]=0.0;

            /* reset spp valid satellite flags */
            if (ssat) ssat[sat-1].vs=0;

            /* exclude satellites with large residuals */
            if (!vsat[idx]) continue;
            if (!(sys=satsys(sat,NULL))) continue;
            
            /* reject duplicated observation data */
            if (i<n-1&&i<MAXOBS-1&&sat==obs[i+1].sat) {
                satno2id(sat,id); trace(7,"duplicated obs data sat=%s\n",id);
                i++; continue;
            }

            /* excluded satellite */
            if (satexclude(sat,vare[i],svh[i],opt)) continue;
            
            /* geometric distance and elevation mask*/
            if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
            if (satazel(pos,e,azel+i*2)<opt->elmin) continue;

            /* test SNR mask */
            if (!snrmask(0,obs+i,azel+i*2,opt)) continue;

            /* ionospheric correction */
            if (!ionocorr(obs+i,time,nav,sat,pos,azel+i*2,opt,&dion,&vion)) {
                continue;
            }
            if ((freq=sat2freq(sat,obs[i].code[fr],nav))==0.0) continue;
            /* convert from L1 to current frequency */
            dion*=SQR(FREQL1/freq); vion*=SQR(SQR(FREQL1/freq));
        
            /* tropospheric correction */
            if (!tropcorr(time,nav,pos,azel+i*2,opt->tropopt,&dtrp,&vtrp)) {
                continue;
            }

            /* pseudorange with code bias correction */
            if ((P=prange(k,obs+i,nav,opt,&vmeas,&dcb))==0.0) continue;
            
            /* pseudorange residual (v=z-h(x)) */
            v[nv]=(P-(r-CLIGHT*dts[i*2]+dion+dtrp))-dtr;
            if (ssat) trace(8,"sat=%d: v=%.3f P=%.3f r=%.3f dts=%.6f dion=%.3f dtrp=%.3f\n",ssat[i].id,v[nv],P,r,dts[i*2],dion,dtrp);

            for (j=0;j<nx;j++) {
                H[j+nv*nx]=(j==IC(0,opt)?1.0:0.0);
            }    
            
            /* design matrix H */
            if (spp_tc_flag) {
                vmulMat3(1.0,e,Cen,Hpp);
                vmvskew(1.0,Hpp,lever_n,Hpa);

                for (j=0;j<3;j++)   H[j+nv*nx]    =Hpa[j];
                for (j=0;j+6<9;j++) H[(j+6)+nv*nx]=Hpp[j];                
            }
            else {
                for (j=0;j<3;j++)   H[j+nv*nx]    =-e[j];
            }

            /* time system offset and receiver bias correction */
            if      (sys==SYS_GLO) {ic=IC(1,opt); v[nv]-=rtk->x[ic]; H[ic+nv*nx]=1.0; }
            else if (sys==SYS_GAL) {ic=IC(2,opt); v[nv]-=rtk->x[ic]; H[ic+nv*nx]=1.0; }
            else if (sys==SYS_CMP) {ic=IC(3,opt); v[nv]-=rtk->x[ic]; H[ic+nv*nx]=1.0; }
            else if (sys==SYS_IRN) {ic=IC(4,opt); v[nv]-=rtk->x[ic]; H[ic+nv*nx]=1.0; }
            else if (sys==SYS_QZS) {ic=IC(5,opt); v[nv]-=rtk->x[ic]; H[ic+nv*nx]=1.0; }
            
            vsat[idx]=1; resp[idx]=v[nv]; (*ns)++; sati[nv]=sat;

            /* update solution status */
            if (ssat) update_ssat(ssat+(sat-1),opt,1,sat,fr,rs+i*6,rr,r,azel+i*2,resp[idx],dtr,dts[i*2],dtrp,dion,0.0,NULL,NULL,dcb);
            
            /* variance of pseudorange error */
            var[nv]=vare[i]+vmeas+vion+vtrp;
            if (ssat) {
                var[nv++]+=varerr_spp(opt,&ssat[i],&obs[i],azel[1+i*2],sys);
            }    
            else {
                var[nv++]+=varerr_spp(opt,NULL,&obs[i],azel[1+i*2],sys);
            }         
            trace(8,"sat=%2d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n",obs[i].sat,
                azel[i*2]*R2D,azel[1+i*2]*R2D,resp[idx],sqrt(var[nv-1]));
        }        
    }

    return nv;
}
/* outlier rejection for spp ---------------------------------------------------------*/
extern int outrej_spp(int nv, int *nv_code, int nx, int nx_code, double thres, double *v, double *H, double *var,
                      const ssat_t *ssat, const int *sati, const int *vi, int *vsat, int it, int *clock_idx) 
{
    double mean,std,dv;
    int j,k,m=0,mask[NX-3]={0},nv_code_=0;
    double *v_,*H_,*var_,sum=0.0;

    if (nv_code) nv_code_=*nv_code;

    v_=mat(nv,1); H_=mat(nv,nx); var_=mat(nv,1);

    mean=calexp(v,nv); std=calstd(v,nv);

    for (j=0;j<nv;j++) {  
        /* remove auxiliary quantities that prevent least squares rank deficiency */
        if (fabs(v[j])<1E-4) dv=0.0;
        /* standardized residuals */
        else dv=fabs(v[j]-mean)/std;     

        /* threshold for outlier detection (only for code) */     
        if ((nv_code_>0&&j>=nv_code_)||dv<thres) {
            v_[m]=v[j];
            for (k=0;k<nx;k++) {
                H_[k+m*nx]=H[k+j*nx];
            }  
            var_[m++]=var[j];
        }
        else {
            if (vsat&&vi) vsat[vi[j]]=0;
            if (nv_code) (*nv_code)--;
            if (ssat) trace(6,"iteration(%2d), outlier rejected(spp) sat=%s, res=%13.4f, dv=%10.4f, thres=%5.2f, el=%4.1f\n",it+1,ssat[sati[j]-1].id,
            fabs(v[j]),fabs(v[j]-mean)/std,thres,ssat[sati[j]-1].azel[1]*R2D);
            continue;
        }               
    }
    /* update the number of valid observations */
    nv=m;
    for (j=0;j<nv;j++) {
        v[j]=v_[j];
        for (k=0;k<nx;k++) {
            H[k+j*nx]=H_[k+j*nx];
        }  
        var[j]=var_[j];           
    }
    if (clock_idx) {
        /* satellite system availability check. avoid measurement matrix H rank deficiency caused by satellite removal. */
        for (j=4;j<nx_code;j++) {
            sum=0.0;
            for (k=0;k<nv;k++) {
                sum+=fabs(H[j+k*nx]);
            }
            if (sum>0.0) {
                if (clock_idx[1]==j) mask[1]=1; /* GLONASS */
                else if (clock_idx[2]==j) mask[2]=1; /* Galileo */
                else if (clock_idx[3]==j) mask[3]=1; /* BDS */
                else if (clock_idx[4]==j) mask[4]=1; /* IRN */
                else if (clock_idx[5]==j) mask[5]=1; /* QZS */
            }
        }
        /* constraint to avoid rank deficiency of H matrix */
        for (k=0;k<nx_code-4;k++) {
            if ((clock_idx[1]==k+4)&&mask[1]) continue; /* GLONASS */
            else if ((clock_idx[2]==k+4)&&mask[2]) continue; /* Galileo */
            else if ((clock_idx[3]==k+4)&&mask[3]) continue; /* BDS */
            else if ((clock_idx[4]==k+4)&&mask[4]) continue; /* IRN */
            else if ((clock_idx[5]==k+4)&&mask[5]) continue; /* QZS */

            v[nv]=0.0;
            for (j=0;j<nx;j++) {
                if ((clock_idx[1]==k+4)&&!mask[1]) H[j+nv*nx]=j==clock_idx[1]?1.0:0.0; 
                else if ((clock_idx[2]==k+4)&&!mask[2]) H[j+nv*nx]=j==clock_idx[2]?1.0:0.0;    
                else if ((clock_idx[3]==k+4)&&!mask[3]) H[j+nv*nx]=j==clock_idx[3]?1.0:0.0;  
                else if ((clock_idx[4]==k+4)&&!mask[4]) H[j+nv*nx]=j==clock_idx[4]?1.0:0.0;   
                else if ((clock_idx[5]==k+4)&&!mask[5]) H[j+nv*nx]=j==clock_idx[5]?1.0:0.0;  
            }
            var[nv++]=0.01;
        }        
    }


    free(v_); free(H_); free(var_);

    return nv;
}


/* validate solution ---------------------------------------------------------*/
extern int valsol(sol_t *sol, const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, double *P, int nv, int nx)
{
    double azels[MAXOBS*2],dop[4],vv,*vP;
    int i,ns;
    
    trace(3,"valsol  : n=%d nv=%d\n",n,nv);
    
    vP=mat(1,nv);

    /* large GDOP check (only for L1 satellites) */
    for (i=ns=0;i<n;i++) {
        if (!vsat[i]) continue;
        azels[  ns*2]=azel[  i*2];
        azels[1+ns*2]=azel[1+i*2];
        ns++;
    }
    /* compute DOP */
    dops(ns,azels,opt->elmin,dop);
    matcpy(sol->dop,dop,4,1);

    /* chi-square validation of residuals */
    matmul("TN",1,nv,nv,v,P,vP,1.0,0.0);
    matmul("NN",1,nv,1,vP,v,&vv,1.0,0.0);

#if 1    
    /* vv=dot(v,v,nv); */   
    if (nv>nx&&vv>chisqr[nv-nx-1]) {
        trace(7,"spp error: large chi-square error ns=%d vv=%.1f threshold=%.1f\n",nv,vv,chisqr[nv-nx-1]);
        free(vP);
        return 0; /* threshold too strict for all use cases, report error but continue on */
    }

    if (dop[0]<=0.0||dop[0]>MAX_GDOP) {
        trace(7,"gdop error nv=%d gdop=%.1f\n",nv,dop[0]);
        free(vP);
        return 0;
    }
#endif

    free(vP);

    return 1;
}

/* time update of position for spp/ins tc mode */
static void udpos_spp(rtk_t *rtk)
{
    prcopt_t *opt=&rtk->opt;
    ins_t *ins=&rtk->ins;
    int i,j,k,nx=6,*ix; /* nx=6, pos/vel */
    double p_ins[6],Cne[9],Cen[9],pos[3],Q[9]={0.0},Qe[9]={0.0},Qv[9]={0.0};
    double rr[3]={0.0},vel[3]={0.0};
    double *F,*P,*FP,*x,*xp,tt=rtk->interval,var=0.0;

    if (GINS_TC==opt->GI_mode||GINS_STC==opt->GI_mode) {
        /* convert INS solutions (pos and vel) to GNSS antenna center */
        ins2gnss(opt,ins,p_ins,6);
        pos2ecef(p_ins,rtk->ru);

        xyz2enu(p_ins,Cne);
        DCMT(Cne,Cen);
        Mat3mulv(1.0,Cen,p_ins+3,rtk->ru+3);        

        if (GINS_STC==opt->GI_mode) {
            /* nominal larger variance for position/velocity */
            /* for (i=0;i<3;i++) initx(rtk,rtk->ru[i],VAR_POS,i);
            for (i=3;i<6;i++) initx(rtk,rtk->ru[i],VAR_VEL,i); */  

            /* INS error state corresponding to covariance ------------
               transform local enu covariance to xyz-ecef covariance */
            for (i=0;i<3;i++) Q[i+i*3]=rtk->lcgins.P[(i+6)+(i+6)*rtk->lcgins.nx];
            covecef(p_ins,Q,Qe);
            for (i=0;i<3;i++) Q[i+i*3]=rtk->lcgins.P[(i+3)+(i+3)*rtk->lcgins.nx];
            covecef(p_ins,Q,Qv);
            for (i=0;i<3;i++) initx(rtk,rtk->ru[i],MAX(Qe[i+i*3],100.0),i);
            for (i=3;i<6;i++) initx(rtk,rtk->ru[i],MAX(Qv[(i-3)+(i-3)*3],10.0),i);  
            
            return;
        }

        /* NOTE: for tightly coupled mode, reset ins related state after feedback */
        if (GINS_TC==opt->GI_mode) reset_instat(rtk);

        return;  
    }
    else { 
#if 0
        /* white noise model */
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);    
#else
        /* spp based on kf (CV mode) */
        if (norm(rtk->x,3)<=0.0) { /* initialize position for first epoch */
            for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
            for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
            trace(7,"udpos_spp   : reset spp kf position due to gnss outage!\n");
            return;
        }

        /* check variance of estimated position */
        for (i=0;i<3;i++) var+=rtk->P[i+i*rtk->nx];
        var/=3.0;
        if (var>VAR_POS) {
            /* reset position with large variance */
            for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
            for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
            trace(7,"udpos_spp   : reset spp position due to large variance: var=%.3f\n",var);
            return;
        }

        ix=imat(6,1); /* generate valid state index (pos and vel) */
        for (i=0;i<nx;i++) {
            if (i<6||(rtk->x[i]!=0.0&&rtk->P[i+i*NX]>0.0)) ix[k++]=i;
        }
        /* state transition of position/velocity */
        F=eye(nx); P=mat(nx,nx); FP=mat(nx,nx); x=mat(nx,1); xp=mat(nx,1);
        
        /* improved CV model, which calculates the average velocity based on the velocities of the previous and current epochs */
        for (i=0;i<3;i++) {
            F[(i+3)+i*nx]=tt/2.0;
        }

        /* copy state */
        for (i=0;i<nx;i++) {
            x[i]=rtk->x[ix[i]];
            for (j=0;j<nx;j++) {
                P[i+j*nx]=rtk->P[ix[i]+ix[j]*NX];
            }
        }

        /* x=F*x, P=F*P*F'+Q */   
        matmul("NN",nx,nx,1,F,x,xp,1.0,0.0);  /* x=F*x */
        for (i=0;i<3;i++) xp[i]+=rtk->sol.rr_old[i+3]/2.0*tt; /* velocity-related input */
        matmul("NN",nx,nx,nx,F,P,FP,1.0,0.0); /* FP=F*P */
        matmul("NT",nx,nx,nx,FP,F,P,1.0,0.0); /* P=FP*F'*/ 
        
        /* update state */
        for (i=0;i<nx;i++) {
            rtk->x[ix[i]]=xp[i];
            for (j=0;j<nx;j++) {
                rtk->P[ix[i]+ix[j]*NX]=P[i+j*nx];
            }
        }

        /* process noise added to only velocity  P=P+Q */
        Q[0]=Q[4]=SQR(opt->prn[3])*fabs(tt);
        Q[8]=SQR(opt->prn[4])*fabs(tt);
        ecef2pos(rtk->x,pos);
        covecef(pos,Q,Qv);
        for (i=0;i<3;i++) for (j=0;j<3;j++) {
            rtk->P[j+i*NX]+=Qv[j+i*3];
            rtk->P[(j+3)+(i+3)*NX]+=Qv[j+i*3]; /* P=F*P*F'+Q*/
        }

        free(ix); free(F); free(P); free(FP); free(x); free(xp);
#endif

    }

}

/* time update of clock for spp/ins tc mode */
static void udclk_spp(rtk_t *rtk)
{
    prcopt_t *opt=&rtk->opt;
    double dtr;
    int i,ic,sys=rtk->opt.navsys,dclk=GNSISB_WN;
    trace(3,"udclk_spp:\n");

    /* initialize GPS clock (white noise) */
	ic=IC(0,opt);
    dtr=rtk->sol.dtr[0];
    if (fabs(dtr)<1.0e-10) dtr=1.0e-10;
    initx(rtk,CLIGHT*dtr,VAR_CLK,ic);


    /* initialize GPS drift (random walk) */
    /* NOTE: the receiver clock drift was modeled as a white noise model, and random walks did not perform well */
    ic=IC(6,opt);
    dtr=rtk->sol.dtr[6];
    if (dclk==GNSISB_WN) {
        /* white noise process */
        if (fabs(dtr)<1.0e-10) dtr=1.0e-10;
        initx(rtk,CLIGHT*dtr,VAR_CLKD,ic);
        } 
    else if (dclk==GNSISB_RW) {
        /* random walk process */
        if (rtk->x[ic]==0.0) {
            if (fabs(dtr)<1.0e-10) dtr=1.0e-10;
            initx(rtk,CLIGHT*dtr,VAR_CLKD,ic);
        }  
        else {
            rtk->P[ic+ic*rtk->nx]+=SQR(1e-2)*fabs(rtk->tt);
        }
    }   

    /* multi system clock initialization */
    for (i=1;i<NSYS;i++) {
        if (!(sys&SYS_GLO)&&i==1) continue;
        if (!(sys&SYS_GAL)&&i==2) continue;
        if (!(sys&SYS_CMP)&&i==3) continue;
        if (!(sys&SYS_IRN)&&i==4) continue;
        if (!(sys&SYS_QZS)&&i==5) continue;

        dtr=rtk->sol.dtr[i];
        ic=IC(i,opt);

        if (opt->sysisb==GNSISB_CT) {
            /* constant */
            if (rtk->x[ic]==0.0) {
                if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
                initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
            }
        }
        else if (opt->sysisb==GNSISB_RW) {
            /* random walk process */
            if (rtk->x[ic]==0.0) {
                if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
                initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
            }  
            else {
                rtk->P[ic+ic*rtk->nx]+=SQR(1e-3)*fabs(rtk->tt);
            }
        }
        else if (opt->sysisb==GNSISB_WN) {
            /* white noise process */
            if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
            initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
        }
    }
}

/* time update of states */
static void udstate_spp(rtk_t *rtk)
{
    /* time update of position */
    udpos_spp(rtk);

    /* time update of clock*/
    udclk_spp(rtk);
}


/* update solution status */
static void update_stat(rtk_t *rtk, int n, int stat)
{
    const prcopt_t *opt=&rtk->opt;
    ins_t *ins=&rtk->ins;
    sol_t *sol=&rtk->sol;
    int nx=rtk->nx,i;
    
    /* if GNSS/INS integration solution is available, reset GNSS outage count to 0 */
    if (rtk->outage<=MAX_OUTIME) rtk->outage=0;
    sol->ns=n;
    sol->stat=stat;
    
    /* update ins state */
    if (GINS_TC==opt->GI_mode) update_instat(opt,ins,rtk->P,sol,nx);
    else { /* spp kf */
        for (i=0;i<6;i++) sol->rr[i]=rtk->x[i];
    }

    /* store clock and isb */
    sol->dtr[0]=rtk->x[IC(0,opt)]/CLIGHT; /* GPS clock (s) */
    sol->dtr[1]=rtk->x[IC(1,opt)]/CLIGHT; /* GLO-GPS */
    sol->dtr[2]=rtk->x[IC(2,opt)]/CLIGHT; /* GAL-GPS */
    sol->dtr[3]=rtk->x[IC(3,opt)]/CLIGHT; /* BDS-GPS */
    sol->dtr[4]=rtk->x[IC(4,opt)]/CLIGHT; /* IRNSS-GPS */
    sol->dtr[5]=rtk->x[IC(5,opt)]/CLIGHT; /* QZSS-GPS */
    sol->dtr[6]=rtk->x[IC(6,opt)]/CLIGHT; /* GPS drift */
}

/* range rate residuals ------------------------------------------------------*/
static int resdop_filter(rtk_t *rtk, const obsd_t *obs, int iter, int m, int n, const double *rs, const double *dts,
                  const nav_t *nav, const double *rr, const double *x,
                  const double *azel, const int *vsat, double *v,
                  double *H, double *var, int mode, int *vi)
{
    ins_t *ins=&rtk->ins;
    prcopt_t *opt=&rtk->opt;
    double freq,rate,pos[3],a[3],e[3],vs[3],cosel,factor,delta;
    double Hva[3],Hvv[3],Cne[9],Cen[9],Cbn[9],wbie[3],wbeb[3],temp1[3],temp2[3],v_temp[9];
    int i,j,k,nv=0,sys,fr,nx=(GINS_TC==opt->GI_mode)?rtk->nx:((SPP_KF==mode)?NX:((SPP_LS_CD==mode)?n:4)),nf=opt->mfspp?opt->nf:1,idx; /* spp_ls_cd/spp+ins tc/spp kf */
    int spp_tc_flag=(GINS_TC==opt->GI_mode);
    
    trace(3,"resdop  : m=%d\n",m);
    
    /* check if there is enough Doppler observations */
    if (m<=0) return 0;

    /* the sign of Doppler observations is determined based on pseudorange variation between adjacent epochs */
    factor=rtk->dopsgn;

    ecef2pos(rr,pos); xyz2enu(pos,Cne); 
    DCMT(Cne,Cen);    

    for (k=0;k<nf;k++) {
        for (i=0;i<m&&i<MAXOBS;i++) {    
            sys=satsys(obs[i].sat,NULL); fr=sys2freid(sys,k,opt);
            freq=sat2freq(obs[i].sat,obs[i].code[fr],nav);

            idx=i+m*k;

            /* check if the observation is valid */
            if (obs[i].D[fr]==0.0||freq==0.0||!vsat[idx]||norm(rs+3+i*6,3)<=0.0) {
                continue;
            }

            /* calculate the line-of-sight vector based on satellite-to-ground distance in spp_LS_cd mode */
            if (SPP_LS_CD==mode) {
                geodist(rs+i*6,rr,e);
            }
            else {
                /* LOS (line-of-sight) vector in ECEF */
                cosel=cos(azel[1+i*2]);
                a[0]=sin(azel[i*2])*cosel;
                a[1]=cos(azel[i*2])*cosel;
                a[2]=sin(azel[1+i*2]);
                matmul("TN",3,3,1,Cne,a,e,1.0,0.0);                
            }

            if (iter>0) {
                /* test elevation mask */
                if (azel[1+i*2]<opt->elmin) continue;

                /* test SNR mask */
                if (!snrmask(0,obs+i,azel+i*2,opt)) continue;                
            }

            /* satellite velocity relative to receiver in ECEF */
            for (j=0;j<3;j++) {
                vs[j]=rs[j+3+i*6]-x[j];
            }
            /* range rate with earth rotation correction */
            rate=dot3(vs,e); /* neglect earth rotation correction (maximum quantity is mm/s ) */
            /* rate=dot3(vs,e)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                        rs[3+i*6]*rr[1]-rs[  i*6]*x[1]); */
            /* trace(7,"delta=%.7f\n",OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                        rs[3+i*6]*rr[1]-rs[  i*6]*x[1])); */                                    
            
            /* range rate residual (v=z-h(x)) (m/s) */
            v[nv]=(factor*obs[i].D[fr]*CLIGHT/freq-(rate-CLIGHT*dts[1+i*2]))-x[3];
            
            /* design matrix */
            if (GINS_TC==opt->GI_mode){
                vmulMat3(1.0,e,Cen,Hvv);

                DCMT(ins->Cnb,Cbn);
                Mat3mulv(1.0,Cbn,ins->eth.wnie,wbie);
                Mat3add2(ins->wbib,1.0,wbie,-1.0,wbeb);
                vskewmv(1.0,wbeb,ins->lever,temp1);
                Mat3mulv(1.0,ins->Cnb,temp1,temp2);
                
                vskew(1.0,temp2,v_temp);
                vmulMat3(1.0,Hvv,v_temp,Hva);

                for (j=0;j<nx;j++)  H[j+nv*nx]=(j==IC(6,opt)?1.0:0.0);
                for (j=0;j<3;j++)   H[j+nv*nx]=Hva[j];
                for (j=0;j+3<6;j++) H[(j+3)+nv*nx]=Hvv[j];            
            }
            else if (SPP_KF==mode) {
                for (j=0;j<nx;j++) {
                    H[j+nv*nx]=0.0;
                    if (j>=3&&j<6) H[j+nv*nx]=-e[j-3];
                    if (j==nx-1)   H[j+nv*nx]=1.0;
                }
            }
            else if (SPP_LS_CD==mode) {
                for (j=0;j<nx;j++)  {
                   H[j+nv*nx]=0.0; 
                   if (j>=nx-4&&j<nx-1) H[j+nv*nx]=-e[j-(nx-4)];
                   if (j==nx-1) H[j+nv*nx]=1.0;
                }
            }
            else {
                for (j=0;j<nx;j++)  H[j+nv*nx]=((j<3)?-e[j]:1.0);
            }

            if (vi) vi[nv]=idx;   

            /* TODO: stochastic model of doppler observations */
            var[nv++]=1e-1*varerr_spp(opt,NULL,&obs[i],azel[1+i*2],sys);
        }        
    }

    return nv;
}

/* estimate receiver velocity ------------------------------------------------*/
extern int estvel(rtk_t *rtk, const obsd_t *obs, int n, const double *rs, const double *dts,
                   const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                   const double *azel, int *vsat)
{
    double x[4]={0},dx[4],Q[16],*v,*H,*var,*P;
    double sig,err=opt->err[4],thres=2.0; /* Doppler error (Hz) */
    int i,j,k,nv,nx=4,stat=0,nf=opt->nf,info=0,LS_mode=opt->respp?Robust_RES:Robust_OFF;
    int *vi;
    
    v=mat(n*nf,1); H=mat(n*nf,nx); var=mat(n*nf,1); P=mat(n*nf,n*nf);
    vi=imat(n*nf,1);
    
    for (i=0;i<MAXITR;i++) {
        /* range rate residuals (m/s) */
        if ((nv=resdop_filter(rtk,obs,1,n,4,rs,dts,nav,sol->rr,x,azel,vsat,v,H,var,SPP_LS_D,vi))<4) {
            stat=0; 
            break;
        }

#if 1 
        /* outlier reject based on normal distribution*/       
        if (i>=2&&nv>=nx) {
            nv=outrej_spp(nv,NULL,nx,0,thres,v,H,var,NULL,NULL,vi,vsat,i,NULL);
        }
#endif

        /* weight by variance (lsq uses sqrt of weight) */
        diag_Cov(nv,var,P,diag_wei);

        /* least square estimation */
        if ((info=lsq_roubst(H,v,P,nx,nv,dx,Q,Robust_RES))) {
            trace(7,"dopple estvel lsq info=%d\n",info); break;
        }
        
        for (j=0;j<nx;j++) x[j]+=dx[j];
        
        if (norm(dx,nx)<1E-6) {
            trace(3,"estvel : vx=%.3f vy=%.3f vz=%.3f, n=%d\n",x[0],x[1],x[2],n);
            matcpy(sol->rr+3,x,3,1);
            sol->dtr[6]=x[3]/CLIGHT; /* clock drift */
            sol->qv[0]=(float)Q[0];  /* xx */
            sol->qv[1]=(float)Q[5];  /* yy */
            sol->qv[2]=(float)Q[10]; /* zz */
            sol->qv[3]=(float)Q[1];  /* xy */
            sol->qv[4]=(float)Q[6];  /* yz */
            sol->qv[5]=(float)Q[2];  /* zx */
            stat=1;
            break;
        }
    }

    free(v); free(H); free(var); free(P); free(vi);
    return stat;
}

extern int spp_sys(const prcopt_t *popt, int *clock_idx)
{
    int nx=3;

    /* GPS is included by default */
    clock_idx[0]=nx++;

    if (popt->navsys&SYS_GLO) {
        clock_idx[1]=nx++;
    }
    if (popt->navsys&SYS_GAL) {
        clock_idx[2]=nx++;
    }
    if (popt->navsys&SYS_CMP) {
        clock_idx[3]=nx++;
    }
    if (popt->navsys&SYS_IRN) {
        clock_idx[4]=nx++;
    }
    if (popt->navsys&SYS_QZS) {
        clock_idx[5]=nx++;
    }

    return nx;

}

/* determine the max obsat  */
extern int maxobsat(const int *ns, int nf)
{
    int j,max,*ns_;

    ns_=imat(nf,1); 
    for (j=0;j<nf;j++) ns_[j]=ns[j];

    /* determine the max obsat */
    for (j=1;j<nf;j++) {
        if (ns_[0]<ns_[j]) {
            max=ns_[0]; ns_[0]=ns_[j]; ns_[j]=max;
        }
    }
    max=ns_[0];

    free(ns_);

    return max;
}

static int dopvel_cons(rtk_t *rtk, int nx, int nx_code, const double *x_pre, const double *vel_pre, const double *x, double *v, double *H, double *var)
{
    int j,k,nv=0,interval=rtk->interval;

    for (j=0;j<3;j++) { 
        /* relative displacement = average velocity * sampling interval */
        v[nv]=(x_pre[j]+(vel_pre[j]+x[nx_code+j])/2.0*interval)-x[j];  
        for (k=0;k<nx;k++) H[k+nv*nx]=0.0;
        for (k=0;k<nx;k++) {
            if (j==k) {
                H[k+nv*nx]=1.0; 
                H[k+nx_code+nv*nx]=-interval/2.0; 
            }
        } 
        var[nv++]=5;
    }  

    return nv;
}

static void save_postat(const double *x, const double *Q, int nx, sol_t *sol, const int *mask, const int *clock_idx)
{
    int j;

    sol->type=0; /* not used ? */
    sol->dtr[0]=x[clock_idx[0]]/CLIGHT; /* receiver clock bias (s) */
    if (mask[1]) sol->dtr[1]=x[clock_idx[1]]/CLIGHT; /* GLO-GPS time offset (s) */
    if (mask[2]) sol->dtr[2]=x[clock_idx[2]]/CLIGHT; /* GAL-GPS time offset (s) */
    if (mask[3]) sol->dtr[3]=x[clock_idx[3]]/CLIGHT; /* BDS-GPS time offset (s) */
    if (mask[4]) sol->dtr[4]=x[clock_idx[4]]/CLIGHT; /* IRN-GPS time offset (s) */
    if (mask[5]) sol->dtr[5]=x[clock_idx[5]]/CLIGHT; /* QZS-GPS time offset (s) */               
    for (j=0;j<3;j++) sol->rr[j]=x[j];
    /* if (GINS_OFF==opt->GI_mode) for (j=0;j<3;j++) sol->rr[j+3]=0.0; */
    for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*nx];
    sol->qr[3]=(float)Q[1];    /* cov xy */
    sol->qr[4]=(float)Q[2+nx]; /* cov yz */
    sol->qr[5]=(float)Q[2];    /* cov zx */

    sol->age=sol->ratio=sol->ADOP=0.0; 
}

/* save vel state for spp based on pseudorange and doppler */
static void save_velstat(const double *x, const double *Q, int nx, sol_t *sol)
{
    int j;

    sol->dtr[6]=x[nx-1];
    for (j=0;j<3;j++) sol->rr_old[j+3]=sol->rr[j+3];
    for (j=0;j<3;j++) sol->rr[j+3]=x[nx-4+j];
    for (j=0;j<3;j++) sol->qv[j]=(float)Q[(j+nx-4)+(j+nx-4)*nx];
    sol->qv[3]=(float)Q[(nx-3)+(nx-4)*nx];    /* cov xy */
    sol->qv[4]=(float)Q[(nx-2)+(nx-3)*nx];    /* cov yz */
    sol->qv[5]=(float)Q[(nx-2)+(nx-4)*nx];    /* cov zx */
}

/* estimate receiver position ------------------------------------------------*/
extern int estpos(rtk_t *rtk, const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, ssat_t *ssat, sol_t *sol, double *azel,
                  int *vsat, double *resp)
{
    double x[NX]={0},dx[NX],Q[NX*NX],x_pre[3]={0.0},*v,*H,*var,sig;
    double *P,*R,interval=rtk->interval,thres=3.0,zupt_time,velvar_thres=5.0,dpos[3]={0.0},ave_vel[3]={0.0}; 
    double *xp,*Pp,vx[4];
    int i,j,k,m,info,nx=0,nx_code=0,stat=SOLQ_NONE,LS_mode=opt->respp?Robust_RES:Robust_OFF,mode,nv=0,*sati,*vi,nf=opt->mfspp?(IONOOPT_IFLC==opt->ionoopt?1:opt->nf):1,ns[nf];
    int max_sat,mask[NX-3]={0},clock_idx[NX-3]={0},spp_mode=opt->spp_mode,dop_cons_flag=0,spp_kf_flag=0;
    int nv_code=0,nv_dop=0,nv_cons=0,spp_tc_flag=0,f=(SPP_LS_CD==spp_mode)?2:1; /* spp/ins tc flag */
    
    trace(8,"estpos  : n=%d\n",n);
    
    /* init pvt status */
    if (GINS_TC==opt->GI_mode&&PMODE_SINGLE==opt->mode) spp_tc_flag=1; 
    else if (GINS_STC==opt->GI_mode&&PMODE_SINGLE==opt->mode) { spp_mode=SPP_KF; spp_kf_flag=1;}
    else if (SPP_KF==spp_mode) spp_kf_flag=1;

    /* number of spp parameters */
    nx=nx_code=spp_sys(opt,clock_idx);
    /* if the sign of doppler observations is no initialized, use SPP_LS_C mode */
    if (!rtk->dopsgn) spp_mode=SPP_LS_C;
    if (SPP_LS_CD==spp_mode) nx+=4; /* add dopple estimation parameters, ecef velocity and clock drift */

    v=mat((f*n+5)*nf,1); H=zeros((f*n+5)*nf,nx); var=mat((f*n+5)*nf,1); P=mat((f*n+5)*nf,(f*n+5)*nf);
    sati=imat((f*n+5)*nf,1); vi=imat((f*n+5)*nf,1);        
    
    /* use the previous epoch position as the initial position */
    for (i=0;i<3;i++) x[i]=sol->rr[i];
    /* dopple velocity constraint flag */
    if (SPP_LS_CD<=spp_mode&&!rtk->outage&&norm(x,3)>=0.0&&norm(sol->rr+3,3)>=0.0) {
        dop_cons_flag=1; 
        matcpy(x_pre,x,3,1);
    } 

    for (i=0;i<MAXITR;i++) {
        /* pseudorange residuals (m) */
        nv_code=nv=rescode(i,obs,n,nx_code,rs,dts,vare,svh,nav,x,opt,ssat,v,H,var,nx,mask,clock_idx,azel,vsat,resp,ns,sati,vi);    

        /* range rate residuals (m/s) */
        if (SPP_LS_CD==spp_mode) {
            nv+=resdop_filter(rtk,obs,i,n,nx,rs,dts,nav,x,x+(nx-4),azel,vsat,v+nv,H+nv*nx,var+nv,SPP_LS_CD,vi+nv);
        }

        /* outlier recject based on standard normal distribution */
        if (i>=2&&nv>=nx) {
            nv=outrej_spp(nv,&nv_code,nx,nx_code,thres,v,H,var,ssat,sati,vi,vsat,i,clock_idx);
        }

        /* dopple velocity constraint */
        if (SPP_LS_CD==spp_mode&&dop_cons_flag&&i>2) {
            /* if the variance of the velocity estimate is too large, doppler constraints should not be used */
            if ((Q[(nx-4)+(nx-4)*nx]+Q[(nx-3)+(nx-3)*nx]+Q[(nx-2)+(nx-2)*nx])/3.0>velvar_thres) {
                dop_cons_flag=0;
            }
            if (dop_cons_flag) {
                nv+=dopvel_cons(rtk,nx,nx_code,x_pre,sol->rr+3,x,v+nv,H+nv*nx,var+nv);
            }
        }

        /* determine the max obsat of all frequency */
        max_sat=maxobsat(ns,nf);

        if (max_sat<nx_code||nv_code<nx_code||nv<nx) {
            trace(7,"spp lack of valid sats, nv=%d\n",max_sat); break;
        }

        /* weight by variance */
        diag_Cov(nv,var,P,diag_wei);

        /* least square estimation */
        if ((info=lsq_roubst(H,v,P,nx,nv,dx,Q,LS_mode))) {
            trace(7,"spp lsq error info=%d\n!",info); break;
        }

        /* tracefilter(12,TRAE_R|TRAE_H|TRAE_v|TRAE_xpre,nx,nv,P,H,NULL,NULL,v,dx,NULL); */

        for (j=0;j<nx;j++) {
            x[j]+=dx[j];
        }
        if (norm(dx,nx)<1E-4) {
            sol->time=timeadd(obs[0].time,-x[3]/CLIGHT); /* receiver time */
            sol->ns=(uint8_t)ns[0];
            save_postat(x,Q,nx,sol,mask,clock_idx);
            /* spp based on pseudorange and doppler */
            if (SPP_LS_CD==spp_mode) save_velstat(x,Q,nx,sol);
            /* for GINS mode, may spp is ok but ppk is not ok */
            if (isGNSS(opt)) rtk->outage=0; /* reset outage counter */
            
            /* validate solution */
            if ((stat=valsol(sol,azel,vsat,n,opt,v,P,nv,nx))) {
                sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
                /* save receiver clock (m) */
                for (j=0;j<n;j++) if (ssat) ssat[obs[j].sat-1].cdtr[0]=x[3];
            }
            
            /* free memory */
            free(v); free(H); free(var); free(P); free(sati); free(vi);
            
            if (spp_tc_flag||spp_kf_flag) break;
            else return stat;
        }
    }

    /* if the iteration exceeds the limit or the solution fails, the solution fails flag is returned */
    if (i>=MAXITR||(SOLQ_NONE==stat&&!spp_tc_flag)) {
        if (i>=MAXITR) trace(7,"spp: iteration over limit i=%d, norm(dx)=%.6f!\n",i,norm(dx,nx));              
        free(v); free(H); free(var); free(P); free(sati); free(vi);  
        return SOLQ_NONE;
    }
    
    /* spp/ins tc integration */
    if (rtk&&(spp_tc_flag||spp_kf_flag)) {
        nx=(spp_tc_flag)?rtk->nx:NX; /* number of states parameters */
        mode=rtk->opt.filter; /* fusion filter mode */
        /* initialize receiver velocity and clock drift */
        if (spp_kf_flag) estvel(rtk,obs,n,rs,dts,nav,opt,sol,azel,vsat);

        /* if GNSS solution fails, do not enable GNSS/INS integration mode */
        if (stat) {
            /* initialization, consider motion constraints (NHC/ZUPT/ZIHR) */
            nv=2*n*nf+3+1; /* ZUPT+ZIHR */
            xp=zeros(nx,1); Pp=zeros(nx,nx); v=mat(nv,1); H=mat(nv,nx); var=mat(nv,1); R=zeros(nv,nv); sati=imat(2*n,1);

            /* time update of ekf states*/
            udstate_spp(rtk);

            /* copy states */
            matcpy(xp,rtk->x,nx,1); matcpy(Pp,rtk->P,nx,nx);

            /* initialize velocity-related state parameters */
            for (i=0;i<3;i++) {
                if (spp_tc_flag) vx[i]=rtk->ru[i+3]; /* spp/in tc */
                else vx[i]=rtk->x[i+3];              /* spp */    
            }
            vx[3]=rtk->x[IC(6,opt)];

            /* prefit residuals (v=z-h(x))*/
            nv=rescode_filter(rtk,obs,n,rs,dts,vare,svh,nav,(spp_tc_flag?rtk->ru:xp),opt,ssat,v,H,var,azel,vsat,resp,ns,sati);               
            /* doppler obs */
            if (rtk->dopsgn) {
                nv_dop=resdop_filter(rtk,obs,1,n,-1,rs,dts,nav,(spp_tc_flag?rtk->ru:xp),vx,azel,vsat,v+nv,H+nv*nx,var+nv,(spp_tc_flag?GINS_TC:SPP_KF),NULL);                
            }
            
            /* motion constraints (nhc/zupt/zihr) */
            if (spp_tc_flag) {  
                nv_cons=motion_meas(rtk,opt,H,v,var,nv+nv_dop,nx);                  
            }      

            /* measurement noise covariance matrix R */
            diag_Cov(nv+nv_dop+nv_cons,var,R,diag_var);

            /* kalman filter measurement update */
            if ((info=filter_gins(rtk,xp,Pp,H,v,R,nx,(nv+nv_dop+nv_cons),(spp_tc_flag?KF_GINS:KF_GNSS),mode))) {
                trace(7,"SPP/INS filter error (info=%d)\n",info);
                free(xp); free(Pp); free(v); free(H); free(R); free(var); free(sati);
                return SOLQ_NONE;
            };                

            /* save receiver clock drift (m/s) */
            for (j=0;j<n;j++) if (ssat) ssat[obs[j].sat-1].cdtr[1]=xp[IC(6,opt)];

            /* tracefilter(12,TRAE_R|TRAE_H|TRAE_Ppre|TRAE_Pp|TRAE_v|TRAE_xpre|TRAE_xp,nx,nv+nv_dop+nv_cons,R,H,rtk->P,Pp,v,rtk->x,xp); */

            /* updates states */
            matcpy(rtk->x,xp,nx,1); matcpy(rtk->P,Pp,nx,nx);

            /* reset cross-covariance */
            /* init_crosscov(rtk,rtk->ins.nx,rtk->nx); */

            /* ins feedback correction and state reset */
            if (spp_tc_flag) {
                ins_fedback(rtk,xp);
                reset_instat(rtk);
            }

            /* update solution status */
            update_stat(rtk,nv,SOLQ_SINGLE);    
                            
            free(xp); free(Pp); free(v); free(H); free(R); free(var); free(sati);
            return stat;
        }
        /* motion constraints (nhc/zupt/zihr) */
        else if (GINS_TC==opt->GI_mode&&is_motionconstraints(opt)) {
            rtk->outage++;
            motion_constraints(rtk,opt);
            return SOLQ_CONS;
        }
        else {
            return SOLQ_NONE;  
        }
    }
}
/* RAIM FDE (failure detection and exclusion) -------------------------------*/
static int raim_fde(const obsd_t *obs, int n, const double *rs,
                    const double *dts, const double *vare, const int *svh,
                    const nav_t *nav, const prcopt_t *opt, ssat_t *ssat, 
                    sol_t *sol, double *azel, int *vsat, double *resp)
{
    obsd_t *obs_e;
    sol_t sol_e={{0}};
    char tstr[32],name[16];
    double *rs_e,*dts_e,*vare_e,*azel_e,*resp_e,rms_e,rms=100.0;
    int i,j,k,nvsat,stat=0,*svh_e,*vsat_e,sat=0;
    
    trace(3,"raim_fde: %s n=%2d\n",time_str(obs[0].time,0),n);
    
    if (!(obs_e=(obsd_t *)malloc(sizeof(obsd_t)*n))) return 0;
    rs_e = mat(6,n); dts_e = mat(2,n); vare_e=mat(1,n); azel_e=zeros(2,n);
    svh_e=imat(1,n); vsat_e=imat(1,n); resp_e=mat(1,n); 
    
    for (i=0;i<n;i++) {       
        /* satellite exclusion */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            obs_e[k]=obs[j];
            matcpy(rs_e +6*k,rs +6*j,6,1);
            matcpy(dts_e+2*k,dts+2*j,2,1);
            vare_e[k]=vare[j];
            svh_e[k++]=svh[j];
        }
        /* estimate receiver position without a satellite */
        if (!estpos(NULL,obs_e,n-1,rs_e,dts_e,vare_e,svh_e,nav,opt,ssat,&sol_e,azel_e,
                    vsat_e,resp_e)) {
            trace(3,"raim_fde: exsat=%2d (%s)\n",obs[i].sat);
            continue;
        }
        for (j=nvsat=0,rms_e=0.0;j<n-1;j++) {
            if (!vsat_e[j]) continue;
            rms_e+=SQR(resp_e[j]);
            nvsat++;
        }
        if (nvsat<5) {
            trace(3,"raim_fde: exsat=%2d lack of satellites nvsat=%2d\n",
                  obs[i].sat,nvsat);
            continue;
        }
        rms_e=sqrt(rms_e/nvsat);
        
        trace(3,"raim_fde: exsat=%2d rms=%8.3f\n",obs[i].sat,rms_e);
        
        if (rms_e>rms) continue;
        
        /* save result */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            matcpy(azel+2*j,azel_e+2*k,2,1);
            vsat[j]=vsat_e[k];
            resp[j]=resp_e[k++];
        }
        stat=1;
        *sol=sol_e;
        sat=obs[i].sat;
        rms=rms_e;
        vsat[i]=0;
    }
    if (stat) {
        time2str(obs[0].time,tstr,2); satno2id(sat,name);
        trace(2,"%s: %s excluded by raim\n",tstr+11,name);
    }
    free(obs_e);
    free(rs_e ); free(dts_e ); free(vare_e); free(azel_e);
    free(svh_e); free(vsat_e); free(resp_e);
    return stat;
}

/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          prcopt_t *opt    I   processing options
*          sol_t  *sol      IO  solution
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
*          ssat_t *ssat     IO  satellite status              (NULL: no output)
*          char   *msg      O   error message for error exit
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int pntpos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, ssat_t *ssat)
{
    prcopt_t opt_=*opt;
    double *rs,*dts,*var,*resp,*azel_;
    int i,j,stat,*vsat,svh[MAXOBS],sys,fr,nf=opt->mfspp?(IONOOPT_IFLC==opt_.ionoopt?1:opt_.nf):1;
    
    trace(3,"pntpos  : tobs=%s n=%d\n",time_str(obs[0].time,3),n);
    
    /* NOTE: for GNSS/INS TC integration, the INS solution is set to the initial solution status */
    if (GINS_TC==opt_.GI_mode) sol->stat=SOLQ_INS;
    else sol->stat=SOLQ_NONE; 
    
    if (n<=0) {
        rtk->outage++;
        /* NOTE: reset spp position if gnss outage in spp_kf mode */
        if (SPP_KF==opt_.spp_mode) for (i=0;i<6;i++) rtk->x[i]=0.0;
        /* if the number of available satellites is 0, output INS solution */
        if (GINS_TC==opt_.GI_mode) update_instat(&rtk->opt,&rtk->ins,rtk->P,sol,rtk->nx);                  
        trace(7,"no observation data");
        return 0;
    }
    sol->time=obs[0].time;
    
    rs=mat(n,6); dts=mat(n,2); var=mat(n*nf,2); azel_=zeros(n,2); resp=mat(n*nf,1);
    vsat=imat(n*nf,1);
    
    /* init ssat struct */
    if (rtk&&ssat) init_ssatpar(rtk,obs,n,SPP_ssat,SOLQ_NONE);
    
    if (opt_.mode!=PMODE_SINGLE||opt_.GI_mode!=GINS_OFF) { /* for precise positioning */
        opt_.spp_mode=SPP_LS_C; /* TODO */
        opt_.sateph=EPHOPT_BRDC;
        opt_.ionoopt=IONOOPT_BRDC;
        opt_.tropopt=TROPOPT_SAAS;
    }
    /* satellite positions, velocities and clocks */
    satposs(sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);
    
    /* estimate receiver position and time with pseudorange */
    stat=estpos(rtk,obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp);

    /* spp/ins tc mode and GNSS unavailable, output INS solution */
    if (!stat&&rtk) {
        rtk->outage++;
        /* NOTE: reset spp position if gnss outage in spp_kf mode */
        if (SPP_KF==opt_.spp_mode) for (i=0;i<6;i++) rtk->x[i]=0.0;
        /* in the GNSS framework, only TC is processed; LC/STC is processed in lc_gins */
        if (GINS_TC==opt->GI_mode) { /* for ppd/ppk/ppp-ins tc mode */
            if (is_motionconstraints(opt)) {
                motion_constraints(rtk,opt);
            }
            else {
                sol->stat=SOLQ_INS;
                update_instat(&rtk->opt,&rtk->ins,rtk->P,sol,rtk->nx);                 
            }          
        }
        return SOLQ_NONE;
    }

    /* RAIM FDE */
    if (!stat&&n>=6&&opt->posopt[4]) {
        stat=raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp);
    }

    /* estimate receiver velocity with Doppler */
    if (isGNSS(opt)&&(SPP_LS_C==opt->spp_mode)&&stat&&rtk->dopsgn) {
        estvel(rtk,obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);
    }

    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];
    }

    free(rs); free(dts); free(var); free(azel_); free(resp); free(vsat);

    return stat;
}
