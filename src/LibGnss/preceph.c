/*------------------------------------------------------------------------------
* preceph.c : precise ephemeris and clock functions
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* references :
*     [1] S.Hilla, The Extended Standard Product 3 Orbit Format (SP3-c),
*         12 February, 2007
*     [2] J.Ray, W.Gurtner, RINEX Extensions to Handle Clock Information,
*         27 August, 1998
*     [3] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*     [4] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*         Space Technology Library, 2004
*     [5] S.Hilla, The Extended Standard Product 3 Orbit Format (SP3-d),
*         February 21, 2016
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2009/01/18 1.0  new
*           2009/01/31 1.1  fix bug on numerical error to read sp3a ephemeris
*           2009/05/15 1.2  support glonass,galileo,qzs
*           2009/12/11 1.3  support wild-card expansion of file path
*           2010/07/21 1.4  added api:
*                               eci2ecef(),sunmoonpos(),peph2pos(),satantoff(),
*                               readdcb()
*                           changed api:
*                               readsp3()
*                           deleted api:
*                               eph2posp()
*           2010/09/09 1.5  fix problem when precise clock outage
*           2011/01/23 1.6  support qzss satellite code
*           2011/09/12 1.7  fix problem on precise clock outage
*                           move sunmmonpos() to rtkcmn.c
*           2011/12/01 1.8  modify api readsp3()
*                           precede later ephemeris if ephemeris is NULL
*                           move eci2ecef() to rtkcmn.c
*           2013/05/08 1.9  fix bug on computing std-dev of precise clocks
*           2013/11/20 1.10 modify option for api readsp3()
*           2014/04/03 1.11 accept extenstion including sp3,eph,SP3,EPH
*           2014/05/23 1.12 add function to read sp3 velocity records
*                           change api: satantoff()
*           2014/08/31 1.13 add member cov and vco in peph_t sturct
*           2014/10/13 1.14 fix bug on clock error variance in peph2pos()
*           2015/05/10 1.15 add api readfcb()
*                           modify api readdcb()
*           2017/04/11 1.16 fix bug on antenna offset correction in peph2pos()
*           2020/11/30 1.17 support SP3-d [5] to accept more than 85 satellites
*                           support NavIC/IRNSS in API peph2pos()
*                           LC defined GPS/QZS L1-L2, GLO G1-G2, GAL E1-E5b,
*                            BDS B1I-B2I and IRN L5-S for API satantoff()
*                           fix bug on reading SP3 file extension
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define NMAX        10              /* order of polynomial interpolation */
#define MAXDTE      900.0           /* max time difference to ephem time (s) */
#define EXTERR_CLK  1E-3            /* extrapolation error for clock (m/s) */
#define EXTERR_EPH  5E-7            /* extrapolation error for ephem (m/s^2) */
#define MAX_BIAS_SYS 5              /* # of constellations supported */

/* table to translate code to code bias table index  */
static int8_t code_Obias_ix[MAX_BIAS_SYS][MAXCODE];
/* initialize code bias lookup table -------------------------------------------
*       -1 = code not supported
*        0 = reference code (0 bias)
*        1-3 = table index for code
* ----------------------------------------------------------------------------*/
static void init_bias_ix(void) {
    int i,j,k;

    for (i=0;i<MAX_BIAS_SYS;i++) for (j=0;j<MAXCODE;j++) {
        code_Obias_ix[i][j]=-1;
    }
        
    /*OSB*/
    /* GPS */
    /* satellite signals with accurately DCB corrections based on .BSX file */
    code_Obias_ix[0][CODE_L1C]=0;
    code_Obias_ix[0][CODE_L1W]=1;
    code_Obias_ix[0][CODE_L2W]=2;
    code_Obias_ix[0][CODE_L2L]=3;
    code_Obias_ix[0][CODE_L2S]=4;
    code_Obias_ix[0][CODE_L2X]=5;
    code_Obias_ix[0][CODE_L2C]=6;
    code_Obias_ix[0][CODE_L5Q]=7;
    code_Obias_ix[0][CODE_L5X]=8;

    /* satellite signals with approximately DCB corrections based on .BSX file */
    code_Obias_ix[0][CODE_L1S]=1;
    code_Obias_ix[0][CODE_L1L]=1;
    code_Obias_ix[0][CODE_L1X]=1;
    code_Obias_ix[0][CODE_L1P]=1;
    code_Obias_ix[0][CODE_L1Y]=1;
    code_Obias_ix[0][CODE_L1M]=1;

    code_Obias_ix[0][CODE_L2D]=2;
    code_Obias_ix[0][CODE_L2P]=2;
    code_Obias_ix[0][CODE_L2Y]=2;
    code_Obias_ix[0][CODE_L2M]=2;
    code_Obias_ix[0][CODE_L5I]=6;

    /* GLONASS */
    code_Obias_ix[1][CODE_L1C]=0;
    code_Obias_ix[1][CODE_L1P]=1;
    code_Obias_ix[1][CODE_L2C]=2;
    code_Obias_ix[1][CODE_L2P]=3;

    /* Galileo */
    /* satellite signals with accurately DCB corrections based on .BSX file */
    code_Obias_ix[2][CODE_L1C]=0;
    code_Obias_ix[2][CODE_L1X]=1;
    code_Obias_ix[2][CODE_L5Q]=2;
    code_Obias_ix[2][CODE_L5X]=3;
    code_Obias_ix[2][CODE_L6C]=4;
    code_Obias_ix[2][CODE_L7Q]=5;
    code_Obias_ix[2][CODE_L7X]=6;
    code_Obias_ix[2][CODE_L8Q]=7;   
    code_Obias_ix[2][CODE_L8X]=8;

    /* satellite signals with approximately DCB corrections based on .BSX file */
    code_Obias_ix[2][CODE_L1A]=1;
    code_Obias_ix[2][CODE_L1B]=1;
    code_Obias_ix[2][CODE_L1Z]=1;
    code_Obias_ix[2][CODE_L5I]=3;
    code_Obias_ix[2][CODE_L7I]=6;
    code_Obias_ix[2][CODE_L8I]=8;
    code_Obias_ix[2][CODE_L6A]=4;
    code_Obias_ix[2][CODE_L6B]=4;
    code_Obias_ix[2][CODE_L6X]=4;
    code_Obias_ix[2][CODE_L6Z]=4;

    /* Beidou */
    /* satellite signals with accurately DCB corrections based on .BSX file */
    code_Obias_ix[3][CODE_L2I]=0;
    code_Obias_ix[3][CODE_L6I]=1;
    code_Obias_ix[3][CODE_L7I]=2;
    code_Obias_ix[3][CODE_L1D]=3;
    code_Obias_ix[3][CODE_L1P]=4;
    code_Obias_ix[3][CODE_L1X]=5;
    code_Obias_ix[3][CODE_L5D]=6;
    code_Obias_ix[3][CODE_L5P]=7;
    code_Obias_ix[3][CODE_L5X]=8;
    code_Obias_ix[3][CODE_L7Z]=9;
    code_Obias_ix[3][CODE_L8X]=10;

    /* satellite signals with approximately DCB corrections based on .BSX file */
    code_Obias_ix[3][CODE_L2Q]=0;
    code_Obias_ix[3][CODE_L2X]=0;
    code_Obias_ix[3][CODE_L6Q]=1;
    code_Obias_ix[3][CODE_L6X]=1;
    code_Obias_ix[3][CODE_L7Q]=2;
    code_Obias_ix[3][CODE_L7X]=2;
    code_Obias_ix[3][CODE_L7D]=7;
    code_Obias_ix[3][CODE_L7P]=7;
    code_Obias_ix[3][CODE_L8D]=8;
    code_Obias_ix[3][CODE_L8P]=8;

    /* QZSS */
    /* satellite signals with accurately DCB corrections based on .BSX file */
    code_Obias_ix[4][CODE_L1C]=0;
    code_Obias_ix[4][CODE_L1X]=1;
    code_Obias_ix[4][CODE_L2L]=2;
    code_Obias_ix[4][CODE_L2X]=3;
    code_Obias_ix[4][CODE_L5Q]=4;
    code_Obias_ix[4][CODE_L5X]=5;
  
}

/* satellite code to satellite system ----------------------------------------*/
static int code2sys(char code)
{
    if (code=='G'||code==' ') return SYS_GPS;
    if (code=='R') return SYS_GLO;
    if (code=='E') return SYS_GAL; /* SP3-d */
    if (code=='J') return SYS_QZS; /* SP3-d */
    if (code=='C') return SYS_CMP; /* SP3-d */
    if (code=='I') return SYS_IRN; /* SP3-d */
    if (code=='L') return SYS_LEO; /* SP3-d */
    return SYS_NONE;
}
/* read SP3 header -----------------------------------------------------------*/
static int readsp3h(FILE *fp, gtime_t *time, char *type, int *sats,
                    double *bfact, char *tsys)
{
    int i=0,j,k=0,ns=0,sys,prn;
    char buff[1024];

    trace(3,"readsp3h:\n");

    for (i=0;;i++) {
        if (!fgets(buff,sizeof(buff),fp)) break;
        
        if (i==0) {
            *type=buff[2];
            if (str2time(buff,3,28,time)) return 0;
        }
        else if (!strncmp(buff,"+ ",2)) { /* satellite id */
            if (ns==0) {
                ns=(int)str2num(buff,3,3);
            }
            for (j=0;j<17&&k<ns;j++) {
                sys=code2sys(buff[9+3*j]);
                prn=(int)str2num(buff,10+3*j,2);
                if (k<MAXSAT) sats[k++]=satno(sys,prn);
            }
        }
        else if (!strncmp(buff,"++",2)) { /* orbit accuracy */
            continue;
        }
        else if (!strncmp(buff,"%c",2)&&tsys[0]=='\0') { /* time system */
            strncpy(tsys,buff+9,3); tsys[3]='\0';
        }
        else if (!strncmp(buff,"%f",2)&&bfact[0]==0.0) { /* fp base number */
            bfact[0]=str2num(buff, 3,10);
            bfact[1]=str2num(buff,14,12);
        }
        else if (!strncmp(buff,"%i",2)) {
            continue;
        }
        else if (!strncmp(buff,"/*",2)) { /* comment */
            continue;
        }
        else if (!strncmp(buff,"* ",2)) { /* first record */
            /* roll back file pointer */
            fseek(fp,0,SEEK_SET);
            break;
        }
    }
    return ns;
}
/* add precise ephemeris -----------------------------------------------------*/
static int addpeph(nav_t *nav, peph_t *peph)
{
    peph_t *nav_peph;

    if (nav->ne>=nav->nemax) {
        nav->nemax+=256;
        if (!(nav_peph=(peph_t *)realloc(nav->peph,sizeof(peph_t)*nav->nemax))) {
            trace(1,"readsp3b malloc error n=%d\n",nav->nemax);
            free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
            return 0;
        }
        nav->peph=nav_peph;
    }
    nav->peph[nav->ne++]=*peph;
    return 1;
}
/* read SP3 body -------------------------------------------------------------*/
static void readsp3b(FILE *fp, char type, int *sats, int ns, double *bfact,
                     char *tsys, int index, int opt, nav_t *nav)
{
    peph_t peph;
    gtime_t time;
    double val,std,base;
    int i,j,sat,sys,prn,n=ns*(type=='P'?1:2),pred_o,pred_c,v;
    char buff[1024];

    trace(3,"readsp3b: type=%c ns=%d index=%d opt=%d\n",type,ns,index,opt);

    while (fgets(buff,sizeof(buff),fp)) {

        if (!strncmp(buff,"EOF",3)) break;

        if (buff[0]!='*'||str2time(buff,3,28,&time)) {
            /* trace(2,"sp3 invalid epoch %31.31s\n",buff); */
            continue;
        }
        if (!strcmp(tsys,"UTC")) time=utc2gpst(time); /* utc->gpst */
        peph.time =time;
        peph.index=index;

        for (i=0;i<MAXSAT;i++) {
            for (j=0;j<4;j++) {
                peph.pos[i][j]=0.0;
                peph.std[i][j]=0.0f;
                peph.vel[i][j]=0.0;
                peph.vst[i][j]=0.0f;
            }
            for (j=0;j<3;j++) {
                peph.cov[i][j]=0.0f;
                peph.vco[i][j]=0.0f;
            }
        }
        for (i=pred_o=pred_c=v=0;i<n&&fgets(buff,sizeof(buff),fp);i++) {

            if (strlen(buff)<4||(buff[0]!='P'&&buff[0]!='V')) continue;

            sys=buff[1]==' '?SYS_GPS:code2sys(buff[1]);
            prn=(int)str2num(buff,2,2);
            if      (sys==SYS_SBS) prn+=100;
            else if (sys==SYS_QZS) prn+=192; /* extension to sp3-c */

            if (!(sat=satno(sys,prn))) continue;

            if (buff[0]=='P') {
                pred_c=strlen(buff)>=76&&buff[75]=='P';
                pred_o=strlen(buff)>=80&&buff[79]=='P';
            }
            for (j=0;j<4;j++) {

                /* read option for predicted value */
                if (j< 3&&(opt&1)&& pred_o) continue;
                if (j< 3&&(opt&2)&&!pred_o) continue;
                if (j==3&&(opt&1)&& pred_c) continue;
                if (j==3&&(opt&2)&&!pred_c) continue;

                val=str2num(buff, 4+j*14,14);
                std=str2num(buff,61+j* 3,j<3?2:3);

                if (buff[0]=='P') { /* position */
                    if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
                        peph.pos[sat-1][j]=val*(j<3?1000.0:1E-6);
                        v=1; /* valid epoch */
                    }
                    if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
                        peph.std[sat-1][j]=(float)(pow(base,std)*(j<3?1E-3:1E-12));
                    }
                }
                else if (v) { /* velocity */
                    if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
                        peph.vel[sat-1][j]=val*(j<3?0.1:1E-10);
                    }
                    if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
                        peph.vst[sat-1][j]=(float)(pow(base,std)*(j<3?1E-7:1E-16));
                    }
                }
            }
        }
        if (v) {
            if (!addpeph(nav,&peph)) return;
        }
    }
}
/* compare precise ephemeris -------------------------------------------------*/
static int cmppeph(const void *p1, const void *p2)
{
    peph_t *q1=(peph_t *)p1,*q2=(peph_t *)p2;
    double tt=timediff(q1->time,q2->time);
    return tt<-1E-9?-1:(tt>1E-9?1:q1->index-q2->index);
}
/* combine precise ephemeris -------------------------------------------------*/
static void combpeph(nav_t *nav, int opt)
{
    int i,j,k,m;

    trace(3,"combpeph: ne=%d\n",nav->ne);

    qsort(nav->peph,nav->ne,sizeof(peph_t),cmppeph);

    if (opt&4) return;

    for (i=0,j=1;j<nav->ne;j++) {

        if (fabs(timediff(nav->peph[i].time,nav->peph[j].time))<1E-9) {

            for (k=0;k<MAXSAT;k++) {
                if (norm(nav->peph[j].pos[k],4)<=0.0) continue;
                for (m=0;m<4;m++) nav->peph[i].pos[k][m]=nav->peph[j].pos[k][m];
                for (m=0;m<4;m++) nav->peph[i].std[k][m]=nav->peph[j].std[k][m];
                for (m=0;m<4;m++) nav->peph[i].vel[k][m]=nav->peph[j].vel[k][m];
                for (m=0;m<4;m++) nav->peph[i].vst[k][m]=nav->peph[j].vst[k][m];
            }
        }
        else if (++i<j) nav->peph[i]=nav->peph[j];
    }
    nav->ne=i+1;

    trace(3,"combpeph: ne=%d\n",nav->ne);
}
/* read sp3 precise ephemeris file ---------------------------------------------
* read sp3 precise ephemeris/clock files and set them to navigation data
* args   : char   *file       I   sp3-c precise ephemeris file
*                                 (wind-card * is expanded)
*          nav_t  *nav        IO  navigation data
*          int    opt         I   options (1: only observed + 2: only predicted +
*                                 4: not combined)
* return : none
* notes  : see ref [1]
*          precise ephemeris is appended and combined
*          nav->peph and nav->ne must by properly initialized before calling the
*          function
*          only files with extensions of .sp3, .SP3, .eph* and .EPH* are read
*-----------------------------------------------------------------------------*/
extern void readsp3(const char *file, nav_t *nav, int opt)
{
    FILE *fp;
    gtime_t time={0};
    double bfact[2]={0};
    int i,j,n,ns,sats[MAXSAT]={0};
    char *efiles[MAXEXFILE],*ext,type=' ',tsys[4]="";

    trace(3,"readpephs: file=%s\n",file);

    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return;
        }
    }
    /* expand wild card in file path */
    n=expath(file,efiles,MAXEXFILE);

    for (i=j=0;i<n;i++) {
        if (!(ext=strrchr(efiles[i],'.'))) continue;

        if (!(fp=fopen(efiles[i],"r"))) {
            trace(1,"sp3 file open error %s!\n",efiles[i]);
            continue;
        }
        /* read sp3 header */
        ns=readsp3h(fp,&time,&type,sats,bfact,tsys);

        /* read sp3 body */
        readsp3b(fp,type,sats,ns,bfact,tsys,j++,opt,nav);

        fclose(fp);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);

    /* combine precise ephemeris */
    if (nav->ne>0) combpeph(nav,opt);
}
/* read satellite antenna parameters -------------------------------------------
* read satellite antenna parameters
* args   : char   *file       I   antenna parameter file
*          gtime_t time       I   time
*          nav_t  *nav        IO  navigation data
* return : status (1:ok,0:error)
* notes  : only support antex format for the antenna parameter file
*-----------------------------------------------------------------------------*/
extern int readsap(const char *file, gtime_t time, nav_t *nav)
{
    spcvs_t pcvs={0};
    rpcvs_t pcvr={0};
    spcv_t pcv0={0},*pcv;
    int i;

    trace(3,"readsap : file=%s time=%s\n",file,time_str(time,0));

    if (!readpcv(file,&pcvs,&pcvr)) return 0;

    for (i=0;i<MAXSAT;i++) {
        pcv=searchspcv(i+1,time,&pcvs);
        nav->spcvs[i]=pcv?*pcv:pcv0;
    }
    free(pcvs.pcv);
    return 1;
}
/* read DCB parameters from DCB file -------------------------------------------
*    - supports satellite and receiver biases
*-----------------------------------------------------------------------------*/
static int readdcbf(const char *file, nav_t *nav, const sta_t *sta)
{
    FILE *fp;
    double cbias;
    char buff[256],str1[32],str2[32]="";
    int i,j,sat,type=0;

    trace(3,"readdcbf: file=%s\n",file);

    if (!(fp=fopen(file,"r"))) {
        trace(1,"GPS dcb file open error: %s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {

        if      (strstr(buff,"DIFFERENTIAL (P1-C1) CODE BIASES")) type=1;
        else if (strstr(buff,"DIFFERENTIAL (P2-C2) CODE BIASES")) type=2;

        str2[0]='\0';
        if (!type||sscanf(buff,"%31s %31s",str1,str2)<1) continue;

        if ((cbias=str2num(buff,26,9))==0.0) continue;

        if (sta&&(!strcmp(str1,"G")||!strcmp(str1,"R"))) { /* receiver DCB */
            for (i=0;i<MAXRCV;i++) {
                if (!strcmp(sta[i].name,str2)) break;
            }
            if (i<MAXRCV) {
                j=!strcmp(str1,"G")?0:1;
                nav->rbias[i][j][type-1]=cbias*1E-9*CLIGHT; /* ns -> m */
            }
        }
        else if ((sat=satid2no(str1))) { /* satellite dcb */
            nav->cbias[sat-1][type-1]=cbias*1E-9*CLIGHT; /* ns -> m */
        }
    }
    fclose(fp);

    return 1;
}
/* satellite system to index */
extern int sys2ix(int sys)
{
    switch (sys) {
        case SYS_GPS: return 0;
        case SYS_SBS: return 0;
        case SYS_GLO: return 1;
        case SYS_GAL: return 2;
        case SYS_CMP: return 3;
        case SYS_QZS: return 4;
        case SYS_IRN: return 5;
    }
    return 0;
}
/* translate satellite system and code to code bias table index ----------------
*       -1 = code not supported
*        0 = reference code (0 bias)
*        1-3 = table index for code
* ----------------------------------------------------------------------------*/
extern int code2bias_ix(int sys, int code) {
    int sys_ix;

    sys_ix=sys2ix(sys);
    if (sys_ix<MAX_BIAS_SYS)
        return code_Obias_ix[sys_ix][code];
    else
        return 0;
}
/* read DCB parameters from BIA or BSX file ------------------------------------
*    - supports satellite code biases only
*-----------------------------------------------------------------------------*/
static int readbiaf(const prcopt_t *prcopt, const char *file, nav_t *nav)
{
    FILE *fp;
    double cbias;
    char buff[256],bias[6]="",svn[6]="",prn[6]="",obs1[6]="",obs2[6],*Osbias;
    int i,id,sat,freq,code1,code2,bias_ix1,bias_ix2,sys,type;
    double code_Dbias[MAXSTA][MAXCODE],gamma1,alpha,beta;

    trace(3,"readbiaf: file=%s\n",file);

    if (!(fp=fopen(file,"r"))) {
        trace(1,"MEGX DCB file open error: %s!\n",file);
        return 0;
    }
    if (!(Osbias=(char *)malloc(sizeof(char)*256))) {
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        if (sscanf(buff,"%4s %5s %4s %4s %4s",bias,svn,prn,obs1,obs2)<5) continue;
        if (obs1[0]!='C') continue;  /* skip phase biases for now */
        if ((cbias=str2num(buff,70,21))==0.0) continue;
        sat=satid2no(prn);
        sys=satsys(sat,NULL);
        
        if (!(code1=obs2code(&obs1[1]))) continue; /* skip if code not valid */
        bias_ix1=code2bias_ix(sys,code1);
        if (strcmp(bias,"OSB")==0) {
            nav->obias_flag=OPT_OSB;
            nav->obias[sat-1][bias_ix1]=cbias*1E-9*CLIGHT; /* ns -> m */
        }
        else if (strcmp(bias,"DSB")==0) {
            /* differential signal bias */
            nav->obias_flag=OPT_DCB;
            if (!(code2=obs2code(&obs2[1]))) continue; /* skip if code not valid */
            bias_ix2=code2bias_ix(sys,code2);
            /* observed signal bias */
            strncpy(Osbias,buff+24,10);
            if (sys==SYS_GPS) {
                if (strstr(Osbias, "C1C  C1W")) type=0;
                else if (strstr(Osbias, "C2C  C2W")) type=1;
                else if (strstr(Osbias, "C2W  C2S")) type=2;
                else if (strstr(Osbias, "C2W  C2L")) type=3;
                else if (strstr(Osbias, "C2W  C2X")) type=4;
                else if (strstr(Osbias, "C1C  C2W")) type=5;
                else if (strstr(Osbias, "C1C  C5Q")) type=6;
                else if (strstr(Osbias, "C1C  C5X")) type=7;
                else if (strstr(Osbias, "C1W  C2W")) type=8;
            }
            else if (sys==SYS_GLO) {
                if (strstr(Osbias, "C1C  C1P")) type=0;
                else if (strstr(Osbias, "C2C  C2P")) type=1;
                else if (strstr(Osbias, "C1C  C2C")) type=2;
                else if (strstr(Osbias, "C1C  C2P")) type=3;
                else if (strstr(Osbias, "C1P  C2P")) type=3;
            }
            else if (sys==SYS_GAL) {
                if (strstr(Osbias, "C1C  C5Q")) type=0;
                else if (strstr(Osbias, "C1C  C6C")) type=1;
                else if (strstr(Osbias, "C1C  C7Q")) type=2;
                else if (strstr(Osbias, "C1C  C8Q")) type=3;
                else if (strstr(Osbias, "C1X  C5X")) type=4;
                else if (strstr(Osbias, "C1X  C7X")) type=5;
                else if (strstr(Osbias, "C1X  C8X")) type=6;
            }
            else if (sys==SYS_CMP) {
                if (strstr(Osbias, "C2I  C7I")) type=0;
                else if (strstr(Osbias, "C2I  C6I")) type=1;
                else if (strstr(Osbias, "C1X  C5X")) type=2;
                else if (strstr(Osbias, "C1P  C5P")) type=3;
                else if (strstr(Osbias, "C1D  C5D")) type=4;
                else if (strstr(Osbias, "C1X  C6I")) type=5;
                else if (strstr(Osbias, "C1P  C6I")) type=6;
                else if (strstr(Osbias, "C1D  C6I")) type=7;
                else if (strstr(Osbias, "C1X  C7Z")) type=8;
                else if (strstr(Osbias, "C1X  C8X")) type=9;
                
            }
            else if (sys==SYS_QZS) {
                if (strstr(Osbias, "C1C  C1X")) type=0;
                else if (strstr(Osbias, "C1C  C2L")) type=1;
                else if (strstr(Osbias, "C1C  C5X")) type=2;
                else if (strstr(Osbias, "C1C  C5Q")) type=3;
                else if (strstr(Osbias, "C1X  C2X")) type=4;
                else if (strstr(Osbias, "C1X  C5X")) type=5;
            }
            code_Dbias[sat][type]=cbias*1E-9*CLIGHT;
        }
    }
    /* NOTE: convert DCB to OSB correction */
    if (OPT_DCB==nav->obias_flag) {
        for (i=1;i<=MAXSTA;i++) {
            sys=satsys(i,&id);
            /*ref:C1W-C2W*/
            if (sys==SYS_GPS) {
                gamma1=SQR(FREQL1/FREQL2);
                alpha=gamma1/(gamma1-1);beta=-1/(gamma1-1);
                nav->obias[i-1][0]=(beta*code_Dbias[i][8]+code_Dbias[i][0]);    /* C1C beta*(C1W-C2W)+C1C-C1W */
                nav->obias[i-1][1]=(beta*code_Dbias[i][8]);                     /* C1W beta*(C1W-C2W) */
                nav->obias[i-1][2]=(-alpha*code_Dbias[i][8]);                   /* C2W -alpha*(C1W-C2W) */
                nav->obias[i-1][3]=(-alpha*code_Dbias[i][8]-code_Dbias[i][3]);  /* C2L -alpha*(C1W-C2W)+C2L-C2W */
                nav->obias[i-1][4]=(-alpha*code_Dbias[i][8]-code_Dbias[i][2]);  /* C2S -alpha*(C1W-C2W)+C2S-C2W */
                nav->obias[i-1][5]=(-alpha*code_Dbias[i][8]-code_Dbias[i][4]);  /* C2X -alpha*(C1W-C2W)+C2X-C2W */
                nav->obias[i-1][6]=(-alpha*code_Dbias[i][8]+code_Dbias[i][1]);  /* C2C -alpha*(C1W-C2W)+C2C-C2W */
                nav->obias[i-1][7]=-(alpha*(code_Dbias[i][6]-code_Dbias[i][0])+beta*(code_Dbias[i][6]-code_Dbias[i][5])); /* C5Q -(alpha*(C1W-C5Q)+beta*(C2W-C5Q)) */
                nav->obias[i-1][8]=-(alpha*(code_Dbias[i][7]-code_Dbias[i][0])+beta*(code_Dbias[i][7]-code_Dbias[i][5])); /* C5X -(alpha*(C1W-C5X)+beta*(C2W-C5X)) */
            }
            /*ref:C1X-C5X*/
            if (sys==SYS_GAL) {
                /*NOTE:C1C/C5Q/C6C/C7Q/C8Q is inaccurate correction*/
                gamma1=SQR(FREQL1/FREQL5);
                alpha=gamma1/(gamma1-1);beta=-1/(gamma1-1);
                nav->obias[i-1][0]=(beta*code_Dbias[i][4]);                     /* C1C beta*(C1X-C5X)+C1C-C5X  miss C1C-C5X*/
                nav->obias[i-1][1]=(beta*code_Dbias[i][4]);                     /* C1X beta*(C1X-C5X)*/
                nav->obias[i-1][2]=(-alpha*code_Dbias[i][4]);                   /* C5Q -alpha*(C1X-C5X)+C5Q-C5X  miss C5Q-C5X*/
                nav->obias[i-1][3]=(-alpha*code_Dbias[i][4]);                   /* C5X -alpha*(C1X-C5X)*/
                nav->obias[i-1][4]=-(alpha*code_Dbias[i][1]+beta*(code_Dbias[i][1]-code_Dbias[i][4])); /* C6C -(alpha*(C1X-C6C)+beta*(C5X-C6C))*/
                nav->obias[i-1][5]=-(alpha*code_Dbias[i][2]+beta*(code_Dbias[i][2]-code_Dbias[i][4])); /* C7Q -(alpha*(C1X-C7Q)+beta*(C5X-C7Q))*/
                nav->obias[i-1][6]=-(alpha*code_Dbias[i][5]+beta*(code_Dbias[i][5]-code_Dbias[i][4])); /* C7X -(alpha*(C1X-C7X)+beta*(C5X-C7X))*/
                nav->obias[i-1][7]=-(alpha*code_Dbias[i][3]+beta*(code_Dbias[i][3]-code_Dbias[i][4])); /* C8Q -(alpha*(C1X-C8Q)+beta*(C5X-C8Q))*/
                nav->obias[i-1][8]=-(alpha*code_Dbias[i][6]+beta*(code_Dbias[i][6]-code_Dbias[i][4])); /* C8X -(alpha*(C1X-C8X)+beta*(C5X-C8X))*/
            }
            
            if (sys==SYS_CMP) {
                /*ref:B3I,C6I*/
                /* Broadcast ephemeris DCB correction */
                if (id<19) {
                    nav->bds_tgd[id-1][0]=code_Dbias[i][1]; /* C2I C2I-C6I */
                    nav->bds_tgd[id-1][1]=0.0;                 /* C6I */
                    nav->bds_tgd[id-1][2]=(code_Dbias[i][1]-code_Dbias[i][0]); /* C7I C7I-C6I=(C2I-C6I)-(C2I-C7I) */ 
                }
                else {
                    nav->bds_tgd[id-1][0]=code_Dbias[i][1]; /* C2I C2I-C6I */
                    nav->bds_tgd[id-1][1]=0.0;              /* C6I */
                    nav->bds_tgd[id-1][2]=0.0;              /* no C7I  */ 
                    nav->bds_tgd[id-1][3]=code_Dbias[i][6]; /* C1P C1P-C6I */ 
                    nav->bds_tgd[id-1][4]=code_Dbias[i][5]; /* C1X C7I-C6I */ 
                    nav->bds_tgd[id-1][5]=(code_Dbias[i][6]-code_Dbias[i][3]); /* C5P C5P-C6I=(C1P-C6I)-(C1P-C5P) */ 
                    nav->bds_tgd[id-1][6]=(code_Dbias[i][5]-code_Dbias[i][2]); /* C5X C5X-C6I=(C1X-C6I)-(C1X-C5X) */ 
                    nav->bds_tgd[id-1][7]=(code_Dbias[i][5]-code_Dbias[i][8]); /* C7Z C7Z-C6I=(C1X-C6I)-(C1X-C7Z) */ 
                    nav->bds_tgd[id-1][8]=(code_Dbias[i][5]-code_Dbias[i][9]); /* C8X C8X-C6I=(C1X-C6I)-(C1X-C8X) */
                    nav->bds_tgd[id-1][9]=code_Dbias[i][7]; /* C1D C1D-C6I */ 
                    nav->bds_tgd[id-1][10]=(code_Dbias[i][7]-code_Dbias[i][4]); /* C5D C5D-C6I=(C1D-C6I)-(C1D-C5D) */ 
                }                    
                /*ref:C2I-C6I*/
                /* Precise ephemeris DCB correction */
                if (EPHOPT_PREC==prcopt->sateph) {
                    gamma1=SQR(FREQ1_CMP/FREQ3_CMP);
                    alpha=gamma1/(gamma1-1);beta=-1/(gamma1-1);                    
                    if (id<19) {
                        nav->obias[i-1][0]=beta*code_Dbias[i][1];       /* C2I beta*(C2I-C6I) */
                        nav->obias[i-1][1]=-alpha*code_Dbias[i][1];     /* C6I -alpha*(C2I-C6I)*/
                        nav->obias[i-1][2]=-(alpha*code_Dbias[i][0]+beta*(code_Dbias[i][0]-code_Dbias[i][1]));     /* C7I -(alpha*(C2I-C7I)+beta*(C6I-C7I)) */ 
                    }
                    else {
                        nav->obias[i-1][0]=beta*code_Dbias[i][1];       /* C2I beta*(C2I-C6I) */
                        nav->obias[i-1][1]=-alpha*code_Dbias[i][1];     /* C6I -alpha*(C2I-C6I)*/
                        nav->obias[i-1][2]=0;                           /* no C7I  */ 
                        nav->obias[i-1][3]=-(alpha*(code_Dbias[i][1]-code_Dbias[i][7])-beta*code_Dbias[i][7]); /* C1D -(alpha*(C2I-C1D)+beta*(C6I-C1D))*/ 
                        nav->obias[i-1][4]=-(alpha*(code_Dbias[i][1]-code_Dbias[i][6])-beta*code_Dbias[i][6]); /* C1P  -(alpha*(C2I-C1P)+beta*(C6I-C1P))*/ 
                        nav->obias[i-1][5]=-(alpha*(code_Dbias[i][1]-code_Dbias[i][5])-beta*code_Dbias[i][5]); /* C1X  -(alpha*(C2I-C1X)+beta*(C6I-C1X))*/ 
                        nav->obias[i-1][6]=-(alpha*(code_Dbias[i][1]-code_Dbias[i][7]+code_Dbias[i][4])+beta*(code_Dbias[i][4]-code_Dbias[i][7])); /* C5D -(alpha*(C2I-C5D)+beta*(C6I-C5D))*/ 
                        nav->obias[i-1][7]=-(alpha*(code_Dbias[i][1]-code_Dbias[i][6]+code_Dbias[i][3])+beta*(code_Dbias[i][3]-code_Dbias[i][6])); /* C5P -(alpha*(C2I-C5P)+beta*(C6I-C5P))*/ 
                        nav->obias[i-1][8]=-(alpha*(code_Dbias[i][1]-code_Dbias[i][5]+code_Dbias[i][2])+beta*(code_Dbias[i][2]-code_Dbias[i][5]));  /* C5X -(alpha*(C2I-C5X)+beta*(C6I-C5X))*/ 
                        nav->obias[i-1][9]=-(alpha*(code_Dbias[i][1]-code_Dbias[i][5]+code_Dbias[i][8])+beta*(code_Dbias[i][8]-code_Dbias[i][5])); /* C7Z -(alpha*(C2I-C7Z)+beta*(C6I-C7Z))*/ 
                        nav->obias[i-1][10]=-(alpha*(code_Dbias[i][1]-code_Dbias[i][5]+code_Dbias[i][9])+beta*(code_Dbias[i][9]-code_Dbias[i][5])); /* C8X -(alpha*(C2I-C8X)+beta*(C6I-C8X))*/
                    }                      
                }                
            }
            /*ref:C1C-C2X?*/
            if (sys==SYS_QZS) {
                gamma1=SQR(FREQL1/FREQL2);
                alpha=gamma1/(gamma1-1);beta=-1/(gamma1-1);
                nav->obias[i-1][0]=(beta*(code_Dbias[i][0]+code_Dbias[i][4]));                  /* C1C beta*(C1C-C2X) */
                nav->obias[i-1][1]=(beta*(code_Dbias[i][0]+code_Dbias[i][4])-code_Dbias[i][0]); /* C1X beta*(C1C-C2X)+C1X-C1C */
                nav->obias[i-1][2]=(-alpha*(code_Dbias[i][0]+code_Dbias[i][4])-(code_Dbias[i][0]-code_Dbias[i][1]+code_Dbias[i][4])); /* C2L -alpha*(C1C-C2X)+C2L-C2X */
                nav->obias[i-1][3]=-(alpha*(code_Dbias[i][0]+code_Dbias[i][4]));                /* C2X -alpha*(C1C-C2X) */
                nav->obias[i-1][4]=-(alpha*code_Dbias[i][3]+beta*(code_Dbias[i][3]-code_Dbias[i][0]-code_Dbias[i][4]));     /* C5Q -(alpha*(C1C-C5Q)beta*(C2X-C5Q)) */
                nav->obias[i-1][5]=-(alpha*(code_Dbias[i][0]+code_Dbias[i][5])+beta*(code_Dbias[i][5]-code_Dbias[i][4]));   /* C5X -(alpha*(C1C-C5X)beta*(C2X-C5X)) */
            }
        }
    }
    fclose(fp);

    return 1;
}
/* store tgd corrections to nav->obias*/
extern int tgdarrge(const prcopt_t *opt, nav_t *nav) 
{   
    int i,j,sys,id;
    double gamma1,alpha,beta,b1,b2;

    /* if has DCB file, using DCB correction */
    /* NOTE:for BDS, DCB/OSB corrections are only valid up to C46, so satellites after C46 use TGD corrections */
    if (nav->obias_flag>0&&!(opt->navsys&SYS_CMP)) {
        return 0;
    }

    /* init translation table from code to table column */
    init_bias_ix(); 

    /* if no DCB and OSB file, using TGD correction*/
    if (0==nav->obias_flag||(opt->navsys&SYS_CMP)) {
        if (0==nav->obias_flag) nav->obias_flag=OPT_TGD;
        for (i=1;i<=MAXSTA;i++) {
            sys=satsys(i,&id);
            /*ref:C1W-C2W*/
            if (sys==SYS_GPS&&OPT_TGD==nav->obias_flag) {
                b1=gettgd(i,nav,0); /*DCB beta*(C1-C2)*/
                gamma1=SQR(FREQL1/FREQL2);
                nav->obias[i-1][0]=b1;                     /* C1C beta*(C1-C2) */
                nav->obias[i-1][1]=b1;                     /* C1W beta*(C1-C2) */
                nav->obias[i-1][2]=gamma1*b1;              /* C2W -alpha*(C1-C2) */
                nav->obias[i-1][3]=gamma1*b1;;             /* C2L -alpha*(C1-C2) */
                nav->obias[i-1][4]=gamma1*b1;;             /* C2S -alpha*(C1-C2) */
                nav->obias[i-1][5]=gamma1*b1;;             /* C2X -alpha*(C1-C2) */
                nav->obias[i-1][8]=gamma1*b1;;             /* C2C -alpha*(C1-C2) */
            }
            /*ref:C1X-C5X*/
            if (sys==SYS_GAL&&OPT_TGD==nav->obias_flag) {
                b1=gettgd(i,nav,0);  /*DCB beta*(C1-C5)*/
                b2=gettgd(i,nav,1);  /*DCB beta*(C1-C7)*/
                /* gamma1=SQR(FREQL1/FREQL5); */
                gamma1=SQR(FREQL1/FREQE5b);
                alpha=gamma1/(gamma1-1);beta=-1/(gamma1-1);
                nav->obias[i-1][0]=b1;                          /* C1C beta*(C1-C5) */
                nav->obias[i-1][1]=b1;                          /* C1X beta*(C1-C5) */
                nav->obias[i-1][2]=gamma1*b1;                   /* C5Q -alpha*(C1-C5) */
                nav->obias[i-1][3]=gamma1*b1;                   /* C5X -alpha*(C1-C5)*/
                nav->obias[i-1][5]=-(alpha*b2+beta*(b2-b1));    /* C7Q -(alpha*(C1-C7)+beta*(C5-C7))*/
                nav->obias[i-1][6]=-(alpha*b2+beta*(b2-b1));    /* C7X -(alpha*(C1-C7)+beta*(C5-C7))*/
            }
            /* if DCB products are not provided, use TGD to correct BDS broadcast ephemeris. */
            if (sys==SYS_CMP&&OPT_DCB!=nav->obias_flag) {
                b1=gettgd(i,nav,0); /*DCB C2I-C6I*/
                b2=gettgd(i,nav,1); /*DCB C7I-C6I*/
                /*ref:B3I,C6I*/
                if (id<19) {
                    nav->bds_tgd[id-1][0]=b1;      /* C2I C2I-C6I */
                    nav->bds_tgd[id-1][1]=0.0;     /* C6I */
                    nav->bds_tgd[id-1][2]=b2;      /* C7I C7I-C6I=(C2I-C6I)-(C2I-C7I) */ 
                }
                else {
                    nav->bds_tgd[id-1][0]=b1;      /* C2I C2I-C6I */
                    nav->bds_tgd[id-1][1]=0.0;     /* C6I */
                    nav->bds_tgd[id-1][2]=0.0;     /* no C7I  */  
                }                         
            }
            /*ref:C1C-C2X?*/
            if (sys==SYS_QZS&&OPT_TGD==nav->obias_flag) {
                b1=gettgd(i,nav,0);     /*DCB C2-C2*/
                gamma1=SQR(FREQL1/FREQL2);
                nav->obias[i-1][0]=b1;                  /* C1C beta*(C1-C2) */
                nav->obias[i-1][1]=b1;                  /* C1X beta*(C1-C2) */
                nav->obias[i-1][2]=gamma1*b1;           /* C2L -alpha*(C1-C2) */
                nav->obias[i-1][3]=gamma1*b1;           /* C2X -alpha*(C1-C2) */
            }
        }
        return 1;
    }
 
}
/* read DCB parameters ---------------------------------------------------------
* read differential code bias (DCB) parameters
* args   : char   *file       I   DCB parameters file (wild-card * expanded)
*          nav_t  *nav        IO  navigation data
*          sta_t  *sta        I   station info data to import receiver DCB
*                                 (NULL: no use)
* return : status (1:ok,0:error)
* notes  : supports DCB, BIA, and BSX file formats
         : currently only support P1-P2, P1-C1 bias in DCB file
         : currently only supports satellite biases in BIA/BSX files
*-----------------------------------------------------------------------------*/
extern int readdcb(const prcopt_t *prcopt, const char *file, nav_t *nav, const sta_t *sta)
{
    int i,j,k,n,dcb_ok=0;
    char *efiles[MAXEXFILE]={0};

    trace(3,"readdcb : file=%s\n",file);

    /* init translation table from code to table column */
    init_bias_ix();  

    /* init code biases (DCB/OSB, BDS broadcast ephemeris DCB correction is handled separately ) */
    for (i=0;i<MAXRCV;i++) for (j=0;j<MAX_CODE_BIAS_FREQS;j++) {
        nav->cbias[i][j]=0.0;
        for (k=0;k<MAX_CODE_BIASES;k++) {
            nav->rbias[i][j][k]=0.0;
        }  
    }
    for (i=0;i<MAXSTA;i++) {
        for (j=0;j<MAXCODE;j++) {
            nav->obias[i][j]=0.0;
        }        
    }
    for (i=0;i<=NSATCMP;i++) {
        for (j=0;j<MAXCODE;j++) {
            nav->bds_tgd[i][j]=0.0;
        }
    }
        
    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return 0;
        }
    }
    n=expath(file,efiles,MAXEXFILE);

    for (i=0;i<n;i++) {
        if (strstr(efiles[i],".BIA")||strstr(efiles[i],".bia")||
            strstr(efiles[i],".BSX")||strstr(efiles[i],".bsx"))
            dcb_ok=readbiaf(prcopt,efiles[i],nav);
        else if (strstr(efiles[i],".DCB")||strstr(efiles[i],".dcb"))
            dcb_ok=readdcbf(efiles[i],nav,sta);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);

    return dcb_ok;
}
/* polynomial interpolation by Neville's algorithm ---------------------------*/
static double interppol(const double *x, double *y, int n)
{
    int i,j;

    for (j=1;j<n;j++) {
        for (i=0;i<n-j;i++) {
            y[i]=(x[i+j]*y[i]-x[i]*y[i+1])/(x[i+j]-x[i]);
        }
    }
    return y[0];
}
/* satellite position by precise ephemeris -----------------------------------*/
static int pephpos(gtime_t time, int sat, const nav_t *nav, double *rs,
                   double *dts, double *vare, double *varc)
{
    double t[NMAX+1],p[3][NMAX+1],c[2],*pos,std=0.0,s[3],sinl,cosl;
    int i,j,k,index;

    trace(3,"pephpos : time=%s sat=%2d\n",time_str(time,3),sat);

    rs[0]=rs[1]=rs[2]=dts[0]=0.0;

    if (nav->ne<NMAX+1||
        timediff(time,nav->peph[0].time)<-MAXDTE||
        timediff(time,nav->peph[nav->ne-1].time)>MAXDTE) {
        trace(2,"no prec ephem %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    /* binary search */
    for (i=0,j=nav->ne-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->peph[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;

    /* polynomial interpolation for orbit */
    i=index-(NMAX+1)/2;
    if (i<0) i=0; else if (i+NMAX>=nav->ne) i=nav->ne-NMAX-1;

    for (j=0;j<=NMAX;j++) {
        t[j]=timediff(nav->peph[i+j].time,time);
        if (norm(nav->peph[i+j].pos[sat-1],3)<=0.0) {
            trace(2,"prec ephem outage %s sat=%2d\n",time_str(time,0),sat);
            return 0;
        }
    }
    for (j=0;j<=NMAX;j++) {
        pos=nav->peph[i+j].pos[sat-1];
        /* correction for earth rotation ver.2.4.0 
        This is just to obtain the satellite position at the time of signal transmission */
        sinl=sin(OMGE*t[j]);
        cosl=cos(OMGE*t[j]);
        p[0][j]=cosl*pos[0]-sinl*pos[1];
        p[1][j]=sinl*pos[0]+cosl*pos[1];
        p[2][j]=pos[2];
    }
    for (i=0;i<3;i++) {
        rs[i]=interppol(t,p[i],NMAX+1);
    }
    if (vare) {
        for (i=0;i<3;i++) s[i]=nav->peph[index].std[sat-1][i];
        std=norm(s,3);

        /* extrapolation error for orbit */
        if      (t[0   ]>0.0) std+=EXTERR_EPH*SQR(t[0   ])/2.0;
        else if (t[NMAX]<0.0) std+=EXTERR_EPH*SQR(t[NMAX])/2.0;
        *vare=SQR(std);
    }
    /* linear interpolation for clock */
    t[0]=timediff(time,nav->peph[index  ].time);
    t[1]=timediff(time,nav->peph[index+1].time);
    c[0]=nav->peph[index  ].pos[sat-1][3];
    c[1]=nav->peph[index+1].pos[sat-1][3];

    if (t[0]<=0.0) {
        if ((dts[0]=c[0])!=0.0) {
            std=nav->peph[index].std[sat-1][3]*CLIGHT-EXTERR_CLK*t[0];
        }
    }
    else if (t[1]>=0.0) {
        if ((dts[0]=c[1])!=0.0) {
            std=nav->peph[index+1].std[sat-1][3]*CLIGHT+EXTERR_CLK*t[1];
        }
    }
    else if (c[0]!=0.0&&c[1]!=0.0) {
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=t[0]<-t[1]?0:1;
        std=nav->peph[index+i].std[sat-1][3]+EXTERR_CLK*fabs(t[i]);
    }
    else {
        dts[0]=0.0;
    }
    if (varc) *varc=SQR(std);
    return 1;
}
/* satellite clock by precise clock ------------------------------------------*/
extern int pephclk(gtime_t time, int sat, const nav_t *nav, double *dts,
                   double *varc)
{
    double t[2],c[2],std;
    int i,j,k,index;

    trace(3,"pephclk : time=%s sat=%2d\n",time_str(time,3),sat);

    if (nav->nc<2||
        timediff(time,nav->pclk[0].time)<-MAXDTE||
        timediff(time,nav->pclk[nav->nc-1].time)>MAXDTE) {
        trace(2,"no prec clock %s sat=%2d\n",time_str(time,0),sat);
        return 1;
    }
    /* binary search */
    for (i=0,j=nav->nc-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->pclk[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;

    /* linear interpolation for clock */
    t[0]=timediff(time,nav->pclk[index  ].time);
    t[1]=timediff(time,nav->pclk[index+1].time);
    c[0]=nav->pclk[index  ].clk[sat-1][0];
    c[1]=nav->pclk[index+1].clk[sat-1][0];

    if (t[0]<=0.0) {
        if ((dts[0]=c[0])==0.0) return 0;
        std=nav->pclk[index].std[sat-1][0]*CLIGHT-EXTERR_CLK*t[0];
    }
    else if (t[1]>=0.0) {
        if ((dts[0]=c[1])==0.0) return 0;
        std=nav->pclk[index+1].std[sat-1][0]*CLIGHT+EXTERR_CLK*t[1];
    }
    else if (c[0]!=0.0&&c[1]!=0.0) {
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=t[0]<-t[1]?0:1;
        std=nav->pclk[index+i].std[sat-1][0]*CLIGHT+EXTERR_CLK*fabs(t[i]);
    }
    else {
        trace(2,"prec clock outage %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    if (varc) *varc=SQR(std);
    return 1;
}
/* satellite antenna phase center offset ---------------------------------------
* compute satellite antenna phase center offset in ecef
* args   : gtime_t time       I   time (gpst)
*          double *rs         I   satellite position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          double *dant       O   satellite antenna phase center offset (ecef)
*                                 {dx,dy,dz} (m) (iono-free LC value)
* return : none
* notes  : iono-free LC frequencies defined as follows:
*            GPS/QZSS : L1-L2
*            GLONASS  : G1-G2
*            Galileo  : E1-E5a
*            BDS      : B1I-B3I
*            NavIC    : L5-S
*-----------------------------------------------------------------------------*/
extern void satantoff(gtime_t time, const double *rs, int fr, int sat, const nav_t *nav,
                      double *dant)
{
    const spcv_t *spcv=nav->spcvs+sat-1;
    double ex[3],ey[3],ez[3],es[3],r[3],rsun[3],gmst,erpv[5]={0},freq[2];
    double C1,C2,dant1,dant2;
    int i,sys,j,k;

    trace(3,"satantoff: time=%s sat=%2d\n",time_str(time,3),sat);

    dant[0]=dant[1]=dant[2]=0.0;

    /* sun position in ecef */
    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,&gmst);

    /* unit vectors of satellite fixed coordinates */
    for (i=0;i<3;i++) r[i]=-rs[i];
    if (!normv3(r,ez)) return;
    for (i=0;i<3;i++) r[i]=rsun[i]-rs[i];
    if (!normv3(r,es)) return;
    cross3(ez,es,r);
    if (!normv3(r,ey)) return;
    cross3(ey,ez,ex);

    if (fr>=0&&norm(spcv->off[fr],3)>0) {
        for (i=0;i<3;i++) dant[i]=spcv->off[fr][0]*ex[i]+spcv->off[fr][1]*ey[i]+spcv->off[fr][2]*ez[i];
    }
    else {
        /* if there is no satellite PCO correction information, use L1/L2 iono-free LC coefficients */
        sys=satsys(sat,NULL);
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2 */
            j=0,k=1;
            freq[0]=FREQL1;
            freq[1]=FREQL2;
        }
        else if (sys==SYS_GLO) { /* G1-G2 */
            j=0;k=1;
            freq[0]=sat2freq(sat,CODE_L1C,nav);
            freq[1]=sat2freq(sat,CODE_L2C,nav);
        }
        else if (sys==SYS_GAL) { /* E1-E5a */
            j=0;k=2;
            freq[0]=FREQL1;
            freq[1]=FREQL5;
            /* freq[1]=FREQE5b; */
        }
        else if (sys==SYS_CMP) { /* B1I-B3I */
            j=0;k=2;
            freq[0]=FREQ1_CMP;
            freq[1]=FREQ3_CMP;
        }
        else if (sys==SYS_QZS) { /* L1-L2 */
            j=0;k=1;
            freq[0]=FREQL1;
            freq[1]=FREQL2;
        }
        else if (sys==SYS_IRN) { /* L1-L2 */
            j=0;k=1;
            freq[0]=FREQL5;
            freq[1]=FREQs;
        }
        else return;

        C1= SQR(freq[0])/(SQR(freq[0])-SQR(freq[1]));
        C2=-SQR(freq[1])/(SQR(freq[0])-SQR(freq[1]));

        /* iono-free LC */
        for (i=0;i<3;i++) {
            dant1=spcv->off[j][0]*ex[i]+spcv->off[j][1]*ey[i]+spcv->off[j][2]*ez[i];
            dant2=spcv->off[k][0]*ex[i]+spcv->off[k][1]*ey[i]+spcv->off[k][2]*ez[i];
            dant[i]=C1*dant1+C2*dant2;
        }
    }
}
/* satellite position/clock by precise ephemeris/clock -------------------------
* compute satellite position/clock with precise ephemeris/clock
* args   : gtime_t time       I   time (gpst)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          int    opt         I   sat position option
*                                 (0: center of mass, 1: antenna phase center)
*          double *rs         O   sat position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          double *dts        O   sat clock {bias,drift} (s|s/s)
*          double *var        IO  sat position and clock error variance (m)
*                                 (NULL: no output)
* return : status (1:ok,0:error or data outage)
* notes  : clock includes relativistic correction but does not contain code bias
*          before calling the function, nav->peph, nav->ne, nav->pclk and
*          nav->nc must be set by calling readsp3(), readrnx() or readrnxt()
*          if precise clocks are not set, clocks in sp3 are used instead
*-----------------------------------------------------------------------------*/
extern int peph2pos(gtime_t time, int sat, const nav_t *nav, int opt,
                    double *rs, double *dts, double *var)
{
    gtime_t time_tt;
    double rss[3],rst[3],dtss[1],dtst[1],dant[3]={0},vare=0.0,varc=0.0,tt=1E-3;
    int i;

    trace(3,"peph2pos: time=%s sat=%2d opt=%d\n",time_str(time,3),sat,opt);

    if (sat<=0||MAXSAT<sat) return 0;

    /* satellite position and clock bias */
    if (!pephpos(time,sat,nav,rss,dtss,&vare,&varc)||
        !pephclk(time,sat,nav,dtss,&varc)) return 0;

    time_tt=timeadd(time,tt);
    if (!pephpos(time_tt,sat,nav,rst,dtst,NULL,NULL)||
        !pephclk(time_tt,sat,nav,dtst,NULL)) return 0;

    /* satellite antenna offset correction */
    /* if (opt) {
        satantoff(time,rss,-1,sat,nav,dant);
    } */
    for (i=0;i<3;i++) {
        rs[i  ]=rss[i]+dant[i];
        rs[i+3]=(rst[i]-rss[i])/tt;
    }
    /* relativistic effect correction */
    if (dtss[0]!=0.0) {
        dts[0]=dtss[0]-2.0*dot3(rs,rs+3)/CLIGHT/CLIGHT;
        dts[1]=(dtst[0]-dtss[0])/tt;
    }
    else { /* no precise clock */
        dts[0]=dts[1]=0.0;
    }
    if (var) *var=vare+varc;

    return 1;
}
