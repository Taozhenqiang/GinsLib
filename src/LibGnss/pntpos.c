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

#define SQR(x)      ((x)*(x))
#define MAX(x,y)    ((x)>=(y)?(x):(y))

#define QZSDT /* enable GPS-QZS time offset estimation */
#ifdef QZSDT
#define NX          (4+5)       /* # of estimated parameters */
#else
#define NX          (4+4)       /* # of estimated parameters */
#endif
#define ERR_ION     5.0         /* ionospheric delay Std (m) */
#define ERR_TROP    3.0         /* tropspheric delay Std (m) */
#define ERR_SAAS    0.3         /* Saastamoinen model error Std (m) */
#define ERR_BRDCI   0.2         /* broadcast ionosphere model error factor */
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */
#define MAX_GDOP   30          /* max gdop for valid solution  */

/* add by tzq */
#define VAR_CLK     SQR(60.0) /* init variance receiver clock (m^2) */
#define VAR_CLKD    SQR(30.0) /* init variance receiver clock drift (m^2/s^2) */

/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NP(opt)     ((opt)->GI_mode==GINS_TC?15:3)
#define NC(opt)     (((opt)->GI_mode==GINS_TC&&(opt)->mode==PMODE_SINGLE)?7:0) /* clock (0-5,GPS,GLONASS,Galileo,BDS,IRNSS,QZSS); clock drift (6:GPS)*/

/* state variable index */
#define IC(s,opt)   (NP(opt)+(s))

/* pseudorange measurement error variance ------------------------------------*/
extern double varerr_spp(const prcopt_t *opt, const ssat_t *ssat, const obsd_t *obs, double el, int sys)
{
    int fr2[2]={0};
    double fact=1.0,varr,snr_rover,maxsnr_rover=0.0;

    fr2[0]=sys2freid(sys,0,opt);
    fr2[1]=sys2freid(sys,1,opt);
    if (ssat) {
        if (IONOOPT_BRDC==opt->ionoopt) maxsnr_rover=ssat->maxsnr_rover[fr2[0]];
        else if (IONOOPT_IFLC==opt->ionoopt) maxsnr_rover=(ssat->maxsnr_rover[fr2[0]]+ssat->maxsnr_rover[fr2[1]])/2.0;
        else maxsnr_rover=ssat->maxsnr_rover[fr2[0]];        
    }

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
    /* var = R^2*(a^2 + (b^2/sin(el) + c^2*(10^(0.1*(snr_max-snr_rover)))) + (d*rcv_std)^2) */
    varr=SQR(opt->err[1])+SQR(opt->err[2])/sin(el);
    if (opt->err[6]>0.0) {  /* if snr term not zero */
        snr_rover=(ssat)?SNR_UNIT*ssat->snr_rover[0]:opt->err[5];
        /* varr+=SQR(opt->err[6])*pow(10,0.1*MAX(opt->err[5]-snr_rover,0)); */
        varr+=SQR(opt->err[6])*pow(10,0.1*MAX(maxsnr_rover-snr_rover,0));
    }

    varr*=SQR(opt->eratio[0]);
    if (opt->err[7]>0.0) {
        varr+=SQR(opt->err[7]*0.01*(1<<(obs->Pstd[0]+5)));  /* 0.01*2^(n+5) m */
    }

    if (opt->ionoopt==IONOOPT_IFLC) varr*=SQR(3.0); /* iono-free */
    return SQR(fact)*varr;
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
static int snrmask(const obsd_t *obs, const double *azel, const prcopt_t *opt)
{
    int sys,fr;
    char id[4];
    fr=sys2freid(sys,0,opt);

    if (testsnr(0,0,azel[1],obs->SNR[fr]*SNR_UNIT,&opt->snrmask)) {
        satno2id(obs->sat,id);
        trace(7,"SNR check failed: %4s, el=%4.1f, SNR=%5.1f\n",id,azel[1]*R2D,obs->SNR[fr]*SNR_UNIT);
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
static double prange(const obsd_t *obs, const nav_t *nav, const prcopt_t *opt,
                     double *var, double *dcb)
{
    double P1,P2,gamma,b1=0.0,freq1=0.0,freq2=0.0;
    int sat,sys,f2,bias_ix[2],flag=0,fr2[2];

    sat=obs->sat;
    sys=satsys(sat,NULL);
    fr2[0]=sys2freid(sys,0,opt);
    fr2[1]=sys2freid(sys,1,opt);
    P1=obs->P[fr2[0]];
    P2=obs->P[fr2[1]];
    *var=0.0;*dcb=0.0;
    
    if (P1==0.0||(opt->ionoopt==IONOOPT_IFLC&&P2==0.0)) return 0.0;
    bias_ix[0]=code2bias_ix(sys,obs->code[fr2[0]]);  /* L1 code bias */
    bias_ix[1]=code2bias_ix(sys,obs->code[fr2[1]]);

    /* obias_flag 1:DCB product, 2:OSB product */
    if (OPT_OSB==nav->obias_flag&&opt->sateph==EPHOPT_BRDC&&sys==SYS_CMP) { 
        flag=0;
    }
    else if (nav->obias_flag>0) {
        if (bias_ix[0]>=0) {
            P1-=nav->obias[sat-1][bias_ix[0]];
            *dcb=nav->obias[sat-1][bias_ix[0]];
        }
        if (bias_ix[1]>=0) P2-=nav->obias[sat-1][bias_ix[1]];
        flag=1;
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
            if      (!flag&&obs->code[fr2[0]]==CODE_L2I) {P1-=gettgd(sat,nav,0);} /* TGD_B1I */
            else if (!flag&&obs->code[fr2[1]]==CODE_L7I) {P2-=gettgd(sat,nav,1);} /* TGD_B2I */
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
            if      (!flag&&obs->code[fr2[0]]==CODE_L2I) b1=gettgd(sat,nav,0); /* TGD_B1I */
            else if (!flag&&obs->code[fr2[0]]==CODE_L7I) b1=gettgd(sat,nav,1); /* TGD_B2I */
            *dcb=b1;
            return P1-b1;
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
extern int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var)
{
    int err=0;

    trace(9,"ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
          time_str(time,3),ionoopt,sat,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* SBAS ionosphere model */
    if (ionoopt==IONOOPT_SBAS) {
        if (sbsioncorr(time,nav,pos,azel,ion,var)) return 1;
        err=1;
    }
    /* IONEX TEC model */
    if (ionoopt==IONOOPT_TEC) {
        if (iontec(time,nav,pos,azel,1,ion,var)) return 1;
        err=1;
    }
    /* QZSS broadcast ionosphere model */
    if (ionoopt==IONOOPT_QZS&&norm(nav->ion_qzs,8)>0.0) {
        *ion=ionmodel(time,nav->ion_qzs,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* GPS broadcast ionosphere model */
    if (ionoopt==IONOOPT_BRDC||ionoopt==IONOOPT_EST||err==1) {
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
static int rescode(int iter, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   ssat_t *ssat, double *v, double *H, double *var,
                   double *azel, int *vsat, double *resp, int *ns, int *sati, int *vi)
{
    gtime_t time;
    double r,freq,dion=0.0,dtrp=0.0,vmeas,vion=0.0,vtrp=0.0,rr[3],pos[3],e[3],P,dtr=0.0,dcb=0.0;
    int i,j,nv=0,sat,sys,mask[NX-3]={0},fr=0;

    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[3];
    
    ecef2pos(rr,pos);
    trace(8,"rescode: rr=%.3f %.3f %.3f\n",rr[0], rr[1], rr[2]);

    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        time=obs[i].time; sat=obs[i].sat;
        fr=sys2freid(sys,0,opt);
        azel[i*2]=azel[1+i*2]=resp[i]=0.0;

        if (!iter) vsat[i]=0;

        /* reset spp valid satellite flags */
        if (ssat) ssat[sat-1].vs=0;

        /* exclude satellites with large residuals */
        if (iter&&!vsat[i]) continue;
        if (!(sys=satsys(sat,NULL))) continue;
        
        /* reject duplicated observation data */
        if (i<n-1&&i<MAXOBS-1&&sat==obs[i+1].sat) {
            trace(7,"duplicated obs data %s sat=%d\n",time_str(time,3),sat);
            i++; continue;
        }
        
        /* excluded satellite */
        if (satexclude(sat,vare[i],svh[i],opt)) continue;
        
        /* geometric distance and elevation mask*/
        if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
        if (satazel(pos,e,azel+i*2)<opt->elmin) continue;
        
        if (iter>0) {        

            /* test SNR mask */
            if (!snrmask(obs+i,azel+i*2,opt)) continue;

            /* ionospheric correction */
            if (!ionocorr(time,nav,sat,pos,azel+i*2,opt->ionoopt,&dion,&vion)) {
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
        if ((P=prange(obs+i,nav,opt,&vmeas,&dcb))==0.0) continue;
        
        /* pseudorange residual */
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
        if (ssat) trace(8,"sat=%d: v=%.3f P=%.3f r=%.3f dts=%.6f dion=%.3f dtrp=%.3f\n",ssat[sat-1].id,v[nv],P,r,dts[i*2],dion,dtrp);
        
        /* design matrix */
        for (j=0;j<NX;j++) {
            H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0);
        }                

        /* time system offset and receiver bias correction */
        if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; mask[1]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; mask[2]=1;}
        else if (sys==SYS_CMP) {v[nv]-=x[6]; H[6+nv*NX]=1.0; mask[3]=1;}
        else if (sys==SYS_IRN) {v[nv]-=x[7]; H[7+nv*NX]=1.0; mask[4]=1;}
#ifdef QZSDT
        else if (sys==SYS_QZS) {v[nv]-=x[8]; H[8+nv*NX]=1.0; mask[5]=1;}
#endif
        else mask[0]=1;
        
        vsat[i]=1; resp[i]=v[nv]; (*ns)++; sati[nv]=sat; vi[nv]=i;

        /* update solution status */
        if (ssat) update_ssat(ssat+(sat-1),opt,1,sat,fr,rs+i*6,rr,azel+i*2,resp[i],dtr,dts[i*2],dtrp,dion,0.0,NULL,NULL,dcb);

        /* variance of pseudorange error */
        var[nv]=vare[i]+vmeas+vion+vtrp;
        if (ssat) {
            var[nv++]+=varerr_spp(opt,ssat+(sat-1),&obs[i],azel[1+i*2],sys); 
        }           
        else {
            var[nv++]+=varerr_spp(opt,NULL,&obs[i],azel[1+i*2],sys);
        }         
        trace(8,"sat=%2d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n",obs[i].sat,
              azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i],sqrt(var[nv-1]));
    }
    /* constraint to avoid rank-deficient */
    for (i=0;i<NX-3;i++) {
        if (mask[i]) continue;
        v[nv]=0.0;
        for (j=0;j<NX;j++) H[j+nv*NX]=j==i+3?1.0:0.0;
        var[nv++]=0.01;
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
    int i,j,nx=rtk->nx,nv=0,sat,sys,fr=0,ic;

    for (i=0;i<3;i++) rr[i]=x[i];
    ic=IC(0,opt);     dtr=rtk->x[ic];

    xyz2enu(ins->pos,Cne); DCMT(Cne,Cen);
    Mat3mulv(1.0,ins->Cnb,ins->lever,lever_n);
    
    ecef2pos(rr,pos);
    trace(8,"rescode: rr=%.3f %.3f %.3f\n",rr[0], rr[1], rr[2]);
    
    for (i=*ns=0;i<n&&i<MAXOBS;i++) { 
        time=obs[i].time; sat=obs[i].sat;
        fr=sys2freid(sys,0,opt);
        azel[i*2]=azel[1+i*2]=resp[i]=0.0;

        /* reset spp valid satellite flags */
        if (ssat) ssat[sat-1].vs=0;

        /* exclude satellites with large residuals */
        if (!vsat[i]) continue;
        if (!(sys=satsys(sat,NULL))) continue;
        
        /* reject duplicated observation data */
        if (i<n-1&&i<MAXOBS-1&&sat==obs[i+1].sat) {
            trace(7,"duplicated obs data %s sat=%d\n",time_str(time,3),sat);
            i++;
            continue;
        }

        /* excluded satellite */
        if (satexclude(sat,vare[i],svh[i],opt)) continue;
        
        /* geometric distance and elevation mask*/
        if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
        if (satazel(pos,e,azel+i*2)<opt->elmin) continue;

        /* test SNR mask */
        if (!snrmask(obs+i,azel+i*2,opt)) continue;

        /* ionospheric correction */
        if (!ionocorr(time,nav,sat,pos,azel+i*2,opt->ionoopt,&dion,&vion)) {
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
        if ((P=prange(obs+i,nav,opt,&vmeas,&dcb))==0.0) continue;
        
        /* pseudorange residual */
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
        if (ssat) trace(8,"sat=%d: v=%.3f P=%.3f r=%.3f dts=%.6f dion=%.3f dtrp=%.3f\n",
            ssat[i].id,v[nv],P,r,dts[i*2],dion,dtrp);
        
        /* design matrix H */
        vmulMat3(1.0,e,Cen,Hpp);
        vmvskew(1.0,Hpp,lever_n,Hpa);

        for (j=0;j<nx;j++) {
            H[j+nv*nx]=(j==IC(0,opt)?1.0:0.0);
        }    
        
        for (j=0;j<3;j++)   H[j+nv*nx]    =Hpa[j];
        for (j=0;j+6<9;j++) H[(j+6)+nv*nx]=Hpp[j];

        /* time system offset and receiver bias correction */
        if      (sys==SYS_GLO) {ic=IC(1,opt); v[nv]-=rtk->x[ic]; H[ic+nv*nx]=1.0; }
        else if (sys==SYS_GAL) {ic=IC(2,opt); v[nv]-=rtk->x[ic]; H[ic+nv*nx]=1.0; }
        else if (sys==SYS_CMP) {ic=IC(3,opt); v[nv]-=rtk->x[ic]; H[ic+nv*nx]=1.0; }
        else if (sys==SYS_IRN) {ic=IC(4,opt); v[nv]-=rtk->x[ic]; H[ic+nv*nx]=1.0; }
#ifdef QZSDT
        else if (sys==SYS_QZS) {ic=IC(5,opt); v[nv]-=rtk->x[ic]; H[ic+nv*nx]=1.0; }
#endif
        
        vsat[i]=1; resp[i]=v[nv]; (*ns)++; sati[nv]=sat;

        /* update solution status */
        if (ssat) update_ssat(ssat+(sat-1),opt,1,sat,fr,rs+i*6,rr,azel+i*2,resp[i],dtr,dts[i*2],dtrp,dion,0.0,NULL,NULL,dcb);
        
        /* variance of pseudorange error */
        var[nv]=vare[i]+vmeas+vion+vtrp;
        if (ssat) {
            var[nv++]+=varerr_spp(opt,&ssat[i],&obs[i],azel[1+i*2],sys);
        }    
        else {
            var[nv++]+=varerr_spp(opt,NULL,&obs[i],azel[1+i*2],sys);
        }         
        trace(8,"sat=%2d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n",obs[i].sat,
              azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i],sqrt(var[nv-1]));
    }

    return nv;
}
/* outlier rejection for spp ---------------------------------------------------------*/
extern int outrej_spp(int nv, int nx, double thres, double *v, double *H, double *var,
                      const ssat_t *ssat, const int *sati, const int *vi, int *vsat, int it) 
{
    double mean,std,dv;
    int j,k,m;
    double *v_,*H_,*var_;

    v_=mat(nv,1); H_=mat(nv,nx); var_=mat(nv,1);

    mean=calexp(v,nv);
    std=calstd(v,nv);
    m=0;
    for (j=0;j<nv;j++) {  

        /* remove auxiliary quantities that prevent least squares rank deficiency */
        if (fabs(v[j])<1E-4) dv=0.0;
        /* standardized residuals */
        else dv=fabs(v[j]-mean)/std;     

        /* threshold for outlier detection */     
        if (dv<thres) {
            v_[m]=v[j];
            for (k=0;k<nx;k++) {
                H_[k+m*nx]=H[k+j*nx];
            }  
            var_[m++]=var[j];
        }
        else {
            if (vsat&&vi) vsat[vi[j]]=0;
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

    free(v_); free(H_); free(var_);

    return nv;
}


/* validate solution ---------------------------------------------------------*/
static int valsol(const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, double *P, int nv, int nx)
{
    double azels[MAXOBS*2],dop[4],vv,*vP;
    int i,ns;
    
    trace(3,"valsol  : n=%d nv=%d\n",n,nv);
    
    vP=mat(1,nv);

    /* chi-square validation of residuals */
    matmul("TN",1,nv,nv,v,P,vP,1.0,0.0);
    matmul("NN",1,nv,1,vP,v,&vv,1.0,0.0);
    /* vv=dot(v,v,nv); */
    if (nv>nx&&vv>chisqr[nv-nx-1]) {
        trace(7,"spp error: large chi-square error ns=%d vv=%.1f threshold=%.1f\n",nv,vv,chisqr[nv-nx-1]);
        free(vP);
        return 0; /* threshold too strict for all use cases, report error but continue on */
    }

    /* large GDOP check */
    for (i=ns=0;i<n;i++) {
        if (!vsat[i]) continue;
        azels[  ns*2]=azel[  i*2];
        azels[1+ns*2]=azel[1+i*2];
        ns++;
    }
    dops(ns,azels,opt->elmin,dop);
    if (dop[0]<=0.0||dop[0]>MAX_GDOP) {
        trace(7,"gdop error nv=%d gdop=%.1f\n",nv,dop[0]);
        free(vP);
        return 0;
    }
    
    free(vP);
    return 1;
}

/* time update of position */
static void udpos_spp(rtk_t *rtk)
{
    int i;
    double p_ins[6],Cne[9],Cen[9];
    ins_t *ins=&rtk->ins;

    /* convert INS solutions to GNSS center */
    ins2gnss(rtk,p_ins,6);
    pos2ecef(p_ins,rtk->ru);

    xyz2enu(p_ins,Cne);
    DCMT(Cne,Cen);
    Mat3mulv(1.0,Cen,p_ins+3,rtk->ru+3);

    /* reset ins related state */
    for(i=0;i<ins->nx;i++) rtk->x[i]=0.0;
}

/* time update of clock*/
static void udclk_spp(rtk_t *rtk)
{
    prcopt_t *opt=&rtk->opt;
    double dtr;
    int i,ic,sys=rtk->opt.navsys;
    trace(3,"udclk_spp:\n");

    /* initialize GPS clock (white noise) */
	ic=IC(0,opt);
    dtr=rtk->sol.dtr[0];
    if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
    initx(rtk,CLIGHT*dtr,VAR_CLK,ic);


    /* initialize GPS drift (random walk) */
    ic=IC(6,opt);
    /* if (rtk->x[ic]==0.0) { */
        dtr=rtk->sol.dtr[6];
        /* dtr=1.0e-10; */        /* initial clock drift */
        initx(rtk,CLIGHT*dtr,VAR_CLKD,ic);
    /* }  
    else {
        rtk->P[ic+ic*rtk->nx]+=SQR(1e-3)*fabs(rtk->tt);
    } */    

    /* multi system clock error initialization */
    for (i=1;i<NSYS;i++) {
        if (!(sys&SYS_GLO)&&i==1) continue;
        if (!(sys&SYS_GAL)&&i==2) continue;
        if (!(sys&SYS_CMP)&&i==3) continue;
        if (!(sys&SYS_IRN)&&i==4) continue;
        if (!(sys&SYS_QZS)&&i==5) continue;

        dtr=rtk->sol.dtr[i];
        ic=IC(i,opt);

        if (opt->sysisb==GNSISB_CT) {
            // constant
            if (rtk->x[ic]==0.0) {
                if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
                initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
            }
        }
        else if (opt->sysisb==GNSISB_RW) {
            // random walk process
            if (rtk->x[ic]==0.0) {
                if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
                initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
            }  
            else {
                rtk->P[ic+ic*rtk->nx]+=SQR(1e-3)*fabs(rtk->tt);
            }
        }
        else if (opt->sysisb==GNSISB_WN) {
            //white noise process
            if (fabs(dtr)<1.0e-16) dtr=1.0e-16;
            initx(rtk,CLIGHT*dtr,VAR_CLK,ic);
        }
    }
}

/* time update of states*/
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
    int nx=rtk->nx;
    
    /* if GNSS/INS integration solution is available, reset GNSS outage count to 0 */
    if (rtk->outage<=MAX_OUTIME) rtk->outage=0;
    sol->ns=n;
    sol->stat=stat;
    
    /* update ins state */
    update_instat(ins,rtk->P,sol,nx);

    /* store clock and isb */
    rtk->sol.dtr[0]=rtk->x[IC(0,opt)]/CLIGHT; /* GPS clock (s) */
    rtk->sol.dtr[1]=rtk->x[IC(1,opt)]/CLIGHT; /* GLO-GPS */
    rtk->sol.dtr[2]=rtk->x[IC(2,opt)]/CLIGHT; /* GAL-GPS */
    rtk->sol.dtr[3]=rtk->x[IC(3,opt)]/CLIGHT; /* BDS-GPS */
    rtk->sol.dtr[4]=rtk->x[IC(4,opt)]/CLIGHT; /* IRNSS-GPS */
    rtk->sol.dtr[5]=rtk->x[IC(5,opt)]/CLIGHT; /* QZSS-GPS */
    rtk->sol.dtr[6]=rtk->x[IC(6,opt)]/CLIGHT; /* GPS drift */
}

/* range rate residuals ------------------------------------------------------*/
static int resdop_filter(rtk_t *rtk, const obsd_t *obs, int n, const double *rs, const double *dts,
                  const nav_t *nav, const double *rr, const double *x,
                  const double *azel, const int *vsat, double *v,
                  double *H, double *var, int flag)
{
    ins_t *ins=&rtk->ins;
    double freq,rate,pos[3],a[3],e[3],vs[3],cosel,factor;
    double Hva[3],Hvv[3],Cne[9],Cen[9],Cbn[9],wbie[3],wbeb[3],temp1[3],temp2[3],v_temp[9];
    int i,j,nv=0,sys,fr,nx=(flag&&(GINS_TC==rtk->opt.GI_mode))?rtk->nx:4;
    
    trace(3,"resdop  : n=%d\n",n);
    
    ecef2pos(rr,pos); xyz2enu(pos,Cne);
    DCMT(Cne,Cen);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        sys=satsys(obs[i].sat,NULL);
        fr=sys2freid(sys,0,&rtk->opt);
        freq=sat2freq(obs[i].sat,obs[i].code[fr],nav);
        
        if (obs[i].D[fr]==0.0||freq==0.0||!vsat[i]||norm(rs+3+i*6,3)<=0.0) {
            continue;
        }
        /* LOS (line-of-sight) vector in ECEF */
        cosel=cos(azel[1+i*2]);
        a[0]=sin(azel[i*2])*cosel;
        a[1]=cos(azel[i*2])*cosel;
        a[2]=sin(azel[1+i*2]);
        matmul("TN",3,3,1,Cne,a,e,1.0,0.0);
        
        /* satellite velocity relative to receiver in ECEF */
        for (j=0;j<3;j++) {
            vs[j]=rs[j+3+i*6]-x[j];
        }
        /* range rate with earth rotation correction */
        rate=dot3(vs,e)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                     rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);        
        
        /* range rate residual (m/s) */
        factor=-1;
        if (fabs(-obs[i].D[fr]*CLIGHT/freq-(rate+x[3]-CLIGHT*dts[1+i*2]))>fabs(-obs[i].D[fr]*CLIGHT/freq)) factor=1;

        v[nv]=(factor*obs[i].D[fr]*CLIGHT/freq-(rate+x[3]-CLIGHT*dts[1+i*2]));
        
        /* design matrix */
        if (flag&&(GINS_TC==rtk->opt.GI_mode)){
            vmulMat3(1.0,e,Cen,Hvv);

            DCMT(ins->Cnb,Cbn);
            Mat3mulv(1.0,Cbn,ins->eth.wnie,wbie);
            Mat3add2(ins->wbib,1.0,wbie,-1.0,wbeb);
            vskewmv(1.0,wbeb,ins->lever,temp1);
            Mat3mulv(1.0,ins->Cnb,temp1,temp2);
            
            vskew(1.0,temp2,v_temp);
            vmulMat3(1.0,Hvv,v_temp,Hva);

            for (j=0;j<nx;j++)  H[j+nv*nx]=(j==IC(6,&rtk->opt)?1.0:0.0);
            for (j=0;j<3;j++)   H[j+nv*nx]=Hva[j];
            for (j=0;j+3<6;j++) H[(j+3)+nv*nx]=Hvv[j];            
        }
        else {
            for (j=0;j<nx;j++)  H[j+nv*nx]=((j<3)?-e[j]:1.0);
        }

        /* TODO: stochastic model of Doppler observations */
        var[nv++]=3*varerr_spp(&rtk->opt,NULL,&obs[i],azel[1+i*2],sys);
    }
    return nv;
}

/* estimate receiver velocity ------------------------------------------------*/
extern int estvel(rtk_t *rtk, const obsd_t *obs, int n, const double *rs, const double *dts,
                   const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                   const double *azel, const int *vsat)
{
    double x[4]={0},dx[4],Q[16],*v,*H,*var;
    double sig,err=opt->err[4]; /* Doppler error (Hz) */
    int i,j,k,nv,stat;
    
    v=mat(n,1); H=mat(4,n); var=mat(n,1);
    
    for (i=0;i<MAXITR;i++) {
        
        /* range rate residuals (m/s) */
        if ((nv=resdop_filter(rtk,obs,n,rs,dts,nav,sol->rr,x,azel,vsat,v,H,var,0))<4) {
            stat=0;
            break;
        }

        /* weight by variance (lsq uses sqrt of weight) */
        for (j=0;j<nv;j++) {
            sig=sqrt(var[j]);
            v[j]/=sig;
            for (k=0;k<4;k++) H[k+j*4]/=sig;
        }

        /* least square estimation */
        if (lsq(H,v,4,nv,dx,Q)) break;
        
        for (j=0;j<4;j++) x[j]+=dx[j];
        
        if (norm(dx,4)<1E-6) {
            trace(3,"estvel : vx=%.3f vy=%.3f vz=%.3f, n=%d\n",x[0],x[1],x[2],n);
            matcpy(sol->rr+3,x,3,1);
            sol->dtr[6]=x[3]/CLIGHT;        /* clock drift */
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
    free(v); free(H); free(var);
    return stat;
}

/* estimate receiver position ------------------------------------------------*/
extern int estpos(rtk_t *rtk, const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, ssat_t *ssat, sol_t *sol, double *azel,
                  int *vsat, double *resp)
{
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
    double *P,*R,thres=2.0,zupt_time;
    double *xp,*Pp,vc[4];
    int i,j,k,it,m,info,stat=SOLQ_NONE,mode,nv=0,nv_dop=0,nv_cons=0,ns,*sati,*vi;
    int tc_flag=0; /* spp/ins tc flag */
    
    trace(8,"estpos  : n=%d\n",n);
    
    v=mat(n+5,1); H=mat(n+5,NX); var=mat(n+5,1); P=mat(n+5,n+5);
    sati=imat(n+5,1); vi=imat(n+5,1);
    
    for (i=0;i<3;i++) x[i]=sol->rr[i];

    for (i=0;i<MAXITR;i++) {

        /* pseudorange residuals (m) */
        nv=rescode(i,obs,n,rs,dts,vare,svh,nav,x,opt,ssat,v,H,var,azel,vsat,resp,&ns,sati,vi); 
        
        /* trace(12,"H=\n"); tracemat(12,H,n+5,NX,9,4,0); */   

        /* outlier recject based on standard normal distribution */
         if (i>=2&&nv>=NX) {
            nv=outrej_spp(nv,NX,thres,v,H,var,ssat,sati,vi,vsat,i);
        }

        if (nv<NX) {
            trace(7,"spp lack of valid sats ns=%d\n",nv);
            break;
        }

        /* trace(12,"H=\n"); tracemat(12,H,nv,NX,9,4,0);
        trace(12,"v=\n"); tracemat(12,v,nv,1,9,4,0); */

        /* weight by variance */
        for (j=0;j<nv;j++) {  
            for (k=0;k<nv;k++) {
                P[k+j*nv]=0.0;  
                if (k==j) P[k+j*nv]=1.0/var[j];
            }        
        }
        /* trace(12,"P=\n"); tracemat(12,P,nv,nv,9,4,0); */

        /* least square estimation */
        if ((info=lsq_roubst(H,v,P,NX,nv,dx,Q,Robust_OFF))) {
            trace(7,"spp lsq error info=%d\n!",info);
            break;
        }
        for (j=0;j<NX;j++) {
            x[j]+=dx[j];
        }
        if (norm(dx,NX)<1E-4) {

            sol->type=0;
            sol->time=timeadd(obs[0].time,-x[3]/CLIGHT);
            sol->dtr[0]=x[3]/CLIGHT; /* receiver clock bias (s) */
            sol->dtr[1]=x[4]/CLIGHT; /* GLO-GPS time offset (s) */
            sol->dtr[2]=x[5]/CLIGHT; /* GAL-GPS time offset (s) */
            sol->dtr[3]=x[6]/CLIGHT; /* BDS-GPS time offset (s) */
            sol->dtr[4]=x[7]/CLIGHT; /* IRN-GPS time offset (s) */
#ifdef QZSDT
            sol->dtr[5]=x[8]/CLIGHT; /* QZS-GPS time offset (s) */
#endif
            for (j=0;j<3;j++) sol->rr[j]=x[j];
            if (GINS_OFF==opt->GI_mode) for (j=0;j<3;j++) sol->rr[j+3]=0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];
            sol->qr[3]=(float)Q[1];    /* cov xy */
            sol->qr[4]=(float)Q[2+NX]; /* cov yz */
            sol->qr[5]=(float)Q[2];    /* cov zx */
            sol->ns=(uint8_t)ns;
            sol->age=sol->ratio=sol->ADOP=0.0;
            
            /* validate solution */
            if ((stat=valsol(azel,vsat,n,opt,v,P,nv,NX))) {
                sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
                /* save receiver clock (m) */
                for (j=0;j<n;j++) if (ssat) ssat[obs[j].sat-1].cdtr=x[3];
            }

            /* free memory */
            free(v);  free(H);    free(var); 
            free(P);  free(sati); free(vi);

            if (GINS_TC==opt->GI_mode&&PMODE_SINGLE==opt->mode) { tc_flag=1; break; }
            else return stat;
        }
    }

    /* if the iteration exceeds the limit or the solution fails, the solution fails flag is returned */
    if (i>=MAXITR||(SOLQ_NONE==stat&&!tc_flag)) {
        if (i>=MAXITR) trace(7,"spp: iteration over limit i=%d!\n",i);       
        
        /* free memory */
        free(v);  free(H);    free(var); 
        free(P);  free(sati); free(vi);  

        return SOLQ_NONE;
    }
    
    /* spp/ins TC mode */
    if (tc_flag) {

        mode=rtk->opt.filter;
        /* detected vehicle stationary time (s)*/
        zupt_time=rtk->ins.zupt.count*rtk->ins.interval*rtk->ins.nn;

        /* if GNSS solution fails, do not enable GNSS/INS integration mode */
        if (stat) {

            /* initialization, consider motion constraints (NHC/ZUPT/ZIHR) */
            xp=zeros(rtk->nx,1); Pp=zeros(rtk->nx,rtk->nx);
            nv=2*n+4;
            v=mat(nv,1); H=mat(nv,rtk->nx); var=mat(nv,1); R=zeros(nv,nv); sati=imat(2*n,1);  

            /* initialize clock drift */
            estvel(rtk,obs,n,rs,dts,nav,opt,sol,azel,vsat);

            /* time update of ekf states*/
            udstate_spp(rtk);
            for (i=0;i<3;i++) vc[i]=rtk->ru[i+3]; vc[3]=rtk->x[IC(6,opt)];
            /* copy states */
            matcpy(xp,rtk->x,rtk->nx,1);
            matcpy(Pp,rtk->P,rtk->nx,rtk->nx);
            /* trace(12,"P_pre=\n"); tracemat(12,Pp,rtk->nx,rtk->nx,9,4,0); */

            /* prefit residuals */
            nv=rescode_filter(rtk,obs,n,rs,dts,vare,svh,nav,rtk->ru,opt,ssat,v,H,var,azel,vsat,resp,&ns,sati);
            nv_dop=resdop_filter(rtk,obs,n,rs,dts,nav,rtk->ru,vc,azel,vsat,v+nv,H+nv*rtk->nx,var+nv,1);
            
            /* NOTE the vehicle is considered stationary only when the zero speed detection is passed, 
            the stationary state is greater than 1s and the calculated vehicle speed is less than 0.1m/s*/
            if (opt->constraint[1]&&zupt_time>1.0&&(norm(rtk->sol.rr+3,3)>0&&norm(rtk->sol.rr+3,3)<0.1)) { /* zupt*/
                nv_cons=motion_update(rtk,H,v,var,nv+nv_dop,rtk->nx,CONS_ZUPT);
                sol->iFlag=SOLF_ZUPT; /* zupt flag */
            }
            else if (opt->constraint[0]) { /* nhc */
                nv_cons=motion_update(rtk,H,v,var,nv+nv_dop,rtk->nx,CONS_NHC);        
            }
            if (opt->constraint[2]&&zupt_time>1.0&&(norm(rtk->sol.rr+3,3)>0&&norm(rtk->sol.rr+3,3)<0.1)) { /* zihr */
                nv_cons+=motion_update(rtk,H,v,var,nv+nv_dop+nv_cons,rtk->nx,CONS_ZIHR);
            }            

            /* measurement noise covariance matrix */
            for (i=0;i<(nv+nv_dop+nv_cons);i++) {
                for (j=0;j<(nv+nv_dop+nv_cons);j++) { 
                    if (i==j) R[j+i*(nv+nv_dop+nv_cons)]=var[i]; 
                }                
            }

            /* trace(12,"v=\n"); tracemat(12,v,nv+nv_dop+nv_cons,1,9,4,0);
            trace(12,"H=\n"); tracemat(12,H,nv+nv_dop+nv_cons,rtk->nx,9,4,0);
            trace(12,"Rn=\n"); tracemat(12,R,nv+nv_dop+nv_cons,nv+nv_dop+nv_cons,9,4,0); */

            /* kalman filter measurement update */
            if ((info=filter_gins(rtk,xp,Pp,H,v,R,rtk->nx,(nv+nv_dop+nv_cons),KF_GINS,mode))) {
                trace(7,"SPP/INS filter error (info=%d)\n",info);
                free(xp); free(Pp); free(v); free(H); free(R); free(var); free(sati);
                return SOLQ_NONE;
            };

            /* updates states */
            matcpy(rtk->x,xp,rtk->nx,1);
            matcpy(rtk->P,Pp,rtk->nx,rtk->nx);
            /* trace(12,"Pp=\n"); tracemat(12,rtk->P,rtk->nx,rtk->nx,9,4,0); */

            /* reset cross-covariance */
            /* init_crosscov(rtk,rtk->ins.nx,rtk->nx); */
            /* trace(12,"Pp=\n"); tracemat(12,rtk->P,rtk->nx,rtk->nx,9,4,0); */

            /* ins feedback correction */
            ins_fedback(rtk,xp);

            /* update solution status */
            update_stat(rtk,nv,SOLQ_SINGLE);    
                            
            free(xp); free(Pp); free(v); free(H); free(R); free(var); free(sati);
            return stat;
        }
        /* motion constraints (nhc/zupt) */
        else if (opt->constraint[0]||opt->constraint[1]||opt->constraint[2]) {
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
        sol_e.eventime = sol->eventime;
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
    double *rs,*dts,*var,*azel_,*resp;
    int i,j,stat,vsat[MAXOBS]={0},svh[MAXOBS],sys,fr;
    
    trace(3,"pntpos  : tobs=%s n=%d\n",time_str(obs[0].time,3),n);
    
    /* NOTE: for GNSS/INS integration, the INS solution is set to the initial state. */
    if (GINS_OFF==opt_.GI_mode) sol->stat=SOLQ_NONE;
    else sol->stat=SOLQ_INS; 
    
    if (n<=0) {
        /* if the number of available satellites is 0, output INS solution */
        if (SOLQ_INS==sol->stat) {
            rtk->outage++;
            if (GINS_LC==opt_.GI_mode||GINS_STC==opt_.GI_mode) {
                update_instat(&rtk->ins,rtk->lcgins.P,&rtk->lcgins.sol,rtk->ins.nx);
            }
            else if (GINS_TC==opt_.GI_mode) {
                update_instat(&rtk->ins,rtk->P,sol,rtk->nx);                  
            }
        }
        trace(7,"no observation data");
        return 0;
    }
    sol->time=obs[0].time;
    sol->eventime=obs[0].eventime;
    
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel_=zeros(2,n); resp=mat(1,n);
    
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            satno2id(i+1,ssat[i].id);
            sys=satsys(i+1,NULL); fr=sys2freid(sys,0,opt);
            ssat[i].sys=sys;
            ssat[i].vs=0;  /* initialize spp valid satellite flag */
            ssat[i].azel[0]=ssat[i].azel[1]=0.0;
            ssat[i].resp[fr]=ssat[i].resc[fr]=0.0;
            ssat[i].snr_rover[fr]=ssat[i].snr_base[fr]=0;
        }
        for (i=0;i<n;i++) {
            sys=satsys(obs[i].sat,NULL); 
            for (j=0;j<opt->nf;j++) {
                fr=sys2freid(sys,j,opt);
                ssat[obs[i].sat-1].snr_rover[fr]=obs[i].SNR[fr];
                ssat[obs[i].sat-1].maxsnr_rover[fr]=MAX((SNR_UNIT*obs[i].SNR[fr]),ssat[obs[i].sat-1].maxsnr_rover[fr]);                
            }         
        }
    }
    
    if (opt_.mode!=PMODE_SINGLE) { /* for precise positioning */
        opt_.sateph=EPHOPT_BRDC;
        opt_.ionoopt=IONOOPT_BRDC;
        opt_.tropopt=TROPOPT_SAAS;
    }
    /* satellite positions, velocities and clocks */
    satposs(sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);
    
    /* estimate receiver position and time with pseudorange */
    stat=estpos(rtk,obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp);

    /* output solution azel */
    outsolazel(rtk,obs,n);

    /* TC mode and GNSS unavailable, output INS solution */
    if (!stat&&GINS_TC==opt->GI_mode) {
        rtk->outage++;
        sol->stat=SOLQ_INS;
        update_instat(&rtk->ins,rtk->P,sol,rtk->nx);
        return SOLQ_NONE;
    }

    /* RAIM FDE */
    if (!stat&&n>=6&&opt->posopt[4]) {
        stat=raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp);
    }

    /* estimate receiver velocity with Doppler */
    if (stat&&(GINS_OFF==rtk->opt.GI_mode)) {
        estvel(rtk,obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);
    }
    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];
    }

    free(rs); free(dts); free(var); free(azel_); free(resp);
    return stat;
}
