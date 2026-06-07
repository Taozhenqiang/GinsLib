/*------------------------------------------------------------------------------
*lc_gins.c : GNSS/INS loose coupled integration
 *-----------------------------------------------------------------------------*/

#include "rtklib.h"

extern int solflags(sol_t *sol)
{
    return (SOLQ_NONE!=sol->stat);
}

extern int isGNSS(const prcopt_t *popt)
{
    return (GINS_OFF==popt->GI_mode);
}

extern int isGINS(const prcopt_t *popt)
{
    return (GINS_LC==popt->GI_mode||GINS_TC==popt->GI_mode||GINS_STC==popt->GI_mode);
}

extern int isGINS_LC(const prcopt_t *popt)
{
    return (GINS_LC==popt->GI_mode||GINS_STC==popt->GI_mode||PMODE_LC_POS==popt->mode);
}

extern int is_motionconstraints(const prcopt_t *popt)
{
    return (popt->constraint[0]||popt->constraint[1]||popt->constraint[2]);
}

extern int isNHC(const prcopt_t *popt)
{
    return (popt->constraint[0]);
}

/* reset INS state parameters */
extern void reset_instat(rtk_t *rtk)
{
    int i;
    for (i=0;i<rtk->ins.nx;i++) rtk->x[i]=0.0;
}

/* Jacobian matrix for pos/vel/att */
static void Jacobi_avp(rtk_t *rtk, int nx, int nv, double *H)
{
    ins_t *ins=&rtk->ins;
    int i,j;
    double lever_n[3],lever_nx[9],*I3;

    I3=eye(3);

    Mat3mulv(1.0,ins->Cnb,ins->lever,lever_n);
    vskew(1.0,lever_n,lever_nx);

    /* measurement matrix H */
    for (i=0;i<nv;i++){
        for (j=0;j<nx;j++){
            if (j<3)        H[j+i*nx]=lever_nx[j+i*3];
            if (j>=6&&j<9)  H[j+i*nx]=I3[(j-6)+i*3];
        }
    }

    free(I3);
}

/* update INS solution state */
extern void update_instat(const prcopt_t *popt, ins_t *ins, double *P, sol_t *sol, int nx)
{
    int i,j;
    double re[3],ve[3],Qa[9],Qvn[9],Qv[9],Qrn[9],Qbg[9],Qba[9],Qr[9],Cne[9],Cen[9];
    double p_gnss[6];
    double sgn=(SOLTYPE_BACKWARD==popt->reverse?-1.0:1.0); 

    /* solution status */
    sol->time=ins->time;
    
    /* convert ins pos/vel to GNSS pos/vel */
    /* NOTE: in backward mode, the GNSS/INS velocity and gyroscope bias have opposite signs to the actual values */
    if (OUTPOS_GNSS==popt->outpos) {
        ins2gnss(popt,ins,p_gnss,6);

        pos2ecef(p_gnss,re);
        xyz2enu(p_gnss,Cne);
        DCMT(Cne,Cen);
        Mat3mulv(1.0,Cen,p_gnss+3,ve);  
    }
    else {
        pos2ecef(ins->pos,re);
        xyz2enu(ins->pos,Cne);
        DCMT(Cne,Cen);
        Mat3mulv(1.0,Cen,ins->vel,ve);
    }

    /* save GNSS/INS pos/vel/att/bg/ba */
    for (i=0;i<3;i++){
        sol->rr [i]=re[i];
        sol->vel[i]=ve[i]*sgn;
        sol->att[i]=ins->att[i]*R2D;
        sol->bg [i]=ins->bg [i]*R2D*3600*sgn;
        sol->ba [i]=ins->ba [i]*1E5; 
    }
    for (i=0;i<4;i++) sol->qnb[i]=ins->qnb[i];

    /* converts the yaw from counterclockwise to clockwise */
    if (sol->att[2]<=0) sol->att[2]=-sol->att[2];
    else sol->att[2]=360.0-sol->att[2];

    for (i=0;i<3;i++){
        for (j=0;j<3;j++){
            Qa[j+i*3]=P[j+i*nx];            /* rad^2 */
            Qvn[j+i*3]=P[(j+3)+(i+3)*nx];   /* m^2/s^2 */
            Qrn[j+i*3]=P[(j+6)+(i+6)*nx];   /* m^2 */
            Qbg[j+i*3]=P[(j+9)+(i+9)*nx];   /* rad^2/s^2*/
            Qba[j+i*3]=P[(j+12)+(i+12)*nx]; /* g^2*/
        }
    }
    
    /* cov of local frame to ecef frame */
    covecef(ins->pos,Qrn,Qr);
    covecef(ins->pos,Qvn,Qv);
    covtosol_att(Qa,sol);
    covtosol_vel(Qv,sol);
    covtosol(Qr,sol);
    covtosol_bga(Qbg,Qba,sol);
}

/* update GNSS/INS LC solution state */
extern void update_lcstat(rtk_t *rtk, int stat){

    ins_t *ins=&rtk->ins;
    sol_t *sol=&rtk->lcgins.sol;
    int i,j,nx=rtk->lcgins.nx;

    /* update ins state */
    update_instat(&rtk->opt,ins,rtk->lcgins.P,sol,nx);

    /* solution status */
    sol->stat=stat;
    if (SOLQ_INS==stat||SOLQ_CONS==stat) {
        sol->time=ins->time;
        sol->ns=0;
        sol->ratio=0.0;
        for (i=0;i<4;i++) sol->dop[i]=0.0;
    }
    else {
        /* if GNSS/INS integration solution is available, reset GNSS outage count to 0 */
        if (rtk->outage<=ins->max_outime) rtk->outage=0;
        sol->time=rtk->sol.time;
        sol->ns=rtk->sol.ns;
        sol->ratio=rtk->sol.ratio;
        matcpy(sol->dop,rtk->sol.dop,4,1);   
    }
}

/* check GNSS solution quality */
static int quality_check(rtk_t *rtk)
{
    int i;
    double posvar=0.0,posvar_thres=10.0;

    /* calc average position variance, will skip LC if too high */
    for (i=0;i<3;i++) posvar+=rtk->P[i+i*rtk->nx];
    posvar/=3.0; 

    if (posvar>posvar_thres||rtk->sol.ns<4) {
        trace(7,"GNSS solution quality check failed(LC), posvar=%f, ns=%d\n",posvar,rtk->sol.ns);
        return 0;
    }
    else return 1;
}

/* GNSS/INS loosely coupled integration measurement (H/v/R) */
static void LCI_meas(rtk_t *rtk, double *H, double *v, double *var, int nx, int nv)
{
    ins_t *ins=&rtk->ins;
    sol_t *sol=&rtk->lcgins.sol;
    prcopt_t *popt=&rtk->opt;
    int i;
    double p_ins[3],p_gnss[3],iFrp[9],dp[3],Re[9],Rn[9];

    earth_update(popt,ins->pos,ins->vel,&ins->eth);
    matcpy(iFrp,ins->eth.Frp,3,3);
    matinv(iFrp,3);

    /* lever arm correction to convert INS position to GNSS position */
    ins2gnss(popt,ins,p_ins,3);
    ecef2pos(rtk->sol.rr,p_gnss);

    /* measurement vector */
    for (i=0;i<3;i++) dp[i]=p_ins[i]-p_gnss[i];
    Mat3mulv(1.0,iFrp,dp,v);

    /* measurement matrix H */
    Jacobi_avp(rtk,nx,nv,H);
    
    /* initialize measurement variance */
    soltocov(&rtk->sol,Re); 
    covenu(ins->pos,Re,Rn); covtodiag(Rn,3);
    for (i=0;i<nv;i++) var[i]=Rn[i+i*nv];  
}

/* GNSS/INS loosely coupled integration */
extern int lc_gins(rtk_t *rtk)
{
    ins_t *ins=&rtk->ins;
    sol_t *sol=&rtk->lcgins.sol;
    prcopt_t *popt=&rtk->opt;
    int nx=rtk->lcgins.nx,nv,nv_cons=0,info,gnss_stat=0,stat=rtk->sol.stat,mode=popt->lcfilter;
    double *I3,*x,*P,*xp,*Pp,*v,*H,*var,*R;

    /* check GNSS solution quality */
    gnss_stat=(solflags(&rtk->sol)&&quality_check(rtk));

    /* check GNSS status and output INS navigation information if GNSS is unavailable */
    if (!gnss_stat&&(!popt->constraint[0]&&!popt->constraint[1])) {
        rtk->outage++;
        sol->stat=SOLQ_INS;
        update_instat(popt,ins,rtk->lcgins.P,sol,nx);
        return 1;
    }

    /* number of GNSS pos measurements */
    nv=(gnss_stat)?3:0;

    /* initialize heap memory, consider NHC/ZUPT/ZIHR constraints (max num=ZUPT+ZIHR=4) */
    x=zeros(nx,1); P=zeros(nx,nx); xp=zeros(nx,1); Pp=zeros(nx,nx);
    v=zeros(nv+4,1); H=zeros(nv+4,nx); var=mat(nv+4,1); R=zeros(nv+4,nv+4);

    /* initialize states */
    matcpy(P,rtk->lcgins.P,nx,nx);

    /* if GNSS is available, don't using GNSS/INS LC */
    if (gnss_stat) {
        LCI_meas(rtk,H,v,var,nx,nv);
    }
    else stat=SOLQ_CONS;

    /* motion constraints */
    nv_cons=motion_meas(rtk,popt,H,v,var,nv,nx);

    /* measurement noise covariance matrix R*/
    diag_Cov(nv+nv_cons,var,R,diag_var);

    /* measurement update of ekf states */
    if ((info=filter_(rtk,x,P,H,v,R,nx,nv+nv_cons,xp,Pp,mode))) {
        trace(2,"lc_gins filter error info=%d\n",info);
        stat=SOLQ_INS;
    }   
    /* tracefilter(12,TRAE_R|TRAE_H|TRAE_Ppre|TRAE_Pp|TRAE_v|TRAE_xpre|TRAE_xp,nx,nv+nv_cons,R,H,P,Pp,v,x,xp); */

    /* update state covariance matrix */
    matcpy(rtk->lcgins.P,Pp,nx,nx);

    /* INS feedback correction */
    ins_fedback(rtk,xp);

    /* save solution status */
    update_lcstat(rtk,stat);

    /* free heap memory */
    free(x); free(P); free(xp); free(Pp);
    free(v); free(H); free(var); free(R);

    return 1;
}