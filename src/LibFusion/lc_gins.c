/*------------------------------------------------------------------------------
*lc_gins.c : GNSS/INS loose coupled integration
 *-----------------------------------------------------------------------------*/

#include "rtklib.h"

/* Integrated navigation initialization */
extern void gins_init(rtk_t *rtk, const prcopt_t *popt)
{
    ins_t *ins=&rtk->ins;
    int i,nx;
    double P[15];

    ins_init(ins,popt);

    if (GINS_LC==popt->GI_mode||GINS_STC==popt->GI_mode)
    {
        sol_t sol0={{0}};
        nx=15; 

        rtk->lcgins.x=zeros(nx,1);
        rtk->lcgins.P=zeros(nx,nx);
        rtk->lcgins.nx=nx;
        rtk->lcgins.sol=sol0;

        for (i=0;i<nx;i++)
        {
            if (i<3)              P[i]=popt->init_att_unc[i]*popt->init_att_unc[i];
            else if (i>=3&&i<6)   P[i]=popt->init_vel_unc[i-3]*popt->init_vel_unc[i-3];
            else if (i>=6&&i<9)   P[i]=popt->init_pos_unc[i-6]*popt->init_pos_unc[i-6];
            else if (i>=9&&i<12)  P[i]=popt->init_bg_unc*popt->init_bg_unc;
            else                  P[i]=popt->init_ba_unc*popt->init_ba_unc;
        }
        for (i=0;i<nx;i++) rtk->lcgins.P[i+i*nx]=P[i];
        /* trace(12,"P=\n"); tracemat(12,rtk->lcgins.P,nx,nx,9,4,0); */
    } 
    if (GINS_TC==popt->GI_mode)
    {
        nx=15; 
        for (i=0;i<nx;i++)
        {
            if (i<3)              P[i]=popt->init_att_unc[i]*popt->init_att_unc[i];
            else if (i>=3&&i<6)   P[i]=popt->init_vel_unc[i-3]*popt->init_vel_unc[i-3];
            else if (i>=6&&i<9)   P[i]=popt->init_pos_unc[i-6]*popt->init_pos_unc[i-6];
            else if (i>=9&&i<12)  P[i]=popt->init_bg_unc*popt->init_bg_unc;
            else                  P[i]=popt->init_ba_unc*popt->init_ba_unc;
        }
        for (i=0;i<nx;i++) rtk->P[i+i*rtk->nx]=P[i];
    }     
}

/* update cross-covariance -------------------------------------------
*args  :  rtk_t    *rtk   IO   rtk structure
*return:none
*-----------------------------------------------------------------------------*/
extern void update_crosscov(rtk_t *rtk)
{
    int i,j,ns=rtk->ins.nx,nx=rtk->nx;
    double *P_IG,*P_IG_,*P_GI_;

    P_IG=zeros(ns,nx-ns);P_IG_=zeros(ns,nx-ns);P_GI_=zeros(ns,nx-ns);

    pmatcpy(P_IG,ns,nx-ns,0,0,ns,nx-ns,rtk->P,nx,nx,0,ns,ns,nx);
    /* trace(12,"P_IG(k-1)=\n"); tracemat(12,P_IG,ns,nx-ns,13,6,0); */

    matmul("NN",ns,ns,nx-ns,rtk->ins.Phi,P_IG,P_IG_,1.0,0.0);
    /* MatirxT(P_IG_,P_GI_,ns,nx-ns); */
    matmul("TT",nx-ns,ns,ns,P_IG,rtk->ins.Phi,P_GI_,1.0,0.0);

    pmatcpy(rtk->P,nx,nx,0,ns,ns,nx,P_IG_,ns,nx-ns,0,0,ns,nx-ns);
    pmatcpy(rtk->P,nx,nx,ns,0,nx,ns,P_GI_,nx-ns,ns,0,0,nx-ns,ns);

    /* trace(12,"P_IG(k)=\n"); tracemat(12,P_IG_,ns,nx-ns,13,6,0);
    trace(12,"P_GI(k)=\n"); tracemat(12,P_GI_,nx-ns,ns,13,6,0); */

    free(P_IG);free(P_IG_);free(P_GI_);
}

/* INS time update*/
extern int ins_update(rtk_t *rtk)
{
    ins_t *ins=&rtk->ins;
    int nx=ins->nx;
    double *FP=zeros(nx,nx),*GQ=zeros(nx,nx),*P=zeros(nx,nx);
    double sec=ins->time.sec,interval=ins->interval*ins->nn,dt;

    /* NOTE: calculate the difference between the current time and the time update time */
    dt=fabs(sec-round((sec+interval/2.0)/ins->discretime)*ins->discretime);

    if (dt>(interval+ins->dttol)/2.0) return 0;

    /* update the state transition matrix Phi */
    phi_update(&rtk->ins,&rtk->opt);

    if (GINS_TC==rtk->opt.GI_mode)  pmatcpy(P,nx,nx,0,0,nx,nx,rtk->P,rtk->nx,rtk->nx,0,0,nx,nx);
    else matcpy(P,rtk->lcgins.P,nx,nx);

    /* if (GINS_TC==rtk->opt.GI_mode) trace(12,"Pk-1=\n"); tracemat(12,rtk->P,rtk->nx,rtk->nx,9,2,0); */
    /* trace(12,"Pk-1=\n"); tracemat(12,P,nx,nx,9,2,0); */

    /* time update */
    matmul("NN",nx,nx,nx,ins->Phi,P,FP,1.0,0.0);          /* FP=F*P */
    matmul("NT",nx,nx,nx,FP,ins->Phi,P,1.0,0.0);          /* FPF=FP*F' */

    matmul("NN",nx,nx,nx,ins->G,ins->Q,GQ,1.0,0.0);      /* GQ=G*Q */
    matmul("NT",nx,nx,nx,GQ,ins->G,P,1.0,1.0);           /* FPF=F*P*F'+G*Q*G' */
   
    if (GINS_TC==rtk->opt.GI_mode)  pmatcpy(rtk->P,rtk->nx,rtk->nx,0,0,nx,nx,P,nx,nx,0,0,nx,nx);
    else matcpy(rtk->lcgins.P,P,nx,nx);

    /* NOTE: update GNSS/INS cross-covariance!!! */
    if (GINS_TC==rtk->opt.GI_mode) update_crosscov(rtk);

    /* if (GINS_TC==rtk->opt.GI_mode) trace(12,"P_pre=\n"); tracemat(12,rtk->P,rtk->nx,rtk->nx,9,2,0); */
    /* trace(12,"P_pre=\n"); tracemat(12,P,nx,nx,9,2,0); */

    free(FP); free(GQ); free(P);
}

/* GNSS/INS loosely coupled integration */
extern int lc_gins(rtk_t *rtk)
{
    ins_t *ins=&rtk->ins;
    sol_t *sol=&rtk->lcgins.sol;
    prcopt_t *popt=&rtk->opt;
    int i,j,nx=rtk->lcgins.nx,nv=3,nv_cons=0,info,stat=rtk->sol.stat,mode=rtk->opt.filter;
    double p_ins[3],p_gnss[3],iFrp[9],dp[3],zupt_time;
    double lever_n[3],lever_nx[9],Re[9],Rn[9];
    double *I3,*x,*P,*xp,*Pp,*v,*H,*var,*R;

    /* check GNSS status and output INS navigation information if GNSS is unavailable */
    if (SOLQ_INS==rtk->sol.stat) {
        rtk->outage++;
        sol->stat=SOLQ_INS;
        update_instat(popt,ins,rtk->lcgins.P,sol,nx);
        return 1;
    }

    /* detected vehicle stationary time span (s)*/
    zupt_time=ins->zupt.count*ins->interval*ins->nn;

    /* initialize heap memory, consider NHC/ZUPT/ZIHR constraints */
    I3=eye(3); x=zeros(nx,1); P=zeros(nx,nx); xp=zeros(nx,1); Pp=zeros(nx,nx);
    v=zeros(nv+4,1); H=zeros(nv+4,nx); var=mat(nv+4,1); R=zeros(nv+4,nv+4);

    /* initialize states */
    matcpy(P,rtk->lcgins.P,nx,nx);
    matcpy(iFrp,ins->eth.Frp,3,3);
    /* trace(12,"P_pre=\n"); tracemat(12,P,nx,nx,9,4,0); */
    earth_update(popt,ins->pos,ins->vel,&ins->eth);
    matinv(iFrp,3);

    /* lever arm correction to convert INS position to GNSS position */
    ins2gnss(popt,&rtk->ins,p_ins,3);
    ecef2pos(rtk->sol.rr,p_gnss);

    Mat3mulv(1.0,ins->Cnb,ins->lever,lever_n);
    vskew(1.0,lever_n,lever_nx);

    /* measurement vector */
    for (i=0;i<3;i++) dp[i]=p_ins[i]-p_gnss[i];
    Mat3mulv(1.0,iFrp,dp,v);

    /* measurement matrix H */
    for (i=0;i<nv;i++){
        for (j=0;j<nx;j++){
            if (j<3)        H[j+i*nx]=lever_nx[j+i*3];
            if (j>=6&&j<9)  H[j+i*nx]=I3[(j-6)+i*3];
        }
    }
    
    /* initialize measurement variance */
    soltocov(&rtk->sol,Re);
    covenu(ins->pos,Re,Rn);
    for (i=0;i<nv;i++) var[i]=Rn[i+i*nv];

    /* motion constraints */
    /* NOTE: the vehicle is considered stationary only when the zero speed detection is passed, 
    the stationary state is greater than 1s and the calculated vehicle speed is less than 0.1m/s */
    if (popt->constraint[1]&&zupt_time>1.0&&(norm(rtk->sol.rr+3,3)>0&&norm(rtk->sol.rr+3,3)<0.1)) { /* zupt*/
        nv_cons=motion_update(rtk,H,v,var,nv,nx,CONS_ZUPT);
        sol->iFlag=SOLF_ZUPT; /* zupt flag */
    }
    else if (popt->constraint[0]) { /* nhc */
        nv_cons=motion_update(rtk,H,v,var,nv,nx,CONS_NHC);      
    }
    if (popt->constraint[2]&&zupt_time>1.0&&(norm(rtk->sol.rr+3,3)>0&&norm(rtk->sol.rr+3,3)<0.1)) { /* zihr */
        nv_cons+=motion_update(rtk,H,v,var,nv+nv_cons,nx,CONS_ZIHR);
    }     

    /* measurement noise covariance matrix R*/
    diag_Cov(nv+nv_cons,var,R,diag_var);
    /* trace(12,"H=\n"); tracemat(12,H,nv,nx,15,10,0); */
    /* trace(12,"Rn=\n"); tracemat(12,Rn,3,3,9,4,0); */

    /* measurement update of ekf states */
    if ((info=filter_(rtk,x,P,H,v,R,nx,nv+nv_cons,xp,Pp,mode))) {
        trace(2,"lc_gins (%d) filter error info=%d\n",i+1,info);
        stat=SOLQ_NONE;
    }   
    /* trace(12,"Pk=\n"); tracemat(12,Pp,nx,nx,9,4,0); */

    /* update state covariance matrix */
    matcpy(rtk->lcgins.P,Pp,nx,nx);

    /* INS feedback correction */
    ins_fedback(rtk,xp);

    /* save solution status */
    update_lcstat(rtk,stat);

    /* free heap memory */
    free(I3);free(x); free(P); free(xp); free(Pp);
    free(v); free(H); free(var); free(R);

    return 1;
}

extern void conpad_fedback(ins_t *ins, double *dw, double *dv, double *da_con, double *dv_rot, double *dv_pad)
{
    double temp1[3],temp2[3];

    /* single sample + previous*/
    if (1==ins->nn){
        matcpy(dw,ins->dw,3,1);
        matcpy(dv,ins->dv,3,1);
        /* rotation error compensation term */
        vskewmv(1.0/2.0,ins->dw,ins->dv,dv_rot);

        /* paddling error compensation term */
        vskewmv(1.0/12.0,ins->p1dw,ins->dv,temp1);
        vskewmv(1.0/12.0,ins->p1dv,ins->dw,temp2);
        vnadd(3,temp1,1.0,temp2,1.0,dv_pad);

        /* second-order cone error compensation term */
        vskewmv(1.0/12.0,ins->p1dw,ins->dw,temp1);
        vnadd(3,ins->dw,1.0,temp1,1.0,da_con);
    }
    /* double sample */
    else if (2==ins->nn){
        vnadd(3,ins->dw,1.0,ins->n1dw,1.0,dw);
        vnadd(3,ins->dv,1.0,ins->n1dv,1.0,dv);

        /* rotation error compensation term */
        vskewmv(1.0/2.0,dw,dv,dv_rot);

        /* paddling error compensation term */
        vskewmv(2.0/3.0,ins->dw,ins->n1dv,temp1);
        vskewmv(2.0/3.0,ins->dv,ins->n1dw,temp2);
        vnadd(3,temp1,1.0,temp2,1.0,dv_pad);

        /* second-order cone error compensation term */
        vskewmv(2.0/3.0,ins->dw,ins->n1dw,temp1);
        vnadd(3,dw,1.0,temp1,1.0,da_con);        
    }
}

/* INS error feedback correction*/
extern void imu_fedback(ins_t *ins, imud_t *imu)
{
    /* correct gyroscope and accelerometer zero bias */
    vnadd(3,imu[0].dw,1.0,ins->bg,-ins->interval,ins->dw);
    vnadd(3,imu[0].dv,1.0,ins->ba,-ins->interval,ins->dv);        

    if (2==ins->nn){
        vnadd(3,imu[1].dw,1.0,ins->bg,-ins->interval,ins->n1dw);
        vnadd(3,imu[1].dv,1.0,ins->ba,-ins->interval,ins->n1dv);        
    }
}

/* convert psi error state to phi error state ------------------------------
*args   : ins_t *ins           I   ins structure
*         const double *dr     I   dblh pos error  (3x1)
*         double *dx           IO  ins error state (15x1)
*return : none
*-------------------------------------------------------------------------------*/
extern void  psi2phi_corr(ins_t *ins, const double *dr, double *dx)
{
    double d_ceta[3],dv[3];
    int i;

    /* equivalent rotation vector phi_nc */
    d_ceta[0]=-dr[0]; d_ceta[1]=dr[1]*cos(ins->pos[0]); d_ceta[2]=dr[1]*sin(ins->pos[0]);

    /* convert Psi attitude misalignment angle (psi_cn') to Phi attitude misalignment angle (phi_nn') , phi_nn'=phi_nc+psi_cn' */     
    for (i=0;i<3;i++) dx[i]=dx[i]+d_ceta[i];

    /* convert Psi velocity error state to Phi velocity error state */
    vskewmv(1.0,d_ceta,ins->vel,dv);
    for (i=0;i<3;i++) dx[i+3]=dx[i+3]-dv[i];

}

/* INS error feedback correction */
extern void ins_fedback(rtk_t *rtk, double *dx)
{
    prcopt_t *popt=&rtk->opt;
    ins_t *ins=&rtk->ins;
    int i;
    double dr[3],phi[9],Cnn_[9];
    double *I3=eye(3),Cnb[9]={0.0};
    double qnn_[4],qn_b[4],phi_nn_[3];

    /* convert dxyz to dblh */
    earth_update(popt,ins->pos,ins->vel,&ins->eth);
    Mat3mulv(1.0,ins->eth.Frp,dx+6,dr);

    /* NOTE: convert psi error state to phi error state */
    if (ERR_PSI==popt->err_model) psi2phi_corr(ins,dr,dx);

    /* qnb=qnn_°qn_b */
    for (i=0;i<4;i++) qn_b[i]=ins->qnb[i];
    for (i=0;i<3;i++) phi_nn_[i]=dx[i];
    
    rv2quat(1.0,phi_nn_,qnn_);
    quatmul(qnn_,qn_b,ins->qnb); 
    qnbnorm(ins->qnb); /* normalize qnb */
    qnb2Cnb(ins->qnb,ins->Cnb);
    Cnb2att(ins->Cnb,ins->att);

    /* Cnb=(I+[phi x])Cn'b */
    /* for (i=0;i<9;i++) Cnb[i]=ins->Cnb[i];

    vskew(1.0,dx,phi);
    Mat3add2(I3,1.0,phi,1.0,Cnn_);
    Mat3mul2(1.0,Cnn_,Cnb,ins->Cnb); 
    Cnb2att(ins->Cnb,ins->att);
    att2qnb(ins->att,ins->qnb); */

    for (i=0;i<ins->nx;i++){
        if (i>=3&&i<6)      ins->vel[i-3]-=dx[i];
        if (i>=6&&i<9)      ins->pos[i-6]-=dr[i-6];
        if (i>=9&&i<12)     ins->bg[i-9] +=dx[i];
        if (i>=12&&i<15)    ins->ba[i-12]+=dx[i];
    }        

    for (i=0;i<3;i++){
        /* update previous epoch pos/vel by kf updated state */
        ins->p1pos[i]=ins->pos[i];         
        ins->p1vel[i]=ins->vel[i];       
    }

    free(I3);
}

/* INS error feedback correction*/
extern void ins_fedback_fix(rtk_t *rtk, double *dx)
{
    prcopt_t *popt=&rtk->opt;
    ins_t *ins=&rtk->ins;
    int i;
    double dr[3],phi[9],Cnn_[9];
    double *I3=eye(3),Cnb[9],Cnb_[9];
    double qnn_[4],qn_b[4],qnb[4],phi_nn_[3];

    /* convert dxyz to dblh */
    earth_update(popt,ins->pos,ins->vel,&ins->eth);
    Mat3mulv(1.0,ins->eth.Frp,dx+6,dr);

    /* NOTE: convert psi error state to phi error state */
    if (ERR_PSI==popt->err_model) psi2phi_corr(ins,dr,dx);

    /* qnb=qnn_°qn_b */
    for (i=0;i<4;i++) qn_b[i]=ins->qnb[i];
    for (i=0;i<3;i++) phi_nn_[i]=dx[i];
    
    rv2quat(1.0,phi_nn_,qnn_);
    quatmul(qnn_,qn_b,qnb); 
    qnbnorm(qnb); /* normalize qnb */
    qnb2Cnb(qnb,Cnb);
    Cnb2att(Cnb,ins->xa);

    /* Cnb=(I+[phi x])Cn'b */
    /* for (i=0;i<9;i++) Cnb=ins->Cnb[i];

    vskew(1.0,dx,phi);
    Mat3add2(I3,1.0,phi,1.0,Cnn_);
    Mat3mul2(1.0,Cnn_,Cnb,Cnb_); 
    Cnb2att(Cnb_,ins->xa); */

    for (i=0;i<ins->nx;i++){
        if (i>=3&&i<6)      ins->xa[i]=ins->vel[i-3]-dx[i];
        if (i>=6&&i<9)      ins->xa[i]=ins->pos[i-6]-dr[i-6];
        if (i>=9&&i<12)     ins->xa[i]=ins->bg[i-9]+dx[i];
        if (i>=12&&i<15)    ins->xa[i]=ins->ba[i-12]+dx[i];
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
        Mat3mulv(1.0*sgn,Cen,ins->vel,ve);
    }

    for (i=0;i<3;i++){
        sol->rr [i]=re[i];
        sol->vel[i]=ve[i];
        sol->att[i]=ins->att[i]*R2D;
        sol->bg [i]=ins->bg [i]*R2D*3600*sgn;
        sol->ba [i]=ins->ba [i]*1E5; 
    }

    /* converts the yaw from counterclockwise to clockwise */
    if (sol->att[2]<=0) sol->att[2]=-sol->att[2];
    else sol->att[2]=360.0-sol->att[2];

    for (i=0;i<3;i++){
        for (j=0;j<3;j++){
            Qa[j+i*3]=P[j+i*nx]*R2D*R2D;    /* deg^2 */
            Qvn[j+i*3]=P[(j+3)+(i+3)*nx];   /* m^2/s^2 */
            Qrn[j+i*3]=P[(j+6)+(i+6)*nx];   /* m^2 */
            Qbg[j+i*3]=P[(j+9)+(j+9)*nx]*(R2D*3600)*(R2D*3600);  /* deg^2/h^2*/
            Qba[j+i*3]=P[(j+12)+(j+12)*nx]*1E5*1E5;              /* ug^2*/
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
    if (stat!=SOLQ_NONE) sol->stat=stat;
    if (PMODE_KINEMA==rtk->opt.mode) sol->ratio=rtk->sol.ratio;
    if (SOLQ_INS==stat) {
        sol->time=ins->time;
        sol->ns=0;
    }
    else {
        /* if GNSS/INS integration solution is available, reset GNSS outage count to 0 */
        if (rtk->outage<=MAX_OUTIME) rtk->outage=0;
        sol->time=rtk->sol.time;
        sol->ns=rtk->sol.ns;      
    }

}