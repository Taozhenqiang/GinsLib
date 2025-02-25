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

    if (GINS_LC==popt->GI_mode)
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
        for (i=0;i<nx;i++) rtk->P[i+i*nx]=P[i];
    }     
}

/* INS time update*/
extern int ins_update(rtk_t *rtk)
{
    ins_t *ins=&rtk->ins;
    int nx=ins->nx;
    double *FP=zeros(nx,nx),*GQ=zeros(nx,nx),*P=zeros(nx,nx);
    double sec=ins->time.sec,interval=ins->interval*ins->nn,dt;

    dt=fabs(sec-round(sec*1e1)/1e1);

    if (dt>(interval+ins->dttol)/2.0) return 0;

    phi_update(&rtk->ins);

    if (GINS_TC==rtk->opt.GI_mode)  matcpy(P,rtk->P,nx,nx);
    else matcpy(P,rtk->lcgins.P,nx,nx);
    trace(12,"Pk-1=\n"); tracemat(12,P,nx,nx,9,4,0);

    matmul("NN",nx,nx,nx,ins->Phi,P,FP,1.0,0.0);          /* FP=F*P */
    matmul("NT",nx,nx,nx,FP,ins->Phi,P,1.0,0.0);          /* FPF=FP*F' */

    matmul("NN",nx,nx,nx,ins->G,ins->Q,GQ,1.0,0.0);      /* GQ=G*Q */
    matmul("NT",nx,nx,nx,GQ,ins->G,P,1.0,1.0);           /* FPF=F*P*F'+G*Q*G' */
   
    if (GINS_TC==rtk->opt.GI_mode)  matcpy(rtk->P,P,nx,nx);
    else matcpy(rtk->lcgins.P,P,nx,nx);
    trace(12,"P_pre=\n"); tracemat(12,P,rtk->nx,rtk->nx,9,4,0);

    free(FP); free(GQ); free(P);
}

/* GNSS/INS loosely coupled integration */
extern void lc_gins(rtk_t *rtk){

    ins_t *ins=&rtk->ins;
    int i,j,nx=rtk->lcgins.nx,nv=3,info,stat=SOLQ_LC;
    double p_ins[3],p_gnss[3],iFpv[9],dp[3];
    double lever_n[3],lever_nx[9],Re[9];
    double *I3,*x,*P,*xp,*Pp,*v,*H,*Rn;

    I3=eye(3);x=zeros(nx,1); P=zeros(nx,nx); xp=zeros(nx,1); Pp=zeros(nx,nx);
    v=zeros(nv,1); H=zeros(nv,nx); Rn=mat(nv,nv);

    matcpy(P,rtk->lcgins.P,nx,nx);
    matcpy(iFpv,ins->eth.Fpv,3,3);
    /* trace(12,"P_pre=\n"); tracemat(12,P,nx,nx,9,4,0); */
    earth_update(ins->pos,ins->vel,&ins->eth);
    matinv(iFpv,3);

    ins2gnss(rtk,p_ins,3);
    ecef2pos(rtk->sol.rr,p_gnss);  

    Mat3mulv(1.0,ins->Cnb,ins->lever,lever_n);
    vskew(1.0,lever_n,lever_nx);

    for (i=0;i<3;i++) dp[i]=p_ins[i]-p_gnss[i];
    Mat3mulv(1.0,iFpv,dp,v);

    for (i=0;i<nv;i++){
        for (j=0;j<nx;j++){
            if (j<3)        H[j+i*nx]=lever_nx[j+i*3];
            if (j>=6&&j<9)  H[j+i*nx]=I3[j-6+i*3];
        }
    }
    /* trace(12,"H=\n"); tracemat(12,H,nv,nx,15,10,0); */
    soltocov(&rtk->sol,Re);
    covenu(ins->pos,Re,Rn);
    /* trace(12,"Rn=\n"); tracemat(12,Rn,3,3,9,4,0); */

    /* measurement update of ekf states */
    if ((info=filter_(x,P,H,v,Rn,nx,nv,xp,Pp))) {
        trace(2,"lc_gins (%d) filter error info=%d\n",i+1,info);
        stat=SOLQ_NONE;
    }   
    /* trace(12,"Pk=\n"); tracemat(12,Pp,nx,nx,9,4,0); */
    matcpy(rtk->lcgins.P,Pp,nx,nx);
    ins_fedback(rtk,xp);

    update_lcstat(rtk,stat);

    free(I3);free(x); free(P); free(xp); free(Pp);
    free(v); free(H); free(Rn);
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

/* INS error feedback correction*/
extern void ins_fedback(rtk_t *rtk, double *dx)
{
    ins_t *ins=&rtk->ins;
    int i;
    double dr[3],phi[9],Cnn_[9];
    double *I3=eye(3),*Cnb=zeros(3,3);

    matcpy(Cnb,ins->Cnb,3,3);

    vskew(1.0,dx,phi);
    Mat3add2(I3,1.0,phi,1.0,Cnn_);
    Mat3mul2(1.0,Cnn_,Cnb,ins->Cnb); /*Cnb=(I+[phi x])Cn'b*/
    Cnb2att(ins->Cnb,ins->att);

    earth_update(ins->pos,ins->vel,&ins->eth);
    Mat3mulv(1.0,ins->eth.Fpv,dx+6,dr);

    for (i=0;i<ins->nx;i++){
        if (i>=3&&i<6)      ins->vel[i-3]-=dx[i];
        if (i>=6&&i<9)      ins->pos[i-6]-=dr[i-6];
        if (i>=9&&i<12)     ins->bg[i-9] +=dx[i];
        if (i>=12&&i<15)    ins->ba[i-12]+=dx[i];
    }

    for (i=0;i<3;i++){
        ins->p1pos[i]=ins->pos[i];         
        ins->p1vel[i]=ins->vel[i];       
    }

    free(I3);free(Cnb);
}

/* update GNSS/INS LC solution status */
extern void update_lcstat(rtk_t *rtk, int stat){

    ins_t *ins=&rtk->ins;
    sol_t *sol=&rtk->lcgins.sol;
    int i,j,nx=rtk->lcgins.nx;
    double re[3],ve[3],Qa[9],Qvn[9],Qv[9],Qrn[9],Qbg[9],Qba[9],Qr[9],Cne[9],Cen[9];

    if (stat!=SOLQ_NONE) sol->stat=stat;

    sol->time=rtk->sol.time;
    sol->ns=rtk->sol.ns;

    pos2ecef(ins->pos,re);
    xyz2enu(ins->pos,Cne);
    DCMT(Cne,Cen);
    Mat3mulv(1.0,Cen,ins->vel,ve);

    for (i=0;i<3;i++){
        sol->rr [i]=re[i];
        sol->vel[i]=ve[i];
        sol->att[i]=ins->att[i]*R2D;
        sol->bg [i]=ins->bg [i]*R2D*3600;
        sol->ba [i]=ins->ba [i]*1E5; 
    }

    /* Converts the yaw from counterclockwise to clockwise */
    sol->att[2]=-sol->att[2];

    for (i=0;i<3;i++){
        for (j=0;j<3;j++){
            Qa[j+i*3]=rtk->lcgins.P[j+i*nx]*R2D*R2D;    /* deg^2 */
            Qvn[j+i*3]=rtk->lcgins.P[(j+3)+(i+3)*nx];   /* m^2/s^2 */
            Qrn[j+i*3]=rtk->lcgins.P[(j+6)+(i+6)*nx];   /* m^2 */
            Qbg[j+i*3]=rtk->lcgins.P[(j+9)+(j+9)*nx]*(R2D*3600)*(R2D*3600);  /* deg^2/h^2*/
            Qba[j+i*3]=rtk->lcgins.P[(j+12)+(j+12)*nx]*1E5*1E5;              /* ug^2*/
        }
    }
    
    /* cov of local frame to ecef frame */
    covecef(ins->pos,Qrn,Qr);
    covecef(ins->pos,Qvn,Qv);
    covtosol_att(Qa,sol);
    covtosol(Qr,sol);
    covtosol_vel(Qv,sol);
    covtosol_bga(Qbg,Qba,sol);
}