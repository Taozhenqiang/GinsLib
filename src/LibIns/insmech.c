/*-----------------
 * ins mechanization and error model
------------------*/
#include "rtklib.h"

/* initialize earth related parameters -----------------------------------------------*/
extern void earth_init(const double *pos, const double *vel, eth_t *eth)
{   
    int i;
    double sinB2,sin2B2,secB,Rmh,Rnh,temp1[3];
    sinB2=sin(pos[0])*sin(pos[0]);
    sin2B2=sin(2.0*pos[0])*sin(2.0*pos[0]);
    secB=1.0/cos(pos[0]);

    eth->g0=G0;
    eth->Rl=RE_WGS84;
    eth->alpha=FE_WGS84;
    eth->Rs=eth->Rl*(1.0-eth->alpha);
    eth->e1=sqrt(eth->Rl*eth->Rl-eth->Rs*eth->Rs)/eth->Rl;
    eth->e2=sqrt(eth->Rl*eth->Rl-eth->Rs*eth->Rs)/eth->Rs;
    eth->wie=OMGE;
    /* gravity related parameters, ref psins */
    eth->beta[0]=5.2790414E-3;
    eth->beta[1]=2.32718E-5;
    eth->beta[2]=3.086E-6;
    eth->beta[3]=8.08E-9;

    eth->RN=eth->Rl/sqrt(1.0-eth->e1*eth->e1*sinB2);
    eth->RM=eth->RN*(1.0-eth->e1*eth->e1)/(1.0-eth->e1*eth->e1*sinB2);
    Rmh=eth->RM+pos[2];Rnh=eth->RN+pos[2];

    /* the default is forward mechanization */
    eth->wnie[0]=0.0; eth->wnie[1]=eth->wie*cos(pos[0]); eth->wnie[2]=eth->wie*sin(pos[0]);
    eth->wnen[0]=-vel[1]/Rmh; 
    eth->wnen[1]= vel[0]/Rnh; 
    eth->wnen[2]= vel[0]*sin(pos[0])/(Rnh*cos(pos[0])); 
    vnadd(3,eth->wnie,1.0,eth->wnen,1.0,eth->wnin);
    for (i=0;i<9;i++)
    {
        eth->F1 [i]=0.0; 
        eth->F2 [i]=0.0;
        eth->F3 [i]=0.0;
        eth->F4 [i]=0.0;
        eth->Fav[i]=0.0;
        eth->Frr[i]=0.0;        
        eth->Frp[i]=0.0;
    }

    eth->g=eth->g0*(1+eth->beta[0]*sinB2-eth->beta[1]*sin2B2)-eth->beta[2]*pos[2];
    eth->gn[0]=0.0; eth->gn[1]=0.0; eth->gn[2]=-1.0*eth->g;

    vnadd(3,eth->wnie,2.0,eth->wnen,1.0,eth->wnien);
    vskewmv(1.0,eth->wnien,vel,temp1);
    /* harmful acceleration */
    vnadd(3,eth->gn,1.0,temp1,-1.0,eth->gcc);

}

/* update earth related parameters -----------------------------------------------*/
extern void earth_update(const prcopt_t *popt, const double *pos, const double *vel, eth_t *eth)
{ 
    int i;
    double sinB,cosB,cos2B,tanB,sin2B,sinB2,sin2B2,secB,secB2,Rmh,Rnh,Rmh2,Rnh2,temp1[3];
    sinB=sin(pos[0]);sin2B=sin(2.0*pos[0]);sinB2=sinB*sinB;sin2B2=sin2B*sin2B;
    cosB=cos(pos[0]);cos2B=cos(2.0*pos[0]);tanB=tan(pos[0]);
    secB=1.0/cosB;secB2=secB*secB;

    eth->RN=eth->Rl/sqrt(1.0-eth->e1*eth->e1*sinB2);
    eth->RM=eth->RN*(1.0-eth->e1*eth->e1)/(1.0-eth->e1*eth->e1*sinB2);
    Rmh=eth->RM+pos[2];Rnh=eth->RN+pos[2];
    Rmh2=Rmh*Rmh;Rnh2=Rnh*Rnh;

    if (MECH_FORWARD==popt->reverse) {
        eth->wnie[0]=0.0; eth->wnie[1]=eth->wie*cos(pos[0]); eth->wnie[2]=eth->wie*sin(pos[0]);       
    }
    /* backward mechanization negates the wie symbol */
    else if (MECH_BACKWARD==popt->reverse) {
        eth->wnie[0]=0.0; eth->wnie[1]=-eth->wie*cos(pos[0]); eth->wnie[2]=-eth->wie*sin(pos[0]);  
    }

    /* wnen */
    eth->wnen[0]=-vel[1]/Rmh; 
    eth->wnen[1]= vel[0]/Rnh; 
    eth->wnen[2]= vel[0]*sin(pos[0])/(Rnh*cos(pos[0])); 
    /* wnin */
    vnadd(3,eth->wnie,1.0,eth->wnen,1.0,eth->wnin);

    for (i=0;i<9;i++)
    {
        eth->F1 [i]=0.0; 
        eth->F2 [i]=0.0;
        eth->F3 [i]=0.0;
        eth->F4 [i]=0.0;
        eth->Fav[i]=0.0;
        eth->Frr[i]=0.0;        
        eth->Frp[i]=0.0;
    }
    eth->Frp[1]=1.0/Rmh;               
    eth->Frp[3]=secB/Rnh;               
    eth->Frp[8]=1.0;

    /* refer psins */
    eth->g=eth->g0*(1+eth->beta[0]*sinB2-eth->beta[1]*sin2B2)-eth->beta[2]*pos[2];
    eth->gn[0]=0.0; eth->gn[1]=0.0; eth->gn[2]=-eth->g;

    vnadd(3,eth->wnie,2.0,eth->wnen,1.0,eth->wnien);
    vskewmv(1.0,eth->wnien,vel,temp1);
    /* harmful acceleration */
    vnadd(3,eth->gn,1.0,temp1,-1.0,eth->gcc);       

    /* state transition matrix related parameters */
    eth->F1 [3]=-eth->wie*sinB;      eth->F1 [6]=eth->wie*cosB;
    eth->F2 [2]=vel[1]/Rmh2;         eth->F2 [5]=-vel[0]/Rnh2;           eth->F2 [6]=vel[0]*secB2/Rnh; eth->F2[8]=-vel[0]*tanB/Rnh2;
    eth->F3 [3]=-2.0*eth->beta[3]*pos[2]*cos2B;                          
    eth->F3 [5]=-eth->beta[3]*sin2B; 
    eth->F3 [6]=-eth->g0*sin2B*(eth->beta[0]-4.0*eth->beta[1]*cos2B);    eth->F3 [8]=eth->beta[2];
    eth->F4 [0]=-eth->g0/Rnh;        eth->F4 [4]=-eth->g0/Rmh;           eth->F4 [8]=-2.0*eth->g0/(sqrt(Rmh*Rnh)+pos[2]);
    eth->Fav[1]=-1.0/Rmh;            eth->Fav[3]=1.0/Rnh;                eth->Fav[6]=tanB/Rnh;    
    eth->Frr[0]=vel[2]/Rnh-vel[1]*tanB/Rmh;                              eth->Frr[1]=vel[0]*tanB/Rmh;
    eth->Frr[2]=-vel[0]/Rnh;         eth->Frr[4]=vel[2]/Rmh;             eth->Frr[5]=-vel[1]/Rmh;
}

/* initialize ins related parameters -----------------------------------------------*/
extern int ins_init(ins_t *ins, const prcopt_t *popt)
{
    double install_angle[3]={0.0};
    int i,nx;

    ins->nx=GINS_NX;nx=ins->nx;
    ins->F=zeros(nx,nx); ins->Phi=zeros(nx,nx);
    ins->G=zeros(nx,nx); ins->Q=zeros(nx,nx);

    ins->time.sec=0.0; ins->time.time=0.0;
    ins->max_outime=popt->max_outime;
    ins->interval=1.0/popt->insample;
    ins->nn=popt->nn;
    ins->dttol=ins->interval/1e3;
    ins->discretime=(popt->insample%10)==0?1e-1:((popt->insample%25)==0?25*ins->interval:ins->interval);
    
    /* accelerometer and gyroscope velocity random walk and angle random walk and bias drive noise */
    ins->corr_time=popt->corr_time; 
    if (IMU_BIAS_MANUAL==popt->init_bias_type) {
        matcpy(ins->init_gyro_bias,popt->init_gyro_bias,3,1); /* rad/s */
        matcpy(ins->init_acce_bias,popt->init_acce_bias,3,1); /* mg */

        for (i=0;i<3;i++) ins->init_acce_bias[i]*=1e-3*G0; /* convert mg to m/s^2 */
        ins->bias_flag=1;        
    }

    ins->psd_gyro=popt->psd_gyro;
    ins->psd_acce=popt->psd_acce;
    ins->psd_bg=popt->psd_bg;
    ins->psd_ba=popt->psd_ba;

    ins->odo=NULL;

    /* init zupt configuration options (sliding window length and detection threshold) */
    ins->zupt.window=popt->insample;
    ins->zupt.gthres=popt->zupt_gthres;

    /* init process noise covariance matrix */
    for (i=0;i<GINS_NX;i++){
        ins->xa[i]=0.0;
        
        if (i<3)             ins->Q[i+i*nx]=ins->psd_gyro*ins->discretime;
        else if(i>=3&&i<6)   ins->Q[i+i*nx]=ins->psd_acce*ins->discretime;
        else if(i>=9&&i<12)  ins->Q[i+i*nx]=ins->psd_bg*ins->discretime;
        else if(i>=12&&i<15) ins->Q[i+i*nx]=ins->psd_ba*ins->discretime;
    }
    /* trace(12,"Q=\n"); tracemat(12,ins->Q,nx,nx,20,16); */ /*ok*/

    for (i=0;i<3;i++){
        /* initialize the lever and motion constraint information */
        install_angle[i]=popt->install_angle[i];
        ins->lever_nhc[i]=popt->lever_nhc[i];
        ins->lever[i]=popt->lever[i];

        /* INS basic parameters */
        ins->dw[i]=0.0;
        ins->dv[i]=0.0;
        ins->p1dv[i]=0.0;
        ins->p1dw[i]=0.0;
        ins->p1vel[i]=0.0;
        ins->p2vel[i]=0.0;
        ins->bg[i]=0.0;
        ins->ba[i]=0.0;
        ins->wbib[i]=0.0;
        ins->fb[i]=0.0;
        ins->fn[i]=0.0;

        ins->p1pos[i]=0.0;
        ins->pos[i]=0.0;
        ins->vel[i]=0.0;
        ins->body_vel[i]=0.0;
        ins->nhc_vel[i]=0.0;
        ins->att[i]=0.0;       
    }

    /* initialize the posture matrix and installation angle matrix */
    att2qnb(ins->att,ins->qnb);
    att2Cnb(ins->att,ins->Cnb);
    att2Cnb(install_angle,ins->Cvb);
}

/* update INS previous related parameters  --------------------------
*
*args   : ins_t *ins     I    ins struct
*return : none
*-------------------------------------------------------------------------------*/
extern void update_ins(ins_t *ins)
{
    int i;

    for (i=0;i<3;i++)
    {
        ins->p1dw[i]=ins->dw[i];
        ins->p1dv[i]=ins->dv[i];
        ins->p2vel[i]=ins->p1vel[i]; 
        /* update previous epoch pos/vel by ins predicted information */       
        ins->p1vel[i]=ins->vel[i];
        ins->p1pos[i]=ins->pos[i];
    }
}

/* rotation error, padding error and second-order cone error compensation */
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
    int i;

    /* if init bias flag is false, then set it to zero */
    if (!ins->bias_flag) {
        for (i=0;i<3;i++) ins->init_gyro_bias[i]=0.0;
        for (i=0;i<3;i++) ins->init_acce_bias[i]=0.0;
    }

    /* correct gyroscope and accelerometer zero bias (includes initial bias and estimated residual bias ) */
    for (i=0;i<3;i++) {
        ins->dw[i]=imu[0].dw[i]-(ins->init_gyro_bias[i]+ins->bg[i])*ins->interval;
        ins->dv[i]=imu[0].dv[i]-(ins->init_acce_bias[i]+ins->ba[i])*ins->interval;
    }

    if (2==ins->nn){
        for (i=0;i<3;i++) {
            ins->n1dw[i]=imu[1].dw[i]-(ins->init_gyro_bias[i]+ins->bg[i])*ins->interval;
            ins->n1dv[i]=imu[1].dv[i]-(ins->init_acce_bias[i]+ins->ba[i])*ins->interval;            
        }    
    }
}

/* ins mechanization -----------------------------------------------*/
extern void ins_mech(ins_t *ins, imud_t *imu, const prcopt_t *popt) 
{
    int i;
    double vel_m[3],pos_m[3],temp1[3],temp2[3],dv_rot[3],dv_pad[3],da_con[3],interval,*I3;
    double rv[3],rs[9],Irs[9],Ce[9],dw[3],dv[3],delta_vb[3],delta_vn[3];
    double Cnnk[9],Cnkn[9],Cbbk[9],Cnb[9],Cbn[9]={0.0};
    double qnb[4],qnkn[4],qbbk[4];

    I3=eye(3);

    /* sample interval of imu */
    interval=ins->interval*ins->nn;

    /* bias correction for gyroscopes and accelerometers */
    imu_fedback(ins,imu);

    /* initialize the attitude matrix and quaternion */
    for (i=0;i<9;i++) Cnb[i]=ins->Cnb[i];
    for (i=0;i<4;i++) qnb[i]=ins->qnb[i];

    /* extrapolate velocity and position at k-1/2 */
    vnadd(3,ins->p1vel,3.0/2.0,ins->p2vel,-1.0/2.0,vel_m);
    vnmul(3,vel_m,interval/2.0,temp1);

    Mat3mulv(1.0,ins->eth.Frp,temp1,temp2);
    vnadd(3,ins->p1pos,1.0,temp2,1.0,pos_m);

    /* update earth parameters at k-1/2 */
    earth_update(popt,pos_m,vel_m,&ins->eth);

    /* coning error, rotation error and paddling error compensation */
    conpad_fedback(ins,dw,dv,da_con,dv_rot,dv_pad);

    vnmul(3,dw,1.0/interval,ins->wbib);
    vnmul(3,dv,1.0/interval,ins->fb);

    if (ATT_DCM==popt->att_type) Mat3mulv(1.0,Cnb,ins->fb,ins->fn);        
    else quatmulv(1.0,qnb,ins->fb,ins->fn);

    /* forward mechanization, sequentially performing velocity update, position update, and attitude update */
    if (MECH_FORWARD==popt->reverse) {
        /* [I-1/2*(zetax)*Cnb(k-1)] of velocity update */
        vnmul(3,ins->eth.wnin,interval,rv);
        vskew(1.0,rv,rs);  
        Mat3add2(I3,1.0,rs,-1.0/2.0,Irs);
        if (ATT_QUAT==popt->att_type) qnb2Cnb(qnb,Cnb);  
        Mat3mul2(1.0,Irs,Cnb,Ce);

        /* rotation and paddling error compensation items of velocity update */
        for (i=0;i<3;i++){
            delta_vb[i]=dv[i]+dv_rot[i]+dv_pad[i];
        }
        Mat3mulv(1.0,Ce,delta_vb,delta_vn);

        /* NOTE: velocity update */
        for (i=0;i<3;i++) {
            ins->vel[i]=ins->p1vel[i]+(ins->eth.gcc[i]*interval)+delta_vn[i];
        }

        /* NOTE: position update */
        Mat3mulv(1.0/2.0*interval,ins->eth.Frp,ins->p1vel,temp1);
        Mat3mulv(1.0/2.0*interval,ins->eth.Frp,ins->vel,temp2);
        for (i=0;i<3;i++) {
            ins->pos[i]=ins->p1pos[i]+(temp1[i]+temp2[i]);
        }

        /* NOTE: attitude update */
        if (ATT_DCM==popt->att_type) {
            rv2DCM(interval,ins->eth.wnin,Cnnk);
            DCMT(Cnnk,Cnkn);
            rv2DCM(1.0,da_con,Cbbk);
            Mat3mul3(Cnkn,Cnb,Cbbk,ins->Cnb);
            Cnb2att(ins->Cnb,ins->att);  
            att2qnb(ins->att,ins->qnb);      
        }
        else {
            rv2quat(-1.0*interval,ins->eth.wnin,qnkn);
            rv2quat(1.0,da_con,qbbk);
            quatmul3(qnkn,qnb,qbbk,ins->qnb);
            qnbnorm(ins->qnb);
            qnb2Cnb(ins->qnb,ins->Cnb);
            Cnb2att(ins->Cnb,ins->att);
        }        
    }
    /* backward mechanization, sequentially performing attitude update, velocity update, and position update */
    else if (MECH_BACKWARD==popt->reverse) {
        /* NOTE: attitude update */
        if (ATT_DCM==popt->att_type) {
            rv2DCM(interval,ins->eth.wnin,Cnnk);
            DCMT(Cnnk,Cnkn);
            rv2DCM(1.0,da_con,Cbbk);
            Mat3mul3(Cnkn,Cnb,Cbbk,ins->Cnb);
            Cnb2att(ins->Cnb,ins->att);  
            att2qnb(ins->att,ins->qnb);      
        }
        else {
            rv2quat(-1.0*interval,ins->eth.wnin,qnkn);
            rv2quat(1.0,da_con,qbbk);
            quatmul3(qnkn,qnb,qbbk,ins->qnb);
            qnbnorm(ins->qnb);
            qnb2Cnb(ins->qnb,ins->Cnb);
            Cnb2att(ins->Cnb,ins->att);
        } 

        /* [I-1/2*(zetax)*Cnb(k-1)] of velocity update */
        vnmul(3,ins->eth.wnin,interval,rv);
        vskew(1.0,rv,rs);  
        Mat3add2(I3,1.0,rs,-1.0/2.0,Irs);
        if (ATT_QUAT==popt->att_type) qnb2Cnb(qnb,Cnb);  
        Mat3mul2(1.0,Irs,Cnb,Ce);

        /* rotation and paddling error compensation items of velocity update */
        for (i=0;i<3;i++)
        {
            delta_vb[i]=dv[i]+dv_rot[i]+dv_pad[i];
        }
        Mat3mulv(1.0,Ce,delta_vb,delta_vn);

        /* NOTE: velocity update */
        for (i=0;i<3;i++) {
            ins->vel[i]=ins->p1vel[i]+(ins->eth.gcc[i]*interval)+delta_vn[i];
        }

        /* NOTE: position update */
        Mat3mulv(1.0/2.0*interval,ins->eth.Frp,ins->p1vel,temp1);
        Mat3mulv(1.0/2.0*interval,ins->eth.Frp,ins->vel,temp2);
        for (i=0;i<3;i++) {
            ins->pos[i]=ins->p1pos[i]+(temp1[i]+temp2[i]);
        }
    }

    /* update ins velocty in b frame*/
    DCMT(ins->Cnb,Cbn);
    Mat3mulv(1.0,Cbn,ins->vel,ins->body_vel);

    /* update INS previous related parameters  */
    update_ins(ins);

    free(I3);
}

/* get non-zero element index */
static int get_nozeroidx(const double *x, const double *P, int nx, int offset, int *ix)
{
    int i,k=0;

    for (i=offset;i<nx;i++) {
        if (x[i]!=0.0&&P[i+i*nx]!=0.0) ix[k++]=i;
    }

    return k;
}

/* update cross-covariance -------------------------------------------
*args  :  rtk_t    *rtk   IO   rtk structure
*return:none
*-----------------------------------------------------------------------------*/
extern void update_crosscov(rtk_t *rtk)
{
    int i,j,k,ns=rtk->ins.nx,nx=rtk->nx,*ix=NULL;
    double *P_IG=NULL,*P_IG_=NULL,*P_GI_=NULL;

#if 1
    ix=imat(nx-ns,1);
    k=get_nozeroidx(rtk->x,rtk->P,nx,ns,ix);

    /* no non-zero element */
    if (k<=0) {
        free(ix); return;    
    }

    P_IG=zeros(ns,k); P_IG_=zeros(ns,k);

    for (i=0;i<ns;i++) for (j=0;j<k;j++) P_IG[j+i*k]=rtk->P[ix[j]+i*nx];

    matmul("NN",ns,ns,k,rtk->ins.Phi,P_IG,P_IG_,1.0,0.0);

    for (i=0;i<ns;i++) for (j=0;j<k;j++) rtk->P[ix[j]+i*nx]=P_IG_[j+i*k]; /* P_IG */
    for (i=0;i<k;i++) for (j=0;j<ns;j++) rtk->P[j+ix[i]*nx]=P_IG_[i+j*k]; /* P_GI */

    /* trace(12,"P_IG(k)=\n"); tracemat(12,P_IG_,ns,k,13,6); */

    free(ix); free(P_IG); free(P_IG_);
#else
    P_IG=zeros(ns,nx-ns);P_IG_=zeros(ns,nx-ns);P_GI_=zeros(ns,nx-ns);

    pmatcpy(P_IG,ns,nx-ns,0,0,ns,nx-ns,rtk->P,nx,nx,0,ns,ns,nx);
    /* trace(12,"P_IG(k-1)=\n"); tracemat(12,P_IG,ns,nx-ns,13,6); */

    matmul("NN",ns,ns,nx-ns,rtk->ins.Phi,P_IG,P_IG_,1.0,0.0);
    /* MatirxT(P_IG_,P_GI_,ns,nx-ns); */
    matmul("TT",nx-ns,ns,ns,P_IG,rtk->ins.Phi,P_GI_,1.0,0.0);

    pmatcpy(rtk->P,nx,nx,0,ns,ns,nx,P_IG_,ns,nx-ns,0,0,ns,nx-ns);
    pmatcpy(rtk->P,nx,nx,ns,0,nx,ns,P_GI_,nx-ns,ns,0,0,nx-ns,ns);

    free(P_IG);free(P_IG_);free(P_GI_);
#endif
}

/* update INS state transition matrix F and noise driving matrix G */
extern void phi_update(ins_t *ins, const prcopt_t *popt)
{
    int i,j,nx,k;
    double Faa[9],Fap[9],Far[9],Fva[9],Fvv[9],Fvp[9],Fvr[9];
    double Fvv1[9],Fvv2[9],F12[9],Fvp1[9],*Fg,*I,*I3;
    double Frr[9];

    nx=ins->nx;
    Fg=zeros(3,3);I=eye(nx);I3=eye(3);

    /* initialize F and G */
    for (i=0;i<nx;i++){
        for (j=0;j<nx;j++){
            ins->F[j+i*nx]=0.0;
            ins->G[j+i*nx]=0.0;
            if (i==j) ins->Phi[j+i*nx]=1.0; else ins->Phi[j+i*nx]=0.0;
        }
    }

    /* update earth parameters */
    earth_update(popt,ins->pos,ins->vel,&ins->eth);

    /* NOTE: Phi angle error model */
    if (ERR_PHI==popt->err_model) {
        vskew(-1.0,ins->eth.wnin,Faa);
        Mat3add2(ins->eth.F1,1.0,ins->eth.F2,1.0,Fap);
        Mat3mul2(1.0,Fap,ins->eth.Frp,Far);

        vskew(1.0,ins->fn,Fva);
        vskewmat3(1.0,ins->vel,ins->eth.Fav,Fvv1);
        vskew(-1.0,ins->eth.wnien,Fvv2);
        Mat3add2(Fvv1,1.0,Fvv2,1.0,Fvv);
        Mat3add2(ins->eth.F1,2.0,ins->eth.F2,1.0,F12);
        vskewmat3(1.0,ins->vel,F12,Fvp1);
        Mat3add2(Fvp1,1.0,ins->eth.F3,1.0,Fvp);
        Mat3mul2(1.0,Fvp,ins->eth.Frp,Fvr);

        if (ins->corr_time>0)
        {
            Fg[0]=-1.0/ins->corr_time;Fg[4]=-1.0/ins->corr_time;Fg[8]=-1.0/ins->corr_time; 
        }
        
        for (i=0;i<nx;i++)
        {
            if (i<3) {
                k=i;
                for (j=0;j<3;j++)
                {
                    ins->F[j+i*nx]=Faa[j+k*3];
                    ins->G[j+i*nx]=-ins->Cnb[j+k*3];
                }
                for (j=3;j<6;j++)
                {
                    ins->F[j+i*nx]=ins->eth.Fav[j-3+k*3];
                }
                for (j=6;j<9;j++)
                {
                    ins->F[j+i*nx]=Far[j-6+k*3];
                }
                for (j=9;j<12;j++)
                {
                    ins->F[j+i*nx]=-ins->Cnb[j-9+k*3];
                } 
            }     
            if (i>=3&&i<6) {
                k=i-3;
                for (j=0;j<3;j++)
                {
                    ins->F[j+i*nx]=Fva[j+k*3];
                }
                for (j=3;j<6;j++)
                {
                    ins->F[j+i*nx]=Fvv[j-3+k*3];
                    ins->G[j+i*nx]=ins->Cnb[j-3+k*3];
                }
                for (j=6;j<9;j++)
                {
                    ins->F[j+i*nx]=Fvr[j-6+k*3];
                }
                for (j=12;j<15;j++)
                {
                    ins->F[j+i*nx]=ins->Cnb[j-12+k*3];
                } 
            } 
            if (i>=6&&i<9) {
                k=i-6;
                for (j=3;j<6;j++)
                {
                    ins->F[j+i*nx]=I3[j-3+k*3];
                }
                for (j=6;j<9;j++)
                {
                    ins->F[j+i*nx]=ins->eth.Frr[j-6+k*3];
                }
            } 
            if (i>=9&&i<12) {
                k=i-9;
                for (j=9;j<12;j++)
                {
                    ins->F[j+i*nx]=Fg[j-9+k*3];
                    if (i==j) ins->G[j+i*nx]=1.0;
                } 
            } 
            if (i>=12&&i<15) {
                k=i-12;
                for (j=12;j<15;j++)
                {
                    ins->F[j+i*nx]=Fg[j-12+k*3];
                    if (i==j) ins->G[j+i*nx]=1.0;
                } 
            }                                                   
        }        
    }
    /* NOTE: Psi angle error model */
    else if (ERR_PSI==popt->err_model) {
        vskew(-1.0,ins->eth.wnin,Faa); /* att */

        vskew(1.0,ins->fn,Fva);
        vskew(-1.0,ins->eth.wnien,Fvv2); /* vel */
       
        vskew(-1.0,ins->eth.wnen,Frr);   /* pos */

        if (ins->corr_time>0)
        {
            Fg[0]=-1.0/ins->corr_time;Fg[4]=-1.0/ins->corr_time;Fg[8]=-1.0/ins->corr_time; 
        }
        
        for (i=0;i<nx;i++)
        {
            if (i<3) {
                k=i;
                for (j=0;j<3;j++)
                {
                    ins->F[j+i*nx]=Faa[j+k*3];
                    ins->G[j+i*nx]=-ins->Cnb[j+k*3];
                }
                for (j=9;j<12;j++)
                {
                    ins->F[j+i*nx]=-ins->Cnb[j-9+k*3];
                } 
            }     
            if (i>=3&&i<6) {
                k=i-3;
                for (j=0;j<3;j++)
                {
                    ins->F[j+i*nx]=Fva[j+k*3];
                }
                for (j=3;j<6;j++)
                {
                    ins->F[j+i*nx]=Fvv2[j-3+k*3];
                    ins->G[j+i*nx]=ins->Cnb[j-3+k*3];
                }
                for (j=6;j<9;j++)
                {
                    ins->F[j+i*nx]=ins->eth.F4[j-6+k*3];
                }
                for (j=12;j<15;j++)
                {
                    ins->F[j+i*nx]=ins->Cnb[j-12+k*3];
                } 
            } 
            if (i>=6&&i<9) {
                k=i-6;
                for (j=3;j<6;j++)
                {
                    ins->F[j+i*nx]=I3[j-3+k*3];
                }
                for (j=6;j<9;j++)
                {
                    ins->F[j+i*nx]=Frr[j-6+k*3];
                }
            } 
            if (i>=9&&i<12) {
                k=i-9;
                for (j=9;j<12;j++)
                {
                    ins->F[j+i*nx]=Fg[j-9+k*3];
                    if (i==j) ins->G[j+i*nx]=1.0;
                } 
            } 
            if (i>=12&&i<15) {
                k=i-12;
                for (j=12;j<15;j++)
                {
                    ins->F[j+i*nx]=Fg[j-12+k*3];
                    if (i==j) ins->G[j+i*nx]=1.0;
                } 
            }                                                   
        }                
    }

    /* discretization of the state transition matrix Phi */
    matmul("NN",nx,nx,nx,ins->F,I,ins->Phi,ins->discretime,1.0);

    /* trace(12,"F=\n"); tracemat(12,ins->F,nx,nx,20,16);
    trace(12,"Phi=\n"); tracemat(12,ins->Phi,nx,nx,20,16);
    trace(12,"G=\n"); tracemat(12,ins->G,nx,nx,9,4); */ /*ok*/

    free(Fg);free(I);free(I3);
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

    if (GINS_TC==rtk->opt.GI_mode) pmatcpy(P,nx,nx,0,0,nx,nx,rtk->P,rtk->nx,rtk->nx,0,0,nx,nx);
    else matcpy(P,rtk->lcgins.P,nx,nx);

    /* if (GINS_TC==rtk->opt.GI_mode) trace(12,"Pk-1=\n"); tracemat(12,rtk->P,rtk->nx,rtk->nx,9,2); */
    /* trace(12,"Pk-1=\n"); tracemat(12,P,nx,nx,9,2); */

    /* time update */
    matmul("NN",nx,nx,nx,ins->Phi,P,FP,1.0,0.0);          /* FP=F*P */
    matmul("NT",nx,nx,nx,FP,ins->Phi,P,1.0,0.0);          /* FPF=FP*F' */

    matmul("NN",nx,nx,nx,ins->G,ins->Q,GQ,1.0,0.0);      /* GQ=G*Q */
    matmul("NT",nx,nx,nx,GQ,ins->G,P,1.0,1.0);           /* FPF=F*P*F'+G*Q*G' */
   
    if (GINS_TC==rtk->opt.GI_mode) pmatcpy(rtk->P,rtk->nx,rtk->nx,0,0,nx,nx,P,nx,nx,0,0,nx,nx);
    else matcpy(rtk->lcgins.P,P,nx,nx);

    /* NOTE: update GNSS/INS cross-covariance!!! */
    if (GINS_TC==rtk->opt.GI_mode&&rtk->outage>0) init_crosscov(rtk,nx,rtk->nx);
    if (GINS_TC==rtk->opt.GI_mode&&rtk->outage==0) update_crosscov(rtk);

    /* if (GINS_TC==rtk->opt.GI_mode) trace(12,"P_pre=\n"); tracemat(12,rtk->P,rtk->nx,rtk->nx,9,2); */
    /* trace(12,"P_pre=\n"); tracemat(12,P,nx,nx,9,2); */

    free(FP); free(GQ); free(P);
}

/* for backward processing mode, the sign of the INS velocity and gyroscope bias is inverted */
extern void pos_reverse(const prcopt_t *popt, ins_t *ins, int *reverse_flag)
{
    int i;

    if (SOLTYPE_COMBINED<=popt->soltype&&SOLTYPE_BACKWARD==popt->reverse) {
        for (i=0;i<3;i++) {
            ins->vel[i]=-ins->vel[i];
            ins->bg[i]=-ins->bg[i];
        }
        *reverse_flag=1;
    }
}