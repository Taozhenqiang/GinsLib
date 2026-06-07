/*----------------------
 * INS/GNSS integration common functions
 *----------------------*/

#include "rtklib.h"

/* integrated navigation initialization */
extern void gins_init(rtk_t *rtk, const prcopt_t *popt)
{
    ins_t *ins=&rtk->ins;
    int i,nx;
    double P[GINS_NX]={0.0};

    /* init INS struct */
    ins_init(ins,popt);

    /* loose coupled and semi-tight coupled mode */
    if (GINS_LC==popt->GI_mode||GINS_STC==popt->GI_mode) {
        sol_t sol0={{0}};
        nx=GINS_NX; 

        rtk->lcgins.x=zeros(nx,1);
        rtk->lcgins.P=zeros(nx,nx);
        rtk->lcgins.nx=nx;
        rtk->lcgins.sol=sol0;

        for (i=0;i<nx;i++){
            if (i<3)              P[i]=popt->init_att_unc[i]*popt->init_att_unc[i];
            else if (i>=3&&i<6)   P[i]=popt->init_vel_unc[i-3]*popt->init_vel_unc[i-3];
            else if (i>=6&&i<9)   P[i]=popt->init_pos_unc[i-6]*popt->init_pos_unc[i-6];
            else if (i>=9&&i<12) {
                if (IMU_BIAS_MANUAL==popt->init_biasunc_type) P[i]=popt->init_bg_unc*popt->init_bg_unc;
                else  P[i]=(sqrt(ins->psd_bg)*1E2)*(sqrt(ins->psd_bg)*1E2); /* auto */
            } 
            else if (i>=12&&i<=15) {
                if (IMU_BIAS_MANUAL==popt->init_biasunc_type) P[i]=popt->init_ba_unc*popt->init_ba_unc;    
                else  P[i]=(sqrt(ins->psd_ba)*1E2)*(sqrt(ins->psd_ba)*1E2); /* auto */
            }            
        }
        for (i=0;i<nx;i++) rtk->lcgins.P[i+i*nx]=P[i];
        /* trace(12,"P=\n"); tracemat(12,rtk->lcgins.P,nx,nx,9,4); */
    } 
    /* tight coupled mode */
    if (GINS_TC==popt->GI_mode) {
        nx=GINS_NX; 
        for (i=0;i<nx;i++){
            if (i<3)              P[i]=popt->init_att_unc[i]*popt->init_att_unc[i];
            else if (i>=3&&i<6)   P[i]=popt->init_vel_unc[i-3]*popt->init_vel_unc[i-3];
            else if (i>=6&&i<9)   P[i]=popt->init_pos_unc[i-6]*popt->init_pos_unc[i-6];
            else if (i>=9&&i<12) {
                if (IMU_BIAS_MANUAL==popt->init_biasunc_type)  P[i]=popt->init_bg_unc*popt->init_bg_unc;
                else  P[i]=(sqrt(ins->psd_bg)*1E2)*(sqrt(ins->psd_bg)*1E2); /* auto */
            }
            else if (i>=12&&i<=15) {
                if (IMU_BIAS_MANUAL==popt->init_biasunc_type) P[i]=popt->init_ba_unc*popt->init_ba_unc;
                else  P[i]=(sqrt(ins->psd_ba)*1E2)*(sqrt(ins->psd_ba)*1E2); /* auto */
            }       
        }
        for (i=0;i<nx;i++) rtk->P[i+i*rtk->nx]=P[i];
    }     
}

/* Convert INS solutions to GNSS center */
extern void ins2gnss(const prcopt_t *popt, ins_t *ins, double *pv_g, int n)
{
    int i;
    double F1[9],lever_n[3],Cbn[9],wbie[3],wbeb[3],temp[3],d_v[3],pv[6];

    Mat3mul2(1.0,ins->eth.Frp,ins->Cnb,F1);
    Mat3mulv(1.0,F1,ins->lever,lever_n);
    vnadd(3,ins->pos,1.0,lever_n,1.0,pv);       

    if (n==6) {
        DCMT(ins->Cnb,Cbn);
        Mat3mulv(1.0,Cbn,ins->eth.wnie,wbie);
        Mat3add2(ins->wbib,1.0,wbie,-1.0,wbeb);
        vskewmv(1.0,wbeb,ins->lever,temp);
        Mat3mulv(1.0,ins->Cnb,temp,d_v);
        vnadd(3,ins->vel,1.0,d_v,1.0,pv+3);        
    }

    for (i=0;i<n;i++)
    {
        pv_g[i]=pv[i];
    }
}

extern void insfix2gnss(ins_t *ins, double *pv_g, int n)
{
    int i;
    double F1[9],lever_n[3],Cbn[9],wbie[3],wbeb[3],temp[3],d_v[3],pv[6];

    Mat3mul2(1.0,ins->eth.Frp,ins->Cnb,F1);
    Mat3mulv(1.0,F1,ins->lever,lever_n);
    vnadd(3,ins->xa+6,1.0,lever_n,1.0,pv);       

    if (n==6) {
        DCMT(ins->Cnb,Cbn);
        Mat3mulv(1.0,Cbn,ins->eth.wnie,wbie);
        Mat3add2(ins->wbib,1.0,wbie,-1.0,wbeb);
        vskewmv(1.0,wbeb,ins->lever,temp);
        Mat3mulv(1.0,ins->Cnb,temp,d_v);
        vnadd(3,ins->xa+3,1.0,d_v,1.0,pv+3);        
    }

    for (i=0;i<n;i++)
    {
        pv_g[i]=pv[i];
    }
}

/* convert GNSS solutions to INS center */
extern void gnss2ins(rtk_t *rtk, double *pv_g, double *pv_i, int mode)
{
    ins_t *ins=&rtk->ins;
    int i;
    double F1[9],lever_n[3],Cbn[9],wbie[3],wbeb[3],temp[3],d_v[3];

    /* convert GNSS position to INS position */
    if (mode==1) {
        /* r_gnss=r_ins+Fpv*Cnb*lever_b */
        Mat3mul2(1.0,ins->eth.Frp,ins->Cnb,F1);
        Mat3mulv(1.0,F1,ins->lever,lever_n);

        vnadd(3,pv_g,1.0,lever_n,-1.0,pv_i);
    }
    /* convert GNSS velocity to INS velocity */
    else if (mode==2) {
        /* vn_gnss=vn_ins+Cnb([wbeb x]lever_b)*/
        DCMT(ins->Cnb,Cbn);
        Mat3mulv(1.0,Cbn,ins->eth.wnie,wbie);
        Mat3add2(ins->wbib,1.0,wbie,-1.0,wbeb);
        vskewmv(1.0,wbeb,ins->lever,temp);
        Mat3mulv(1.0,ins->Cnb,temp,d_v);

        vnadd(3,pv_g,1.0,d_v,-1.0,pv_i);      
    }
}

/* convert psi error state to phi error state ------------------------------
*args   : ins_t *ins           I   ins structure
*         const double *dr     I   dblh pos error  (3x1)
*         double *dx           IO  ins error state (15x1)
*return : none
*-------------------------------------------------------------------------------*/
extern void psi2phi_corr(ins_t *ins, const double *dr, double *dx)
{
    double d_ceta[3],dv[3];
    int i;

    /* equivalent rotation vector phi_nc */
    d_ceta[0]=-dr[0]; d_ceta[1]=dr[1]*cos(ins->pos[0]); d_ceta[2]=dr[1]*sin(ins->pos[0]);

    /* NOTE: convert Psi attitude misalignment angle (psi_cn') to Phi attitude misalignment angle (phi_nn') , phi_nn'=phi_nc+psi_cn' */     
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

    /* convert denu to dblh */
    earth_update(popt,ins->pos,ins->vel,&ins->eth);
    Mat3mulv(1.0,ins->eth.Frp,dx+6,dr);

    /* NOTE: convert psi error state to phi error state */
    if (ERR_PSI==popt->err_model) psi2phi_corr(ins,dr,dx);

    /* quaternion-based attitude feedback correction, qnb=qnn_°qn_b */
    for (i=0;i<4;i++) qn_b[i]=ins->qnb[i];
    for (i=0;i<3;i++) phi_nn_[i]=dx[i];
    
    rv2quat(1.0,phi_nn_,qnn_);
    quatmul(qnn_,qn_b,ins->qnb); 
    qnbnorm(ins->qnb); /* normalize qnb */
    qnb2Cnb(ins->qnb,ins->Cnb);
    Cnb2att(ins->Cnb,ins->att);

    /* DCM-based attitude feedback correction, Cnb=(I+[phi x])Cn'b */
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

    /* update previous epoch pos/vel by kf updated state */
    for (i=0;i<3;i++){
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

    /* convert denu to dblh */
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