/*-----------------------
* motion constraint model
* ---------------------*/

#include "rtklib.h"

static double init_gyro_bias[3]; /* init gyro bias (rad/s) */
static double init_acce_bias[3]; /* init acce bias (m/s^2) */

/* zero speed detection */
extern void zerovel_detect(rtk_t *rtk, imud_t *imu)
{
    ins_t *ins=&rtk->ins;
    zupt_t *zupt=&ins->zupt;
    int i,j,nn=ins->nn,window=(int)20/(ins->interval*nn); /* 20 s static imu data */
    double N=zupt->window,dw[3]={0.0},ndw;

    /* for multi-sample mode, the average value of the gyroscope is used */
    for (j=0;j<3;j++) for (i=0;i<nn;i++) dw[j]+=imu[i].dw[j];
    for (j=0;j<3;j++) dw[j]/=nn;
    ndw=norm(dw,3);

    /* try IMU static bias saving only before INS alignment is complete */
    if (!rtk->align) {
        if (zupt->count&&IMU_BIAS_AUTO==rtk->opt.init_bias_type&&!ins->bias_flag) {
            for (i=0;i<3;i++) for (j=0;j<nn;j++) {
                init_gyro_bias[i]+=imu[j].dw[i];
                init_acce_bias[i]+=imu[j].dv[i];
            }
            if (zupt->count>=window) {
                for (i=0;i<3;i++) {
                    ins->init_gyro_bias[i]=init_gyro_bias[i]/(window*nn)/ins->interval; /* rad/s */
                    ins->init_acce_bias[i]=init_acce_bias[i]/(window*nn)/ins->interval; /* m/s^2*/
                }
                ins->init_acce_bias[2]-=G0; /* remove gravity from z-axis accelerometer bias */
                /* NOTE: init bias of accelerometer not enabled */
                for (i=0;i<3;i++) ins->init_acce_bias[i]=0.0;
                ins->bias_flag=1;
            }
        }
        else {
            for (i=0;i<3;i++) init_gyro_bias[i]=init_acce_bias[i]=0.0;
        }        
    }
    
    /* initialization */
    if (!zupt->iimu) {
        zupt->Gm=ndw;
        zupt->Gd=0.0;
    }
    else {
        zupt->Gm=(N-1.0)/N*zupt->old_Gm+1.0/N*ndw;
        if (!zupt->Gd) zupt->Gd=fabs(zupt->Gm-ndw);
        else zupt->Gd=(N-1.0)/N*zupt->old_Gd+1.0/N*fabs(zupt->Gm-ndw);
    }
    zupt->iimu++;

    /* determine whether the vehicle is stationary */
    if (zupt->Gd&&zupt->Gd<zupt->gthres) zupt->count++;
    else zupt->count=0;

    /* save detection statistics for the current epoch */
    zupt->old_Gm=zupt->Gm;
    zupt->old_Gd=zupt->Gd;
}

/* motion measurement */
extern int motion_meas(rtk_t *rtk, const prcopt_t *popt, double *H, double *v, double *var, int nv, int nx)
{
    ins_t *ins=&rtk->ins;
    sol_t *sol=(GINS_TC==popt->GI_mode)?&rtk->sol:&rtk->lcgins.sol;
    double zupt_time=0.0,vel=0.0;
    int nv_cons=0;

    /* detected vehicle stationary time (s) and GNSS velocity */
    zupt_time=ins->zupt.count*ins->interval*ins->nn;
    vel=norm(rtk->sol.rr+3,3);

    /* NOTE: the vehicle is considered stationary only when the zero speed detection is passed, 
    the stationary state is greater than 1s and the calculated vehicle speed is less than 0.1m/s */
    if (popt->constraint[1]&&zupt_time>1.0&&(vel>0&&vel<0.1)) { /* zupt */
        nv_cons=motion_update(rtk,H,v,var,nv,nx,CONS_ZUPT);
        sol->iFlag=SOLF_ZUPT; /* zupt flag */
    }
    else if (popt->constraint[0]) { /* nhc */
        nv_cons=motion_update(rtk,H,v,var,nv,nx,CONS_NHC);        
    }
    if (popt->constraint[2]&&zupt_time>1.0&&(vel>0&&vel<0.1)) { /* zihr */
        nv_cons+=motion_update(rtk,H,v,var,nv+nv_cons,nx,CONS_ZIHR);
    }   
    
    return nv_cons;
}

/* motion constraints */
extern void motion_constraints(rtk_t *rtk, const prcopt_t *popt) 
{
    ins_t *ins=&rtk->ins;
    sol_t *sol=(GINS_TC==popt->GI_mode)?&rtk->sol:&rtk->lcgins.sol;
    int i,j,nx=ins->nx,nv=4,info; /* max nv:ZUPT+ZIHR=4 */
    double *xp,*Pp,*H,*v,*var,*R;

    /* initializing memory */
    xp=zeros(nx,1); Pp=zeros(nx,nx); R=zeros(nv,nv);
    H=mat(nv,nx); v=mat(nv,1); var=mat(nv,1);

    /* copy the covariance matrix */
    if (GINS_LC==popt->GI_mode||GINS_STC==popt->GI_mode) matcpy(Pp,rtk->lcgins.P,nx,nx);
    else if (GINS_TC==popt->GI_mode) pmatcpy(Pp,nx,nx,0,0,nx,nx,rtk->P,rtk->nx,rtk->nx,0,0,nx,nx); 

    nv=motion_meas(rtk,popt,H,v,var,0,nx);
    
    /* measurement noise covariance matrix */
    diag_Cov(nv,var,R,diag_var);

    /* measurement update */
    if ((info=filter_gins(rtk,xp,Pp,H,v,R,nx,nv,KF_GINS,Robust_OFF))) {
        trace(7,"motion_constraints: filter_gins error info=%d\n",info);
        sol->stat=SOLQ_INS;
        /* update solution status */
        update_instat(popt,ins,Pp,sol,nx);
        
        free(xp); free(Pp); free(H); free(v); free(var);
        return;
    }
    /* update solution status */
    sol->stat=SOLQ_CONS;

    /* update the covariance matrix */
    if (GINS_LC==popt->GI_mode||GINS_STC==popt->GI_mode) matcpy(rtk->lcgins.P,Pp,nx,nx);
    else if (GINS_TC==popt->GI_mode) pmatcpy(rtk->P,rtk->nx,rtk->nx,0,0,nx,nx,Pp,nx,nx,0,0,nx,nx); 
    
    /* ins feedback correction */
    ins_fedback(rtk,xp);

    /* update solution status */
    update_instat(popt,ins,Pp,sol,nx);

    free(xp); free(Pp); free(H); free(v); free(var);
}

/* zupt/nhc update */
extern int motion_update(rtk_t *rtk, double *H, double *v, double *var, int nv, int nx, int mode)
{
    ins_t *ins=&rtk->ins;
    int i,j,k,inv=0,k2[2]={0,2},k3[3]={0,1,2};
    int odo_flag=0;
    double Cbn[9]={0.0},Cvn[9]={0.0},Cne[9]={0.0},lever_v[9]={0.0},vel_v[9]={0.0},Ha_bg[3]={0.0},att[3]={0.0},yaw=0.0;
    double ref_veln[3]={0.0},ref_velv[3]={0.0};

    DCMT(ins->Cnb,Cbn);
    Mat3mul2(1.0,ins->Cvb,Cbn,Cvn);
    /* vehichle velocity of v frame */
    Mat3mulv(1.0,Cvn,ins->vel,ins->nhc_vel);
    Mat3mvskew(-1.0,Cvn,ins->vel,vel_v);
    Mat3mvskew(-1.0,ins->Cvb,ins->lever_nhc,lever_v);

    /* sim odo constraint */
    if (ODO_SIM==rtk->opt.odopt&&ins->odo) {
        xyz2enu(ins->pos,Cne);
        Mat3mulv(1.0,Cne,ins->odo->odo_vel,ref_veln);
        Mat3mulv(1.0,Cvn,ref_veln,ref_velv);

        /* ins forward velocity - odo forward velocity */
        ins->nhc_vel[1]-=ref_velv[1];
        odo_flag=1;

        trace(7,"motion_update: odo_flag=%d\n",odo_flag);
    }

    /* H of ZIHR */
    if (CONS_ZIHR==mode) {
        Ha_bg[0]=-sin(ins->att[1])/cos(ins->att[0])*rtk->interval; 
        Ha_bg[1]=0.0; 
        Ha_bg[2]=cos(ins->att[1])/cos(ins->att[0])*rtk->interval;        
    }

    /* converts the yaw from clockwise to counterclockwise */
    if (GINS_LC==rtk->opt.GI_mode) for (i=0;i<3;i++) att[i]=rtk->lcgins.sol.att[i];
    else for (i=0;i<3;i++) att[i]=rtk->sol.att[i];

    if (att[2]<=180) yaw=-att[2]*D2R;
    else yaw=(360.0-att[2])*D2R;

    /* determine constraint model */
    if (CONS_NHC==mode) {
        /* sim odo constraint */
        if (odo_flag) inv=3;
        else inv=2;
        trace(8,"nhc_constraints: v=\n");tracemat(8,ins->nhc_vel,3,1,9,4);
    }
    else if (CONS_ZUPT==mode) {
        inv=3;
        trace(8,"zupt_constraints: v=\n");tracemat(8,ins->nhc_vel,3,1,9,4);
    }
    else if (CONS_ZIHR==mode) {
        inv=1;
        trace(8,"zihr_constraints: yaw=\n");tracemat(8,ins->att,3,1,9,4);
    }

    for (i=0;i<inv;i++) {
        if (CONS_NHC==mode) if (odo_flag) k=k3[i]; else k=k2[i];
        else if (CONS_ZUPT==mode) k=k3[i];
        else if (CONS_ZIHR==mode) k=0;

        for (j=0;j<nx;j++) {
            /* if H is NULL, return 0 */
            if (H) H[j+nv*nx]=0.0;
            else return 0;
            
            /* update the measurement coefficient of NHC/ZUPT/ZIHR in H */
            if (CONS_NHC==mode||CONS_ZUPT==mode) {
                if (j<3)             H[j+nv*nx]=vel_v[j+k*3];
                else if (j>=3&&j<6)  H[j+nv*nx]=Cvn[(j-3)+k*3];
                else if (j>=9&&j<12) H[j+nv*nx]=lever_v[(j-9)+k*3];                
            }
            else if (CONS_ZIHR==mode) {
                if (j>=9&&j<12) H[j+nv*nx]=Ha_bg[j-9];
            }
        }

        if (CONS_NHC==mode||CONS_ZUPT==mode) {
            v[nv]=ins->nhc_vel[k];
            var[nv]=0.1; /* variance of the constraint, can be adjusted */  
        }
        else if (CONS_ZIHR==mode) {
            v[nv]=ins->att[2]-yaw;       /* the difference in yaw between the last GNSS update and the current INS update */
            var[nv]=(0.1*D2R)*(0.1*D2R); /* variance of the constraint (rad), can be adjusted */          
        }
        nv++;
    }

    return inv;
}