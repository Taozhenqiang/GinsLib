/*----------------------
 * ins init alignment 
-----------------------*/
#include "rtklib.h"

/* check if the GNSS status meets the INS alignment requirements */
static int gnss_aid_insalign(rtk_t *rtk)
{
    prcopt_t *popt=&rtk->opt;
    sol_t *sol=&rtk->sol;
    int flag=0;

    if (PMODE_SINGLE==popt->mode&&sol->ns>=6) flag=1; /* spp mode */
    else if (PMODE_DGPS==popt->mode&&sol->ns>=6) flag=1; /* DGPS mode */
    else if (PMODE_KINEMA<=popt->mode&&popt->mode<=PMODE_FIXED) { /* rtk mode */
        /* AR is fixed and ns is greater than 6 */
        if (popt->artype>OFF&&sol->stat==SOLQ_FIX&&sol->ns>=6) flag=1;
        else if (OFF==popt->artype&&sol->ns>=6) flag=1;
    }
    else if (PMODE_PPP_KINEMA<=popt->mode&&popt->mode<=PMODE_PPP_FIXED) { /* ppp mode */
        if (sol->ns>=6) flag=1;
    }

    return flag;

}

/* initialize INS position, velocity and attitude */
extern void init_inspva(ins_t *ins, const double *pos, const double *vel, const double *att) 
{
    int i;

    for (i=0;i<3;i++)
    {
        ins->p1pos[i]=pos[i];
        ins->pos[i]=pos[i];
        ins->p1vel[i]=vel[i];
        ins->p2vel[i]=vel[i];
        ins->vel[i]=vel[i];
        if (att) ins->att[i]=att[i];
    }
    if (att) {
        att2qnb(ins->att,ins->qnb);
        att2Cnb(ins->att,ins->Cnb);         
    }

    earth_init(ins->pos,ins->vel,&ins->eth);  
}

static void save_tdcp_att(rtk_t *rtk, const double *att)
{
    rtk->sol.pitch=att[0]*R2D;
    /* yaw (counterclockwise to clockwise, deg) */
    rtk->sol.yaw=att[2]*R2D;
    if (rtk->sol.yaw<0) rtk->sol.yaw*=-1;
    else rtk->sol.yaw=360.0-rtk->sol.yaw;
}

/* ins initial alignment -------------------------------------------*/
extern int ins_align(rtk_t *rtk, obsd_t *obs, int n, nav_t *nav, const prcopt_t *opt, int vel_flag)
{
    /* NOTE: The structure copy is a shallow copy! */
    prcopt_t popt=*opt;
    rtk_t  rtk_={0}; 
    ins_t *ins=&rtk->ins;
    int i,j,init_flag=0;
    double att[3]={0.0},pos[3]={0.0},vn[3]={0.0};

    popt.GI_mode=GINS_OFF;  /* set to GINS_OFF mode */
    /* NOTE: initialize rtk_ instead of assigning rtk to rtk_ to avoid shallow copying of the structure */
    if (SYNC_YES==rtk->upte) { 
        rtkinit(&rtk_,&popt,NULL); 
        init_flag=1; 
    }

    /* for GNSS/INS LC with pos file ,set alignment type to manual */
    if (PMODE_LC_POS==popt.mode) popt.alingetype=INSALI_MANUAL;

    /* manual alignment */
    if (!rtk->align&&INSALI_MANUAL==popt.alingetype&&popt.ts.time) {   
        if ((fabs(timediff(ins->time,popt.ts))-rtk->ins.dttol)<=rtk->ins.nn*rtk->ins.interval/2.0) {
            /* initialize ins position, velocity and attitude */
            init_inspva(ins,popt.initpos,popt.initvel,popt.initatt); 
            if (init_flag) rtkfree(&rtk_);  

            trace(12,"INS initial alignment completed: %s!\n",Debug_Glo.chTime); 
            showerr("INS initial alignment completed: %s!",Debug_Glo.chTime); 
            return 1;            
        }
        else if (timediff(ins->time,popt.ts)>0) {
            if (init_flag) rtkfree(&rtk_);  
            showmsg("warning : start time is smaller than GNSS/INS matching time!\n"); 
            return 0;
        }
    }

    /* velocity vector-assisted alignment and INS re-initialization require GNSS/INS time synchronization */
    if (SYNC_YES==rtk->upte&&vel_flag) {
        /* velocity vector assisted yaw initialization based on tdcp */
        if (!rtk->align&&INSALI_VELTOR==popt.alingetype)  {
            /* initialize INS position using GNSS solution */
            if (!rtkpos(&rtk_,obs,n,nav)||!gnss_aid_insalign(&rtk_)) {
                if (init_flag) rtkfree(&rtk_); trace(7,"rtkpos error: GNSS unavailable during INS align!\n");
                return 0;
            }
            /* initialize position (from GNSS) and velocity (from tdcp) */
            ecef2pos(rtk_.sol.rr,pos);        
            ecef2enu(pos,rtk->sol.rr+3,vn);
            
            /* initialize pitch and yaw using the velocity in the n frame */
            att[0]=atan2(vn[2],sqrt(vn[0]*vn[0]+vn[1]*vn[1])); /* pitch angle */
            att[2]=-atan2(vn[0],vn[1]);                        /* yaw angle */

            /* output the initial alignment solution status */
            save_tdcp_att(rtk,att);
            outsolstat(rtk,nav);

            /* initialize ins position, velocity and attitude ,consider lever arm correction */
            gnss2ins(rtk,pos,ins->pos,1);
            /* reverse the velocity vector if the solution type is backward */
            if (SOLTYPE_BACKWARD==popt.reverse) {
                for (i=0;i<3;i++) vn[i]=-vn[i];
            }
            gnss2ins(rtk,vn,ins->vel,2);
            init_inspva(ins,ins->pos,ins->vel,att);

            trace(12,"INS initial alignment completed: %s!\n",Debug_Glo.chTime); 
            showerr("INS initial alignment completed: %s!",Debug_Glo.chTime); 
            if (init_flag) rtkfree(&rtk_); /* NOTE: free the rtk_ structure!!! */

            return 1;           
        }

        /* reinitialize INS in the event of a long-term GNSS outage */
        if (rtk->align&&rtk->outage>ins->max_outime&&!outsim.valid_flag) {
            /* initialize INS position using GNSS solution */
            if (!rtkpos(&rtk_,obs,n,nav)||!gnss_aid_insalign(&rtk_)) {
                if (init_flag) rtkfree(&rtk_); trace(7,"rtkpos error: GNSS unavailable during INS reinitialization!\n");
                return rtk->align?1:0;
            }
            /* if GNSS becomes available after a long outage, reset the GNSS outage count */
            rtk->outage=0;

            /* reinitialize ins position and velocity, consider lever arm correction */
            ecef2pos(rtk_.sol.rr,pos);        
            ecef2enu(pos,rtk->sol.rr+3,vn);
            gnss2ins(rtk,pos,ins->pos,1);
            gnss2ins(rtk,vn,ins->vel,2);
            init_inspva(ins,ins->pos,ins->vel,NULL);

            trace(12,"INS reinitialization completed: %s!\n",Debug_Glo.chTime);
            showerr("INS reinitialization completed: %s!",Debug_Glo.chTime);   
            if (init_flag) rtkfree(&rtk_); /* NOTE: free the rtk_ structure!!! */
            return 1;     
        }  
        /* NOTE: free the rtk_ structure!!! */
        if (init_flag) rtkfree(&rtk_);   
        return rtk->align?1:0; 
    }
    else {
        /* NOTE: free the rtk_ structure!!! */
        if (init_flag) rtkfree(&rtk_);   
        return rtk->align?1:0; 
    }
}