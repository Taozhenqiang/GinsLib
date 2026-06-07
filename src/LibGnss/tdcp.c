/*--------------
 * TDCP module
 *-----------------------------*/

#include "rtklib.h"

/* TDCP-assisted motion alignment */
extern int tdcp_vel(rtk_t *rtk, rtk_t *rtk_main, const obsd_t *obs, const obsd_t *obs_old, int n, int n_old, const nav_t *nav, const prcopt_t *opt)
{
    prcopt_t opt_=*opt;
    sol_t sol={0},sol_old={0};
    double *rs,*rs_old,*dts,*dts_old,*vare,*vare_old,*resp,*resp_old,*azel,*azel_old;
    double rr[3],rr_old[3],r,dr[3],dr_old[3],er=0.0,er_old=0.0,e[3],e_old[3],freq,thres_ouj=2.0,thres=3.0;
    double sgn=(SOLTYPE_BACKWARD==opt->reverse?-1.0:1.0);
    double *v,*H,*var,*P,dx[4]={0},Q[4*4];
    int align=rtk_main->align;
    int sat[MAXSAT],ir_old[MAXSAT],ir[MAXSAT];
    int stat=0,stat_old=0,i,j,k,m,nf=rtk->opt.nf,sys,fr,nx=4,nv=0,vnv[MAXFREQ]={0},max_vnv=0,info,pos_flag=1,tdcp_flag=1,vel_flag=0,mode=Robust_RES;
    int vsat[MAXOBS]={0},vsat_old[MAXOBS]={0},svh[MAXOBS]={0},svh_old[MAXOBS]={0};

    /* check if there are valid observations */
    if (!n||!n_old) return 0;

    /* initialize rtk_tdcp parameters */
    rtk->interval=rtk_main->interval; rtk->dopsgn=rtk_main->dopsgn;

    /* initialize rtk TDCP velocity */
    for (i=0;i<3;i++) rtk->sol.tdcp_vel[i]=0.0; 

    /* initializing memory */
    rs=mat(n,6); rs_old=mat(n_old,6); dts=mat(n,2); dts_old=mat(n_old,2);
    vare=mat(n,1); vare_old=mat(n_old,1); resp=mat(n,1); resp_old=mat(n_old,1);
    azel=zeros(n,2); azel_old=zeros(n_old,2);
    v=mat(nf*n,1);  H=mat(nf*n,nx); var=mat(nf*n,1); P=zeros(nf*n,nf*n); 

    /* reset receiver velocity */
    for (i=0;i<3;i++) rtk->sol.rr[i+3]=0.0;

    /* configured in spp mode */
    if (opt_.mode!=PMODE_SINGLE||opt_.GI_mode!=GINS_OFF) {
        opt_.GI_mode=GINS_OFF; opt_.mode=PMODE_SINGLE; opt_.spp_mode=SPP_LS_C; /* TOdo */
        opt_.sateph =EPHOPT_BRDC; opt_.ionoopt=IONOOPT_BRDC; opt_.tropopt=TROPOPT_SAAS;
        rtk->opt=opt_;
    }

    /* satellite positons, velocities and clocks of current and previous epoch */
    satposs(obs[0].time,obs,n,nav,opt_.sateph,rs,dts,vare,svh);
    satposs(obs_old[0].time,obs_old,n_old,nav,opt_.sateph,rs_old,dts_old,vare_old,svh_old);

    /* receiver positions of current and previous epoch */
    stat=estpos(rtk,obs,n,rs,dts,vare,svh,nav,&opt_,NULL,&sol,azel,vsat,resp);
    stat_old=estpos(rtk,obs_old,n_old,rs_old,dts_old,vare_old,svh_old,nav,&opt_,NULL,&sol_old,azel_old,vsat_old,resp_old);

    /* check solution status */
    if (!stat||!stat_old) {
        pos_flag=tdcp_flag=0; trace(7,"tdcp_vel: estpos error stat=%d, stat_old=%d\n",stat,stat_old); 
    }
    else {
        /* receiver position in the previous epoch and the current epoch in spp mode */
        rtk->sol.stat=stat;
        for (i=0;i<3;i++) {
            rr[i]=rtk->sol.rr[i]=sol.rr[i];
            rr_old[i]=sol_old.rr[i];
        }
    }

    /* if GNSS outage, tdcp fails */
    if (rtk->interval&&timediff(obs[0].time,obs_old[0].time)>rtk->interval) {
        tdcp_flag=0; trace(7,"tdcp_vel: time difference between current and previous epoch is too large, tt=%.2f\n",timediff(obs[0].time,obs_old[0].time));
    }

    /* epoch-to-epoch average velocity estimation based on tdcp */
    if (tdcp_flag) {
        /* check whether a cycle slip occurs in the current epoch observation */
        init_ssatpar(rtk,NULL,n,ssat_slip,SOLQ_NONE);

        /* detect cycle slip by LLI/geometry-free/Melbourne-Wubbena linear combination */
        detslp_ll_ppp(rtk,obs,n);
        detslp_gf_ppp(rtk,obs,n,nav);
        detslp_mw_ppp(rtk,obs,n,nav);

        /* select the common satellites between the previous epoch and the current epoch */
        for (i=0,j=0,k=0;i<n_old&&j<n;i++,j++)
        {
            sys=satsys(obs[j].sat,NULL);   
            if      (obs_old[i].sat<obs[j].sat) j--;
            else if (obs_old[i].sat>obs[j].sat) i--;
            else {
                /* exclude satellites that are not involved in the solution or have unhealthy ephemeris */
                if (satexclude(obs_old[i].sat,vare_old[i],svh_old[i],&opt_)||satexclude(obs[j].sat,vare[j],svh[j],&opt_)) continue;      
                /* exclude satellites with large residuals*/
                if (!vsat_old[i]||!vsat[j]) continue;
                /* exclude satellites that have cycle slips */
                for (m=0;m<nv;m++) {
                    fr=sys2freid(sys,m,&opt_);
                    if (rtk->ssat[obs[j].sat-1].slip[fr]) break; 
                }
                if (m<nv) continue;

                /* save common satellites idx of the previous epoch and the current epoch */
                sat[k]=obs[j].sat;ir_old[k]=i;ir[k++]=j;
            }
        }

        /* construct the error equation with v and H */
        for (m=0;m<nf;m++) {    
            fr=sys2freid(sys,m,&opt_);
            for (i=0;i<k;i++) {
                er=er_old=0.0;
                sys=satsys(sat[i],NULL);            
                freq=sat2freq(sat[i],obs[ir[i]].code[fr],nav);

                /* excluding satellites with missing observations */
                if (obs[ir[i]].L[fr]==0.0||obs_old[ir_old[i]].L[fr]==0.0) continue;

                if ((r=geodist(rs_old+ir_old[i]*6,rr_old,e_old))<=0.0) continue;
                if ((r=geodist(rs+ir[i]*6,rr,e))<=0.0) continue;

                vnadd(3,rs_old+ir_old[i]*6,1.0,rr_old,-1.0,dr_old);
                vnadd(3,rs+ir[i]*6,1.0,rr_old,-1.0,dr);
                
                for (j=0;j<3;j++) {
                    er+=e[j]*dr[j];
                    er_old+=e_old[j]*dr_old[j];
                }
                /* observation vector */
                v[nv]=CLIGHT/freq*(obs[ir[i]].L[fr]-obs_old[ir_old[i]].L[fr])+CLIGHT*(dts[ir[i]*2]-dts_old[ir_old[i]*2])-(er-er_old);

                /* design matrix */   
                for (j=0;j<4;j++) {
                    H[j+nv*4]=j<3?-e[j]:(j==3?1.0:0.0);
                }   

                /* determine the variance of the observations */
                var[nv++]=varerr_spp(&opt_,NULL,&obs[ir[i]],azel[1+ir[i]*2],sys)+vare[ir[i]];
                /* valid satellite observations at the current frequency */
                vnv[m]++;
            }  
        }         

        /* determine the maximum number of satellites available on a single frequency */
        max_vnv=maxobsat(vnv,nf);

        if (max_vnv<nx) {
            tdcp_flag=0;trace(7,"tdcp_vel: not enough valid satellites nv=%d\n",nv);
        }
        else {
            /* mode=(max_vnv>nx)?Robust_OFF:Robust_RES; */
            nv=outrej_spp(nv,NULL,4,0,thres_ouj,v,H,var,NULL,NULL,NULL,vsat,0,NULL);

            /* calculate the weight matrix */
            diag_Cov(nv,var,P,diag_wei);

            /* least square estimation */
            if ((info=lsq_roubst(H,v,P,4,nv,dx,Q,mode))) {
                tdcp_flag=0;trace(7,"tdcp lsq error info=%d\n!",info);
            }

            /* tracefilter(12,TRAE_R|TRAE_H|TRAE_Ppre|TRAE_v|TRAE_xpre,4,nv,P,H,rtk->P,NULL,v,dx,NULL); */

            /* calculate the posterior residuals */
            matmul("NN",nv,nx,1,H,dx,v,1.0,-1.0);

            /* validate solution */
            if (valsol(&sol,azel,vsat,n,opt,v,P,nv,nx)==0) {
                tdcp_flag=0;trace(7,"validation of tdcp solution failed nv=%d\n",nv);
            }
            
            if (tdcp_flag) {
                /* update the receiver velocity with tdcp */
                for (i=0;i<3;i++) rtk->sol.rr[i+3]=sgn*dx[i]/rtk->interval; 
                vel_flag=1;

                /* save tdcp velocity */ 
                matcpy(rtk->sol.tdcp_vel,rtk->sol.rr+3,3,1);
                /* outsolstat(rtk,nav); */                   
            }           
        }      
    }

    /* if tdcp fails, velocity estimation is performed using Doppler observations */
    if (stat&&!vel_flag&&rtk->dopsgn&&estvel(rtk,obs,n,rs,dts,nav,&opt_,&sol,azel,vsat)) {
        vel_flag=1;
        matcpy(rtk->sol.rr+3,sol.rr+3,3,1); /* update the receiver velocity with dopple */
    }

    /* if tdcp and dopple fail, velocity estimation is performed using the spp position difference between epochs */
    if (pos_flag&&!vel_flag) {
        vel_flag=1;
        for (i=0;i<3;i++) rtk->sol.rr[i+3]=sgn*(rr[i]-rr_old[i])/rtk->interval; /* update the receiver velocity with dpos */
    }

    /* only the vehicle velocity exceeds the threshold, the alignment is considered complete */
    if (!align&&vel_flag&&norm(rtk->sol.rr+3,3)<thres) {
        vel_flag=0;
    }

    /* copy TDCP estimated velocity to rtk_main struct */
    if (norm(rtk->sol.rr+3,3)>0.0) matcpy(rtk_main->sol.rr+3,rtk->sol.rr+3,3,1); 

    /* freeing up memory */
    free(rs);   free(rs_old);   free(dts);  free(dts_old);
    free(vare); free(vare_old); free(resp); free(resp_old);
    free(azel); free(azel_old);
    free(v);    free(H);        free(var);  free(P);        

    return vel_flag;
}