/*--------------------------
 * ginsio.c - I/O functions for GINS
 *--------------------------*/

#include "rtklib.h"

/* get pos obs for LC mode */
extern void getpos(pos_t pos, sol_t *sol, int ipos)
{
    posd_t *posd=NULL;
    int i;

    if (ipos<0||ipos>=pos.n) {
        trace(7,"getpos: posfile reaches the end!\n");
        return;
    }

    posd=&pos.data[ipos];

    sol->stat=SOLQ_FLOAT;
    sol->time=posd->time;
    for (i=0;i<3;i++){
        sol->rr [i]=posd->pos[i];
        sol->vel[i]=posd->vel[i];
        sol->qr [i]=posd->qr[i];
        sol->qv [i]=posd->qv[i];
    }
    for (i=3;i<6;i++) sol->qr[i]=posd->qr[i];
}

/* get odo vel based on gins time */
extern int getodovel(odo_t odo, ins_t *ins, gtime_t gins_time, int iodo)
{
    odod_t *odod=NULL;
    int i;

    if (iodo<0||iodo>=odo.n) {
        trace(7,"getodovel: odo file reaches the end!\n");
        return odo.n;
    }

    for (i=iodo;i<odo.n;i++) {
        if (fabs(timediff(gins_time,odo.data[i].time))<=ins->interval/2.0) {
            ins->odo=&odo.data[i];
            break;
        }
    }

    if (i>=odo.n) {
        ins->odo=NULL;
        return iodo; 
    }

    return i;
}

/* add pos data ------------------------------------------------------*/
static int addposdata(pos_t *pos, const posd_t *data)
{
    posd_t *pos_data;
    
    if (pos->nmax<=pos->n) {
        if (pos->nmax<=0) pos->nmax=NMAXPOS; else pos->nmax*=2;
        if (!(pos_data=(posd_t *)realloc(pos->data,sizeof(posd_t)*pos->nmax))) {
            trace(1,"addposdata: malloc error n=%dx%d\n",sizeof(posd_t),pos->nmax);
            free(pos->data); pos->data=NULL; pos->n=pos->nmax=0;
            return -1;
        }
        pos->data=pos_data;
    }
    pos->data[pos->n++]=*data;
    return 1;
}

/* read pos result data -----------------------------------------------*/
extern int readpos(const char *file, const prcopt_t *popt, pos_t *poss, int gps_week)
{
    FILE *fp;
    posd_t posd={0};
    char buff[1024];
    double week,sec,data[12]={0.0},pos[3]={0.0},vel[3]={0.0},Qr[9]={0.0},Qr2[9]={0.0},Qv[9]={0.0},Qv2[9]={0.0},temp;
    int i,vel_flag=0,pos_flag=0,stat=0;

    if (!(fp=fopen(file,"r"))) {
        trace(7,"Error : pos file open failed %s",file);
        showerr("Error : pos file open failed %s",file);
        return 0;
    }

    poss->nmax=poss->n=0;

    while (fgets(buff,sizeof(buff),fp)) 
    {
        if (strchr(buff,'%')) continue;

        /* replace spaces with commas */
        repspace(buff);
        if (sscanf(buff,"%lf,%lf,%lf,%lf,%lf,%lf,%lf",&sec,data,data+1,data+2,data+3,data+4,data+5)<7) continue;
        /* week,sec,pos(llh/xyz),Q,ns,var_pos(xyz) for GINSLIB/RTKLIB */
        else if (sscanf(buff,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                &week,&sec,data,data+1,data+2,data+3,data+4,data+5,data+6,data+7,data+8,data+9,data+10)==13&&data[3]<10) pos_flag=1; 
        /* week,sec,pos(llh/xyz),vel(xyz),var_pos(xyz),var_vel(xyz) */
        else if (sscanf(buff,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                &week,&sec,data,data+1,data+2,data+3,data+4,data+5,data+6,data+7,data+8,data+9,data+10,data+11)==14) vel_flag=1;
        /* week,sec,pos(llh/xyz),vel(xyz),var_pos(xyz),var_vel(xyz) */
        else if (sscanf(buff,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                &sec,data,data+1,data+2,data+3,data+4,data+5,data+6,data+7,data+8,data+9,data+10,data+11)==13) { week=gps_week; vel_flag=1; }       
        /* week,sec,pos(llh/xyz),var_pos(xyz) */
        else if (sscanf(buff,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",&week,&sec,data,data+1,data+2,data+3,data+4,data+5)==8) {}
        /* sec,pos(llh/xyz),var_pos(xyz) */
        else if (sscanf(buff,"%lf,%lf,%lf,%lf,%lf,%lf,%lf",&sec,data,data+1,data+2,data+3,data+4,data+5)==7) week=gps_week;
        else continue; 

        posd.time=gpst2time(week,sec);
        
        /* GINSLIB/RTKLIB format */
        if (pos_flag) {
            for (i=0;i<3;i++) {
                pos[i]=data[i];
                Qr[i+i*3]=data[i+5]*data[i+5]; 
            }  
            Qr[1]=Qr[3]=powstd(data[8]);
            Qr[5]=Qr[7]=powstd(data[9]);
            Qr[2]=Qr[6]=powstd(data[10]);
        }
        else { /* other format */
            for (i=0;i<3;i++) {
                if (vel_flag) {
                    pos[i]=data[i];
                    vel[i]=data[i+3];
                    Qr [i+i*3]=data[i+6]*data[i+6];
                    Qv [i+i*3]=data[i+9]*data[i+9];
                }
                else {
                    pos[i]=data[i];
                    Qr[i+i*3]=data[i+3]*data[i+3];
                }
            }            
        }

        if (POSF_XYZ==popt->postype) {
            for (i=0;i<3;i++) {
                posd.pos[i]=pos[i];
                posd.qr [i]=Qr[i+i*3];
            }
            posd.qr[3]=Qr[1]; posd.qr[4]=Qr[5]; posd.qr[5]=Qr[2];
        }
        else if (POSF_NED==popt->postype) {
            for (i=0;i<2;i++) pos[i]*=D2R;
            pos2ecef(pos,posd.pos);
            /* convert covariance matrix from ned to enu */
            temp=Qr[0]; Qr[0]=Qr[4]; Qr[4]=temp;
            covecef(posd.pos,Qr,Qr2);
            for (i=0;i<3;i++) posd.qr[i]=Qr2[i+i*3];
        }
        else if (POSF_ENU==popt->postype) {
            for (i=0;i<2;i++) pos[i]*=D2R;
            pos2ecef(pos,posd.pos); 
            covecef(posd.pos,Qr,Qr2);
            for (i=0;i<3;i++) posd.qr[i]=Qr2[i+i*3];
        }

        stat=addposdata(poss,&posd);
    }

    fclose(fp);

    return stat;
}

/* add odo data */
static int addododata(odo_t *odo, const odod_t *data)
{
    odod_t *odo_data;

    if (odo->nmax<=odo->n) {
        if (odo->nmax<=0) odo->nmax=NMAXODO; else odo->nmax*=2;
        if (!(odo_data=(odod_t *)realloc(odo->data,sizeof(odod_t)*odo->nmax))) {
            trace(1,"addododata: malloc error n=%dx%d\n",sizeof(odod_t),odo->nmax);
            free(odo->data); odo->data=NULL; odo->n=odo->nmax=0;
            return -1;
        }
        odo->data=odo_data;
    }
    odo->data[odo->n++]=*data;
    return 1;
}

/* read odo data -----------------------------------------------*/
extern int readodo(gtime_t ts, gtime_t te, const char *file, const prcopt_t *popt, odo_t *odo)
{
    FILE *fp;
    odod_t odod;
    gtime_t time;
    char buff[256];
    int stat=0;
    double week,sec,data[3]={0.0},interval=1/popt->insample;

    if (!(fp=fopen(file,"r"))) {
        trace(7,"Error : odo file open failed %s",file);
        showerr("Error : odo file open failed %s",file);
        return 0;
    }

    while (fgets(buff,sizeof(buff),fp)) {

        if (strchr(buff,'%')) continue;

        repspace(buff); /* replace spaces with commas */
        if (sscanf(buff,"%lf,%lf,%lf,%lf,%lf",&week,&sec,data,data+1,data+2)!=5) continue;

        /* integer second filtering */
        if (fmod(sec,1.0)>interval/2.0) continue;

        time=gpst2time(week,sec);
        /* screen data by time */
        if ((ts.time!=0&&timediff(time,ts)<0.0)||(te.time!=0&&timediff(time,te)>interval/2.0)) continue;

        odod.time=time;
        matcpy(odod.odo_vel,data,3,1); /* odo velocity in m/s2 */

        stat=addododata(odo,&odod);
    }

    fclose(fp);

    return stat;
}

/* free pos data -----------------------------------------------------*/
extern void freepos(pos_t *pos)
{
    trace(3,"freepos:\n");

    free(pos->data); pos->data=NULL; pos->n =pos->nmax =0;
}

/* free odo data ----------------------------------------------------*/
extern void freeodo(odo_t *odo)
{
    trace(3,"freeodo:\n");

    free(odo->data); odo->data=NULL; odo->n=odo->nmax=0;
}