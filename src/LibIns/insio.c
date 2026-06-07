/* -----------------
 * insio.c - I/O functions for INS
 * ----------------- */

#include "rtklib.h"

/* get imu data --------------------------------------------------------------
* get imu data to imu[] from imus.data[iimu]
*
*args   : const prcopt_t *popt  I   processing options
*         imud_t *imu            O   imu data
*         imu_t imus             I   imu data
*         int iimu               I   index of imu data
*         const int nn           I   number of imu data to copy
*return : none
*-----------------------------------------------------------------------------*/
extern void getimu(const prcopt_t *popt, imud_t *imu, imu_t imus, int iimu, const int nn)
{
    int i,j;

    /* single sample + previous */   
    if (1==nn) {
        imu[0]=imus.data[iimu];     
    }
    /* double sample */
    else if (2==nn) {
        if (MECH_FORWARD==popt->reverse) {
            for (i=0;i<nn;i++) imu[i]=imus.data[iimu+i];            
        }
        else if (MECH_BACKWARD==popt->reverse) {
            for (i=0;i<nn;i++) imu[i]=imus.data[iimu-i];
        }
    }

    /* backward mechanization negates the sign of the gyroscope output */
    if (MECH_BACKWARD==popt->reverse) {
        for (i=0;i<nn;i++) {
            for (j=0;j<3;j++) {
                imu[i].dw[j]=-imu[i].dw[j];
            }
        }
    }
}

/* add imu data ------------------------------------------------------*/
static int addimudata(imu_t *imu, const imud_t *data)
{
    imud_t *imu_data;

    if (imu->nmax<=imu->n) {
        if (imu->nmax<=0) imu->nmax=NINCIMU; else imu->nmax*=2;
        if (!(imu_data=(imud_t *)realloc(imu->data,sizeof(imud_t)*imu->nmax))) {
            trace(1,"addimudata: malloc error n=%dx%d\n",sizeof(imud_t),imu->nmax);
            free(imu->data); imu->data=NULL; imu->n=imu->nmax=0;
            return -1;
        }
        imu->data=imu_data;
    }
    imu->data[imu->n++]=*data;
    return 1;
}

/* read imu data -----------------------------------------------*/
extern int readimu(gtime_t ts, gtime_t te, const char *file, const prcopt_t *popt, imu_t *imu, int gps_week)
{
    FILE *fp;
    imud_t imud;
    gtime_t time;
    int i,stat=0;
    char buff[256];
    double week,sec,data[6]={0.0},factor=1.0,dw[3],dv[3],Cvb[9]={0.0};

    if (ts.time!=0) ts.time-=1;

    /* calculate the rotation matrix from b' frame to v frame */
    att2Cnb(popt->rotation_angle,Cvb);

    if (!(fp=fopen(file,"r")))
    {
        trace(7,"Error: IMU file open failed: %s!\n",file);
        showerr("Error: IMU file open failed: %s!\n",file);
        return 0;
    }

    imu->data=NULL; imu->n=imu->nmax=0;

    /* convert rate to incremental measurement */
    if (IMUT_RATE==popt->imudatype)
    {
        factor=1.0/popt->insample;
    }

    while (fgets(buff,sizeof(buff),fp))
    {
        /* replace spaces with commas */
        repspace(buff);
        if (sscanf(buff,"%lf,%lf,%lf,%lf,%lf,%lf,%lf",&sec,data,data+1,data+2,data+3,data+4,data+5)<7) continue;
        else if (sscanf(buff,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",&week,&sec,data,data+1,data+2,data+3,data+4,data+5)==8) {}
        else if (sscanf(buff,"%lf,%lf,%lf,%lf,%lf,%lf,%lf",&sec,data,data+1,data+2,data+3,data+4,data+5)==7) week=gps_week;
        else continue; 
        imud.time=gpst2time(week,sec);

        /* screen data by time */
        if ((ts.time!=0&&timediff(imud.time,ts)<0.0)||(te.time!=0&&timediff(imud.time,te)>0.5/popt->insample)) continue;

        if (strstr(popt->imu_order,"AgGd")!=NULL) {
            for (i=0;i<6;i++)
            {
                if (i<3) {
                    dv[i]  =factor*data[i];
                    data[i]=0.0;
                }     
                else if (i<6) {
                    dw[i-3]=factor*data[i]*D2R; 
                    data[i]=0.0;
                }
            }            
        }
        else if (strstr(popt->imu_order,"AgGr")!=NULL) {
            for (i=0;i<6;i++)
            {
                if (i<3) {
                    dv[i]  =factor*data[i];
                    data[i]=0.0;
                }     
                else if (i<6) {
                    dw[i-3]=factor*data[i]; 
                    data[i]=0.0;
                }
            }            
        }
        else if (strstr(popt->imu_order,"GdAg")!=NULL) {
            for (i=0;i<6;i++)
            {
                if (i<3) {
                    dw[i]  =factor*data[i]*D2R;
                    data[i]=0.0;
                }     
                else if (i<6) {
                    dv[i-3]=factor*data[i]; 
                    data[i]=0.0;
                }
            }            
        }
        else if (strstr(popt->imu_order,"GrAg")!=NULL) {
            for (i=0;i<6;i++)
            {
                if (i<3) {
                    dw[i]  =factor*data[i];
                    data[i]=0.0;
                }     
                else if (i<6) {
                    dv[i-3]=factor*data[i]; 
                    data[i]=0.0;
                }
            }            
        } 

        /* Body frame adjustment, normalized to FRU frame */
        if (BODYF_FRD==popt->bodyframe) {
            /* gyroscope data */
            for (i=0;i<3;i++) data[i]=dw[i];
            dw[0]=data[1];
            dw[1]=data[0];
            dw[2]=-data[2];
            /* accelerometer data */
            for (i=0;i<3;i++) data[i]=dv[i];
            dv[0]=data[1];
            dv[1]=data[0];
            dv[2]=-data[2];
            for (i=0;i<3;i++) data[i]=0.0;
        }

        /* rotation angle compensation */
        Mat3mulv(1.0,Cvb,dw,imud.dw);
        Mat3mulv(1.0,Cvb,dv,imud.dv);

        if (norm(imud.dw,3)<=0.0||norm(imud.dv,3)<=0.0) {
            showerr("warning: IMU measurement output is zero: %s, week=%.0f, sec=%.4f!",file,week,sec); 
            trace(7,"warning: IMU measurement output is zero: %s, week=%.0f, sec=%.4f!\n",file,week,sec);
        }
        stat=addimudata(imu,&imud);
    }

    fclose(fp);

    return stat;
}

/* free imu data -----------------------------------------------------*/
extern void freeimu(imu_t *imu)
{
    trace(3,"freeimu:\n");

    free(imu->data); imu->data=NULL; imu->n =imu->nmax =0;
}