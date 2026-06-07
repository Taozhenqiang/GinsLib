/*------------------------------------------------------------------------------
* inspure.c : pure inertial navigation solution
*
*          Copyright (C) 2007-2016 by T.TAKASU, All rights reserved.
*
* version : $Revision: 1.1 $ $Date: 2024/12/03 21:55:16 $
* history : 2024/12/03  1.0 new
*-----------------------------------------------------------------------------*/
#include <stdarg.h>
#include "rtklib.h"

#define PROGNAME    "inspure"          /* program name */
#define MAXFILE     1                  /* max number of input files */

/* solution option to field separator ----------------------------------------*/
static const char *opt2sep(const solopt_t *opt)
{
    if (!*opt->sep) return " ";
    else if (!strcmp(opt->sep,"\\t")) return "\t";
    return opt->sep;
}

/* write header to output file -----------------------------------------------*/
static int outhead(const char *outfile, const imu_t imu, const solopt_t *sopt)
{
    FILE *fp=stdout;
    int w1,w2;
    double t1,t2;
    char s2[32],s3[32];
    gtime_t ts,te;

    trace(3,"outhead: outfile=%s\n",outfile);
    const char *sep=opt2sep(sopt);

    if (*outfile) {
        createdir(outfile);

        if (!(fp=fopen(outfile,"wb"))) {
            trace(1,"error : open output file %s",outfile);
            return 0;
        }
    }

    fprintf(fp,"%s program   : RTKLIB ver.%s %s\n",COMMENTH,VER_RTKLIB,PATCH_LEVEL);

    ts=imu.data[0].time;
    te=imu.data[imu.n-1].time;
    t1=time2gpst(ts,&w1);
    t2=time2gpst(te,&w2);

    time2str(ts,s2,1);
    time2str(te,s3,1);

    fprintf(fp,"%s obs start : %s GPST (week%04d %8.1fs)\n",COMMENTH,s2,w1,t1);
    fprintf(fp,"%s obs end   : %s GPST (week%04d %8.1fs)\n",COMMENTH,s3,w2,t2);

    fprintf(fp,"%s\n",COMMENTH);

    fprintf(fp,"%s%10s%s%10s%s%10s%s%10s%s%10s%s%10s%s%10s%s%10s%s%10s%s%10s\n",
                        COMMENTH,"GPST",sep,"x-ecef(m)",sep,"y-ecef(m)",sep,"z-ecef(m)",sep,"ve(m/s)",sep,
                       "vn(m/s)",sep,"vu(m/s)",sep,"pitch(deg)",sep,"roll(deg)",sep,"yaw(deg)");

    if (*outfile) fclose(fp);
} 

/* open output file for append -----------------------------------------------*/
static FILE *openfile(const char *outfile)
{
    trace(3,"openfile: outfile=%s\n",outfile);

    return !*outfile?stdout:fopen(outfile,"ab");
}

/* output solution body --------------------------------------------------------
* output solution body to file
* args   : FILE   *fp       I   output file pointer
*          sol_t  *sol      I   solution
*          solopt_t *opt    I   solution options
* return : none
*-----------------------------------------------------------------------------*/
static void ioutsol(FILE *fp, ins_t *ins, const prcopt_t *popt, const solopt_t *opt)
{
    const char *sep=opt2sep(opt); 
    int week;
    double pos[3],tow;

    pos2ecef(ins->pos,pos);

    tow=time2gpst(ins->time,&week);

    fprintf(fp,"%4d%s%14.4f%s",week,sep,tow,sep);
    /* fprintf(fp,"%14.4f%s%14.4f%s%14.4f%s%14.4f%s%14.4f%s%14.4f%s%14.4f%s%14.4f%s%14.4f%s%14.5f%s%14.5f%s%14.5f\n",
               pos[0],sep,pos[1],sep,pos[2],sep,ins->vel[0],sep,ins->vel[1],sep,ins->vel[2],
               sep,ins->att[0]*R2D,sep,ins->att[1]*R2D,sep,ins->att[2]*R2D,sep,ins->dv[0]/ins->interval,sep,ins->dv[1]/ins->interval,sep,ins->dv[2]/ins->interval); */
    fprintf(fp,"%14.10f%s%14.10f%s%14.6f%s%14.4f%s%14.4f%s%14.4f%s%14.4f%s%14.4f%s%14.4f\n",
               ins->pos[0],sep,ins->pos[1],sep,ins->pos[2],sep,ins->vel[0],sep,ins->vel[1],sep,ins->vel[2],
               sep,ins->att[0]*R2D,sep,ins->att[1]*R2D,sep,ins->att[2]*R2D);
}

/* pure inertial navigation */
extern int inspure(gtime_t ts, gtime_t te, const prcopt_t *popt, const solopt_t *sopt, const char *infile, const char *outfile)
{
    int i,n;
    double sec=0.0,thres=0.0;

    ins_t *ins = (ins_t *)malloc(sizeof(ins_t));

    readimu(ts,te,infile,popt,&imus,0);

    outhead(outfile,imus,sopt);
    FILE *fp=openfile(outfile);

    ins_init(ins,popt);

    n=imus.n;

    for (i=0;i<n;i++) {
        /* DebugTime(imus.data[i].time,436804,2188); */
        sec=imus.data[i].time.sec;
        thres=sec>0.5?(1-sec):sec;
        earth_init(ins->pos,ins->vel,&ins->eth);
        ins_mech(ins,&imus.data[i],popt);
        /* update_ins(ins); */
        ioutsol(fp,ins,popt,sopt);
        /* if (thres<=/2.0)
        {
            ioutsol(fp,ins,popt,sopt);            
        } */

    }
    
    free(ins);
}

/* inspure main -------------------------------------------------------------*/
int main(int argc, char **argv)
{
    prcopt_t prcopt=prcopt_default;
    solopt_t solopt=solopt_default;
    filopt_t filopt={""};
    gtime_t ts={0},te={0};
    double tint=0.0,es[]={2000,1,1,0,0,0},ee[]={2000,12,31,23,59,59},pos[3];
    int i,j,n,ret;
    const char *infile,*outfile="",*p;

    prcopt.mode  =PMODE_INSPURE;
    solopt.timef=0;
    sprintf(solopt.prog ,"%s ver.%s %s",PROGNAME,VER_RTKLIB,PATCH_LEVEL);

    /* load options from configuration file */
    for (i=1;i<argc;i++) {
        if (!strcmp(argv[i],"-k")&&i+1<argc) {
            resetsysopts();
            if (!loadopts(argv[++i],sysopts)) return EXIT_FAILURE;
            getsysopts(&prcopt,&solopt,&filopt);
        }
    }
    for (i=1,n=0;i<argc;i++) {
        if      (!strcmp(argv[i],"-i")&&i+1<argc) infile=argv[++i];
        else if (!strcmp(argv[i],"-y")&&i+1<argc) solopt.sstat=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-x")&&i+1<argc) solopt.trace=atoi(argv[++i]);        
        else if (!strcmp(argv[i],"-o")&&i+1<argc) outfile=argv[++i];
        else if (!strcmp(argv[i],"-ts")&&i+1<argc) {
            sscanf(argv[++i],"%lf/%lf/%lf %lf:%lf:%lf",es,es+1,es+2,es+3,es+4,es+5);
            ts=epoch2time(es);
        }
        else if (!strcmp(argv[i],"-te")&&i+1<argc) {
            sscanf(argv[++i],"%lf/%lf/%lf %lf:%lf:%lf",ee,ee+1,ee+2,ee+3,ee+4,ee+5);
            te=epoch2time(ee);
        }
    }
    ret=inspure(ts,te,&prcopt,&solopt,infile,outfile);

    if (!ret) fprintf(stderr,"%40s\r","");
    return ret?EXIT_FAILURE:0;
}
