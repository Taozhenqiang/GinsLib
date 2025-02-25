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
    sprintf(filopt.trace,"%s.trace",PROGNAME);

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
