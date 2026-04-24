/*------------------------------------------------------------------------------
* rnx2rtkp.c : read rinex obs/nav files and compute receiver positions
*
*          Copyright (C) 2007-2016 by T.TAKASU, All rights reserved.
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:55:16 $
* history : 2007/01/16  1.0 new
*           2007/03/15  1.1 add library mode
*           2007/05/08  1.2 separate from postpos.c
*           2009/01/20  1.3 support rtklib 2.2.0 api
*           2009/12/12  1.4 support glonass
*                           add option -h, -a, -l, -x
*           2010/01/28  1.5 add option -k
*           2010/08/12  1.6 add option -y implementation (2.4.0_p1)
*           2014/01/27  1.7 fix bug on default output time format
*           2015/05/15  1.8 -r or -l options for fixed or ppp-fixed mode
*           2015/06/12  1.9 output patch level in header
*           2016/09/07  1.10 add option -sys
*           2024/09/28  modified by tzq
*-----------------------------------------------------------------------------*/
#include <stdarg.h>
#include "rtklib.h"

#define PROGNAME    "rnx2rtkp"          /* program name */
#define MAXFILE     16                  /* max number of input files */

/* help text -----------------------------------------------------------------*/
static const char *help[]={
"",
" usage: rnx2rtkp [option]... file file [...]",
"",
" Read RINEX OBS/NAV/GNAV/HNAV/CLK, SP3, SBAS message log files and compute ",
" receiver (rover) positions and output position solutions.",
" The first RINEX OBS file shall contain receiver (rover) observations. For the",
" relative mode, the second RINEX OBS file shall contain reference",
" (base station) receiver observations. At least one RINEX NAV/GNAV/HNAV",
" file shall be included in input files. To use SP3 precise ephemeris, specify",
" the path in the files. The extension of the SP3 file shall be .sp3 or .eph.",
" All of the input file paths can include wild-cards (*). To avoid command",
" line deployment of wild-cards, use \"...\" for paths with wild-cards.",
" Command line options are as follows ([]:default). A maximum number of",
" input files is currently set to 16. With -k option, the",
" processing options are input from the configuration file. In this case,",
" command line options precede options in the configuration file.",
"",
" -?        print help",
" -k file   input options from configuration file [off]",
" -o file   set output file [stdout]",
" -ts ds ts start day/time (ds=y/m/d ts=h:m:s) [obs start time]",
" -te de te end day/time   (de=y/m/d te=h:m:s) [obs end time]",
" -ti tint  time interval (sec) [all]",
" -p mode   mode (0:single,1:dgps,2:kinematic,3:static,4:static-start,",
"                 5:moving-base,6:fixed,7:ppp-kinematic,8:ppp-static,9:ppp-fixed) [2]",
" -m mask   elevation mask angle (deg) [15]",
" -sys s[,s...] nav system(s) (s=G:GPS,R:GLO,E:GAL,J:QZS,C:BDS,I:IRN) [G|R]",
" -f freq   number of frequencies for relative mode (1:L1,2:L1+L2,3:L1+L2+L5) [2]",
" -v thres  validation threshold for integer ambiguity (0.0:no AR) [3.0]",
" -b        backward solutions [off]",
" -c        forward/backward combined solutions [off]",
" -i        instantaneous integer ambiguity resolution [off]",
" -h        fix and hold for integer ambiguity resolution [off]",
" -bl bl,std     baseline distance and stdev",
" -e        output x/y/z-ecef position [latitude/longitude/height]",
" -a        output e/n/u-baseline [latitude/longitude/height]",
" -n        output NMEA-0183 GGA sentence [off]",
" -g        output latitude/longitude in the form of ddd mm ss.ss' [ddd.ddd]",
" -t        output time in the form of yyyy/mm/dd hh:mm:ss.ss [sssss.ss]",
" -u        output time in utc [gpst]",
" -d col    number of decimals in time [3]",
" -s sep    field separator [' ']",
" -r x y z  reference (base) receiver ecef pos (m) [average of single pos]",
"           rover receiver ecef pos (m) for fixed or ppp-fixed mode",
" -l lat lon hgt reference (base) receiver latitude/longitude/height (deg/m)",
"           rover latitude/longitude/height for fixed or ppp-fixed mode",
" -y level  output solution status (0:off,1:states,2:residuals) [0]",
" -x level  debug trace level (0:off) [0]",
" --version display release version",
};

/* print help ----------------------------------------------------------------*/
static void printhelp(void)
{
    int i;
    for (i=0;i<(int)(sizeof(help)/sizeof(*help));i++) fprintf(stderr,"%s\n",help[i]);
    exit(0);
}
/* rnx2rtkp main -------------------------------------------------------------*/
int main(int argc, char **argv)
{
    prcopt_t prcopt=prcopt_default;
    solopt_t solopt=solopt_default;
    filopt_t filopt={""};
    gtime_t ts={0},te={0};
    double tint=0.0,es[]={2000,1,1,0,0,0},ee[]={2000,12,31,23,59,59},pos[3],offset[3]={0.0};
    int i,j,n,stat,pos2kml=0,pcolor=0;
    char *refile="",*outfile="",*p;

    clock_t start_time, end_time;
    double total_time;

    prcopt.mode  =PMODE_KINEMA;
    prcopt.navsys=0;
    prcopt.refpos=POSOPT_SINGLE;
    prcopt.glomodear=GLO_ARMODE_ON;
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

    ts=prcopt.ts;te=prcopt.te;

    for (i=1,n=0;i<argc;i++) {
        if (!strcmp(argv[i],"-ts")&&i+2<argc) {
            sscanf(argv[++i],"%lf/%lf/%lf",es,es+1,es+2);
            sscanf(argv[++i],"%lf:%lf:%lf",es+3,es+4,es+5);
            ts=epoch2time(es);
        }
        else if (!strcmp(argv[i],"-te")&&i+2<argc) {
            sscanf(argv[++i],"%lf/%lf/%lf",ee,ee+1,ee+2);
            sscanf(argv[++i],"%lf:%lf:%lf",ee+3,ee+4,ee+5);
            te=epoch2time(ee);
        }
        else if (!strcmp(argv[i],"-ti")&&i+1<argc) tint=atof(argv[++i]);
        else if (!strcmp(argv[i],"-k")&&i+1<argc)  {++i; continue;}
        else if (!strcmp(argv[i],"-p")&&i+1<argc)    prcopt.mode=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-f")&&i+1<argc)    prcopt.nf=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-gins")&&i+1<argc) prcopt.GI_mode=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-ion")&&i+1<argc)  prcopt.ionoopt=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-tro")&&i+1<argc)  prcopt.tropopt=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-eph")&&i+1<argc)  prcopt.sateph=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-flt")&&i+1<argc)  prcopt.filter=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-amb")&&i+1<argc)  prcopt.modear=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-ambt")&&i+1<argc) prcopt.artype=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-solt")&&i+1<argc) prcopt.soltype=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-kml"))            pos2kml=1;
        else if (!strcmp(argv[i],"-pcolor")&&i+1<argc) pcolor=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-ref")&&i+1<argc)  refile=argv[++i];
        else if (!strcmp(argv[i],"-imu")&&i+4<argc) {
            prcopt.psd_gyro=atof(argv[++i]);
            prcopt.psd_acce=atof(argv[++i]);
            prcopt.psd_bg=atof(argv[++i]);
            prcopt.psd_ba=atof(argv[++i]);
        }
        else if (!strcmp(argv[i],"-sys")&&i+1<argc) {
            prcopt.navsys=SYS_NONE;
            for (p=argv[++i];*p;p++) {
                switch (*p) {
                    case 'G': prcopt.navsys|=SYS_GPS;break;
                    case 'R': prcopt.navsys|=SYS_GLO;break;
                    case 'E': prcopt.navsys|=SYS_GAL;break;
                    case 'C': prcopt.navsys|=SYS_CMP;break;
                    case 'J': prcopt.navsys|=SYS_QZS;break;
                    case 'I': prcopt.navsys|=SYS_IRN;break;
                }
            }
        }
        else if (!strcmp(argv[i],"-m")&&i+1<argc) prcopt.elmin=atof(argv[++i])*D2R;
        else if (!strcmp(argv[i],"-v")&&i+1<argc) prcopt.thresar[0]=atof(argv[++i]);
        else if (!strcmp(argv[i],"-s")&&i+1<argc) strcpy(solopt.sep,argv[++i]);
        else if (!strcmp(argv[i],"-d")&&i+1<argc) solopt.timeu=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-b")) prcopt.soltype=1;
        else if (!strcmp(argv[i],"-c")) prcopt.soltype=2;
        else if (!strcmp(argv[i],"-h")) prcopt.modear=3;
        else if (!strcmp(argv[i],"-t")) solopt.timef=1;
        else if (!strcmp(argv[i],"-u")) solopt.times=TIMES_UTC;
        else if (!strcmp(argv[i],"-e")) solopt.posf=SOLF_XYZ;
        else if (!strcmp(argv[i],"-a")) solopt.posf=SOLF_ENU;
        else if (!strcmp(argv[i],"-n")) solopt.posf=SOLF_NMEA;
        else if (!strcmp(argv[i],"-g")) solopt.degf=1;
        else if (!strcmp(argv[i],"-bl")&&i+2<argc) {
            for (j=0;j<2;j++) prcopt.baseline[j]=atof(argv[++i]);
        }
        else if (!strcmp(argv[i],"-r")&&i+3<argc) {
            prcopt.refpos=prcopt.rovpos=POSOPT_POS_XYZ;
            for (j=0;j<3;j++) prcopt.rb[j]=atof(argv[++i]);
            matcpy(prcopt.ru,prcopt.rb,3,1);
        }
        else if (!strcmp(argv[i],"-l")&&i+3<argc) {
            prcopt.refpos=prcopt.rovpos=POSOPT_POS_LLH;
            for (j=0;j<3;j++) pos[j]=atof(argv[++i]);
            for (j=0;j<2;j++) pos[j]*=D2R;
            pos2ecef(pos,prcopt.rb);
            matcpy(prcopt.ru,prcopt.rb,3,1);
        }
        else if (!strcmp(argv[i],"-y")&&i+1<argc) solopt.sstat=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-x")&&i+1<argc) solopt.trace=atoi(argv[++i]);
        else if (!strcmp(argv[i], "--version")) {
            fprintf(stderr, "rnx2rtkp RTKLIB %s %s\n", VER_RTKLIB, PATCH_LEVEL);
            exit(0);
        }
        else if (*argv[i]=='-') printhelp();
    }
    if (!prcopt.navsys) {
        prcopt.navsys=SYS_GPS|SYS_GLO;
        showerr("No navigation system specified, defaulting to GPS+GLO!\n");
    }

    start_time=clock(); /* start time */
    
    /* post-processing positioning */
    stat=execses(ts,te,tint,&prcopt,&solopt,&filopt);

    /* convert solution to kml file */
    if (pos2kml) {
        convkml(filopt.sol,refile,outfile,ts,te,&solopt,tint,0,offset,2,pcolor,1,0,0);   
    }    

    end_time=clock(); /* end time */
    total_time=(double)(end_time-start_time)/CLOCKS_PER_SEC;
    showerr("Total processing time of program: %f sec!", total_time); 

    if (!stat) fprintf(stderr,"%40s\r","");
    return stat?EXIT_FAILURE:0;
}
