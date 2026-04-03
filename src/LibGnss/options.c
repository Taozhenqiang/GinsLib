/*------------------------------------------------------------------------------
* options.c : options functions
*
*          Copyright (C) 2010-2020 by T.TAKASU, All rights reserved.
*
* version : $Revision:$ $Date:$
* history : 2010/07/20  1.1  moved from postpos.c
*                            added api:
*                                searchopt(),str2opt(),opt2str(),opt2buf(),
*                                loadopts(),saveopts(),resetsysopts(),
*                                getsysopts(),setsysopts()
*           2010/09/11  1.2  add options
*                                pos2-elmaskhold,pos1->snrmaskena
*                                pos1-snrmask1,2,3
*           2013/03/11  1.3  add pos1-posopt1,2,3,4,5,pos2-syncsol
*                                misc-rnxopt1,2,pos1-snrmask_r,_b,_L1,_L2,_L5
*           2014/10/21  1.4  add pos2-bdsarmode
*           2015/02/20  1.4  add ppp-fixed as pos1-posmode option
*           2015/05/10  1.5  add pos2-arthres1,2,3,4
*           2015/05/31  1.6  add pos2-armaxiter, pos1-posopt6
*                            add selection precise for pos1-pospot3
*           2015/11/26  1.7  modify pos1-frequency 4:l1+l2+l5+l6 -> l1+l5
*           2015/12/05  1.8  add misc-pppopt
*           2016/06/10  1.9  add ant2-maxaveep,ant2-initrst
*           2016/07/31  1.10 add out-outsingle,out-maxsolstd
*           2017/06/14  1.11 add out-outvel
*           2020/11/30  1.12 change options pos1-frequency, pos1-ionoopt,
*                             pos1-tropopt, pos1-sateph, pos1-navsys,
*                             pos2-gloarmode,
*-----------------------------------------------------------------------------*/
#define _POSIX_C_SOURCE 199506
#include "rtklib.h"

/* system options buffer -----------------------------------------------------*/
static prcopt_t prcopt_;
static solopt_t solopt_;
static filopt_t filopt_;
static double elmask_,elmaskar_,elmaskhold_;
static double antpos_[2][3];
static char exsats_[1024];
static char snrmask_[NFREQ][1024];
static char time_[2][1024];
static char fre_[MAXSYS][1024];
static char constraint_[1024];
static char rotation_angle_[1024];
static char install_angle_[1024];
static char initpose_[3][1024];
static char initunc_[3][1024];
static char lever_[1024];
static char lever_nhc_[1024];
static char stat_[statopt];

/* system options table ------------------------------------------------------*/
#define SWTOPT  "0:off,1:on"
#define GIOPT   "0:off,1:LC,2:TC,3:STC"
#define POSDOPT "0:ECEF,1:NED,2:ENU"
#define MODOPT  "0:single,1:dgps,2:kinematic,3:static,4:tight,5:static-start,6:movingbase,7:fixed,8:ppp-kine,9:ppp-static,10:ppp-fixed,11:LC-pos,12:PINS"
#define FILOPT  "0:KF,1:Robust_INO,2:Robust_RES,3:Robust_Chi,4:Robust_ST,5:Robust_MST"
#define MESOPT  "0:IGG3,1:Huber,2:MCKF"
#define TYPOPT  "0:forward,1:backward,2:combined,3:combined-nophasereset"
#define IONOPT  "0:off,1:brdc,2:sbas,3:dual-freq,4:est-stec,5:ionex-tec,6:qzs-brdc"
#define TRPOPT  "0:off,1:saas,2:sbas,3:est-ztd,4:est-ztdgrad"
#define EPHOPT  "0:brdc,1:precise,2:brdc+sbas,3:brdc+ssrapc,4:brdc+ssrcom"
#define NAVOPT  "1:gps+2:sbas+4:glo+8:gal+16:qzs+32:bds+64:navic"
#define GAROPT  "0:off,1:on,2:autocal,3:fix-and-hold"
#define WEIGHTOPT "0:elevation,1:snr"
#define SOLOPT  "0:llh,1:xyz,2:enu,3:nmea"
#define TSYOPT  "0:gpst,1:utc,2:jst"
#define TFTOPT  "0:tow,1:hms"
#define DFTOPT  "0:deg,1:dms"
#define HGTOPT  "0:ellipsoidal,1:geodetic"
#define GEOOPT  "0:internal,1:egm96,2:egm08_2.5,3:egm08_1,4:gsi2000"
#define STSOPT  "0:off,1:state,2:residual"
#define ARMOPT  "0:off,1:continuous,2:instantaneous,3:fix-and-hold"
#define STAOPT  "0:single,1:vrs"
#define POSOPT  "0:llh,1:xyz,2:single,3:posfile,4:rinexhead,5:rtcm,6:raw"
#define TIDEOPT "0:off,1:on,2:otl"
#define PHWOPT  "0:off,1:on,2:precise"

EXPORT opt_t sysopts[]={
    {"pos1-GINS",       3,  (void *)&prcopt_.GI_mode,    GIOPT  },
    {"pos1-postype",    3,  (void *)&prcopt_.postype,    POSDOPT},
    {"pos1-week",       0,  (void *)&prcopt_.week,       ""     },
    {"pos1-mfspp",      3,  (void *)&prcopt_.mfspp,      SWTOPT },
    {"pos1-respp",      3,  (void *)&prcopt_.respp,      SWTOPT },
    {"pos1-cdspp",      3,  (void *)&prcopt_.cdspp,      SWTOPT },
    {"pos1-posmode",    3,  (void *)&prcopt_.mode,       MODOPT },
    {"pos1-frequency",  0,  (void *)&prcopt_.nf,         ""     },    
    {"pos1-ts",         2,  (void *)&time_[0],           ""     },
    {"pos1-te",         2,  (void *)&time_[1],           ""     },
    {"pos1-GPSfreid",   2,  (void *)&fre_[0],            ""     },
    {"pos1-GLOfreid",   2,  (void *)&fre_[1],            ""     },
    {"pos1-GALfreid",   2,  (void *)&fre_[2],            ""     },
    {"pos1-BDSfreid",   2,  (void *)&fre_[3],            ""     },    
    {"pos1-QZSfreid",   2,  (void *)&fre_[4],            ""     },
    {"pos1-soltype",    3,  (void *)&prcopt_.soltype,    TYPOPT },
    {"pos1-elmask",     1,  (void *)&elmask_,            "deg"  },
    {"pos1-snrmask_r",  3,  (void *)&prcopt_.snrmask.ena[0],SWTOPT},
    {"pos1-snrmask_b",  3,  (void *)&prcopt_.snrmask.ena[1],SWTOPT},
    {"pos1-snrmask_L1", 2,  (void *)snrmask_[0],         ""     },
    {"pos1-snrmask_L2", 2,  (void *)snrmask_[1],         ""     },
    {"pos1-snrmask_L5", 2,  (void *)snrmask_[2],         ""     },
    {"pos1-dynamics",   3,  (void *)&prcopt_.dynamics,   SWTOPT },
    {"pos1-tidecorr",   0,  (void *)&prcopt_.tidecorr,   ""     },
    {"pos1-sysisb",     0,  (void *)&prcopt_.sysisb,     ""     },
    {"pos1-ionoopt",    3,  (void *)&prcopt_.ionoopt,    IONOPT },
    {"pos1-ionoise",    0,  (void *)&prcopt_.ionoise,    ""     },
    {"pos1-tropopt",    3,  (void *)&prcopt_.tropopt,    TRPOPT },
    {"pos1-sateph",     3,  (void *)&prcopt_.sateph,     EPHOPT },
    {"pos1-posopt1",    3,  (void *)&prcopt_.posopt[0],  SWTOPT },
    {"pos1-posopt2",    3,  (void *)&prcopt_.posopt[1],  SWTOPT },
    {"pos1-posopt3",    3,  (void *)&prcopt_.posopt[2],  PHWOPT },
    {"pos1-posopt4",    3,  (void *)&prcopt_.posopt[3],  SWTOPT },
    {"pos1-posopt5",    3,  (void *)&prcopt_.posopt[4],  SWTOPT },
    {"pos1-posopt6",    3,  (void *)&prcopt_.posopt[5],  SWTOPT },
    {"pos1-exclsats",   2,  (void *)exsats_,             "prn ..."},
    {"pos1-bds2",       3,  (void *)&prcopt_.bdsflag[0], SWTOPT },
    {"pos1-bds3",       3,  (void *)&prcopt_.bdsflag[1], SWTOPT },
    {"pos1-navsys",     0,  (void *)&prcopt_.navsys,     NAVOPT },
    {"pos1-filter",     3,  (void *)&prcopt_.filter,     FILOPT },
    {"pos1-M_robust",   3,  (void *)&prcopt_.M_robust,   MESOPT },

    {"pos2-artype",     0,  (void *)&prcopt_.artype,     ""     },    
    {"pos2-armode",     3,  (void *)&prcopt_.modear,     ARMOPT },
    {"pos2-gloarmode",  3,  (void *)&prcopt_.glomodear,  GAROPT },
    {"pos2-bdsarmode",  3,  (void *)&prcopt_.bdsmodear,  SWTOPT },
    {"pos2-arfilter",   3,  (void *)&prcopt_.arfilter,   SWTOPT },
    {"pos2-arthres",    1,  (void *)&prcopt_.thresar[0], ""     },
    {"pos2-arthresmin", 1,  (void *)&prcopt_.thresar[5], ""     },
    {"pos2-arthresmax", 1,  (void *)&prcopt_.thresar[6], ""     },
    {"pos2-arthres1",   1,  (void *)&prcopt_.thresar[1], ""     },
    {"pos2-arthres2",   1,  (void *)&prcopt_.thresar[2], ""     },
    {"pos2-arthres3",   1,  (void *)&prcopt_.thresar[3], ""     },
    {"pos2-arthres4",   1,  (void *)&prcopt_.thresar[4], ""     },
    {"pos2-varholdamb", 1,  (void *)&prcopt_.varholdamb, "cyc^2"},
    {"pos2-gainholdamb",1,  (void *)&prcopt_.gainholdamb,""     },
    {"pos2-arlockcnt",  0,  (void *)&prcopt_.minlock,    ""     },
    {"pos2-minfixsats", 0,  (void *)&prcopt_.minfixsats, ""     },
    {"pos2-minholdsats",0,  (void *)&prcopt_.minholdsats,""     },
    {"pos2-mindropsats",0,  (void *)&prcopt_.mindropsats,""     },
    {"pos2-arelmask",   1,  (void *)&elmaskar_,          "deg"  },
    {"pos2-arminfix",   0,  (void *)&prcopt_.minfix,     ""     },
    {"pos2-armaxiter",  0,  (void *)&prcopt_.armaxiter,  ""     },
    {"pos2-elmaskhold", 1,  (void *)&elmaskhold_,        "deg"  },
    {"pos2-aroutcnt",   0,  (void *)&prcopt_.maxout,     ""     },
    {"pos2-maxage",     1,  (void *)&prcopt_.maxtdiff,   "s"    },
    {"pos2-slipthres",  1,  (void *)&prcopt_.thresslip,  "m"    },
    {"pos2-dopthres",   1,  (void *)&prcopt_.thresdop,   "m"    },
    {"pos2-rejphase",   1,  (void *)&prcopt_.maxinno[0], "m"    },
    {"pos2-rejcode",    1,  (void *)&prcopt_.maxinno[1], "m"    },
    {"pos2-niter",      0,  (void *)&prcopt_.niter,      ""     },
    {"pos2-baselen",    1,  (void *)&prcopt_.baseline[0],"m"    },
    {"pos2-basesig",    1,  (void *)&prcopt_.baseline[1],"m"    },

    {"ins-type",        2,  (void *)&filopt_.ins_type,   ""     }, 
    {"ins-dataorder",   2,  (void *)&prcopt_.imu_order,  ""     },
    {"ins-bodyframe",   0,  (void *)&prcopt_.bodyframe,  ""     },
    {"ins-imudatype",   0,  (void *)&prcopt_.imudatype,  ""     },
    {"ins-nnts",        0,  (void *)&prcopt_.nn,         ""     },
    {"ins-insample",    0,  (void *)&prcopt_.insample,   ""     },
    {"ins-aligntype",   0,  (void *)&prcopt_.alingetype,  ""    },
    {"ins-attupdtype",  0,  (void *)&prcopt_.att_type,    ""    },  
    {"ins-errmodel",    0,  (void *)&prcopt_.err_model,   ""    },  
    {"ins-constraints", 2,  (void *)&constraint_,        ""     },
    {"ins-install_angle",2, (void *)&install_angle_,      ""    },
    {"ins-nhc_lever",   2,  (void *)&lever_nhc_,          ""    },
    {"ins-zupt_gthres", 1,  (void *)&prcopt_.zupt_gthres, ""    },
    {"ins-rotaion_angle",2, (void *)&rotation_angle_,     ""    },
    {"ins-initpos",     2,  (void *)&initpose_[0],       ""     },
    {"ins-initvel",     2,  (void *)&initpose_[1],       ""     },
    {"ins-initatt",     2,  (void *)&initpose_[2],       ""     },
    {"ins-lever",       2,  (void *)&lever_,             ""     },
    {"ins-init_pos_unc",2,  (void *)&initunc_[0],       ""      },
    {"ins-init_vel_unc",2,  (void *)&initunc_[1],       ""      },
    {"ins-init_att_unc",2,  (void *)&initunc_[2],       ""      },
    {"ins-init_bg_unc", 1,  (void *)&prcopt_.init_bg_unc,""     },
    {"ins-init_ba_unc", 1,  (void *)&prcopt_.init_ba_unc,""     },
    {"ins-corr_time",   1,  (void *)&prcopt_.corr_time,  ""     },     
    {"ins-psd_gyro",    1,  (void *)&prcopt_.psd_gyro,   ""     },
    {"ins-psd_acce",    1,  (void *)&prcopt_.psd_acce,   ""     }, 
    {"ins-psd_bg",      1,  (void *)&prcopt_.psd_bg,     ""     },
    {"ins-psd_ba",      1,  (void *)&prcopt_.psd_ba,     ""     },             
    
    {"out-solformat",   3,  (void *)&solopt_.posf,       SOLOPT },
    {"out-outhead",     3,  (void *)&solopt_.outhead,    SWTOPT },
    {"out-outopt",      3,  (void *)&solopt_.outopt,     SWTOPT },
    {"out-outpos",      0,  (void *)&prcopt_.outpos,     ""     },
    {"out-outvel",      3,  (void *)&solopt_.outvel,     SWTOPT },
    {"out-outatt",      3,  (void *)&solopt_.outatt,     SWTOPT },
    {"out-outbga",      3,  (void *)&solopt_.outbga,     SWTOPT },
    {"out-outiFlag",    3,  (void *)&solopt_.outiflag,   SWTOPT },
    {"out-timesys",     3,  (void *)&solopt_.times,      TSYOPT },
    {"out-timeform",    3,  (void *)&solopt_.timef,      TFTOPT },
    {"out-timendec",    0,  (void *)&solopt_.timeu,      ""     },
    {"out-degform",     3,  (void *)&solopt_.degf,       DFTOPT },
    {"out-fieldsep",    2,  (void *)&solopt_.sep,        ""     },
    {"out-outsingle",   3,  (void *)&prcopt_.outsingle,  SWTOPT },
    {"out-maxsolstd",   1,  (void *)&solopt_.maxsolstd,  "m"    },
    {"out-height",      3,  (void *)&solopt_.height,     HGTOPT },
    {"out-geoid",       3,  (void *)&solopt_.geoid,      GEOOPT },
    {"out-nmeaintv1",   1,  (void *)&solopt_.nmeaintv[0],"s"    },
    {"out-nmeaintv2",   1,  (void *)&solopt_.nmeaintv[1],"s"    },
    {"out-outstat",     3,  (void *)&solopt_.sstat,      STSOPT },
    {"out-outazel",     3,  (void *)&solopt_.azel,       SWTOPT },
    {"out-outsatdop",   3,  (void *)&solopt_.satdop,     SWTOPT },
    {"out-outipos",     3,  (void *)&solopt_.ipos,       SWTOPT },
    {"out-statopt",     2,  (void *)&stat_,              ""     },

    {"stats-eratio1",   1,  (void *)&prcopt_.eratio[0],  ""     },
    {"stats-eratio2",   1,  (void *)&prcopt_.eratio[1],  ""     },
    {"stats-eratio5",   1,  (void *)&prcopt_.eratio[2],  ""     },
    {"stats-errphase",  1,  (void *)&prcopt_.err[1],     "m"    },
    {"stats-errphaseel",1,  (void *)&prcopt_.err[2],     "m"    },
    {"stats-errphasebl",1,  (void *)&prcopt_.err[3],     "m/10km"},
    {"stats-errdoppler",1,  (void *)&prcopt_.err[4],     "Hz"   },
    {"stats-snrmax",    1,  (void *)&prcopt_.err[5],     "dB.Hz"},
    {"stats-errsnr",    1,  (void *)&prcopt_.err[6],     "m"    },
    {"stats-errrcv",    1,  (void *)&prcopt_.err[7],     " "    },
    {"stats-stdbias",   1,  (void *)&prcopt_.std[0],     "m"    },
    {"stats-stdiono",   1,  (void *)&prcopt_.std[1],     "m"    },
    {"stats-stdtrop",   1,  (void *)&prcopt_.std[2],     "m"    },
    {"stats-prnaccelh", 1,  (void *)&prcopt_.prn[3],     "m/s^2"},
    {"stats-prnaccelv", 1,  (void *)&prcopt_.prn[4],     "m/s^2"},
    {"stats-prnbias",   1,  (void *)&prcopt_.prn[0],     "m"    },
    {"stats-prniono",   1,  (void *)&prcopt_.prn[1],     "m"    },
    {"stats-prntrop",   1,  (void *)&prcopt_.prn[2],     "m"    },
    {"stats-prnpos",    1,  (void *)&prcopt_.prn[5],     "m"    },
    {"stats-clkstab",   1,  (void *)&prcopt_.sclkstab,   "s/s"  },
    
    {"ant1-postype",    3,  (void *)&prcopt_.rovpos,     POSOPT },
    {"ant1-pos1",       1,  (void *)&antpos_[0][0],      "deg|m"},
    {"ant1-pos2",       1,  (void *)&antpos_[0][1],      "deg|m"},
    {"ant1-pos3",       1,  (void *)&antpos_[0][2],      "m|m"  },
    {"ant1-anttype",    2,  (void *)prcopt_.anttype[0],  ""     },
    {"ant1-antdele",    1,  (void *)&prcopt_.antdel[0][0],"m"   },
    {"ant1-antdeln",    1,  (void *)&prcopt_.antdel[0][1],"m"   },
    {"ant1-antdelu",    1,  (void *)&prcopt_.antdel[0][2],"m"   },
    
    {"ant2-statype",    3,  (void *)&prcopt_.statype,    STAOPT },
    {"ant2-postype",    3,  (void *)&prcopt_.refpos,     POSOPT },
    {"ant2-pos1",       1,  (void *)&antpos_[1][0],      "deg|m"},
    {"ant2-pos2",       1,  (void *)&antpos_[1][1],      "deg|m"},
    {"ant2-pos3",       1,  (void *)&antpos_[1][2],      "m|m"  },
    {"ant2-anttype",    2,  (void *)prcopt_.anttype[1],  ""     },
    {"ant2-antdele",    1,  (void *)&prcopt_.antdel[1][0],"m"   },
    {"ant2-antdeln",    1,  (void *)&prcopt_.antdel[1][1],"m"   },
    {"ant2-antdelu",    1,  (void *)&prcopt_.antdel[1][2],"m"   },
    {"ant2-maxaveep",   0,  (void *)&prcopt_.maxaveep    ,""    },
    {"ant2-initrst",    3,  (void *)&prcopt_.initrst,    SWTOPT },
    
    {"misc-timeinterp", 3,  (void *)&prcopt_.intpref,    SWTOPT },
    {"misc-sbasatsel",  0,  (void *)&prcopt_.sbassatsel, "0:all"},
    {"misc-rnxopt1",    2,  (void *)prcopt_.rnxopt[0],   ""     },
    {"misc-rnxopt2",    2,  (void *)prcopt_.rnxopt[1],   ""     },
    {"misc-pppopt",     2,  (void *)prcopt_.pppopt,      ""     },
   
    {"file-obsufile",   2,  (void *)&filopt_.obs_u,      ""     },
    {"file-obsbfile",   2,  (void *)&filopt_.obs_b,      ""     },
    {"file-navfile",    2,  (void *)&filopt_.nav,        ""     },
    {"file-sp3file",    2,  (void *)&filopt_.sp3,        ""     },
    {"file-clkfile",    2,  (void *)&filopt_.clk,        ""     },
    {"file-imufile",    2,  (void *)&filopt_.imu,        ""     }, 
    {"file-posfile",    2,  (void *)&filopt_.pos,        ""     },    
    {"file-antfile",    2,  (void *)&filopt_.antp,       ""     },
    {"file-mgexdcbfile",2,  (void *)&filopt_.mgex_dcb,   ""     },
    {"file-staposfile", 2,  (void *)&filopt_.stapos,     ""     },
    {"file-geoidfile",  2,  (void *)&filopt_.geoid,      ""     },
    {"file-ionofile",   2,  (void *)&filopt_.iono,       ""     },
    {"file-dcbfile",    2,  (void *)&filopt_.dcb,        ""     },
    {"file-eopfile",    2,  (void *)&filopt_.eop,        ""     },
    {"file-blqfile",    2,  (void *)&filopt_.blq,        ""     },
    {"file-tempdir",    2,  (void *)&filopt_.tempdir,    ""     },
    {"file-geexefile",  2,  (void *)&filopt_.geexe,      ""     },
    
    {"",0,NULL,""} /* terminator */
};
/* discard space characters at tail ------------------------------------------*/
static void chop(char *str)
{
    char *p;
    if ((p=strchr(str,'#'))) *p='\0'; /* comment */
    for (p=str+strlen(str)-1;p>=str&&!isgraph((int)*p);p--) *p='\0';
}
/* enum to string ------------------------------------------------------------*/
static int enum2str(char *s, const char *comment, int val)
{
    char str[32],*p,*q;
    int n;
    
    n=sprintf(str,"%d:",val);
    if (!(p=strstr(comment,str))) {
        return sprintf(s,"%d",val);
    }
    if (!(q=strchr(p+n,','))&&!(q=strchr(p+n,')'))) {
        strcpy(s,p+n);
        return (int)strlen(p+n);
    }
    strncpy(s,p+n,q-p-n); s[q-p-n]='\0';
    return (int)(q-p-n);
}
/* String to enum ------------------------------------------------------------
 * Note if str is empty then the first comment digit is returned.
 */
static int str2enum(const char *str, const char *comment, int *val) {
    for (const char *p = comment;; p++) {
        p=strstr(p, str);
        if (!p) break;
        size_t i=p-comment;
        if (i<1) continue;
        if (comment[--i]!=':') continue;
        /* Search for preceding digits */
        size_t j=i;
        while (j>0) {
            char c=comment[j-1];
            if (c<'0' || c>'9') break;
            j--;
        }
        if (j==i) continue; /* No digits found */
        return sscanf(comment+j,"%d",val)==1;
    }
    char s[32];
    snprintf(s,sizeof(s),"%.30s:",str);
    const char *p=strstr(comment,s);
    if (p) { /* Number */
        return sscanf(p,"%d",val)==1;
    }
    return 0;
}
/* search option ---------------------------------------------------------------
* search option record
* args   : char   *name     I  option name
*          opt_t  *opts     I  options table
*                              (terminated with table[i].name="")
* return : option record (NULL: not found)
*-----------------------------------------------------------------------------*/
extern opt_t *searchopt(const char *name, const opt_t *opts)
{
    int i;
    
    trace(3,"searchopt: name=%s\n",name);
    
    for (i=0;*opts[i].name;i++) {
        if (strstr(opts[i].name,name)) return (opt_t *)(opts+i);
    }
    return NULL;
}
/* string to option value ------------------------------------------------------
* convert string to option value
* args   : opt_t  *opt      O  option
*          char   *str      I  option value string
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int str2opt(opt_t *opt, const char *str)
{
    switch (opt->format) {
        case 0: *(int    *)opt->var=atoi(str); break;
        case 1: *(double *)opt->var=atof(str); break;
        case 2: strcpy((char *)opt->var,str);  break;
        case 3: return str2enum(str,opt->comment,(int *)opt->var);
        default: return 0;
    }
    return 1;
}
/* option value to string ------------------------------------------------------
* convert option value to string
* args   : opt_t  *opt      I  option
*          char   *str      O  option value string
* return : length of output string
*-----------------------------------------------------------------------------*/
extern int opt2str(const opt_t *opt, char *str)
{
    char *p=str;
    
    trace(3,"opt2str : name=%s\n",opt->name);
    
    switch (opt->format) {
        case 0: p+=sprintf(p,"%d"   ,*(int   *)opt->var); break;
        case 1: p+=sprintf(p,"%.15g",*(double*)opt->var); break;
        case 2: p+=sprintf(p,"%s"   , (char  *)opt->var); break;
        case 3: p+=enum2str(p,opt->comment,*(int *)opt->var); break;
    }
    return (int)(p-str);
}
/* option to string -------------------------------------------------------------
* convert option to string (keyword=value # comment)
* args   : opt_t  *opt      I  option
*          char   *buff     O  option string
* return : length of output string
*-----------------------------------------------------------------------------*/
extern int opt2buf(const opt_t *opt, char *buff)
{
    char *p=buff;
    int n;
    
    trace(3,"opt2buf : name=%s\n",opt->name);
    
    p+=sprintf(p,"%-18s =",opt->name);
    p+=opt2str(opt,p);
    if (*opt->comment) {
        if ((n=(int)(buff+30-p))>0) p+=sprintf(p,"%*s",n,"");
        p+=sprintf(p," # (%s)",opt->comment);
    }
    return (int)(p-buff);
}
/* load options ----------------------------------------------------------------
* load options from file
* args   : char   *file     I  options file
*          opt_t  *opts     IO options table
*                              (terminated with table[i].name="")
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int loadopts(const char *file, opt_t *opts)
{
    FILE *fp;
    opt_t *opt;
    char buff[2048],*p;
    int n=0;
    
    trace(3,"loadopts: file=%s\n",file);
    
    if (!(fp=fopen(file,"r"))) {
        trace(1,"loadopts: options file open error (%s)\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        n++;
        chop(buff);
        
        if (buff[0]=='\0') continue;
        
        if (!(p=strstr(buff,"="))) {
            fprintf(stderr,"invalid option %s (%s:%d)\n",buff,file,n);
            continue;
        }
        *p++='\0';
        chop(buff);
        if (!(opt=searchopt(buff,opts))) continue;
        
        if (!str2opt(opt,p)) {
            fprintf(stderr,"invalid option value %s (%s:%d)\n",buff,file,n);
            continue;
        }
    }
    fclose(fp);
    
    return 1;
}
/* save options to file --------------------------------------------------------
* save options to file
* args   : char   *file     I  options file
*          char   *mode     I  write mode ("w":overwrite,"a":append);
*          char   *comment  I  header comment (NULL: no comment)
*          opt_t  *opts     I  options table
*                              (terminated with table[i].name="")
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int saveopts(const char *file, const char *mode, const char *comment,
                    const opt_t *opts)
{
    FILE *fp;
    char buff[2048];
    int i;
    
    trace(3,"saveopts: file=%s mode=%s\n",file,mode);
    
    if (!(fp=fopen(file,mode))) {
        trace(1,"saveopts: options file open error (%s)\n",file);
        return 0;
    }
    if (comment) fprintf(fp,"# %s\n\n",comment);
    
    for (i=0;*opts[i].name;i++) {
        opt2buf(opts+i,buff);
        fprintf(fp,"%s\n",buff);
    }
    fclose(fp);
    return 1;
}
/* system options buffer to options ------------------------------------------*/
static void buff2sysopts(void)
{
    double es[6],pos[3],*rr;
    char buff[1024],*p,*q,*id,*sep=NULL,sol_path[1024]={};
    int i,j,sat,ps;

    /* start time */
    sscanf(time_[0],"%lf/%lf/%lf %lf:%lf:%lf",es,es+1,es+2,es+3,es+4,es+5); prcopt_.ts=epoch2time(es);
    sscanf(time_[1],"%lf/%lf/%lf %lf:%lf:%lf",es,es+1,es+2,es+3,es+4,es+5); prcopt_.te=epoch2time(es);

    prcopt_.elmin     =elmask_    *D2R;
    prcopt_.elmaskar  =elmaskar_  *D2R;
    prcopt_.elmaskhold=elmaskhold_*D2R;
    
    /* receiver position */
    for (i=0;i<2;i++) {
        ps=i==0?prcopt_.rovpos:prcopt_.refpos;
        rr=i==0?prcopt_.ru:prcopt_.rb;
        
        if (ps==POSOPT_POS_LLH) { /* lat/lon/hgt */
            pos[0]=antpos_[i][0]*D2R;
            pos[1]=antpos_[i][1]*D2R;
            pos[2]=antpos_[i][2];
            pos2ecef(pos,rr);
        }
        else if (ps==POSOPT_POS_XYZ) { /* xyz-ecef */
            rr[0]=antpos_[i][0];
            rr[1]=antpos_[i][1];
            rr[2]=antpos_[i][2];
        }
    }
    /* excluded satellites */
    for (i=0;i<MAXSAT;i++) prcopt_.exsats[i]=0;
    if (exsats_[0]!='\0') {
        strcpy(buff,exsats_);
        for (p=strtok_r(buff," ",&q);p;p=strtok_r(NULL," ",&q)) {
            if (*p=='+') id=p+1; else id=p;
            if (!(sat=satid2no(id))) continue;
            prcopt_.exsats[sat-1]=*p=='+'?2:1;
        }
    }
    /* exclude BDS2 */
    if(1==prcopt_.bdsflag[0]) {
        strcpy(buff,BDS2);
        for (p=strtok_r(buff," ",&q);p;p=strtok_r(NULL," ",&q)) {
            if (*p=='+') id=p+1; else id=p;
            if (!(sat=satid2no(id))) continue;
            prcopt_.exsats[sat-1]=*p=='+'?2:1;
        }
    }
    /* exclude BDS3 */
    if(1==prcopt_.bdsflag[1]) {
        strcpy(buff,BDS3);
        for (p=strtok_r(buff," ",&q);p;p=strtok_r(NULL," ",&q)) {
            if (*p=='+') id=p+1; else id=p;
            if (!(sat=satid2no(id))) continue;
            prcopt_.exsats[sat-1]=*p=='+'?2:1;
        }
    }
    /* snrmask, currently, only triple-frequency is supported */
    for (i=0;i<3;i++) {
        for (j=0;j<9;j++) prcopt_.snrmask.mask[i][j]=0.0;
        strcpy(buff,snrmask_[i]);
        for (p=strtok_r(buff,",",&q),j=0;p&&j<9;p=strtok_r(NULL,",",&q)) {
            prcopt_.snrmask.mask[i][j++]=atof(p);
        }
    }
    /* frequency id */
    for (i=0;i<MAXSYS;i++) {
        for (j=0;j<MAXFREQ;j++) prcopt_.fre[i][j]=0;
        strcpy(buff,fre_[i]);
        for (p=strtok_r(buff,",",&q),j=0;p&&j<MAXFREQ;p=strtok_r(NULL,",",&q)) {
            prcopt_.fre[i][j++]=atoi(p);
        }
    }
    /* motion constraints options*/ 
    for (j=0;j<3;j++) prcopt_.constraint[j]=0.0;
    strcpy(buff,constraint_);
    for (p=strtok_r(buff,",",&q),j=0;p&&j<3;p=strtok_r(NULL,",",&q)) {
        prcopt_.constraint[j++]=atoi(p);
    }

    /* ins installation angle */
    for (j=0;j<3;j++) prcopt_.install_angle[j]=0.0;
    strcpy(buff,install_angle_);
    for (p=strtok_r(buff,",",&q),j=0;p&&j<3;p=strtok_r(NULL,",",&q)) {
        prcopt_.install_angle[j++]=atof(p)*D2R;
    }

    /* nhc lever */
    for (j=0;j<3;j++) prcopt_.lever_nhc[j]=0.0;
    strcpy(buff,lever_nhc_);
    for (p=strtok_r(buff,",",&q),j=0;p&&j<3;p=strtok_r(NULL,",",&q)) {
        prcopt_.lever_nhc[j++]=atof(p);
    }

    /* rotation angle */
    for (j=0;j<3;j++) prcopt_.rotation_angle[j]=0.0;
    strcpy(buff,rotation_angle_);
    for (p=strtok_r(buff,",",&q),j=0;p&&j<3;p=strtok_r(NULL,",",&q)) {
        prcopt_.rotation_angle[j++]=atof(p)*D2R;
    }

    /* init ins position (ecef frame[X,Y,Z] (m)/local frame [lat,lon,h] (deg,m) )*/
    for (j=0;j<3;j++) prcopt_.initpos[j]=0.0;
    strcpy(buff,initpose_[0]);
    for (p=strtok_r(buff,",",&q),j=0;p&&j<3;p=strtok_r(NULL,",",&q)) {
        pos[j++]=atof(p);
    }
    /* local frame [lat,lon,h] (deg,m) */
    if (fabs(pos[0])<=90&&fabs(pos[1])<=180) {
        for (j=0;j<2;j++) pos[j]*=D2R;
        matcpy(prcopt_.initpos,pos,3,1);
    }
    else { /* ecef frame[X,Y,Z] (m) */
        ecef2pos(pos,prcopt_.initpos);        
    }

    /* init ins velocity (n [E,N,U] (m/s) )*/
    for (j=0;j<3;j++) prcopt_.initvel[j]=0.0;
    strcpy(buff,initpose_[1]);
    for (p=strtok_r(buff,",",&q),j=0;p&&j<3;p=strtok_r(NULL,",",&q)) {
        prcopt_.initvel[j++]=atof(p);
    }  

    /* init ins attitude ([pitch,roll,yaw] (deg) )*/
    for (j=0;j<3;j++) prcopt_.initatt[j]=0.0;
    strcpy(buff,initpose_[2]);
    for (p=strtok_r(buff,",",&q),j=0;p&&j<3;p=strtok_r(NULL,",",&q)) {
        prcopt_.initatt[j++]=atof(p)*D2R;
    }
    if (prcopt_.initatt[2]>PI) prcopt_.initatt[2]-=2*PI; /* yaw is in range [-pi,pi] */
    /* counterclockwise is positive */
    prcopt_.initatt[2]=-prcopt_.initatt[2];

    /* init ins position std */
    for (j=0;j<3;j++) prcopt_.init_pos_unc[j]=0.0;
    strcpy(buff,initunc_[0]);
    for (p=strtok_r(buff,",",&q),j=0;p&&j<3;p=strtok_r(NULL,",",&q)) {
        prcopt_.init_pos_unc[j++]=atof(p);
    }

    /* init ins velocity std (m/s) */
    for (j=0;j<3;j++) prcopt_.init_vel_unc[j]=0.0;
    strcpy(buff,initunc_[1]);
    for (p=strtok_r(buff,",",&q),j=0;p&&j<3;p=strtok_r(NULL,",",&q)) {
        prcopt_.init_vel_unc[j++]=atof(p);
    }

    /* init ins attitude std (rad) */
    for (j=0;j<3;j++) prcopt_.init_att_unc[j]=0.0;
    strcpy(buff,initunc_[2]);
    for (p=strtok_r(buff,",",&q),j=0;p&&j<3;p=strtok_r(NULL,",",&q)) {
        prcopt_.init_att_unc[j++]=atof(p)*D2R;
    }

    /* lever */
    for (j=0;j<3;j++) prcopt_.lever[j]=0.0;
    strcpy(buff,lever_);
    for (p=strtok_r(buff,",",&q),j=0;p&&j<3;p=strtok_r(NULL,",",&q)) {
        prcopt_.lever[j++]=atof(p);
    }

    /* solution status opt*/
    strcpy(buff,stat_);
    for (p=strtok_r(buff,",",&q),j=0;p&&j<statopt;p=strtok_r(NULL,",",&q)) {
        solopt_.stato[j++]=atoi(p);
    }   

    /* path splicing */
    #if 0
    if (*filopt_.sol_path) {
        strcpy(sol_path,filopt_.sol_path);
        /* separator */
        sep=(strrchr(sol_path,'/'))?"/":"\\";
        strcat(sol_path,sep);

        /* path splicing */
        if (*filopt_.obs_u)    { strcpy(buff,sol_path); strcat(buff,filopt_.obs_u);    strcpy(filopt_.obs_u,buff); buff[0]='\0'; }
        if (*filopt_.obs_b)    { strcpy(buff,sol_path); strcat(buff,filopt_.obs_b);    strcpy(filopt_.obs_b,buff); buff[0]='\0'; }
        if (*filopt_.nav)      { strcpy(buff,sol_path); strcat(buff,filopt_.nav);      strcpy(filopt_.nav,buff); buff[0]='\0'; }
        if (*filopt_.sp3)      { strcpy(buff,sol_path); strcat(buff,filopt_.sp3);      strcpy(filopt_.sp3,buff); buff[0]='\0'; }
        if (*filopt_.clk)      { strcpy(buff,sol_path); strcat(buff,filopt_.clk);      strcpy(filopt_.clk,buff); buff[0]='\0'; }
        if (*filopt_.imu)      { strcpy(buff,sol_path); strcat(buff,filopt_.imu);      strcpy(filopt_.imu,buff); buff[0]='\0'; }
        if (*filopt_.pos)      { strcpy(buff,sol_path); strcat(buff,filopt_.pos);      strcpy(filopt_.pos,buff); buff[0]='\0'; }
        if (*filopt_.antp)     { strcpy(buff,sol_path); strcat(buff,filopt_.antp);     strcpy(filopt_.antp,buff); buff[0]='\0'; }
        if (*filopt_.mgex_dcb) { strcpy(buff,sol_path); strcat(buff,filopt_.mgex_dcb); strcpy(filopt_.mgex_dcb,buff); buff[0]='\0'; }
        if (*filopt_.stapos)   { strcpy(buff,sol_path); strcat(buff,filopt_.stapos);   strcpy(filopt_.stapos,buff); buff[0]='\0'; }
        if (*filopt_.blq)      { strcpy(buff,sol_path); strcat(buff,filopt_.blq);      strcpy(filopt_.blq,buff); buff[0]='\0'; }
        if (*filopt_.trace)    { strcpy(buff,sol_path); strcat(buff,filopt_.trace);    strcpy(filopt_.trace,buff); buff[0]='\0'; }
    }
    #endif
    /* Guard number of frequencies */
    if (prcopt_.nf>NFREQ) {
        fprintf(stderr,"Number of frequencies %d limited to %d, rebuild with NFREQ=%d\n",prcopt_.nf, NFREQ, prcopt_.nf);
        prcopt_.nf=NFREQ;
    }

}
/* options to system options buffer ------------------------------------------*/
static void sysopts2buff(void)
{
    double pos[3],*rr;
    char id[8],*p;
    int i,j,sat,ps;
    
    elmask_    =prcopt_.elmin     *R2D;
    elmaskar_  =prcopt_.elmaskar  *R2D;
    elmaskhold_=prcopt_.elmaskhold*R2D;
    
    for (i=0;i<2;i++) {
        ps=i==0?prcopt_.rovpos:prcopt_.refpos;
        rr=i==0?prcopt_.ru:prcopt_.rb;
        
        if (ps==POSOPT_POS_LLH) {
            ecef2pos(rr,pos);
            antpos_[i][0]=pos[0]*R2D;
            antpos_[i][1]=pos[1]*R2D;
            antpos_[i][2]=pos[2];
        } else if (ps==POSOPT_POS_XYZ) {
            antpos_[i][0] = rr[0];
            antpos_[i][1] = rr[1];
            antpos_[i][2] = rr[2];
        }
    }
    /* excluded satellites */
    exsats_[0]='\0';
    for (sat=1,p=exsats_;sat<=MAXSAT&&p-exsats_<(int)sizeof(exsats_)-32;sat++) {
        if (prcopt_.exsats[sat-1]) {
            satno2id(sat,id);
            p+=sprintf(p,"%s%s%s",p==exsats_?"":" ",
                       prcopt_.exsats[sat-1]==2?"+":"",id);
        }
    }
    /* snrmask */
    for (i=0;i<NFREQ;i++) {
        snrmask_[i][0]='\0';
        p=snrmask_[i];
        for (j=0;j<9;j++) {
            p+=sprintf(p,"%s%.0f",j>0?",":"",prcopt_.snrmask.mask[i][j]);
        }
    }
    /* number of frequency (4:L1+L5) TODO ???? */
    /*if (prcopt_.nf==3&&prcopt_.freqopt==1) {
        prcopt_.nf=4;
        prcopt_.freqopt=0;
    }*/
}
/* reset system options to default ---------------------------------------------
* reset system options to default
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
extern void resetsysopts(void)
{
    int i,j;
    
    trace(3,"resetsysopts:\n");
    
    prcopt_=prcopt_default;
    solopt_=solopt_default;
    filopt_.obs_u  [0]='\0';
    filopt_.obs_b  [0]='\0';
    filopt_.nav    [0]='\0';
    filopt_.sp3    [0]='\0';
    filopt_.clk    [0]='\0';
    filopt_.imu    [0]='\0';
    filopt_.pos    [0]='\0';
    filopt_.mgex_dcb[0]='\0';
    filopt_.antp   [0]='\0';
    filopt_.stapos [0]='\0';
    filopt_.geoid  [0]='\0';
    filopt_.dcb    [0]='\0';
    filopt_.eop    [0]='\0';
    filopt_.blq    [0]='\0';
    elmask_=15.0;
    elmaskar_=0.0;
    elmaskhold_=0.0;
    for (i=0;i<2;i++) for (j=0;j<3;j++) {
        antpos_[i][j]=0.0;
    }
    exsats_[0] ='\0';
}
/* get system options ----------------------------------------------------------
* get system options
* args   : prcopt_t *popt   IO processing options (NULL: no output)
*          solopt_t *sopt   IO solution options   (NULL: no output)
*          folopt_t *fopt   IO file options       (NULL: no output)
* return : none
* notes  : to load system options, use loadopts() before calling the function
*-----------------------------------------------------------------------------*/
extern void getsysopts(prcopt_t *popt, solopt_t *sopt, filopt_t *fopt)
{
    trace(3,"getsysopts:\n");
    
    buff2sysopts();
    if (popt) *popt=prcopt_;
    if (sopt) *sopt=solopt_;
    if (fopt) *fopt=filopt_;
}
/* set system options ----------------------------------------------------------
* set system options
* args   : prcopt_t *prcopt I  processing options (NULL: default)
*          solopt_t *solopt I  solution options   (NULL: default)
*          filopt_t *filopt I  file options       (NULL: default)
* return : none
* notes  : to save system options, use saveopts() after calling the function
*-----------------------------------------------------------------------------*/
extern void setsysopts(const prcopt_t *prcopt, const solopt_t *solopt,
                       const filopt_t *filopt)
{
    trace(3,"setsysopts:\n");
    
    resetsysopts();
    if (prcopt) prcopt_=*prcopt;
    if (solopt) solopt_=*solopt;
    if (filopt) filopt_=*filopt;
    sysopts2buff();
}