/*------------------------------------------------------------------------------
* convkml.c : google earth kml converter
*
*          Copyright (C) 2007-2017 by T.TAKASU, All rights reserved.
*
* references :
*     [1] Open Geospatial Consortium Inc., OGC 07-147r2, OGC(R) KML, 2008-04-14
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/20  1.0  new
*           2007/03/15  1.1  modify color sequence
*           2007/04/03  1.2  add geodetic height option
*                            support input of NMEA GGA sentence
*                            delete altitude info for track
*                            add time stamp option
*                            separate readsol.c file
*           2009/01/19  1.3  fix bug on display mark with by-q-flag option
*           2010/05/10  1.4  support api readsolt() change
*           2010/08/14  1.5  fix bug on readsolt() (2.4.0_p3)
*           2017/06/10  1.6  support wild-card in input file
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants -----------------------------------------------------------------*/

#define SIZP     0.2            /* mark size of rover positions */
#define SIZR     0.3            /* mark size of reference position */
#define TINT     60.0           /* time label interval (sec) */


static ref_t ref={0};
static err_t err={0};

static const char *head1="<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
static const char *head2="<kml xmlns=\"http://earth.google.com/kml/2.1\">";
static const char *mark="http://maps.google.com/mapfiles/kml/pal2/icon18.png";

/* output track --------------------------------------------------------------*/
static void outtrack(FILE *f, const solbuf_t *solbuf, const char *color,
                     int outalt, int outtime)
{
    double pos[3];
    int i;
    
    fprintf(f,"<Placemark>\n");
    fprintf(f,"<name>Rover Track</name>\n");
    fprintf(f,"<Style>\n");
    fprintf(f,"<LineStyle>\n");
    fprintf(f,"<color>%s</color>\n",color);
    fprintf(f,"</LineStyle>\n");
    fprintf(f,"</Style>\n");
    fprintf(f,"<LineString>\n");
    if (outalt) fprintf(f,"<altitudeMode>absolute</altitudeMode>\n");
    fprintf(f,"<coordinates>\n");
    for (i=0;i<solbuf->n;i++) {
        ecef2pos(solbuf->data[i].rr,pos);
        if      (outalt==0) pos[2]=0.0;
        else if (outalt==2) pos[2]-=geoidh(pos);
        fprintf(f,"%13.9f,%12.9f,%5.3f\n",pos[1]*R2D,pos[0]*R2D,pos[2]);
    }
    fprintf(f,"</coordinates>\n");
    fprintf(f,"</LineString>\n");
    fprintf(f,"</Placemark>\n");
}
/* output point --------------------------------------------------------------*/
static void outpoint(FILE *fp, gtime_t time, const solbuf_t *solbuf, const double *pos,
                     const char *label, int style, int outalt, int outtime)
{
    const err_t *err=NULL;
    const sol_t *sol=NULL;
    double ep[6],alt=0.0;
    char str[256]="", name_str[256]="", desc_str[4096]="";
    int i;
    double pos_ref[3], pos_deg[3], vel_enu[3], att_deg[3], pos_sig[3], vel_sig[3], att_sig[3];
    double pos_err[3],vel_err[3]; 
    
    fprintf(fp,"<Placemark>\n");

    if (solbuf) {
        sprintf(name_str,"%.3f",time2gpst(time,NULL));
        fprintf(fp,"<name>%s</name>\n",name_str);
    } else if (*label) {
        fprintf(fp,"<name>%s</name>\n",label);
    }

    fprintf(fp,"<Snippet maxLines=\"0\"></Snippet>\n");

    if (solbuf) {
        err=solbuf->err;
        /* convert position to degree-minute-second format */
        pos_deg[0]=pos[0]*R2D;  /* latitude */
        pos_deg[1]=pos[1]*R2D;  /* longitude */
        pos_deg[2]=pos[2];      /* elevation */
        
        /* find corresponding sol data */
        for (i=0;i<solbuf->n;i++) {
            if (fabs(timediff(solbuf->data[i].time,time))<0.01) {
                sol=&solbuf->data[i];
                break;
            }
        }
        
        if (sol) {
            /* compute velocity of ENU frame */
            if (norm(sol->rr+3,3)>0) {
                ecef2pos(pos,pos_ref);
                ecef2enu(pos_ref,sol->rr+3,vel_enu);
            } else {
                vel_enu[0]=vel_enu[1]=vel_enu[2]=0.0;
            }
            
            /* attitude */
            if (sol->att[2]>0) {
                att_deg[0]=sol->att[0];  /* pitch */
                att_deg[1]=sol->att[1];  /* roll */
                att_deg[2]=sol->att[2];  /* yaw */
            } else {
                att_deg[0]=att_deg[1]=att_deg[2]=0.0;
            }

            /* error analysis */
            ecef2enu(pos,err->data[i].pos,pos_err);
            ecef2enu(pos,err->data[i].vel,vel_err);
            
            /* position accuracy */
            pos_sig[0]=sqrt(sol->qr[0]);  /* E */
            pos_sig[1]=sqrt(sol->qr[1]);  /* N */
            pos_sig[2]=sqrt(sol->qr[2]);  /* U */
            
            /* velocity accuracy */
            if (norm(sol->qv,6)>0) {
                vel_sig[0]=sqrt(sol->qv[0]);  /* E */
                vel_sig[1]=sqrt(sol->qv[1]);  /* N */
                vel_sig[2]=sqrt(sol->qv[2]);  /* U */
            } else {
                vel_sig[0]=vel_sig[1]=vel_sig[2]=0.0;
            }
            
            /* attitude accuracy */
            if (norm(sol->qa,6)>0) {
                att_sig[0]=sqrt(sol->qa[0]);  /* pitch */
                att_sig[1]=sqrt(sol->qa[1]);  /* roll */
                att_sig[2]=sqrt(sol->qa[2]);  /* yaw */
            } else {
                att_sig[0]=att_sig[1]=att_sig[2]=0.0;
            }
            
            /* generate description table */
            time2epoch(time,ep);
            sprintf(desc_str,
                "<![CDATA[<B>Epoch:%.0f - Q%d</B><BR><BR>\n"
                "<TABLE border=\"1\" width=\"100\" Align=\"center\">\n"
                "<TR ALIGN=RIGHT>\n"
                "<TR ALIGN=RIGHT><TD ALIGN=LEFT>Time</TD><TD>%04.0f/%02.0f/%02.0f</TD><TD>%02.0f:%02.0f:%05.2f</TD><TD>%.0f</TD><TD>%.2f</TD></TR>\n"
                "<TR ALIGN=RIGHT><TD ALIGN=LEFT>Position</TD><TD>%.0f %.0f %.6f</TD><TD>%.0f %.0f %.6f</TD><TD>%.3f</TD><TD>(DMS,m)</TD></TR>\n"
                "<TR ALIGN=RIGHT><TD ALIGN=LEFT>Velocity</TD><TD>%.3f</TD><TD>%.3f</TD><TD>%.3f</TD><TD>(m/s)</TD></TR>\n"
                "<TR ALIGN=RIGHT><TD ALIGN=LEFT>Attitude</TD><TD>%.5f</TD><TD>%.5f</TD><TD>%.5f</TD><TD>(deg)</TD></TR>\n"
                "<TR ALIGN=RIGHT><TD ALIGN=LEFT>Poserr</TD><TD>%.3f</TD><TD>%.3f</TD><TD>%.3f</TD><TD>(m)</TD></TR>\n"
                "<TR ALIGN=RIGHT><TD ALIGN=LEFT>Velerr</TD><TD>%.3f</TD><TD>%.3f</TD><TD>%.3f</TD><TD>(m/s)</TD></TR>\n"
                "<TR ALIGN=RIGHT><TD ALIGN=LEFT>Atterr</TD><TD>%.5f</TD><TD>%.5f</TD><TD>%.5f</TD><TD>(deg)</TD></TR>\n"
                "<TR ALIGN=RIGHT><TD ALIGN=LEFT>PosSig</TD><TD>%.3f</TD><TD>%.3f</TD><TD>%.3f</TD><TD>(m)</TD></TR>\n"
                "<TR ALIGN=RIGHT><TD ALIGN=LEFT>VelSig</TD><TD>%.3f</TD><TD>%.3f</TD><TD>%.3f</TD><TD>(m/s)</TD></TR>\n"
                "<TR ALIGN=RIGHT><TD ALIGN=LEFT>AttSig</TD><TD>%.5f</TD><TD>%.5f</TD><TD>%.5f</TD><TD>(deg)</TD></TR>\n"
                "<TR ALIGN=RIGHT><TD ALIGN=LEFT>Quality</TD><TD>Q%d</TD><TD>%s</TD><TD>%.2f</TD><TD>%s</TD></TR>\n"
                "<TR ALIGN=RIGHT><TD ALIGN=LEFT>PDOP</TD><TD>%.2f</TD><TD>Age</TD><TD>%.4f</TD><TD>(s)</TD></TR>\n"
                "<TR ALIGN=RIGHT><TD ALIGN=LEFT>NSAT</TD><TD>%2d</TD></TR>\n"
                "</TABLE>]]>",
                time2gpst(time,NULL), sol->stat,
                ep[0], ep[1], ep[2], ep[3], ep[4], ep[5], time2gpst(time,NULL), time2doy(time),
                floor(fabs(pos_deg[0])), floor(fmod(fabs(pos_deg[0])*60,60)), fmod(fabs(pos_deg[0])*3600,60),
                floor(fabs(pos_deg[1])), floor(fmod(fabs(pos_deg[1])*60,60)), fmod(fabs(pos_deg[1])*3600,60),
                pos_deg[2],
                vel_enu[0], vel_enu[1], vel_enu[2],
                att_deg[0], att_deg[1], att_deg[2],
                norm(pos_err,2),pos_err[2],norm(pos_err,3),
                norm(vel_err,2),vel_err[2],norm(vel_err,3),
                err->data[i].att[0],err->data[i].att[1],err->data[i].att[2],
                pos_sig[0], pos_sig[1], pos_sig[2],
                vel_sig[0], vel_sig[1], vel_sig[2],
                att_sig[0], att_sig[1], att_sig[2],
                sol->stat, (SOLQ_SINGLE==sol->stat)?"SPP":((SOLQ_DGPS==sol->stat)?"DGNSS":((SOLQ_FLOAT==sol->stat)?"AMB_FLOAT":((SOLQ_FIX==sol->stat)?"AMB_FIX":((SOLQ_PPP==sol->stat)?"PPP":"NONE")))), sol->ratio,
                (sol->iFlag==1)?"IMU_ZUPT":"IMU_NONE",
                sol->dop[1], sol->age,
                sol->ns);
            fprintf(fp,"<description>%s</description>\n",desc_str);
        }
    }

    fprintf(fp,"<styleUrl>#P%d</styleUrl>\n",style);

    if (outtime) {
        if      (outtime==2) time=gpst2utc(time);
        else if (outtime==3) time=timeadd(gpst2utc(time),9*3600.0);
        time2epoch(time,ep);
        sprintf(str,"%04.0f-%02.0f-%02.0fT%02.0f:%02.0f:%05.2fZ",
                ep[0],ep[1],ep[2],ep[3],ep[4],ep[5]);
        fprintf(fp,"<TimeStamp><when>%s</when></TimeStamp>\n",str);
    }

    fprintf(fp,"<Point>\n");
    if (outalt) {
        fprintf(fp,"<extrude>1</extrude>\n");
        fprintf(fp,"<altitudeMode>absolute</altitudeMode>\n");
        alt=pos[2]-(outalt==2?geoidh(pos):0.0);
    }
    fprintf(fp,"<coordinates>%13.9f,%12.9f,%5.3f</coordinates>\n",pos[1]*R2D,
            pos[0]*R2D,alt);
    fprintf(fp,"</Point>\n");

    fprintf(fp,"</Placemark>\n");
}

/* save kml file -------------------------------------------------------------*/
static int savekml(const char *file, const solbuf_t *solbuf, int tcolor,
                   int pcolor, int outalt, int outtime)
{
    FILE *fp;
    double pos[3];
    int i,qcolor[]={0,1,2,5,4,3,0};
    char *color[]={
        "ffffffff","ff008800","ff00aaff","ff0000ff","ff00ffff","ffff00ff"
    };
    if (!(fp=fopen(file,"w"))) {
        fprintf(stderr,"file open error : %s\n",file);
        return 0;
    }
    fprintf(fp,"%s\n%s\n",head1,head2);
    fprintf(fp,"<Document>\n");
    for (i=0;i<6;i++) {
        fprintf(fp,"<Style id=\"P%d\">\n",i);
        fprintf(fp,"  <IconStyle>\n");
        fprintf(fp,"    <color>%s</color>\n",color[i]);
        fprintf(fp,"    <scale>%.1f</scale>\n",i==0?SIZR:SIZP);
        fprintf(fp,"    <Icon><href>%s</href></Icon>\n",mark);
        fprintf(fp,"  </IconStyle>\n");
        fprintf(fp,"  <LabelStyle><scale>%d</scale></LabelStyle>\n",0); /* Control whether Google Earth displays coordinate point labels (generally set to off). */
        fprintf(fp,"  <BalloonStyle>\n");
        fprintf(fp,"     <text><![CDATA[<b><font color=%s size=%d>$[name]</font></b><br>$[description]</font><br/>]]></text>\n","#cc0000",3);
        fprintf(fp,"     <bgColor>%s</bgColor>\n","ffd5f3fa");
        fprintf(fp,"  </BalloonStyle>\n");
        fprintf(fp,"</Style>\n");
    }
    if (tcolor>0) {
        outtrack(fp,solbuf,color[tcolor-1],outalt,outtime);
    }
    if (pcolor>0) {
        fprintf(fp,"<Folder>\n");
        fprintf(fp,"  <name>Rover Position</name>\n");
        for (i=0;i<solbuf->n;i++) {
            ecef2pos(solbuf->data[i].rr,pos);
            outpoint(fp,solbuf->data[i].time,solbuf,pos,"",
                     pcolor==5?qcolor[solbuf->data[i].stat]:pcolor-1,outalt,outtime);
        }
        fprintf(fp,"</Folder>\n");
    }
    if (norm(solbuf->rb,3)>0.0) {
        ecef2pos(solbuf->rb,pos);
        outpoint(fp,solbuf->data[0].time,NULL,pos,"Reference Position",0,outalt,0);
    }
    fprintf(fp,"</Document>\n");
    fprintf(fp,"</kml>\n");
    fclose(fp);
    return 1;
}
/* convert to google earth kml file --------------------------------------------
* convert solutions to google earth kml file
* args   : char   *infile   I   input solutions file (wild-card (*) is expanded)
*          char   *outfile  I   output google earth kml file ("":<infile>.kml)
*          gtime_t ts,te    I   start/end time (gpst)
*          int    tint      I   time interval (s) (0.0:all)
*          int    qflg      I   quality flag (0:all)
*          double *offset   I   add offset {east,north,up} (m)
*          int    tcolor    I   track color
*                               (0:none,1:white,2:green,3:orange,4:red,5:yellow)
*          int    pcolor    I   point color
*                               (0:none,1:white,2:green,3:orange,4:red,5:by qflag)
*          int    outalt    I   output altitude (0:off,1:elipsoidal,2:geodetic)
*          int    outtime   I   output time (0:off,1:gpst,2:utc,3:jst)
* return : status (0:ok,-1:file read,-2:file format,-3:no data,-4:file write)
* notes  : see ref [1] for google earth kml file format
*-----------------------------------------------------------------------------*/
extern int convkml(const char *infile, const char *refile, const char *outfile, gtime_t ts,
                   gtime_t te, solopt_t *sopt, double tint, int qflg, double *offset,
                   int tcolor, int pcolor, int outerr, int outalt, int outtime)
{
    solbuf_t solbuf={0};
    double rr[3]={0},pos[3],dr[3];
    int i,j,nfile,stat,stat2,r;
    char *p,file[1024],errfile[1024],*files[MAXEXFILE]={0};
    
    trace(3,"convkml : infile=%s outfile=%s\n",infile,outfile);
    
    /* expand wild-card of infile */
    for (i=0;i<MAXEXFILE;i++) {
        if (!(files[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(files[i]);
            return -4;
        }
    }
    if ((nfile=expath(infile,files,MAXEXFILE))<=0) {
        for (i=0;i<MAXEXFILE;i++) free(files[i]);
        return -3;
    }
    if (!*outfile) {
        if ((p=strrchr(infile,'.'))) {
            strncpy(file,infile,p-infile);
            strcpy(file+(p-infile),".kml");
        }
        else sprintf(file,"%s.kml",infile);
    }
    else strcpy(file,outfile);

    /* output error file name */
    if (p=strrchr(infile,'.')) {
        strncpy(errfile,infile,p-infile);
        strcpy(errfile+(p-infile),".err");
    }
    else sprintf(errfile,"%s.err",infile);
    
    /* read solution file */
    stat=readsolt((const char **)files,nfile,ts,te,tint,qflg,&solbuf);

    /* read reference file */
    if (refile&&outerr) {
        if (readref(refile,&ref)<0) {
            trace(7,"convkml: readref error %s\n",refile);
            return -1;
        }  
        /* navigation error analysis */
        err_analysis(&ref,sopt,&solbuf,&err);
        solbuf.err=&err; /* assign err to solbuf */
        /* save error file */
        saverr(&solbuf,&err,errfile);
    }
    
    for (i=0;i<MAXEXFILE;i++) free(files[i]);
    
    if (!stat) {
        return -1;
    }
    /* mean position */
    for (i=0;i<3;i++) {
        for (j=0;j<solbuf.n;j++) rr[i]+=solbuf.data[j].rr[i];
        rr[i]/=solbuf.n;
    }
    /* add offset */
    ecef2pos(rr,pos);
    enu2ecef(pos,offset,dr);
    for (i=0;i<solbuf.n;i++) {
        for (j=0;j<3;j++) solbuf.data[i].rr[j]+=dr[j];
    }
    if (norm(solbuf.rb,3)>0.0) {
        for (i=0;i<3;i++) solbuf.rb[i]+=dr[i];
    }
    /* save kml file */
    r=savekml(file,&solbuf,tcolor,pcolor,outalt,outtime)?0:-4;
    freesolbuf(&solbuf);
    return r;
}
