/*------------------------------------------------------------------------------
*ins.c : ins common functions
 *-----------------------------------------------------------------------------*/

#include "rtklib.h"

static imu_t imus={0};          /* imu data */

extern void pos2sol(pos_t pos, sol_t *sol, int ipos)
{
    posd_t *posd=&pos.data[ipos];
    int i;

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

/* copy imu data --------------------------------------------------------------
* copy imu data to imu[] from imus.data[iimu]
*
*args   : const prcopt_t *popt  I   processing options
*         imud_t *imu            O   imu data
*         imu_t imus             I   imu data
*         int iimu               I   index of imu data
*         const int nn           I   number of imu data to copy
*return : none
*-----------------------------------------------------------------------------*/
extern void imucpy(const prcopt_t *popt, imud_t *imu, imu_t imus, int iimu, const int nn)
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

/* transform attitude to direction cosine matirx (DCM) --------------------------
*
*args   : double *att      I   attitude {pitch,roll,yaw} (rad)
*         double *Cnb      O   direction cosine matirx form b frame to n frame
*return : none
*-------------------------------------------------------------------------------*/
extern void att2Cnb(const double *att, double *Cnb)
{
    double sp=sin(att[0]),sr=sin(att[1]),sy=sin(att[2]);
    double cp=cos(att[0]),cr=cos(att[1]),cy=cos(att[2]);

    Cnb[0]= cr*cy-sp*sr*sy; Cnb[1]=-cp*sy; Cnb[2]= sr*cy+sp*cr*sy;
    Cnb[3]= cr*sy+sp*sr*cy; Cnb[4]= cp*cy; Cnb[5]= sr*sy-sp*cr*cy;
    Cnb[6]=-cp*sr;          Cnb[7]= sp;    Cnb[8]= cp*cr;
}

/* transform direction cosine matirx(DCM) to attitude --------------------------
*
*args   : double *Cnb      I   direction cosine matirx form b frame to n frame
*         double *att      O   attitude {pitch [-pi/2,pi/2],roll [-pi,pi],yaw [0,2*pi]} (rad)
*return : none
*-------------------------------------------------------------------------------*/
extern void Cnb2att(const double *Cnb, double *att)
{
    double roll=0.0,pitch=0.0,yaw=0.0;

    pitch=atan(Cnb[7]/sqrt(Cnb[6]*Cnb[6]+Cnb[8]*Cnb[8]));
    roll =-atan2(Cnb[6],Cnb[8]);
    yaw  =-atan2(Cnb[1],Cnb[4]);

    att[0]=pitch; att[1]=roll; att[2]=yaw;
}

/* transform attitude to attitude Quaternion qnb --------------------------
*
*args   : double *att      I   attitude {pitch,roll,yaw} (rad)
*         double *qnb      O   quaternion form b frame to n frame
*return : none
*-------------------------------------------------------------------------------*/
extern void att2qnb(const double *att, double *qnb)
{
    double s2p=sin(att[0]/2),s2r=sin(att[1]/2),s2y=sin(att[2]/2);
    double c2p=cos(att[0]/2),c2r=cos(att[1]/2),c2y=cos(att[2]/2);

    qnb[0]= c2p*c2r*c2y-s2p*s2r*s2y; /* qnb[0] = q0 */
    qnb[1]= s2p*c2r*c2y-c2p*s2r*s2y; /* qnb[1] = q1 */
    qnb[2]= s2p*c2r*s2y+c2p*s2r*c2y; /* qnb[2] = q2 */
    qnb[3]= c2p*c2r*s2y+s2p*s2r*c2y; /* qnb[3] = q3 */
}

/* transform attitude Quaternion qnb to DCM --------------------------
*
*args   : double *qnb      I   quaternion form b frame to n frame
*         double *Cnb      O   direction cosine matirx form b frame to n frame
*return : none
*-------------------------------------------------------------------------------*/
extern void qnb2Cnb(const double *qnb, double *Cnb)
{
    double q0=qnb[0],q1=qnb[1],q2=qnb[2],q3=qnb[3];

    Cnb[0]= q0*q0+q1*q1-q2*q2-q3*q3; Cnb[1]= 2.0*(q1*q2-q0*q3);       Cnb[2]= 2.0*(q1*q3+q0*q2);
    Cnb[3]= 2.0*(q1*q2+q0*q3);       Cnb[4]= q0*q0-q1*q1+q2*q2-q3*q3; Cnb[5]= 2.0*(q2*q3-q0*q1);
    Cnb[6]= 2.0*(q1*q3-q0*q2);       Cnb[7]= 2.0*(q2*q3+q0*q1);       Cnb[8]= q0*q0-q1*q1-q2*q2+q3*q3;
}

/* quaternion normalization --------------------------
*
*args   : double *qnb_      IO   quaternion form b frame to n frame
*return : none
*-------------------------------------------------------------------------------*/
extern void qnbnorm(double *qnb_)
{
    double qnb_m=norm(qnb_,4);

    if (qnb_m<=0.0) {
        qnb_[0]=1.0; qnb_[1]=0.0; qnb_[2]=0.0; qnb_[3]=0.0;
    }
    else {
        qnb_[0]=qnb_[0]/qnb_m; qnb_[1]=qnb_[1]/qnb_m; qnb_[2]=qnb_[2]/qnb_m; qnb_[3]=qnb_[3]/qnb_m;
    }
}

/* quaternion conjugation --------------------------
*
*args   : double *q1      I   quaternion 1
*         double *q2      O   quaternion q2=q1*
*return : none
*-------------------------------------------------------------------------------*/
extern void quatconj(const double *q1, double *q2)
{
    q2[0]=q1[0]; q2[1]=-q1[1]; q2[2]=-q1[2]; q2[3]=-q1[3];
}

/* two quaternion multiplication --------------------------
*
*args   : double *q1      I   quaternion 1
*         double *q2      I   quaternion 2
*         double *q      O   quaternion q=q1°q2
*return : none
*-------------------------------------------------------------------------------*/
extern void quatmul(const double *q1, const double *q2, double *q)
{
    q[0]=q1[0]*q2[0]-q1[1]*q2[1]-q1[2]*q2[2]-q1[3]*q2[3];
    q[1]=q1[0]*q2[1]+q1[1]*q2[0]+q1[2]*q2[3]-q1[3]*q2[2];
    q[2]=q1[0]*q2[2]-q1[1]*q2[3]+q1[2]*q2[0]+q1[3]*q2[1];
    q[3]=q1[0]*q2[3]+q1[1]*q2[2]-q1[2]*q2[1]+q1[3]*q2[0];
}

/* three quaternion multiplication --------------------------
*
*args   : double *q1      I   quaternion 1
*         double *q2      I   quaternion 2
*         double *q3      I   quaternion 3
*         double *q      O   quaternion q=q1°q2°q3
*return : none
*-------------------------------------------------------------------------------*/
extern void quatmul3(const double *q1, const double *q2, const double *q3, double *q)
{
    double q12[4];

    quatmul(q1,q2,q12); /* q12=q1°q2 */
    quatmul(q12,q3,q);  /* q=q12°q3 */
}

/* three-dimensional vector transformation based on quaternion --------------------------
*
*args   : double *f       I   vector coefficient
*         double *q       I   quaternion
*         double *vi      I   vector (3x1)
*         double *vo      O   vector (3x1) vo=q°vi°q*
*return : none
*-------------------------------------------------------------------------------*/
extern void quatmulv(const double f, const double *q, const double *vi, double *vo)
{
    double v[4],q2[4],qv1[4]={0.0},qv2[4]={0.0};

    /* v = [0, vi] */
    v[0]=0.0;  v[1]=vi[0]; v[2]=vi[1]; v[3]=vi[2];
    /* q2 = [q0, -q1,-q2,-q3] */
    q2[0]=q[0]; q2[1]=-q[1]; q2[2]=-q[2]; q2[3]=-q[3];

    quatmul(q,v,qv1);    /* qv1=q°vi */
    quatmul(qv1,q2,qv2); /* qv2=q°vi°q* */

    vo[0]=f*qv2[1]; vo[1]=f*qv2[2]; vo[2]=f*qv2[3]; /* qv2 = [0,vo] */
}

/* DCM Transpose --------------------------
*
*args   : double *Cxy    I   DCM y to x (3x3)
*         double *Cyx    O   DCM x to y (3x3)
*return : none
*-------------------------------------------------------------------------------*/
extern void DCMT(const double *Cxy, double *Cyx)
{
    Cyx[0]=Cxy[0]; Cyx[1]=Cxy[3]; Cyx[2]=Cxy[6];
    Cyx[3]=Cxy[1]; Cyx[4]=Cxy[4]; Cyx[5]=Cxy[7];
    Cyx[6]=Cxy[2]; Cyx[7]=Cxy[5]; Cyx[8]=Cxy[8];
}

/* Skew-symmetric matrix of the vector --------------------------
*
*args   : double  f      I   matrix coefficient
*         double *v      I   vector (3x1)
*         double *vx     O   skew symmetric matrix (3x3)
*return : none
*-------------------------------------------------------------------------------*/
extern void vskew(double f, const double *v, double *vx)
{
    int i;

    vx[0]=0;     vx[1]=-v[2]; vx[2]= v[1];
    vx[3]= v[2]; vx[4]=0;     vx[5]=-v[0];
    vx[6]=-v[1]; vx[7]= v[0]; vx[8]=0;

    for (i=0;i<9;i++)
    {
        vx[i]*=f;
    }
}

/* Multiplying two matrixs --------------------------
*
*args   : double  f        I   matrix coefficient
*         double *mat1     I   input  matrix1 (3x3)
*         double *mat2     I   input  matrix2 (3x3)
          double *mat      O   output matrix  (3x3)
*return : none
*-------------------------------------------------------------------------------*/
extern void Mat3mul2(double f, const double *mat1, const double *mat2, double *mat)
{
    int i;

    mat[0]=mat1[0]*mat2[0]+mat1[1]*mat2[3]+mat1[2]*mat2[6];
    mat[1]=mat1[0]*mat2[1]+mat1[1]*mat2[4]+mat1[2]*mat2[7];
    mat[2]=mat1[0]*mat2[2]+mat1[1]*mat2[5]+mat1[2]*mat2[8];

    mat[3]=mat1[3]*mat2[0]+mat1[4]*mat2[3]+mat1[5]*mat2[6];
    mat[4]=mat1[3]*mat2[1]+mat1[4]*mat2[4]+mat1[5]*mat2[7];
    mat[5]=mat1[3]*mat2[2]+mat1[4]*mat2[5]+mat1[5]*mat2[8];

    mat[6]=mat1[6]*mat2[0]+mat1[7]*mat2[3]+mat1[8]*mat2[6];
    mat[7]=mat1[6]*mat2[1]+mat1[7]*mat2[4]+mat1[8]*mat2[7];
    mat[8]=mat1[6]*mat2[2]+mat1[7]*mat2[5]+mat1[8]*mat2[8];    

    for (i=0;i<9;i++)
    {
        mat[i]*=f;
    }
}

/* transform rotation vector to direction cosine matirx(DCM) --------------------------
*
*args   : double  f       I   rotation vector coefficient
*         double *rv      I   rotation vector (3x1)
*         double *DCM     O   DCM (3x3)
*return : none
*-------------------------------------------------------------------------------*/
extern void rv2DCM(double f, const double *rv, double *dcm)
{
    int i;
    double rv1[3],rv_m,temp1[9],temp2[9],*I3;

    I3=eye(3);
    for (i=0;i<3;i++)
    {
        rv1[i]=f*rv[i];
    }
    rv_m=norm(rv1,3);
    vskew(1.0,rv1,temp1); Mat3mul2(1.0,temp1,temp1,temp2);

    for (i=0;i<9;i++)
    {
        dcm[i]=I3[i]+sin(rv_m)/rv_m*temp1[i]+(1.0-cos(rv_m))/(rv_m*rv_m)*temp2[i];
    }

    free(I3); 
}

/* transform rotation vector to attitude quaternion --------------------------
*
*args   : double  f       I   rotation vector coefficient
*         double *rv      I   rotation vector (3x1)
*         double *q       O   quaternion (4x1)
*return : none
*-------------------------------------------------------------------------------*/
extern void rv2quat(double f, const double *rv, double *q)
{
    int i;
    double rv1[3],rv_m,n_rv[3]={0.0};

    for (i=0;i<3;i++)
    {
        rv1[i]=f*rv[i];
    }
    rv_m=norm(rv1,3);

    for (int i=0;i<3;i++) {
        n_rv[i]=rv1[i]/rv_m;
    }

    q[0]=cos(rv_m/2.0);
    q[1]=sin(rv_m/2.0)*n_rv[0];
    q[2]=sin(rv_m/2.0)*n_rv[1];
    q[3]=sin(rv_m/2.0)*n_rv[2];
}

/* transform quaternion to rotation vector --------------------------
*
*args   : double *q       I   quaternion (4x1)
*         double *rv      O   rotation vector (3x1)
*return : none
*-------------------------------------------------------------------------------*/
extern void quat2rv(const double *q, double *rv)
{
    int i;
    double rv_2=0.0,srv_2=0.0;

    rv_2=acos(q[0]);
    srv_2=2.0*rv_2/sin(rv_2);

    for (i=0;i<3;i++) {
        rv[i]=srv_2*q[i+1];
    }

}
/* nx1 vector multiply the number --------------------------
*
*args   : int     n        I   size of the vector (3x3)
*         double *v1       I   input  vector1 (nx1)
*         double  f1       I   coefficient of vector1
          double *vo       O   output vector  (nx1)
*return : none
*-------------------------------------------------------------------------------*/
extern void vnmul(const int n, const double *v1, double f1, double *vo)
{  
    int i;
     
    for (i=0;i<n;i++)
    {
        vo[i]=f1*v1[i];
    }
}

/* nx1 vector add the nx1 vector --------------------------
*
*args   : int     n        I   size of the vector
*         double *v1       I   input  vector1 (nx1)
*         double  f1       I   coefficient of vector1
*         double *v2       I   input  vector2 (nx1)
*         double  f2       I   coefficient of vector2
          double *vo       O   output vector  (nx1)
*return : none
*-------------------------------------------------------------------------------*/
extern void vnadd(const int n, const double *v1, double f1, const double *v2, double f2, double *vo)
{  
    int i;
     
    for (i=0;i<n;i++)
    {
        vo[i]=f1*v1[i]+f2*v2[i];
    }
}

/* 3x3 matrix add 3x3 matrix --------------------------
*
*args   : double *mat1      I   input  matrix1 (3x3)
*         double  f1        I   coefficient of matrix1
*         double *mat2      I   input  matrix2 (3x3)
*         double  f2        I   coefficient of matrix2
          double *mat3      O   output matrix  (3x3)
*return : none
*-------------------------------------------------------------------------------*/
extern void Mat3add2(const double *mat1, double f1, const double *mat2, double f2, double *mat3)
{
    int i;

    for (i=0;i<9;i++)
    {
        mat3[i]=f1*mat1[i]+f2*mat2[i];
    }
}

/* 3x3 matrix multiply the vector --------------------------
*
*args   : double  f        I   matrix coefficient
*         double *mat      I   input  matrix (3x3)
*         double *vi       I   input  vector (3x1)
          double *vo       O   output vector (3x1)
*return : none
*-------------------------------------------------------------------------------*/
extern void Mat3mulv(double f, const double *mat, const double *vi, double *vo)
{
    int i;

    vo[0]=mat[0]*vi[0]+mat[1]*vi[1]+mat[2]*vi[2];
    vo[1]=mat[3]*vi[0]+mat[4]*vi[1]+mat[5]*vi[2];
    vo[2]=mat[6]*vi[0]+mat[7]*vi[1]+mat[8]*vi[2];

    for (i=0;i<3;i++)
    {
        vo[i]*=f;
    }
}

/* 1x3 vector multiply 3x3 matrix  --------------------------
*
*args   : double  f        I   matrix coefficient
*         double *vi       I   input  vector (3x1)
*         double *mat      I   input  matrix (3x3)
          double *vo       O   output vector (3x1)
*return : none
*-------------------------------------------------------------------------------*/
extern void vmulMat3(double f, const double *vi, const double *mat, double *vo)
{
    int i;

    vo[0]=vi[0]*mat[0]+vi[1]*mat[3]+vi[2]*mat[6];
    vo[1]=vi[0]*mat[1]+vi[1]*mat[4]+vi[2]*mat[7];
    vo[2]=vi[0]*mat[2]+vi[1]*mat[5]+vi[2]*mat[8];

    for (i=0;i<3;i++)
    {
        vo[i]*=f;
    }
}

/* Skew-symmetric matrix multiply the vector --------------------------
*
*args   : double  f       I   Skew-symmetric matrix coefficient
*         double *v1      I   vector (3x1)
*         double *v2      I   vector (3x1)
*         double *vx      O   vector (3x1)
*return : none
*-------------------------------------------------------------------------------*/
extern void vskewmv(double f, const double *v1, const double *v2, double *vx)
{
    double mat[9]={0};

    mat[0]=0;            mat[1]=-f*v1[2];   mat[2]= f*v1[1];
    mat[3]= f*v1[2];     mat[4]=0;          mat[5]=-f*v1[0];
    mat[6]=-f*v1[1];     mat[7]= f*v1[0];   mat[8]=0; 
    Mat3mulv(1.0,mat,v2,vx);
}

/* vector multiply the skew-symmetric matrix --------------------------
*
*args   : double  f       I   Skew-symmetric matrix coefficient
*         double *v1      I   vector (1x3)
*         double *v2      I   vector (3x1)
*         double *vx      O   vector (3x1)
*return : none
*-------------------------------------------------------------------------------*/
extern void vmvskew(double f, const double *v1, const double *v2, double *vx)
{
    double mat[9]={0};

    mat[0]=0;            mat[1]=-f*v2[2];   mat[2]= f*v2[1];
    mat[3]= f*v2[2];     mat[4]=0;          mat[5]=-f*v2[0];
    mat[6]=-f*v2[1];     mat[7]= f*v2[0];   mat[8]=0; 
    vmulMat3(1.0,v1,mat,vx);
}

/* Skew-symmetric matrix multiply the 3x3 matrix --------------------------
*
*args   : double  f       I   Skew-symmetric matrix coefficient
*         double *v1      I   vector (3x1)
*         double *mat1    I   vector (3x3)
*         double *mat2    O   vector (3x3)
*return : none
*-------------------------------------------------------------------------------*/
extern void vskewmat3(double f, const double *v1, const double *mat1, double *mat2)
{
    double mat[9]={0};

    vskew(f,v1,mat);
    Mat3mul2(1.0,mat,mat1,mat2);
}

/* 3x3 matrix multiply the skew-symmetric matrix --------------------------
*
*args   : double  f       I   Skew-symmetric matrix coefficient
*         double *mat1    I   vector (3x3)
*         double *v1      I   vector (3x1)
*         double *mat2    O   vector (3x3)
*return : none
*-------------------------------------------------------------------------------*/
extern void Mat3mvskew(double f, const double *mat1, const double *v1, double *mat2)
{
    double mat[9]={0};

    vskew(f,v1,mat);
    Mat3mul2(1.0,mat1,mat,mat2);
}

/* Multiplying three matrixs --------------------------
*
*args   : double *mat1     I   input  matrix1 (3x3)
*         double *mat2     I   input  matrix2 (3x3)
          double *mat3     I   input  matrix3 (3x3)
          double *mat      O   output matrix  (3x3)
*return : none
*-------------------------------------------------------------------------------*/
extern void Mat3mul3(const double *mat1, const double *mat2, const double *mat3, double *mat)
{
    double temp[9]={0};
    Mat3mul2(1.0,mat1,mat2,temp);
    Mat3mul2(1.0,temp,mat3,mat);
}

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
            trace(7,"warning: IMU measurement output is zero: %s, week=%.0f, sec=%.4f\n!",file,week,sec);
        }
        stat=addimudata(imu,&imud);
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

/* free imu data -----------------------------------------------------*/
extern void freeimu(imu_t *imu)
{
    trace(3,"freeimu:\n");

    free(imu->data); imu->data=NULL; imu->n =imu->nmax =0;
}

/* initialize earth related parameters -----------------------------------------------*/
extern void earth_init(const double *pos, const double *vel, eth_t *eth)
{   
    int i;
    double sinB2,sin2B2,secB,Rmh,Rnh,temp1[3];
    sinB2=sin(pos[0])*sin(pos[0]);
    sin2B2=sin(2.0*pos[0])*sin(2.0*pos[0]);
    secB=1.0/cos(pos[0]);

    eth->g0=9.780327;
    eth->Rl=RE_WGS84;
    eth->alpha=FE_WGS84;
    eth->Rs=eth->Rl*(1.0-eth->alpha);
    eth->e1=sqrt(eth->Rl*eth->Rl-eth->Rs*eth->Rs)/eth->Rl;
    eth->e2=sqrt(eth->Rl*eth->Rl-eth->Rs*eth->Rs)/eth->Rs;
    eth->wie=OMGE;
    /* gravity related parameters, ref psins */
    eth->beta[0]=5.2790414E-3;
    eth->beta[1]=2.32718E-5;
    eth->beta[2]=3.086E-6;
    eth->beta[3]=8.08E-9;

    eth->RN=eth->Rl/sqrt(1.0-eth->e1*eth->e1*sinB2);
    eth->RM=eth->RN*(1.0-eth->e1*eth->e1)/(1.0-eth->e1*eth->e1*sinB2);
    Rmh=eth->RM+pos[2];Rnh=eth->RN+pos[2];

    /* the default is forward mechanization */
    eth->wnie[0]=0.0; eth->wnie[1]=eth->wie*cos(pos[0]); eth->wnie[2]=eth->wie*sin(pos[0]);
    eth->wnen[0]=-vel[1]/Rmh; 
    eth->wnen[1]= vel[0]/Rnh; 
    eth->wnen[2]= vel[0]*sin(pos[0])/(Rnh*cos(pos[0])); 
    vnadd(3,eth->wnie,1.0,eth->wnen,1.0,eth->wnin);
    for (i=0;i<9;i++)
    {
        eth->F1 [i]=0.0; 
        eth->F2 [i]=0.0;
        eth->F3 [i]=0.0;
        eth->F4 [i]=0.0;
        eth->Fav[i]=0.0;
        eth->Frr[i]=0.0;        
        eth->Frp[i]=0.0;
    }

    eth->g=eth->g0*(1+eth->beta[0]*sinB2-eth->beta[1]*sin2B2)-eth->beta[2]*pos[2];
    eth->gn[0]=0.0; eth->gn[1]=0.0; eth->gn[2]=-1.0*eth->g;

    vnadd(3,eth->wnie,2.0,eth->wnen,1.0,eth->wnien);
    vskewmv(1.0,eth->wnien,vel,temp1);
    /* harmful acceleration */
    vnadd(3,eth->gn,1.0,temp1,-1.0,eth->gcc);

}

/* update earth related parameters -----------------------------------------------*/
extern void earth_update(const prcopt_t *popt, const double *pos, const double *vel, eth_t *eth)
{ 
    int i;
    double sinB,cosB,cos2B,tanB,sin2B,sinB2,sin2B2,secB,secB2,Rmh,Rnh,Rmh2,Rnh2,temp1[3];
    sinB=sin(pos[0]);sin2B=sin(2.0*pos[0]);sinB2=sinB*sinB;sin2B2=sin2B*sin2B;
    cosB=cos(pos[0]);cos2B=cos(2.0*pos[0]);tanB=tan(pos[0]);
    secB=1.0/cosB;secB2=secB*secB;

    eth->RN=eth->Rl/sqrt(1.0-eth->e1*eth->e1*sinB2);
    eth->RM=eth->RN*(1.0-eth->e1*eth->e1)/(1.0-eth->e1*eth->e1*sinB2);
    Rmh=eth->RM+pos[2];Rnh=eth->RN+pos[2];
    Rmh2=Rmh*Rmh;Rnh2=Rnh*Rnh;

    if (MECH_FORWARD==popt->reverse) {
        eth->wnie[0]=0.0; eth->wnie[1]=eth->wie*cos(pos[0]); eth->wnie[2]=eth->wie*sin(pos[0]);       
    }
    /* backward mechanization negates the wie symbol */
    else if (MECH_BACKWARD==popt->reverse) {
        eth->wnie[0]=0.0; eth->wnie[1]=-eth->wie*cos(pos[0]); eth->wnie[2]=-eth->wie*sin(pos[0]);  
    }

    /* wnen */
    eth->wnen[0]=-vel[1]/Rmh; 
    eth->wnen[1]= vel[0]/Rnh; 
    eth->wnen[2]= vel[0]*sin(pos[0])/(Rnh*cos(pos[0])); 
    /* wnin */
    vnadd(3,eth->wnie,1.0,eth->wnen,1.0,eth->wnin);

    for (i=0;i<9;i++)
    {
        eth->F1 [i]=0.0; 
        eth->F2 [i]=0.0;
        eth->F3 [i]=0.0;
        eth->F4 [i]=0.0;
        eth->Fav[i]=0.0;
        eth->Frr[i]=0.0;        
        eth->Frp[i]=0.0;
    }
    eth->Frp[1]=1.0/Rmh;               
    eth->Frp[3]=secB/Rnh;               
    eth->Frp[8]=1.0;

    eth->g=eth->g0*(1+eth->beta[0]*sinB2-eth->beta[1]*sin2B2)-eth->beta[2]*pos[2];
    eth->gn[0]=0.0; eth->gn[1]=0.0; eth->gn[2]=-eth->g;

    vnadd(3,eth->wnie,2.0,eth->wnen,1.0,eth->wnien);
    vskewmv(1.0,eth->wnien,vel,temp1);
    /* harmful acceleration */
    vnadd(3,eth->gn,1.0,temp1,-1.0,eth->gcc);       

    /* state transition matrix related parameters */
    eth->F1 [3]=-eth->wie*sinB;      eth->F1 [6]=eth->wie*cosB;
    eth->F2 [2]=vel[1]/Rmh2;         eth->F2 [5]=-vel[0]/Rnh2;           eth->F2 [6]=vel[0]*secB2/Rnh; eth->F2[8]=-vel[0]*tanB/Rnh2;
    eth->F3 [3]=-2.0*eth->beta[3]*pos[2]*cos2B;                          
    eth->F3 [5]=-eth->beta[3]*sin2B; 
    eth->F3 [6]=-eth->g0*sin2B*(eth->beta[0]-4.0*eth->beta[1]*cos2B);    eth->F3 [8]=eth->beta[2];
    eth->F4 [0]=-eth->g0/Rnh;        eth->F4 [4]=-eth->g0/Rmh;           eth->F4 [8]=-2.0*eth->g0/(sqrt(Rmh*Rnh)+pos[2]);
    eth->Fav[1]=-1.0/Rmh;            eth->Fav[3]=1.0/Rnh;                eth->Fav[6]=tanB/Rnh;    
    eth->Frr[0]=vel[2]/Rnh-vel[1]*tanB/Rmh;                              eth->Frr[1]=vel[0]*tanB/Rmh;
    eth->Frr[2]=-vel[0]/Rnh;         eth->Frr[4]=vel[2]/Rmh;             eth->Frr[5]=-vel[1]/Rmh;
}

/* initialize ins related parameters -----------------------------------------------*/
extern int ins_init(ins_t *ins, const prcopt_t *popt)
{
    double install_angle[3]={0.0};
    int i,nx;

    ins->nx=15;nx=ins->nx;
    ins->F=zeros(nx,nx); ins->Phi=zeros(nx,nx);
    ins->G=zeros(nx,nx); ins->Q=zeros(nx,nx);

    ins->time.sec=0.0; ins->time.time=0.0;
    ins->interval=1.0/popt->insample;
    ins->nn=popt->nn;
    ins->dttol=ins->interval/1e3;
    ins->discretime=(popt->insample%10)==0?1e-1:((popt->insample%25)==0?25*ins->interval:ins->interval);
    
    /* accelerometer and gyroscope velocity random walk and angle random walk and bias drive noise */
    ins->corr_time=popt->corr_time; 
    ins->psd_gyro=popt->psd_gyro;
    ins->psd_acce=popt->psd_acce;
    ins->psd_bg=popt->psd_bg;
    ins->psd_ba=popt->psd_ba;

    /* initialize zupt configuration options (sliding window length and detection threshold) */
    ins->zupt.window=popt->insample;
    ins->zupt.gthres=popt->zupt_gthres;

    for (i=0;i<15;i++)
    {
        ins->xa[i]=0.0;
        
        if (i<3)             ins->Q[i+i*nx]=ins->psd_gyro*ins->discretime;
        else if(i>=3&&i<6)   ins->Q[i+i*nx]=ins->psd_acce*ins->discretime;
        else if(i>=9&&i<12)  ins->Q[i+i*nx]=ins->psd_bg*ins->discretime;
        else if(i>=12&&i<15) ins->Q[i+i*nx]=ins->psd_ba*ins->discretime;
    }
    /* trace(12,"Q=\n"); tracemat(12,ins->Q,nx,nx,20,16,0); */ /*ok*/

    for (i=0;i<3;i++)
    {
        /* initialize the lever and motion constraint information */
        install_angle[i]=popt->install_angle[i];
        ins->lever_nhc[i]=popt->lever_nhc[i];
        ins->lever[i]=popt->lever[i];

        /* INS basic parameters */
        ins->dw[i]=0.0;
        ins->dv[i]=0.0;
        ins->p1dv[i]=0.0;
        ins->p1dw[i]=0.0;
        ins->p1vel[i]=0.0;
        ins->p2vel[i]=0.0;
        ins->bg[i]=0.0;
        ins->ba[i]=0.0;
        ins->wbib[i]=0.0;
        ins->fb[i]=0.0;
        ins->fn[i]=0.0;

        ins->p1pos[i]=0.0;
        ins->pos[i]=0.0;
        ins->vel[i]=0.0;
        ins->body_vel[i]=0.0;
        ins->nhc_vel[i]=0.0;
        ins->att[i]=0.0;       
    }

    /* initialize the posture matrix and installation angle matrix */
    att2qnb(ins->att,ins->qnb);
    att2Cnb(ins->att,ins->Cnb);
    att2Cnb(install_angle,ins->Cvb);
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

/* ins initial alignment -------------------------------------------*/
extern int ins_align(rtk_t *rtk, obsd_t *obs, obsd_t *obs_old, int n, int n_old, nav_t *nav, imud_t *imu, const prcopt_t *opt)
{
    /* NOTE: The structure copy is a shallow copy! */
    prcopt_t popt=*opt;
    rtk_t  rtk_={0}; 
    ins_t *ins=&rtk->ins;
    int i,j,vel_flag,nr=0,nr_old=0;
    double att[3]={0.0},pos[3]={0.0},vn[3]={0.0};

    popt.GI_mode=GINS_OFF;  /* set to GINS_OFF mode */
    rtkinit(&rtk_,&popt,NULL); /* note: initialize rtk_ instead of assigning rtk to rtk_ to avoid shallow copying of the structure. */

    /* initialize rtk_ TDCP velocity */
    for (i=0;i<3;i++) rtk->sol.tdcp_vel[i]=rtk_.sol.tdcp_vel[i]=0.0; 

    /* for GNSS/INS LC with pos file ,set alignment type to manual */
    if (PMODE_LC_POS==popt.mode) popt.alingetype=INSALI_MANUAL;

    /* determine the number of satellites of rover in the current epoch and the previous epoch */
    for (i=0;i<n;i++) if (obs[i].rcv==1) nr++;
    for (i=0;i<n_old;i++) if (obs_old[i].rcv==1) nr_old++;

    /* manual alignment */
    if (!rtk->align&&INSALI_MANUAL==popt.alingetype&&popt.ts.time)
    {   
        if ((fabs(timediff(ins->time,popt.ts))-rtk->ins.dttol)<=rtk->ins.nn*rtk->ins.interval/2.0) {
            /* initialize ins position, velocity and attitude */
            init_inspva(ins,popt.initpos,popt.initvel,popt.initatt); 

            trace(12,"INS initial alignment completed: %s!\n",Debug_Glo.chTime); 
            showerr("INS initial alignment completed: %s!",Debug_Glo.chTime); 
            return 1;            
        }
        else if (timediff(ins->time,popt.ts)>0) {
            showmsg("warning : start time is smaller than GNSS/INS matching time!\n"); 
            return 0;
        }
    }

    /* velocity vector assisted yaw initialization based on tdcp */
    if (!rtk->align&&INSALI_VELTOR==popt.alingetype&&obs_old[0].time.time&&SYNC_YES==rtk->upte) {

        if (vel_flag=tdcp_vel(rtk,obs,obs_old,nr,nr_old,nav,&popt)) {
            /* initialize INS position using GNSS solution */
            if (!rtkpos(&rtk_,obs,n,nav)) {
                trace(7,"rtkpos error: GNSS unavailable during INS align!\n");
                return 0;
            }
            /* initialize position (from GNSS) and velocity (from tdcp) */
            ecef2pos(rtk_.sol.rr,pos);        
            ecef2enu(pos,rtk->sol.rr+3,vn);
            
            /* initialize pitch and yaw using the velocity in the n frame */
            att[0]=atan2(vn[2],sqrt(vn[0]*vn[0]+vn[1]*vn[1])); /* pitch angle */
            att[2]=-atan2(vn[0],vn[1]);                        /* yaw angle */
            /* yaw (counterclockwise to clockwise, deg) */
            rtk->sol.yaw=att[2]*R2D;
            if (rtk->sol.yaw<0) rtk->sol.yaw*=-1;
            else rtk->sol.yaw=360.0-rtk->sol.yaw;
            /* output the initial alignment solution status */
            /* outsolstat(rtk,nav); */

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
            return 1;            
        }
    }

    /* reinitialize INS in the event of a long-term GNSS outage */
    if (rtk->align&&rtk->outage>MAX_OUTIME&&SYNC_YES==rtk->upte) {

        if (vel_flag=tdcp_vel(rtk,obs,obs_old,nr,nr_old,nav,&popt)) {
            /* initialize INS position using GNSS solution */
            if (!rtkpos(&rtk_,obs,n,nav)) {
                trace(7,"rtkpos error: GNSS unavailable during INS reinitialization!\n");
                return 1;
            }
            /* if GNSS becomes available after a long interruption, set the GNSS interruption count to zero */
            rtk->outage=0;

            /* reinitialize ins position and velocity, consider lever arm correction */
            ecef2pos(rtk_.sol.rr,pos);        
            ecef2enu(pos,rtk->sol.rr+3,vn);
            gnss2ins(rtk,pos,ins->pos,1);
            gnss2ins(rtk,vn,ins->vel,2);
            init_inspva(ins,ins->pos,ins->vel,NULL);

            /* for (i=3;i<6;i++) {
                for (j=0;j<rtk->lcgins.nx;j++) rtk->lcgins.P[i+j*rtk->lcgins.nx]=0.0;
                for (j=0;j<rtk->lcgins.nx;j++) rtk->lcgins.P[j+i*rtk->lcgins.nx]=0.0; 
                rtk->lcgins.P[i+i*rtk->lcgins.nx]=1.0;

                for (j=0;j<rtk->lcgins.nx;j++) rtk->lcgins.P[(i+3)+j*rtk->lcgins.nx]=0.0;
                for (j=0;j<rtk->lcgins.nx;j++) rtk->lcgins.P[j+(i+3)*rtk->lcgins.nx]=0.0;
                rtk->lcgins.P[(i+3)+(i+3)*rtk->lcgins.nx]=10.0;
            } */
           /* for (i=0;i<3;i++) {
                for (j=0;j<rtk->nx;j++) rtk->P[(i+3)+j*rtk->nx]=0.0;
                for (j=0;j<rtk->nx;j++) rtk->P[j+(i+3)*rtk->nx]=0.0; 
                rtk->P[(i+3)+(i+3)*rtk->nx]=1.0;

                for (j=0;j<rtk->nx;j++) rtk->P[(i+6)+j*rtk->nx]=0.0;
                for (j=0;j<rtk->nx;j++) rtk->P[j+(i+6)*rtk->nx]=0.0;
                rtk->P[(i+6)+(i+6)*rtk->nx]=100.0;
            } */
            /* trace(12,"P=\n"); tracemat(12,rtk->P,rtk->nx,rtk->nx,16,8,0); */
            /* trace(12,"P=\n"); tracemat(12,rtk->lcgins.P,rtk->lcgins.nx,rtk->lcgins.nx,16,8,0); */

            trace(12,"INS reinitialization completed: %s!\n",Debug_Glo.chTime);
            showerr("INS reinitialization completed: %s!",Debug_Glo.chTime);             
        } 
        return 1;
    }

    return (rtk->outage>MAX_OUTIME)?1:0;
}

/* TDCP-assisted motion alignment */
extern int tdcp_vel(rtk_t *rtk, const obsd_t *obs, const obsd_t *obs_old, int n, int n_old, const nav_t *nav, const prcopt_t *opt)
{
    prcopt_t opt_=*opt;
    sol_t sol={0},sol_old={0};
    double *rs,*rs_old,*dts,*dts_old,*vare,*vare_old,*resp,*resp_old,*azel,*azel_old;
    double rr[3],rr_old[3],r,dr[3],dr_old[3],er=0.0,er_old=0.0,e[3],e_old[3],freq,thres_ouj=2.0,thres=3.0;
    double sgn=(SOLTYPE_BACKWARD==opt->reverse?-1.0:1.0);
    double *v,*H,*var,*P,dx[4]={0},Q[4*4];
    int sat[MAXSAT],ir_old[MAXSAT],ir[MAXSAT];
    int stat=0,stat_old=0,i,j,k,m,nf=rtk->opt.nf,sys,fr,nx=4,nv=0,vnv[MAXFREQ]={0},max_vnv=0,info,flag=1,vel_flag=1,mode=Robust_RES;
    int vsat[MAXOBS]={0},vsat_old[MAXOBS]={0},svh[MAXOBS]={0},svh_old[MAXOBS]={0};

    if (!n||!n_old) return 0;

    /* initializing memory*/
    rs=mat(n,6);    rs_old=mat(n_old,6);   dts=mat(n,2);   dts_old=mat(n_old,2);
    vare=mat(n,1);  vare_old=mat(n_old,1); resp=mat(n,1);  resp_old=mat(n_old,1);
    azel=zeros(n,2);azel_old=zeros(n_old,2);
    v=mat(nf*n,1);  H=mat(nf*n,nx); var=mat(nf*n,1); P=zeros(nf*n,nf*n); 

    /* reset receiver velocity */
    for (i=0;i<3;i++) rtk->sol.rr[i+3]=0.0;

    /* configured in spp mode */
    if (opt_.mode!=PMODE_SINGLE||opt_.GI_mode!=GINS_OFF) {
        opt_.GI_mode=GINS_OFF;     opt_.mode=PMODE_SINGLE;
        opt_.sateph =EPHOPT_BRDC;
        opt_.ionoopt=IONOOPT_BRDC; opt_.tropopt=TROPOPT_SAAS;
    }

    /* satellite positons, velocities and clocks of current and previous epoch */
    satposs(obs[0].time,obs,n,nav,opt_.sateph,rs,dts,vare,svh);
    satposs(obs_old[0].time,obs_old,n_old,nav,opt_.sateph,rs_old,dts_old,vare_old,svh_old);

    /* satellite positons, velocities and clocks of current and previous epoch */
    stat=estpos(rtk,obs,n,rs,dts,vare,svh,nav,&opt_,NULL,&sol,azel,vsat,resp);
    stat_old=estpos(rtk,obs_old,n_old,rs_old,dts_old,vare_old,svh_old,nav,&opt_,NULL,&sol_old,azel_old,vsat_old,resp_old);

    /* check solution status */
    if (!stat||!stat_old) {
        flag=0; trace(7,"tdcp_vel: estpos error stat=%d, stat_old=%d\n",stat,stat_old); 
    }
    else {
        /* receiver position in the previous epoch and the current epoch in spp mode */
        for (i=0;i<3;i++) {
            rr[i]=sol.rr[i];
            rr_old[i]=sol_old.rr[i];
        }
    }

    /* if GNSS outage, tdcp fails */
    if (rtk->interval&&timediff(obs[0].time,obs_old[0].time)>rtk->interval) {
        flag=0; trace(7,"tdcp_vel: time difference between current and previous epoch is too large tt=%.2f\n",timediff(obs[0].time,obs_old[0].time));
    }

    /* epoch-to-epoch average velocity estimation based on tdcp */
    if (flag) {
        /* check whether a cycle slip occurs in the current epoch observation */
        init_ssatpar(rtk,NULL,n,RTK_slip,SOLQ_NONE);

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
                if (satexclude(obs_old[i].sat,vare_old[i],svh_old[i],&opt_)||satexclude(obs[j].sat,vare[i],svh[j],&opt_)) continue;      
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
        for (i=0;i<nf;i++) {
            if (!i) max_vnv=vnv[i];
            else if (vnv[i]>max_vnv) max_vnv=vnv[i];
        }

        if (max_vnv<nx) {
            vel_flag=0;trace(7,"tdcp_vel: not enough valid satellites nv=%d\n",nv);
        }
        else {
            /* mode=(max_vnv>nx)?Robust_OFF:Robust_RES; */
            nv=outrej_spp(nv,4,thres_ouj,v,H,var,NULL,NULL,NULL,vsat,0);

            /* calculate the weight matrix */
            diag_Cov(nv,var,P,diag_wei);

            /* trace(12,"TDCP H=\n");tracemat(12,H,nv,nx,9,4,0);
            trace(12,"TDCP v=\n");tracemat(12,v,nv,1,9,4,0);
            trace(12,"TDCP P=\n");tracemat(12,P,nv,nv,9,4,0); */

            /* least square estimation */
            if ((info=lsq_roubst(H,v,P,4,nv,dx,Q,mode))) {
                vel_flag=0;trace(7,"tdcp lsq error info=%d\n!",info);
            }

            /* calculate the posterior residuals */
            matmul("NN",nv,nx,1,H,dx,v,1.0,-1.0);

            /* validate solution */
            if (valsol(&sol,azel,vsat,n,opt,v,P,nv,nx)==0) {
                vel_flag=0;trace(7,"validation of tdcp solution failed nv=%d\n",nv);
            }
            
            if (vel_flag) {
                for (i=0;i<3;i++) rtk->sol.rr[i+3]=sgn*dx[i]/rtk->interval;  
                matcpy(rtk->sol.tdcp_vel,rtk->sol.rr+3,3,1);
                /* outsolstat(rtk,nav); */                   
            }           
        }      
    }

    /* if tdcp fails, velocity estimation is performed using Doppler observations */
    if (flag&&!vel_flag&&rtk->dopsgn&&estvel(rtk,obs,n,rs,dts,nav,&opt_,&sol,azel,vsat)) {
        vel_flag=1;
        matcpy(rtk->sol.rr+3,sol.rr+3,3,1);
    }

    /* if tdcp and dopple fail, velocity estimation is performed using the spp position difference between epochs */
    if (flag&&!vel_flag) {
        vel_flag=1;
        /* calculate the average velocity between epochs */
        for (i=0;i<3;i++) rtk->sol.rr[i+3]=sgn*(rr[i]-rr_old[i])/rtk->interval;
    }

    /* when the vehicle velocity exceeds the threshold, the alignment is considered complete */
    if (!rtk->align&&vel_flag&&norm(rtk->sol.rr+3,3)<thres) {
        vel_flag=0;
    }

    /* freeing up memory */
    free(rs);   free(rs_old);   free(dts);  free(dts_old);
    free(vare); free(vare_old); free(resp); free(resp_old);
    free(azel); free(azel_old);
    free(v);    free(H);        free(var);  free(P);        

    return flag?vel_flag:0;
}

/* zero speed detection*/
extern void zerovel_detect(rtk_t *rtk, imud_t *imu)
{
    zupt_t *zupt=&rtk->ins.zupt;
    int i,j,nn=rtk->ins.nn;
    double N=zupt->window,dw[3]={0.0},ndw;

    /* for multi-sample mode, the average value of the gyroscope is used */
    for (j=0;j<3;j++) for (i=0;i<nn;i++) dw[j]+=imu[i].dw[j];
    for (j=0;j<3;j++) dw[j]/=nn;
    ndw=norm(dw,3);
    
    /* initialization */
    if (!zupt->iimu) {
        zupt->Gm=ndw;
        zupt->Gd=0.0;
    }
    if (zupt->iimu) {
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

/* motion constraints */
extern void motion_constraints(rtk_t *rtk, const prcopt_t *popt) 
{
    ins_t *ins=&rtk->ins;
    sol_t *sol=(GINS_TC==popt->GI_mode)?&rtk->sol:&rtk->lcgins.sol;
    int i,j,nx=ins->nx,nv=4,info;
    double zupt_time,vel,*xp,*Pp,*H,*v,*var,*R;

    /* detected vehicle stationary time (s) and GNSS velocity */
    zupt_time=ins->zupt.count*ins->interval*ins->nn;
    vel=norm(rtk->sol.rr+3,3);

    /* initializing memory */
    xp=zeros(nx,1); Pp=zeros(nx,nx); R=zeros(nv,nv);
    H=mat(nv,nx); v=mat(nv,1); var=mat(nv,1);

    /* copy the covariance matrix */
    if (GINS_LC==popt->GI_mode||GINS_STC==popt->GI_mode) matcpy(Pp,rtk->lcgins.P,nx,nx);
    else if (GINS_TC==popt->GI_mode) pmatcpy(Pp,nx,nx,0,0,nx,nx,rtk->P,rtk->nx,rtk->nx,0,0,nx,nx); 

    /* NOTE: the vehicle is considered stationary only when the zero speed detection is passed, 
    the stationary state is greater than 1s and the calculated vehicle speed is less than 0.1m/s */
    if (popt->constraint[1]&&zupt_time>1.0&&(vel>0&&vel<0.1)) { /* zupt */
        nv=motion_update(rtk,H,v,var,0,nx,CONS_ZUPT);
        sol->iFlag=SOLF_ZUPT; /* zupt flag */
    }
    else if (popt->constraint[0]) { /* nhc */
        nv=motion_update(rtk,H,v,var,0,nx,CONS_NHC);        
    }
    if (popt->constraint[2]&&zupt_time>1.0&&(vel>0&&vel<0.1)) { /* zihr */
        nv=motion_update(rtk,H,v,var,nv,nx,CONS_ZIHR);
    }
    
    /* measurement noise covariance matrix */
    diag_Cov(nv,var,R,diag_var);

    /* measurement update */
    if ((info=filter_gins(rtk,xp,Pp,H,v,R,nx,nv,KF_GINS,Robust_OFF))) {
        trace(7,"motion_constraints: filter_gins error info=%d\n",info);
        sol->stat=SOLQ_INS;
        /* update solution status */
        update_instat(popt,ins,Pp,sol,nx);
        
        free(xp); free(Pp); free(H); free(v); free(var);
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
    double Cbn[9]={0.0},Cvn[9]={0.0},lever_v[9]={0.0},vel_v[9]={0.0},Ha_bg[3]={0.0},att[3]={0.0},yaw=0.0;

    DCMT(ins->Cnb,Cbn);
    Mat3mul2(1.0,ins->Cvb,Cbn,Cvn);
    /* vehichle velocity of v frame */
    Mat3mulv(1.0,Cvn,ins->vel,ins->nhc_vel);
    Mat3mvskew(-1.0,Cvn,ins->vel,vel_v);
    Mat3mvskew(-1.0,ins->Cvb,ins->lever_nhc,lever_v);

    /* H of ZIHR */
    Ha_bg[0]=-sin(ins->att[1])/cos(ins->att[0])*rtk->interval; 
    Ha_bg[1]=0.0; 
    Ha_bg[2]=cos(ins->att[1])/cos(ins->att[0])*rtk->interval;

    /* converts the yaw from clockwise to counterclockwise */
    if (GINS_LC==rtk->opt.GI_mode) for (i=0;i<3;i++) att[i]=rtk->lcgins.sol.att[i];
    else for (i=0;i<3;i++) att[i]=rtk->sol.att[i];

    if (att[2]<=180) yaw=-att[2]*D2R;
    else yaw=(360.0-att[2])*D2R;

    /* determine constraint model */
    if (CONS_NHC==mode) {
        inv=2;
        trace(12,"nhc_constraints: v=\n");tracemat(12,ins->nhc_vel,3,1,9,4,0);
    }
    else if (CONS_ZUPT==mode) {
        inv=3;
        trace(12,"zupt_constraints: v=\n");tracemat(12,ins->nhc_vel,3,1,9,4,0);
    }
    else if (CONS_ZIHR==mode) {
        inv=1;
        trace(12,"zihr_constraints: yaw=\n");tracemat(12,ins->att,3,1,9,4,0);
    }

    for (i=0;i<inv;i++) {
        if (CONS_NHC==mode) k=k2[i];
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

/* for backward processing mode, the sign of the INS velocity and gyroscope bias is inverted */
extern void pos_reverse(const prcopt_t *popt, ins_t *ins, int *reverse_flag)
{
    int i;

    if (SOLTYPE_COMBINED<=popt->soltype&&SOLTYPE_BACKWARD==popt->reverse) {
        for (i=0;i<3;i++) {
            ins->vel[i]=-ins->vel[i];
            ins->bg[i]=-ins->bg[i];
        }
        *reverse_flag=1;
    }
}

/* ins mechanization -----------------------------------------------*/
extern void ins_mech(ins_t *ins, imud_t *imu, const prcopt_t *popt) 
{
    int i;
    double vel_m[3],pos_m[3],temp1[3],temp2[3],dv_rot[3],dv_pad[3],da_con[3],interval,*I3;
    double rv[3],rs[9],Irs[9],Ce[9],dw[3],dv[3],delta_vb[3],delta_vn[3];
    double Cnnk[9],Cnkn[9],Cbbk[9],Cnb[9],Cbn[9]={0.0};
    double qnb[4],qnkn[4],qbbk[4];

    I3=eye(3);

    /* sample interval of imu */
    interval=ins->interval*ins->nn;

    /* bias correction for gyroscopes and accelerometers */
    imu_fedback(ins,imu);

    /* initialize the attitude matrix and quaternion */
    for (i=0;i<9;i++) Cnb[i]=ins->Cnb[i];
    for (i=0;i<4;i++) qnb[i]=ins->qnb[i];

    /* extrapolate velocity and position at k-1/2 */
    vnadd(3,ins->p1vel,3.0/2.0,ins->p2vel,-1.0/2.0,vel_m);
    vnmul(3,vel_m,interval/2.0,temp1);

    Mat3mulv(1.0,ins->eth.Frp,temp1,temp2);
    vnadd(3,ins->p1pos,1.0,temp2,1.0,pos_m);

    /* update earth parameters at k-1/2 */
    earth_update(popt,pos_m,vel_m,&ins->eth);

    /* coning error, rotation error and paddling error compensation */
    conpad_fedback(ins,dw,dv,da_con,dv_rot,dv_pad);

    vnmul(3,dw,1.0/interval,ins->wbib);
    vnmul(3,dv,1.0/interval,ins->fb);

    if (ATT_DCM==popt->att_type) Mat3mulv(1.0,Cnb,ins->fb,ins->fn);        
    else quatmulv(1.0,qnb,ins->fb,ins->fn);

    /* forward mechanization, sequentially performing velocity update, position update, and attitude update */
    if (MECH_FORWARD==popt->reverse) {

        /* [I-1/2*(zetax)*Cnb(k-1)] of velocity update */
        vnmul(3,ins->eth.wnin,interval,rv);
        vskew(1.0,rv,rs);  
        Mat3add2(I3,1.0,rs,-1.0/2.0,Irs);
        if (ATT_QUAT==popt->att_type) qnb2Cnb(qnb,Cnb);  
        Mat3mul2(1.0,Irs,Cnb,Ce);

        /* rotation and paddling error compensation items of velocity update */
        for (i=0;i<3;i++)
        {
            delta_vb[i]=dv[i]+dv_rot[i]+dv_pad[i];
        }
        Mat3mulv(1.0,Ce,delta_vb,delta_vn);

        /* NOTE: velocity update */
        for (i=0;i<3;i++) {
            ins->vel[i]=ins->p1vel[i]+(ins->eth.gcc[i]*interval)+delta_vn[i];
        }

        /* NOTE: position update */
        Mat3mulv(1.0/2.0*interval,ins->eth.Frp,ins->p1vel,temp1);
        Mat3mulv(1.0/2.0*interval,ins->eth.Frp,ins->vel,temp2);
        for (i=0;i<3;i++) {
            ins->pos[i]=ins->p1pos[i]+(temp1[i]+temp2[i]);
        }

        /* NOTE: attitude update */
        if (ATT_DCM==popt->att_type) {
            rv2DCM(interval,ins->eth.wnin,Cnnk);
            DCMT(Cnnk,Cnkn);
            rv2DCM(1.0,da_con,Cbbk);
            Mat3mul3(Cnkn,Cnb,Cbbk,ins->Cnb);
            Cnb2att(ins->Cnb,ins->att);  
            att2qnb(ins->att,ins->qnb);      
        }
        else {
            rv2quat(-1.0*interval,ins->eth.wnin,qnkn);
            rv2quat(1.0,da_con,qbbk);
            quatmul3(qnkn,qnb,qbbk,ins->qnb);
            qnbnorm(ins->qnb);
            qnb2Cnb(ins->qnb,ins->Cnb);
            Cnb2att(ins->Cnb,ins->att);
        }        
    }
    /* backward mechanization, sequentially performing attitude update, velocity update, and position update */
    else if (MECH_BACKWARD==popt->reverse) {

        /* NOTE: attitude update */
        if (ATT_DCM==popt->att_type) {
            rv2DCM(interval,ins->eth.wnin,Cnnk);
            DCMT(Cnnk,Cnkn);
            rv2DCM(1.0,da_con,Cbbk);
            Mat3mul3(Cnkn,Cnb,Cbbk,ins->Cnb);
            Cnb2att(ins->Cnb,ins->att);  
            att2qnb(ins->att,ins->qnb);      
        }
        else {
            rv2quat(-1.0*interval,ins->eth.wnin,qnkn);
            rv2quat(1.0,da_con,qbbk);
            quatmul3(qnkn,qnb,qbbk,ins->qnb);
            qnbnorm(ins->qnb);
            qnb2Cnb(ins->qnb,ins->Cnb);
            Cnb2att(ins->Cnb,ins->att);
        } 

        /* [I-1/2*(zetax)*Cnb(k-1)] of velocity update */
        vnmul(3,ins->eth.wnin,interval,rv);
        vskew(1.0,rv,rs);  
        Mat3add2(I3,1.0,rs,-1.0/2.0,Irs);
        if (ATT_QUAT==popt->att_type) qnb2Cnb(qnb,Cnb);  
        Mat3mul2(1.0,Irs,Cnb,Ce);

        /* rotation and paddling error compensation items of velocity update */
        for (i=0;i<3;i++)
        {
            delta_vb[i]=dv[i]+dv_rot[i]+dv_pad[i];
        }
        Mat3mulv(1.0,Ce,delta_vb,delta_vn);

        /* NOTE: velocity update */
        for (i=0;i<3;i++) {
            ins->vel[i]=ins->p1vel[i]+(ins->eth.gcc[i]*interval)+delta_vn[i];
        }

        /* NOTE: position update */
        Mat3mulv(1.0/2.0*interval,ins->eth.Frp,ins->p1vel,temp1);
        Mat3mulv(1.0/2.0*interval,ins->eth.Frp,ins->vel,temp2);
        for (i=0;i<3;i++) {
            ins->pos[i]=ins->p1pos[i]+(temp1[i]+temp2[i]);
        }
    }

    /* update ins velocty in b frame*/
    DCMT(ins->Cnb,Cbn);
    Mat3mulv(1.0,Cbn,ins->vel,ins->body_vel);

    /* update INS previous related parameters  */
    update_ins(ins);

    free(I3);
}

/* update INS state transition matrix F and noise driving matrix G */
extern void phi_update(ins_t *ins, const prcopt_t *popt)
{
    int i,j,nx,k;
    double Faa[9],Fap[9],Far[9],Fva[9],Fvv[9],Fvp[9],Fvr[9];
    double Fvv1[9],Fvv2[9],F12[9],Fvp1[9],*Fg,*I,*I3;
    double Frr[9];

    nx=ins->nx;
    Fg=zeros(3,3);I=eye(nx);I3=eye(3);

    /* initialize F and G */
    for (i=0;i<nx;i++){
        for (j=0;j<nx;j++){
            ins->F[j+i*nx]=0.0;
            ins->G[j+i*nx]=0.0;
            if (i==j) ins->Phi[j+i*nx]=1.0; else ins->Phi[j+i*nx]=0.0;
        }
    }

    earth_update(popt,ins->pos,ins->vel,&ins->eth);

    /* NOTE: Phi angle error model */
    if (ERR_PHI==popt->err_model) {

        vskew(-1.0,ins->eth.wnin,Faa);
        Mat3add2(ins->eth.F1,1.0,ins->eth.F2,1.0,Fap);
        Mat3mul2(1.0,Fap,ins->eth.Frp,Far);

        vskew(1.0,ins->fn,Fva);
        vskewmat3(1.0,ins->vel,ins->eth.Fav,Fvv1);
        vskew(-1.0,ins->eth.wnien,Fvv2);
        Mat3add2(Fvv1,1.0,Fvv2,1.0,Fvv);
        Mat3add2(ins->eth.F1,2.0,ins->eth.F2,1.0,F12);
        vskewmat3(1.0,ins->vel,F12,Fvp1);
        Mat3add2(Fvp1,1.0,ins->eth.F3,1.0,Fvp);
        Mat3mul2(1.0,Fvp,ins->eth.Frp,Fvr);

        if (ins->corr_time>0)
        {
            Fg[0]=-1.0/ins->corr_time;Fg[4]=-1.0/ins->corr_time;Fg[8]=-1.0/ins->corr_time; 
        }
        
        for (i=0;i<15;i++)
        {
            if (i<3) {
                k=i;
                for (j=0;j<3;j++)
                {
                    ins->F[j+i*nx]=Faa[j+k*3];
                    ins->G[j+i*nx]=-ins->Cnb[j+k*3];
                }
                for (j=3;j<6;j++)
                {
                    ins->F[j+i*nx]=ins->eth.Fav[j-3+k*3];
                }
                for (j=6;j<9;j++)
                {
                    ins->F[j+i*nx]=Far[j-6+k*3];
                }
                for (j=9;j<12;j++)
                {
                    ins->F[j+i*nx]=-ins->Cnb[j-9+k*3];
                } 
            }     
            if (i>=3&&i<6) {
                k=i-3;
                for (j=0;j<3;j++)
                {
                    ins->F[j+i*nx]=Fva[j+k*3];
                }
                for (j=3;j<6;j++)
                {
                    ins->F[j+i*nx]=Fvv[j-3+k*3];
                    ins->G[j+i*nx]=ins->Cnb[j-3+k*3];
                }
                for (j=6;j<9;j++)
                {
                    ins->F[j+i*nx]=Fvr[j-6+k*3];
                }
                for (j=12;j<15;j++)
                {
                    ins->F[j+i*nx]=ins->Cnb[j-12+k*3];
                } 
            } 
            if (i>=6&&i<9) {
                k=i-6;
                for (j=3;j<6;j++)
                {
                    ins->F[j+i*nx]=I3[j-3+k*3];
                }
                for (j=6;j<9;j++)
                {
                    ins->F[j+i*nx]=ins->eth.Frr[j-6+k*3];
                }
            } 
            if (i>=9&&i<12) {
                k=i-9;
                for (j=9;j<12;j++)
                {
                    ins->F[j+i*nx]=Fg[j-9+k*3];
                    if (i==j) ins->G[j+i*nx]=1.0;
                } 
            } 
            if (i>=12&&i<15) {
                k=i-12;
                for (j=12;j<15;j++)
                {
                    ins->F[j+i*nx]=Fg[j-12+k*3];
                    if (i==j) ins->G[j+i*nx]=1.0;
                } 
            }                                                   
        }        
    }
    /* NOTE: Psi angle error model */
    else if (ERR_PSI==popt->err_model) {

        vskew(-1.0,ins->eth.wnin,Faa); /* att */

        vskew(1.0,ins->fn,Fva);
        vskew(-1.0,ins->eth.wnien,Fvv2); /* vel */
       
        vskew(-1.0,ins->eth.wnen,Frr);   /* pos */

        if (ins->corr_time>0)
        {
            Fg[0]=-1.0/ins->corr_time;Fg[4]=-1.0/ins->corr_time;Fg[8]=-1.0/ins->corr_time; 
        }
        
        for (i=0;i<15;i++)
        {
            if (i<3) {
                k=i;
                for (j=0;j<3;j++)
                {
                    ins->F[j+i*nx]=Faa[j+k*3];
                    ins->G[j+i*nx]=-ins->Cnb[j+k*3];
                }
                for (j=9;j<12;j++)
                {
                    ins->F[j+i*nx]=-ins->Cnb[j-9+k*3];
                } 
            }     
            if (i>=3&&i<6) {
                k=i-3;
                for (j=0;j<3;j++)
                {
                    ins->F[j+i*nx]=Fva[j+k*3];
                }
                for (j=3;j<6;j++)
                {
                    ins->F[j+i*nx]=Fvv2[j-3+k*3];
                    ins->G[j+i*nx]=ins->Cnb[j-3+k*3];
                }
                for (j=6;j<9;j++)
                {
                    ins->F[j+i*nx]=ins->eth.F4[j-6+k*3];
                }
                for (j=12;j<15;j++)
                {
                    ins->F[j+i*nx]=ins->Cnb[j-12+k*3];
                } 
            } 
            if (i>=6&&i<9) {
                k=i-6;
                for (j=3;j<6;j++)
                {
                    ins->F[j+i*nx]=I3[j-3+k*3];
                }
                for (j=6;j<9;j++)
                {
                    ins->F[j+i*nx]=Frr[j-6+k*3];
                }
            } 
            if (i>=9&&i<12) {
                k=i-9;
                for (j=9;j<12;j++)
                {
                    ins->F[j+i*nx]=Fg[j-9+k*3];
                    if (i==j) ins->G[j+i*nx]=1.0;
                } 
            } 
            if (i>=12&&i<15) {
                k=i-12;
                for (j=12;j<15;j++)
                {
                    ins->F[j+i*nx]=Fg[j-12+k*3];
                    if (i==j) ins->G[j+i*nx]=1.0;
                } 
            }                                                   
        }                
    }

    /* discretization of the state transition matrix Phi */
    matmul("NN",nx,nx,nx,ins->F,I,ins->Phi,ins->discretime,1.0);
    /* trace(12,"F=\n"); tracemat(12,ins->F,nx,nx,20,16,0);
    trace(12,"Phi=\n"); tracemat(12,ins->Phi,nx,nx,20,16,0);
    trace(12,"G=\n"); tracemat(12,ins->G,nx,nx,9,4,0); */ /*ok*/

    free(Fg);free(I);free(I3);
}


/* Convert INS solutions to GNSS center */
extern void ins2gnss(const prcopt_t *popt, ins_t *ins, double *pv_g, int n)
{
    int i;
    double F1[9],lever_n[3],Cbn[9],wbie[3],wbeb[3],temp[3],d_v[3],pv[6];

    Mat3mul2(1.0,ins->eth.Frp,ins->Cnb,F1);
    Mat3mulv(1.0,F1,ins->lever,lever_n);
    vnadd(3,ins->pos,1.0,lever_n,1.0,pv);       

    if (n==6) {
        DCMT(ins->Cnb,Cbn);
        Mat3mulv(1.0,Cbn,ins->eth.wnie,wbie);
        Mat3add2(ins->wbib,1.0,wbie,-1.0,wbeb);
        vskewmv(1.0,wbeb,ins->lever,temp);
        Mat3mulv(1.0,ins->Cnb,temp,d_v);
        vnadd(3,ins->vel,1.0,d_v,1.0,pv+3);        
    }

    for (i=0;i<n;i++)
    {
        pv_g[i]=pv[i];
    }
}

extern void insfix2gnss(ins_t *ins, double *pv_g, int n)
{
    int i;
    double F1[9],lever_n[3],Cbn[9],wbie[3],wbeb[3],temp[3],d_v[3],pv[6];

    Mat3mul2(1.0,ins->eth.Frp,ins->Cnb,F1);
    Mat3mulv(1.0,F1,ins->lever,lever_n);
    vnadd(3,ins->xa+6,1.0,lever_n,1.0,pv);       

    if (n==6) {
        DCMT(ins->Cnb,Cbn);
        Mat3mulv(1.0,Cbn,ins->eth.wnie,wbie);
        Mat3add2(ins->wbib,1.0,wbie,-1.0,wbeb);
        vskewmv(1.0,wbeb,ins->lever,temp);
        Mat3mulv(1.0,ins->Cnb,temp,d_v);
        vnadd(3,ins->xa+3,1.0,d_v,1.0,pv+3);        
    }

    for (i=0;i<n;i++)
    {
        pv_g[i]=pv[i];
    }
}

/* convert GNSS solutions to INS center */
extern void gnss2ins(rtk_t *rtk, double *pv_g, double *pv_i, int mode)
{
    ins_t *ins=&rtk->ins;
    int i;
    double F1[9],lever_n[3],Cbn[9],wbie[3],wbeb[3],temp[3],d_v[3];

    /* convert GNSS position to INS position */
    if (mode==1) {
        /* r_gnss=r_ins+Fpv*Cnb*lever_b */
        Mat3mul2(1.0,ins->eth.Frp,ins->Cnb,F1);
        Mat3mulv(1.0,F1,ins->lever,lever_n);

        vnadd(3,pv_g,1.0,lever_n,-1.0,pv_i);
    }
    /* convert GNSS velocity to INS velocity */
    else if (mode==2) {
        /* vn_gnss=vn_ins+Cnb([wbeb x]lever_b)*/
        DCMT(ins->Cnb,Cbn);
        Mat3mulv(1.0,Cbn,ins->eth.wnie,wbie);
        Mat3add2(ins->wbib,1.0,wbie,-1.0,wbeb);
        vskewmv(1.0,wbeb,ins->lever,temp);
        Mat3mulv(1.0,ins->Cnb,temp,d_v);

        vnadd(3,pv_g,1.0,d_v,-1.0,pv_i);      
    }
}

/* update INS previous related parameters  --------------------------
*
*args   : ins_t *ins     I    ins struct
*return : none
*-------------------------------------------------------------------------------*/
extern void update_ins(ins_t *ins)
{
    int i;

    for (i=0;i<3;i++)
    {
        ins->p1dw[i]=ins->dw[i];
        ins->p1dv[i]=ins->dv[i];
        ins->p2vel[i]=ins->p1vel[i]; 
        /* update previous epoch pos/vel by ins predicted information */       
        ins->p1vel[i]=ins->vel[i];
        ins->p1pos[i]=ins->pos[i];
    }
}

/* pure Inertial Navigation */
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

    for (i=0;i<n;i++)
    {
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
                   

