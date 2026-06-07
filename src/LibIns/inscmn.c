/*------------------------------------------------------------------------------
*ins.c : ins common functions
 *-----------------------------------------------------------------------------*/

#include "rtklib.h"

static imu_t imus={0};          /* imu data */

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
    double rv1[3],rv_m,s=0.0;

    for (i=0;i<3;i++)
    {
        rv1[i]=f*rv[i];
    }
    rv_m=norm(rv1,3);
    
    /* if the rotation vector is small, use the Taylor expansion (ref:psins rv2q.m )*/
    if (rv_m<1e-4) {
        /* cos(n/2)=1-n2/8+n4/384; sin(n/2)/n=1/2-n2/48+n4/3840 */
        q[0]=1.0-rv_m*rv_m*(1.0/8.0-rv_m*rv_m/384.0);
        s=1.0/2.0-rv_m*rv_m*(1.0/48.0-rv_m*rv_m/3840.0);
    }
    else {
        q[0]=cos(rv_m/2.0);
        s=sin(rv_m/2.0)/rv_m;
    }

    q[1]=s*rv1[0];
    q[2]=s*rv1[1];
    q[3]=s*rv1[2];
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
                   

