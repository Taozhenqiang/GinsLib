/*------------------------------------------------------------------------------
* lambda.c : integer ambiguity resolution
*
*          Copyright (C) 2007-2008 by T.TAKASU, All rights reserved.
*
* reference :
*     [1] P.J.G.Teunissen, The least-square ambiguity decorrelation adjustment:
*         a method for fast GPS ambiguity estimation, J.Geodesy, Vol.70, 65-82,
*         1995
*     [2] X.-W.Chang, X.Yang, T.Zhou, MLAMBDA: A modified LAMBDA method for
*         integer least-squares estimation, J.Geodesy, Vol.79, 552-565, 2005
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/13 1.0 new
*           2015/05/31 1.1 add api lambda_reduction(), lambda_search()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants/macros ----------------------------------------------------------*/

#define LOOPMAX     10000           /* maximum count of search loop */

#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define ROUND(x)    (floor((x)+0.5))
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)

/* LD factorization (Q=L'*diag(D)*L) -----------------------------------------*/
static int LD(int n, const double *Q, double *L, double *D)
{
    int i,j,k,info=0;
    double a,*A=mat(n,n);
    
    memcpy(A,Q,sizeof(double)*n*n);
    for (i=n-1;i>=0;i--) {
        if ((D[i]=A[i+i*n])<=0.0) {info=-1; break;}
        a=sqrt(D[i]);
        for (j=0;j<=i;j++) L[i+j*n]=A[i+j*n]/a;
        for (j=0;j<=i-1;j++) for (k=0;k<=j;k++) A[j+k*n]-=L[i+k*n]*L[i+j*n];
        for (j=0;j<=i;j++) L[i+j*n]/=L[i+i*n];
    }
    /* trace(12,"Q=\n"); tracemat(12,Q,n,n,10,5,0);
    trace(12,"D=\n"); tracemat(12,D,n,1,10,5,0);
    trace(12,"L=\n"); tracemat(12,L,n,n,10,5,0); */
    free(A);
    if (info) fprintf(stderr,"%s : LD factorization error\n",__FILE__);
    return info;
}
/* integer gauss transformation ----------------------------------------------*/
static void gauss(int n, double *L, double *Z, int i, int j)
{
    int k,mu;
    
    if ((mu=(int)ROUND(L[i+j*n]))!=0) {
        for (k=i;k<n;k++) L[k+n*j]-=(double)mu*L[k+i*n];
        for (k=0;k<n;k++) Z[k+n*j]-=(double)mu*Z[k+i*n];
    }
}
/* permutations --------------------------------------------------------------*/
static void perm(int n, double *L, double *D, int j, double del, double *Z)
{
    int k;
    double eta,lam,a0,a1;
    
    eta=D[j]/del;
    lam=D[j+1]*L[j+1+j*n]/del;
    D[j]=eta*D[j+1]; D[j+1]=del;
    for (k=0;k<=j-1;k++) {
        a0=L[j+k*n]; a1=L[j+1+k*n];
        L[j+k*n]=-L[j+1+j*n]*a0+a1;
        L[j+1+k*n]=eta*a0+lam*a1;
    }
    L[j+1+j*n]=lam;
    for (k=j+2;k<n;k++) SWAP(L[k+j*n],L[k+(j+1)*n]);
    for (k=0;k<n;k++) SWAP(Z[k+j*n],Z[k+(j+1)*n]);
}
/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
static void reduction(int n, double *L, double *D, double *Z)
{
    int i,j,k;
    double del;
    
    j=n-2; k=n-2;
    while (j>=0) {
        if (j<=k) for (i=j+1;i<n;i++) gauss(n,L,Z,i,j);
        del=D[j]+L[j+1+j*n]*L[j+1+j*n]*D[j+1];
        if (del+1E-6<D[j+1]) { /* compared considering numerical error */
            perm(n,L,D,j,del,Z);
            k=j; j=n-2;
        }
        else j--;
    }
}
/* modified lambda (mlambda) search (ref. [2]) -------------------------------
* args   : n      I  number of float parameters
*          m      I  number of candidate fixed solution
           L,D    I  transformed covariance matrix
           zs     I  transformed double-diff float solutions
           zn     O  transformed double-diff fixed solutions
           s      O  sum of residuals for fixed solutions                    */
static int search(int n, int m, const double *L, const double *D,
                  const double *zs, double *zn, double *s)
{
    int i,j,k,l,c,nn=0,imax=0;
    double newdist,maxdist=1E99,y;
    double *S=zeros(n,n),*dist=mat(n,1),*zb=mat(n,1),*z=mat(n,1),*step=mat(n,1);
    
    k=n-1; dist[k]=0.0;
    zb[k]=zs[k];
    z[k]=ROUND(zb[k]);
    y=zb[k]-z[k];
    step[k]=SGN(y);  /* step towards closest integer */
    for (c=0;c<LOOPMAX;c++) {
        newdist=dist[k]+y*y/D[k];  /* newdist=sum(((z(j)-zb(j))^2/d(j))) */
        if (newdist<maxdist) {
            /* Case 1: move down */
            if (k!=0) {
                dist[--k]=newdist;
                for (i=0;i<=k;i++)
                    S[k+i*n]=S[k+1+i*n]+(z[k+1]-zb[k+1])*L[k+1+i*n];
                zb[k]=zs[k]+S[k+k*n];
                z[k]=ROUND(zb[k]); /* next valid integer */
                y=zb[k]-z[k];
                step[k]=SGN(y);
            }
            /* Case 2: store the found candidate and try next valid integer */
            else {
                /* store the first m initial points, execute only once */
                if (nn<m) {  
                    if (nn==0||newdist>s[imax]) imax=nn;
                    for (i=0;i<n;i++) zn[i+nn*n]=z[i];
                    s[nn++]=newdist;
                    if (nn==m-1) { trace(12,"initial integer candidates Zn:\n"); tracemat(12,zn,nn,n,10,5,0); }
                }
                else {
                    if (newdist<s[imax]) {
                        for (i=0;i<n;i++) zn[i+imax*n]=z[i];
                        s[imax]=newdist;
                        for (i=imax=0;i<m;i++) if (s[imax]<s[i]) imax=i;
                    }
                    maxdist=s[imax];
                }
                z[0]+=step[0]; /* next valid integer */
                y=zb[0]-z[0];
                step[0]=-step[0]-SGN(step[0]);
            }
        }
        /* Case 3: exit or move up */
        else {
            if (k==n-1) break;
            else {
                k++;  /* move up */
                z[k]+=step[k];  /* next valid integer */
                y=zb[k]-z[k];
                step[k]=-step[k]-SGN(step[k]);
            }
        }
    }

    /* sort the searched set of integer candidates based on quadratic residuals (s) */
    for (i=0;i<m-1;i++) { 
        for (j=i+1;j<m;j++) {
            if (s[i]<s[j]) continue;
            SWAP(s[i],s[j]);
            for (l=0;l<n;l++) SWAP(zn[l+i*n],zn[l+j*n]);
        }
    }
    free(S); free(dist); free(zb); free(z); free(step);
    
    if (c>=LOOPMAX) {
        /* fprintf(stderr,"%s : search loop count overflow\n",__FILE__); */
        trace(12,"search loop count overflow\n");
        /* trace(12,"not ok, integer candidates Zn:\n"); tracemat(12,zn,m,n,10,5,0); */
        return -2;
    }
    else {
        /* trace(12,"ok, integer candidates Zn:\n"); tracemat(12,zn,m,n,10,5,0); */        
    }



    return 0;
}
/* lambda/mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).
* args   : int    n      I  number of float parameters
*          int    m      I  number of candidate fixed solutions
*          double *a     I  float parameters (n x 1) (double-diff phase biases)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
* notes  : matrix stored by column-major order (fortran convension)
*-----------------------------------------------------------------------------*/
extern int lambda(rtk_t *rtk, int n, int m, const double *a, const double *Q, double *F,
                  double *s, const int *ix, const int *ixf, int *low_ix)
{
    int info,i,j;
    double *L,*D,*Z,*z,*E,*zQ,*Qz,*dQz;
    
    if (n<=0||m<=0) return -1;
    L=zeros(n,n); D=mat(n,1); Z=eye(n); z=mat(n,1); E=mat(n,m); zQ=mat(n,n); Qz=mat(n,n); dQz=mat(n,1);
    
    /* LD (lower diagonal) factorization (Q=L'*diag(D)*L) */
    if (!(info=LD(n,Q,L,D))) {
        
        /* lambda decorrelation (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) */
        reduction(n,L,D,Z);

        /* decorrelated ambiguity covariance matrix */
        matmul("NN",n,n,n,Z,Q,zQ,1.0,0.0);
        matmul("NT",n,n,n,zQ,Z,Qz,1.0,0.0);
        /* trace(12,"Qz=\n");tracemat(12,Qz,n,n,10,5,0); */
        /* trace(12,"Z'=\n"); tracemat(12,Z,n,n,10,5,0);
        trace(12,"Q=\n"); tracemat(12,Q,n,n,10,5,0);
        trace(12,"L=\n"); tracemat(12,L,n,n,10,5,0);
        trace(12,"D=\n"); tracemat(12,D,n,1,10,5,0); */
        
        /* the ambiguity covariance after decorrelation is used as a screening criterion for PAR */
        for (i=0;i<n;i++) dQz[i]=Qz[i+i*n]; /* dQz=diag(Qz) */
        for (i=1,j=0;i<n;i++) {
            if (dQz[i]>dQz[j]) j=i;
        }

        /* index of maximum variance */
        /* low_ix[0]=ix[2*j+1]-(ixf[2*j]*MAXSAT+rtk->na); low_ix[1]=ixf[2*j+1]; */ 

        matmul("NN",n,n,1,Z,a,z,1.0,0.0);
        /* trace(12,"z=\n"); tracemat(12,z,n,1,7,2,0); */
        /* matmul("TN",n,1,n,Z,a,z); */ /* z=Z'*a */
        
        /* mlambda search 
            z = transformed double-diff phase biases
            L,D = transformed covariance matrix */
        if (!(info=search(n,m,L,D,z,E,s))) {  /* returns 0 if no error */
            
            /* transform the fixed integer ambiguity to the original space, F=Z'\E */
            /* the fixed solution here is the row vector, so F(mxn)=(Z'\E)'=E'*Z^-1 */
            info=solve("T",Z,E,n,m,F); 
        }
    }
    free(L); free(D); free(Z); free(z); free(E); free(zQ); free(Qz); free(dQz);
    return info;
}
/* lambda reduction ------------------------------------------------------------
* reduction by lambda (ref [1]) for integer least square
* args   : int    n      I  number of float parameters
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *Z     O  lambda reduction matrix (n x n)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_reduction(int n, const double *Q, double *Z)
{
    double *L,*D;
    int i,j,info;
    
    if (n<=0) return -1;
    
    L=zeros(n,n); D=mat(n,1);
    
    for (i=0;i<n;i++) for (j=0;j<n;j++) {
        Z[i+j*n]=i==j?1.0:0.0;
    }
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D);
        return info;
    }
    /* lambda reduction */
    reduction(n,L,D,Z);
     
    free(L); free(D);
    return 0;
}
/* mlambda search --------------------------------------------------------------
* search by  mlambda (ref [2]) for integer least square
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_search(int n, int m, const double *a, const double *Q,
                         double *F, double *s)
{
    double *L,*D;
    int info;
    
    if (n<=0||m<=0) return -1;
    
    L=zeros(n,n); D=mat(n,1);
    
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D);
        return info;
    }
    /* mlambda search */
    info=search(n,m,L,D,a,F,s);
    
    free(L); free(D);
    return info;
}

/* calculate the matrix determinant */
extern double determinant(const double* Qb, int n) 
{
    double *matrix=mat(n,n);
    double det=1.0,temp;
    int i,j,k;

    /* 复制输入矩阵到工作矩阵 */
    matcpy(matrix,Qb,n,n);

    if (matrix==NULL||n<=0) {
        fprintf(stderr, "Error: Invalid matrix or dimension.\n");
        return -1;  
    }

    /* 使用高斯消元法计算行列式 */
    for (i=0;i<n;i++) {
        /* 找到当前列最大的元素（进行列主元选择）*/
        int max_row=i;
        for (j=i+1;j<n;j++) {
            if (fabs(matrix[j*n+i])>fabs(matrix[max_row*n+i])) {
                max_row=j;
            }
        }

        /* 如果主元为0，则行列式为0 */
        if (matrix[max_row*n+i]==0) {
            return 0.0;
        }

        /* 行交换（交换当前行和最大元素所在的行）*/
        if (i!=max_row) {
            for (k=0;k<n;k++) {
                temp=matrix[i*n+k];
                matrix[i*n+k]=matrix[max_row*n+k];
                matrix[max_row*n+k]=temp;
            }
            det=-det;  /* 每次交换行，行列式符号改变 */
        }

        /* 消去当前列下方元素 */
        for (j=i+1;j<n;j++) {
            double factor=matrix[j*n+i]/matrix[i*n+i];
            for (k=i;k<n;k++) {
                matrix[j*n+k]-=matrix[i*n+k]*factor;
            }
        }

        /* 乘上对角线元素 */
        det*=matrix[i*n+i];
    }

    free(matrix);

    return det;
}

extern int amb_BIE_qc(rtk_t *rtk, const double *Qab, const double *Qb, const double *y, const double *b, 
                     int na, int nb, int num_candidate, int mode, double *b_BIE, double *temp) 
{
    int i,j,vnum,nx=rtk->nx;
    double gamma,sum_p,delta_b,res_da=0.0;
    double *QaIb,*IQb,*Qa,*db,*da;

    IQb=mat(nb,nb); QaIb=mat(na,nb); Qa=mat(na,na); db=mat(nb,1); da=mat(na,1);

    /* quality control for BIE */
    if (BIE_amb_rec==mode) {
        matcpy(IQb,Qb,nb,nb);
        if (matinv(IQb,nb)) {
            trace(12,"resamb_LAMBDA: floating ambiguity covariance matrix inversion error\n");
            return 0;
        }        
    }

    /* QaIb=Qab*Qb^-1 */
    for (j=0;j<nb;j++) db[j]=y[j]-b[j];
    if (BIE_amb_rec==mode) {
        matmul("NN",na,nb,nb,Qab,IQb,QaIb,1.0,0.0);
        for (i=0;i<na;i++) Qa[i+i*na]=rtk->P[i+i*nx];
        /* da=Qab*Qb^-1*(b0-b)*/
        matmul("NN",na,nb,1,QaIb,db,da,1.0,0.0);  
        res_da=quadratic(da,Qa,na);      
    }
    else res_da=0.0;

    /* overall quality control, test threshold */
    gamma=rtk->sol.thres*rtk->sol.thres*(quadratic(db,Qb,nb)+res_da);
    vnum=num_candidate;
    while (1) {
        for (i=0,sum_p=0.0;i<vnum;i++) {
            for (j=0;j<nb;j++) db[j]=y[j]-b[j+nb*i];
            if (BIE_amb_rec==mode) {
                matmul("NN",na,nb,1,QaIb,db,da,1.0,0.0);  
                res_da=quadratic(da,Qa,na);              
            }
            sum_p+=(quadratic(db,Qb,nb)+res_da);
        }
        if (sum_p>gamma*vnum) {
            vnum--;
            continue;
        }
        else break;
    }

    /* BIE soluiton */
    for (i=0,sum_p=0.0;i<vnum;i++) {
        for (j=0;j<nb;j++) db[j]=y[j]-b[j+nb*i];
        if (BIE_amb_rec==mode) {
            matmul("NN",na,nb,1,QaIb,db,da,1.0,0.0);  
            res_da=quadratic(da,Qa,na);              
        }
        sum_p+=exp(-0.5*(quadratic(db,Qb,nb)+res_da));
    }
    for (i=0;i<vnum;i++) {
        for (j=0;j<nb;j++) db[j]=y[j]-b[j+nb*i];
        if (BIE_amb_rec==mode) {
            matmul("NN",na,nb,1,QaIb,db,da,1.0,0.0);  
            res_da=quadratic(da,Qa,na);              
        }
        for (j=0;j<nb;j++) b_BIE[j]+=b[j+nb*i]*exp(-0.5*(quadratic(db,Qb,nb)+res_da))/sum_p;
    }
    /* round BIE solution to integer */
    for (i=0;i<nb;i++) {
        delta_b=fabs(b_BIE[i]-round(b_BIE[i]));
        if (delta_b<0.1) b_BIE[i]=round(b_BIE[i]);
    }

    /* residual quadratic of BIE solution */
    for (i=0;i<nb;i++) db[i]=y[i]-b_BIE[i]; 
    if (BIE_amb_rec==mode) {
        matmul("NN",na,nb,1,QaIb,db,da,1.0,0.0);  
        res_da=quadratic(da,Qa,na);              
    }
    temp[0]=quadratic(db,Qb,nb)+res_da;

    /* residual quadratic of ILS solution */
    for (i=0;i<nb;i++) db[i]=y[i]-b[i]; 
    if (BIE_amb_rec==mode) {
        matmul("NN",na,nb,1,QaIb,db,da,1.0,0.0);  
        res_da=quadratic(da,Qa,na);              
    }  
    temp[1]=quadratic(db,Qb,nb)+res_da;

    trace(12,"BIE/chi=%.4f BIE/ILS=%.3f BIE=\n",temp[0]/chisqr[na+nb-1],temp[0]/temp[1]); tracemat(12,b_BIE,1,nb,7,2,0); 

    free(QaIb); free(IQb); free(Qa); free(db); free(da);

    return 1;
}