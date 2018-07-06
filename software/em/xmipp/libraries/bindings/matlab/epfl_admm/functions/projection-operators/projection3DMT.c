// g = projection3D(cK, r1,r2, Ti, Tp, LT, supp)
//
// 3D projection in parallel-beam geometry.
//
// cK:    Input image.
// r1, r2:projection coordinate.
// Ti:    sampling step in image domain Ti[3].
// Tp:    sampling step in the projection domain Tp[2].
// LT:    lookup table.
// numberOfSample: number Of sample in lookup table.
// supp : supprot of basis function in the projection domain
#include <pthread.h> /* for threading */
#include <mex.h>
#include <matrix.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static int NTHREAD = 30;
double *LT; 
double supp, step;
double cx1, cx2, cx3, cy1, cy2;
int numberOfSample, numberOfView;


typedef struct
{
    double *Ck, *proj, *r1, *r2, *Tp, *Ti;
    int    index       ;
    int    *N,*M       ;
} tdata;



// function to read lookup table
double readLT(double u, double v) {
    double r = sqrt(u*u+v*v);
    if (r > supp) {
        return 0;
    } else {
        r = r/step ;
        int rmin = (int)(r)    ;
        rmin     = (rmin>r ? (rmin-1):rmin);
        int rmax = rmin+1    ;
        double p = r-rmin    ;
        return p * LT[rmax] + (1-p) * LT[rmin];
        //int rround = round(r);
        //return LT[rmin];
    }
}

int myFloor(double u){
    
    int x = (int) u;
    if (x<=u)
        return x  ;
    else
        return x+1;
    
    }

void *computeProjection(void * t_data)
{
    tdata *th_data = (tdata *) t_data;
    double *Ck, *proj, *r1, *r2;
    double *Tp, *Ti;
    int    index       ;
    int    *N, *M      ;
    
    Ck  = (*th_data).Ck;
    proj= (*th_data).proj;
    r1  = (*th_data).r1;
    r2  = (*th_data).r2;
    Tp  = (*th_data).Tp;
    Ti  = (*th_data).Ti;
    index= (*th_data).index;
    N  = (*th_data).N;
    M  = (*th_data).M;

    int nx1, nx2, nx3, ny1, ny2;
    int i;
    double x1n,x2n,x3n;
    double x1nr1,x1nr2;
    double x2nr1,x2nr2;
    double x3nr1,x3nr2;
    double y1,y2;
    for (i=index; i<numberOfView; i=i+NTHREAD){
       
        
        for (nx3=0; nx3<N[2]; nx3++) {
            // initial computation
            x3n  = (nx3-cx3)*Ti[2];
            x3nr1= x3n*r1[2+i*3];
            x3nr2= x3n*r2[2+i*3];
         
            for (nx1=0; nx1<N[1]; nx1++) {
                // initial computation
                x1n  = (nx1-cx1)*Ti[1];
                x1nr1= x1n*r1[i*3]  ;
                x1nr2= x1n*r2[i*3]  ;
                
                for (nx2=0; nx2<N[0]; nx2++) {
                    // initial computation
                    x2n  = (nx2-cx2)*Ti[0];
                    x2nr1= x2n*r1[1+i*3];
                    x2nr2= x2n*r2[1+i*3];
                    
                    y1 = x1nr1+x2nr1+x3nr1;
                    y2 = x1nr2+x2nr2+x3nr2;
                    
                    
                    
                    int ny1min = (int)(((y1-supp*Ti[0])/Tp[1])+cy1);
                    ny1min = ((ny1min)>(((y1-supp*Ti[0])/Tp[1])+cy1) ? ny1min:(ny1min+1));
                    int ny1max = (int)(((y1+supp*Ti[0])/Tp[1])+cy1);
                    ny1max = (ny1max>(((y1+supp*Ti[0])/Tp[1])+cy1) ? (ny1max-1):ny1max)  ;
                    int ny2min = (int)(((y2-supp*Ti[0])/Tp[0])+cy2);
                    ny2min = (ny2min>(((y2-supp*Ti[0])/Tp[0])+cy2) ? ny2min:(ny2min+1))  ;
                    int ny2max = (int)(((y2+supp*Ti[0])/Tp[0])+cy2);
                    ny2max = (ny2max>(((y2+supp*Ti[0])/Tp[0])+cy2) ? (ny2max-1):ny2max)  ;
                    
                    
                    if (ny1min<0)
                        ny1min = 0;
                    if (ny1max>(M[1]-1))
                        ny1max = M[1]-1;
                    if (ny2min<0)
                        ny2min = 0;
                    if (ny2max>(M[0]-1))
                        ny2max = M[0]-1;
                    
                    for (ny1=ny1min; ny1<=ny1max; ny1++) {
                        for (ny2=ny2min; ny2<=ny2max; ny2++) {
                                proj[i*((int) M[0]*(int) M[1])+ny1*((int) M[0])+ny2] += Ti[0]*Ck[(nx3*N[1]+nx1)*N[0]+nx2]* readLT((Tp[1]*(ny1-cy1)-y1)/Ti[0], (Tp[0]*(ny2-cy2)-y2)/Ti[0]);
                        }
                    }
                }
            }
        }
        
    }
    
}

// mex file to implent the forward projection
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mwSize *N;
    const mwSize *N1;
    const mwSize *N2;
    int M[2];
    double *Ck, *proj, *r1, *r2;
    double *Tp;
    double *Ti;
    // Check for proper number of input arguments
    if (nrhs!=7) mexErrMsgTxt("There must be exactly 7 input arguments.");
    if (nlhs!=1) mexErrMsgTxt("There must be exactly 1 output argument.");
    
    // Check for proper type of input arguments
    if (!mxIsDouble(prhs[0]) | mxIsComplex(prhs[0])) mexErrMsgTxt("input image must be a real-valued double-precision array.");
    if (mxGetNumberOfDimensions(prhs[0])!=3) mexErrMsgTxt("input image must be a 3D array.");
    if (mxGetNumberOfElements(prhs[0])==0) mexErrMsgTxt("input image cannot be empty.");
    N = mxGetDimensions(prhs[0]);
    
    N1  = mxGetDimensions(prhs[1]);
    if (N1[0]!=3) mexErrMsgTxt("input r1 must be 2D vector whose first dimension is 3.");

    N2  = mxGetDimensions(prhs[2]);
    if (N2[0]!=3) mexErrMsgTxt("input r2 must be 2D vector whose first dimension is 3.");

    //mexPrintf("%d %d",N1[0],N1[1]);
    if (N2[1]!=N1[1]) mexErrMsgTxt("input r1 and r2 must be 2D vector whose second dimension should be the same.");

    numberOfView = (int) N1[1];
     
    if (!mxIsDouble(prhs[3]) | mxGetNumberOfElements(prhs[3])!=3) mexErrMsgTxt("Ti must be a double-precision 3D vector.");
    Ti = mxGetPr(prhs[3]);
    
    if (!mxIsDouble(prhs[4]) | mxGetNumberOfElements(prhs[4])!=2) mexErrMsgTxt("Tp must be a double-precision 2D vector.");
    Tp = mxGetPr(prhs[4]);
    
    if (!mxIsDouble(prhs[5]) | mxIsComplex(prhs[5])) mexErrMsgTxt("lookup table must be a real-valued double-precision array.");
    if (mxGetNumberOfDimensions(prhs[5])!=2) mexErrMsgTxt("lookup table must be a 1D array.");
    numberOfSample = mxGetNumberOfElements(prhs[5]);
    if (numberOfSample==0) mexErrMsgTxt("lookup table cannot be empty.");
    
    if (!mxIsDouble(prhs[6]) | mxGetNumberOfElements(prhs[6])!=1 | prhs[6]<=0) mexErrMsgTxt("supp must be a positive double-precision scalar.");
    supp = mxGetScalar(prhs[6])             ;
    step = supp/((double)(numberOfSample-1));
    
    double p = sqrt(N[0]*N[0]*Ti[0]+N[1]*N[1]*Ti[1]+N[2]*N[2]*Ti[2]);
    M[0] = ceil((p+2*supp*Ti[0])/Tp[0]);
    M[1] = ceil((p+2*supp*Ti[0])/Tp[1]);
    if (M[0] % 2 !=0)
        M[0] = M[0]+1;
    
    if (M[1] % 2 !=0)
        M[1] = M[1]+1;
   
    // Allocate output variable
    const mwSize Mprm[3] = {(mwSize) M[0], (mwSize) M[1], (mwSize) numberOfView };
    plhs[0] = mxCreateNumericArray(3, Mprm, mxDOUBLE_CLASS, mxREAL);
    if (!plhs[0]) mexErrMsgTxt("Could not allocate memory for output variable.");
    
    
    // Center of object
    cx1 = ((double) N[1]-1)/2;
    cx2 = ((double) N[0]-1)/2;
    cx3 = ((double) N[2]-1)/2;
    
    // Center of image
    cy1 = ((double) M[1]-1)/2;
    cy2 = ((double) M[0]-1)/2;
    
    // Tomographic projection
    Ck = mxGetPr(prhs[0]);
    r1 = mxGetPr(prhs[1]);
    r2 = mxGetPr(prhs[2]);
    LT = mxGetPr(prhs[5]);
    proj = mxGetPr(plhs[0]);
    
    
    
    // handle threads
    pthread_t threads[NTHREAD];
    
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
    tdata *data[NTHREAD];

    
    int t;
    for (t=0 ; t<NTHREAD ; t++){
        data[t] = (tdata*)  mxCalloc( 1,sizeof(tdata) );
        (*data[t]).Ck     =  Ck  ;
        (*data[t]).proj   =  proj;
        (*data[t]).r1     =  r1  ;
        (*data[t]).r2     =  r2  ;
        (*data[t]).Tp     =  Tp  ;
        (*data[t]).Ti     =  Ti  ;
        (*data[t]).index  =  t   ;
        (*data[t]).N      =  (int *)N;
        (*data[t]).M      =  (int *)M;
  
    }
    
    int rc;
    void * status;
    
    for (t=0; t < NTHREAD ; t++)
    {
        rc = pthread_create(&threads[t], &attr, computeProjection,(void*) data[t]);
    }
    
    
    // mexPrintf("threads created : waiting for end\n");
    // wait for the threads to finish
    
    pthread_attr_destroy(&attr);
    
    
   
    
    for ( t=0; t < NTHREAD; t++)
    {
        rc = pthread_join(threads[t], &status);
        //      if (rc)
        //     {
        //         mexPrintf("error %i joining thread %i\n",rc,t);
        //     } else {
        //         mexPrintf("thread %i joined successfully\n",t);
        //    }
        
    }
    
}