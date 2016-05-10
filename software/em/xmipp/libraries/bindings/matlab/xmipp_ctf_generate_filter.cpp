#include <mex.h>
#include <data/ctf.h>

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /* check for proper number of arguments */
    if(nrhs!=15)
        mexErrMsgTxt("15 inputs required.");
    if(nlhs!=1)
        mexErrMsgTxt("1 output required.");

    /* Get parameters */
    CTFDescription ctf;
    int Xdim=(int)mxGetScalar(prhs[0]);
    ctf.Tm=mxGetScalar(prhs[1]);
    ctf.K=mxGetScalar(prhs[2]);
    ctf.kV=mxGetScalar(prhs[3]);
    ctf.DeltafU=mxGetScalar(prhs[4]);
    ctf.DeltafV=mxGetScalar(prhs[5]);
    ctf.azimuthal_angle=mxGetScalar(prhs[6]);
    ctf.Cs=mxGetScalar(prhs[7]);
    ctf.Ca=mxGetScalar(prhs[8]);
    ctf.espr=mxGetScalar(prhs[9]);
    ctf.ispr=mxGetScalar(prhs[10]);
    ctf.alpha=mxGetScalar(prhs[11]);
    ctf.DeltaF=mxGetScalar(prhs[12]);
    ctf.DeltaR=mxGetScalar(prhs[13]);
    ctf.Q0=mxGetScalar(prhs[14]);
    
    ctf.produceSideInfo();

    /*  set the output pointer to the output matrix */
    mwSize dims[2];
    dims[0] = Xdim;
    dims[1] = Xdim;

    MultidimArray<double> ctfImage;
    ctfImage.initZeros(Xdim,Xdim);
    ctf.generateCTF(Xdim, Xdim, ctfImage);
    
    plhs[0]=mxCreateNumericArray((mwSize) 2, dims, mxDOUBLE_CLASS, mxREAL);
    double *ptrMatlab=mxGetPr(plhs[0]);
    memcpy(ptrMatlab,MULTIDIM_ARRAY(ctfImage),Xdim*Xdim*sizeof(double));
}
