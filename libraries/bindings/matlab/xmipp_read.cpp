#include <mex.h>
#include <data/xmipp_image.h>

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /* check for proper number of arguments */
    if(nrhs!=1)
        mexErrMsgTxt("1 input required.");
    if(nlhs!=1)
        mexErrMsgTxt("1 output required.");

    /* Get parameters */
    char fnImg[1024];
    mxGetString(prhs[0],fnImg,mxGetN(prhs[0])+1);

    /* Read image */
    Image<double> I;
    I.read(fnImg);
    MultidimArray<double> &mI=I();

    /*  set the output pointer to the output matrix */
    mwSize dims[4];
    dims[0] = YSIZE(mI);
    dims[1] = XSIZE(mI);
    dims[2] = ZSIZE(mI);
    dims[3] = NSIZE(mI);

    plhs[0]=mxCreateNumericArray((mwSize) mI.getDim(), dims, mxDOUBLE_CLASS, mxREAL);

    double *ptrMatlab=mxGetPr(plhs[0]);
    double *ptrC=MULTIDIM_ARRAY(mI);

    for (size_t n=0; n<NSIZE(mI); ++n)
        for (size_t z=0; z<ZSIZE(mI); ++z)
        {
            for (size_t y=0; y<YSIZE(mI); ++y)
                for (size_t x=0; x<XSIZE(mI); ++x, ptrC++)
                    ptrMatlab[x*YSIZE(mI)+y]=*ptrC;
            ptrMatlab+=YXSIZE(mI);
        }
}
