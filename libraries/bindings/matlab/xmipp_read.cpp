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
  dims[0]=XSIZE(mI);
  dims[1]=YSIZE(mI);
  dims[2]=ZSIZE(mI);
  dims[3]=NSIZE(mI);

  plhs[0]=mxCreateNumericArray((mwSize) mI.getDim(), dims, mxDOUBLE_CLASS, mxREAL);
  double *Cp=mxGetPr(plhs[0]);
  
  memcpy(Cp,MULTIDIM_ARRAY(mI),MULTIDIM_SIZE(mI)*sizeof(double));
}
