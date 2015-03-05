#include <mex.h>
#include <data/xmipp_image.h>

/* the gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
	/* Get parameters */
	char fnImg[1024];
	mxGetString(prhs[1],fnImg,mxGetN(prhs[1])+1);

	/* Let's extract the size of input data */
	const mwSize *dims;
 	dims = mxGetDimensions(prhs[0]);
	int number_of_dims = mxGetNumberOfDimensions(prhs[0]);
	size_t X,Y,Z,N;
        X=dims[1];//matlab oddities
        Y=dims[0];
        Z=N=1;
	if(number_of_dims>2)	
	   Z=dims[2];
	if(number_of_dims>3)	
	   N=dims[3];

	/* Let's get the pointer to the input data (the numbers are converted to double) */
    	double *ptrMatlab=mxGetPr(prhs[0]);
    	Image<double> I;
    	MultidimArray<double> &mI=I();
    	mI.resizeNoCopy(N,Z,Y,X);
    	double *ptrC=MULTIDIM_ARRAY(mI);


	for (size_t n=0; n<NSIZE(mI); ++n)
            for (size_t z=0; z<ZSIZE(mI); ++z)
            {
            	for (size_t x=0; x<X; ++x)
                    for (size_t y=0; y<Y; ++y, ptrMatlab++)
                       ptrC[y*XSIZE(mI)+x]=*ptrMatlab;
            ptrC+=YXSIZE(mI);
           }
    I.write(fnImg);
}
