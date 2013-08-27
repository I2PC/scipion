#include <mex.h>
#include <data/metadata.h>

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  /* check for proper number of arguments */
  if(nrhs!=1) 
    mexErrMsgTxt("1 input required.");
  if(nlhs!=3)
    mexErrMsgTxt("3 outputs are required.");
  
  /* Get parameters */
  char nmaDir[1024];
  mxGetString(prhs[0],nmaDir,mxGetN(prhs[0])+1);
  
  /* Read images */
  MetaData mdImages;
  mdImages.read(((String)nmaDir)+"/images.xmd");
  mdImages.removeDisabled();
  int nImgs=(int)mdImages.size();

  /* Read modes */
  MetaData mdModes;
  mdModes.read(((String)nmaDir)+"/modes.xmd");
  mdModes.removeDisabled();
  int nModes=(int)mdModes.size();

  // Allocate output
  mwSize dims[1];
  dims[0]=(mwSize)nImgs;
  plhs[0]=mxCreateCellArray((mwSize)1, dims);
  plhs[1]=mxCreateDoubleMatrix((mwSize)nImgs, (mwSize)nModes, mxREAL);
  double *ptrNMADistplacements=mxGetPr(plhs[1]);
  plhs[2]=mxCreateDoubleMatrix((mwSize)nImgs, (mwSize)1, mxREAL);
  double *ptrCost=mxGetPr(plhs[2]);

  // Fill output
  int i=0;
  String fnImg;
  std::vector<double> lambda;
  FOR_ALL_OBJECTS_IN_METADATA(mdImages)
  {
	  mdImages.getValue(MDL_IMAGE,fnImg,__iter.objId);
	  mxSetCell(plhs[0], i, mxCreateString(fnImg.c_str()));
	  mdImages.getValue(MDL_NMA,lambda,__iter.objId);
	  for (int j=0; j<nModes; ++j)
		  ptrNMADistplacements[j*nImgs+i]=lambda[j]; // x*Ydim+y
	  mdImages.getValue(MDL_COST,*ptrCost,__iter.objId);
	  i++;
	  ptrCost++;
  }
}
