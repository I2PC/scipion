#include <mex.h>
#include <data/metadata.h>

/* the gateway function */
/* xmipp_nma_save_cluster(NMAdirectory,clusterName,inCluster) */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  /* check for proper number of arguments */
  if(nrhs!=3)
    mexErrMsgTxt("3 input required.");
  
  /* Get parameters */
  char nmaDir[1024];
  mxGetString(prhs[0],nmaDir,mxGetN(prhs[0])+1);
  char clusterName[256];
  mxGetString(prhs[1],clusterName,mxGetN(prhs[1])+1);
  double *ptrInCluster = mxGetPr(prhs[2]);
  
  /* Read images */
  MetaData mdImages, mdImagesOut;
  mdImages.read(((String)nmaDir)+"/images.xmd");
  mdImages.removeDisabled();

  // Fill output
  MDRow row;
  FOR_ALL_OBJECTS_IN_METADATA(mdImages)
  {
	  if (*ptrInCluster!=0)
	  {
		  mdImages.getRow(row,__iter.objId);
		  size_t id=mdImagesOut.addObject();
		  mdImagesOut.setRow(row,id);
	  }
	  ptrInCluster++;
  }
  mdImagesOut.write(((String)nmaDir)+"/images_"+clusterName+".xmd");
}
