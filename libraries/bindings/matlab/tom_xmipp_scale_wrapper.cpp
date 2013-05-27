/*=================================================================
 *
 * tom_xmipp_scale is a wrapper to xmipp_scale
 *
 * The calling syntax is:
 *
 *		im_out = tom_xmipp_scale_wrapper(image,outsize,gridding)
 *
 * Electron Tomography toolbox of the
 * Max-Planck-Institute for Biochemistry
 * Dept. Molecular Structural Biology
 * 82152 Martinsried, Germany
 * http://www.biochem.mpg.de
 *
 * and
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 *
 * created: 30/08/2007
 * by: Andreas Korinek
 *
 *=================================================================*/

/*xmipp includes */
#include "xmipp_image.h"
#include "tom_xmipp_helpers.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    int outsize[3];
    const double *p_outsize=mxGetPr(prhs[1]);
    outsize[0]=(int)p_outsize[0];
    outsize[1]=(int)p_outsize[1];
    outsize[2]=(int)p_outsize[2];
    
    bool gridding = (bool)mxGetScalar(prhs[2]);
    mwSize ndims = mxGetNumberOfDimensions(prhs[0]);
    
    if (ndims == 2)
    {
        Image<double> image;
        getMatrix2D(prhs[0],image());
		selfScaleToSize(BSPLINE3, image(), outsize[0], outsize[1]);
        setMatrix2D(image(),plhs[0]);
    }
    else 
    {
        Image<double> volume;
        getMatrix3D(prhs[0],volume());
		selfScaleToSize(BSPLINE3, volume(), outsize[0], outsize[1], outsize[2]);
        setMatrix3D(volume(),plhs[0]);
    }
}	

