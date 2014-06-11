/*=================================================================
 *
 * tom_xmipp_scale_pyramid is a wrapper to xmipp_scale_pyramid
 *
 * The calling syntax is:
 *
 *		psd = tom_xmipp_scale_pyramid_wrapper(image,operation,levels)
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
#include "tom_xmipp_helpers.h"
#include <data/xmipp_image.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    /*true: expand, false:reduce*/
    bool operation = (bool)mxGetScalar(prhs[1]);
    int levels = (int)mxGetScalar(prhs[2]);
    mwSize ndims = mxGetNumberOfDimensions(prhs[0]);
    
    MultidimArray<double> result;
    if (ndims == 2)
    {
        Image<double> img;
        getMatrix2D(prhs[0],img());

        if (operation == true)
        {
        	pyramidExpand(BSPLINE3,result,img(),levels);
        }
        else 
        {
        	pyramidReduce(BSPLINE3,result,img(),levels);
        }

        setMatrix2D(result,plhs[0]);
    }
    else 
    {
        Image<double> vol;
        getMatrix3D(prhs[0],vol());
       
        if (operation == true)
        {
        	pyramidExpand(BSPLINE3,result,vol(),levels);
        }
        else
        {
        	pyramidReduce(BSPLINE3,result,vol(),levels);
        }
        
        setMatrix3D(result,plhs[0]);
    }
}	

