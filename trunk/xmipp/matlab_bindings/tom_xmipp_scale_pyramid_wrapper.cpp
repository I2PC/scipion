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
#include "image.h"
#include "volume.h"
#include "tom_xmipp_helpers.h"

/*Matlab includes*/
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    /*true: expand, false:reduce*/
    bool operation = (bool)mxGetScalar(prhs[1]);
    int levels = (int)mxGetScalar(prhs[2]);
    mwSize ndims = mxGetNumberOfDimensions(prhs[0]);
    
    if (ndims == 2)
    {
        ImageXmipp img;
        getMatrix2D(prhs[0],img());

        float Xoff, Yoff;
        img.get_originOffsets(Xoff, Yoff);
        float scale_factor = (float)(pow(2.0, levels));
        
        Matrix2D<double> result;

        if (operation == true)
        {
            img().pyramidExpand(result, levels);
            img.set_originOffsets(Xoff*scale_factor, Yoff*scale_factor);
        }
        else 
        {
            img().pyramidReduce(result, levels);
            img.set_originOffsets(Xoff / scale_factor, Yoff / scale_factor);
        }

        setMatrix2D(result,plhs[0]);
    }
    else 
    {
        Volume vol;
        getMatrix3D(prhs[0],vol());
        Matrix3D<double> result; 
       
        if (operation == true)
        {
            vol().pyramidExpand(result, levels);
        }
        else
        {
            vol().pyramidReduce(result, levels);
        }
        
        setMatrix3D(result,plhs[0]);
    }
}	

