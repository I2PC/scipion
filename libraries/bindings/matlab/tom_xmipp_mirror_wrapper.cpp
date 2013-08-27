/*=================================================================
 *
 * tom_xmipp_mirror is a wrapper to xmipp_mirror
 *
 * The calling syntax is:
 *
 *		im_out = tom_xmipp_mirror_wrapper(image,[x,y,z])
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
   
    bool flipx;
    bool flipy;
    bool flipz;
    const double *p_flip=mxGetPr(prhs[1]);
    flipx=(bool)p_flip[0];
    flipy=(bool)p_flip[1];
    flipz=(bool)p_flip[2];

    mwSize ndims = mxGetNumberOfDimensions(prhs[0]);
    
    if (ndims == 2)
    {
        Image<double> img;
        getMatrix2D(prhs[0],img());
        if (flipx) img().selfReverseX();
        if (flipy) img().selfReverseY();
        setMatrix2D(img(),plhs[0]);     
    }
    else
    {
        Image<double> vol;
        getMatrix3D(prhs[0],vol());
        if (flipx) vol().selfReverseX();
        if (flipy) vol().selfReverseY();
        if (flipz) vol().selfReverseZ();
        setMatrix3D(vol(),plhs[0]);     
    }
}
