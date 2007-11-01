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
#include "image.h"
#include "volume.h"
#include "tom_xmipp_helpers.h"
#include "gridding.h"

/*Matlab includes*/
#include "mex.h"
#include "matrix.h"

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
        Image image;
        getMatrix2D(prhs[0],image());
        
        if (gridding)
        {
            Matrix2D<double> A(3, 3);
            A.initIdentity();
            KaiserBessel kb;
            Matrix2D<double> Maux;
            produceReverseGriddingMatrix2D(image(),Maux,kb);
            DIRECT_MAT_ELEM(A, 0, 0) = (double) outsize[0] / (double) XSIZE(image());
            DIRECT_MAT_ELEM(A, 1, 1) = (double) outsize[1] / (double) YSIZE(image());
            applyGeometryReverseGridding(image(), A, Maux, kb, IS_NOT_INV, WRAP, outsize[0], outsize[1]);
        }
        else
        {
            image().selfScaleToSizeBSpline(3, outsize[0], outsize[1]);
	    }
        setMatrix2D(image(),plhs[0]);
    }
    else 
    {
        Volume volume;
        getMatrix3D(prhs[0],volume());

        if (gridding)
        {
            Matrix2D<double> B(4, 4);
            B.initIdentity();
            KaiserBessel kb;
            Matrix3D<double> Maux;
            produceReverseGriddingMatrix3D(volume(),Maux,kb);
            DIRECT_MAT_ELEM(B, 0, 0) = (double) outsize[0] / (double) XSIZE(volume());
            DIRECT_MAT_ELEM(B, 1, 1) = (double) outsize[1] / (double) YSIZE(volume());
            DIRECT_MAT_ELEM(B, 2, 2) = (double) outsize[2] / (double) ZSIZE(volume());
            applyGeometryReverseGridding(volume(), B, Maux, kb, IS_NOT_INV, WRAP, outsize[0], outsize[1], outsize[2]);
        }
        else
        {
            volume().selfScaleToSizeBSpline(3, outsize[0], outsize[1], outsize[2]);
	    }
        setMatrix3D(volume(),plhs[0]);
    }

}	

