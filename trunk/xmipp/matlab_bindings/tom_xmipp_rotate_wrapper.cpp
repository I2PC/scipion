/*=================================================================
 *
 * tom_xmipp_rotate is a wrapper to xmipp_rotate
 *
 * The calling syntax is:
 *
 *		im_out = tom_xmipp_rotate_wrapper(image,angs,axis,mode,gridding,wrap)
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
   
    double angs[3];
    const double *p_angs=mxGetPr(prhs[1]);
    angs[0]=(double)p_angs[0];
    angs[1]=(double)p_angs[1];
    angs[2]=(double)p_angs[2];

    Matrix1D<double> axis;
    getMatrix1D(prhs[2],axis);

    Matrix2D<double> A3D, A2D;
    
    bool gridding = (bool)mxGetScalar(prhs[4]);
    bool wrap = (bool)mxGetScalar(prhs[5]);
    mwSize ndims = mxGetNumberOfDimensions(prhs[0]);
    

    /*mode: 1 euler, 2 align_with_Z, 3 axis, 4 tform*/
    switch ((int)mxGetScalar(prhs[3])) {
        case 1:
            A3D = Euler_rotation3DMatrix(angs[0], angs[1], angs[2]);
            break;
        case 2:
            A3D = alignWithZ(axis);
            break;
        case 3:
            A3D = rotation3DMatrix(angs[0], axis);
            A2D = A3D;
            A2D.window(0, 0, 2, 2);
            break;
        case 4:
            getMatrix2D(prhs[6],A2D);
            A3D = A2D;
            break;
    }
    
    if (ndims == 2)
    {
        Image img, img_out;
        getMatrix2D(prhs[0],img());
        if (XSIZE(A2D) != 0)
        {
            if (gridding)
            {
                KaiserBessel kb;
                produceReverseGriddingMatrix2D(img(),img_out(),kb);
                applyGeometryReverseGridding(img(), A2D, img_out(), kb, IS_NOT_INV, wrap);
                setMatrix2D(img(),plhs[0]);
            }
            else
            {
                applyGeometryBSpline(img_out(), A2D, img(), 3, IS_NOT_INV, wrap);
                setMatrix2D(img_out(),plhs[0]);
            }
        }
        else 
        {
            setMatrix2D(img(),plhs[0]);
        }
    }
    else
    {
        Volume vol, vol_out;
        getMatrix3D(prhs[0],vol());
        if (gridding)
        {
            KaiserBessel kb;
            produceReverseGriddingMatrix3D(vol(),vol_out(),kb);
            applyGeometryReverseGridding(vol(), A3D, vol_out(), kb, IS_NOT_INV, wrap);
            setMatrix3D(vol(),plhs[0]);
        }
        else
        {
            applyGeometryBSpline(vol_out(), A3D, vol(), 3, IS_NOT_INV, wrap);
            setMatrix3D(vol_out(),plhs[0]);
        }
    }
    
}
