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
#include "tom_xmipp_helpers.h"
#include <data/xmipp_image.h>

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
    
    bool wrap = (bool)mxGetScalar(prhs[5]);
    mwSize ndims = mxGetNumberOfDimensions(prhs[0]);
    

    /*mode: 1 euler, 2 align_with_Z, 3 axis, 4 tform*/
    switch ((int)mxGetScalar(prhs[3])) {
        case 1:
            Euler_angles2matrix(angs[0], angs[1], angs[2], A3D, true);
            break;
        case 2:
            alignWithZ(axis,A3D);
            break;
        case 3:
            rotation3DMatrix(angs[0], axis, A3D);
            A2D.resize(3,3);
            A2D(0,0)=A3D(0,0); A2D(0,1)=A3D(0,1); A2D(0,2)=A3D(0,2);
            A2D(1,0)=A3D(1,0); A2D(1,1)=A3D(1,1); A2D(1,2)=A3D(1,2);
            A2D(2,0)=A3D(2,0); A2D(2,1)=A3D(2,1); A2D(2,2)=A3D(2,2);
            break;
        case 4:
            getMatrix2D(prhs[6],A2D);
            A3D = A2D;
            break;
    }
    
    if (ndims == 2)
    {
        Image<double> img, img_out;
        getMatrix2D(prhs[0],img());
        if (MAT_XSIZE(A2D) != 0)
        {
			applyGeometry(BSPLINE3,img_out(), img(), A2D, IS_NOT_INV, wrap);
			setMatrix2D(img_out(),plhs[0]);
        }
        else 
        {
            setMatrix2D(img(),plhs[0]);
        }
    }
    else
    {
        Image<double> vol, vol_out;
        getMatrix3D(prhs[0],vol());
		applyGeometry(BSPLINE3, vol_out(), vol(), A3D, IS_NOT_INV, wrap);
		setMatrix3D(vol_out(),plhs[0]);
    }
}
