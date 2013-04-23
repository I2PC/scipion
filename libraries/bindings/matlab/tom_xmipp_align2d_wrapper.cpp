/*=================================================================
 *
 * tom_xmipp_align2d is a wrapper to xmipp_align2d
 *
 * The calling syntax is:
 *
 *		im_out = tom_xmipp_align2d_wrapper(image,outsize,gridding)
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
#include "polar.h"
#include "filters.h"
#include "tom_xmipp_helpers.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    Image<double> image;
    getMatrix2D(prhs[0],image());
    
    MultidimArray<double> ref;
    getMatrix2D(prhs[1],ref);
    
    int mode = (int) mxGetScalar(prhs[2]);
    
    float max_shift = (float) mxGetScalar(prhs[3]);
    float max_rot = (float) mxGetScalar(prhs[4]);
    
    /*mode: 1: only trans, 2:only rot, 3:complete*/
    double shiftX=0, shiftY=0,rot=0,scale=1;
    bool flip;
    CorrelationAux aux;
    RotationalCorrelationAux auxRot;
    Polar<std::complex<double> > polarFourierImage, polarFourierRef;
    Polar_fftw_plans *plans;
    Matrix2D<double> M(3,3);
    M.initIdentity();
    try {
        switch (mode)
        {
            case 1:
            	bestShift(ref,image(),shiftX,shiftY,aux,NULL,max_shift);
            	M(0,2)=shiftX;
            	M(1,2)=shiftY;
                break;
            case 2:
                normalizedPolarFourierTransform(image(), polarFourierImage, true, XSIZE(image()) / 5, XSIZE(image()) / 2-2, plans, 1);
                normalizedPolarFourierTransform(ref,     polarFourierRef,   true, XSIZE(image()) / 5, XSIZE(image()) / 2-2, plans, 1);
                rot = best_rotation(polarFourierRef,polarFourierImage,auxRot);
                rotation2DMatrix(rot,M,true);
                break;
            case 3:
            	alignImages(ref,image(),M,true);
            	transformationMatrix2Parameters2D(M, flip, scale, shiftX, shiftY, rot);
        }
    }
    catch (XmippError Xe)
    {
        mexErrMsgTxt(Xe.msg.c_str());
    }
    
    /*fetch the results*/
    const char *field_names[] = {"Xoff","Yoff","Psi","Tform"};
    mwSize dims[2] = {1, 1};
    plhs[0] = mxCreateStructArray(2, dims, NUMBER_OF_FIELDS, field_names); 
    
    mxArray *field1 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field1) = shiftX - XSIZE(image())/2+1;
    mxSetField(plhs[0],0,"Xoff",field1);
    
    mxArray *field2 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field2) = shiftY- YSIZE(image())/2+1;
    mxSetField(plhs[0],0,"Yoff",field2);
    
    mxArray *field3 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field3) = rot;
    mxSetField(plhs[0],0,"Psi",field3);

    mxArray *field4;
    setMatrix2D(M,field4);
    mxSetField(plhs[0],0,"Tform",field4);
}
