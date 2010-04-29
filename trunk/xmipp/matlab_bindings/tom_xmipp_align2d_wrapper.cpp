/*=================================================================
 *
 * tom_xmipp_align2d is a wrapper to xmipp_align2d
 *
 * The calling syntax is:
 *
 *      im_out = tom_xmipp_align2d_wrapper(image,outsize,gridding)
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
#include "tom_xmipp_helpers.h"
#include "reconstruction/align2d.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{

    Prog_align2d_prm alignParams;

    ImageXmipp image;
    getMatrix2D(prhs[0],image());

    Matrix2D<double> ref;
    getMatrix2D(prhs[1],ref);

    int mode = (int) mxGetScalar(prhs[2]);

    float max_shift = (float) mxGetScalar(prhs[3]);
    float max_rot = (float) mxGetScalar(prhs[4]);
    float psi_interval = (float) mxGetScalar(prhs[5]);
    float Rin = (float) mxGetScalar(prhs[6]);
    float Rout = (float) mxGetScalar(prhs[7]);
    double outside = (double) mxGetScalar(prhs[8]);

    bool success;
    /*mode: 1: only trans, 2:only rot, 3:complete*/
    try
    {
        switch (mode)
        {
        case 1:
            success = alignParams.align_trans(image,ref,max_shift,outside);
            break;
        case 2:
            success = alignParams.align_rot(image,ref,max_rot,Rin,Rout,outside);
            break;
        case 3:
            success = alignParams.align_complete_search(image,ref,max_shift,max_rot,psi_interval,Rin,Rout,outside);
            break;
        }
    }
    catch (Xmipp_error Xe)
    {
        mexErrMsgTxt(Xe.msg.c_str());
    }

    /*fetch the results*/
    const char *field_names[] = {"Xoff","Yoff","Psi","Tform","success"};
    mwSize dims[2] = {1, 1};
    plhs[0] = mxCreateStructArray(2, dims, NUMBER_OF_FIELDS, field_names);

    mxArray *field1 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field1) = image.Xoff() - XSIZE(image())/2+1;
    mxSetField(plhs[0],0,"Xoff",field1);

    mxArray *field2 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field2) = image.Yoff() - YSIZE(image())/2+1;
    mxSetField(plhs[0],0,"Yoff",field2);

    mxArray *field3 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field3) = image.Psi();
    mxSetField(plhs[0],0,"Psi",field3);

    mxArray *field4;
    setMatrix2D(image.get_transformation_matrix(),field4);
    mxSetField(plhs[0],0,"Tform",field4);

    mxArray *field5 = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(field5) = (double) success;
    mxSetField(plhs[0],0,"success",field5);

}
