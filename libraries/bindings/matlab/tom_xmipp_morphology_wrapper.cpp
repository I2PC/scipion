/*=================================================================
 *
 * tom_xmipp_morphology is a wrapper to xmipp_morphology
 *
 * The calling syntax is:
 *
 *		im_out = tom_xmipp_morphology_wrapper(image,operation,size,count,neig)
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
#include <data/morphology.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    #define DILATION 1
    #define EROSION  2
    #define OPENING  3
    #define CLOSING  4
     
    int operation = (int) mxGetScalar(prhs[1]);
    int size = (int) mxGetScalar(prhs[2]);
    int count = (int) mxGetScalar(prhs[3]);
    int neig = (int) mxGetScalar(prhs[4]);
    
    mwSize ndims = mxGetNumberOfDimensions(prhs[0]);
    
    if (ndims == 2)
    {
        Image<double> img, img_out;
        getMatrix2D(prhs[0],img());
        img_out() = img();

        switch (operation)
        {
        case DILATION:
            dilate2D(img(),img_out(),neig,count,size);
            break;
        case EROSION:
            erode2D(img(),img_out(),neig,count,size);
            break;
        case OPENING:
            opening2D(img(),img_out(),neig,count,size);
            break;
        case CLOSING:
            closing2D(img(),img_out(),neig,count,size);
            break;
        }
        
        setMatrix2D(img_out(),plhs[0]);     
    }
    else
    {
        Image<double> vol, vol_out;
        getMatrix3D(prhs[0],vol());
        vol_out = vol();
        
        switch (operation)
        {
        case DILATION:
            dilate3D(vol(),vol_out(),neig,count,size);
            break;
        case EROSION:
            erode3D(vol(),vol_out(),neig,count,size);
            break;
        case OPENING:
            opening3D(vol(),vol_out(),neig,count,size);
            break;
        case CLOSING:
            closing3D(vol(),vol_out(),neig,count,size);
            break;
        }
        
        setMatrix3D(vol_out(),plhs[0]);     
    }
    
}
