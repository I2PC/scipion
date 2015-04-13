/*=================================================================
 *
 * tom_xmipp_normalize is a wrapper to xmipp_normalize
 *
 * The calling syntax is:
 *
 *		image_out = tom_xmipp_normalize_wrapper(image,mask,method)
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
 * created: 24/08/2007
 * by: Andreas Korinek & Carlos Oscar Sorzano
 *
 *=================================================================*/

/*xmipp includes */
#include "tom_xmipp_helpers.h"
#include <data/xmipp_image.h>
#include <data/normalize.h>

#define NONE 0
#define OLDXMIPP 1
#define NEAR_OLDXMIPP 2
#define NEWXMIPP 3
#define MICHAEL 4
#define NEWXMIPP2 5
#define RANDOM 6
#define RAMP 7
#define NEIGHBOUR 8

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{

    Image<double> I;
    getMatrix2D(prhs[0],I());
    
    MultidimArray<int> mask;
    getMatrix2D(prhs[1],mask);
        
    switch((int)mxGetScalar(prhs[2]))
    {
    case OLDXMIPP:
        normalize_OldXmipp(I());
        break;
    case NEAR_OLDXMIPP:
        normalize_Near_OldXmipp(I(), mask);
        break;
    case NEWXMIPP:
        normalize_NewXmipp(I(), mask);
        break;
    case NEWXMIPP2:
        normalize_NewXmipp2(I(), mask);
        break;
    case RAMP:
        normalize_ramp(I(), mask);
        break;
    /*case NEIGHBOUR:
        normalize_remove_neighbours(img, bg_mask, thresh_neigh);
        break;*/
    case MICHAEL:
        normalize_Michael(I(), mask);
        break;
    /*case RANDOM:
        a = rnd_unif(a0, aF);
        b = rnd_unif(b0, bF);

        FOR_ALL_ELEMENTS_IN_MATRIX2D((*I)())
        (*I)(i, j) = a * (*I)(i, j) + b;

        break;*/
    }
    
    setMatrix2D(I(),plhs[0]);
}	

