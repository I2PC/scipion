/*=================================================================
 *
 * tom_xmipp_volume_segment is a wrapper to xmipp_volume_segment
 *
 * The calling syntax is:
 *
 *		im_out = tom_xmipp_volume_segment_wrapper(volume,sampling_rate,mass,masstype,en_threshold,threshold,wang_radius,probabilistic)
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
#include "tom_xmipp_helpers.h"
#include "volume_segment.h"

#define VOXEL_MASS 1
#define DALTON_MASS 2
#define AA_MASS 3

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    ProgVolumeSegment prm;
    Image<double> vol;
    getMatrix3D(prhs[0],vol());
    prm.V = vol;
    prm.fn_mask = "";
    prm.sampling_rate = (double) mxGetScalar(prhs[1]); 
    double sampling_rate3 = prm.sampling_rate * prm.sampling_rate * prm.sampling_rate;

    double mass = (double) mxGetScalar(prhs[2]); 
    int masstype = (int) mxGetScalar(prhs[3]);  

    switch (masstype)
    {
        case VOXEL_MASS:
            prm.voxel_mass = mass;
            break;
        case DALTON_MASS:
            prm.voxel_mass = mass * 1.207 / sampling_rate3;
            break;
        case AA_MASS:
            prm.voxel_mass = mass * 110 * 1.207 / sampling_rate3;
    }
 
    int en_threshold = (int) mxGetScalar(prhs[4]);
    if (en_threshold)
        prm.en_threshold = true;
    else
        prm.en_threshold = false;
         
    prm.threshold = (double) mxGetScalar(prhs[5]); 
    prm.wang_radius = (int) mxGetScalar(prhs[6]);     
    int do_prob = (int) mxGetScalar(prhs[7]);
    if (do_prob)
        prm.do_prob = true;
    else
        prm.do_prob = false;
    
    Image<double> mask;
    try
    {
        prm.segment(mask);
    }
    catch (XmippError Xe)
    {
       mexErrMsgTxt(Xe.msg.c_str());
    }

    setMatrix3D(mask(),plhs[0]);
    
}
