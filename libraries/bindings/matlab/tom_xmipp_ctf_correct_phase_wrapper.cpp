/*=================================================================
 *
 * tom_xmipp_ctf_correct_phase is a wrapper to xmipp_ctf_correct_phase
 *
 * The calling syntax is:
 *
 *		im_corr = tom_xmipp_ctf_correct_phase(image,ctfmodel,method,epsilon)
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
 * created: 22/10/2007
 * by: Andreas Korinek
 *
 *=================================================================*/

/*xmipp includes */
#include "ctf_phase_flip.h"
#include "tom_xmipp_helpers.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    MultidimArray<double> I;
    getMatrix2D(prhs[0],I);
    
    ProgCTFPhaseFlipping prm;
    CTFDescription ctf;

    ctf.K = (double) mxGetScalar(prhs[3]);
    ctf.Tm = (double) mxGetScalar(prhs[4]);
    ctf.kV = (double) mxGetScalar(prhs[5]);
    ctf.DeltafU = (double) mxGetScalar(prhs[6]);
    ctf.DeltafV = (double) mxGetScalar(prhs[7]);
    ctf.azimuthal_angle = (double) mxGetScalar(prhs[8]);
    ctf.Cs = (double) mxGetScalar(prhs[9]);
    ctf.Ca = (double) mxGetScalar(prhs[10]);
    ctf.espr = (double) mxGetScalar(prhs[11]);
    ctf.ispr = (double) mxGetScalar(prhs[12]);
    ctf.alpha = (double) mxGetScalar(prhs[13]);
    ctf.DeltaF = (double) mxGetScalar(prhs[14]);
    ctf.DeltaR = (double) mxGetScalar(prhs[15]);
    ctf.Q0 = (double) mxGetScalar(prhs[16]);
    
    try {
    	actualPhaseFlip(I,ctf);
    }
    catch (XmippError Xe)
    {
        mexErrMsgTxt(Xe.msg.c_str());
    }

    setMatrix2D(I,plhs[0]);
}	
