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
#include "ctf_correct_phase.h"
#include "tom_xmipp_helpers.h"

/*Matlab includes*/
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    Matrix2D<double> I;
    getMatrix2D(prhs[0],I);
    Matrix2D< complex<double> > fft;
    FourierTransform(I, fft);
    
    CorrectPhaseParams prm;
    prm.epsilon = (double) mxGetScalar(prhs[2]);
    prm.method = (int) mxGetScalar(prhs[1]);
    
    prm.ctf.FilterBand = CTF;
    prm.ctf.ctf.enable_CTFnoise = false;

    prm.ctf.ctf.K = (double) mxGetScalar(prhs[3]);
    prm.ctf.ctf.Tm = (double) mxGetScalar(prhs[4]);
    prm.ctf.ctf.kV = (double) mxGetScalar(prhs[5]);
    prm.ctf.ctf.DeltafU = (double) mxGetScalar(prhs[6]);
    prm.ctf.ctf.DeltafV = (double) mxGetScalar(prhs[7]);
    prm.ctf.ctf.azimuthal_angle = (double) mxGetScalar(prhs[8]);    
    prm.ctf.ctf.Cs = (double) mxGetScalar(prhs[9]);
    prm.ctf.ctf.Ca = (double) mxGetScalar(prhs[10]);
    prm.ctf.ctf.espr = (double) mxGetScalar(prhs[11]);
    prm.ctf.ctf.ispr = (double) mxGetScalar(prhs[12]);
    prm.ctf.ctf.alpha = (double) mxGetScalar(prhs[13]);
    prm.ctf.ctf.DeltaF = (double) mxGetScalar(prhs[14]);
    prm.ctf.ctf.DeltaR = (double) mxGetScalar(prhs[15]);
    prm.ctf.ctf.Q0 = (double) mxGetScalar(prhs[16]);
    prm.ctf.ctf.Produce_Side_Info();
    prm.correct(fft);
    InverseFourierTransform(fft, I);
    setMatrix2D(I,plhs[0]);
}	
