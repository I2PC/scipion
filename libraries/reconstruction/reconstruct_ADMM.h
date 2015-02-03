/***************************************************************************
 *
 * Authors:        Masih Nilchian (masih.nilchian@epfl.ch)
 *                 Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
#ifndef _PROG_RECONSTRUCT_ADMM_HH
#define _PROG_RECONSTRUCT_ADMM_HH

#include <data/xmipp_program.h>
#include <data/matrix1d.h>
#include <data/xmipp_fftw.h>
#include <data/ctf.h>
#include <data/mask.h>
#include <data/symmetries.h>

/**@defgroup ReconstructADMMProgram Reconstruct Alternating Direction Method of Multipliers
   @ingroup ReconsLibrary */
//@{

class AdmmKernel
{
public:
	double supp;
	double projectionStep;
	double autocorrStep;
	double alpha;
	Matrix1D<double> projectionProfile; // Look-up table
	MultidimArray< std::complex<double> > FourierProjectionAutocorr;
	FourierTransformer transformer;
	MultidimArray<double> projectionAutocorrWithCTF;

	void initializeKernel(double alpha, double a, double _projectionStep);

	double projectionValueAt(double u, double v);

	void convolveKernelWithItself(double _autocorrStep);

	void applyCTFToKernelAutocorrelation(CTFDescription &ctf, double Ts, MultidimArray<double> &autocorrelationWithCTF);

	void getKernelAutocorrelation(MultidimArray<double> &autocorrelation);

	/* Direction=x,y,z, mode=L (Lx), T (L^T x), 2 (L^TLx) */
	void computeGradient(MultidimArray<double> &gradient, char direction, bool adjoint=false);

	/** Compute the kernel in 3D space */
	void computeKernel3D(MultidimArray<double> &kernel);
};

class ProgReconsADMM: public XmippProgram
{
public:
	FileName fnIn; // Set of input images
	FileName fnRoot; // Rootname for output files
	FileName fnFirst; // First reconstruction filename
	FileName fnHtb; // Htb filename
	FileName fnHtKH; // HtKH filename
	String kernelShape;
	double a;     // Kaiser-Bessel maximum radius
	double alpha; // Kaiser-Bessel shape
	double mu; // Augmented Lagrange
	double lambda; // Total variation regularization
	double lambda1; // Tikhonov regularization

	double Ti; // Downsampling factor for the volume
	double Tp; // Downsampling factor for the projections
	bool useWeights; // Use weights if available in the input metadata
	bool useCTF; // Use CTF if available in the input metadata
	int Ncgiter; // Number of conjugate gradient iterations
	int Nadmmiter; // Number of ADMM iterations
	bool positivity; // Positivity constraint
	bool applyMask;
	double Ts; // Sampling rate
	Mask mask; // Mask
	String symmetry;
	bool saveIntermediate;
	size_t Nprocs;
	size_t rank;
public:
	ProgReconsADMM();
    void defineParams();
    void readParams();
    void show();
    void run();

    void produceSideInfo();

    /** Project the volume V onto P using rot,tilt,psi */
    void project(double rot, double tilt, double psi, MultidimArray<double> &P, bool adjoint=false, double weight=1.);

    /** Project the volume V onto P using r1 and r2 as the coordinate system */
    void project(const Matrix1D<double> &r1, const Matrix1D<double> &r2, MultidimArray<double> &P, bool adjoint=false, double weight=1.);

    /** H^t*b
     */
    void constructHtb();

    /** Compute H^t*K*H.
     * H is the projection operator, K is the CTF operator
     */
    void computeHtKH(MultidimArray<double> &kernelV);

    /** Add regularization to the kernel */
    void addRegularizationTerms();

    /** Apply kernel 3D */
    void applyKernel3D(MultidimArray<double> &x, MultidimArray<double> &AtAx);

    /** Apply Lt filter */
    void applyLFilter(MultidimArray< std::complex<double> > &fourierL, bool adjoint=false);
    void applyLtFilter(MultidimArray< std::complex<double> > &fourierL, MultidimArray<double> &u, MultidimArray<double> &d);

    /** Apply conjugate gradient */
    void applyConjugateGradient();

    /** POCS projection */
    void doPOCSProjection();

    /** Update u and d */
    void updateUD();

    /** Convert from coefficients to voxel volume */
    void produceVolume();

    // Share a volume among nodes
    virtual void shareVolume(MultidimArray<double> &V) {}

	// Redefine how to synchronize
	virtual void synchronize() {}

public:
	AdmmKernel           kernel;
	Image<double>        CHtb; // First reconstructed volume
	Image<double>        Ck, Vk; // Reconstructed volume
	MetaData             mdIn; // Set of images and angles
	MultidimArray<std::complex<double> > fourierKernelV;
	MultidimArray<double> paddedx;
	FourierTransformer    transformerPaddedx, transformerL;

	MultidimArray<double> ux, uy, uz, dx, dy, dz, ud;
	MultidimArray< std::complex<double> > fourierLx, fourierLy, fourierLz;
	SymList SL;
};

//@}
#endif
