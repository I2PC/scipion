/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#ifndef _FOURIER_FILTER_HH
#define _FOURIER_FILTER_HH

#include <data/ctf.h>
#include <data/xmipp_fftw.h>
#include <data/filters.h>

/**@defgroup FourierMasks Masks in Fourier space
   @ingroup ReconsLibrary */
//@{
/** Filter class for Fourier space.

   Example of use for highpass filtering
   @code
      Image<double> I;
      I.read("image.xmp");
      FourierFilter Filter;
      Filter.FilterBand=HIGHPASS;
      Filter.w1=w_cutoff;
      Filter.raised_w=slope;
      I().setXmippOrigin();
      Filter.applyMaskSpace(I());
      I.write("filtered_image.xmp");
   @endcode
   
   Example of use for wedge filtering
   @code
        Image<double> V;
        V.read("1rux64.vol");
        V().setXmippOrigin();
        FourierFilter Filter;
        Filter.FilterBand=WEDGE;
        Filter.FilterShape=WEDGE;
        Filter.t1=-60;
        Filter.t2= 60;
        Filter.rot=Filter.tilt=Filter.psi=0;
        Filter.do_generate_3dmask=true;
        Filter.generateMask(V());
        Filter.applyMaskSpace(V());
   @endcode

   For volumes you the mask is computed on the fly and
   in this way memory is saved (unless do_generate_3dmask == true).
*/
class FourierFilter: public XmippFilter
{
public:
#define RAISED_COSINE 1
    /** Shape of the decay in the filter.
       Valid types are RAISED_COSINE. */
    int FilterShape;

#define LOWPASS       1
#define HIGHPASS      2
#define BANDPASS      3
#define STOPBAND      4
#define CTF           5
#define WEDGE         7
#define GAUSSIAN      8
#define CONE          9
#define CTFPOS       10
#define BFACTOR      11
#define REALGAUSSIAN 12
#define SPARSIFY     13
    /** Pass band. LOWPASS, HIGHPASS, BANDPASS, STOPBAND, CTF, CTFPOS,
       WEDGE, CONE, GAUSSIAN, FROM_FILE, REALGAUSSIAN, BFACTOR, SPARSIFY */
    int FilterBand;

    /** Cut frequency for Low and High pass filters, first freq for bandpass.
        Normalized to 1/2*/
    double w1;

    /** Second frequency for bandpass and stopband. Normalized to 1/2 */
    double w2;

    /** Wedge and cone filter parameters */
    double t1, t2,rot,tilt,psi;

    /** Percentage of coefficients to throw */
    double percentage;

    /** Filename in which store the mask (valid only for fourier masks) */
    FileName maskFn;

    /** Pixels around the central frequency for the raised cosine */
    double raised_w;

    /** CTF parameters. */
    CTFDescription ctf;
    
    /** Correct phase before applying CTF */
    bool do_correct_phase;

    /** Flag to generate 3D mask */
    bool do_generate_3dmask;

public:
    /** Define parameters */
    static void defineParams(XmippProgram * program);

    /** Read parameters from command line.
        If a CTF description file is provided it is read. */
    void readParams(XmippProgram * program);

    /** Process one image */
    void apply(MultidimArray<double> &img);

    /** Empty constructor */
    FourierFilter();

    /** Clear */
    void init();

    /** Show. */
    void show();

    /** Compute the mask value at a given frequency.
        The frequency must be normalized so that the maximum frequency
        in each direction is 0.5 */
    double maskValue(const Matrix1D<double> &w);

    /** Generate nD mask. */
    void generateMask(MultidimArray<double> &v);

    /** Apply mask in real space. */
    void applyMaskSpace(MultidimArray<double> &v);

    /** Apply mask in Fourier space.
     * The image remains in Fourier space.
     */
    void applyMaskFourierSpace(const MultidimArray<double> &v, MultidimArray<std::complex<double> > &V);

    /** Get the power of the nD mask. */
    double maskPower();
    
    /** Correct phase */
    void correctPhase();
public:
    // Auxiliary vector for representing frequency values
    Matrix1D<double> w;

    // Auxiliary mask for the filter in 3D
    MultidimArray<int> maskFourier;

    // Auxiliary mask for the filter in 3D
    MultidimArray<double> maskFourierd;

    // Transformer
    FourierTransformer transformer;

    // Auxiliary variables for sparsify
    MultidimArray<double> vMag, vMagSorted;
};

/** Fast access to bandpass filter.
 * Frequencies are normalized to 0.5 */
void bandpassFilter(MultidimArray<double> &img, double w1, double w2, double raised_w);

/** Fast access to Gaussian filter.
 * Frequencies are normalized to 0.5 */
void gaussianFilter(MultidimArray<double> &img, double w1);

/** Fast access to real gaussian filter.
 * Sigma is in pixel units.
 */
void realGaussianFilter(MultidimArray<double> &img, double sigma);
//@}
#endif
