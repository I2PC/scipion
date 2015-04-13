/***************************************************************************
 *
 * Author:     Sjors H.W. Scheres               (scheres@cnb.csic.es)
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

#ifndef VOLUME_CORRECT_BFACTOR_H
#define VOLUME_CORRECT_BFACTOR_H

#include <data/geometry.h>
#include <data/xmipp_fft.h>
#include <data/xmipp_fftw.h>
#include <data/xmipp_image.h>
#include <data/xmipp_program.h>

/** @defgroup Correct Bfactor
  * @ingroup ReconsLibrary
  */

/** correct_bfactor parameters. */
class ProgVolumeCorrectBfactor: public XmippMetadataProgram
{
protected:

    /** Low and high resolution limits for automated fit */
    double fit_minres, fit_maxres, apply_maxres;
    /** Pixels size in Angstrom */
    double sampling_rate;

    /** Mode for B-factor correction
     * Three modes are provided: 
     *  1. automated fit according to Rosenthal and Henderson (2003)
     *  2. fit according to a reference map
     *  3. ad-hoc correction with a user-defined B-factor
     */
#define BFACTOR_AUTO 1
#define BFACTOR_REF 2
#define BFACTOR_ADHOC 3
#define ALLPOINTS_REF 4

    int mode;

    /** X-size of the original volume or image */
    int xsize;

    /** Reference map file name (for mode BFACTOR_REF) */
    FileName fn_ref;

    /** Ad hoc B-factor (for mode BFACTOR_ADHOC) */
    double adhocB;

    /** Filename for FSC curve */
    FileName fn_fsc;

    /** Functions from XmippMetadataProgram */

    void show();
    void defineParams();
    void readParams();
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

public:

    /** Constructor */
    ProgVolumeCorrectBfactor();

    /** Destructor */
    ~ProgVolumeCorrectBfactor()
    { }
    ;

    /** Make Guinier plot (logarithm of spherically averaged amplitudes versus resolution in 1/A^2)
     * @ingroup Correct Bfactor
     */
    void make_guinier_plot(MultidimArray< std::complex< double > > &m1,
                           std::vector<fit_point2D> &guinier);

    /** Apply B-factor
     * @ingroup Correct Bfactor
     * This will apply the following operation to the FT1: exp (- bfactor / 4 * d *d)
     */
    void apply_bfactor(MultidimArray< std::complex< double > > &FT1,
                       double bfactor);

    /** Apply B-factor
     * @ingroup Correct Bfactor
     * This will adjust the power spectrum according to the difference Guinier
     */
    void apply_allpoints(MultidimArray< std::complex< double > > &FT1,
                         std::vector<fit_point2D> &guinier_diff);

    /** Read FSC file and convert to signal-to-noise weights
     * @ingroup Correct Bfactor
     */
    void get_snr_weights(std::vector<double> &snr);

    /** Apply SNR weights to Fourier Transform of a volume
     * @ingroup Correct Bfactor
     */
    void apply_snr_weights(MultidimArray< std::complex< double > > &FT1,
                           std::vector<double> &snr);

    /** Write Guinier plot in a textfile
     * @ingroup Correct Bfactor
     */
    void write_guinierfile(const FileName& fn_guinier,
                           std::vector<fit_point2D> &guinierin,
                           std::vector<fit_point2D> &guinierweighted,
                           std::vector<fit_point2D> &guiniernew,
                           double intercept,
                           std::vector<fit_point2D> &guinierref);

    /** B-factor correction (sharpening)
     * @ingroup FourierOperations
     */
    void bfactor_correction(MultidimArray< double > &m1, const FileName &fn_guinier);
};
//@}
#endif
