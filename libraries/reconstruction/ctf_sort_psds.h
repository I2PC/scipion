/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Slavica Jonic (slavica.jonic@impmc.jussieu.fr)
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
#ifndef _PROG_CTF_SORT_PSDS_HH
#  define _PROG_CTF_SORT_PSDS_HH

#include <data/xmipp_funcs.h>
#include <data/metadata.h>
#include <data/multidim_array.h>
#include <data/xmipp_program.h>

/**@defgroup SortPSD psd_Sort (Sort visualization of the PSD)
   @ingroup ReconsLibrary */
//@{

/** Results of an evaluation of a PSD */
class PSDEvaluation {
public:
	double defocusU;
	double defocusV;
	double PSDcorrelation90;
	double firstZeroRatio;
	double firstZeroAvg;
	double maxFreq;
	double firstZeroDisagreement;
	double beating;
	double maxDampingAtBorder;
	double PSDradialIntegral;
	double fittingScore;
	double fittingCorr13;
	double PSDVariance;
	double PSDPC1Variance;
	double PSDPCRunsTest;
	double histogramNormality;
};

/* Sort PSD Program Parameters ------------------------------------------ */
/** Parameter class for the project program */
class ProgPSDSort: public XmippMetadataProgram
{
public:
    /// Bandpass filter low frequency (in Fourier space, max 0.5)
    double filter_w1;

    /// Bandpass filter high frequency (in Fourier space, max 0.5)
    double filter_w2;

    /// Decay width (raised cosine)
    double decay_width;

    /// Lower frequency for the mask (in Fourier space, max 0.5)
    double mask_w1;

    /// Higher frequency for the mask (in Fourier space, max 0.5)
    double mask_w2;

public:
    /** Empty constructor */
    ProgPSDSort();

    /** Use micrographs for iterating in the metadata */
    void defineLabelParam();

    /** Read from a command line. */
    void readParams();

    /** Define parameters. */
    void defineParams();

    /** Show parameters. */
    void show();

    /** Process micrograph */
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);
};
//@}
#endif
