/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
 *
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

#ifndef SCORE_MICROGRAPH_H_
#define SCORE_MICROGRAPH_H_

#include <data/xmipp_program.h>
#include "ctf_estimate_from_micrograph.h"
#include "ctf_estimate_from_psd.h"
#include "ctf_sort_psds.h"



/**@defgroup ScoreMicrograph score_micrograph (Evaluates the score of a micrograph)
   @ingroup ReconsLibrary */
//@{
/* Score Micrograph Program Parameters ------------------------------------------ */
/** Parameter class for the project program */
class ProgScoreMicrograph: public XmippProgram
{
public:

	//First we assume that the particle size is the same for the X and Y axis
	int particleSize;

	ProgCTFEstimateFromMicrograph prmEstimateCTFFromMicrograph;

    /// Parameters for adjust_CTF program
    ProgCTFEstimateFromPSD  prmEstimateCTFFromPSD;

    /// Micrograph filename
    FileName                fn_micrograph;
    /// Dimension of micrograph pieces
    int                     pieceDim;
    /** Overlap among pieces (0=No overlap, 1=Full overlap */
    double                  overlap;
    /** Skip borders.
     * The number of pieces around the border that must be skipped. */
    int                     skipBorders;
    /** Number of pieces (Nsubpiece x Nsubpiece) for the piece averaging */
    int                     Nsubpiece;
    /** Bootstrap N */
    int                     bootstrapN;
    /// Estimate a CTF for each PSD
    bool 					estimate_ctf;
    /// Defocus range
    double               defocus_range;

    ProgCTFEstimateFromMicrograph::TPSD_mode PSDEstimator_mode;

    ProgCTFEstimateFromMicrograph::TPSD_mode psd_mode;

    ProgPSDSort prmPSDSort;



public:
    /** Read from a command line.*/
    void readParams();

    /** Define parameters. */
    void defineParams();

    void run();
    /** Process one image */
    //void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

};

//@}


#endif /* SCORE_MICROGRAPH_H_ */
