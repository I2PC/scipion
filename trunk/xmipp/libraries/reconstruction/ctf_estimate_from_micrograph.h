/***************************************************************************
 *
 * Authors:     Javier Angel Velazquez Muriel (javi@cnb.csic.es)
 *              Carlos Oscar Sanchez Sorzano
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

#ifndef _PROG_CTF_ESTIMATE_FROM_MICROGRAPH
#define _PROG_CTF_ESTIMATE_FROM_MICROGRAPH

#include "ctf_estimate_from_psd.h"
#include "ctf_estimate_psd_with_arma.h"

/**@defgroup AssignCTF ctf_estimate_from_micrograph (CTF estimation from a micrograph)
   @ingroup ReconsLibrary
   This program assign different CTFs to the particles in a micrograph */
//@{

/** Assign CTF parameters. */
class ProgCTFEstimateFromMicrograph: public XmippProgram
{
public:
    typedef enum {ARMA, Periodogram} TPSDEstimator_mode;
    typedef enum {OnePerMicrograph, OnePerRegion, OnePerParticle} TPSD_mode;

public:
    /// Parameters for adjust_CTF program
    ProgCTFEstimateFromPSD  prmEstimateCTFFromPSD;
    /// Parameters for ARMA
    ARMA_parameters         ARMA_prm;
    /// Position file
    FileName                fn_pos;
    /// Micrograph filename
    FileName                fn_micrograph;
    /// Output rootname
    FileName                fn_root;
    /// Partition mode
    TPSD_mode               psd_mode;
    /// Dimension of micrograph pieces
    int                     pieceDim;
    /** Overlap among pieces (0=No overlap, 1=Full overlap */
    double                  overlap;
    /** Number of pieces (Nsubpiece x Nsubpiece) for the piece averaging */
    int                     Nsubpiece;
    /** PSDEstimator_mode */
    TPSDEstimator_mode      PSDEstimator_mode;
    /** Bootstrap N */
    int                     bootstrapN;
    /// Estimate a CTF for each PSD
    bool 					estimate_ctf;
public:
    /** Read parameters */
    void readParams();

    /** Define parameters */
    void defineParams();

    /** PSD averaging within a piece.
        Compute the PSD of a piece by subdividing it in smaller pieces and
        averaging their PSDs. The piece will be cut into 3x3 overlapping
        pieces of size N/2 x N/2.*/
    void PSD_piece_by_averaging(MultidimArray<double> &piece,
                                MultidimArray<double> &psd);

    /// Process the whole thing
    void run();
};

/** Fast estimate enhanced PSD.
 *  Set downsampling to 2 for halving the image size. */
void fastEstimateEnhancedPSD(const FileName &fnMicrograph, double downsampling, MultidimArray<double> &enhancedPSD);

//@}
#endif
