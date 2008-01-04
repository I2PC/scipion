/***************************************************************************
 *
 * Authors:     Javier Angel Velazquez Muriel (javi@cnb.uam.es)
 *              Carlos Oscar Sánchez Sorzano
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#ifndef _PROG_ASSIGN_CTF
#define _PROG_ASSIGN_CTF

#include "adjust_ctf.h"
#include "sparma.h"

/**@defgroup AssignCTF ctf_estimate_from_micrograph (CTF estimation from a micrograph)
   @ingroup ReconsLibraryPrograms
   This program assign different CTFs to the particles in a micrograph */
//@{

/** Assign CTF parameters. */
class Prog_assign_CTF_prm
{
public:
    typedef enum {ARMA, Periodogram} TPSD_mode;

public:
    /// Parameters for adjust_CTF program
    Adjust_CTF_Parameters   adjust_CTF_prm;
    /// Parameters for ARMA
    ARMA_parameters         ARMA_prm;
    /// Reversed endian
    bool                  reversed;
    /// X dimension of micrograph pieces
    int                     N_horizontal;
    /// Y dimension of micrograph pieces
    int                     N_vertical;
    /// Selfile with all projections
    FileName                selfile_fn;
    /// Position file
    FileName                picked_fn;
    /// PSD files root
    FileName                PSDfn_root;
    /// Micrograph filename
    FileName                image_fn;
    /** the center of the windows in which the CTF is computed
        are the particles (stored at the .pos file) instead of
        a regular grid. By default this is false.
    */
    bool                    compute_at_particle;
    /** Perform averaging. */
    bool                    micrograph_averaging;
    /** Perform averaging. */
    bool                    piece_averaging;
    /** Number of pieces (Nside_piece x Nside_piece) for the piece averaging */
    int                     Nside_piece;
    /** PSD_mode */
    TPSD_mode               PSD_mode;
    /** Don't adjust CTF.
        That is, only compute the PSD */
    bool                    dont_adjust_CTF;
public:
    /** Selfile mode.
        If image_fn finishes in .sel, then the selfile_mode is on,
        otherwise is off */
    bool                    selfile_mode;

public:
    /** Read parameters from file.
        If do_not_read_files is TRUE then all FileNames parameters
        are not read */
    void read(const FileName &fn_prm, bool do_not_read_files = false);

    /** Write parameters to file.
        The directory is an option used in the grid. */
    void write(const FileName &fn_prm, const std::string & directory = "");

    /** PSD averaging within a piece.
        Compute the PSD of a piece by subdividing it in smaller pieces and
        averaging their PSDs. The piece will be cut into 3x3 overlapping
        pieces of size N/2 x N/2.*/
    void PSD_piece_by_averaging(Matrix2D<double> &piece,
                                Matrix2D<double> &psd);

    /// Process the whole thing
    void process();
};
//@}
#endif
