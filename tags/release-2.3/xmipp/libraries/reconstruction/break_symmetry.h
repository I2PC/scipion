/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.uam.es (2004)
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

#include <data/args.h>
#include <data/funcs.h>
#include <data/selfile.h>
#include <data/docfile.h>
#include <data/image.h>
#include <data/volume.h>
#include <data/filters.h>
#include <data/mask.h>

#include "symmetries.h"
#include "projection.h"
#include "directions.h"
#include "reconstruct_wbp.h"
#include "symmetrize.h"

#include <vector>

/**@defgroup BreakSym break_symmetry (Break symmetry of a set of reference volumes)
   @ingroup ReconsLibraryPrograms */
//@{
/** Break_Sym parameters. */
class Prog_Break_Sym_prm
{

public:
    // Various filenames
    FileName fn_sel, fn_vol, fn_sym, fn_root, fn_mask, fn_iter;
    // Selfile with experimental images or reference volumes
    SelFile SF, SFvol;
    // Symmetry list
    SymList SL;
    // Verbosity flag
    int verb;
    // dimension
    int dim;
    // Number of volumes to process
    int Nvols;
    // Iteration numbering
    int Niter, istart;
    // Convergence check
    double eps;
    // Radius for masking of volume
    double mask_radius;
    // File handler for the history file
    std::ofstream fh_hist;
    // Reference volumes
    std::vector<Matrix3D<double> > vols;
    // Mask
    VolumeXmipp mask;

public:

    /// Read additional arguments for 3D-process from command line
    void read(int argc, char **argv);

    /// Usage
    void usage();

    /// Show
    void show();

    /// Project the reference volume in evenly sampled directions
    void process_one_image(ImageXmipp &img, int &opt_vol, int &opt_sym, double &maxcorr);

    /// Process various images
    void process_selfile(SelFile &SF, std::vector<SelFile> &SFout, double &avecorr);

    /// reconstruction by (weighted ART) or WBP
    void reconstruction(int argc, char **argv, SelFile &SF, int iter, int volno);

};
//@}
