/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2004)
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

#ifndef _ML_REFINE3D_H
#define _ML_REFINE3D_H

#include <data/fft.h>
#include <data/args.h>
#include <data/funcs.h>
#include <data/metadata.h>
#include <data/image.h>
#include <data/filters.h>
#include <data/mask.h>
#include <data/morphology.h>
#include <data/grids.h>
#include <data/blobs.h>
#include <data/symmetries.h>
#include <data/sampling.h>
#include <data/projection.h>
#include "directions.h"
#include "reconstruct_art.h"
// FIXME:
//#include "reconstruct_fourier.h"
#include "ml_align2d.h"
#include "mlf_align2d.h"
#include "symmetrize.h"
#include "volume_segment.h"

#include <vector>

/**@defgroup Refine3d ml_refine3d (Maximum likelihood 3D refinement)
   @ingroup ReconsLibrary */
//@{
/** Refine3d parameters. */
class ProgRefine3D: public XmippProgram
{

public:
    // Filename for reference volume, symmetry file and output rootname
    FileName fn_sel, fn_vol, fn_sym, fn_root, fn_solv, fn_iter, fn_symmask;
    // Selfile with reference volumes
    MetaData SFvol;
    // Number of volumes to refine
    int Nvols;
    // Iteration numbers
    int istart, Niter;
    // Convergence check
    double eps;
    // Angular sampling interval (degree)
    double angular;
    /// File handler for the history file
    std::ofstream fh_hist;
    // Use fourier-interpolation instead of WLS-ART for reconstruction in ML
    bool reconstruct_fourier;
    // Low-pass filter digital frequency
    double lowpass;
    // For user-provided tilt range
    double tilt_range0, tilt_rangeF;
    // Do not use a starting volume in wlsART reconstruction
    bool wlsart_no_start;
    // Threshold for flooding-like solvent mask
    double threshold_solvent;
    // Do Wang/Terwilliger-like probabilistic solvent flattening
    bool do_prob_solvent;
    // Do separate_object-like selection of largest connected volume in
    // binarized solvent mask
    bool do_deblob_solvent;
    // Dilate (binarized) solvent mask?
    int dilate_solvent;
    // Use MLF mode
    bool fourier_mode;
    // Flag to skip reconstruction
    bool skip_reconstruction;
    // Perturb angles of reference projections
    bool do_perturb;
    // sampling object
    Sampling mysampling;
    // Symmetry setup
    int symmetry, sym_order;
    // Number of reference projections per 3D model
    int nr_projections;

    // A pointer to the 2D alignment and classification program
    ProgML2D * ml2d;

private:
    ///Helper function to show to screen and file history
    void showToStream(std::ostream &out);

public:

    /// Empty constructor, call the constructor of ProgML2D
    /// with the ML3D flag set to true
    ProgRefine3D(bool fourier = false);
    /** Destructor */
    ~ProgRefine3D();

    /// Read additional arguments for 3D-process from command line
    void readParams();

    /// Define the parameters accepted
    void defineParams();

    /// MLF Usage
    void MLF_usage();

    /// Show
    void show();

    //Call produceSideInfo of ML2D and
    // Fill sampling and create DFlib
    virtual void produceSideInfo();

    /// Project the reference volume in evenly sampled directions
    void projectReferenceVolume(MetaData &SFlib, int rank = 0, int size = 1) ;

    /// (For mpi-version only:) calculate noise averages and write to disc
    void makeNoiseImages(std::vector<Image<double>  > &Iref) ;

    /// reconstruction by (weighted ART) or Fourier interpolation
    void reconstruction(int argc, char **argv,
                        int iter, int volno, int noise = 0);

    /// Calculate 3D SSNR according to Unser ea. (2005)
    void calculate3DSSNR(MultidimArray<double> &spectral_signal, int iter);

    /// After reconstruction update reference volume selfile
    void remakeSFvol(int iter, bool rewrite = false, bool include_noise = false) ;

    /// Merge MLalign2D classification selfiles into volume classes
    void concatenateSelfiles(int iter);

    /// Maksing, filtering etc. of the volume
    void postProcessVolumes(int argc, char **argv) ;

    /// Convergency check
    bool checkConvergence(int iter) ;



    ///Provides implementation of the run function
    void run();
};
//@}
#endif
