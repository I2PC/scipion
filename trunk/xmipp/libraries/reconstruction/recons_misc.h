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

#ifndef _RECONS_MISC_HH
#define _RECONS_MISC_HH

#include <data/xmipp_funcs.h>
#include <data/metadata.h>
#include <data/symmetries.h>
#include <data/projection.h>
#include <data/grids.h>
#include "basic_art.h"



/** Build from a Selection File and a Symmetry List.
    The result is stored in the Recons_info array which should point
    to NULL when it is not initialized. */
void buildReconsInfo(MetaData &selfile,
                     const FileName &fn_ctf, const SymList &SL, ReconsInfo * &IMG_Inf,
                     bool do_not_use_symproj);

/** Sort projections orthogonally.
   This function sorts a number of images given by numIMG, whose information
   about their Euler angles are in IMG_inf, into an ordered list which
   gives the indexes. First an image is chosen randomly from the whole
   set. Then all images are compared to the first one, and the most
   perpendicular one is chosen. The remaining set of images are compared
   to this two images and the most perpendicular one to the former two
   is chosen, and so on until no image is left in the set.

   If the result in ordered list is 4, 70, 54, 203, 1, 0, ... it means
   that the first image is the number 4, then goes the 70, then the 54, ...

   If N!=-1 then the product is done only with the last N images. A very
   useful value is N=2*/
void sortPerpendicular(int numIMG, ReconsInfo *IMG_Inf,
                       MultidimArray<int> &ordered_list, int N = 2);

/** No projection sorting at all.
    This function directly returns the same order as in the selection file */
void noSort(int numIMG, MultidimArray<int> &ordered_list);

/** Randomize the projections.
   This function sorts randomly a number of images given by numIMG. */
void sortRandomly(int numIMG, MultidimArray<int> &ordered_list);

/** Update residual vector for WLS ART */
void updateResidualVector(BasicARTParameters &prm, GridVolume &vol_basis,
                          double &kappa, double &pow_residual_vol, double &pow_residual_imgs);

/* ------------------------------------------------------------------------- */
/** Variability structure */
class VariabilityClass
{
public:
    typedef enum {VAR_none, VAR_measuring, VAR_analyzing} t_VAR_status;
    t_VAR_status VAR_state;
    int Zoutput_volume_size;
    int Youtput_volume_size;
    int Xoutput_volume_size;
    BasicARTParameters *prm;

    /// Vector of training vectors
    std::vector < MultidimArray<double> > VA;

    /// Number of updates so far
    int N;

    /// Constructor
    VariabilityClass(BasicARTParameters *_prm,
                     int _Zoutput_volume_size, int _Youtput_volume_size,
                     int _Xoutput_volume_size);

    /** Start a new ART iteration. */
    void newIteration();

    /** Update data with a new volume.
        The update volume is set to zeros after this function */
    void newUpdateVolume(GridVolume *ptr_vol_out, Projection &read_proj);

    /** Finish analysis. */
    void finishAnalysis();
};

/* ------------------------------------------------------------------------- */
/** POCS structure */
class POCSClass
{
public:
    typedef enum {POCS_measuring, POCS_use, POCS_lowering, POCS_N_measure,
                  POCS_N_use} t_POCS_status;
    t_POCS_status POCS_state;
    double POCS_avg;
    double POCS_stddev;
    double POCS_min;
    double POCS_max;
    double POCS_mean_error;
    double POCS_max_error;
    double POCS_global_mean_error;
    int POCS_freq;
    int POCS_i;
    int POCS_vec_i;
    int POCS_used;
    int POCS_N;
    int Zoutput_volume_size;
    int Youtput_volume_size;
    int Xoutput_volume_size;
    bool apply_POCS;
    MultidimArray<double> POCS_errors;
    BasicARTParameters *prm;

    /// Constructor
    POCSClass(BasicARTParameters *_prm,
              int _Zoutput_volume_size, int _Youtput_volume_size,
              int _Xoutput_volume_size);

    /// Start New ART iteration
    void newIteration();

    /// Start new Projection
    void newProjection();

    /// Apply
    void apply(GridVolume &vol_basis, int it, int images);
};


#endif
