/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
#ifndef _PROG_ART_HH
#  define _PROG_ART_HH

#include <iostream>

#include "basic_art.h"

/**@name Plain ART
   This module provides the Extra paramter structure (empty in fact)
   to run the ART process with single particles. The main function of this
   module is the single step ART process.
*/
//@{
/* ART parameters ---------------------------------------------------------- */
/** Plain ART parameters.
    This class contains extra information for: Plain ART.
    In this case nothing. Notice that mainly this classes are created to
    support other versions of ART (crystals, ...) */
class Plain_ART_Parameters {
    /* User parameters ..................................................... */
    /* Side information .................................................... */
public:
    /// cout << eprm;
    friend ostream & operator << (ostream &o,
       const Plain_ART_Parameters &eprm);

    /// Produce Plain side information from the Extra parameters
    void produce_Side_Info(const Basic_ART_Parameters &prm,
       GridVolume &vol_basis0);
};

/** Run a single step of ART.
    An ART iteration is compound of as many steps as projections,
    this function runs a single step of the process. In ART the
    pointer to the output volume must point to the same vol_in,
    while in SIRT it should point to a second volume. The read projection
    must be provided to the algorithm but the rest of projections
    are output by the routine. The mean error is also an output.
    numIMG is a normalizing factor to be used in SIRT, if you are
    running pure ART then this factor should be 1.

    The symmetry matrix from which the view is derived must be given in
    sym_no. In fact, it is not used in this version of ART, but it is
    needed for the crystal counterpart. */
void ART_single_step(GridVolume &vol_in, GridVolume *vol_out,
   Basic_ART_Parameters &prm, Plain_ART_Parameters &eprm,
   Projection &theo_proj, Projection &read_proj,
   int sym_no,
   Projection &diff_proj, Projection &corr_proj, Projection &alig_proj,
   double &mean_error, int numIMG, double lamdba, int act_proj,
   const FileName &fn_ctf);

/** Update residual vector for WLS ART */
void update_residual_vector(Basic_ART_Parameters &prm, GridVolume &vol_basis,
			    double &kappa, double &pow_residual_vol, double &pow_residual_imgs);

/** Finish iterations.
    For WLS: delete residual images
    Else: do nothing. */
void finish_ART_iterations(const Basic_ART_Parameters &prm,
   const Plain_ART_Parameters &eprm, GridVolume &vol_basis);

/** Force the {\it trial} volume to be symmetric. So far only implemented
    for crystals.*/
void apply_symmetry(GridVolume &vol_in, GridVolume *vol_out,
                    const Plain_ART_Parameters &eprm,
		    int grid_type);
//@}
#endif
