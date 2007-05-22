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
#ifndef _PROG_ART_CRYSTAL_HH
#  define _PROG_ART_CRYSTAL_HH

#include <data/matrix1d.h>
#include <data/matrix2d.h>
#include <data/matrix3d.h>
#include <data/image.h>

#include "reconstruct_art.h"
#include "symmetries.h"

/**@name ART+crystal */
//@{
/** ART+crystal parameters.
    Here only those specific parameters for crystals are found, the rest of
    parameters common with normal ART should be looked up in
    \Ref{Prog_Art_Parameters} */
class Crystal_ART_Parameters
{
public:
    /* User parameters ...................................................... */
    /**@name User parameters */
    //@{
    /// First lattice vector module (user supplies it in arbitrary units
    /// but it is divided by sampling as soon as it is initialized)
    double a_mag;
    /// Second lattice vector module (user supplies it in arbitrary units
    /// but it is divided by sampling as soon as it is initialized)
    double b_mag;
    /// angle from a to b (degrees)
    double ang_a2b_deg;
    /// angle from x axis to a (degrees)
    double ang_x2a_deg;
    /// First lattice vector  (voxel units)
    matrix1D<double> avox;
    /// Second lattice vector (voxel units)
    matrix1D<double> bvox;
    /// Fill space, repeat unit cell all over the space
    bool fill_space;
    //@}

    /* Side information ..................................................... */
    /**@name Side information */
    //@{
    /**@name Crystal lattice related */
    //@{
    /// First lattice vector (BCC units)
    matrix1D<double> a;
    /// Second lattice vector (BCC units)
    matrix1D<double> b;
    /// First lattice vector approximated to integer numbers (BCC units)
    matrix1D<double> aint;
    /// Second lattice vector approximated to integer numbers (BCC units)
    matrix1D<double> bint;
    /// ai=aint/2 as double numbers
    matrix1D<double> ai;
    /// bi=aint/2 as double numbers
    matrix1D<double> bi;
    /// space_group
    int space_group;
    //@}
    /** Unit cell mask.
        This mask is defined in the BCC space
        and it represent a parallelogram defined by points (-a-b)/2, (-a+b)/2,
        (a-b)/2, (a+b)/2. The reconstruction will be only performed for
        the basis inside this mask. Be careful that this is a 2D mask for a 3D
        reconstruction. */
    matrix2D<int> unit_cell_mask;
    //@}
public:
    /// cout << crystal_prm;
    friend ostream & operator << (ostream &o,
                                  const Crystal_ART_Parameters &eprm);

    /** Read special parameters from a command line.
        This function reads the parameters from a command line
        defined by argc and argv. An exception might be thrown by any
        of the internal conversions, this would mean that there is
        an error in the command line and you might show a usage message. */
    void read(int argc, char **argv, Basic_ART_Parameters &prm);

    /** Usage message.
        This function shows the way of introducing these parameters. */
    void usage_more();

    /** Produce Initial and Side information for ART+crystal.
        This function computes from the ART+crystal parameters things like
        the integer<-->float lattice vector matrices, the unit cell mask.
        The reconstructed volume size can be modified such that the unit cell
        mask fits.
    */
    void produce_Side_Info(Basic_ART_Parameters &prm,
                           GridVolume &vol_basis0);
};

/** Compute integer lattice vectors and passing matrix.
    Given two real lattice vectors, this function returns the corresponding
    integer vectors and the matrix V which passes from a to aint. All
    vectors are supposed to be 2x1, so V is 2x2. Before computing
    \begin{verbatim}
       aint = Vinv*a; a=V*aint;
       bint = Vinv*b; b=V*bint;
    \end{verbatim} */
void compute_integer_lattice(const matrix1D<double> &a,
                             const matrix1D<double> &b,
                             double a_mag, double b_mag,
                             double ang_a2b_deg,
                             matrix1D<double> &aint,
                             matrix1D<double> &bint,
                             matrix2D<double> &V,
                             int space_group);

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
    sym_no. */
void ART_single_step(GridVolume &vol_in, GridVolume *vol_out,
                     Basic_ART_Parameters &prm, Crystal_ART_Parameters &eprm,
                     Projection &theo_proj, Projection &read_proj,
                     int sym_no,
                     Projection &diff_proj, Projection &corr_proj, Projection &align_proj,
                     double &mean_error, int numIMG, double lambda, int act_proj,
                     const FileName &fn_ctf);

/* Finish ART iterations.
   Expand output volume to fill space if necessary. */
void finish_ART_iterations(const Basic_ART_Parameters &prm,
                           const Crystal_ART_Parameters &eprm, GridVolume &vol_basis);

/** Expand basis values to fill space.
    Copy basis values as a crystal in order to fill the whole space determined
    by the output volume. */
void expand_to_fill_space(const Basic_ART_Parameters &prm,
                          const Crystal_ART_Parameters &eprm, GridVolume &vol);

/** Force the {\it trial} volume to be symmetric. */
void apply_symmetry(GridVolume &vol_in, GridVolume *vol_out,
                    const Crystal_ART_Parameters &eprm,
                    int grid_type);

//@}

#endif
