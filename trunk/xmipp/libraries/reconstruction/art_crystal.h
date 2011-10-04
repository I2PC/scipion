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
#ifndef _PROG_ART_CRYSTAL_HH
#  define _PROG_ART_CRYSTAL_HH

#include <data/xmipp_image.h>
#include <data/matrix1d.h>
#include <data/matrix2d.h>
#include <data/symmetries.h>

#include "base_art_recons.h"

/**@defgroup ARTCrystal art_crystal (ART for crystals)
   @ingroup ReconsLibrary */
//@{
/** ART+crystal parameters.
    Here only those specific parameters for crystals are found, the rest of
    parameters common with normal ART should be looked up in
    \ref BasicARTParameters */

class CrystalARTRecons: public ARTReconsBase
{
public:
    /* User parameters ...................................................... */
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
    Matrix1D<double> avox;
    /// Second lattice vector (voxel units)
    Matrix1D<double> bvox;
    /// Fill space, repeat unit cell all over the space
    bool fill_space;

    /* Side information ..................................................... */
    /// First lattice vector (BCC units)
    Matrix1D<double> a;
    /// Second lattice vector (BCC units)
    Matrix1D<double> b;
    /// First lattice vector approximated to integer numbers (BCC units)
    Matrix1D<double> aint;
    /// Second lattice vector approximated to integer numbers (BCC units)
    Matrix1D<double> bint;
    /// ai=aint/2 as double numbers
    Matrix1D<double> ai;
    /// bi=aint/2 as double numbers
    Matrix1D<double> bi;
    /// space_group
    int space_group;
    /** Unit cell mask.
        This mask is defined in the BCC space
        and it represent a parallelogram defined by points (-a-b)/2, (-a+b)/2,
        (a-b)/2, (a+b)/2. The reconstruction will be only performed for
        the basis inside this mask. Be careful that this is a 2D mask for a 3D
        reconstruction. */
    MultidimArray<int> unit_cell_mask;

public:

    static void defineParams(XmippProgram * program, const char* prefix=NULL, const char* comment=NULL);

    /** Read special parameters from command line. */
    void readParams(XmippProgram * program);

    /// std::cout << crystal_prm;
    void print(std::ostream &o) const;


    /** Produce Initial and Side information for ART+crystal.
        This function computes from the ART+crystal parameters things like
        the integer<-->float lattice vector matrices, the unit cell mask.
        The reconstructed volume size can be modified such that the unit cell
        mask fits.
    */
    void produceSideInfo(GridVolume &vol_basis0);

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
    void singleStep(GridVolume &vol_in, GridVolume *vol_out,
                    Projection &theo_proj, Projection &read_proj,
                    int sym_no,
                    Projection &diff_proj, Projection &corr_proj, Projection &alig_proj,
                    double &mean_error, int numIMG, double lambda, int act_proj,
                    const FileName &fn_ctf, const MultidimArray<int> *maskPtr,
                    bool refine);

    /* Finish ART iterations.
       Expand output volume to fill space if necessary. */
    void postIterations(GridVolume &vol_basis);

    /** Force the trial volume to be symmetric. */
    void applySymmetry(GridVolume &vol_in, GridVolume *vol_out,int grid_type);
};

/** Compute integer lattice vectors and passing matrix.
    Given two real lattice vectors, this function returns the corresponding
    integer vectors and the matrix V which passes from a to aint. All
    vectors are supposed to be 2x1, so V is 2x2. Before computing
    @code
       aint = Vinv*a; a=V*aint;
       bint = Vinv*b; b=V*bint;
    @endcode */
void computeIntegerLattice(const Matrix1D<double> &a,
                             const Matrix1D<double> &b,
                             double a_mag, double b_mag,
                             double ang_a2b_deg,
                             Matrix1D<double> &aint,
                             Matrix1D<double> &bint,
                             Matrix2D<double> &V,
                             int space_group);

/** Expand basis values to fill space.
    Copy basis values as a crystal in order to fill the whole space determined
    by the output volume. */
void expandToFillSpace(const BasicARTParameters &prm,
                          const CrystalARTRecons &eprm, GridVolume &vol);

//@}

#endif
