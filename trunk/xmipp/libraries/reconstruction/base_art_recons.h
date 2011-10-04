/***************************************************************************
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Joaquin Oton (joton@cnb.csic.es)
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

#ifndef BASE_ART_RECONS_H_
#define BASE_ART_RECONS_H_

#include "basic_art.h"
#include <data/xmipp_program.h>

/**@defgroup common ART Reconstruction stuff
   @ingroup ReconsLibrary
    The main difference between ART applied to different cases (single
    particles, crystals, ...) is the single step applied to each case.
    Most of the tasks in the ART are common to all ART processes. All
    these common tasks as well as the common parameters are comprised
    in the ARTReconsBase class. These common tasks are based on the existence
    of a BasicARTParameters class containing all the specific information
    for the ART process.
*/
//@{

/* ART Reconstruction  ------------------------------------------ */
/** ART Base Reconstruction.
    This class contains all the basic routines needed about the ART reconstruction
    process. This is the class used in order to reconstruct single particles.
    Any other type of reconstruction must inherit from this.
    */
class ARTReconsBase
{
public:
    BasicARTParameters artPrm;

    virtual ~ARTReconsBase()
    {}

    static void defineParams(XmippProgram * program, bool mpiMode = false)
    {
        BasicARTParameters::defineParams(program, mpiMode);
    }

    /* --- Virtual methods to be implemented by children --- */

    /* Read params from command line
     */
    virtual void readParams(XmippProgram * program);

    /* Show reconstruction information
     */
    virtual void print(std::ostream &o)const;

    /// Produce Plain side information from the Class parameters
    virtual void preIterations(GridVolume &vol_basis0, int level = FULL, int rank = -1);

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
    virtual void singleStep(GridVolume &vol_in, GridVolume *vol_out,
                            Projection &theo_proj, Projection &read_proj,
                            int sym_no,
                            Projection &diff_proj, Projection &corr_proj, Projection &alig_proj,
                            double &mean_error, int numIMG, double lambda, int act_proj,
                            const FileName &fn_ctf, const MultidimArray<int> *maskPtr,
                            bool refine);

    /** Finish iterations.
        For WLS: delete residual images
        Else: do nothing. */
    virtual void postIterations(GridVolume &vol_basis);

    /** Force the trial volume to be symmetric. So far only implemented
        for crystals.*/
    virtual void applySymmetry(GridVolume &vol_in, GridVolume *vol_out,int grid_type);

    /* --- Methods that do not have to be implemented by children --- */

    /** Write first part of ART history.
        This function writes all ART parameters, projection angles, symmetry
        matrices, and grid structure in the history handler provided by
        the side information. At the end the routine makes a call to the
        operator << of the Extra_ART_Parameters to show the specific part
        of the History.

        BasicARTParameters is not constant since things are written in
        \ref BasicARTParameters::fh_hist.*/
    void initHistory(const GridVolume &vol_basis0);

    /** Perform all ART iterations.
        This function performs the iterations according to the ART parameters,
        it needs the side information to be fully computed. It throws
        a lot of information to the screen and to the history file (side.fh_hist),
        specially this one must exist.

        The GridVolume must have as input an initial guess for the solution,
        and when this function finishes, it contains the final solution volume
        in basis.

        The rank is the identification number of the process running this function.
        If it is -1, the function is run in sequential mode. If it is 0, then
        it is the root process.

        See the \ref BasicARTParameters for more information
        about how to generate the iterations.

        This method is not virtual as it should be common for all ARTRecons classes.
    */
    void iterations(GridVolume &vol_basis, int rank = -1);

    friend std::ostream & operator<< (std::ostream &o, const ARTReconsBase& artRecons);

};



class SinPartARTRecons : public ARTReconsBase
{
    pthread_t *th_ids;

public:
    SinPartARTRecons()
    {}

    virtual ~SinPartARTRecons()
    {}

    void preIterations(GridVolume &vol_basis0, int level = FULL, int rank = -1);

    virtual void singleStep(GridVolume &vol_in, GridVolume *vol_out,
                            Projection &theo_proj, Projection &read_proj,
                            int sym_no,
                            Projection &diff_proj, Projection &corr_proj, Projection &alig_proj,
                            double &mean_error, int numIMG, double lambda, int act_proj,
                            const FileName &fn_ctf, const MultidimArray<int> *maskPtr,
                            bool refine);

    void postIterations(GridVolume &vol_basis);
}
;




//@}
#endif /* BASE_ART_RECONS_H_ */
