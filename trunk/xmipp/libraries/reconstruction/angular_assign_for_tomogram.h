/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2006)
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
#ifndef _PROG_ANGULAR_PREDICT_TOMOGRAPHY
#define _PROG_ANGULAR_PREDICT_TOMOGRPAHY

#include <data/funcs.h>
#include <data/selfile.h>
#include <data/volume.h>

/**@defgroup AngularPredictTomography angular_assign_for_tomogram (Discrete angular assignment for tomography)
   @ingroup ReconsLibraryPrograms */
//@{
/** Alignment structure */
class AlignmentTomography
{
public:
    double rot;
    double tilt;
    double psi;
    double x;
    double y;
    double corr;
    FileName fn_img;
    FileName fn_mask;
    
    AlignmentTomography();
};

/** Angular Predict parameters. */
class Prog_angular_predict_tomography_prm
{
public:
    /** Filename of the reference volume */
    FileName fn_ref;
    /** Filename of the images */
    FileName fn_sel;
    /** Filename of the masks */
    FileName fn_masksel;
    /** Root filename for the output */
    FileName fn_out;
    /** Maximum rotation change. */
    double max_rot_change;
    /** Max tilt change */
    double max_tilt_change;
    /** Maximum psi change. */
    double max_psi_change;
    /** Rot step */
    double rot_step;
    /** Tilt step */
    double tilt_step;
    /** Psi step */
    double psi_step;
    /** Maximum shift change */
    double max_shift_change;
    /** Shift step */
    double shift_step;
    /** Adjust gray */
    bool adjustGray;
public:
    VolumeXmipp V;
    std::vector<AlignmentTomography> list_of_assigned;
public:
    /// Read argument from command line
    void read(int argc, char **argv);

    /// Show
    void show();

    /// Usage
    void usage();

    /** Produce side info.
        An exception is thrown if any of the files is not found*/
    void produce_side_info();

    /** Predict angles and shift.
        This function searches in the shift-psi space and for each combination
        it correlates with the whole reference set. */
    void predict_angles(int i);

    /** Run the algorithm over all images */
    void run();
};
//@}
#endif
