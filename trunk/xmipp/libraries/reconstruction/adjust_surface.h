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
#ifndef _PROG_ADJUST_SURFACE_HH
#  define _PROG_ADJUST_SURFACE_HH

#include <data/image.h>
#include <data/volume.h>
#include <data/funcs.h>
#include <data/vectorial.h>

/**@name Adjust Surface program */
//@{
/* Adjust Surface Program Parameters --------------------------------------- */
/** Parameter class for the project program */
class Prog_Adjust_Surface_Parameters
{
public:
    /// Filename with the input volume.
    FileName fn_vol;
    /// Filename with the input surface.
    FileName fn_in_surface;
    /// Filename with the output surface.
    FileName fn_out_surface;
    /// Exhaustive search
    bool exhaustive;
    /// Apply best correlation
    bool apply;
#define TOP2BOTTOM 1
#define BOTTOM2TOP 2
    /// It's a phantom
    bool phantom;
    /// Direction. Either TOP2BOTTOM or BOTTOM2TOP
    int direction;
    /// User given ztop
    bool given_ztop;
    /// Minimum height to search the top surface
    int ztop0;
    /// Maximum height to search the top surface
    int ztopF;
    /// Step for ztop
    int ztop_step;
    /// Minimum angle to search the surface
    double angle0;
    /// Maximum angle to search the surface
    double angleF;
    /// Angle step
    double angle_step;
    /// Initial shift X
    double shiftX0;
    /// Initial shift Y
    double shiftY0;
    /// Maximum shift allowed
    double shift_dist;
    /// Shift step
    int  shift_step;
    /// User given zbottom
    bool given_zbottom;
    /// Minimum height to search the bottom surface
    int zbottom0;
    /// Maximum height to search the bottom surface
    int zbottomF;
    /// Step for zbottom
    int zbottom_step;
    /// Initial scale factor in X
    double scaleX0;
    /// Final scale factor in X
    double scaleXF;
    /// ScaleX step
    double scaleX_step;
    /// Initial scale factor in Y
    double scaleY0;
    /// Final scale factor in Y
    double scaleYF;
    /// ScaleY step
    double scaleY_step;
#define CORR_2D      0x2
#define CORR_GRAD    0x4
#define MANUAL_ORDER 0x8
    /// Tell
    int tell;
    /**@name Side parameters */
    //@{
    /// Input volume
    VolumeXmipp V;
    /// Surface to adjust
    ImageXmipp surface;
    /// Min surface value
    double min_val;
    /// Max surface value
    double max_val;
    /// Shift mask
    matrix2D<int> shift_mask;

    /// Gradient of input volume
    Vectorial_matrix3D V_grad;

    /// Working surface volume
    VolumeXmipp Vsurf;
    /// Working surface image
    ImageXmipp wsurface;
    /// Working gradient volume
    Vectorial_matrix3D Vsurf_grad;
    /// Working best combination
    matrix1D<double> p;
    //@}
public:
    /** Read from a command line.
        An exception might be thrown by any of the internal conversions,
        this would mean that there is an error in the command line and you
        might show a usage message.*/
    void read(int argc, char **argv);

    /** Usage message.
        This function shows the way of introducing this parameters. */
    void usage() const;

    /** Produce Side Information.
        Read phantom file and assign ztop and zbottom if they
        are not assigned by user, ie, set it to zdim. An exception is thrown
        if ztop or zbottom exceeds volume limits*/
    void produce_Side_Info();
};

/** Create surface mask.
    From a surface image, create a volume mask with the desired surface */
void create_surface_mask(const Image *surf, const Volume *V, Volume *Vsurf,
                         int direction);

/** Correlate surface and volume (3D).
    This function constructs a mask volume according to the information
    given by mask. Then the correlation is computed of the two volumes
    between the two heights given. The correlation coefficient is returned.

    The direction indicates in which direction the volume mask will be
    extended from the surface, it can take TOP2BOTTOM or BOTTOM2TOP. An
    exception is thrown if non of these directions is specified.

    tell can take the MANUAL_ORDER flag, then the intermidiate volumes
    are written as PPPsurface.vol with the generated volumetric surface
    mask and PPPsign.vol with the differences in the interest planes
    between the input and the mask volumes*/
double correlate_surface_and_volume_3D(const Image *surf, const Volume *V,
                                       Volume *Vsurf, int ktop, int kbottom, int direction, int tell = 0)
;

/** Correlate surface and volume (2D).
    This function constructs a mask volume according to the information
    given by mask. Then the correlation is computed of the two volumes
    between the two heights given. The correlation coefficient is returned.

    The direction indicates in which direction the volume mask will be
    extended from the surface, it can take TOP2BOTTOM or BOTTOM2TOP. An
    exception is thrown if non of these directions is specified.

    tell can take the MANUAL_ORDER flag, then the intermidiate projection
    are written as PPPsurface.img with the projection of the interest region
    and PPPsign.img with the differences between the input surface and
    thr projection*/
double correlate_surface_and_volume_2D(const Image *surf, const Volume *V,
                                       Volume *Vsurf, int ktop, int kbottom, int direction, int tell = 0)
;

/** Correlate surface and volume (gradients).
    The same as above but this time the correlation is perfomed attending
    to the gradient at surface positions inside the volume. */
double correlate_surface_and_volume_gradients(const Image *surf,
        const Vectorial_matrix3D &V_grad, Vectorial_matrix3D & Vsurf_grad,
        int ktop, int kbottom, int direction, int tell = 0)
;

/** Run adjust surface.
    Search for the best fit of the surface within the volume */
void ROUT_adjust_surface(Prog_Adjust_Surface_Parameters &prm);
//@}
#endif
