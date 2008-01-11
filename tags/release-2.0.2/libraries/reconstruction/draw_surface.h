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
#ifndef _PROG_DRAW_SURFACE_HH
#  define _PROG_DRAW_SURFACE_HH

#include <data/funcs.h>
#include <data/image.h>
#include <data/volume.h>

/**@defgroup DrawSurface draw_surface (Draw a surface on a volume)
   @ingroup ReconsLibraryPrograms */
//@{
/* Surface Program Parameters ---------------------------------------------- */
/** Parameter class for the project program */
class Prog_Draw_Surface_Parameters
{
public:
    /// Filename with the input volume.
    FileName fn_in;
    /** Filename with the output volume. Maybe empty and then the input one
        is used */
    FileName fn_out;
    /** Surface to draw */
    FileName fn_surf;
    /** Force surface size */
    bool     enable_adjust;
    /** New forced surface size */
    int      ztop, zbottom;
    /** Colouring density */
    float    color;

    /// Volume where to draw
    VolumeXmipp vol;
    /// Surface
    ImageXmipp surface;
public:
    /** Read from a command line.
        An exception might be thrown by any of the internal conversions,
        this would mean that there is an error in the command line and you
        might show a usage message. */
    void read(int argc, char **argv);

    /** Usage message.
        This function shows the way of introducing this parameters. */
    void usage() const;

    /** Produce Side Information.
        Read phantom file and assign ztop and zbottom if they
        are not assigned by user, ie, set it to zdim. */
    void produce_Side_Info();
};

/** Draw surface on volume.
    The given color is assigned to the volumem voxel specified by the
    surface (V(S(i,j),i,j)=color). An exception is thrown if the
    surface is not exactly of the same shape as the volume XY plane. */
void draw_surface(Volume *vol, const Image *surf, float color);

/** Run draw_surface.
    Produces side information, draw the surface and save output volume*/
void ROUT_draw_surface(Prog_Draw_Surface_Parameters &prm);
//@}
#endif
