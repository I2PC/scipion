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

#include "draw_surface.h"

#include <data/args.h>

/* Read from command line -------------------------------------------------- */
void Prog_Draw_Surface_Parameters::read(int argc, char **argv)
{
    fn_in   = getParameter(argc, argv, "-i");
    fn_surf = getParameter(argc, argv, "-s");
    fn_out  = getParameter(argc, argv, "-o", "");
    enable_adjust = checkParameter(argc, argv, "-ztop");
    if (enable_adjust)
    {
        ztop = textToInteger(getParameter(argc, argv, "-ztop"));
        zbottom = textToInteger(getParameter(argc, argv, "-zbottom"));
    }
    color = textToFloat(getParameter(argc, argv, "-color", "2"));
}

/* Usage ------------------------------------------------------------------- */
void Prog_Draw_Surface_Parameters::usage() const
{
    std::cerr << "Usage: draw_surface [parameters]\n"
    << "   -i <volume>                   : volume where drawing surface\n"
    << "   -s <surface>                  : image with a surface coded\n"
    << "  [-o <output volume>]           : if not given the same as input\n"
    << "  [-ztop <zmin> -zbottom <zmax>] : force surface to be in this range\n"
    << "  [-color <density=2>]           : density value for surface\n";
}

/* Produce side information ------------------------------------------------ */
void Prog_Draw_Surface_Parameters::produce_Side_Info()
{
    vol.read(fn_in);
    surface.read(fn_surf);
    if (enable_adjust)
        surface().range_adjust(ztop, zbottom);
    vol().setXmippOrigin();
    surface().setXmippOrigin();
}

/* Draw surface in volume -------------------------------------------------- */
#define V VOLMATRIX(*vol)
#define S IMGMATRIX(*surf)
void draw_surface(Volume *vol, const Image *surf, float color)
{
    if (XSIZE(V) != XSIZE(S)         || YSIZE(V) != YSIZE(S) ||
        STARTINGX(V) != STARTINGX(S) || STARTINGY(V) != STARTINGY(S))
        REPORT_ERROR(1, "draw_surface: Volume and surface are of different size");
    FOR_ALL_ELEMENTS_IN_MATRIX2D(S)
    if (MAT_ELEM(S, i, j) >= STARTINGZ(V) && MAT_ELEM(S, i, j) <= FINISHINGZ(V))
        VOL_ELEM(V, ROUND(MAT_ELEM(S, i, j)), i, j) = color;
}

/* Routine Draw surface ---------------------------------------------------- */
void ROUT_draw_surface(Prog_Draw_Surface_Parameters &prm)
{
    prm.produce_Side_Info();
    draw_surface(&(prm.vol), &(prm.surface), prm.color);
    if (prm.fn_out != "") prm.vol.write(prm.fn_out);
    else                prm.vol.write(prm.fn_in);
}
