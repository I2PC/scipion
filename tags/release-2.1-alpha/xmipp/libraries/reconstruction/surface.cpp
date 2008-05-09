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

#include "surface.h"

/* Read from command line ================================================== */
void Prog_Surface_Parameters::read(int argc, char **argv)
{
    fnPhantom     = getParameter(argc, argv, "-i", "");
    probe_radius   = textToFloat(getParameter(argc, argv, "-r", "0.5"));
    fn_top         = getParameter(argc, argv, "-top", "");
    fn_bottom      = getParameter(argc, argv, "-bottom", "");
    fn_mask        = getParameter(argc, argv, "-o", "");
    enable_ztop    = checkParameter(argc, argv, "-ztop");
    if (enable_ztop) ztop = textToFloat(getParameter(argc, argv, "-ztop"));
    enable_zbottom = checkParameter(argc, argv, "-zbottom");
    if (enable_zbottom) zbottom = textToFloat(getParameter(argc, argv, "-zbottom"));
    zdim           = textToInteger(getParameter(argc, argv, "-zdim", "0"));

    if (fn_top == "" && fn_bottom == "")
        REPORT_ERROR(1, "Prog_Surface_Parameters::read: No surface given!!!");
    if (fnPhantom == "" && zdim <= 0)
        REPORT_ERROR(1, "Prog_Surface_Parameters::read: No valid Z dimension");
}

/* Usage =================================================================== */
void Prog_Surface_Parameters::usage() const
{
    std::cout << "\nUsage:\n";
    std::cout << "surface\n"
    << "  [-i <Phantom file>]             : Phantom description file\n"
    << "  [-o <volume_mask>]              : Output mask\n"
    << "  [-r <probe_radius=0.5>]         : Probe radius for surface generation\n"
    << "  [-top <top surface file>]       : Output Top surface\n"
    << "  [-bottom <bottom surface file>] : Output Bottom surface\n"
    << "  [-ztop <ztop>]                  : Maximum height for top ray\n"
    << "  [-zbottom <zbottom>]            : Maximum height for bottom ray\n"
    << "  [-zdim <zdim>]                  : Output Z dimension\n";
}

/* Produce side information ================================================ */
void Prog_Surface_Parameters::produce_Side_Info()
{
    if (fnPhantom != "")
    {
        phantom.read(fnPhantom);
        if (zdim == 0)         zdim   = phantom.zdim;
        if (!enable_ztop)    ztop   = phantom.zdim;
        if (!enable_zbottom) zbottom = phantom.zdim;
    }
}

/* Create surface mask ===================================================== */
#define VOL    (*V)()
void create_surface_mask(const Image *top, const Image *bottom, int zdim,
                         Volume *V)
{
    const Matrix2D<double> *surf;
    if (top != NULL && bottom != NULL)
    {
        if (!SAME_SHAPE2D((*top)(), (*bottom)()))
            REPORT_ERROR(1, "create_surface_mask: Top and bottom surfaces hasn't "
                         "got the same shape");
        surf = &((*top)());
    }
    else if (top != NULL) surf = &((*top)());
    else if (bottom != NULL) surf = &((*bottom)());
    else
        REPORT_ERROR(1, "create_surface_mask: Both surfaces are pointing to NULL");

    // Resize V
    VOL.resize(zdim, YSIZE(*surf), XSIZE(*surf));
    STARTINGZ(VOL) = FIRST_XMIPP_INDEX(zdim);
    STARTINGY(VOL) = STARTINGY(*surf);
    STARTINGX(VOL) = STARTINGX(*surf);
    VOL.initConstant(1);

    // Create mask
    int kmin = FIRST_XMIPP_INDEX(zdim);
    int kmax = LAST_XMIPP_INDEX(zdim);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(*surf)
    {
        int k0, kF;
        if (top != NULL)    k0 = FLOOR(IMGPIXEL(*top   , i, j)) + 1;
        else k0 = kmin;
        if (bottom != NULL) kF = CEIL(IMGPIXEL(*bottom, i, j)) - 1;
        else kF = kmax;
        for (int k = k0; k <= kF; k++)
            if (k >= kmin && k <= kmax) VOLVOXEL(*V, k, i, j) = 0;
    }
}

/* Create surface ========================================================== */
void ROUT_surface(Prog_Surface_Parameters &prm)
{
    ImageXmipp top_surface, bottom_surface;

    // Create surfaces from a phantom ---------------------------------------
    if (prm.fnPhantom != "")
    {
        // Create top surface
        if (prm.fn_top != "")
        {
            top_surface.adapt_to_size(prm.zdim, prm.zdim);
            prm.phantom.surface(prm.ztop, prm.probe_radius,
                                NEG_POS, (Image *) &top_surface);
            top_surface.write(prm.fn_top);
        }

        // Create bottom surface
        if (prm.fn_bottom != "")
        {
            bottom_surface.adapt_to_size(prm.zdim, prm.zdim);
            prm.phantom.surface(prm.zbottom, prm.probe_radius,
                                POS_NEG, (Image *) &bottom_surface);
            bottom_surface.write(prm.fn_bottom);
        }
        // Read surfaces --------------------------------------------------------
    }
    else
    {
        if (prm.fn_top != "")
        {
            top_surface.read(prm.fn_top);
            top_surface().setXmippOrigin();
        }
        if (prm.fn_bottom != "")
        {
            bottom_surface.read(prm.fn_bottom);
            bottom_surface().setXmippOrigin();
        }
    }

    // Create volume mask
    if (prm.fn_mask != "")
    {
        VolumeXmipp mask;
        int zdim_to_use;
        if (XSIZE(bottom_surface()) != 0 && XSIZE(top_surface()) != 0)
            create_surface_mask((Image *) &top_surface, (Image *) &bottom_surface,
                                prm.zdim, (Volume *) &mask);
        else if (XSIZE(top_surface()) != 0)
            create_surface_mask((Image *) &top_surface, NULL,
                                prm.zdim, (Volume *) &mask);
        else
            create_surface_mask(NULL, (Image *) &bottom_surface,
                                prm.zdim, (Volume *) &mask);
        mask.write(prm.fn_mask);
    }
}
