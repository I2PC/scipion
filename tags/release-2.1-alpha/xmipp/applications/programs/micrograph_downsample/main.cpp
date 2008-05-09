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

#include <data/micrograph_downsample.h>
#include <data/args.h>
#include <data/xvsmooth.h>

void Usage(const Prog_downsample_prm &prm);

int main(int argc, char **argv)
{
    Prog_downsample_prm prm;
    bool                smooth;
    bool                reversed;

    // Get input parameters -------------------------------------------------
    try
    {
        prm.read(argc, argv);
        smooth         = checkParameter(argc, argv, "-smooth");
        reversed       = checkParameter(argc, argv, "-reverse_endian");
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage(prm);
        exit(1);
    }

    try
    {
        prm.generate_kernel();
        prm.open_input_micrograph();
        prm.create_empty_output_file();
        if (smooth)
        {
            Micrograph Mp;
            Mp.open_micrograph(prm.fn_downsampled, reversed);
            byte rgb[256];
            for (int i = 0; i < 256; i++) rgb[i] = i;
            byte *result = SmoothResize((byte *)(prm.M.array8()),
                                        prm.Xdim, prm.Ydim, prm.Xpdim, prm.Ypdim,
                                        rgb, rgb, rgb, rgb, rgb, rgb, 256);
            for (int i = 0; i < prm.Ypdim; i++)
                for (int j = 0; j < prm.Xpdim; j++)
                    Mp.set_val(j, i, result[i*prm.Xpdim+j]);
            Mp.close_micrograph();
        }
        else prm.Downsample();
        prm.close_input_micrograph();
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
}

/* Usage =================================================================== */
void Usage(const Prog_downsample_prm &prm)
{
    std::cerr << "Purpose: This file allows you to downsample raw images\n"
    << "Usage: downsample [parameters]\n"
    << "   -i <input_file>        : Raw input file, <input_file>.inf\n"
    << "                            must exist\n"
    << "   -o <output_file>       : Must be different from input one\n"
    << "  [-smooth]               : Use Smoothing for downsampling\n"
    ;
    prm.usage();
}
