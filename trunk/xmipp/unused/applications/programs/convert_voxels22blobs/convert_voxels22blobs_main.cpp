/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2000)
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

#include <data/blobs.h>
#include <data/grids.h>

void Usage();

int main(int argc, char *argv[])
{
    FileName fn_in, fn_out;
    bool   voxels_to_blobs;
    struct blobtype blob;
    float lambda;
    float final_error;
    float grid_relative_size;
#define  CC 0
#define FCC 1
#define BCC 2
    int grid_type;
    double R;

// Read the command line ---------------------------------------------------
    try
    {
        if (checkParameter(argc, argv, "-voxels"))
        {
            voxels_to_blobs = true;
            fn_in = getParameter(argc, argv, "-voxels");
            grid_relative_size = textToFloat(getParameter(argc, argv, "-g", "1.41"));
            if (checkParameter(argc, argv, "-FCC")) grid_type = FCC;
            else if (checkParameter(argc, argv, "-CC"))  grid_type = CC;
            else                                      grid_type = BCC;
        }
        else if (checkParameter(argc, argv, "-blobs"))
        {
            voxels_to_blobs = false;
            fn_in = getParameter(argc, argv, "-blobs");
        }
        else
            REPORT_ERROR(ERR_ARG_INCORRECT, "Not recognised input file type");
        fn_out = getParameter(argc, argv, "-o");
        lambda             = textToFloat(getParameter(argc, argv, "-l",    "0.05"));
        final_error        = textToFloat(getParameter(argc, argv, "-final_error", "0.01"));
        blob.radius        = textToFloat(getParameter(argc, argv, "-r",    "2"));
        blob.order         = textToInteger(getParameter(argc, argv, "-m",    "2"));
        blob.alpha         = textToFloat(getParameter(argc, argv, "-a",    "10.4"));
        R                  = textToFloat(getParameter(argc, argv, "-R",    "-1"));
    }
    catch (XmippError &XE)
    {
        std::cout << XE;
        Usage();
        exit(1);
    }
    try
    {
// Really convert ----------------------------------------------------------
        GridVolume  vol_blobs;
        Image<double> vol_voxels;
        if (voxels_to_blobs)
        {
            vol_voxels.read(fn_in);
            vol_voxels().setXmippOrigin();
            voxels2blobs(&(vol_voxels()), blob, vol_blobs, grid_type,
                         grid_relative_size, lambda, NULL, NULL, final_error, false, R);
            vol_blobs.write(fn_out);
        }
        else
        {
            vol_blobs.read(fn_in,"blobs");
            blobs2voxels(vol_blobs, blob, &(vol_voxels()));
            vol_voxels.write(fn_out);
        }
    }
    catch (XmippError XE)
    {
        std::cout << XE;
    }

}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    std::cout << "Usage: Voxels22blobs [Parameters]\n"
    << "   (-voxels | -blobs) <file_in>       : Input file\n"
    << "    -o <file_out>                     : of the opposite type\n"
    << "   [-r <blob radius=2>]               : blob radius\n"
    << "   [-m <blob order=2>]                : blob derivative order\n"
    << "   [-a <blob alpha=10.4>]             : controls smoothness\n"
    << "Only if voxels:\n"
    << "   [-g <grid_relative_size=1.41>]     : size between grid samples\n"
    << "   [-FCC | -CC]                       : by default, BCC grid\n"
    << "   [-l <lambda=0.05>]                 : convergence rate\n"
    << "   [-final_error <error=0.01>]        : minimum change percentage\n"
    << "   [-R <R=-1>]                        : interest radius\n";
}
