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

#include <reconstruction/reconstruct_art.h>
//#include <reconstruction/art_crystal.h>

int main(int argc, char *argv[])
{
// Variables
    Basic_ART_Parameters   art_prm;
    Plain_ART_Parameters   dummy;
//    Crystal_ART_Parameters crystal_art_prm;
    Image<double>        vol_voxels;
    GridVolume             vol_blobs;
    GridVolume             *vol_blobs_var = NULL;
    int                    crystal_mode;

// Read Art Parameters
    try
    {
        art_prm.read(argc, argv);
        // Crystal
        crystal_mode = checkParameter(argc, argv, "-crystal");
        if (crystal_mode)
        	REPORT_ERROR(1,"crystal art temporarily deactivated.");
//        	crystal_art_prm.read(argc, argv, art_prm);
    }
    catch (Xmipp_error &XE)
    {
        std::cout << XE;
        bool usage_more = checkParameter(argc, argv, "-more_help");
        if (usage_more)
        {
            art_prm.usage_more();
            //crystal_art_prm.usage_more();
        }
        else
            art_prm.usage();
        exit(1);
    }

// Call main ART routine
    try
    {
        if (!crystal_mode)
            Basic_ROUT_Art(art_prm, dummy, vol_voxels, vol_blobs);
        else
        	std::cerr << "never get here..." <<std::endl;

            //Basic_ROUT_Art(art_prm, crystal_art_prm, vol_voxels, vol_blobs);
        std::cerr.flush();
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        exit(1);
    }
    exit(0);
}
