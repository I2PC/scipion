/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (2002)
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

#include <data/volume.h>
#include <data/args.h>
#include <data/mask.h>
#include <data/filters.h>
#include <data/geometry.h>

#include <cstdio>

void Usage();

int main(int argc, char **argv)
{
    FileName        fn_input;
    Mask_Params     mask_prm(INT_MASK);
    double          rot0,  rotF,  step_rot;
    double          tilt0, tiltF, step_tilt;
    int             rot_sym;
    bool            centerVolume;

    VolumeXmipp     volume;

    // Read arguments --------------------------------------------------------
    try
    {
        fn_input = get_param(argc, argv, "-i");
        rot_sym = AtoI(get_param(argc, argv, "-rot_sym", "0"));
        centerVolume = check_param(argc, argv, "-center_volume");
        int i;
        if ((i = position_param(argc, argv, "-rot")) != -1)
        {
            if (i + 3 >= argc)
                REPORT_ERROR(1, "findcenter3D: Not enough parameters behind -rot");
            rot0    = AtoF(argv[i+1]);
            rotF    = AtoF(argv[i+2]);
            step_rot = AtoF(argv[i+3]);
        }
        else
        {
            rot0 = 0;
            rotF = 355;
            step_rot = 5;
        }
        if ((i = position_param(argc, argv, "-tilt")) != -1)
        {
            if (i + 3 >= argc)
                REPORT_ERROR(1, "findcenter3D: Not enough parameters behind -tilt");
            tilt0    = AtoF(argv[i+1]);
            tiltF    = AtoF(argv[i+2]);
            step_tilt = AtoF(argv[i+3]);
        }
        else
        {
            tilt0 = 0;
            tiltF = 90;
            step_tilt = 5;
        }
        mask_prm.read(argc, argv);
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        Usage();
        mask_prm.usage();
        exit(1);
    }

    // Find Center and symmetry elements ------------------------------------
    try
    {
        // Read input volumes
        volume.read(fn_input);
        volume().set_Xmipp_origin();
        mask_prm.generate_3Dmask(volume());
        VolumeXmipp volume_sym, volume_aux;

        // Compute center of mass
        matrix1D<double> center_of_mass;
        volume().center_of_mass(center_of_mass, &mask_prm.get_binary_mask3D());
        cout << "Center of mass (X,Y,Z)= " << center_of_mass.transpose() << endl;

        // Move origin to that center of mass
        volume().self_translate(-center_of_mass, DONT_WRAP);
        if (centerVolume) volume.write();

        // Look for the rotational symmetry axis
        if (rot_sym > 1)
        {
            double best_corr = 0, best_rot = rot0 - step_rot, best_tilt = tilt0 - step_tilt;
            matrix2D<double> Euler;
            matrix1D<double> sym_axis(3);
            int maxsteps = FLOOR((rotF - rot0) / step_rot) * FLOOR((tiltF - tilt0) / step_tilt);
            cerr << "Searching symmetry axis ...\n";
            init_progress_bar(maxsteps);
            int i = 0;
            for (double rot = rot0; rot <= rotF; rot += step_rot)
                for (double tilt = tilt0; tilt <= tiltF; tilt += step_tilt)
                {
                    // Compute symmetry axis
                    matrix2D<double> Euler;
                    Euler_angles2matrix(rot, tilt, 0, Euler);
                    matrix1D<double> sym_axis(3);
                    Euler.getRow(2, sym_axis);
                    sym_axis.self_transpose();

                    // Symmetrize along this axis
                    volume_sym() = volume();
                    for (int n = 1; n < rot_sym; n++)
                    {
                        matrix2D<double> sym_matrix;
                        sym_matrix = rot3D_matrix(360.0 / rot_sym * n, sym_axis);
                        apply_geom(volume_aux(), sym_matrix, volume(), IS_NOT_INV,
                                   DONT_WRAP);
                        volume_sym() += volume_aux();
                    }

                    // Measure correlation
                    double corr = correlation_index(volume(), volume_sym(),
                                                    &mask_prm.get_binary_mask3D());
                    if (corr > best_corr)
                    {
                        best_corr = corr;
                        best_rot = rot;
                        best_tilt = tilt;
                    }

                    // progress bar
                    if ((i++) % MAX(maxsteps / 60, 1) == 0) progress_bar(i);
                }
            progress_bar(maxsteps);
            cout << "Symmetry axis (rot,tilt)= " << best_rot << " "
            << best_tilt << " --> ";
            Euler_angles2matrix(best_rot, best_tilt, 0, Euler);
            Euler.getRow(2, sym_axis);
            cout << sym_axis << endl;
        }

    }
    catch (Xmipp_error XE)
    {
        cout << XE;
    }
    exit(0);
} //main

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    cerr << "Purpose:\n";
    cerr << "    Finds the 3D center of mass within a mask\n"
    << "    and a symmetry rotational axis passing through that center\n";

    cerr << "Usage: findcenter3D [options]" << endl
    << "    -i <volume>                         : volume to process\n"
    << "   [-center_volume]                     : save the centered volume\n"
    << "   [-rot_sym <n=0>]                     : order of the rotational axis\n"
    << "   [-rot  <rot0=0>  <rotF=355> <step=5>]: limits for rotational axis\n"
    << "   [-tilt <tilt0=0> <tiltF=90> <step=5>]: limits for rotational axis\n"
    ;
}
