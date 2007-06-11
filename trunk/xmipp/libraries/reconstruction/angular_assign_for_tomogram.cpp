/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.uam.es (2006)
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

#include "angular_assign_for_tomogram.h"
#include "projection.h"

#include <data/args.h>
#include <data/docfile.h>

#include <algorithm>

// Empty constructor =======================================================
Prog_angular_predict_tomography_prm::Prog_angular_predict_tomography_prm()
{
    each_image_produces_an_output = true;
}

// Read arguments ==========================================================
void Prog_angular_predict_tomography_prm::read(int argc, char **argv)
{
    Prog_parameters::read(argc, argv);
    fn_ref = getParameter(argc, argv, "-ref");
    fn_out_ang = getParameter(argc, argv, "-oang");
    max_rot_change = textToFloat(getParameter(argc, argv, "-max_rot_change", "5"));
    max_tilt_change = textToFloat(getParameter(argc, argv, "-max_tilt_change", "2"));
    max_psi_change = textToFloat(getParameter(argc, argv, "-max_psi_change", "5"));
    rot_step = textToFloat(getParameter(argc, argv, "-rot_step", "1"));
    tilt_step = textToFloat(getParameter(argc, argv, "-tilt_step", "1"));
    psi_step = textToFloat(getParameter(argc, argv, "-psi_step", "1"));
    max_shift_change = textToFloat(getParameter(argc, argv, "-max_shift_change", "3"));
    shift_step = textToFloat(getParameter(argc, argv, "-shift_step", "1"));
    onlyX = checkParameter(argc, argv, "-onlyX");
    onlyY = checkParameter(argc, argv, "-onlyY");
    produce_side_info();
}

// Show ====================================================================
void Prog_angular_predict_tomography_prm::show()
{
    Prog_parameters::show();
    cout << "Reference images: " << fn_ref << endl
    << "Ouput angular file: " << fn_out_ang << endl
    << "Max rot change: " << max_rot_change << " step: " << rot_step << endl
    << "Max tilt change: " << max_tilt_change << " step: " << tilt_step << endl
    << "Max psi change: " << max_psi_change << " step: " << psi_step << endl
    << "Max shift change: " << max_shift_change << " step: " << shift_step << endl
    << "OnlyX: " << onlyX << endl
    << "OnlyY: " << onlyY << endl
    ;
}

// usage ===================================================================
void Prog_angular_predict_tomography_prm::usage()
{
    Prog_parameters::usage();
    cerr << "   -ref <volume>             : Reference volume\n"
    << "   -oang <angle file>        : DocFile with output angles\n"
    << "  [-max_rot_change <ang=5>]  : Maximum change allowed in rot\n"
    << "  [-max_tilt_change <ang=2>] : Maximum change allowed in tilt\n"
    << "  [-max_psi_change <ang=5>]  : Maximum change allowed in psi\n"
    << "  [-max_shift_change <r=3>]  : Maximum change allowed in shift\n"
    << "  [-rot_step <ang=1>]        : Rot search step\n"
    << "  [-tilt_step <ang=1>]       : Tilt search step\n"
    << "  [-psi_step <ang=3>]        : Psi search step\n"
    << "  [-shift_step <r=1>]        : Step in shift in pixels\n"
    << "  [-onlyX]                   : Apply shifts only in X\n"
    << "  [-onlyY]                   : Apply shifts only in Y\n"
    ;
}

// Produce side information ================================================
void Prog_angular_predict_tomography_prm::produce_side_info()
{
    V.read(fn_ref);
    V().setXmippOrigin();
}

// Predict shift and psi -----------------------------------------------------
//#define DEBUG
double Prog_angular_predict_tomography_prm::predict_angles(ImageXmipp &I,
        double &assigned_shiftX, double &assigned_shiftY,
        double &assigned_rot, double &assigned_tilt, double &assigned_psi)
{
    vector<Alignment> list_of_alignments;

    for (double rot = -max_rot_change; rot <= max_rot_change; rot += rot_step)
        for (double tilt = I.tilt() - max_tilt_change;
             tilt <= I.tilt() + max_tilt_change;
             tilt += tilt_step)
        {
            // Take a projection from the given direction
            Projection theo;
            project_Volume(V(), theo, YSIZE(V()), XSIZE(V()), rot, tilt, 0);
            double theo_avg, theo_stddev, min_val, max_val;
            theo().computeStats(theo_avg, theo_stddev, min_val, max_val);
            theo() -= theo_avg;

            // Compare it to all possible rotations and shifts of the experimental
            // image
            ImageXmipp Ip;
            Matrix1D<double> shift(2);
            double max_shift_X, max_shift_Y;
            if (onlyX)
            {
                max_shift_Y = 0;
                max_shift_X = max_shift_change;
            }
            else if (onlyY)
            {
                max_shift_X = 0;
                max_shift_Y = max_shift_change;
            }
            else
            {
                max_shift_X = max_shift_Y = max_shift_change;
            }
            for (double x = I.Xoff() - max_shift_X; x <= I.Xoff() + max_shift_X; x += shift_step)
                for (double y = I.Yoff() - max_shift_Y; y <= I.Yoff() + max_shift_Y; y += shift_step)
                {
                    if ((x - I.Xoff())*(x - I.Xoff()) + (y - I.Yoff())*(y - I.Yoff()) > max_shift_change*max_shift_change) continue;
                    bool first_psi = true;
                    Alignment best_psi_alignment;
                    for (double psi = I.psi() - max_psi_change; psi <= I.psi() + max_psi_change; psi += psi_step)
                    {
                        // Shift image if necessary
                        if (x == 0 && y == 0) Ip() = I();
                        else
                        {
                            VECTOR_R2(shift, x, y);
                            I().translate(shift, Ip(), DONT_WRAP);
                        }

                        // Rotate image if necessary
                        // Adding 2 is a trick to avoid that the 0, 90, 180 and 270
                        // are treated in a different way
                        if (psi != 0)
                        {
                            Ip().selfRotate(psi + 2, DONT_WRAP);
                            Ip().selfRotate(-2, DONT_WRAP);
                        }

                        // Compute the correlation index
                        double read_avg, read_stddev;
                        Ip().computeStats(read_avg, read_stddev, min_val, max_val);
                        double correlation_index = 0;
                        FOR_ALL_ELEMENTS_IN_MATRIX2D(Ip())
                        correlation_index += (Ip(i, j) - read_avg) * theo(i, j);
                        correlation_index /= XSIZE(Ip()) * YSIZE(Ip());
                        correlation_index /= read_stddev * theo_stddev;

                        // Keep the value
                        Alignment A;
                        A.rot = rot;
                        A.tilt = tilt;
                        A.psi = psi;
                        A.x = x;
                        A.y = y;
                        A.corr = correlation_index;

                        if (A > best_psi_alignment || first_psi)
                        {
                            first_psi = false;
                            best_psi_alignment = A;
                        }

#ifdef DEBUG
                        cout << A;
                        ImageXmipp save;
                        save() = theo() - theo_avg;
                        save.write("PPPtheo.xmp");
                        save() = Ip() - read_avg;
                        save.write("PPPexp.xmp");
                        multiplyElements((theo() - theo_avg), (Ip() - read_avg), save());
                        save.write("PPPprod.xmp");
                        char c;
                        cin >> c;
#endif
                    }
                    list_of_alignments.push_back(best_psi_alignment);
                }
        }

    sort(list_of_alignments.begin(), list_of_alignments.end());

    // Select best alignment
    int imax = list_of_alignments.size() - 1;
    assigned_rot = list_of_alignments[imax].rot;
    assigned_tilt = list_of_alignments[imax].tilt;
    assigned_psi = list_of_alignments[imax].psi;
    assigned_shiftX = list_of_alignments[imax].x;
    assigned_shiftY = list_of_alignments[imax].y;
    list_of_assigned.push_back(list_of_alignments[imax]);
    image_name.push_back(I.name());
    return list_of_alignments[imax].corr;
}
#undef DEBUG

// Finish processing ---------------------------------------------------------
void Prog_angular_predict_tomography_prm::finish_processing()
{
    // Save predicted angles
    int p = list_of_assigned.size();
    DocFile DF;
    DF.reserve(p + 1);
    DF.append_comment("Predicted_Rot Predicted_Tilt Predicted_Psi Predicted_ShiftX Predicted_ShiftY Corr");
    Matrix1D<double> v(6);
    for (int i = 0; i < p; i++)
    {
        v(0) = list_of_assigned[i].rot;
        v(1) = list_of_assigned[i].tilt;
        v(2) = list_of_assigned[i].psi;
        v(3) = list_of_assigned[i].x;
        v(4) = list_of_assigned[i].y;
        v(5) = list_of_assigned[i].corr;
        DF.append_comment(image_name[i]);
        DF.append_data_line(v);
    }
    DF.write(fn_out_ang);
}
