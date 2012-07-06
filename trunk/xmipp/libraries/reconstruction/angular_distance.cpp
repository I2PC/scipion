/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2002)
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

#include "angular_distance.h"

#include <data/args.h>
#include <data/histogram.h>

// Read arguments ==========================================================
void ProgAngularDistance::readParams()
{
    fn_ang1 = getParam("--ang1");
    fn_ang2 = getParam("--ang2");
    fn_out = getParam("--oroot");
    fn_sym = getParam("--sym");
    check_mirrors = checkParam("--check_mirrors");
    object_rotation = checkParam("--object_rotation");
}

// Show ====================================================================
void ProgAngularDistance::show()
{
    std::cout
    << "Angular docfile 1: " << fn_ang1       << std::endl
    << "Angular docfile 2: " << fn_ang2       << std::endl
    << "Angular output   : " << fn_out    << std::endl
    << "Symmetry file    : " << fn_sym        << std::endl
    << "Check mirrors    : " << check_mirrors << std::endl
    << "Object rotation  : " << object_rotation<<std::endl
    ;
}

// usage ===================================================================
void ProgAngularDistance::defineParams()
{
    addUsageLine("Computes the angular distance between two angle files. The angular distance ");
    addUsageLine("is defined as the average angular distance between the 3 vectors of the ");
    addUsageLine("coordinate system defined by the Euler angles (taking into account any ");
    addUsageLine("possible symmetry). ");
    addParamsLine("   --ang1 <Metadata1>        : Angular document file 1");
    addParamsLine("   --ang2 <Metadata2>        : Angular document file 2");
    addParamsLine("  [--oroot <rootname=\"\">]  : Output rootname");
    addParamsLine("                             : rootname.xmd Angular comparison;");
    addParamsLine("                             : rootname_vec_diff_hist.txt Histogram of the differences in vector directions;");
    addParamsLine("                             : rootname_shift_diff_hist.txt Histogram of the differences in shifts;");
    addParamsLine("                             :+ rootname_rot_diff_hist.txt (verbose>=2) Histogram of the differences in rot;");
    addParamsLine("                             :+ rootname_tilt_diff_hist.txt (verbose>=2) Histogram of the differences in tilt;");
    addParamsLine("                             :+ rootname_psi_diff_hist.txt (verbose>=2) Histogram of the differences in psi;");
    addParamsLine("                             :+ rootname_X_diff_hist.txt (verbose>=2) Histogram of the differences in shiftX;");
    addParamsLine("                             :+ rootname_Y_diff_hist.txt (verbose>=2) Histogram of the differences in shiftY;");
    addParamsLine("  [--sym <symmetry=\"\">]    : Symmetry file if any");
    addParamsLine("                             :+The definition of the symmetry is described at [[transform_symmetrize_v3][transform_symmetrize]]");
    addParamsLine("  [--check_mirrors]          : Check if mirrored projections give better results");
    addParamsLine("  [--object_rotation]        : Use object rotations rather than projection directions");
    addParamsLine("                             : fit (Spider, APMQ)");
}

// Produce side information ================================================
void ProgAngularDistance::produce_side_info()
{
    if (fn_sym != "")
        SL.readSymmetryFile(fn_sym);

    // Check that both docfiles are of the same length
    if (fn_ang1!="" && fn_ang2!="")
    {
        DF1.read(fn_ang1);
        DF2.read(fn_ang2);
        if (DF1.size() != DF2.size())
            REPORT_ERROR(ERR_MD_OBJECTNUMBER,
                         "Angular_distance: Input Docfiles with different number of entries");
    }
}

//#define DEBUG
// Compute distance --------------------------------------------------------
void ProgAngularDistance::run()
{
    produce_side_info();

    MetaData DF_out;
    double angular_distance=0;
    double shift_distance=0;

    MultidimArray<double> rot_diff, tilt_diff, psi_diff, vec_diff,
    X_diff, Y_diff, shift_diff;
    rot_diff.resize(DF1.size());
    tilt_diff.resize(rot_diff);
    psi_diff.resize(rot_diff);
    vec_diff.resize(rot_diff);
    X_diff.resize(rot_diff);
    Y_diff.resize(rot_diff);
    shift_diff.resize(rot_diff);

    // Build output comment
    DF_out.setComment("image rot1 rot2 diff_rot tilt1 tilt2 diff_tilt psi1 psi2 diff_psi ang_dist X1 X2 Xdiff Y1 Y2 Ydiff ShiftDiff");

    int i = 0;
    size_t id;
    FileName fnImg;
    std::vector<double> output;
    output.resize(17,0);
    bool fillOutput=fn_out!="";
    FOR_ALL_OBJECTS_IN_METADATA2(DF1, DF2)
    {
        // Read input data
        double rot1,  tilt1,  psi1;
        double rot2,  tilt2,  psi2;
        double rot2p, tilt2p, psi2p;
        double distp;
        double X1, X2, Y1, Y2;
        DF1.getValue(MDL_IMAGE,fnImg,__iter.objId);

        DF1.getValue(MDL_ANGLE_ROT,rot1,__iter.objId);
        DF1.getValue(MDL_ANGLE_TILT,tilt1,__iter.objId);
        DF1.getValue(MDL_ANGLE_PSI,psi1,__iter.objId);
        DF1.getValue(MDL_SHITF_X,X1,__iter.objId);
        DF1.getValue(MDL_SHITF_Y,Y1,__iter.objId);

        DF2.getValue(MDL_ANGLE_ROT,rot2,__iter2.objId);
        DF2.getValue(MDL_ANGLE_TILT,tilt2,__iter2.objId);
        DF2.getValue(MDL_ANGLE_PSI,psi2,__iter2.objId);
        DF2.getValue(MDL_SHITF_X,X2,__iter2.objId);
        DF2.getValue(MDL_SHITF_Y,Y2,__iter2.objId);

        // Bring both angles to a normalized set
        rot1 = realWRAP(rot1, -180, 180);
        tilt1 = realWRAP(tilt1, -180, 180);
        psi1 = realWRAP(psi1, -180, 180);

        rot2 = realWRAP(rot2, -180, 180);
        tilt2 = realWRAP(tilt2, -180, 180);
        psi2 = realWRAP(psi2, -180, 180);

        // Apply rotations to find the minimum distance angles
        rot2p = rot2;
        tilt2p = tilt2;
        psi2p = psi2;
        distp = SL.computeDistance(rot1, tilt1, psi1, rot2p, tilt2p, psi2p, false,
                                   check_mirrors, object_rotation);
        angular_distance += distp;

        // Compute angular difference
        rot_diff(i) = rot1 - rot2p;
        tilt_diff(i) = tilt1 - tilt2p;
        psi_diff(i) = psi1 - psi2p;
        vec_diff(i) = distp;
        X_diff(i) = X1 - X2;
        Y_diff(i) = Y1 - Y2;
        shift_diff(i) = sqrt(X_diff(i)*X_diff(i)+Y_diff(i)*Y_diff(i));
        shift_distance += shift_diff(i);

        // Fill the output result
        if (fillOutput)
        {
            output[0]=rot1;
            output[1]=rot2p;
            output[2]=rot_diff(i);
            output[3]=tilt1;
            output[4]=tilt2p;
            output[5]=tilt_diff(i);
            output[6]=psi1;
            output[7]=psi2p;
            output[8]=psi_diff(i);
            output[9]=distp;
            output[10]=X1;
            output[11]=X2;
            output[12]=X_diff(i);
            output[13]=Y1;
            output[14]=Y2;
            output[15]=Y_diff(i);
            output[16]=shift_diff(i);

            id = DF_out.addObject();
            DF_out.setValue(MDL_IMAGE,fnImg,id);
            DF_out.setValue(MDL_ANGLE_COMPARISON,output, id);
        }

        i++;
    }
    angular_distance /= i;
    shift_distance /=i;

    if (fillOutput)
    {
        DF_out.write(fn_out + ".xmd");
        Histogram1D hist;
        compute_hist(vec_diff, hist, 0, 180, 180);
        hist.write(fn_out + "_vec_diff_hist.txt");
        compute_hist(shift_diff, hist, 20);
        hist.write(fn_out + "_shift_diff_hist.txt");
        if (verbose==2)
        {
            compute_hist(rot_diff, hist, 100);
            hist.write(fn_out + "_rot_diff_hist.txt");
            compute_hist(tilt_diff, hist, 100);
            hist.write(fn_out + "_tilt_diff_hist.txt");
            compute_hist(psi_diff, hist, 100);
            hist.write(fn_out + "_psi_diff_hist.txt");
            compute_hist(X_diff, hist, 20);
            hist.write(fn_out + "_X_diff_hist.txt");
            compute_hist(Y_diff, hist, 20);
            hist.write(fn_out + "_Y_diff_hist.txt");
        }
    }

    std::cout << "Global angular distance = " << angular_distance << std::endl;
    std::cout << "Global shift   distance = " << shift_distance   << std::endl;
}
