/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.uam.es (2004)
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

#include "directions.h"

// Check whether projection directions are unique =================================
bool directions_are_unique(double rot,  double tilt,
                           double rot2, double tilt2,
                           double rot_limit, double tilt_limit,
                           SymList &SL, bool include_mirrors)
{

    bool are_unique = true;
    double rot2p, tilt2p, psi2p, psi2 = 0.;
    double diff_rot, diff_tilt;
    matrix2D<double>  L(4, 4), R(4, 4);

    for (int isym = 0; isym <= SL.SymsNo(); isym++)
    {

        if (isym == 0)
        {
            rot2p = rot2;
            tilt2p = tilt2;
            psi2p = psi2;
        }
        else
        {
            SL.get_matrices(isym - 1, L, R);
            L.resize(3, 3); // Erase last row and column
            R.resize(3, 3); // as only the relative orientation
            // is useful and not the translation
            Euler_apply_transf(L, R, rot2, tilt2, psi2, rot2p, tilt2p, psi2p);
        }

        diff_rot = ABS(realWRAP(rot - rot2p, -180, 180));
        diff_tilt = ABS(realWRAP(tilt - tilt2p, -180, 180));
        if ((rot_limit - diff_rot) > 1e-3 && (tilt_limit - diff_tilt) > 1e-3) are_unique = false;
        Euler_another_set(rot2p, tilt2p, psi2p, rot2p, tilt2p, psi2p);
        diff_rot = ABS(realWRAP(rot - rot2p, -180, 180));
        diff_tilt = ABS(realWRAP(tilt - tilt2p, -180, 180));
        if ((rot_limit - diff_rot) > 1e-3 && (tilt_limit - diff_tilt) > 1e-3) are_unique = false;
        if (!include_mirrors)
        {
            Euler_up_down(rot2p, tilt2p, psi2p, rot2p, tilt2p, psi2p);
            diff_rot = ABS(realWRAP(rot - rot2p, -180, 180));
            diff_tilt = ABS(realWRAP(tilt - tilt2p, -180, 180));
            if ((rot_limit - diff_rot) > 1e-3 && (tilt_limit - diff_tilt) > 1e-3) are_unique = false;
            Euler_another_set(rot2p, tilt2p, psi2p, rot2p, tilt2p, psi2p);
            diff_rot = ABS(realWRAP(rot - rot2p, -180, 180));
            diff_tilt = ABS(realWRAP(tilt - tilt2p, -180, 180));
            if ((rot_limit - diff_rot) > 1e-3 && (tilt_limit - diff_tilt) > 1e-3) are_unique = false;
        }
    }

    return are_unique;

}

double distance_directions(double rot1, double tilt1,
                           double rot2, double tilt2,
                           bool include_mirrors)
{

    double            rot2p, tilt2p, psi2p, dist, mindist;
    double            diff_rot, diff_tilt;

    diff_rot = ABS(realWRAP(rot1 - rot2, -180, 180));
    diff_tilt = ABS(realWRAP(tilt1 - tilt2, -180, 180));
    dist = sqrt((diff_rot * diff_rot) + (diff_tilt * diff_tilt));

    Euler_another_set(rot2, tilt2, 0., rot2p, tilt2p, psi2p);
    diff_rot = ABS(realWRAP(rot1 - rot2p, -180, 180));
    diff_tilt = ABS(realWRAP(tilt1 - tilt2p, -180, 180));
    dist = MIN(dist, sqrt((diff_rot * diff_rot) + (diff_tilt * diff_tilt)));

    if (include_mirrors)
    {
        Euler_up_down(rot2, tilt2, 0., rot2p, tilt2p, psi2p);
        diff_rot = ABS(realWRAP(rot1 - rot2p, -180, 180));
        diff_tilt = ABS(realWRAP(tilt1 - tilt2p, -180, 180));
        dist = MIN(dist, sqrt((diff_rot * diff_rot) + (diff_tilt * diff_tilt)));

        Euler_another_set(rot2p, tilt2p, 0., rot2p, tilt2p, psi2p);
        diff_rot = ABS(realWRAP(rot1 - rot2p, -180, 180));
        diff_tilt = ABS(realWRAP(tilt1 - tilt2p, -180, 180));
        dist = MIN(dist, sqrt((diff_rot * diff_rot) + (diff_tilt * diff_tilt)));
    }

    return dist;

}
// Fill DF with evenly distributed rot & tilt =================================
void make_even_distribution(DocFile &DF, double sampling,
                            SymList &SL, bool include_mirror)
{

    int rot_nstep, tilt_nstep = ROUND(180. / sampling) + 1;
    double rotp, tiltp, psip, rot_sam, tilt, rot, tilt_sam, psi = 0.;
    bool append;
    Matrix1D<double> dataline(3);
    tilt_sam = (180. / tilt_nstep);

    DF.clear();
    // Create evenly distributed angles
    for (int tilt_step = 0; tilt_step < tilt_nstep; tilt_step++)
    {
        tilt = ((double)tilt_step / (tilt_nstep - 1)) * 180.;
        if (tilt > 0) rot_nstep = CEIL(360. * sin(DEG2RAD(tilt)) / sampling);
        else rot_nstep = 1;
        rot_sam = 360. / (double)rot_nstep;
        for (double rot = 0.; rot < 360.; rot += rot_sam)
        {
            // Check whether by symmetry or mirror the angle has been included already
            append = true;
            DF.go_first_data_line();
            while (!DF.eof())
            {
                if (!directions_are_unique(rot, tilt, DF(0), DF(1), rot_sam, tilt_sam, SL, include_mirror))
                {
                    append = false;
                    break;
                }
                DF.next_data_line();
            }
            if (append)
            {
                dataline(0) = rot;
                dataline(1) = tilt;
                dataline(2) = 0.;
                DF.append_data_line(dataline);
            }
        }
    }


}

void limit_tilt_range(DocFile &DF, double tilt_range0, double tilt_rangeF)
{

    // Select use-provided tilt range
    double tilt;
    DocLine DL;
    DocFile Dt;
    DF.go_first_data_line();
    while (!DF.eof())
    {
        DL = DF.get_current_line();
        tilt = DF(1);
        if (tilt >= tilt_range0 && tilt <= tilt_rangeF)
        {
            Dt.append_line(DL);
        }
        DF.next_data_line();
    }
    DF = Dt;
    Dt.clear();

}


int find_nearest_direction(double rot1, double tilt1,
                           DocFile &DFlib, int col_rot, int col_tilt, SymList &SL)
{

    int               dir, optdir;
    double            dist, mindist;
    double            newrot, newtilt, newpsi;
    matrix2D<double>  L(4, 4), R(4, 4);

    DFlib.go_first_data_line();
    optdir = dir = 1;
    mindist = 9999.;
    while (!DFlib.eof())
    {
        dist = distance_directions(rot1, tilt1, DFlib(col_rot), DFlib(col_tilt), false);
        if (dist < mindist)
        {
            mindist = dist;
            optdir = dir;
        }
        for (int i = 0; i < SL.SymsNo(); i++)
        {
            SL.get_matrices(i, L, R);
            L.resize(3, 3);
            R.resize(3, 3);
            Euler_apply_transf(L, R, rot1, tilt1, 0., newrot, newtilt, newpsi);
            dist = distance_directions(newrot, newtilt, DFlib(col_rot), DFlib(col_tilt), false);
            if (dist < mindist)
            {
                mindist = dist;
                optdir = dir;
            }
        }

        DFlib.next_data_line();
        dir++;
    }

    return optdir;

}

