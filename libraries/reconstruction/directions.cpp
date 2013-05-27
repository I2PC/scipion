/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2004)
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

#include "directions.h"

// Check whether projection directions are unique =================================
bool directions_are_unique(double rot,  double tilt,
                           double rot2, double tilt2,
                           double rot_limit, double tilt_limit,
                           SymList &SL, bool include_mirrors,
                           Matrix2D<double> &Laux, Matrix2D<double> &Raux)
{
    bool are_unique = true;
    double rot2p, tilt2p, psi2p, psi2 = 0.;
    double diff_rot, diff_tilt;

    int isymmax=SL.symsNo();
    for (int isym = 0; isym <= isymmax; isym++)
    {
        if (isym == 0)
        {
            rot2p = rot2;
            tilt2p = tilt2;
            psi2p = psi2;
        }
        else
        {
            SL.getMatrices(isym - 1, Laux, Raux,false);
            Euler_apply_transf(Laux, Raux, rot2, tilt2, psi2, rot2p, tilt2p, psi2p);
        }

        double aux=rot - rot2p;
        diff_rot = fabs(realWRAP(aux, -180, 180));
        aux=tilt - tilt2p;
        diff_tilt = fabs(realWRAP(aux, -180, 180));
        if ((rot_limit - diff_rot) > 1e-3 && (tilt_limit - diff_tilt) > 1e-3)
            are_unique = false;
        Euler_another_set(rot2p, tilt2p, psi2p, rot2p, tilt2p, psi2p);
        aux=rot - rot2p;
        diff_rot = fabs(realWRAP(aux, -180, 180));
        aux=tilt - tilt2p;
        diff_tilt = fabs(realWRAP(aux, -180, 180));
        if ((rot_limit - diff_rot) > 1e-3 && (tilt_limit - diff_tilt) > 1e-3)
            are_unique = false;
        if (!include_mirrors)
        {
            Euler_up_down(rot2p, tilt2p, psi2p, rot2p, tilt2p, psi2p);
            aux=rot - rot2p;
            diff_rot = fabs(realWRAP(aux, -180, 180));
            aux=tilt - tilt2p;
            diff_tilt = fabs(realWRAP(aux, -180, 180));
            if ((rot_limit - diff_rot) > 1e-3 && (tilt_limit - diff_tilt) > 1e-3)
                are_unique = false;
            Euler_another_set(rot2p, tilt2p, psi2p, rot2p, tilt2p, psi2p);
            aux=rot - rot2p;
            diff_rot = fabs(realWRAP(aux, -180, 180));
            aux=tilt - tilt2p;
            diff_tilt = fabs(realWRAP(aux, -180, 180));
            if ((rot_limit - diff_rot) > 1e-3 && (tilt_limit - diff_tilt) > 1e-3)
                are_unique = false;
        }
    }

    return are_unique;
}

double distance_directions(double rot1, double tilt1,
                           double rot2, double tilt2,
                           bool include_mirrors)
{

    double            rot2p, tilt2p, psi2p, dist;
    double            diff_rot, diff_tilt;

    double aux=rot1 - rot2;
    diff_rot = fabs(realWRAP(aux, -180, 180));
    aux=tilt1 - tilt2;
    diff_tilt = fabs(realWRAP(aux, -180, 180));
    dist = sqrt((diff_rot * diff_rot) + (diff_tilt * diff_tilt));

    Euler_another_set(rot2, tilt2, 0., rot2p, tilt2p, psi2p);
    aux=rot1 - rot2p;
    diff_rot = fabs(realWRAP(aux, -180, 180));
    aux=tilt1 - tilt2p;
    diff_tilt = fabs(realWRAP(aux, -180, 180));
    dist = fmin(dist, sqrt((diff_rot * diff_rot) + (diff_tilt * diff_tilt)));

    if (include_mirrors)
    {
        Euler_up_down(rot2, tilt2, 0., rot2p, tilt2p, psi2p);
        aux=rot1 - rot2p;
        diff_rot = fabs(realWRAP(aux, -180, 180));
        aux=tilt1 - tilt2p;
        diff_tilt = fabs(realWRAP(aux, -180, 180));
        dist = fmin(dist, sqrt((diff_rot * diff_rot) + (diff_tilt * diff_tilt)));

        Euler_another_set(rot2p, tilt2p, 0., rot2p, tilt2p, psi2p);
        aux=rot1 - rot2p;
        diff_rot = fabs(realWRAP(aux, -180, 180));
        aux=tilt1 - tilt2p;
        diff_tilt = fabs(realWRAP(aux, -180, 180));
        dist = fmin(dist, sqrt((diff_rot * diff_rot) + (diff_tilt * diff_tilt)));
    }

    return dist;

}
// Fill DF with evenly distributed rot & tilt =================================
void make_even_distribution(std::vector<double> &rotList, std::vector<double> &tiltList,
							double sampling, SymList &SL, bool include_mirror)
{
    int rot_nstep, tilt_nstep = ROUND(180. / sampling) + 1;
    double rot_sam, tilt, tilt_sam;
    bool append;
    tilt_sam = (180. / tilt_nstep);

    // Create evenly distributed angles
    rotList.clear();
    tiltList.clear();
    rotList.reserve(20000);  // Normally there are many less directions than 20000
    tiltList.reserve(20000); // Set to 20000 to avoid resizing
    Matrix2D<double> L(3,3),R(3,3);
    for (int tilt_step = 0; tilt_step < tilt_nstep; tilt_step++)
    {
        tilt = ((double)tilt_step / (tilt_nstep - 1)) * 180.;
        if (tilt > 0)
            rot_nstep = CEIL(360. * sin(DEG2RAD(tilt)) / sampling);
        else
            rot_nstep = 1;
        rot_sam = 360. / (double)rot_nstep;
        for (double rot = 0.; rot < 360.; rot += rot_sam)
        {
            // Check whether by symmetry or mirror the angle has been included already
            append = true;
            size_t imax=rotList.size();
            double *ptrRot=NULL;
            double *ptrTilt=NULL;
            if (imax>0)
            {
            	ptrRot=&rotList[0];
            	ptrTilt=&tiltList[0];
            }
            for (size_t i=0; i<imax; ++i, ++ptrRot, ++ptrTilt)
            {
                if (!directions_are_unique(rot, tilt, *ptrRot, *ptrTilt, rot_sam, tilt_sam, SL,
                                           include_mirror, L, R))
                {
                    append = false;
                    break;
                }
            }
            if (append)
            {
            	rotList.push_back(rot);
            	tiltList.push_back(tilt);
            }
        }
    }
}

void limit_tilt_range(MetaData &DF, double tilt_range0, double tilt_rangeF)
{

    MetaData DFaux;
    DFaux.importObjects(DF, MDValueRange(MDL_ANGLE_TILT, tilt_range0, tilt_rangeF));
    DF = DFaux;
}


int find_nearest_direction(double rot1, double tilt1,
                           std::vector<double> &rotList, std::vector<double> &tiltList,
                           SymList &SL, Matrix2D<double> &Laux, Matrix2D<double> &Raux)
{

    int               optdir;
    double            dist, mindist;
    double            newrot, newtilt, newpsi;

    optdir = -1;
    mindist = 9999.;
    int imax=SL.symsNo();
    size_t nmax=rotList.size();
    double *ptrRot=NULL;
    double *ptrTilt=NULL;
    if (nmax>0)
    {
    	ptrRot=&rotList[0];
    	ptrTilt=&tiltList[0];
    }
    for (size_t n=0; n<nmax; ++n, ++ptrRot, ++ptrTilt)
    {
        dist = distance_directions(rot1, tilt1, *ptrRot, *ptrTilt, false);
        if (dist < mindist)
        {
            mindist = dist;
            optdir = n;
        }

        for (int i = 0; i < imax; i++)
        {
            SL.getMatrices(i, Laux, Raux, false);
            Euler_apply_transf(Laux, Raux, rot1, tilt1, 0., newrot, newtilt, newpsi);
            dist = distance_directions(newrot, newtilt, *ptrRot, *ptrTilt, false);
            if (dist < mindist)
            {
                mindist = dist;
                optdir = n;
            }
        }
    }

    return optdir;
}

