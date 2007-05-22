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

#include "radon.h"

#include <data/geometry.h>

void Radon_Transform(Volume *vol, double rot, double tilt,
                     matrix1D<double> &RT)
{
    matrix3D<double>   Rvol;

    // Align the Radon direction with the Z axis
    if (rot != 0 || tilt != 0) Euler_rotate(VOLMATRIX(*vol), rot, tilt, 0.0F, Rvol);
    else Rvol = VOLMATRIX(*vol);

    // Project onto one line
    RT.init_zeros(Rvol.zdim);
    STARTINGX(RT) = STARTINGZ(Rvol);

    for (int k = STARTINGZ(Rvol); k < FINISHINGZ(Rvol); k++)
        for (int i = STARTINGY(Rvol); i < FINISHINGY(Rvol); i++)
            for (int j = STARTINGX(Rvol); j < FINISHINGX(Rvol); j++)
                VEC_ELEM(RT, k) += VOL_ELEM(Rvol, k, i, j);
}

// The number of voxels in each layer is a double to make easier some
// computations outside
void Local_Radon_Transform(Volume *vol, double rot, double tilt,
                           int label, Volume *vol_label, matrix1D<double> &RT,
                           matrix1D<double> &RT_n)
{
    matrix3D<double>   Rvol;
    matrix3D<double>   Lvol;

    // Align the Radon direction with the Z axis
    if (rot != 0 || tilt != 0)
    {
        Euler_rotate(VOLMATRIX(*vol), rot, tilt, 0.0F, Rvol);
        Euler_rotate(VOLMATRIX(*vol_label), rot, tilt, 0.0F, Lvol);
    }
    else
    {
        Rvol = VOLMATRIX(*vol);
        Lvol = VOLMATRIX(*vol_label);
    }

    // Project onto one line
    RT.init_zeros(Rvol.zdim);
    STARTINGX(RT) = STARTINGZ(Rvol);
    RT_n = RT;

    for (int k = STARTINGZ(Rvol); k < FINISHINGZ(Rvol); k++)
        for (int i = STARTINGY(Rvol); i < FINISHINGY(Rvol); i++)
            for (int j = STARTINGX(Rvol); j < FINISHINGX(Rvol); j++)
                if (VOL_ELEM(Lvol, k, i, j) == label)
                {
                    VEC_ELEM(RT, k) += VOL_ELEM(Rvol, k, i, j);
                    VEC_ELEM(RT_n, k)++;
                }
}

/* Radon Transform of an image --------------------------------------------- */
void Radon_Transform(const matrix2D<double> &I, double rot_step,
                     matrix2D<double> &RT)
{
    matrix2D<double> rot_I;
    RT.init_zeros(CEIL(360.0 / rot_step), XSIZE(I));
    STARTINGX(RT) = STARTINGX(I);
    int l = 0;
    for (double rot = 0; rot < 360; rot += rot_step, l++)
    {
        // Rotate image
        I.rotate(rot, rot_I);
        // Sum by columns
        FOR_ALL_ELEMENTS_IN_MATRIX2D(rot_I) RT(l, j) += rot_I(i, j);
    }
}
