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

#include "projection.h"
#include "geometry.h"

/* Reset =================================================================== */
void Projection::reset(int Ydim, int Xdim)
{
    img.init_zeros(Ydim,Xdim);
    move_origin_to_center();
}

/* Set angles ============================================================== */
void Projection::set_angles(double _rot, double _tilt, double _psi)
{
    set_eulerAngles(_rot,_tilt,_psi);
    Euler_angles2matrix(_rot,_tilt,_psi,euler);
    eulert=euler.transpose();
    euler.getRow(2,direction);
    direction.self_transpose();
}

/* Read ==================================================================== */
void Projection::read(const FileName &fn, const bool &apply_shifts)
{
    ImageXmipp::read(fn,false,false,false,apply_shifts);
    Euler_angles2matrix(rot(),tilt(),psi(),euler);
    eulert=euler.transpose();
    euler.getRow(2,direction);
    direction.self_transpose();
}

/* Assignment ============================================================== */
Projection & Projection::operator = (const Projection &P)
{
    // Esto hay que ponerlo m�s elegantemente accediendo al = del padre
    *(ImageXmipp *)this = * ((ImageXmipp *) &P);
    direction = P.direction;
    euler     = P.euler;
    eulert    = P.eulert;
    return *this;
}

/* Another function for assignment ========================================= */
void Projection::assign (const Projection &P)
{
    *this=P;
}

