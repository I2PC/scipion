/***************************************************************************
 *
 * Authors:     Slavica JONIC (slavica.jonic@impmc.jussieu.fr, slavica.jonic@a3.epfl.ch)
 *              Jean-Noël PIOCHE (jnp95@hotmail.com)
 *
 * Biomedical Imaging Group, EPFL (Lausanne, Suisse).
 * Structures des Assemblages Macromoléculaires, IMPMC UMR 7590 (Paris, France).
 * IUT de Reims-Châlons-Charleville (Reims, France).
 *
 * Last modifications by JNP the 27/05/2009 15:52:45
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

#ifndef __PROJECTION_REAL_SHEARS_H__
#define __PROJECTION_REAL_SHEARS_H__

#include <data/projection.h>

/**@name ProjRealShears Projection library program using Real-Shears */
//@{
/// Structure for holding a volume
class RealShearsInfo
{
public:
    const MultidimArray<double> *volume;
    int Xdim;
    Matrix2D<double> Ac, Acinv;
    MultidimArray<double> Coef_x, Coef_y, Coef_z;

    // Constructor
    RealShearsInfo(const MultidimArray<double> &V);
};

/// Make projection
void projectVolume(RealShearsInfo &Data, Projection &P, int Ydim, int Xdim,
                    double rot, double tilt, double psi,
                    double shiftX=0, double shiftY=0);

//@}
#endif
