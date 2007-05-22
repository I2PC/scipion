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

/* ------------------------------------------------------------------------- */
/* INRIA Library Xmipp front end                                             */
/* ------------------------------------------------------------------------- */
#ifndef _XMIPP_INRIA_HH
#   define _XMIPP_INRIA_HH

#ifdef _HAVE_INRIA
#include <data/matrix3d.h>
#include <data/vectorial.h>

/**@name INRIA Library */
//@{
/**@name Derivatives and gradients */
//@{
/** Compute derivative along a certain direction.
    Valid derivative types are "X", "Y", "Z", "XX", "XY", ... "ZZZ".
    Values in the output volume are the corresponding derivatives at
    each volume position*/
void compute_derivative(const matrix3D<double> &in_vol,
                        matrix3D<double> &out_vol, char *type, double sigma = 1);

/** Gradient of a volume */
void compute_gradient(const matrix3D<double> &in_vol,
                      Vectorial_matrix3D &out_vol, double sigma = 1);
//@}

//@}
#endif
#endif
