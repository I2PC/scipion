/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@mipg.upenn.edu)
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
#ifndef _PROG_EULER_HH
#define _PROG_EULER_HH

#include <data/args.h>
#include <data/docfile.h>
#include <data/funcs.h>
#include <data/matrix2d.h>
#include <data/matrix1d.h>
#include <data/projection.h>
#include <data/geometry.h>

/**@defgroup CrystalEuler Euler (Calculate new tilt after volume compresion)
   @ingroup ReconsLibraryPrograms */
//@{

/** Compute the Euler angles in real space. */
void ROUT_EULER(const double rot,
                const double tilt,
                const double psi,
                const double d00,
                const double d01,
                const double d10,
                const double d11);

/** Compute the Euler angles in Fourier space. */
void ROUT_EULER1(const double rot,
                 const double tilt,
                 const double psi,
                 const double d00,
                 const double d01,
                 const double d10,
                 const double d11);
//@}
#endif
