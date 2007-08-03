/***************************************************************************
 *
 * Authors:     Carlos Oscar Sanchez Sorzano (coss@cnb.uam.es)
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
/*****************************************************************************/
/* INTERACTION WITH PDBs                                                     */
/*****************************************************************************/

#ifndef _XMIPP_PDB_HH
#define _XMIPP_PDB_HH

#include <string>
#include <data/matrix1d.h>

/**@name PDB */
//@{
/** Returns the charge of an atom.
    Returns 0 if the atom is not within the short list (H, C, N, O, S, P, Fe)
    of valid atoms. */
int atomCharge(const string &atom);

/** Returns the radius of an atom.
    Returns 0 if the atom is not within the short list (H, C, N, O, S, P, Fe)
    of valid atoms. */
double atomRadius(const string &atom);

/** Compute the center of mass and limits of a PDB file.
    The limits are referred to the center of mass, i.e., the
    extension of the PDB goes from centerOfMass-limit0 to
    centerOfMass+limitF.
*/
void computePDBgeometry(const std::string &fnPDB,
   Matrix1D<double> &centerOfMass,
   Matrix1D<double> &limit0, Matrix1D<double> &limitF);

/** Apply geometry transformation to an input PDB.
    The result is written in the output PDB. Set centerPDB if you
    want to compute the center of mass first and apply the transformation
    after centering the PDB. */
void applyGeometry(const string &fn_in, const string &fn_out,
    const Matrix2D<double> &A, bool centerPDB=true);
//@}
#endif
