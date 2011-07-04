/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2002)
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
#ifndef _PROG_ANGULAR_GCAR
#define _PROG_ANGULAR_GCAR

#include <data/program.h>
#include <data/matrix2d.h>
#include <data/matrix1d.h>
#include <data/funcs.h>
#include <math.h>
#include <vector>
#include <blas1c.h>
#include <lapackc.h>
#include <arsnsym.h>

/**@defgroup AngularGCAR Globally Convergent Angular Reconstitution
   @ingroup ReconsLibrary */
//@{
/** Angular GCAR parameters. */
class ProgAngularGCAR: public XmippProgram
{
public:
public:
public:
    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Usage
    void defineParams();

    /** Run */
    void run();

    ///  qrand
    void qrand(Matrix2D<double>& q, int K);

	/// qToRot
	void qToRot(const Matrix1D<double>& q, Matrix2D<double> &rotMatrix);

	/// commonLineQ
	void commonLineQ(const Matrix2D<double>& q, int k1, int k2, int nTheta, int& idx1, int& idx2);

	/// clmatrixCheaatQ
	void clmatrixCheatQ(const Matrix2D<double>& q, int nTheta, Matrix2D<int>& clmatrix, Matrix2D<int>&  clcorr);

	void cryoCosmetify(const Matrix2D<double>& PHI, Matrix2D<double>& PHI2, int K, int L);

	void Q2S2(const Matrix2D<double>& Q, Matrix2D<double>& PHIRef, int NTheta);

	void cryoOrientationsSdp(const SparseMatrix2D W, int nEigs, Matrix2D<double> PHI);
};

class EigElement {
public:
	double eigenvalue;
	int    pos;
};

inline bool operator<(const EigElement& e1, const EigElement& e2)
{
	return e1.eigenvalue>e2.eigenvalue; //Orden descendente
}

//@}
#endif
