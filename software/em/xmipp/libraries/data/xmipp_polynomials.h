/***************************************************************************
 * Authors:     Javier Vargas (jvargas@cnb.csic.es)
 *
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

#ifndef POLYNOMIALS_H_
#define POLYNOMIALS_H_

#include "multidim_array.h"
#include "matrix2d.h"

// This class performs all the work related with Polynomials. This is a base class
class Polynomials
{

#define COEFFICIENTS(poly) (poly.fittedCoeffs)

public :
        // Destructor
        virtual ~Polynomials() {}

	//Fitted coefficients
	Matrix1D<double> fittedCoeffs;

public:
	// fitting a surface given by member im using the polynomials
	virtual void fit(const Matrix1D<int> & coef, MultidimArray<double> & im, MultidimArray<double> &weight,
			         MultidimArray<bool> & ROI, int verbose=0)=0;

protected:
	// Create the polynomials
	virtual void create(const Matrix1D<int> & coef)=0;

private:

};

//The class PolyZernikes heritage from Polynomials
//We follow the Zernike implementation explained in: Efficient Cartesian representation of Zernike polynomials in computer memory
//SPIE Vol. 3190 pp. 382
class PolyZernikes: public Polynomials
{

private:
	std::vector<Matrix2D<int> > fMatV;

public:
	//Create not really the polynomials, This function creates a set of coefficient matrix that are efficient
	// to be stored in memory that give us the analytical expression of the polynomials
	//NOTE: take a look to:
	//"Efficient Cartesian representation of Zernike polynomials in computer memory
	//SPIE Vol. 3190 pp. 382 to get more details about the implementation"
	void create(const Matrix1D<int> & coef);
	//This function obtains the Zernike coefficients from the matrix 1D coeff. This array is formed by zeros and ones.
	// If the value of one element of coef is
	void fit(const Matrix1D<int> & coef, MultidimArray<double> & im, MultidimArray<double> &weight,
			 MultidimArray<bool> & ROI, int verbose=0);
	//Gives the zernike polynomio of index zerIndex that it is stored in im
	void zernikePols(const Matrix1D<int> coef, MultidimArray<double> & im, MultidimArray<bool> & ROI, int verbose=0);
};

#endif /* POLYNOMIALS_H_ */
