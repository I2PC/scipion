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
protected :

	//Number of polynomials
	int numPol;
	// Auxiliary and temporal data to store the generated polynomials
	MultidimArray<double> imPols;
	// Pointer to the image to be fitted
	MultidimArray<double> * im;

public:
	// fitting a surface given by member im by the polynomials
	virtual void fit(const Matrix1D<int> & coef, MultidimArray<double> & im)=0;

protected:

	virtual void create(const Matrix1D<int> & coef)=0;

private:

};

//The class PolyZernikes heritage from Polynomials
//We follow the Zernike implementation explained in: Efficient Cartesian representation of Zernike polynomials in computer memory
//SPIE Vol. 3190 pp. 382
class PolyZernikes: public Polynomials
{
private:
	std::vector<Matrix2D< int> > fMatV;

public:
	void create(const Matrix1D<int> & coef);
	void fit(const Matrix1D<int> & coef, MultidimArray<double> & im);
};

#endif /* POLYNOMIALS_H_ */
