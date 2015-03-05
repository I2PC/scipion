/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
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

#ifndef SPHERICAL_HARMONICS_H_
#define SPHERICAL_HARMONICS_H_

#include <alglib/src/specialfunctions.h>

#include "multidim_array.h"
#include "matrix2d.h"
//#include "../../../kk/gsl-1.16/gsl/gsl_sf_gegenbauer.h"

# define pi           3.14159265358979323846  /* pi */
# define sampling  100;

#endif /* SPHERICAL_HARMONICS_H_ */

//Spherical harmonics (SH) as defined by Guan Koay in "A signal transformational
//framework for breaking the noise floor and its applications in MRI".
//he obtained Spherical harmonics are real and orthonormal. Theta is zenith angle and varphi is azimuthal angle, as defined in
// physics space
class PolySphericalHarmonics
{
private:

	bool save;

public:

	void shPols(int degree, MultidimArray<double> & im, int verbose=0);

	inline double shPols(int degree, int m, double theta,double varphi);

private:


	//We store it as double because there is a lot of division between factorials in the
	// SH implementation
	inline  double factorial(size_t n)
	{
		double f=1;

		for (size_t var = 1; var <= n; var++)
			f *= var;

		return f;
	}

public:

	PolySphericalHarmonics(bool doSave=false)
	{
			save = doSave;
	}

};




