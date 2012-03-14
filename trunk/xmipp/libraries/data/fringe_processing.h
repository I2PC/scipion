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

#ifndef FRINGEPROCESSING_H_
#define FRINGEPROCESSING_H_

#include "multidim_array.h"
#include "matrix2d.h"

//This class groups all the fringe processing methods that will be used to process the CTF
//in the non parametric CTF estimation approach
class FringeProcessing
{

public:

	enum FP_TYPE { SIMPLY_OPEN_FRINGES, SIMPLY_CLOSED_FRINGES, COMPLEX_OPEN_FRINGES, COMPLEX_CLOSED_FRINGES };

	//Function to simulate some test fringe patterns
    void simulPattern(MultidimArray<double> & im, enum FP_TYPE type, int xdim, int ydim, double noiseLevel, const double freq, Matrix1D<int> coefs);

	//Function to simulate some test fringe patterns.
	//sd=SPHT(c) computes the quadrture term of c still affected by the
	//direction phase factor. Therefore for a real c=b*cos(phi)
	//sd=SPHT(c)=i*exp(i*dir)*b*sin(phi)
	//Ref: Kieran G. Larkin, Donald J. Bone, and Michael A. Oldfield, "Natural
	//demodulation of two-dimensional fringe patterns. I. General background of the spiral phase quadrature transform," J. Opt. Soc. Am. A 18, 1862-1870 (2001)
    void SPTH(MultidimArray<double> & im, MultidimArray< std::complex <double> > & imProc);

    //Orientation by minimun diference fit. This method computes the orientation using a window size. The default
    //value are wSize=5.
    //A reference to the current method can be found in:
    //Yang, Xia; Yu, Qifeng, and Fu, Sihua. An algorithm for estimating both fringe orientation and fringe density. Optics Communications.
    //2007 Jun 15; 274(2):286-292
    void orMinDer(const MultidimArray<double> & im, MultidimArray<double > & orMap,  MultidimArray<double > & orModMap, int wSize);

    //This function computes the normalized version of the fringe pattern im = a+m*cos(phi) that it is
    //imN = cos(phi) and computes also the modulation map m that it is called imModMap; fmin and fmax are
    //rough estimations of the maximun and minimun frequency of the fringes and num is the number of
    //Gabor filters used.
    //Ref: Juan Antonio Quiroga, Manuel Servin, "Isotropic n-dimensional fringe
    //pattern normalization", Optics Communications, 224, Pages 221-227 (2003)
    void normalize(MultidimArray<double> & im, MultidimArray<double > & imN,  MultidimArray<double > & imModMap, int fmin, int fmax, int num);

protected:

    void direction(const MultidimArray<double> & im, double lambda, int size);
};


#endif /* FRINGEPROCESSING_H_ */
