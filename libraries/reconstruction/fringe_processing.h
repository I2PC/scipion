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

#include <data/multidim_array.h>
#include <data/matrix2d.h>
#include <data/histogram.h>
#include <data/xmipp_fftw.h>

/** @defgroup FringeProcessing Routines to work with fringes
 *  @ingroup ReconsLibrary
 *
 * @{
 */

/** Types of fringes.
 * This class groups all the fringe processing methods that will be used to process the CTF
 * in the non parametric CTF estimation approach
 */
enum FP_TYPE { SIMPLY_OPEN_FRINGES, SIMPLY_CLOSED_FRINGES, COMPLEX_OPEN_FRINGES, COMPLEX_CLOSED_FRINGES, SIMPLY_CLOSED_FRINGES_MOD };

/// Function to simulate some test fringe patterns
void simulPattern(MultidimArray<double> & im, enum FP_TYPE type, int xdim, int ydim, double noiseLevel, const double freq, Matrix1D<int> coefs);

/** Function to simulate some test fringe patterns.
 * sd=SPHT(c) computes the quadrture term of c still affected by the
 * direction phase factor. Therefore for a real c=b*cos(phi)
 * sd=SPHT(c)=i*exp(i*dir)*b*sin(phi)
 * Ref: Kieran G. Larkin, Donald J. Bone, and Michael A. Oldfield, "Natural
 * demodulation of two-dimensional fringe patterns. I. General background of the spiral phase quadrature transform," J. Opt. Soc. Am. A 18, 1862-1870 (2001) */
void SPTH(FourierTransformer &ftrans, MultidimArray<double> & im, MultidimArray< std::complex <double> > & imProc);

//Orientation by minimun diference fit. This method computes the orientation using a window size. The default
//value are wSize=5.
//A reference to the current method can be found in:
//Yang, Xia; Yu, Qifeng, and Fu, Sihua. An algorithm for estimating both fringe orientation and fringe density. Optics Communications.
//2007 Jun 15; 274(2):286-292
void orMinDer(const MultidimArray<double> & im, MultidimArray<double > & orMap,  MultidimArray<double > & orModMap, int wSize, MultidimArray<bool > & ROI);

//This function computes the normalized version of the fringe pattern im = a+m*cos(phi) that it is
//imN = cos(phi) and computes also the modulation map m that it is called imModMap; R and S are
//rough estimations of the number of fringes in the image that give us the fringe frequency and S is the variance of the exponential
//that filter the frequency of the fringes
//Ref: Juan Antonio Quiroga, Manuel Servin, "Isotropic n-dimensional fringe
//pattern normalization", Optics Communications, 224, Pages 221-227 (2003)
void normalize(FourierTransformer &ftrans, MultidimArray<double> & im, MultidimArray<double > & imN,  MultidimArray<double > & imModMap,
		double R, double S, MultidimArray<bool> & ROI);

//This method is similar to the method normalize but it uses a different Fourier Filter H.
void normalizeWB(MultidimArray<double> & im, MultidimArray<double > & imN,  MultidimArray<double > & imModMap, double rmax, double rmin, MultidimArray<bool> & ROI);

//This method is similar to the method normalize and normalizeWB2 but it uses a different Fourier Filter H.
void normalizeWB2(MultidimArray<double> & im, MultidimArray<double > & imN,  MultidimArray<double > & imModMap, double rmax, double rmin, MultidimArray<bool> & ROI);

//This method obtains the phase direction map from the fringe orientation map solving the sign ambiguity problem that exists in the fringe orientation map.
//Once computed the phase direction map the modulating phase can be obtained from the SPTH transform.
//REF: Jesus Villa, Ismael De la Rosa, and Gerardo Miramontes, Juan Antonio Quiroga, Phase recovery from a single fringe pattern using an orientational
//vector-field-regularized estimator J. Opt. Soc. Am. A Vol. 22, No. 12, (2005)
void direction(const MultidimArray<double> & orMap, MultidimArray<double> & qualityMap, double lambda, int size, MultidimArray<double> & dirMap, int x, int y);

//This method performs the phase unwrapping process obtaining the absolute phase from a wrapped phase that has 2pi jumps. The method is an improved version
//of the work presented in:
//Miguel A. Navarro, Julio C. Estrada, M. Servin, Juan A. Quiroga, and Javier Vargas,
//"Fast two-dimensional simultaneous phase unwrapping and low-pass filtering," Opt. Express 20, 2556-2561 (2012)
void unwrapping(const MultidimArray<double> & wrappedPhase, MultidimArray<double> & qualityMap, double lambda, int size, MultidimArray<double> & unwrappedPhase);

//This method performs all the process to demodulate a CTF (im) in order to obtain the phase and the envelope that it is called mod. The method calls the rest of
//the methods for this purpose.  R and S are rough estimations of the number of fringes in the image (R) and S is the variance of the exponential that filter the frequency of the fringes.
//lambda and size is the regularization parameter and the window size where the regularization it is performed. In order to see typical values of these
//parameters see the test: test_fringe_processing_main.cpp
//COMMENTS: verbose is an argument to save intermediate maps
//verbose == 1 saves the normalize map of the fringe pattern
//verbose == 2 saves the orientation map of the fringe pattern
//verbose == 3 saves the direction map
//verbose == 4 saves the wrapped phase map
//verbose == 5 saves all
void demodulate(MultidimArray<double> & im, double lambda, int size, int x, int y, int rmin, int rmax,
		MultidimArray<double> & phase, MultidimArray<double> & mod, Matrix1D<double> & coeffs, int verbose=0);

void demodulate2(MultidimArray<double> & im, double lambda, int size, int x, int y, int rmin, int rmax,
				Matrix1D<double> & coeffs, int verbose=0);

void firsPSDZero(MultidimArray<double> & enhancedPSD, Matrix1D<double> & xPoints,Matrix1D<double> & yPoints, double rmin,
		double rmax, int numAngles, int verbose=0);

void fitEllipse(Matrix1D<double> & xPts, Matrix1D<double> & yPts, double & x0, double & y0, double & majorAxis, double & minorAxis,
		double & ellipseAngle);

void fitEllipse(MultidimArray<double> & normImag, double & x0, double & y0, double & majorAxis, double & minorAxis,
		double & ellipseAngle, size_t & area);

void calculateDefocus(double & defocusU,double & defocusV, double majorAxis, double minorAxis,  double Q0, double lambda, double Cs,
						double imgSize, double Ts);

//@}
#endif /* FRINGEPROCESSING_H_ */
