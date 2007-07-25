/***************************************************************************
 *
 * Authors: Sjors H.W. Scheres (scheres@cnb.uam.es)
 *
 *  This code is strongly based on ideas by Pawel Penczek & Zhengfan
 *  Yang as implemented in SPARX at the University of Texas - Houston 
 *  Medical School
 *
 *  see P. A. Penczek, R. Renka, and H. Schomberg,
 *      J. Opt. Soc. Am. _21_, 449 (2004)
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
#ifndef POLAR_H
#define POLAR_H
#include "funcs.h"
#include "matrix2d.h"
#include "gridding.h"

#define HALF_CIRCLES 0
#define FULL_CIRCLES 1

/** Produces a vector of integers with the information for 2D polar
 *  coordinates generation.
 *  For each ring, this vector contains 3 numbers:
 *  [i*3]:   Radius of the ring
 *  [i*3+1]: Total number of samples of all inner rings
 *  [i*3+2]: Number of samples of this ring.
 */
vector<int> getPolarRingInfo2D(int first_ring, int last_ring, 
			       int step = 1, int mode = FULL_CIRCLES,
			       bool const_sam = false);

/** Produces a vector of doubles with the weights for the 2D polar
 *  coordinates. 
 *  
 *  These weights are ~ 1/radius
 */
vector<double> getPolarWeights2D(vector<int> ring_info, int mode = FULL_CIRCLES);

/** Converts input matrix2D M1 to an output vector, which is in polar
 *  coordinates, according to the ring structur as stored in ring_info.
 * 
 *  Note that the shifts options are only valid for real-space matrices.
 */
template<typename T>
vector<Matrix1D<T> > polarCoordinatesGridding(const Matrix2D<T> &M1, vector<int> ring_info,
					      KaiserBessel &kb, int mode = FULL_CIRCLES,
					      double xoff = 0., double yoff = 0.)
{
    int nsam;
    int nring = ring_info.size()/3;
    double radius, twopi, dphi, phi; 
    double xp, yp;
    vector<Matrix1D<T> > out;
    Matrix1D<T> Mring;
    out.clear();

    // Check no pixel falls outside the image
    if (ABS(xoff) + ring_info[3 * (nring - 1)] > XSIZE(M1)/4 || 
	ABS(yoff) + ring_info[3 * (nring - 1)] > YSIZE(M1)/4)
	REPORT_ERROR(1,"polarCoordinatesGridding: last ring falls outside image");

    if (mode == FULL_CIRCLES)
	twopi = 2*PI;
    else if (mode == HALF_CIRCLES)
	twopi = PI;
    else
	REPORT_ERROR(1,"Incorrect mode for getPolarWeights2D!");

    for (int iring = 0; iring < nring; iring++)
    {
	radius = (double) ring_info[3 * iring];
	nsam = ring_info[3 * iring + 2];
	dphi = twopi / (double)nsam;
	Mring.resize(nsam);
	for (int iphi = 0; iphi < nsam; iphi++)
	{
	    // polar coordinates
	    phi = iphi * dphi;
	    xp = sin(phi) * radius;
	    yp = cos(phi) * radius;
	    xp += xoff;
	    yp += yoff; 
	    Mring(iphi) = (T) interpolatedElementGridding(M1,xp,yp,kb);
	}
	out.push_back(Mring);
    }
    return out;

}

vector<Matrix1D<complex<double> > > fourierTransformRings(
    const vector<Matrix1D<double> > &in, bool conjugated = false);

vector<Matrix1D<double> > inverseFourierTransformRings(
    const vector<Matrix1D<complex<double> > > &in);

/** This function find the optimal rotational angle between images img
 *  and ref. 

 * It is assumed that ref has the corresponding weights applied and
 * that both images are Fourier-transforms of the rings of the images
 * in polar coordinates.
 */
double getBestAnglePolarGridding(const vector<Matrix1D<complex<double> > > &img,
				 const vector<Matrix1D<complex<double> > > &ref,
				 double &maxcorr, bool interpolate = false);
				 


#endif
