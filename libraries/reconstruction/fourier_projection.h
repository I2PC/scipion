/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
 * 02111-1307  USAcd ..
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef _FOURIER_FOURIER_PROJECTION_H
#define _FOURIER_FOURIER_PROJECTION_H

#include <data/xmipp_image.h>
#include <data/xmipp_program.h>
#include <data/filters.h>
#include <data/xmipp_fftw.h>
#include <data/projection.h>

/**@defgroup FourierProjection Fourier projection
   @ingroup ReconsLibrary */
//@{

/** Program class to create projections in Fourier space */
class FourierProjector
{
public:
    /// Padding factor
    double paddingFactor;
    /// Maximum Frequency for pixels
    double maxFrequency;
    /// The order of B-Spline for interpolation
    double BSplineDeg;

public:
    // Auxiliary FFT transformer
    FourierTransformer transformer3D, transformer2D;

    // Volume to project
    MultidimArray<double> *volume;

    // Real and imaginary B-spline coefficients for Fourier of the volume
    MultidimArray< double > VfourierRealCoefs, VfourierImagCoefs;

    // Projection in Fourier space
    MultidimArray< std::complex<double> > projectionFourier;

    // Projection in real space
    Image<double> projection;

    // Original volume size
    int volumeSize;

    // Euler matrix
    Matrix2D<double> E;
public:
    /*
     * The constructor of the class
     */
    FourierProjector(MultidimArray<double> &V, double paddFactor, double maxFreq, int BSplinedegree);
    /**
     * This method gets the volume's Fourier and the Euler's angles as the inputs and interpolates the related projection
     */
    void project(double rot, double tilt, double psi);
private:
    /*
     * This is a private method which provides the values for the class variable
     */
    void produceSideInfo();
};

/*
 * This function gets an object form the FourierProjection class and makes the desired projection in Fourier space
 */
void projectVolume(FourierProjector &projection, Projection &P, int Ydim, int Xdim,
                   double rot, double tilt, double psi);
//@}

#endif
