/***************************************************************************
 *
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
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

#ifndef _PSF_XR_HH
#define _PSF_XR_HH

#include "fftw.h"
#include "matrix2d.h"
#include "image.h"
#include "volume.h"
#include <complex>


/**@defgroup PSFXRSupport X-Ray PSF support classes
 @ingroup DataLibrary */
//@{
/** X-ray PSF class.


 */
class XmippXRPSF
{
public:
    // Current OTF
    Matrix2D< std::complex<double> > OTF;

    /* RX Microscope configuration */
    /// Lens Aperture Radius
    double Rlens;
    /// Object plane on Focus (Reference)
    double Zo;
    /// Object plane
    double Z;
    /// Image plane (CCD position)
    double Zi;
    /* Minimum resolution condition.
     The same for both axis x-y, due to the simmetry of the lens aperture */
    double dxiMax;

public:
    /// Lambda
    double lambda;

    /* RX Microscope configuration */

    /// Focal length in mm
    double Flens;
    /// Number of zones in zone plate
    double Nzp;
    /// Magnification
    double Ms;
    /// Z axis global shift
    double DeltaZo;

    /// object space XY-plane sampling rate
    double dxo;
    /// Image space XY-plane sampling rate
    double dxi;
    /// object space Z sampling rate
    double dzo;



    /** Empty constructor. */
    XmippXRPSF()
    {
        clear();
    }

    /** Read from file.
        An exception is thrown if the file cannot be open.*/
    void read(const FileName &fn);

    /** Write to file.
        An exception is thrown if the file cannot be open.*/
    void write(const FileName &fn);

    /// Usage
    void usage();

    /// Show
    friend std::ostream & operator <<(std::ostream &out, const XmippXRPSF &psf);

    /// Clear.
    void clear();

    /// Produce Side information
    void produceSideInfo();

    /// Apply OTF to an image
    template <typename T>
    void applyOTF(Matrix2D<T> &Im)
    {
        Matrix2D<std::complex<double> > ImFT;
        XmippFftw transformer;

//        Im.setXmippOrigin();
#define DEBUG
#ifdef DEBUG

        ImageXmipp _Im;
        _Im().resize(Im);
        _Im().setXmippOrigin();
        FOR_ALL_ELEMENTS_IN_MATRIX2D(Im)
        _Im(i,j) = abs(Im(i,j));

        _Im.write(("psfxr-Imin.spi"));
#endif

        transformer.FourierTransform(Im, ImFT, false);

        FOR_ALL_ELEMENTS_IN_MATRIX2D(ImFT)
        ImFT(i,j) *= OTF(i,j);

        transformer.inverseFourierTransform();

#ifdef DEBUG

//        FOR_ALL_ELEMENTS_IN_MATRIX2D(OTF)
//			_Im(i,j) = abs(OTF(i,j));
//        _Im.write(("apply-otf.spi"));

        CenterOriginFFT(ImFT, 1);
        _Im().resize(ImFT);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(ImFT)
			_Im(i,j) = abs(ImFT(i,j));

        _Im.write(("psfxr-imft.spi"));

#endif

        //        CenterOriginFFT(Im, 1);

#ifdef DEBUG
        _Im().resize(Im);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(Im)
        _Im(i,j) = abs(Im(i,j));

        _Im.write(("psfxr-imout.spi"));
#endif
    }

    /// Generate OTF image.
    void generateOTF(Matrix2D<double> &Im) ;

    void generateOTF(Matrix2D<std::complex<double> > &Im) ;

};


/// Generate the quadratic phase distribution of a ideal lens
void lensPD(Matrix2D<std::complex<double> > &Im, double Flens, double lambda, double dx, double dy);


/// Generate projection for an X-ray microscope... TBC
void project_xr(XmippXRPSF &psf, VolumeXmipp &vol, ImageXmipp &imOut);



//@}
#endif
