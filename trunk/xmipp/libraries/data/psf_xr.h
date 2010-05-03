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


/// This enum defines which method should be used to
/// correct the limitation due to Nyquist limit in diffraction
enum psfxrAdjust
{
	PSFXR_STD, /// Standard mode, nothing changes
	PSFXR_INT, /// Increasing the image size by Interpolating
	PSFXR_ZPAD /// Increasing the image size by Zeropadding
};


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



    /// Size of the input image (object plane size)
    double Nox, Noy;
    /* Minimum resolution condition.
     The same for both axis x-y, due to the simmetry of the lens aperture */
    double dxiMax;
    /// Pixel size in lens plane
    double dxl, dyl;

    /// Parameters to change image size to avoid Nyquist limit
    psfxrAdjust AdjustType;
    double npMin;

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
    /// Size of the image in image plane, to be rescaled if needed
    double Nix, Niy;



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

//#define DEBUG
#ifdef DEBUG

        ImageXmipp _Im;
        _Im().resize(Im);
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Im)
        dMij(_Im(),i,j) = abs(dMij(Im,i,j));

        _Im.write(("psfxr-Imin.spi"));
#endif

        transformer.FourierTransform(Im, ImFT, false);

#ifdef DEBUG

        _Im().resize(ImFT);
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(ImFT)
        dMij(_Im(),i,j) = abs(dMij(ImFT,i,j));
        _Im.write(("psfxr-imft1.spi"));
#endif

        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(ImFT)
        dMij(ImFT,i,j) *= dMij(OTF,i,j);

#ifdef DEBUG

        _Im().resize(ImFT);
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(ImFT)
        dMij(_Im(),i,j) = abs(dMij(ImFT,i,j));
        _Im.write(("psfxr-imft2.spi"));
#endif

        transformer.inverseFourierTransform();

        //        CenterOriginFFT(Im, 1);

#ifdef DEBUG

        _Im().resize(Im);
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Im)
        dMij(_Im(),i,j) = abs(dMij(Im,i,j));
        _Im.write(("psfxr-imout.spi"));
#endif

    }

    /// Generate OTF image.
    void generateOTF(Matrix2D<double> &Im) ;

    void generateOTF(Matrix2D<std::complex<double> > &Im) ;

    /// Calculate if a rescaling of the images is needed due to the limitations of the system parameters
    void adjustParam(ImageXmipp &Im) ;

    void adjustParam(VolumeXmipp &Vol) ;

};


/// Generate the quadratic phase distribution of a ideal lens
void lensPD(Matrix2D<std::complex<double> > &Im, double Flens, double lambda, double dx, double dy);


/// Generate projection for an X-ray microscope... TBC
void project_xr(XmippXRPSF &psf, VolumeXmipp &vol, ImageXmipp &imOut);



//@}
#endif
