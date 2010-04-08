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

#include "matrix2d.h"

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
	/// Object plane on Focus
	double Zo;
	/// Image plane (CCD position)
	double Zi;
	/// Minimum resolution condition
	double dxiMax;

public:
    /// Lambda
    double lambda;

    /* RX Microscope configuration */

    /// Focal length
    double Flens;
    /// Number of zones in zone plate
    double Nzp;
    /// Magnification
    double Ms;

    /// object space XY-plane sampling rate (nm/pixel)
    double dxo;
    /// object space Z sampling rate (nm/pixel)
    double dzo;
    /// Shift

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
    friend std::ostream & operator << (std::ostream &out, const XmippXRPSF &psf);

    /// Clear.
    void clear();

    /// Produce Side information
    void produceSideInfo();

    /// Apply OTF to an image
    void applyOTF(Matrix2D < double > &I) const;

    /** Generate OTF image.
        The sample image is used only to take its dimensions. */
    template <class T>
    void generateOTF(const Matrix2D<T> &sample_image) const
    {
    	generateOTF(YSIZE(sample_image), XSIZE(sample_image));
        STARTINGX(OTF) = STARTINGX(sample_image);
        STARTINGY(OTF) = STARTINGY(sample_image);
    }

    /// Generate OTF image.
    void generateOTF(int Ydim, int Xdim) const;
};
//@}
#endif
