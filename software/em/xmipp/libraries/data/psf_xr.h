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

//#include <complex>
#include "xmipp_image.h"
#include "xmipp_fftw.h"
#include "projection.h"
#include "xmipp_program.h"

/**@defgroup PSFXRSupport X-Ray Microscope PSF class
   @ingroup DataLibrary */
//@{

/** This enum defines which method should be used to
 correct the constraint due to Nyquist limit in diffraction. */
enum PsfxrAdjust
{
    PSFXR_STD, /// Standard mode, image size does not changes
    PSFXR_INT, /// Increasing the image size by Interpolating
    PSFXR_ZPAD /// Increasing the image size by Zeropadding
} ;

// PSF Generation algorithm
/* ANALITIC_ZP is based on  O. Mendoza-Yero et als."PSF analysis of nanometric Fresnel
 * zone plates," EOS Topical Meeting on Diffractive Optics 2010,
 * ISBN 978-3-00-024193-2, 14th-18th February
 * 2010, Koli, Finland. (contact email:omendoza@uji.es)
 */
enum PsfType
{
    IDEAL_FRESNEL_LENS,
    ANALYTIC_ZP
};

/** X-ray PSF class.
 * Here goes how to filter an image with the Point Spread Function of a X-ray microscope optics

 @code

 #include <data/psf_xr.h>


 int main(int argc, char **argv)
 {
  FileName fnPSF, fn_input,fn_output;
  XRayPSF psf;

 /// Read the microscope parameters

 if (checkParameter(argc, argv, "-psf"))
 {
  fnPSF = getParameter(argc, argv, "-psf");
  psf.read(fnPSF);
 }
 else
  psf.clear(); /// Use a default reference microscope


 fn_input = getParameter(argc, argv, "-i");

 psf.produceSideInfo();

 Image<double>   inputVol, imOut;

 inputVol.read(fn_input);
 psf.adjustParam(inputVol);
 project_xr(psf, inputVol, imOut);

 if (checkParameter(argc, argv, "-out"))
  fn_output = getParameter(argc, argv, "-out");
 else
  fn_output = file_name.without_extension().add_extension("out").add_extension("spi");

 imOut.write(fn_output);
 }
 @endcode
 */
class XRayPSF
{
public:

    // Operation modes
    typedef enum
    {
        GENERATE_PSF,
        PSF_FROM_FILE
    } operMode;

    operMode mode;

    /// Define the selected PSF generation algorithm.
    PsfType type;

//    /// Current OTF
//    MultidimArray< std::complex<double> > OTF;
    ImageGeneric psfGen;
    /// 3D PSF read from file
    Image<double>  psfVol;
    /// Working PSF with nonlinear zdim whose slices are the mean PSF for slabs
    MultidimArray<double> PSF;
//    /// Axial intensity
//    MultidimArray<float> axialInt;
    /// Threshold to separate The volume into slabs to use the same PSF
    double slabThr;
    /// Z positions in the original PSF Volume to determine de slabs
    std::vector<int> slabIndex;
    // Transformation Matrix when reading PSF from file
    Matrix2D<double>  T;

    /// Lens shape Mask
    MultidimArray<double> *mask; //TODO: As in threadXrayProject a copy of this class is done, this pointer is never deleted. This class should use
    // applyOTF simultaneously from several threads

    /* RX Microscope configuration */
    /// Lens Aperture Radius
    double Rlens;
    /// Object plane on Focus (Reference)
    double Zo;
    /// Object plane
//    double Z;
    /// Image plane (CCD position)
    double Zi;
    /// Depth of focus. Only for information purposes
    double DoF;

    /* Digital Parameters */
    /// X size of the input image (object plane size)
    size_t Nox;
    /// Y size of the input image (object plane size)
    size_t Noy;
    /// Z size of the input image (object plane size)
    size_t Noz;
    /* Maximum pixel size in image plane (Minimum resolution condition).
     The same for both axis x-y, due to the symmetry of the lens aperture */
    double dxiMax;
    /// Pixel size in X-dim in lens plane
    double dxl;
    /// Pixel size in Y-dim in lens plane
    double dyl;
    /// Z limits around Zo in the psf generation due to Nyquist Limit
    double deltaZMaxX, deltaZMaxY, deltaZMinX, deltaZMinY;

    /// Parameters to change image size to avoid Nyquist limit
    PsfxrAdjust AdjustType;
    /// Minimum diameter size of the microscope pupile in the lens plane, measured in pixels
    double pupileSizeMin;

//    // Fourier Transformer to generate OTF, declared in class to avoid copy output
//    FourierTransformer ftGenOTF;

public:
    /// Lambda of illumination
    double lambda;

    /* RX Microscope configuration */
    /// Focal length
    double Flens;
    /// Number of zones in zone plate
    double Nzp;
    /// Outermost zone width
    double deltaR;
    /// Magnification
    double Ms;
    /// Z axis global shift
    double DeltaZo;


    /* Digital Parameters */
    /// object space XY-plane sampling rate
    double dxo;
    /// Image space XY-plane sampling rate
    double dxi;
    /// object space Z sampling rate
    double dzo;
    /// Size of the image in image plane, to be rescaled if needed
    size_t Nix, Niy;

    /// object space XY-plane sampling rate of the PSF Volume
    double dxoPSF;
    /// object space Z sampling rate of the PSF Volume
    double dzoPSF;

    /// Switch to control verbose mode
    int verbose;

    // Number of threads
    int nThr;

public:
    /** Empty constructor. */
    XRayPSF();

    /* Destructor
     */
    ~XRayPSF();

    /* Initialization of parameters
     */
    void init();

    /// Clear.
    void clear();

    /* Definition of params to be read from command line
     */
    static void defineParams(XmippProgram * program);

    /* Read of params from command line
     */
    void readParams(XmippProgram * program);

    /** Read from file.
        An exception is thrown if the file cannot be open.*/
    void read(const FileName &fn, bool readVolume = true);

    /** Write to file.
        An exception is thrown if the file cannot be open.*/
    void write(const FileName &fn);

    /// Show the microscope parameters
    void show();

    /// Show
    friend std::ostream & operator <<(std::ostream &out, const XRayPSF &psf);

    /// Add focal shift to previously read psf zshift
    void setFocalShift(double zShift);

    /// Produce Side information
    void calculateParams(double _dxo, double _dzo = -1, double threshold = 0.);

    /// Calculate the width of the slabs to reduce computing time and the mean PSF for each
    void reducePSF2Slabs(double threshold);

    /// Apply the OTF to the image, by means of the convolution
    void applyOTF(MultidimArray<double> &Im, const double sliceOffset) const;

    /// Generate the Optical Transfer Function (OTF) for a slice according to Microscope and Im parameters.
    void generateOTF(MultidimArray<std::complex<double> > &OTF, double Zpos) const;

    /// Generate the 3D Point Spread Function (PSF) according to Microscope parameters.
    void generatePSF();

    /// Calculate if a resize of the X-Y plane is needed to avoid the Nyquist Limit
    void adjustParam();
    void adjustParam(MultidimArray<double> &Vol);

protected:
    /// Generate the PSF for a single plane according to a ideal lens.
    void generatePSFIdealLens(MultidimArray<double> &PSFi, double Zpos) const;
};

/// Generate the quadratic phase distribution of a ideal lens
void lensPD(MultidimArray<std::complex<double> > &Im, double Flens, double lambda, double dx, double dy);

//@}
#endif
