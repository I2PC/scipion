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

#include "psf_xr.h"
#include "args.h"
#include "fft.h"
#include "mask.h"

/* Read -------------------------------------------------------------------- */
void XmippXRPSF::read(const FileName &fn)
{
    FILE *fh_param;
    if ((fh_param = fopen(fn.c_str(), "r")) == NULL)
        REPORT_ERROR(1,
                     (std::string)"XmippXROTF::read: There is a problem "
                     "opening the file " + fn);

    try
    {
        lambda = textToFloat(getParameter(fh_param, "lambda", 0, "2.43"))* 1e-9;
        Flens = textToFloat(getParameter(fh_param, "focal_length", 0, "1.4742")) * 1e-3;
        Nzp = textToFloat(getParameter(fh_param, "zones_number", 0, "560"));
        Ms = textToFloat(getParameter(fh_param, "magnification", 0, "2304"));
        dxo = textToFloat(getParameter(fh_param, "sampling_rate", 0, "1")) *1e-9;
        if (checkParameter(fh_param, "z_sampling_rate"))
            dzo = textToFloat(getParameter(fh_param, "z_sampling_rate", 0)) *1e-9;
        else
            dzo = dxo;
        DeltaZo = textToFloat(getParameter(fh_param, "z_axis_shift", 0, "0")) *1e-6;

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE << std::endl;
        REPORT_ERROR(1, (std::string)"There is an error reading " + fn);
    }
    fclose(fh_param);
}

/* Write ------------------------------------------------------------------- */
void XmippXRPSF::write(const FileName &fn)
{
    std::ofstream fh_param;
    fh_param.open(fn.c_str());
    if (!fh_param)
        REPORT_ERROR(1, (std::string)"Xmipp_CTF::write: Cannot open " + fn +
                     " for output");
    fh_param << *this << std::endl;
    fh_param.close();
}

/* Usage ------------------------------------------------------------------- */
void XmippXRPSF::usage()
{
    std::cerr << "  [lambda=<lambda=2.43>]            : Wavelength in nm\n"
    << "  [focal_length=<Flens=1.4742>]     : Focal length in mm\n"
    << "  [zones_number=<Nzp=560>]          : zone plate number\n"
    << "  [magnification=<Ms=2304>]         : Microscope magnification\n"
    << "  [sampling_rate=<dxo=1>]           : Object XY plane resolution in nm/pixel\n"
    << "  [z_sampling_rate=<dzo=dxo>]       : Object Z axis resolution in nm/pixel\n"
    << "  [z_axis_shift=<DeltaZo=0>]        : Z axis shift in um\n";
}

/* Show -------------------------------------------------------------------- */
std::ostream & operator <<(std::ostream &out, const XmippXRPSF &psf)
{

    out << "lambda=               " << psf.lambda * 1e9 << " nm" << std::endl
    << "focal_length=         " << psf.Flens * 1e3 << " mm" << std::endl
    << "zones_number=         " << psf.Nzp << std::endl
    << "magnification=        " << psf.Ms << std::endl
    << "sampling_rate=        " << psf.dxo * 1e9 << " nm" << std::endl
    << "z_sampling_rate=      " << psf.dzo * 1e9 << " nm" << std::endl
    << "DeltaZo=              " << psf.DeltaZo * 1e6 << " um" << std::endl
    << std::endl
    << "Lens Radius=          " << psf.Rlens * 1e6 << " um" << std::endl
    << "Zo=                   " << psf.Zo * 1e3 << " mm" << std::endl
    << "Zi=                   " << psf.Zi * 1e3 << " mm" << std::endl
    << "Minimum Resolution=   " << psf.dxiMax * 1e6 << " um" << std::endl
    ;

    return out;
}

/* Default values ---------------------------------------------------------- */
void XmippXRPSF::clear()
{
    lambda = 2.43e-9;
    Flens = 1.4742e-3;
    Nzp = 560;
    Ms = 2304;
    dzo = dxo = 1e-9;
    DeltaZo = 0;
}

/* Produce Side Information ------------------------------------------------ */
void XmippXRPSF::produceSideInfo()
{
    /// Calculation of the rest of microscope parameters
    Rlens = sqrt(Nzp * lambda * Flens);
    Zo = (1 + 1 / Ms) * Flens;
    Zi = Zo * Ms;
    dxi = dxo * Ms;
    dxiMax = lambda * Zi / (2 * Rlens);

    //        Z = 0.99999*Zo;
    Z = Zo;
}

/* Apply the OTF to an image ----------------------------------------------- */
//void XmippXRPSF::applyOTF(Matrix2D<double> &I) const
//{
//}

/* Generate the Intensity PSF for a specific XR microscope configuration     ------------- */
/* Generate OTF Image ------------------------------------------------------ */
//#define DEBUG
void XmippXRPSF::generateOTF(Matrix2D<std::complex<double> > &Im)
{
    Matrix2D<double> I2(Im.ydim,Im.xdim);
    generateOTF(I2);
}

void XmippXRPSF::generateOTF(Matrix2D<double> &Im)
{
#define DEBUG
    /// REMEMBER TO INCLUDE AND/OR ANALYZE THE MINIMUM RESOLUTION CONDITION !!!! ////


    double Nx = Im.xdim, Ny = Im.ydim;
    double focalEquiv = 1/(1/(Z + DeltaZo) - 1/Zo); // inverse of defocus = 1/Z - 1/Zo
    double dxl = lambda*Zi / (Nx * dxi); // Pixel X-size en the plane of lens aperture
    double dyl = lambda*Zi / (Ny * dxi); // Pixel Y-size en the plane of lens aperture

#ifdef DEBUG

    std::cout << std::endl;
    std::cout << "XmippXRPSF::GenerateOTF - Parameters:" << std::endl;
    std::cout << "Nx = " << Nx << "   Ny = " << Ny << std::endl;
    std::cout << "dxl = " << dxl << "   dyl = " << dyl << std::endl;
    std::cout << "Equivalent focal = " << focalEquiv << std::endl;
    std::cout << "Discrete X-Radius in pixels = " << Rlens / dxl  << std::endl;
    std::cout << "Discrete Y-Radius in pixels = " << Rlens / dyl  << std::endl;
    std::cout << std::endl;
#endif

    Matrix2D< std::complex<double> > PSFi;
    XmippFftw transformer;
    //    Mask_Params mask_prm; TODO do we have to include masks using this method?

    OTF.resize(Ny,Nx);
    lensPD(OTF, focalEquiv, lambda, dxl, dyl);

    //    OTF.window(-128,-128,127,255,10);


    FOR_ALL_ELEMENTS_IN_MATRIX2D(OTF)
    {
        if (sqrt(double(i*i)*dyl*dyl + double(j*j)*dxl*dxl) > Rlens)
            OTF(i,j)=0;
        //        if (sqrt(double(i*i)+ double(j*j)) > 64) OTF(i,j)=0;
    }

#ifdef DEBUG
    ImageXmipp _Im;
    _Im().resize(OTF);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(OTF)
    _Im(i,j) = abs(OTF(i,j));
    _Im.write("psfxr-lens.spi");
#endif


    transformer.FourierTransform(OTF, PSFi, false);
    //    CenterOriginFFT(PSFi, 1);
    double norm=0;

    //     FOR_ALL_ELEMENTS_IN_MATRIX2D(LPDFT)
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(PSFi)
    {
        PSFi.data[n] = abs(PSFi.data[n]);
        PSFi.data[n] *= PSFi.data[n];
        norm +=  PSFi.data[n].real();
    }
    PSFi /= norm;

    transformer.inverseFourierTransform();

#ifdef DEBUG

    CenterOriginFFT(OTF,1);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(OTF)
    _Im(i,j) = abs(OTF(i,j));
    _Im.write("psfxr-otf.spi");
    CenterOriginFFT(OTF,0);

    CenterOriginFFT(PSFi,1);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(PSFi)
    _Im(i,j) = abs(PSFi(i,j));
    _Im.write("psfxr-psfi.spi");
#endif
}

/* Generate the quadratic phase distribution of an ideal lens ------------- */
void lensPD(Matrix2D<std::complex<double> > &Im, double Flens, double lambda, double dx, double dy)
{

    double Lx0 = Im.xdim * dx, Ly0 = Im.ydim * dy, x, y, phase;

    Im.setXmippOrigin();

    FOR_ALL_ELEMENTS_IN_MATRIX2D(Im)
    {
        /// For indices in standard fashion
        //   x = (double) j * dx + (dx - Lx0) / 2;
        //   y = (double) i * dx + (dx - Ly0) / 2;

        x = (double) j * dx + (dx) / 2;
        y = (double) i * dy + (dy) / 2;

        phase = (-PI / lambda / Flens) * (x * x + y * y);

        (Im(i, j)).real() = cos(phase);
        (Im(i, j)).imag() = sin(phase);

    }
}

void project_xr(XmippXRPSF &psf, VolumeXmipp &vol, ImageXmipp &imOut)
{
    imOut() = Matrix2D<double> (vol().ydim, vol().xdim);
    imOut().initZeros();
//    imOut()+= 1;
    imOut().setXmippOrigin();

    Matrix2D<double> imTemp(imOut()), intExp(imOut());
    intExp.initZeros();

    vol().setXmippOrigin();

#define DEBUG
#ifdef DEBUG
    ImageXmipp _Im(imOut);
    _Im().setXmippOrigin();

#endif

    for (int k=((vol()).zinit); k<=((vol()).zinit + (vol()).zdim - 1); k++)
    {
        FOR_ALL_ELEMENTS_IN_MATRIX2D(intExp)
        {
            intExp(i, j) = intExp(i, j) + vol(k, i, j);
//            imTemp(i, j) = (intExp(i,j)*psf.dzo*(-1))*vol(k,i,j)*psf.dzo;
            imTemp(i, j) = exp(intExp(i,j)*(-1))*vol(k,i,j);
            _Im(i,j) = imTemp(i,j);
        }
        psf.Z = psf.Zo - k*psf.dzo;
        psf.generateOTF(imTemp);
        psf.applyOTF(imTemp);
        imOut() += imTemp;

#ifdef DEBUG
//        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(imTemp)
//        _Im(i,j) = abs(imTemp(i,j));

        _Im.write("psfxr-imTemp.spi");
#endif
    }

    imOut() = 1-imOut();

}


#undef DEBUG

