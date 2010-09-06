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
void XRayPSF::read(const FileName &fn)
{
    FILE *fh_param;

    if (fn == "")
    {
        clear();
        return;
    }
    if ((fh_param = fopen(fn.c_str(), "r")) == NULL)
        REPORT_ERROR(ERR_IO_NOTOPEN,
                     (std::string)"XmippXROTF::read: There is a problem "
                     "opening the file " + fn);

    try
    {
        lambda = textToFloat(getParameter(fh_param, "lambda", 0, "2.43"))* 1e-9;
        deltaR = textToFloat(getParameter(fh_param, "outer_zone_width", 0, "40")) * 1e-9;
        Nzp = textToFloat(getParameter(fh_param, "zones_number", 0, "560"));
        Ms = textToFloat(getParameter(fh_param, "magnification", 0, "2304"));
        dxo = textToFloat(getParameter(fh_param, "sampling_rate", 0, "1")) *1e-9;
        if (checkParameter(fh_param, "z_sampling_rate"))
            dzo = textToFloat(getParameter(fh_param, "z_sampling_rate", 0)) *1e-9;
        else
            dzo = dxo;
        DeltaZo = textToFloat(getParameter(fh_param, "z_axis_shift", 0, "0")) *1e-6;

    }
    catch (XmippError XE)
    {
        std::cout << XE << std::endl;
        REPORT_ERROR(ERR_IO_NOREAD, (std::string)"There is an error reading " + fn);
    }
    fclose(fh_param);
}

/* Write ------------------------------------------------------------------- */
void XRayPSF::write(const FileName &fn)
{
    std::ofstream fh_param;
    fh_param.open(fn.c_str());
    if (!fh_param)
        REPORT_ERROR(ERR_IO_NOTOPEN, (std::string)"Xmipp_CTF::write: Cannot open " + fn +
                     " for output");
    fh_param << *this << std::endl;
    fh_param.close();
}

/* Usage ------------------------------------------------------------------- */
void XRayPSF::usage()
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
std::ostream & operator <<(std::ostream &out, const XRayPSF &psf)
{

    out     << std::endl
    << "XmippXRPSF:: X-Ray Microscope parameters:" << std::endl
    << "lambda=               " << psf.lambda * 1e9 << " nm" << std::endl
    << "zones_number=         " << psf.Nzp << std::endl
    << "outer_zone_width=     " << psf.deltaR * 1e9 << " nm" << std::endl
    << "magnification=        " << psf.Ms << std::endl
    << "sampling_rate=        " << psf.dxo * 1e9 << " nm" << std::endl
    << "z_sampling_rate=      " << psf.dzo * 1e9 << " nm" << std::endl
    << "z_axis_shift=         " << psf.DeltaZo * 1e6 << " um" << std::endl
    << std::endl
    << "focal_length=         " << psf.Flens * 1e3 << " mm" << std::endl
    << "Lens Radius=          " << psf.Rlens * 1e6 << " um" << std::endl
    << "Zo=                   " << psf.Zo * 1e3 << " mm" << std::endl
    << "Zi=                   " << psf.Zi * 1e3 << " mm" << std::endl
    << "Depth_of_Focus=       " << psf.DoF * 1e6 << " um" << std::endl
    << std::endl
    << "dxi=                  " << psf.dxi * 1e6 << " um" << std::endl
    << "dxiMax=               " << psf.dxiMax * 1e6 << " um" << std::endl
    ;

    return out;
}

/* Default values ---------------------------------------------------------- */
void XRayPSF::clear()
{
    lambda = 2.43e-9;
    //    Flens = 1.4742e-3;
    deltaR = 40e-9;
    Nzp = 560;
    Ms = 2304;
    dzo = dxo = 1e-9;
    DeltaZo = 0;
    pupileSizeMin = 5;
}

/* Produce Side Information ------------------------------------------------ */
void XRayPSF::produceSideInfo()
{
    /// Calculation of the rest of microscope parameters
    Flens = (4*Nzp*deltaR*deltaR)/lambda;
    //    Rlens = sqrt(Nzp * lambda * Flens);
    Rlens = 4 * Nzp * deltaR;
    Zo = (1 + 1 / Ms) * Flens;
    Zi = Zo * Ms;
    dxi = dxo * Ms;
    dxiMax = lambda * Zi / (2 * Rlens);
    DoF = 4*deltaR*deltaR/lambda;

    //        Z = 0.99999*Zo;
    //    Z = Zo;

    if (verbose)
        std::cout << *this << std::endl;
}

/* Apply the OTF to an image ----------------------------------------------- */
//void XRayPSF::applyOTF(MultidimArray<double> &I) const
//{
//}

/* Generate the Intensity PSF for a specific XR microscope configuration     ------------- */
/* Generate OTF Image ------------------------------------------------------ */

void XRayPSF::generateOTF(MultidimArray<std::complex<double> > &Im)
{
    MultidimArray<double> I2(Im.ydim,Im.xdim);
    generateOTF(I2);
}

void XRayPSF::generateOTF(MultidimArray<double> &Im)
{
    //#define DEBUG
    /// REMEMBER TO INCLUDE AND/OR ANALYZE THE MINIMUM RESOLUTION CONDITION !!!! ////


    //    double Nx = Im.xdim, Ny = Im.ydim;
    double focalEquiv = 1/(1/Z - 1/Zo); // inverse of defocus = 1/Z - 1/Zo
    //    double dxl = lambda*Zi / (Nx * dxi); // Pixel X-size en the plane of lens aperture
    //    double dyl = lambda*Zi / (Ny * dxi); // Pixel Y-size en the plane of lens aperture


    //#ifdef DEBUG
    //
    //    std::cout << std::endl;
    //    std::cout << "XRayPSF::GenerateOTF - Parameters:" << std::endl;
    //    std::cout << "Nx = " << Nx << "   Ny = " << Ny << std::endl;
    //    std::cout << "dxl = " << dxl << "   dyl = " << dyl << std::endl;
    //    std::cout << "Equivalent focal = " << focalEquiv << std::endl;
    //    std::cout << "Discrete X-Radius in pixels = " << Rlens / dxl  << std::endl;
    //    std::cout << "Discrete Y-Radius in pixels = " << Rlens / dyl  << std::endl;
    //    std::cout << std::endl;
    //#endif

    MultidimArray< std::complex<double> > OTFTemp, PSFi;
    FourierTransformer transformer;
    //    Mask_Params mask_prm; TODO do we have to include masks using this method?

    OTFTemp.resize(Niy,Nix);
    lensPD(OTFTemp, focalEquiv, lambda, dxl, dyl);

    //    OTFTemp.window(-128,-128,127,255,10);


    FOR_ALL_ELEMENTS_IN_ARRAY2D(OTFTemp)
    {
        if (sqrt(double(i*i)*dyl*dyl + double(j*j)*dxl*dxl) > Rlens)
            OTFTemp(i,j)=0;
        //        if (sqrt(double(i*i)+ double(j*j)) > 64) OTFTemp(i,j)=0;
    }

#ifdef DEBUG
    Image<double> _Im;
    _Im().resize(OTFTemp);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(OTFTemp)
    dAij(_Im(),i,j) = abs(dAij(OTFTemp,i,j));
    _Im.write("psfxr-lens.spi");
#endif


    transformer.FourierTransform(OTFTemp, PSFi, false);
    //    CenterOriginFFT(PSFi, 1);
    double norm=0;

    //     FOR_ALL_ELEMENTS_IN_ARRAY2D(LPDFT)
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(PSFi)
    {
        PSFi.data[n] = abs(PSFi.data[n]);
        PSFi.data[n] *= PSFi.data[n];
        norm +=  PSFi.data[n].real();
    }
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(PSFi)
    {
        PSFi.data[n] /= norm;
    }

    transformer.inverseFourierTransform();

    OTF.resize(OTFTemp.ydim, OTFTemp.xdim/2+1);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(OTF)
    dAij(OTF,i,j) = dAij(OTFTemp,i,j);




#ifdef DEBUG

    //    CenterOriginFFT(OTFTemp,1);
    _Im().resize(OTFTemp);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(OTFTemp)
    dAij(_Im(),i,j) = abs(dAij(OTFTemp,i,j));
    _Im.write("psfxr-otf1.spi");
    //    CenterOriginFFT(OTFTemp,0);
    _Im().resize(OTF);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(OTF)
    dAij(_Im(),i,j) = abs(dAij(OTF,i,j));
    _Im.write("psfxr-otf2.spi");

    //    CenterOriginFFT(PSFi,1);
    _Im().resize(PSFi);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(PSFi)
    dAij(_Im(),i,j) = abs(dAij(PSFi,i,j));
    _Im.write("psfxr-psfi.spi");
#endif
}

void XRayPSF::adjustParam(Image<double> &vol)
{

    produceSideInfo();

    Nox = vol().xdim;
    Noy = vol().ydim;

    dxl = lambda*Zi / (Nox * dxi); // Pixel X-size en the plane of lens aperture
    dyl = lambda*Zi / (Noy * dxi); // Pixel Y-size en the plane of lens aperture

    deltaZMaxX = Zo*((Rlens*2*dxl)/(Rlens*2*dxl - Zo*lambda) - 1);
    deltaZMaxY = Zo*((Rlens*2*dyl)/(Rlens*2*dyl - Zo*lambda) - 1);
    deltaZMinX = Zo*((Rlens*2*dxl)/(Rlens*2*dxl + Zo*lambda) - 1);
    deltaZMinY = Zo*((Rlens*2*dyl)/(Rlens*2*dyl + Zo*lambda) - 1);


    if (verbose)
    {
        std::cout << std::endl;
        std::cout << "XmippXRPSF::Param adjust:" << std::endl;
        std::cout << "(Nox,Noy,Nz) = (" << Nox << "," << Noy << "," << vol().zdim << ")" << std::endl;
        std::cout << "Larger volume Z plane = " << (ABS(DeltaZo)+vol().zdim/2*dzo)*1e6 << " um" << std::endl;
        std::cout << "Larger discrete Z plane (x,y) = (" << deltaZMaxX*1e6 << ", " << deltaZMaxY*1e6 << ") um" << std::endl;
    }

    if (dxi>dxiMax) /// Lens Radius in pixels higher than image
    {
        Nix = ceil(Nox * dxi/dxiMax);
        Niy = ceil(Noy * dxi/dxiMax);

        dxi *= Nox/Nix;

        AdjustType = PSFXR_INT;

        if (verbose)
            std::cout << "XmippXRPSF: Image plane sampling too small: increasing resolution" << std::endl;

    }
    else
    {
        Nix = Nox;
        Niy = Noy;
        AdjustType = PSFXR_STD;

        if (dxi < pupileSizeMin/Nox * dxiMax)
        {
            Nix = ceil(pupileSizeMin * dxiMax/dxi);
            AdjustType = PSFXR_ZPAD;
        }
        if (dxi < pupileSizeMin/Noy * dxiMax)
        {
            Niy = ceil(pupileSizeMin * dxiMax/dxi);
            AdjustType = PSFXR_ZPAD;
        }
        if (DeltaZo + vol().zdim/2*dzo > deltaZMaxX)
        {
            Nix = XMIPP_MAX(Nix,ceil(Zi*Rlens*2*ABS(DeltaZo+vol().zdim/2*dzo)/(Zo*dxi*(Zo+DeltaZo+vol().zdim/2*dzo))));
            AdjustType = PSFXR_ZPAD;
        }
        if (DeltaZo - vol().zdim/2*dzo < deltaZMinX)
        {
            Nix = XMIPP_MAX(Nix,ceil(Zi*Rlens*2*ABS(DeltaZo-vol().zdim/2*dzo)/(Zo*dxi*(Zo+DeltaZo-vol().zdim/2*dzo))));
            AdjustType = PSFXR_ZPAD;
        }
        if (DeltaZo + vol().zdim/2*dzo > deltaZMaxY)
        {
            Niy = XMIPP_MAX(Niy,ceil(Zi*Rlens*2*ABS(DeltaZo+vol().zdim/2*dzo)/(Zo*dxi*(Zo+DeltaZo+vol().zdim/2*dzo))));
            AdjustType = PSFXR_ZPAD;
        }
        if ( DeltaZo - vol().zdim/2*dzo < deltaZMinY)
        {
            Niy = XMIPP_MAX(Niy,ceil(Zi*Rlens*2*ABS(DeltaZo-vol().zdim/2*dzo)/(Zo*dxi*(Zo+DeltaZo-vol().zdim/2*dzo))));
            AdjustType = PSFXR_ZPAD;
        }
    }

    dxl = lambda*Zi / (Nix * dxi); // Pixel X-size en the plane of lens aperture
    dyl = lambda*Zi / (Niy * dxi); // Pixel Y-size en the plane of lens aperture

    if (verbose)
    {
        if (AdjustType!=PSFXR_STD)
        {
            if (AdjustType==PSFXR_ZPAD)
                std::cout << "XmippXRPSF: Image plane size too small: increasing size" << std::endl;
            std::cout << "New slice dimensions:  " <<  "(Nix, Niy) = (" << Nix << ", " << Niy << ")" << std::endl;
        }
        std::cout << "(dxl, dyl) = (" << dxl*1e6 << ", " << dyl*1e6 << ") um" << std::endl;
        std::cout << "Pupile Diameter in pixels (NpX, NpY) = (" << ceil(2*Rlens / dxl)
        << ", " << ceil(2*Rlens / dyl) << ")"  << std::endl;
        std::cout << std::endl;
    }


}

/* Generate the quadratic phase distribution of an ideal lens ------------- */
void lensPD(MultidimArray<std::complex<double> > &Im, double Flens, double lambda, double dx, double dy)
{

    double Lx0 = Im.xdim * dx, Ly0 = Im.ydim * dy, x, y, phase;

    Im.setXmippOrigin();

    FOR_ALL_ELEMENTS_IN_ARRAY2D(Im)
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
/// TODO: func description
void project_xr(XRayPSF &psf, Image<double> &vol, Image<double> &imOut, int idxSlice)
{

    // psf.adjustParam(vol);

    imOut() = MultidimArray<double> (psf.Niy, psf.Nix);
    imOut().initZeros();
    imOut().setXmippOrigin();

    MultidimArray<double> imTemp(psf.Noy, psf.Nox), intExp(psf.Noy, psf.Nox), imTempSc(imOut()), *imTempP;
    intExp.initZeros();
    imTemp.initZeros();
    imTempSc.initZeros();
    intExp.setXmippOrigin();
    imTemp.setXmippOrigin();
    imTempSc.setXmippOrigin();

    vol().setXmippOrigin();

    //#define DEBUG
#ifdef DEBUG

    Image<double> _Im(imOut);
#endif

        init_progress_bar(vol().zdim-1);

    for (int k=((vol()).zinit); k<=((vol()).zinit + (vol()).zdim - 1); k++)
    {
        FOR_ALL_ELEMENTS_IN_ARRAY2D(intExp)
        {
            intExp(i, j) = intExp(i, j) + vol(k, i, j);
            imTemp(i, j) = (exp(-intExp(i,j)*psf.dzo))*vol(k,i,j)*psf.dzo;
            //            imTemp(i, j) = 1./(exp(intExp(i,j)))*vol(k,i,j);
        }
#ifdef DEBUG
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(imTemp)
        dAij(_Im(),i,j) = dAij(imTemp,i,j);
        _Im.write("psfxr-imTemp.spi");
#endif

        psf.Z = psf.Zo + psf.DeltaZo - k*psf.dzo + imOut.Zoff();

        switch (psf.AdjustType)
        {
        case PSFXR_INT:
            imTempP = &imTempSc;
            scaleToSize(LINEAR,*imTempP,imTemp,psf.Nix,psf.Niy);
            //          imTemp.scaleToSize(psf.Niy, psf.Nix, *imTempP);
            break;

        case PSFXR_STD:
            imTempP = &imTemp;
            break;

        case PSFXR_ZPAD:
            //    (*imTempSc).resize(imTemp);
            imTempSc = imTemp;
            imTempSc.window(-ROUND(psf.Niy/2)+1,-ROUND(psf.Nix/2)+1,ROUND(psf.Niy/2)-1,ROUND(psf.Nix/2)-1);
            imTempP = &imTempSc;
            break;
        }

        psf.generateOTF(*imTempP);




#ifdef DEBUG

        _Im().resize(intExp);
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(intExp)
        dAij(_Im(),i,j) = dAij(intExp,i,j);
        _Im.write("psfxr-intExp.spi");
        _Im().resize(*imTempP);
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(*imTempP)
        dAij(_Im(),i,j) = dAij(*imTempP,i,j);
        _Im.write("psfxr-imTempEsc.spi");
#endif

        psf.applyOTF(*imTempP);

#ifdef DEBUG

        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(*imTempP)
        dAij(_Im(),i,j) = dAij(*imTempP,i,j);
        _Im.write("psfxr-imTempEsc2.spi");
#endif

        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(*imTempP)
        dAij(imOut(),i,j) += dAij(*imTempP,i,j);

        //        imOut.write("psfxr-imout.spi");

        if (idxSlice > -1)
            progress_bar((idxSlice - 1)*vol().zdim + k - vol().zinit);
    }

    imOut() = 1-imOut();
    //            imOut.write("psfxr-imout.spi");


    switch (psf.AdjustType)
    {
    case PSFXR_INT:
        //    imOut().selfScaleToSize(psf.Noy, psf.Nox);
        selfScaleToSize(LINEAR,imOut(), psf.Nox, psf.Noy);
        break;

    case PSFXR_ZPAD:
        imOut().window(-ROUND(psf.Noy/2)+1,-ROUND(psf.Nox/2)+1,ROUND(psf.Noy/2)-1,ROUND(psf.Nox/2)-1);
        break;
    }
    imOut().setXmippOrigin();
    //    imOut.write("psfxr-imout2.spi");

}

#undef DEBUG

