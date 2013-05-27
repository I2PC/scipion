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


XRayPSF::XRayPSF()
{
    init();
}

XRayPSF::~XRayPSF()
{
    clear();
}

/* Definition of command line parameters
 */
void XRayPSF::defineParams(XmippProgram * program)
{
    //    program->addParamsLine(" [-file <psf_param_file>]       : Read X-ray psf parameters from file.");
    program->addParamsLine(" [--lambda <v=2.43>]             : X-ray wavelength (nm).");
    program->addParamsLine(" [--out_width <deltaR=40>]       : Outermost zone width of the X-ray Fresnel lens (nm).");
    program->addParamsLine(" [--zones <N=560>]               : Number of zones of the X-ray Fresnel lens.");
    program->addParamsLine(" [--mag <Ms=2304>]               : Magnification of the X-ray microscope.");
    program->addParamsLine(" [--sampling <dxy=10> <dz=dxy>]  : Sampling rate in X-Y plane and Z axis (nm).");
    program->addParamsLine(" [--zshift <deltaZ=0>]           : Longitudinal displacement along Z axis in microns.");
    program->addParamsLine(" [--size <x> <y=x> <z=x>]        : Size of the X-ray PSF volume.");
    program->addParamsLine(" [--type  <lens_type=ideal>]     : Lens type to generate the PSF.");
    program->addParamsLine("           where <lens_type>");
    program->addParamsLine("                  ideal          : Ideal phase Fresnel lens.");
    program->addParamsLine("                  zp             : Fresnel Zone Plate lens.");

}

/* Read params
 */
void XRayPSF::readParams(XmippProgram * program)
{
    lambda  = program->getDoubleParam("--lambda")*1e-9;
    deltaR  = program->getDoubleParam("--out_width")*1e-9;
    Nzp     = program->getDoubleParam("--zones");
    Ms      = program->getDoubleParam("--mag");
    dxoPSF  = program->getDoubleParam("--sampling",0)*1e-9;
    dzoPSF  = STR_EQUAL(program->getParam("--sampling", 1), "dxy") ? dxoPSF : program->getDoubleParam("--sampling", 1)*1e-9;
    DeltaZo = program->getDoubleParam("--zshift")*1e-6;
    Nox     = program->getIntParam("--size",0);
    Noy     = STR_EQUAL(program->getParam("--size", 1),"x") ? Nox : program->getIntParam("--size", 1);
    Noz     = STR_EQUAL(program->getParam("--size", 2),"x") ? Nox : program->getIntParam("--size", 2);
    type    = STR_EQUAL(program->getParam("--type"), "zp") ? ANALYTIC_ZP : IDEAL_FRESNEL_LENS;

    verbose = program->getIntParam("-v");
}

/* Read -------------------------------------------------------------------- */
void XRayPSF::read(const FileName &fn, bool readVolume)
{

    if (fn == "")
    {
        clear();
        return;
    }

    if (fn.isMetaData())
    {
        MetaData MD;
        MD.read(fn);
        size_t id = MD.firstObject();

        FileName fnPSF;
        if ( MD.getValue(MDL_IMAGE,fnPSF, id) && readVolume )
        {
            mode = PSF_FROM_FILE;
            psfGen.readMapped(fn.getDir()+fnPSF);
            ArrayDim aDim;
            psfGen.getDimensions(aDim);
            //            MULTIDIM_ARRAY(psfVol).getDimensions(aDim);
            Nox = aDim.xdim;
            Noy = aDim.ydim;
            Noz = aDim.zdim;

            psfGen().setXmippOrigin();

            //            MULTIDIM_ARRAY(psfVol).setXmippOrigin();
        }
        else
        {
            std::vector<double> dimV;
            if (!MD.getValue(MDL_CTF_DIMENSIONS,dimV, id))
                REPORT_ERROR(ERR_ARG_MISSING, MDL::label2Str(MDL_CTF_DIMENSIONS) + " argument not present.");
            Nox = (int)dimV[0];
            Noy = (int)dimV[1];
            Noz = (int)dimV[2];
        }
        String typeS;
        if (!MD.getValue(MDL_CTF_XRAY_LENS_TYPE, typeS, id))
            REPORT_ERROR(ERR_ARG_MISSING, MDL::label2Str(MDL_CTF_XRAY_LENS_TYPE) + " argument not present.");
        if (typeS == "ZP")
            type = ANALYTIC_ZP;
        else
            type = IDEAL_FRESNEL_LENS;

        if (!MD.getValue(MDL_CTF_LAMBDA,lambda, id))
            REPORT_ERROR(ERR_ARG_MISSING, MDL::label2Str(MDL_CTF_LAMBDA) + " argument not present.");
        lambda *= 1e-9;
        if (!MD.getValue(MDL_CTF_XRAY_OUTER_ZONE_WIDTH,deltaR, id))
            REPORT_ERROR(ERR_ARG_MISSING, MDL::label2Str(MDL_CTF_XRAY_OUTER_ZONE_WIDTH) + " argument not present.");
        deltaR *= 1e-9;
        if (!MD.getValue(MDL_CTF_XRAY_ZONES_NUMBER,Nzp, id))
            REPORT_ERROR(ERR_ARG_MISSING, MDL::label2Str(MDL_CTF_XRAY_ZONES_NUMBER) + " argument not present.");
        if (!MD.getValue(MDL_MAGNIFICATION,Ms, id))
            REPORT_ERROR(ERR_ARG_MISSING, MDL::label2Str(MDL_MAGNIFICATION) + " argument not present.");
        if (!MD.getValue(MDL_CTF_SAMPLING_RATE,dxoPSF, id))
            REPORT_ERROR(ERR_ARG_MISSING, MDL::label2Str(MDL_CTF_SAMPLING_RATE) + " argument not present.");
        dxoPSF *= 1e-9;
        if (!MD.getValue(MDL_CTF_SAMPLING_RATE_Z,dzoPSF, id))
            REPORT_ERROR(ERR_ARG_MISSING, MDL::label2Str(MDL_CTF_SAMPLING_RATE_Z) + " argument not present.");
        dzoPSF *= 1e-9;
        if (!MD.getValue(MDL_CTF_LONGITUDINAL_DISPLACEMENT,DeltaZo, id))
            REPORT_ERROR(ERR_ARG_MISSING, MDL::label2Str(MDL_CTF_LONGITUDINAL_DISPLACEMENT) + " argument not present.");
        DeltaZo *= 1e-6;
    }
    else
    {
        FILE *fh_param;
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
            dxoPSF = textToFloat(getParameter(fh_param, "sampling_rate", 0, "1")) *1e-9;
            if (checkParameter(fh_param, "z_sampling_rate"))
                dzoPSF = textToFloat(getParameter(fh_param, "z_sampling_rate", 0)) *1e-9;
            else
                dzoPSF = dxoPSF;
            DeltaZo = textToFloat(getParameter(fh_param, "z_axis_shift", 0, "0")) *1e-6;

            Nox = textToInteger(getParameter(fh_param, "x_dim", 0, "0"));

        }
        catch (XmippError &XE)
        {
            std::cout << XE << std::endl;
            REPORT_ERROR(ERR_IO_NOREAD, (std::string)"There is an error reading " + fn);
        }
        fclose(fh_param);
    }
}

/* Write ------------------------------------------------------------------- */
void XRayPSF::write(const FileName &fn)
{
    MetaData MD;
    MD.setColumnFormat(false);
    size_t id = MD.addObject();

    FileName fnPSF = fn.withoutExtension().addExtension("vol");

    MD.setValue(MDL_IMAGE, fnPSF, id);
    String typeS = (type==ANALYTIC_ZP)? "ZP": "ideal";
    MD.setValue(MDL_CTF_XRAY_LENS_TYPE, typeS, id);
    MD.setValue(MDL_CTF_LAMBDA,lambda*1e9, id);
    MD.setValue(MDL_CTF_XRAY_OUTER_ZONE_WIDTH,deltaR*1e9, id);
    MD.setValue(MDL_CTF_XRAY_ZONES_NUMBER,Nzp, id);
    MD.setValue(MDL_MAGNIFICATION,Ms, id);
    MD.setValue(MDL_CTF_SAMPLING_RATE,dxoPSF*1e9, id);
    MD.setValue(MDL_CTF_SAMPLING_RATE_Z,dzoPSF*1e9, id);
    MD.setValue(MDL_CTF_LONGITUDINAL_DISPLACEMENT,DeltaZo*1e6, id);
    std::vector<double> dimV(3);
    dimV[0] = Nox;
    dimV[1] = Noy;
    dimV[2] = Noz;
    MD.setValue(MDL_CTF_DIMENSIONS,dimV, id);

    MD.write(fn.withoutExtension().addExtension("xmd"));
    //    PSF->write(fnPSF);
    psfVol.write(fnPSF);
}

/* Show -------------------------------------------------------------------- */
std::ostream & operator <<(std::ostream &out, const XRayPSF &psf)
{

    out     << std::endl
    << "--------------------------------------" << std::endl
    << "XrayPSF:: X-Ray Microscope parameters:" << std::endl
    << "--------------------------------------" << std::endl
    << "type =             " << ((psf.type==ANALYTIC_ZP)? "ZP": "ideal") << std::endl
    << "lambda =           " << psf.lambda * 1e9 << " nm" << std::endl
    << "zones_number =     " << psf.Nzp << std::endl
    << "outer_zone_width = " << psf.deltaR * 1e9 << " nm" << std::endl
    << "magnification =    " << psf.Ms << std::endl
    << "sampling_rate =    " << psf.dxo * 1e9 << " nm" << std::endl
    << "z_sampling_rate =  " << psf.dzo * 1e9 << " nm" << std::endl
    << "z_axis_shift =     " << psf.DeltaZo * 1e6 << " um" << std::endl
    << std::endl
    << "focal_length =     " << psf.Flens * 1e3 << " mm" << std::endl
    << "Lens Radius =      " << psf.Rlens * 1e6 << " um" << std::endl
    << "Zo =               " << psf.Zo * 1e3 << " mm" << std::endl
    << "Zi =               " << psf.Zi * 1e3 << " mm" << std::endl
    << "Depth_of_Focus =   " << psf.DoF * 1e6 << " um" << std::endl
    << std::endl
    << "dxi =              " << psf.dxi * 1e6 << " um" << std::endl
    << "dxiMax =           " << psf.dxiMax * 1e6 << " um" << std::endl
    ;

    return out;
}

/* Show the microscope parameters------------------------------------------- */
void XRayPSF::show()
{
    std::cout << *this << std::endl;
}

/* Initialization of parameters -------------------------------------------  */
void XRayPSF::init()
{
    lambda = 2.43e-9;
    //    Flens = 1.4742e-3;
    deltaR = 40e-9;
    Nzp = 560;
    Ms = 2304;
    dzo = dxo = 1e-9;
    DeltaZo = 0;
    pupileSizeMin = 5;

    mode = GENERATE_PSF;
    type = IDEAL_FRESNEL_LENS;

    //    ftGenOTF.setNormalizationSign(FFTW_BACKWARD);
}
/* Default values ---------------------------------------------------------- */
void XRayPSF::clear()
{
    psfVol.clear();
    init();
}

/* Produce Side Information ------------------------------------------------ */
void XRayPSF::calculateParams(double _dxo, double _dzo, double threshold)
{

    dxo = (_dxo > 0)? _dxo : dxoPSF;
    dzo = (_dzo > 0)? _dzo : dxo;

    /// Calculation of the rest of microscope parameters

    Rlens = (4 * Nzp * deltaR)*0.5; // Eq. 9.13 from Soft X-rays and .... David Attwood
    Flens = (2 * Rlens * deltaR)/lambda;

    Zo = (1 + 1 / Ms) * Flens;
    Zi = Zo * Ms;
    dxi = dxo * Ms;
    dxiMax = lambda * Zi / (2 * Rlens);
    DoF = 4*deltaR*deltaR/lambda;

    if(verbose > 0)
    {
        show();
        verbose++; // Show only once the param adjust info
    }

    if (mode == PSF_FROM_FILE)
    {
        T.initIdentity(4);
        double scaleFactor = dxo/dxoPSF;
        //        double scaleFactorZ = dxo/dzoPSF;
        double scaleFactorZ = 1;

        dMij(T, 0, 0) = scaleFactor;
        dMij(T, 1, 1) = scaleFactor;
        dMij(T, 2, 2) = scaleFactorZ;
        /** TODO: The scale in Z may not be necessary to save memory
         *  and should be used the nearest PSF slice for each Z plane */

        MultidimArray<double> &mdaPsfVol = psfVol();

        mdaPsfVol.resize(1, (size_t)(Noz/scaleFactorZ),
                         (size_t)(Noy/scaleFactor),
                         (size_t)(Nox/scaleFactor), false);

        mdaPsfVol.setXmippOrigin();

        applyGeometry(LINEAR, mdaPsfVol, psfGen(), T,
                      IS_INV, DONT_WRAP);

        psfGen.clear(); // Free mem in case image is not mapped

        // Rescale the PSF values to keep the normalization in each X-Y plane
        mdaPsfVol *= scaleFactor*scaleFactor;

        // Reset the sampling values so now psfVol matches input phantom's sampling
        dxoPSF = dxo;
        //        dzoPSF = dzo;
        T.initIdentity(4);

        if (threshold != 0.)
            reducePSF2Slabs(threshold);
    }
    else // Mask is used when creating PSF on demand
    {
        mask = new MultidimArray<double>;
    }
}

void XRayPSF::reducePSF2Slabs(double threshold)
{
    // find max position
    ArrayCoord maxPos;
    MultidimArray<double> &mdaPsfVol = psfVol();

    //MULTIDIM_ARRAY_GENERIC(PSFGen).maxIndex(maxPos);
    mdaPsfVol.maxIndex(maxPos);
    double maxValue = A3D_ELEM(mdaPsfVol,maxPos.z,maxPos.y,maxPos.x);
    threshold *= maxValue;
    slabIndex.clear();

    double maxTemp = maxValue;
    for (int k = 0; k > STARTINGZ(mdaPsfVol); --k)
    {
        double tempV = A3D_ELEM(mdaPsfVol,k,maxPos.y,maxPos.x);
        if ( maxTemp - tempV >  threshold)
        {
            slabIndex.insert(slabIndex.begin(), k);
            maxTemp = tempV;
        }
    }
    slabIndex.insert(slabIndex.begin(), STARTINGZ(mdaPsfVol));

    maxTemp = maxValue;
    for (int k = 0; k < FINISHINGZ(mdaPsfVol); ++k)
    {
        double tempV = A3D_ELEM(mdaPsfVol,k,maxPos.y,maxPos.x);
        if ( maxTemp - tempV >  threshold)
        {
            slabIndex.push_back(k);
            maxTemp = tempV;
        }
    }
    slabIndex.push_back(FINISHINGZ(mdaPsfVol));

} //reducePSF2Slabs

/* Apply the OTF to an image ----------------------------------------------- */
void XRayPSF::applyOTF(MultidimArray<double> &Im, const double sliceOffset) const
{
    //    double Z;
    //    MultidimArray< std::complex<double> > OTF;

    if (mode == GENERATE_PSF)
        REPORT_ERROR(ERR_VALUE_NOTSET, "XrayPSF::applyOTF: a PSF from file must be read.");

    double Z = Zo + DeltaZo + sliceOffset;

    MultidimArray<std::complex<double> > OTF;
    generateOTF(OTF, Z);

    MultidimArray<std::complex<double> > ImFT;
    FourierTransformer FTransAppOTF(FFTW_BACKWARD);

    //#define DEBUG
#ifdef DEBUG

    Image<double> _Im;
    _Im() = Im;

    _Im.write(("psfxr-ImIn.spi"));
#endif

    FTransAppOTF.FourierTransform(Im, ImFT, false);

    ////#define DEBUG
#ifdef DEBUG

    _Im() = Im;

    _Im.write(("psfxr-ImIn_after.spi"));

    //    Image<double> _Im;
    _Im().clear();
    _Im().resizeNoCopy(OTF);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(OTF)
    dAij(_Im(),i,j) = abs(dAij(OTF,i,j));
    _Im.write("psfxr-otf.spi");

    _Im().resizeNoCopy(ImFT);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(ImFT)
    dAij(_Im(),i,j) = abs(dAij(ImFT,i,j));
    _Im.write(("psfxr-imft1.spi"));
#endif

    //    double normSize = MULTIDIM_SIZE(OTF);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(ImFT)
    dAij(ImFT,i,j) *= dAij(OTF,i,j);

#ifdef DEBUG

    _Im().resize(ImFT);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(ImFT)
    dAij(_Im(),i,j) = abs(dAij(ImFT,i,j));
    _Im.write(("psfxr-imft2.spi"));
#endif

    FTransAppOTF.inverseFourierTransform();

#ifdef DEBUG

    _Im() = Im;

    //    _Im().resize(Im);
    //    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Im)
    //    dAij(_Im(),i,j) = abs(dAij(Im,i,j));
    _Im.write(("psfxr-imout.spi"));
#endif
#undef DEBUG
}

void XRayPSF::generateOTF(MultidimArray<std::complex<double> > &OTF, double Zpos) const
{

    MultidimArray<double> PSFi;
    FourierTransformer ftGenOTF(FFTW_BACKWARD);
    switch (mode)
    {
    case GENERATE_PSF:
        {
            if (type == IDEAL_FRESNEL_LENS)
                generatePSFIdealLens(PSFi, Zpos);
            break;
        }
    case PSF_FROM_FILE:
        {
            PSFi.resizeNoCopy(Niy, Nix);
            PSFi.setXmippOrigin();

            const MultidimArray<double> &mPsfVol = MULTIDIM_ARRAY(psfVol);

            double zIndexPSF = (Zo - Zpos)/dzoPSF; // Slice index for Zpos

            if (zIndexPSF < STARTINGZ(mPsfVol) ||
                zIndexPSF > FINISHINGZ(mPsfVol))
                PSFi.initConstant(1/YXSIZE(PSFi));
            else
            {   /* Actually T transform matrix is set to identity, as sampling is the same for both psfVol and phantom
                                 * It is only missing to set Z shift to select the slice */
                dMij(T, 2, 3) = zIndexPSF; // Distance from the focal plane
                applyGeometry(LINEAR, PSFi, mPsfVol, T,
                              IS_INV, DONT_WRAP, dAkij(mPsfVol,0,0,0));
                CenterFFT(PSFi, true);
            }
            break;
        }
    }

    ftGenOTF.FourierTransform(PSFi, OTF);

    //#define DEBUG
#ifdef DEBUG

    Image<double> _Im;
    //    _Im() = PSFi;
    CenterFFT(PSFi, false);
    _Im().alias(PSFi);
    _Im.write("psfxr-psfi.spi");
    //    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(PSFi)
    //    dAij(_Im(),i,j) = abs(dAij(PSFi,i,j));
    //    _Im.write("psfxr-psfi.spi");
    //    CenterOriginFFT(OTFTemp,1);
    //    _Im().resize(OTFTemp);
    //    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(OTFTemp)
    //    dAij(_Im(),i,j) = abs(dAij(OTFTemp,i,j));
    //    _Im.write("psfxr-otf1.spi");
    //    CenterOriginFFT(OTFTemp,0);
    _Im.clear();
    _Im().resize(OTF);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(OTF)
    dAij(_Im(),i,j) = abs(dAij(OTF,i,j));
    _Im.write("psfxr-otf2.spi");

#endif
#undef DEBUG
}

/* Generate the Intensity PSF for a specific XR microscope configuration     ------------- */
void XRayPSF::generatePSF()
{
    calculateParams(dxoPSF, dzoPSF);
    adjustParam();

    //    PSF = new Image<double>(Nox, Noy, Noz);
    //    VOLMATRIX(*PSF).setXmippOrigin();

    MultidimArray<double> &mPsfVol = MULTIDIM_ARRAY(psfVol);
    mPsfVol.resize(1,Noz, Niy, Nix, false);
    mPsfVol.setXmippOrigin();

    MultidimArray<double> PSFi;

    init_progress_bar(ZSIZE(mPsfVol));

    for (int k = STARTINGZ(mPsfVol), k2 = 0; k <= FINISHINGZ(mPsfVol); k++, k2++)
    {
        /* We keep sign of Z, Zo and DeltaZo positives in object space for the sake of simplicity in calculations,
         * it is, they increase opposite to Zdim in Xmipp reference system, so this is the reason to calculate
         * the plane Z by subtracting k*dzoPSF which keeps positive in XMIPP system.
         */

        double Z = Zo + DeltaZo - k*dzoPSF;

        switch (type)
        {
        case IDEAL_FRESNEL_LENS:
            {
                generatePSFIdealLens(PSFi, Z);
                break;
            }
        case ANALYTIC_ZP:
            {
                REPORT_ERROR(ERR_NOT_IMPLEMENTED, "Zone plates algorithm not implemented yet.");
                break;
            }
        }

        CenterFFT(PSFi, false);
        mPsfVol.setSlice(k, PSFi);

        progress_bar(k2+1);
    }

    /* When PSF is generated from IDEAL_FRESNEL_LENS, dimensions are adjusted to avoid subsampling,
     * in that case it is not recommended to crop XY plane to avoid the loss of normalization in PSF plane. */

    //    if (Nix != Nox || Niy != Noy)
    //        mPsfVol.selfWindow(STARTINGZ(mPsfVol),
    //                           FIRST_XMIPP_INDEX(Noy),FIRST_XMIPP_INDEX(Nox),
    //                           FINISHINGZ(mPsfVol),
    //                           FIRST_XMIPP_INDEX(Noy)+Noy-1,FIRST_XMIPP_INDEX(Nox)+Nox-1);
}

/* Generate the PSF for a ideal lens --------------------------------------- */
void XRayPSF::generatePSFIdealLens(MultidimArray<double> &PSFi, double Zpos) const
{
    double focalEquiv = 1/(1/Zpos - 1/Zo); // inverse of defocus = 1/Z - 1/Zo

    MultidimArray< std::complex<double> > OTFTemp(Niy, Nix), PSFiTemp;
    //    Mask_Params mask_prm; TODO do we have to include masks using this class?

    lensPD(OTFTemp, focalEquiv, lambda, dxl, dyl);

    // Apply the Shape mask
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(OTFTemp)
    {
        if (dAi(*mask,n) == 0)
            dAi(OTFTemp,n) = 0;
    }

    //#define DEBUG
#ifdef DEBUG
    Image<double> _Im;
    _Im().resize(OTFTemp);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(OTFTemp)
    dAij(_Im(),i,j) = arg(dAij(OTFTemp,i,j));
    _Im.write("phase_lens.spi");
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(OTFTemp)
    dAij(_Im(),i,j) = abs(dAij(OTFTemp,i,j));
    _Im.write("abs_lens.spi");
#endif
#undef DEBUG

    FourierTransformer transformer(FFTW_BACKWARD);
    transformer.FourierTransform(OTFTemp, PSFiTemp, false);
    double norm=0, iNorm;

    PSFi.resizeNoCopy(PSFiTemp);

#ifdef DEBUG

    _Im.data.resizeNoCopy(PSFiTemp);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(PSFiTemp)
    DIRECT_MULTIDIM_ELEM(_Im.data,n) = abs(DIRECT_MULTIDIM_ELEM(PSFiTemp,n));
    _Im.write("psfitemp.spi");
#endif

    double aux, aux2;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(PSFi)
    {
        aux=abs(DIRECT_MULTIDIM_ELEM(PSFiTemp,n));
        DIRECT_MULTIDIM_ELEM(PSFi,n)=aux2=aux*aux;
        norm += aux2;
    }

    iNorm = 1/norm;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(PSFi)
    DIRECT_MULTIDIM_ELEM(PSFi,n) *= iNorm;

    //#define DEBUG
#ifdef DEBUG

    _Im() = PSFi;
    CenterFFT(_Im(),0);
    _Im.write("generate-psfi.spi");
#endif
#undef DEBUG


#ifdef DEBUG

    transformer.inverseFourierTransform();

    CenterOriginFFT(OTFTemp,0);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(OTFTemp)
    {
        DIRECT_MULTIDIM_ELEM(_Im.data,n) = abs(DIRECT_MULTIDIM_ELEM(OTFTemp,n));
    }
    _Im.write("otftemp.spi");
#endif

#undef DEBUG
}

void XRayPSF::adjustParam(MultidimArray<double> &vol)
{
    Nox = XSIZE(vol);
    Noy = YSIZE(vol);
    Noz = ZSIZE(vol);

    adjustParam();
}

void XRayPSF::adjustParam()
{
    Nix = Nox;
    Niy = Noy;
    AdjustType = PSFXR_STD;

    if (mode == GENERATE_PSF)
    {
        switch (type)
        {
        case IDEAL_FRESNEL_LENS:
            {
                // Reset the image plane sampling value
                dxi = dxo * Ms;

                dxl = lambda*Zi / (Nox * dxi); // Pixel X-size en the plane of lens aperture
                dyl = lambda*Zi / (Noy * dxi); // Pixel Y-size en the plane of lens aperture

                deltaZMaxX = Zo*((Rlens*2*dxl)/(Rlens*2*dxl - Zo*lambda) - 1);
                deltaZMaxY = Zo*((Rlens*2*dyl)/(Rlens*2*dyl - Zo*lambda) - 1);
                deltaZMinX = Zo*((Rlens*2*dxl)/(Rlens*2*dxl + Zo*lambda) - 1);
                deltaZMinY = Zo*((Rlens*2*dyl)/(Rlens*2*dyl + Zo*lambda) - 1);


                if (verbose > 1)
                {
                    std::cout << std::endl;
                    std::cout << "----------------------" << std::endl;
                    std::cout << "XrayPSF::Param adjust:" << std::endl;
                    std::cout << "----------------------" << std::endl;
                    std::cout << "(Nox,Noy,Nz) = (" << Nox << "," << Noy << "," << Noz << ")" << std::endl;
                    std::cout << "Larger volume Z plane to be calculated = " << (ABS(DeltaZo)+(Noz-1)/2*dzo)*1e6 << " um" << std::endl;
                    std::cout << "Larger allowed discrete Z plane (x,y) = (" << deltaZMaxX*1e6 << ", " << deltaZMaxY*1e6 << ") um" << std::endl;
                }

                if (dxi>dxiMax) /// Lens Radius in pixels higher than image
                {
                    Nix = (size_t)ceil(Nox * dxi/dxiMax);
                    Niy = (size_t)ceil(Noy * dxi/dxiMax);

                    dxi *= Nox/Nix;

                    AdjustType = PSFXR_INT;

                    if (verbose > 1)
                        std::cout << "XrayPSF: Image plane sampling too small: increasing resolution" << std::endl;

                }
                else
                {
                    if (dxi < pupileSizeMin/Nox * dxiMax)
                    {
                        Nix = (size_t)ceil(pupileSizeMin * dxiMax/dxi);
                        AdjustType = PSFXR_ZPAD;
                    }
                    if (dxi < pupileSizeMin/Noy * dxiMax)
                    {
                        Niy = (size_t)ceil(pupileSizeMin * dxiMax/dxi);
                        AdjustType = PSFXR_ZPAD;
                    }
                    if (DeltaZo + Noz/2*dzo > deltaZMaxX)
                    {
                        Nix = std::max(Nix,(size_t)ceil(Zi*Rlens*2*fabs(DeltaZo+Noz/2*dzo)/(Zo*dxi*(Zo+DeltaZo+Noz/2*dzo))));
                        AdjustType = PSFXR_ZPAD;
                    }
                    if (DeltaZo - Noz/2*dzo < deltaZMinX)
                    {
                        Nix = std::max(Nix,(size_t)ceil(Zi*Rlens*2*fabs(DeltaZo-Noz/2*dzo)/(Zo*dxi*(Zo+DeltaZo-Noz/2*dzo))));
                        AdjustType = PSFXR_ZPAD;
                    }
                    if (DeltaZo + Noz/2*dzo > deltaZMaxY)
                    {
                        Niy = std::max(Niy,(size_t)ceil(Zi*Rlens*2*fabs(DeltaZo+Noz/2*dzo)/(Zo*dxi*(Zo+DeltaZo+Noz/2*dzo))));
                        AdjustType = PSFXR_ZPAD;
                    }
                    if ( DeltaZo - Noz/2*dzo < deltaZMinY)
                    {
                        Niy = std::max(Niy,(size_t)ceil(Zi*Rlens*2*fabs(DeltaZo-Noz/2*dzo)/(Zo*dxi*(Zo+DeltaZo-Noz/2*dzo))));
                        AdjustType = PSFXR_ZPAD;
                    }
                }

                dxl = lambda*Zi / (Nix * dxi); // Pixel X-size en the plane of lens aperture
                dyl = lambda*Zi / (Niy * dxi); // Pixel Y-size en the plane of lens aperture

                if (verbose > 1)
                {
                    if (AdjustType!=PSFXR_STD)
                    {
                        if (AdjustType==PSFXR_ZPAD)
                            std::cout << "XrayPSF: Image plane size too small: increasing size" << std::endl;
                        std::cout << "New slice dimensions:  " <<  "(Nix, Niy) = (" << Nix << ", " << Niy << ")" << std::endl;
                    }
                    std::cout << "(dxl, dyl) = (" << dxl*1e6 << ", " << dyl*1e6 << ") um" << std::endl;
                    std::cout << "Pupile Diameter in pixels (NpX, NpY) = (" << ceil(2*Rlens / dxl)
                    << ", " << ceil(2*Rlens / dyl) << ")"  << std::endl;
                    std::cout << std::endl;
                }

                // Calculate the mask to be applied when generating PSFIdealLens

                mask->initZeros(Niy,Nix);

                double Rlens2=Rlens*Rlens;
                double auxY = dyl*(1 - (int)Niy);
                double auxX = dxl*(1 - (int)Nix);

                for (size_t i=0; i<YSIZE(*mask); i++)
                {
                    double y = (double) i * dyl + auxY * 0.5;
                    double y2 = y * y;
                    for (size_t j=0; j<XSIZE(*mask); j++)// Circular mask
                    {
                        /// For indices in standard fashion
                        double x = (double) j * dxl + auxX * 0.5;
                        if (x*x + y2 <= Rlens2)
                            dAij(*mask,i,j) = 1;
                    }
                }


                break;
            }
        case ANALYTIC_ZP:
            {
                REPORT_ERROR(ERR_NOT_IMPLEMENTED, "Zone Plate algorithm not implemented yet.");
                break;
            }
        }
    }
    else
    {
        //        OTF.clear(); // Only for each projection
    }

    if (verbose > 1)
        --verbose;
}

/* Generate the quadratic phase distribution of an ideal lens ------------- */
void lensPD(MultidimArray<std::complex<double> > &Im, double Flens, double lambda, double dx, double dy)
{

    double Lx0 = XSIZE(Im)*dx, Ly0 = YSIZE(Im)*dy, x, y, phase;

    Im.setXmippOrigin();

    double K = (-PI / (lambda * Flens));

    for (size_t i=0; i < YSIZE(Im); ++i)
    {
        y = (double) i * dy + (dy - Ly0) * 0.5;
        double y2 =  y * y;

        for (size_t j=0; j < XSIZE(Im); ++j)
        {
            /// For indices in standard fashion
            x = (double) j * dx + (dx - Lx0) *0.5;

            phase = K * (x * x + y2);

            dAij(Im,i,j) = std::complex<double>(cos(phase),sin(phase));
        }
    }
}
#undef DEBUG
