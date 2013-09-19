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
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "ctf.h"
#include "xmipp_fft.h"
#include <math.h>

/* Read -------------------------------------------------------------------- */
void CTFDescription::readFromMdRow(const MDRow &row, bool disable_if_not_K)
{
    row.getValueOrDefault(MDL_CTF_SAMPLING_RATE, Tm, 1);

    if (enable_CTF)
    {
        row.getValueOrDefault(MDL_CTF_VOLTAGE, kV, 100);
        row.getValueOrDefault(MDL_CTF_DEFOCUSU, DeltafU, 0);
        row.getValueOrDefault(MDL_CTF_DEFOCUSV, DeltafV, DeltafU);
        row.getValueOrDefault(MDL_CTF_DEFOCUS_ANGLE, azimuthal_angle, 0);
        row.getValueOrDefault(MDL_CTF_CS, Cs, 0);
        row.getValueOrDefault(MDL_CTF_CA, Ca, 0);
        row.getValueOrDefault(MDL_CTF_ENERGY_LOSS, espr, 0);
        row.getValueOrDefault(MDL_CTF_LENS_STABILITY, ispr, 0);
        row.getValueOrDefault(MDL_CTF_CONVERGENCE_CONE, alpha, 0);
        row.getValueOrDefault(MDL_CTF_LONGITUDINAL_DISPLACEMENT, DeltaF, 0);
        row.getValueOrDefault(MDL_CTF_TRANSVERSAL_DISPLACEMENT, DeltaR, 0);
        row.getValueOrDefault(MDL_CTF_Q0, Q0, 0);
        row.getValueOrDefault(MDL_CTF_K, K, 1);

        if (K == 0 && disable_if_not_K)
            enable_CTF = false;
    }
    if (enable_CTFnoise)
    {
        row.getValueOrDefault(MDL_CTF_BG_GAUSSIAN_K, gaussian_K, 0);
        row.getValueOrDefault(MDL_CTF_BG_GAUSSIAN_SIGMAU, sigmaU, 0);
        row.getValueOrDefault(MDL_CTF_BG_GAUSSIAN_SIGMAV, sigmaV, sigmaU);
        row.getValueOrDefault(MDL_CTF_BG_GAUSSIAN_CU, cU, 0);
        row.getValueOrDefault(MDL_CTF_BG_GAUSSIAN_CV, cV, cU);
        row.getValueOrDefault(MDL_CTF_BG_GAUSSIAN_ANGLE, gaussian_angle, 0);
        row.getValueOrDefault(MDL_CTF_BG_SQRT_K, sqrt_K, 0);
        row.getValueOrDefault(MDL_CTF_BG_SQRT_U, sqU, 0);
        row.getValueOrDefault(MDL_CTF_BG_SQRT_V, sqV, sqU);
        row.getValueOrDefault(MDL_CTF_BG_SQRT_ANGLE, sqrt_angle, 0);
        row.getValueOrDefault(MDL_CTF_BG_BASELINE, base_line, 0);
        row.getValueOrDefault(MDL_CTF_BG_GAUSSIAN2_K, gaussian_K2, 0);
        row.getValueOrDefault(MDL_CTF_BG_GAUSSIAN2_SIGMAU, sigmaU2, 0);
        row.getValueOrDefault(MDL_CTF_BG_GAUSSIAN2_SIGMAV, sigmaV2, sigmaU2);
        row.getValueOrDefault(MDL_CTF_BG_GAUSSIAN2_CU, cU2, 0);
        row.getValueOrDefault(MDL_CTF_BG_GAUSSIAN2_CV, cV2, cU2);
        row.getValueOrDefault(MDL_CTF_BG_GAUSSIAN2_ANGLE, gaussian_angle2, 0);

        if (gaussian_K == 0 && sqrt_K == 0 && base_line == 0 && gaussian_K2 == 0 &&
            disable_if_not_K)
            enable_CTFnoise = false;
    }
}

/* Read -------------------------------------------------------------------- */
void CTFDescription::readFromMetadataRow(const MetaData &md, size_t id, bool disable_if_not_K)
{
    MDRow row;
    md.getRow(row, id);
    readFromMdRow(row, disable_if_not_K);
}

void CTFDescription::read(const FileName &fn, bool disable_if_not_K)
{
    if (fn.isMetaData())
    {
        MetaData md;
        md.read(fn);
        MDRow row;
        md.getRow(row, md.firstObject());
        readFromMdRow(row, disable_if_not_K);
    }
}

/* Write ------------------------------------------------------------------- */
void CTFDescription::write(const FileName &fn)
{

    MDRow row;
    row.setValue(MDL_CTF_SAMPLING_RATE, Tm);
    if (enable_CTF)
    {
        row.setValue(MDL_CTF_VOLTAGE, kV);
        row.setValue(MDL_CTF_DEFOCUSU, DeltafU);
        row.setValue(MDL_CTF_DEFOCUSV, DeltafV);
        row.setValue(MDL_CTF_DEFOCUS_ANGLE, azimuthal_angle);
        row.setValue(MDL_CTF_CS, Cs);
        row.setValue(MDL_CTF_CA, Ca);
        row.setValue(MDL_CTF_ENERGY_LOSS, espr);
        row.setValue(MDL_CTF_LENS_STABILITY, ispr);
        row.setValue(MDL_CTF_CONVERGENCE_CONE, alpha);
        row.setValue(MDL_CTF_LONGITUDINAL_DISPLACEMENT, DeltaF);
        row.setValue(MDL_CTF_TRANSVERSAL_DISPLACEMENT, DeltaR);
        row.setValue(MDL_CTF_Q0, Q0);
        row.setValue(MDL_CTF_K, K);
    }
    if (enable_CTFnoise)
    {
        row.setValue(MDL_CTF_BG_GAUSSIAN_K, gaussian_K);
        row.setValue(MDL_CTF_BG_GAUSSIAN_SIGMAU, sigmaU);
        row.setValue(MDL_CTF_BG_GAUSSIAN_SIGMAV, sigmaV);
        row.setValue(MDL_CTF_BG_GAUSSIAN_CU, cU);
        row.setValue(MDL_CTF_BG_GAUSSIAN_CV, cV);
        row.setValue(MDL_CTF_BG_GAUSSIAN_ANGLE, gaussian_angle);
        row.setValue(MDL_CTF_BG_SQRT_K, sqrt_K);
        row.setValue(MDL_CTF_BG_SQRT_U, sqU);
        row.setValue(MDL_CTF_BG_SQRT_V, sqV);
        row.setValue(MDL_CTF_BG_SQRT_ANGLE, sqrt_angle);
        row.setValue(MDL_CTF_BG_BASELINE, base_line);
        row.setValue(MDL_CTF_BG_GAUSSIAN2_K, gaussian_K2);
        row.setValue(MDL_CTF_BG_GAUSSIAN2_SIGMAU, sigmaU2);
        row.setValue(MDL_CTF_BG_GAUSSIAN2_SIGMAV, sigmaV2);
        row.setValue(MDL_CTF_BG_GAUSSIAN2_CU, cU2);
        row.setValue(MDL_CTF_BG_GAUSSIAN2_CV, cV2);
        row.setValue(MDL_CTF_BG_GAUSSIAN2_ANGLE, gaussian_angle2);
    }
    if (isLocalCTF)
    {
        row.setValue(MDL_CTF_X0, x0);
        row.setValue(MDL_CTF_XF, xF);
        row.setValue(MDL_CTF_Y0, y0);
        row.setValue(MDL_CTF_YF, yF);
    }

    MetaData md;
    md.setColumnFormat(false);
    md.addRow(row);
    md.write(fn);
}

/* Define Params ------------------------------------------------------------------- */
void CTFDescription::defineParams(XmippProgram * program)
{
    program->addParamsLine("== CTF description");
    program->addParamsLine("  [--ctf_similar_to++ <ctfFile>]        : ctfparam file");
    program->addParamsLine("                                        : Parameters from this file are ");
    program->addParamsLine("                                        : overriden by the parameters in the command line");
    program->addParamsLine("  [--sampling_rate <Tm>]                : Angstroms/pixel. Ex: 1.4");
    program->addParamsLine("                                        : This parameter is compulsory if a CTF is needed.");
    program->addParamsLine("  [--voltage <kV>]                      : Accelerating voltage (kV). Ex: 200");
    program->addParamsLine("                                        : This parameter is compulsory if a CTF is needed.");
    program->addParamsLine("     alias --kV;");
    program->addParamsLine("  [--spherical_aberration <Cs>]         : Milimiters. Ex: 5.6");
    program->addParamsLine("                                        : This parameter is compulsory if a CTF is needed.");
    program->addParamsLine("     alias --Cs;");
    program->addParamsLine("  [--defocusU <DeltafU>]                : Defocus in Angstroms (Ex: 2000)");
    program->addParamsLine("  [--defocusV++ <DeltafV>]              : If astigmatic");
    program->addParamsLine("  [--azimuthal_angle++ <ang=0>]         : Angle between X and U (degrees)");
    program->addParamsLine("  [--chromatic_aberration++ <Ca=0>]     : Milimiters. Ex: 2");
    program->addParamsLine("  [--energy_loss++ <espr=0>]            : eV. Ex: 1");
    program->addParamsLine("  [--lens_stability++ <ispr=0>]         : ppm. Ex: 1");
    program->addParamsLine("  [--convergence_cone++ <alpha=0>]      : mrad. Ex: 0.5");
    program->addParamsLine("  [--longitudinal_displace++ <DeltaF=0>]: Angstrom. Ex: 100");
    program->addParamsLine("  [--transversal_displace++ <DeltaR=0>] : Angstrom. Ex: 3");
    program->addParamsLine("  [--Q0++ <Q0=0>]                       : Percentage of cosine (Q0>0)");
    program->addParamsLine("  [--K++ <K=1>]                         : Global gain");
}

/* Read from command line -------------------------------------------------- */
void CTFDescription::readParams(XmippProgram * program)
{
    kV=Tm=Cs=0;
    if (program->checkParam("--ctf_similar_to"))
        read(program->getParam("--ctf_similar_to"));
    if (program->checkParam("--sampling_rate"))
        Tm=program->getDoubleParam("--sampling_rate");
    if (Tm==0)
        REPORT_ERROR(ERR_ARG_MISSING,"--sampling_rate");
    if (program->checkParam("--voltage"))
        kV=program->getDoubleParam("--voltage");
    if (kV==0)
        REPORT_ERROR(ERR_ARG_MISSING,"--voltage");
    if (program->checkParam("--spherical_aberration"))
        Cs=program->getDoubleParam("--spherical_aberration");
    if (Cs==0)
        REPORT_ERROR(ERR_ARG_MISSING,"--spherical_aberration");
    if (program->checkParam("--defocusU"))
        DeltafU=program->getDoubleParam("--defocusU");
    if (program->checkParam("--defocusV"))
        DeltafV=program->getDoubleParam("--defocusV");
    else
        DeltafV=DeltafU;
    if (program->checkParam("--Q0"))
        Q0=program->getDoubleParam("--Q0");
    azimuthal_angle=program->getDoubleParam("--azimuthal_angle");
    Ca=program->getDoubleParam("--chromatic_aberration");
    espr=program->getDoubleParam("--energy_loss");
    ispr=program->getDoubleParam("--lens_stability");
    alpha=program->getDoubleParam("--convergence_cone");
    DeltaF=program->getDoubleParam("--longitudinal_displace");
    DeltaR=program->getDoubleParam("--transversal_displace");
    Q0=program->getDoubleParam("--Q0");
    K=program->getDoubleParam("--K");
}

/* Show -------------------------------------------------------------------- */
std::ostream & operator << (std::ostream &out, const CTFDescription &ctf)
{
    if (ctf.enable_CTF)
    {
        out
        << "sampling_rate=        " << ctf.Tm              << std::endl
        << "voltage=              " << ctf.kV              << std::endl
        << "defocusU=             " << ctf.DeltafU         << std::endl
        << "defocusV=             " << ctf.DeltafV         << std::endl
        << "azimuthal_angle=      " << ctf.azimuthal_angle << std::endl
        << "spherical_aberration= " << ctf.Cs              << std::endl
        << "chromatic_aberration= " << ctf.Ca              << std::endl
        << "energy_loss=          " << ctf.espr            << std::endl
        << "lens_stability=       " << ctf.ispr            << std::endl
        << "convergence_cone=     " << ctf.alpha           << std::endl
        << "longitudinal_displace=" << ctf.DeltaF          << std::endl
        << "transversal_displace= " << ctf.DeltaR          << std::endl
        << "Q0=                   " << ctf.Q0              << std::endl
        << "K=                    " << ctf.K               << std::endl
        ;
    }
    if (ctf.enable_CTFnoise)
    {
        out
        << "gaussian_K=           " << ctf.gaussian_K      << std::endl
        << "sigmaU=               " << ctf.sigmaU          << std::endl
        << "sigmaV=               " << ctf.sigmaV          << std::endl
        << "cU=                   " << ctf.cU              << std::endl
        << "cV=                   " << ctf.cV              << std::endl
        << "gaussian_angle=       " << ctf.gaussian_angle  << std::endl
        << "sqrt_K=               " << ctf.sqrt_K          << std::endl
        << "sqU=                  " << ctf.sqU             << std::endl
        << "sqV=                  " << ctf.sqV             << std::endl
        << "sqrt_angle=           " << ctf.sqrt_angle      << std::endl
        << "base_line=            " << ctf.base_line       << std::endl
        << "gaussian_K2=          " << ctf.gaussian_K2     << std::endl
        << "sigmaU2=              " << ctf.sigmaU2         << std::endl
        << "sigmaV2=              " << ctf.sigmaV2         << std::endl
        << "cU2=                  " << ctf.cU2             << std::endl
        << "cV2=                  " << ctf.cV2             << std::endl
        << "gaussian_angle2=      " << ctf.gaussian_angle2 << std::endl
        ;
    }
    return out;
}

/* Default values ---------------------------------------------------------- */
void CTFDescription::clear()
{
    enable_CTF = true;
    enable_CTFnoise = false;
    isLocalCTF = false;
    clearNoise();
    clearPureCtf();
    y0=x0=xF=yF=0;
}

void CTFDescription::clearNoise()
{
    base_line = 0;
    cU = cV = sigmaU = sigmaV = gaussian_angle = gaussian_K = 0;
    sqU = sqV = sqrt_K = sqrt_angle = 0;
    cU2 = cV2 = sigmaU2 = sigmaV2 = gaussian_angle2 = gaussian_K2 = 0;
    isLocalCTF = false;
}

void CTFDescription::clearPureCtf()
{
    enable_CTF = true;
    enable_CTFnoise = false;
    Tm = 2;
    kV = 100;
    DeltafU = DeltafV = azimuthal_angle = 0;
    Cs = Ca = espr = ispr = alpha = DeltaF = DeltaR = 0;
    K = 1;
    Q0 = 0;
    isLocalCTF = false;
}

/* Produce Side Information ------------------------------------------------ */
void CTFDescription::produceSideInfo()
{
    // Change units
    double local_Cs = Cs * 1e7;
    double local_Ca = Ca * 1e7;
    double local_kV = kV * 1e3;
    double local_ispr = ispr * 1e6;
    rad_azimuth = DEG2RAD(azimuthal_angle);
    rad_gaussian = DEG2RAD(gaussian_angle);
    rad_gaussian2 = DEG2RAD(gaussian_angle2);
    rad_sqrt = DEG2RAD(sqrt_angle);

    // lambda=h/sqrt(2*m*e*kV)
    //    h: Planck constant
    //    m: electron mass
    //    e: electron charge
    lambda=12.2643247/sqrt(local_kV*(1.+0.978466e-6*local_kV)); // See http://en.wikipedia.org/wiki/Electron_diffraction
    //
    defocus_average   = -(DeltafU + DeltafV) * 0.5;
    defocus_deviation = -(DeltafU - DeltafV) * 0.5;
    // Phase shift for spherical aberration
    // X(u)=-PI*deltaf(u)*lambda*u^2+PI/2*Cs*lambda^3*u^4
    // ICE: X(u)=-PI/2*deltaf(u)*lambda*u^2+PI/2*Cs*lambda^3*u^4
    //          = K1*deltaf(u)*u^2         +K2*u^4
    K1 = PI / 2 * 2 * lambda;
    K2 = PI / 2 * local_Cs * lambda * lambda * lambda;

    // Envelope
    // D(u)=Ed(u)*Ealpha(u)
    // Ed(u)=exp(-1/2*PI^2*lambda^2*D^2*u^4)
    // Ealpha(u)=exp(-PI^2*alpha^2*u^2*(Cs*lambda^2*u^2+Deltaf(u))^2)
    // ICE: Eespr(u)=exp(-(1/4*PI*Ca*lambda*espr/kV)^2*u^4/log2)
    // ICE: Eispr(u)=exp(-(1/2*PI*Ca*lambda*ispr)^2*u^4/log2)
    // ICE: EdeltaF(u)=bessj0(PI*DeltaF*lambda*u^2)
    // ICE: EdeltaR(u)=sinc(u*DeltaR)
    // ICE: Ealpha(u)=exp(-PI^2*alpha^2*(Cs*lambda^2*u^3+Deltaf(u)*u)^2)
    // CO: K3=pow(0.25*PI*Ca*lambda*(espr/kV,2)/log(2); Both combines in new K3
    // CO: K4=pow(0.5*PI*Ca*lambda*ispr,2)/log(2);
    K3 = pow(0.25 * PI * local_Ca * lambda * (espr / kV + 2 * local_ispr), 2) / log(2.0);
    K5 = PI * DeltaF * lambda;
    K6 = PI * PI * alpha * alpha;
    K7 = local_Cs * lambda * lambda;
    Ksin = sqrt(1-Q0*Q0);
    Kcos = Q0;
}

/* Precompute values ------------------------------------------------------- */
void CTFDescription::precomputeValues(const MultidimArray<double> &cont_x_freq,
                                      const MultidimArray<double> &cont_y_freq)
{
    precomputedImage.reserve(MULTIDIM_SIZE(cont_x_freq));
    precomputedImageXdim=XSIZE(cont_x_freq);

    FOR_ALL_ELEMENTS_IN_ARRAY2D(cont_x_freq)
    {
        double X=A2D_ELEM(cont_x_freq,i,j);
        double Y=A2D_ELEM(cont_y_freq,i,j);
        precomputeValues(X, Y);
        if (fabs(X) < XMIPP_EQUAL_ACCURACY &&
            fabs(Y) < XMIPP_EQUAL_ACCURACY)
            precomputed.deltaf=0;
        else
            precomputed.deltaf=-1;
        precomputedImage.push_back(precomputed);
    }
}

/* Zero -------------------------------------------------------------------- */
//#define DEBUG
void CTFDescription::zero(int n, const Matrix1D<double> &u, Matrix1D<double> &freq)
{
    double wmax = 1 / (2 * Tm);
    double wstep = wmax / 300;
    int sign_changes = 0;
    double last_ctf = getValuePureNoPrecomputedAt(0,0), ctf;
    double w;
    for (w = 0; w <= wmax; w += wstep)
    {
        V2_BY_CT(freq, u, w);
        ctf = getValuePureNoPrecomputedAt(XX(freq), YY(freq));
        if (SGN(ctf) != SGN(last_ctf))
        {
            sign_changes++;
            if (sign_changes == n)
                break;
        }
        last_ctf = ctf;
    }
    if (sign_changes != n)
    {
        VECTOR_R2(freq, -1, -1);
    }
    else
    {
        // Compute more accurate zero
#ifdef DEBUG
        std::cout << n << " zero: w=" << w << " (" << wmax << ") freq="
        << (u*w).transpose()
        << " last_ctf=" << last_ctf << " ctf=" << ctf << " ";
#endif

        w += ctf * wstep / (last_ctf - ctf);
        V2_BY_CT(freq, u, w);
#ifdef DEBUG

        std::cout << " final w= " << w << " final freq=" << freq.transpose() << std::endl;
#endif

    }
}
#undef DEBUG

/* Apply the CTF to an image ----------------------------------------------- */
void CTFDescription::applyCTF(MultidimArray < std::complex<double> > &FFTI)
{
    Matrix1D<int>    idx(2);
    Matrix1D<double> freq(2);
    if ( ZSIZE(FFTI) > 1 )
        REPORT_ERROR(ERR_MULTIDIM_DIM,"ERROR: Apply_CTF only works on 2D images, not 3D.");

    FOR_ALL_ELEMENTS_IN_ARRAY2D(FFTI)
    {
        XX(idx) = j;
        YY(idx) = i;
        FFT_idx2digfreq(FFTI, idx, freq);
        precomputeValues(XX(freq), YY(freq));
        double ctf = getValueAt();
        FFTI(i, j) *= ctf;
    }
}

/* Get profiles ------------------------------------------------------------ */
void CTFDescription::getProfile(double angle, double fmax, int nsamples,
                                MultidimArray<double> &profiles)
{
    double step = fmax / nsamples;

    profiles.resizeNoCopy(nsamples, 4);
    double sinus = sin(angle);
    double cosinus = cos(angle);
    double f;
    size_t i;

    for (i = 0, f = 0; i < YSIZE(profiles); i++, f += step)
    {
        double fx = f * cosinus;
        double fy = f * sinus;

        // Compute current frequencies.
        precomputeValues(fx, fy);

        // Store values.
        double bgNoise = getValueNoiseAt();
        double ctf = getValuePureAt();
        double E = getValueDampingAt();

        A2D_ELEM(profiles, i, 0) = bgNoise;
        A2D_ELEM(profiles, i, 1) = bgNoise + E * E;
        A2D_ELEM(profiles, i, 2) = bgNoise + ctf * ctf;
        A2D_ELEM(profiles, i, 3) = ctf;
    }
}

/* Get average profiles ----------------------------------------------------- */
void CTFDescription::getAverageProfile(double fmax, int nsamples,
                                       MultidimArray<double> &profiles)
{
    double step = fmax / nsamples;
    profiles.initZeros(nsamples, 4);

    for (double angle = 0.0; angle < 360; angle++)
    {
        double sinus = sin(angle);
        double cosinus = cos(angle);
        double f;
        size_t i;
        for (i = 0, f = 0; i < YSIZE(profiles); i++, f += step)
        {
            double fx = f * cosinus;
            double fy = f * sinus;

            // Compute current frequencies.
            precomputeValues(fx, fy);

            // Store values.
            double bgNoise = getValueNoiseAt();
            double ctf = getValuePureAt();
            double E = getValueDampingAt();

            A2D_ELEM(profiles, i, 0) += bgNoise;
            A2D_ELEM(profiles, i, 1) += bgNoise + E * E;
            A2D_ELEM(profiles, i, 2) += bgNoise + ctf * ctf;
            A2D_ELEM(profiles, i, 3) += ctf;
        }
    }
    profiles*=1.0/360;
}

/* Generate CTF Image ------------------------------------------------------ */
//#define DEBUG
void CTFDescription::generateCTF(int Ydim, int Xdim,
                                 MultidimArray < std::complex<double> > &CTF)
{
    Matrix1D<int>    idx(2);
    Matrix1D<double> freq(2);
    CTF.resize(Ydim, Xdim);
#ifdef DEBUG

    std::cout << "CTF:\n" << *this << std::endl;
#endif

    FOR_ALL_ELEMENTS_IN_ARRAY2D(CTF)
    {
        XX(idx) = j;
        YY(idx) = i;
        FFT_idx2digfreq(CTF, idx, freq);
        digfreq2contfreq(freq, freq, Tm);
        precomputeValues(XX(freq), YY(freq));
        //FIXME: Check why not to use the macro
        CTF(i, j) = getValueAt();
#ifdef DEBUG

        if (i == 0)
            std::cout << i << " " << j << " " << YY(freq) << " " << XX(freq)
            << " " << CTF(i, j) << std::endl;
#endif

    }
}
#undef DEBUG

/* Physical meaning -------------------------------------------------------- */
//#define DEBUG
bool CTFDescription::hasPhysicalMeaning()
{
    bool retval;
    if (enable_CTF)
    {
        precomputeValues(0,0);
        retval =
            K >= 0       && base_line >= 0  &&
            kV >= 50     && kV <= 1000      &&
            espr >= 0    && espr <= 20      &&
            ispr >= 0    && ispr <= 20      &&
            Cs >= 0      && Cs <= 20        &&
            Ca >= 0      && Ca <= 3         &&
            alpha >= 0   && alpha <= 5      &&
            DeltaF >= 0  && DeltaF <= 1000  &&
            DeltaR >= 0  && DeltaR <= 100   &&
            Q0 >= 0      && Q0 <= 0.40      &&
            DeltafU >= 0 && DeltafV >= 0    &&
            getValueAt() >= 0;
#ifdef DEBUG

        if (retval == false)
        {
            std::cout << *this << std::endl;
            std::cout << "K>=0       && base_line>=0  " << (K >= 0       && base_line >= 0) << std::endl
            << "kV>=50     && kV<=1000      " << (kV >= 50     && kV <= 1000)     << std::endl
            << "espr>=0    && espr<=20      " << (espr >= 0    && espr <= 20)     << std::endl
            << "ispr>=0    && ispr<=20      " << (ispr >= 0    && ispr <= 20)     << std::endl
            << "Cs>=0      && Cs<=20        " << (Cs >= 0      && Cs <= 20)       << std::endl
            << "Ca>=0      && Ca<=3         " << (Ca >= 0      && Ca <= 3)        << std::endl
            << "alpha>=0   && alpha<=5      " << (alpha >= 0   && alpha <= 5)     << std::endl
            << "DeltaF>=0  && DeltaF<=1000  " << (DeltaF >= 0  && DeltaF <= 1000) << std::endl
            << "DeltaR>=0  && DeltaR<=100   " << (DeltaR >= 0  && DeltaR <= 100)  << std::endl
            << "Q0>=0      && Q0<=0.4       " << (Q0 >= 0      && Q0 <= 0.4)      << std::endl
            << "DeltafU>=0 && DeltafV>=0    " << (DeltafU >= 0 && DeltafV >= 0)   << std::endl
            << "getValueAt(0,0)>=0       " << (getValueAt() >= 0)         << std::endl
            ;
            std::cout << "getValueAt(0,0)=" << getValueAt(true) << std::endl;
        }
#endif

    }
    else
        retval = true;
    bool retval2;
    if (enable_CTFnoise)
    {
        double min_sigma = XMIPP_MIN(sigmaU, sigmaV);
        double min_c = XMIPP_MIN(cU, cV);
        double min_sigma2 = XMIPP_MIN(sigmaU2, sigmaV2);
        double min_c2 = XMIPP_MIN(cU2, cV2);
        retval2 =
            base_line >= 0       &&
            gaussian_K >= 0      &&
            sigmaU >= 0          && sigmaV >= 0          &&
            sigmaU <= 100e3      && sigmaV <= 100e3       &&
            cU >= 0              && cV >= 0              &&
            sqU >= 0             && sqV >= 0             &&
            sqrt_K >= 0          &&
            gaussian_K2 >= 0     &&
            sigmaU2 >= 0         && sigmaV2 >= 0          &&
            sigmaU2 <= 100e3     && sigmaV2 <= 100e3      &&
            cU2 >= 0             && cV2 >= 0              &&
            gaussian_angle >= 0  && gaussian_angle <= 90  &&
            sqrt_angle >= 0      && sqrt_angle <= 90      &&
            gaussian_angle2 >= 0 && gaussian_angle2 <= 90
            ;
        if (min_sigma > 0)
            retval2 = retval2 && ABS(sigmaU - sigmaV) / min_sigma <= 3;
        if (min_c > 0)
            retval2 = retval2 && ABS(cU - cV) / min_c <= 3;
        if (gaussian_K != 0)
            retval2 = retval2 && (cU * Tm >= 0.01) && (cV * Tm >= 0.01);
        if (min_sigma2 > 0)
            retval2 = retval2 && ABS(sigmaU2 - sigmaV2) / min_sigma2 <= 3;
        if (min_c2 > 0)
            retval2 = retval2 && ABS(cU2 - cV2) / min_c2 <= 3;
        if (gaussian_K2 != 0)
            retval2 = retval2 && (cU2 * Tm >= 0.01) && (cV2 * Tm >= 0.01);
#ifdef DEBUG

        if (retval2 == false)
        {
            std::cout << *this << std::endl;
            std::cout << "base_line>=0       &&        " << (base_line >= 0)           << std::endl
            << "gaussian_K>=0      &&        " << (gaussian_K >= 0)    << std::endl
            << "sigmaU>=0      && sigmaV>=0     " << (sigmaU >= 0      && sigmaV >= 0)   << std::endl
            << "sigmaU<=100e3      && sigmaV<=100e3       " << (sigmaU <= 100e3      && sigmaV <= 100e3)   << std::endl
            << "cU>=0       && cV>=0      " << (cU >= 0       && cV >= 0)   << std::endl
            << "sqU>=0      && sqV>=0      " << (sqU >= 0        && sqV >= 0)   << std::endl
            << "sqrt_K>=0      &&        " << (sqrt_K >= 0)     << std::endl
            << "gaussian_K2>=0     &&        " << (gaussian_K2 >= 0)    << std::endl
            << "sigmaU2>=0      && sigmaV2>=0     " << (sigmaU2 >= 0      && sigmaV2 >= 0)   << std::endl
            << "sigmaU2<=100e3     && sigmaV2<=100e3      " << (sigmaU2 <= 100e3     && sigmaV2 <= 100e3)  << std::endl
            << "cU2>=0       && cV2>=0      " << (cU2 >= 0      && cV2 >= 0)  << std::endl
            << "gaussian_angle>=0  && gaussian_angle<=90  " << (gaussian_angle >= 0  && gaussian_angle <= 90)  << std::endl
            << "sqrt_angle>=0      && sqrt_angle<=90      " << (sqrt_angle >= 0      && sqrt_angle <= 90)      << std::endl
            << "gaussian_angle2>=0 && gaussian_angle2<=90 " << (gaussian_angle2 >= 0 && gaussian_angle2 <= 90) << std::endl
            ;
            if (min_sigma > 0)
                std::cout << "ABS(sigmaU-sigmaV)/min_sigma<=3         " << (ABS(sigmaU - sigmaV) / min_sigma <= 3)     << std::endl;
            if (min_c > 0)
                std::cout << "ABS(cU-cV)/min_c<=3                     " << (ABS(cU - cV) / min_c <= 3)                 << std::endl;
            if (gaussian_K > 0)
                std::cout << "(cU*Tm>=0.01) && (cV*Tm>=0.01)          " << ((cU*Tm >= 0.01) && (cV*Tm >= 0.01))      << std::endl;
            if (min_sigma2 > 0)
                std::cout << "ABS(sigmaU2-sigmaV2)/min_sigma2<=3      " << (ABS(sigmaU2 - sigmaV2) / min_sigma2 <= 3)  << std::endl;
            if (min_c2 > 0)
                std::cout << "ABS(cU2-cV2)/min_c2<=3                  " << (ABS(cU2 - cV2) / min_c2 <= 3)              << std::endl;
            if (gaussian_K2 > 0)
                std::cout << "(cU2*Tm>=0.01) && (cV2*Tm>=0.01)        " << ((cU2*Tm >= 0.01) && (cV2*Tm >= 0.01))    << std::endl;
            std::cout << cV2*Tm << std::endl;
        }
#endif

    }
    else
        retval2 = true;
#ifdef DEBUG
    // std::cout << "Retval= " << retval << " retval2= " << retval2 << std::endl;
#endif

    return retval && retval2;
}
#undef DEBUG

/* Force Physical meaning -------------------------------------------------- */
void CTFDescription::forcePhysicalMeaning()
{
    if (enable_CTF)
    {
        if (K < 0)
            K = 0;
        if (base_line < 0)
            base_line = 0;
        if (kV < 50)
            kV = 50;
        if (kV > 1000)
            kV = 1000;
        if (espr < 0)
            espr = 0;
        if (espr > 20)
            espr = 20;
        if (ispr < 0)
            ispr = 0;
        if (ispr > 20)
            ispr = 20;
        if (Cs < 0)
            Cs = 0;
        if (Cs > 20)
            Cs = 20;
        if (Ca < 0)
            Ca = 0;
        if (Ca > 3)
            Ca = 3;
        if (alpha < 0)
            alpha = 0;
        if (alpha > 5)
            alpha = 5;
        if (DeltaF < 0)
            DeltaF = 0;
        if (DeltaF > 1000)
            DeltaF = 1000;
        if (DeltaR < 0)
            DeltaR = 0;
        if (DeltaR > 1000)
            DeltaR = 1000;
        if (Q0 > 0.40)
            Q0 = 0.40;
        if (Q0 < 0)
            Q0 = 0;
        if (DeltafU < 0)
            DeltafU = 0;
        if (DeltafV < 0)
            DeltafV = 0;
    }
    if (enable_CTFnoise)
    {
        double min_sigma = XMIPP_MIN(sigmaU, sigmaV);
        double min_c = XMIPP_MIN(cU, cV);
        double min_sigma2 = XMIPP_MIN(sigmaU2, sigmaV2);
        double min_c2 = XMIPP_MIN(cU2, cV2);
        if (base_line < 0)
            base_line = 0;
        if (gaussian_K < 0)
            gaussian_K = 0;
        if (sigmaU < 0)
            sigmaU = 0;
        if (sigmaV < 0)
            sigmaV = 0;
        if (sigmaU > 100e3)
            sigmaU = 100e3;
        if (sigmaV > 100e3)
            sigmaV = 100e3;
        if (cU < 0)
            cU = 0;
        if (cV < 0)
            cV = 0;
        if (sqU < 0)
            sqU = 0;
        if (sqV < 0)
            sqV = 0;
        if (sqrt_K < 0)
            sqrt_K = 0;
        if (gaussian_K2 < 0)
            gaussian_K2 = 0;
        if (sigmaU2 < 0)
            sigmaU2 = 0;
        if (sigmaV2 < 0)
            sigmaV2 = 0;
        if (sigmaU2 > 100e3)
            sigmaU2 = 100e3;
        if (sigmaV2 > 100e3)
            sigmaV2 = 100e3;
        if (cU2 < 0)
            cU2 = 0;
        if (cV2 < 0)
            cV2 = 0;
        if (gaussian_angle < 0)
            gaussian_angle = 0;
        if (gaussian_angle > 90)
            gaussian_angle = 90;
        if (sqrt_angle < 0)
            sqrt_angle = 0;
        if (sqrt_angle > 90)
            sqrt_angle = 90;
        if (gaussian_angle2 < 0)
            gaussian_angle2 = 0;
        if (gaussian_angle2 > 90)
            gaussian_angle2 = 90;
        if (min_sigma > 0)
            if (ABS(sigmaU - sigmaV) / min_sigma > 3)
            {
                if (sigmaU < sigmaV)
                    sigmaV = 3.9 * sigmaU;
                else
                    sigmaU = 3.9 * sigmaV;
            }
        if (min_c > 0)
            if (ABS(cU - cV) / min_c > 3)
            {
                if (cU < cV)
                    cV = 3.9 * cU;
                else
                    cU = 3.9 * cV;
            }
        if (gaussian_K != 0)
        {
            if (cU*Tm < 0.01)
                cU = 0.011 / Tm;
            if (cV*Tm < 0.01)
                cV = 0.011 / Tm;
        }
        if (min_sigma2 > 0)
            if (ABS(sigmaU2 - sigmaV2) / min_sigma2 > 3)
            {
                if (sigmaU2 < sigmaV2)
                    sigmaV2 = 3.9 * sigmaU2;
                else
                    sigmaU2 = 3.9 * sigmaV2;
            }
        if (min_c2 > 0)
            if (ABS(cU2 - cV2) / min_c2 > 3)
            {
                if (cU2 < cV2)
                    cV2 = 3.9 * cU2;
                else
                    cU2 = 3.9 * cV2;
            }
        if (gaussian_K2 != 0)
        {
            if (cU2*Tm < 0.01)
                cU2 = 0.011 / Tm;
            if (cV2*Tm < 0.01)
                cV2 = 0.011 / Tm;
        }
    }
}
#undef DEBUG

void generateCTFImageWith2CTFs(const MetaData &MD1, const MetaData &MD2, int Xdim, MultidimArray<double> &imgOut)
{
    CTFDescription CTF1, CTF2;
    CTF1.enable_CTF=true;
    CTF1.enable_CTFnoise=false;
    CTF1.readFromMetadataRow(MD1,MD1.firstObject());
    CTF1.produceSideInfo();

    CTF2.enable_CTF=true;
    CTF2.enable_CTFnoise=false;
    CTF2.readFromMetadataRow(MD2,MD2.firstObject());
    CTF2.produceSideInfo();

    imgOut.initZeros(Xdim,Xdim);
    Matrix1D<int> idx(2);
    Matrix1D<double> freq(2);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(imgOut)
    {
        XX(idx) = j;
        YY(idx) = i;

        // Digital frequency
        FFT_idx2digfreq(imgOut, idx, freq);
        if (XX(freq)>=0 && XX(freq)<0.5)
        {
            digfreq2contfreq(freq, freq, CTF1.Tm);
            CTF1.precomputeValues(XX(freq),YY(freq));
            A2D_ELEM(imgOut,i,j)=CTF1.getValueAt();
        }
        else
        {
            digfreq2contfreq(freq, freq, CTF2.Tm);
            CTF2.precomputeValues(XX(freq),YY(freq));
            A2D_ELEM(imgOut,i,j)=CTF2.getValueAt();
        }
    }
    CenterFFT(imgOut,false);
}

double errorBetween2CTFs( MetaData &MD1,
                          MetaData &MD2,
                          size_t Xdim,
                          double minFreq,
                          double maxFreq)
{

    Matrix1D<double> freq(2); // Frequencies for Fourier plane

    CTFDescription CTF1, CTF2;

    CTF1.enable_CTF=true;
    CTF1.enable_CTFnoise=false;
    CTF1.readFromMetadataRow(MD1,MD1.firstObject());
    CTF1.produceSideInfo();

    CTF2.enable_CTF=true;
    CTF2.enable_CTFnoise=false;
    CTF2.readFromMetadataRow(MD2,MD2.firstObject());
    CTF2.produceSideInfo();

    double iTm=1.0/CTF1.Tm;
    size_t xDim, yDim;
    xDim = yDim = Xdim;
    //#define DEBUG
#ifdef DEBUG

    Image<double> img1(xDim,yDim);
    Image<double> img2(xDim,yDim);
    Image<double> img3(xDim,yDim);

    MultidimArray<double> &dummy1 = img1.data;
    MultidimArray<double> &dummy2 = img2.data;
    MultidimArray<double> &dummy3 = img3.data;
#endif

    double error =0.;
    double _freq=0.;

    for (int i=0; i<(int)yDim; ++i)
    {
        FFT_IDX2DIGFREQ(i, yDim, YY(freq));
        YY(freq) *= iTm;
        for (int j=0; j<(int)xDim; ++j)
        {
            FFT_IDX2DIGFREQ(j, xDim, XX(freq));
            _freq = freq.module();
            if (_freq < minFreq || _freq > maxFreq)
                continue;
            XX(freq) *= iTm;
            //freq *= CTF1.Tm;
            CTF1.precomputeValues(XX(freq),YY(freq));
            CTF2.precomputeValues(XX(freq),YY(freq));
            error += fabs(CTF2.getValuePureWithoutDampingAt()- CTF1.getValuePureWithoutDampingAt());
#ifdef DEBUG

            double a = CTF1.getValuePureWithoutDampingAt();
            double b = CTF2.getValuePureWithoutDampingAt();
            DIRECT_A2D_ELEM(dummy1,i,j)=a;
            DIRECT_A2D_ELEM(dummy2,i,j)=b;
            DIRECT_A2D_ELEM(dummy3,i,j)=fabs(b-a);
#endif

        }
    }
#ifdef DEBUG
    img1.write("/tmp/img1.spi");
    img2.write("/tmp/img2.spi");
    img3.write("/tmp/img3.spi");
#endif
#undef DEBUG
    return error;
}


double errorMaxFreqCTFs( MetaData &MD1,
                         MetaData &MD2)
{
    //This is a quick and extremely dirty solution
    //return angstroms

    CTFDescription CTF1, CTF2;

    CTF1.enable_CTF=true;
    CTF1.enable_CTFnoise=false;
    CTF1.readFromMetadataRow(MD1,MD1.firstObject());
    double defocus1 = (CTF1.DeltafU + CTF1.DeltafV)/2.0;
    CTF1.produceSideInfo();

    CTF2.enable_CTF=true;
    CTF2.enable_CTFnoise=false;
    CTF2.readFromMetadataRow(MD2,MD2.firstObject());
    double defocus2 = (CTF2.DeltafU + CTF2.DeltafV)/2.0;
    CTF2.produceSideInfo();
#ifdef DEBUG

    Matrix1D<double> freq(2); // Frequencies for Fourier plane
    double iTm=1.0/CTF1.Tm;
    for (int i=0; i<256; ++i)
    {
        FFT_IDX2DIGFREQ(i, 512, YY(freq));
        FFT_IDX2DIGFREQ(i, 512, XX(freq));
        freq *= iTm;
        CTF1.precomputeValues(XX(freq),YY(freq));
        CTF2.precomputeValues(XX(freq),YY(freq));
        std::cerr << "DEBUG_ROB: "
        		  << XX(freq) << " "
        		  << CTF1.getValuePureWithoutDampingAt() << " "
        		  << CTF2.getValuePureWithoutDampingAt() << " "
        		  << std::endl;
    }
#endif
    return 1.0/sqrt(PI/(CTF1.K1*abs(defocus1-defocus2)));
}

void generatePSDCTFImage(MultidimArray<double> &img, const MetaData &MD)
{
    img.rangeAdjust(-1,1);
    CenterFFT(img,false);

    CTFDescription CTF;
    CTF.enable_CTF=true;
    CTF.enable_CTFnoise=false;
    CTF.readFromMetadataRow(MD,MD.firstObject());
    CTF.produceSideInfo();

    Matrix1D<int> idx(2);
    Matrix1D<double> freq(2);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(img)
    {
        XX(idx) = j;
        YY(idx) = i;

        // Digital frequency
        FFT_idx2digfreq(img, idx, freq);
        if (XX(freq)>=0 && XX(freq)<0.5)
        {
            digfreq2contfreq(freq, freq, CTF.Tm);
            CTF.precomputeValues(XX(freq),YY(freq));
            double aux=CTF.getValueAt();
            A2D_ELEM(img,i,j)=aux*aux;
        }
    }
    CenterFFT(img,true);
}
