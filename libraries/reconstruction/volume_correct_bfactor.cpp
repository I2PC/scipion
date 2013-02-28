/***************************************************************************
 *
 * Author:     Sjors H.W. Scheres               (scheres@cnb.csic.es)
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

#include "volume_correct_bfactor.h"

void ProgVolumeCorrectBfactor::defineParams()
{
    addUsageLine("Sharpen a volume by applying a negative B-factor");
    addUsageLine("+The high-resolution features are enhanced, thereby correcting ");
    addUsageLine("+the envelope functions of the microscope, detector etc. Three ");
    addUsageLine("+modes are available:");
    addUsageLine("+");
    addUsageLine("+1. An automated mode based on methodology developed by [[http://www.ncbi.nlm.nih.gov/pubmed/14568533][Rosenthal and Henderson]]");
    addUsageLine("+");
    addUsageLine("+2. Based on the fall-off in a reference map (possibly obtained using xmipp_volume_from_pdb)");
    addUsageLine("+");
    addUsageLine("+3. An user-provided (ad hoc) B-factor");
    addSeeAlsoLine("volume_enhance_contrast, volume_from_pdb, resolution_fsc");
    XmippMetadataProgram::defineParams();
    addParamsLine(" --auto                              : Use automated B-factor fit in flat Wilson region");
    addParamsLine("                                     : Note: do not use the automated mode for maps with resolutions");
    addParamsLine("                                     : lower than 12-15 Angstroms!");
    addParamsLine(" or --ref <fn_ref> <mode=bfactorref> : Fit B-factor according to the reference ");
    addParamsLine("         where <mode>");
    addParamsLine("                  bfactorref         : Fit B-factor");
    addParamsLine("                  allpoints          : Adjust power spectrum to reference");
    addParamsLine(" or --adhoc <B>              : Use a user-provided (negative) B-factor");
    addParamsLine("== Specific parameters == ");
    addParamsLine("  --sampling <float>         : Pixel size of the input volume (in Angstroms/pixel) ");
    addParamsLine("  --maxres <float>           : High-resolution limit for B-factor correction (in Ang.) ");
    addParamsLine("  -o <filename>              : Output file Name with corrected volume ");
    addParamsLine(" [--fit_minres+ <f=15>]      : Low-resolution  limit (in Ang.) for fit in --auto or --ref ");
    addParamsLine(" [--fit_maxres+ <f=-1>]      : High-resolution limit (in Ang.) for fit in --auto or --ref,");
    addParamsLine("                             : -1 means maximum resolution ");
    addParamsLine(" [--fsc+ <fscFile=\"\">]     : FSC file produced by xmipp_resolution_fsc");
    addExampleLine("xmipp_correct_bfactor -i volume.vol -o correctedVolume.vol --auto --sampling 1.4 --maxres 10");
    addExampleLine("To plot the Guinier file you may use:",false);
    addExampleLine("gnuplot");
    addExampleLine("set xlabel \"d^(-2)\"");
    addExampleLine("plot \"correctedVolume.vol.guinier\" using 1:2 title \"log F\" with line, \"correctedVolume.vol.guinier\" using 1:4 title \"Corrected log F\"with line");
}

void ProgVolumeCorrectBfactor::readParams()
{
    produces_an_output = true;

    XmippMetadataProgram::readParams();
    if (checkParam("--ref"))
    {
        fn_ref= getParam("--ref");
        String strMode=getParam("--ref",1);
        mode =  strMode=="allpoints" ? ALLPOINTS_REF : BFACTOR_REF;
    }
    else if (checkParam("--adhoc"))
    {
        mode = BFACTOR_ADHOC;
        adhocB = getDoubleParam("--adhoc");
    }
    else if (checkParam("--auto"))
        mode = BFACTOR_AUTO;
    else
        REPORT_ERROR(ERR_DEBUG_IMPOSIBLE, "This should not happens, review program definition");
    sampling_rate = getDoubleParam("--sampling");
    apply_maxres = getDoubleParam("--maxres");
    fit_minres = getDoubleParam("--fit_minres");
    fit_maxres = getDoubleParam("--fit_maxres");
    fn_fsc = getParam("--fsc");

    if (fit_maxres < 0.)
        fit_maxres = apply_maxres;
}

void ProgVolumeCorrectBfactor::show()
{
    XmippMetadataProgram::show();
    std::cout << "Pixel size : " << sampling_rate << " Angstrom" << std::endl;
    std::cout << "Maximum resolution: " << apply_maxres << " Angstrom" << std::endl;
    if (mode == BFACTOR_REF || mode == BFACTOR_AUTO)
    {
        std::cerr<<"Fit within resolutions: " << fit_minres << " - " << fit_maxres << " Angstrom" << std::endl;
    }
    if (mode == BFACTOR_REF)
    {
        std::cout << "Adjust B-factor according to reference "<<fn_ref<<std::endl;
    }
    else if (mode == BFACTOR_ADHOC)
    {
        std::cout << "Apply ad-hoc B-factor of "<< adhocB <<" squared Angstroms" << std::endl;
    }
    else
    {
        std::cout << "Use automated B-factor fit (Rosenthal and Henderson, 2003) " << std::endl;
    }
    if (fn_fsc != "")
        std::cout << "Use signal-to-noise weight based on "<< fn_fsc <<std::endl;
}

void ProgVolumeCorrectBfactor::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    Image<double> vol;
    vol.read(fnImg);
    vol().checkDimensionWithDebug(3,__FILE__,__LINE__);
    FileName fn_guinier = fn_out + ".guinier";
    bfactor_correction(vol(), fn_guinier);
    vol.write(fn_out);
}

ProgVolumeCorrectBfactor::ProgVolumeCorrectBfactor()
{
    fit_minres    = -1.;
    fit_maxres    = -1.;
    apply_maxres  = -1.;
    sampling_rate = 1.;
    mode          = BFACTOR_AUTO;
    xsize         = -1;
    fn_ref        = "";
    fn_fsc        = "";
    adhocB        = 0.;
}

void  ProgVolumeCorrectBfactor::make_guinier_plot(MultidimArray< std::complex< double > > &FT1,
        std::vector<fit_point2D> &guinier)
{
    MultidimArray< int >  radial_count(xsize);
    MultidimArray<double> lnF(xsize);
    Matrix1D<double>      f(3);
    fit_point2D      onepoint;

    lnF.initZeros();
    for (size_t k=0; k<ZSIZE(FT1); k++)
    {
        FFT_IDX2DIGFREQ(k,ZSIZE(FT1),ZZ(f));
        double z2=ZZ(f)*ZZ(f);
        for (size_t i=0; i<YSIZE(FT1); i++)
        {
            FFT_IDX2DIGFREQ(i,YSIZE(FT1),YY(f));
            double y2z2=z2+YY(f)*YY(f);
            for (size_t j=0; j<XSIZE(FT1); j++)
            {
                FFT_IDX2DIGFREQ(j,xsize,XX(f));
                double R2=y2z2+XX(f)*XX(f);
                if (R2>0.25)
                    continue;
                double R=sqrt(R2);
                int idx=ROUND(R*xsize);
                A1D_ELEM(lnF,idx) += abs(dAkij(FT1, k, i, j));
                ++A1D_ELEM(radial_count,idx);
            }
        }
    }

    guinier.clear();
    for (size_t i = 0; i < XSIZE(radial_count); i++)
    {
        double res = (xsize * sampling_rate)/(double)i;
        if (res >= apply_maxres)
        {
            onepoint.x = 1. / (res * res);
            if (lnF(i)>0.)
            {
                onepoint.y = log ( lnF(i) / radial_count(i) );
                if (res <= fit_minres && res >= fit_maxres)
                {
                    onepoint.w = 1.;
                }
                else
                {
                    onepoint.w = 0.;
                }
            }
            else
            {
                onepoint.y = 0.;
                onepoint.w = 0.;
            }
            guinier.push_back(onepoint);
        }
    }
}

void ProgVolumeCorrectBfactor::get_snr_weights(std::vector<double> &snr)
{
    std::ifstream  fh;
    std::string    line;
    float ires, fsc, noise, res;
    int line_no = 0;

    snr.clear();

    fh.open((fn_fsc).c_str(), std::ios::in);
    if (!fh)
        REPORT_ERROR(ERR_IO_NOTEXIST, fn_fsc);

    // Count the number of lines
    fh.peek();
    getline(fh, line);
    while (!fh.eof())
    {
        getline(fh, line);
        sscanf(line.c_str(), "%f %f %f %f", &ires, &fsc, &noise, &res);
        double mysnr = XMIPP_MAX((double)(2*fsc) / (1+fsc), 0.);
        snr.push_back( sqrt(mysnr) );
        line_no++;
        fh.peek();
    }

    fh.close();

    if (line_no != xsize/2)
    {
        std::cerr<<"line_no = "<<line_no <<" neq xsize/2= "<<xsize/2<<std::endl;
        REPORT_ERROR(ERR_VALUE_INCORRECT,"ERROR: invalid FSC file");
    }
}

void  ProgVolumeCorrectBfactor::apply_snr_weights(MultidimArray< std::complex< double > > &FT1,
        std::vector<double> &snr)
{

    Matrix1D<double> f(3);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(FT1)
    {
        FFT_IDX2DIGFREQ(j,xsize,XX(f));
        FFT_IDX2DIGFREQ(i,YSIZE(FT1),YY(f));
        FFT_IDX2DIGFREQ(k,ZSIZE(FT1),ZZ(f));
        double R = f.module();
        if (sampling_rate / R >= apply_maxres)
        {
            int idx=ROUND(R*xsize);
            dAkij(FT1, k, i, j) *= snr[idx];
        }
    }
}

void  ProgVolumeCorrectBfactor::apply_bfactor(MultidimArray< std::complex< double > > &FT1,
                                        double bfactor)
{
    Matrix1D<double> f(3);
    double isampling_rate=1.0/sampling_rate;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(FT1)
    {
        FFT_IDX2DIGFREQ(j,xsize,XX(f));
        FFT_IDX2DIGFREQ(i,YSIZE(FT1),YY(f));
        FFT_IDX2DIGFREQ(k,ZSIZE(FT1),ZZ(f));
        double R = f.module() * isampling_rate;
        if (1./R >= apply_maxres)
        {
            dAkij(FT1, k, i, j) *= exp( -0.25* bfactor  * R * R);
        }
    }
}

void  ProgVolumeCorrectBfactor::apply_allpoints(MultidimArray< std::complex< double > > &FT1,
        std::vector<fit_point2D> &guinier_diff)
{
    Matrix1D<double> f(3);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(FT1)
    {
        FFT_IDX2DIGFREQ(j,xsize,XX(f));
        FFT_IDX2DIGFREQ(i,YSIZE(FT1),YY(f));
        FFT_IDX2DIGFREQ(k,ZSIZE(FT1),ZZ(f));
        double R=f.module();
        if (R>0.5)
            continue;
        size_t idx=lround(R*xsize);
        if (idx < guinier_diff.size() && guinier_diff[idx].w > 0.)
        {
            dAkij(FT1, k, i, j) *= exp( -guinier_diff[idx].y );
        }
    }
}

void  ProgVolumeCorrectBfactor::write_guinierfile(const FileName &fn_guinier,
        std::vector<fit_point2D> &guinierin,
        std::vector<fit_point2D> &guinierweighted,
        std::vector<fit_point2D> &guiniernew,
        double intercept,
        std::vector<fit_point2D> &guinierref)
{
    std::ofstream fh;
    fh.open((fn_guinier).c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR(ERR_IO_NOWRITE, fn_guinier);

    fh << "# 1/d^2   lnF     weighted-lnF   corrected-lnF   (model)\n";
    for (size_t i = 0; i < guinierin.size(); i++)
    {
        fh << (guinierin[i]).x << " " << (guinierin[i]).y << " " << (guinierweighted[i]).y << " " <<(guiniernew[i]).y;
        if (mode==BFACTOR_AUTO)
            fh << " " << intercept;
        else if (mode==BFACTOR_REF || mode==ALLPOINTS_REF)
        {
            fh << " " << (guinierref[i]).y + intercept;
        }
        fh << std::endl;
    }
    fh.close();
}

void ProgVolumeCorrectBfactor::bfactor_correction(MultidimArray< double > &m1,
        const FileName &fn_guinier)
{
    MultidimArray< std::complex< double > > FT1, FT2;
    FourierTransformer transformer;
    double slope, intercept;
    std::vector<fit_point2D>  guinierin, guinierweighted, guinierref, guinierdiff;
    std::vector<double> snr;

    // Transform
    xsize=XSIZE(m1);
    transformer.FourierTransform(m1, FT1, true);
    make_guinier_plot(FT1,guinierin);

    // Get SNR weights and apply to guinier1
    if (fn_fsc != "")
    {
        get_snr_weights(snr);
        apply_snr_weights(FT1,snr);
        make_guinier_plot(FT1,guinierweighted);
    }
    else
        guinierweighted=guinierin;

    if (mode == BFACTOR_AUTO)
    {
        least_squares_line_fit(guinierweighted, slope, intercept);
        std::cerr<<" Fitted slope= "<<slope<<" intercept= "<<intercept<<std::endl;
        adhocB = 4. * slope;
    }
    else if (mode == BFACTOR_REF || mode == ALLPOINTS_REF)
    {
        Image<double> ref;
        ref.read(fn_ref);
        ref().setXmippOrigin();
        transformer.FourierTransform(ref(), FT2, true);
        make_guinier_plot(FT2,guinierref);
        guinierdiff = guinierweighted;
        for (size_t i = 0; i < guinierdiff.size(); i++)
        {
            (guinierdiff[i]).y -= (guinierref[i]).y;
        }
        if (mode == BFACTOR_REF)
        {
            least_squares_line_fit(guinierdiff, slope, intercept);
            std::cerr<<" Fitted slope= "<<slope<<" intercept= "<<intercept<<std::endl;
            adhocB = 4. * slope;
        }
    }

    if (mode == ALLPOINTS_REF)
    {
        // Now apply the allpoints correction
        std::cerr<<"Adjust power spectrum to that of reference "<<std::endl;
        apply_allpoints(FT1,guinierdiff);
    }
    else
    {
        // Now apply the B-factor
        std::cerr<<"Applying B-factor of "<< adhocB << " squared Angstroms"<<std::endl;
        apply_bfactor(FT1,adhocB);
    }
    // Now backtransform
    transformer.inverseFourierTransform(FT1,m1);

    // Write Guinier-plot
    std::vector<fit_point2D>  guiniernew;
    make_guinier_plot(FT1,guiniernew);
    write_guinierfile(fn_guinier, guinierin, guinierweighted, guiniernew, intercept, guinierref);
}
