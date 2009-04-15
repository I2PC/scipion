/***************************************************************************
 *
 * Author:     Sjors H.W. Scheres               (scheres@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include "correct_bfactor.h"

Prog_correct_bfactor_prm::Prog_correct_bfactor_prm()
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
    return;

}

void  Prog_correct_bfactor_prm::make_guinier_plot(Matrix3D< std::complex< double > > &FT1,
                                                  std::vector<fit_point2D> &guinier)
{

    Matrix1D< int >  radial_count(xsize);
    Matrix1D<double> lnF(xsize),f(3);
    fit_point2D      onepoint;

    lnF.initZeros();
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(FT1)
    {
        FFT_IDX2DIGFREQ(j,xsize,XX(f));
        FFT_IDX2DIGFREQ(i,YSIZE(FT1),YY(f));
        FFT_IDX2DIGFREQ(k,ZSIZE(FT1),ZZ(f));
        double R=f.module();
        if (R>0.5) continue;
        int idx=ROUND(R*xsize);
        lnF(idx) += abs(dVkij(FT1, k, i, j));
        radial_count(idx)++;
    }

    guinier.clear();
    for (int i = 0; i < XSIZE(radial_count); i++)
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

void Prog_correct_bfactor_prm::get_snr_weights(std::vector<double> &snr)
{
    std::ifstream  fh;
    std::string    line;
    float ires, fsc, noise, res;
    int line_no = 0;

    snr.clear();

    fh.open((fn_fsc).c_str(), std::ios::in);
    if (!fh)
        REPORT_ERROR(3008, (std::string)"Prog_correct_bfactor_prm: Cannot read file: " + fn_fsc);

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
        REPORT_ERROR(1,"ERROR: invalid FSC file");
    }

}
void  Prog_correct_bfactor_prm::apply_snr_weights(Matrix3D< std::complex< double > > &FT1,
                                                  std::vector<double> &snr)
{

    Matrix1D<double> f(3);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(FT1)
    {
        FFT_IDX2DIGFREQ(j,xsize,XX(f));
        FFT_IDX2DIGFREQ(i,YSIZE(FT1),YY(f));
        FFT_IDX2DIGFREQ(k,ZSIZE(FT1),ZZ(f));
        double R = f.module();
        if (sampling_rate / R >= apply_maxres)
        {
            int idx=ROUND(R*xsize);
            dVkij(FT1, k, i, j) *= snr[idx];
        }
    }
}


void  Prog_correct_bfactor_prm::apply_bfactor(Matrix3D< std::complex< double > > &FT1,
                                              double bfactor)
{

    Matrix1D<double> f(3);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(FT1)
    {
        FFT_IDX2DIGFREQ(j,xsize,XX(f));
        FFT_IDX2DIGFREQ(i,YSIZE(FT1),YY(f));
        FFT_IDX2DIGFREQ(k,ZSIZE(FT1),ZZ(f));
        double R = f.module() / sampling_rate;
        if (1./R >= apply_maxres)
        {
            dVkij(FT1, k, i, j) *= exp( -(bfactor / 4)  * R * R);
        }
    }
}

void  Prog_correct_bfactor_prm::apply_allpoints(Matrix3D< std::complex< double > > &FT1,
                                                std::vector<fit_point2D> &guinier_diff)
{

    Matrix1D<double> f(3);
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(FT1)
    {
        FFT_IDX2DIGFREQ(j,xsize,XX(f));
        FFT_IDX2DIGFREQ(i,YSIZE(FT1),YY(f));
        FFT_IDX2DIGFREQ(k,ZSIZE(FT1),ZZ(f));
        double R=f.module();
        if (R>0.5) continue;
        int idx=ROUND(R*xsize);
        if (idx < guinier_diff.size() && guinier_diff[idx].w > 0.) 
        {
            dVkij(FT1, k, i, j) *= exp( -guinier_diff[idx].y );
        }
    }
}



void  Prog_correct_bfactor_prm::write_guinierfile(FileName fn_guinier, 
                                                  std::vector<fit_point2D> &guinierin,
                                                  std::vector<fit_point2D> &guinierweighted,
                                                  std::vector<fit_point2D> &guiniernew,
                                                  double intercept,
                                                  std::vector<fit_point2D> &guinierref)
{
    std::ofstream fh;
    fh.open((fn_guinier).c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR(3008, (std::string)"Prog_correct_bfactor_prm: Cannot write file: " + fn_guinier);

    fh << "# 1/d^2   lnF     weighted-lnF   corrected-lnF   (model)\n";
    for (int i = 0; i < guinierin.size(); i++)
    {
        fh << (guinierin[i]).x << " " << (guinierin[i]).y << " " << (guinierweighted[i]).y << " " <<(guiniernew[i]).y;
        if (mode==BFACTOR_AUTO)
            fh << " " << intercept;
        else if (mode==BFACTOR_REF || mode==ALLPOINTS_REF)
        {
            fh << " " << (guinierref[i]).y + intercept;
        }
        fh << "\n"<<std::endl;
    }
    fh.close();
    
}

void Prog_correct_bfactor_prm::bfactor_correction(Matrix3D< double > &m1, FileName fn_guinier)
{
    
    Matrix3D< std::complex< double > > FT1, FT2;
    XmippFftw transformer;
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
        VolumeXmipp ref;
        ref.read(fn_ref);
        ref().setXmippOrigin();
        transformer.FourierTransform(ref(), FT2, true);
        make_guinier_plot(FT2,guinierref);
        guinierdiff = guinierweighted;
        for (int i = 0; i < guinierdiff.size(); i++)
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
