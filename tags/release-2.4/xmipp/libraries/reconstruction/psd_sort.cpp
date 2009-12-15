/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Slavica Jonic (slavica.jonic@impmc.jussieu.fr)
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

#include <vector>
#include "psd_sort.h"
#include "psd_enhance.h"
#include "ctf_estimate_from_micrograph.h"
#include <data/args.h>
#include <data/filters.h>

/* Read parameters --------------------------------------------------------- */
void Prog_Sort_PSD_Parameters::read(int argc, char **argv)
{
    fnSel = getParameter(argc,argv,"-i");
    fnOut = getParameter(argc,argv,"-o");
    windowSize = textToInteger(getParameter(argc, argv, "-N", "512"));
    filter_w1 = textToFloat(getParameter(argc, argv, "-f1", "0.05"));
    filter_w2 = textToFloat(getParameter(argc, argv, "-f2", "0.2"));
    decay_width = textToFloat(getParameter(argc, argv, "-decay", "0.02"));
    mask_w1 = textToFloat(getParameter(argc, argv, "-m1", "0.025"));
    mask_w2 = textToFloat(getParameter(argc, argv, "-m2", "0.2"));
}

/* Usage ------------------------------------------------------------------- */
void Prog_Sort_PSD_Parameters::usage() const
{
    std::cerr << "psd_sort\n"
              << "   -i <selfile>             : Selfile with micrographs\n"
              << "   -o <outputfile>          : List of micrographs with its correlations\n"
              << "  [-N <size=512>]           : Size of windows for the PSD estimation\n"
              << "  [-f1 <freq_low=0.05>]     : Low freq. for band pass filtration, max 0.5\n"
              << "  [-f2 <freq_high=0.2>]     : High freq. for band pass filtration, max 0.5\n"
              << "  [-decay <freq_decay=0.02>]: Decay for the transition bands\n"
              << "  [-m1 <freq_low=0.025>]    : Low freq. for mask, max 0.5\n"
              << "  [-m2 <freq_high=0.2>      : High freq. for mask, max 0.5\n"
    ;
}

/* Show -------------------------------------------------------------------- */
void Prog_Sort_PSD_Parameters::show() const
{
    std::cout << "Selfile:      " << fnSel       << std::endl
              << "Output:       " << fnOut       << std::endl
              << "N:            " << windowSize  << std::endl
              << "Filter w1:    " << filter_w1   << std::endl
              << "Filter w2:    " << filter_w2   << std::endl
              << "Filter decay: " << decay_width << std::endl
              << "Mask w1:      " << mask_w1     << std::endl
              << "Mask w2:      " << mask_w2     << std::endl
    ;
}

/* Produce side information ------------------------------------------------ */
void Prog_Sort_PSD_Parameters::produceSideInfo()
{
    SF.read(fnSel);
}

/* Compute Correlation ----------------------------------------------------- */
double Prog_Sort_PSD_Parameters::computeCorrelation(
    const FileName &fnMicrograph) const
{
    std::cerr << "Processing " << fnMicrograph << std::endl;
    FileName fn_root  = fnMicrograph.without_extension();

    // Compute PSD
    Prog_assign_CTF_prm assignCTF;
    assignCTF.image_fn = fnMicrograph;
    assignCTF.selfile_mode = false;
    assignCTF.reversed=false;
    assignCTF.N_horizontal=windowSize;
    assignCTF.N_vertical=assignCTF.N_horizontal;
    assignCTF.compute_at_particle=false;
    assignCTF.micrograph_averaging=true;
    assignCTF.piece_averaging=false;
    assignCTF.Nside_piece=5;
    assignCTF.PSD_mode = Prog_assign_CTF_prm::Periodogram;
    assignCTF.dont_adjust_CTF=true;
    assignCTF.selfile_fn = "";
    assignCTF.picked_fn = "";
    assignCTF.PSDfn_root = fn_root + "_Periodogram";
    assignCTF.process();
    
    // Load the PSD
    ImageXmipp PSD;
    PSD.read(fn_root+"_Periodogramavg.psd");

    // Enhance the PSD
    Prog_Enhance_PSD_Parameters enhancePSD;
    enhancePSD.center = true;
    enhancePSD.take_log = true;
    enhancePSD.filter_w1 = filter_w1;
    enhancePSD.filter_w2 = filter_w2;
    enhancePSD.decay_width = decay_width;
    enhancePSD.mask_w1 = mask_w1;
    enhancePSD.mask_w2 = mask_w2;
    enhancePSD.apply(PSD());
    PSD.write(fn_root+"_Periodogramavg_enhanced.psd");

    // Rotate 90º and compute correlation
    ImageXmipp PSDrotated(PSD);
    PSDrotated().selfRotateBSpline(3,90);
    return correlation_index(PSD(),PSDrotated());
}

/* Run --------------------------------------------------------------------- */
void Prog_Sort_PSD_Parameters::run()
{
    // Compute the correlation of all micrographs
    SF.go_first_ACTIVE();
    correlation.resize(SF.ImgNo());
    std::vector<FileName> filenames;
    FOR_ALL_ELEMENTS_IN_MATRIX1D(correlation)
    {
        FileName fnMicrograph=SF.NextImg();
        filenames.push_back(fnMicrograph);
        correlation(i)=computeCorrelation(fnMicrograph);
    }
    
    // Sort the correlations
    Matrix1D<int> idx=correlation.indexSort();
    
    // Produce output
    std::ofstream fhOut;
    fhOut.open(fnOut.c_str());
    if (!fhOut)
        REPORT_ERROR(1,(std::string)"Cannot open "+fnOut+" for output");
    FOR_ALL_ELEMENTS_IN_MATRIX1D(idx)
    {
        int ii=idx(i)-1;
        fhOut << filenames[ii] << " \t" << correlation(ii) << std::endl;
    }
    fhOut.close();
}
