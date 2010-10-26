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
#include <data/transformations.h>

/* Read parameters --------------------------------------------------------- */
void Prog_Sort_PSD_Parameters::readParams()
{
    fnSel = getParam("-i");
    filter_w1 = getDoubleParam("-f1");
    filter_w2 = getDoubleParam("-f2");
    decay_width = getDoubleParam("-decay");
    mask_w1 = getDoubleParam("-m1");
    mask_w2 = getDoubleParam("-m2");
}

/* Usage ------------------------------------------------------------------- */
void Prog_Sort_PSD_Parameters::defineParams()
{
    addUsageLine("This program evaluates the CTFs and PSDs of a set of micrographs");
    addUsageLine("as produced by the preprocessing micrographs step of the Xmipp protocols");
    addParamsLine("   -i <selfile>              : Selfile with micrographs");
    addParamsLine("==+ Enhancement filter parameters");
    addParamsLine("  [-f1 <freq_low=0.02>]      : Low freq. for band pass filtration, max 0.5");
    addParamsLine("  [-f2 <freq_high=0.2>]      : High freq. for band pass filtration, max 0.5");
    addParamsLine("  [-decay <freq_decay=0.02>] : Decay for the transition bands");
    addParamsLine("  [-m1 <mfreq_low=0.01>]     : Low freq. for mask, max 0.5");
    addParamsLine("  [-m2 <mfreq_high=0.45>]    : High freq. for mask, max 0.5");
}

/* Show -------------------------------------------------------------------- */
void Prog_Sort_PSD_Parameters::show() const
{
    std::cout << "Selfile:      " << fnSel << std::endl
    << "Filter w1:    " << filter_w1 << std::endl
    << "Filter w2:    " << filter_w2 << std::endl
    << "Filter decay: " << decay_width << std::endl
    << "Mask w1:      " << mask_w1 << std::endl
    << "Mask w2:      " << mask_w2 << std::endl;
}

/* Compute Correlation ----------------------------------------------------- */
double Prog_Sort_PSD_Parameters::evaluate(const FileName &fnMicrograph,
    const FileName &fnPSD, const FileName &fnCTF, const FileName &fnCTF2,
    PSDEvaluation &evaluation) const
{
    FileName fnRoot = fnMicrograph.withoutExtension();

    // Read input data
    Image<double> PSD;
    PSD.read(fnPSD);
    CTFDescription CTF1, CTF2;
    CTF1.read(fnCTF);
    CTF1.Produce_Side_Info();
    if (fnCTF2!="")
    {
    	CTF2.read(fnCTF2);
    	CTF2.Produce_Side_Info();
    }

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
    PSD.write(fnRoot + "_Periodogramavg_enhanced.xmp");

    // Evaluate the radial integral
    PSD().setXmippOrigin();
    Matrix1D< int > center_of_rot(2);
    MultidimArray< double > radial_mean;
    MultidimArray<int> radial_count;
    radialAverage(PSD(),center_of_rot,radial_mean,radial_count);
    radial_mean.selfABS();
    evaluation.PSDradialIntegral=radial_mean.sum();

    // Rotate 90 degrees and compute correlation
    Image<double> PSDrotated(PSD);
    selfRotate(BSPLINE3, PSDrotated(), 90);
    evaluation.PSDcorrelation90=correlation_index(PSD(), PSDrotated());

    // Get the fitting score
    MetaData MD;
    MD.read(fnCTF);
    MD.getValue(MDL_CTF_CRITERION_FITTINGSCORE,evaluation.fittingScore);

    // Explore the CTF
    Matrix1D<double> u(2), freqZero1(2), freqZero2(2);
    double wmax=0.5/CTF1.Tm;
    double maxModule=0, minModule=1e38, min;
    double N=0;
    evaluation.maxDampingAtBorder=0;
    evaluation.firstZeroDisagreement=0;
    evaluation.firstZeroAvg=0;
    for (double alpha=0; alpha<=PI; alpha+=PI/180, N++)
    {
    	VECTOR_R2(u,cos(alpha),sin(alpha));

    	// Get the zero in the direction of u
    	CTF1.zero(1, u, freqZero1);
    	double module=1.0/freqZero1.module();
    	maxModule=XMIPP_MAX(maxModule,module);
    	minModule=XMIPP_MIN(minModule,module);
    	evaluation.firstZeroAvg+=module;
    	double wx=wmax*XX(u);
    	double wy=wmax*YY(u);
    	double damping=CTF1.CTFdamping_at(wx,wy);
    	damping=damping*damping;
    	evaluation.maxDampingAtBorder=XMIPP_MAX(evaluation.maxDampingAtBorder,damping);
    	if (fnCTF2!="") {
        	CTF2.zero(1, u, freqZero2);
        	double module2=1.0/freqZero2.module();
        	double diff=ABS(module-module2);
        	evaluation.firstZeroDisagreement=XMIPP_MAX(evaluation.firstZeroDisagreement,diff);
    	}
    }
    evaluation.firstZeroAvg/=N;
    evaluation.firstZeroRatio=maxModule/minModule;
}

/* Run --------------------------------------------------------------------- */
void Prog_Sort_PSD_Parameters::run()
{
    MetaData SF(fnSel);
    PSDEvaluation evaluation;
    init_progress_bar(SF.size());
    int idx=0;
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        FileName fnMicrograph, fnPSD, fnCTF, fnCTF2;
        SF.getValue(MDL_IMAGE,fnMicrograph);
        SF.getValue(MDL_PSD,fnPSD);
        SF.getValue(MDL_CTFMODEL,fnCTF);
        if (SF.containsLabel(MDL_CTFMODEL2))
        	SF.getValue(MDL_CTFMODEL2,fnCTF2);
        evaluate(fnMicrograph, fnPSD, fnCTF, fnCTF2, evaluation);
        SF.setValue(MDL_CTF_CRITERION_DAMPING,evaluation.maxDampingAtBorder);
        SF.setValue(MDL_CTF_CRITERION_FIRSTZEROAVG,evaluation.firstZeroAvg);
        SF.setValue(MDL_CTF_CRITERION_FIRSTZERODISAGREEMENT,evaluation.firstZeroDisagreement);
        SF.setValue(MDL_CTF_CRITERION_FIRSTZERORATIO,evaluation.firstZeroRatio);
        SF.setValue(MDL_CTF_CRITERION_FITTINGSCORE,evaluation.fittingScore);
        SF.setValue(MDL_CTF_CRITERION_PSDCORRELATION90,evaluation.PSDcorrelation90);
        SF.setValue(MDL_CTF_CRITERION_PSDRADIALINTEGRAL,evaluation.PSDradialIntegral);
        progress_bar(++idx);
    }
    SF.write(fnSel);
}
