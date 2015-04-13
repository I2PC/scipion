/***************************************************************************
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
 *
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

#include "classify_evaluate_classes.h"
#include <data/metadata_extension.h>

ClassEvaluation::ClassEvaluation()
{
	FRC_05=DPR_05=0;
	overfitting=false;
}

void evaluateClass(MetaData &MD, ClassEvaluation &eval)
{
    eval.FRC_05=0;
    eval.DPR_05=0;
    if (MD.size()<10)
    	return;

    MetaData MDrandomized;
	std::vector<MetaData> vMD;
    MultidimArray<double> I0, I1, freq, frc, dpr, frc_noise, error_l2;

    // Compute FRC
    MDrandomized.randomize(MD);
    MDrandomized.split(2,vMD,MDL_IMAGE);
    getAverageApplyGeo(vMD[0],I0);
    getAverageApplyGeo(vMD[1],I1);
    I0.setXmippOrigin();
    I1.setXmippOrigin();
    frc_dpr(I0, I1, 1, freq, frc, frc_noise, dpr, error_l2, true);

    // Compute the frequency of FRC=0.5
    int i_05=-1;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(frc)
    if (dAi(frc,i)<0.5)
    {
    	i_05=i;
    	break;
    }
    if (i_05==-1)
    {
    	i_05=XSIZE(frc)-1;
    	eval.overfitting=true;
    }

    // Extract evaluators
    eval.FRC_05=dAi(freq,i_05);
    eval.DPR_05=dAi(dpr,i_05);
}

// Evaluate classes program -------------------------------------------
void ProgEvaluateClass::readParams()
{
	fnClass=getParam("-i");
	fnOut=getParam("-o");
}

void ProgEvaluateClass::defineParams()
{
    addUsageLine("Evaluate the quality of a set of classes");
    addUsageLine("+The program associates to each class a number of quality measures:");
    addUsageLine("+");
    addUsageLine("+FRC_05 is the digital frequency (<0.5) at which the Fourier Ring Correlation drops below 0.5 ");
    addUsageLine("+(to convert from a digital frequency to one measured in Angstrom invert the digital frequency and multiply by the sampling rate).");
    addUsageLine("+");
    addUsageLine("+DPR_05 is the Differential Phase Residual at the frequency of FRC_05");
    addParamsLine(" -i <infile>             : Metadata with the classification (normally the output of CL2D or ML2D)");
    addParamsLine(" [-o <outfile=\"\">]     : Output file");
    addExampleLine("xmipp_classify_evaluate_classes -i 2D/CL2D/run_001/results_level_00_classes.xmd");
}

void ProgEvaluateClass::show()
{
	if (!verbose)
		return;
	std::cout
	<< "Input:  " << fnClass << std::endl
	<< "Output: " << fnOut << std::endl;
}

void ProgEvaluateClass::run()
{
	MetaData MD((String)"classes@"+fnClass), MDclass;
	ClassEvaluation eval;
	if (verbose>0)
		init_progress_bar(MD.size());
	int idx=0;
	FOR_ALL_OBJECTS_IN_METADATA(MD)
	{
		int classNo;
		MD.getValue(MDL_REF,classNo,__iter.objId);
		MDclass.read(formatString("class%06d_images@%s",classNo,fnClass.c_str()));
		evaluateClass(MDclass,eval);
		MD.setValue(MDL_CLASSIFICATION_FRC_05,eval.FRC_05,__iter.objId);
		MD.setValue(MDL_CLASSIFICATION_DPR_05,eval.DPR_05,__iter.objId);
		idx++;
		if (verbose>0)
			progress_bar(idx);
	}
	if (verbose>0)
		progress_bar(MD.size());
	if (fnOut=="")
		fnOut=fnClass;
	MD.write((String)"classes@"+fnOut,MD_APPEND);
}
