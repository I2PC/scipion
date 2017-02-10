/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano coss@cnb.csic.es
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

#include "transform_adjust_image_grey_levels.h"
#include <data/mask.h>
#include <data/numerical_tools.h>

// Empty constructor =======================================================
ProgTransformImageGreyLevels::ProgTransformImageGreyLevels()
{
    produces_a_metadata = true;
    each_image_produces_an_output = true;
    projector = NULL;
}

ProgTransformImageGreyLevels::~ProgTransformImageGreyLevels()
{
	delete projector;
}

// Read arguments ==========================================================
void ProgTransformImageGreyLevels::readParams()
{
	XmippMetadataProgram::readParams();
    fnVol = getParam("--ref");
    maxResol = getDoubleParam("--max_resolution");
    maxA = getDoubleParam("--max_gray_scale");
    maxB = getDoubleParam("--max_gray_shift");
    Ts = getDoubleParam("--sampling");
    Rmax = getIntParam("--Rmax");
    pad = getIntParam("--padding");
}

// Show ====================================================================
void ProgTransformImageGreyLevels::show()
{
    if (!verbose)
        return;
	XmippMetadataProgram::show();
    std::cout
    << "Reference volume:    " << fnVol              << std::endl
    << "Max. Resolution:     " << maxResol           << std::endl
	<< "Max. Gray Scale:     " << maxA               << std::endl
	<< "Max. Gray Shift:     " << maxB               << std::endl
    << "Sampling:            " << Ts                 << std::endl
    << "Max. Radius:         " << Rmax               << std::endl
    << "Padding factor:      " << pad                << std::endl
    ;
}

// usage ===================================================================
void ProgTransformImageGreyLevels::defineParams()
{
    addUsageLine("Make a continuous angular assignment");
	defaultComments["-i"].clear();
	defaultComments["-i"].addComment("Metadata with images and alignment");
	defaultComments["-o"].clear();
	defaultComments["-o"].addComment("Stack of images prepared for 3D reconstruction");
    XmippMetadataProgram::defineParams();
    addParamsLine("   --ref <volume>              : Reference volume");
    addParamsLine("  [--max_resolution <f=4>]     : Maximum resolution (A)");
    addParamsLine("  [--max_gray_scale <a=0.05>]  : Maximum gray scale change");
    addParamsLine("  [--max_gray_shift <b=0.05>]  : Maximum gray shift change as a factor of the image standard deviation");
    addParamsLine("  [--sampling <Ts=1>]          : Sampling rate (A/pixel)");
    addParamsLine("  [--Rmax <R=-1>]              : Maximum radius (px). -1=Half of volume size");
    addParamsLine("  [--padding <p=2>]            : Padding factor");
}

// Produce side information ================================================
void ProgTransformImageGreyLevels::preProcess()
{
    // Read the reference volume
    Image<double> V;
    V.read(fnVol);
    V().setXmippOrigin();
    Xdim=XSIZE(V());

    // Construct mask
    if (Rmax<0)
    	Rmax=Xdim/2;
    Mask mask;
    mask.type = BINARY_CIRCULAR_MASK;
    mask.mode = INNER_MASK;
    mask.R1 = Rmax;
    mask.generate_mask(Xdim,Xdim);
    mask2D=mask.get_binary_mask();
    iMask2Dsum=1.0/mask2D.sum();

    // Construct projector
    projector = new FourierProjector(V(),pad,Ts/maxResol,BSPLINE3);

    // Low pass filter
    filter.FilterBand=LOWPASS;
    filter.w1=Ts/maxResol;
    filter.raised_w=0.02;
}

double transformImageGrayCost(double a, double b, ProgTransformImageGreyLevels *prm)
{
	if (fabs(a-1)>prm->maxA)
		return 1e38;
	if (fabs(b)>prm->maxB*prm->Istddev)
		return 1e38;

    double cost=0;
	const MultidimArray<double> &mP=prm->P();
	const MultidimArray<int> &mMask2D=prm->mask2D;
	MultidimArray<double> &mIfiltered=prm->Ifiltered();
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mMask2D)
	{
//		if (DIRECT_MULTIDIM_ELEM(mMask2D,n))
		{
			double val=fabs((a*DIRECT_MULTIDIM_ELEM(mP,n)+b)-DIRECT_MULTIDIM_ELEM(mIfiltered,n));
			cost+=val*val;
		}
	}
	cost*=prm->iMask2Dsum;
	return cost;
}

double transformImageGrayCostAB(double *x, void *_prm)
{
	return transformImageGrayCost(x[1],x[2],(ProgTransformImageGreyLevels *)_prm);
}

double transformImageGrayCostBA(double *x, void *_prm)
{
	return transformImageGrayCost(x[2],x[1],(ProgTransformImageGreyLevels *)_prm);
}

// Predict =================================================================
//#define DEBUG
void ProgTransformImageGreyLevels::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    rowOut=rowIn;

    // Read input image and initial parameters
    double rot, tilt, psi;
	rowIn.getValue(MDL_ANGLE_ROT,rot);
	rowIn.getValue(MDL_ANGLE_TILT,tilt);
	rowIn.getValue(MDL_ANGLE_PSI,psi);

	double olda=1.0, oldb=0.0;
	if (rowIn.containsLabel(MDL_CONTINUOUS_GRAY_A)){
		rowIn.getValue(MDL_CONTINUOUS_GRAY_A,olda);
		rowIn.getValue(MDL_CONTINUOUS_GRAY_B,oldb);
	}

	I.read(fnImg);
	I().setXmippOrigin();
	Istddev=I().computeStddev();

    Ifiltered()=I();
    filter.applyMaskSpace(Ifiltered());

	projectVolume(*projector, P, (int)XSIZE(I()), (int)XSIZE(I()),  rot, tilt, psi);

#ifdef DEBUG
    Image<double> save;
    save()=P();
    std::cout << "P: "; P().printStats(); std::cout << std::endl;
    save.write("PPPprojection.xmp");
    save()=I();
    std::cout << "I: "; I().printStats(); std::cout << std::endl;
    save.write("PPPexperimental.xmp");
    save()=P()-I();
    std::cout << "P-I init: "; save().printStats(); std::cout << std::endl;
    save.write("PPPdiffInit.xmp");
#endif

    Matrix1D<double> pAB(2), steps(2), pBA(2), p(2);
    steps.initConstant(1);

    // Optimize
	double costAB=-1, costBA=-1;
	int iter;
	try
	{
		costAB=1e38;
	    pAB(0)=olda; // a in I'=a*I+b
	    pAB(1)=oldb; // b in I'=a*I+b
		powellOptimizer(pAB, 1, 2, &transformImageGrayCostAB, this, 0.01, costAB, iter, steps, verbose>=2);
		if (costAB>1e30)
		{
			rowOut.setValue(MDL_ENABLED,-1);
			pAB.initZeros();
			pAB(0)=1; // a in I'=a*I+b
			pAB(1)=0; // b in I'=a*I+b
		}

		costBA=1e38;
	    pBA(0)=oldb; // a in I'=a*I+b
	    pBA(1)=olda; // b in I'=a*I+b
		powellOptimizer(pBA, 1, 2, &transformImageGrayCostBA, this, 0.01, costBA, iter, steps, verbose>=2);
		if (costBA>1e30)
		{
			rowOut.setValue(MDL_ENABLED,-1);
			pBA.initZeros();
			pBA(0)=0; // a in I'=a*I+b
			pBA(1)=1; // b in I'=a*I+b
		}

		// Decide
		if (costAB<costBA)
			p=pAB;
		else
		{
			p(0)=pBA(1);
			p(1)=pBA(0);
		}

		// Apply
		MultidimArray<double> &mI=I();
		double ia=1.0/p(0);
		double b=p(1);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mI)
		{
			if (DIRECT_MULTIDIM_ELEM(mask2D,n) || true)
				DIRECT_MULTIDIM_ELEM(mI,n)=ia*(DIRECT_MULTIDIM_ELEM(mI,n)-b);
			else
				DIRECT_MULTIDIM_ELEM(mI,n)=0.0;
		}
		rowOut.setValue(MDL_IMAGE,fnImgOut);
		I.write(fnImgOut);
	}
	catch (XmippError XE)
	{
		std::cerr << XE << std::endl;
		std::cerr << "Warning: Cannot refine " << fnImg << std::endl;
		rowOut.setValue(MDL_ENABLED,-1);
	}
    rowOut.setValue(MDL_CONTINUOUS_GRAY_A,p(0));
    rowOut.setValue(MDL_CONTINUOUS_GRAY_B,p(1));

#ifdef DEBUG
    save()=I();
    save.write("PPPexperimentalCorrected.xmp");
    std::cout << "I corrected: "; I().printStats(); std::cout << std::endl;
    save()=P()-I();
    std::cout << "P-I final: "; save().printStats(); std::cout << std::endl;
    save.write("PPPdiff.xmp");
    std::cout << fnImgOut << " rewritten\n";
    std::cout << "Press any key" << std::endl;
    char c; std::cin >> c;
#endif
}
#undef DEBUG
