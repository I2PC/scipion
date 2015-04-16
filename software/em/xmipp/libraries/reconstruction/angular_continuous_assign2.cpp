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

#include "angular_continuous_assign2.h"
#include <data/mask.h>
#include <data/numerical_tools.h>

// Empty constructor =======================================================
ProgAngularContinuousAssign2::ProgAngularContinuousAssign2()
{
    produces_a_metadata = true;
    each_image_produces_an_output = true;
    projector = NULL;
}

ProgAngularContinuousAssign2::~ProgAngularContinuousAssign2()
{
	delete projector;
}

// Read arguments ==========================================================
void ProgAngularContinuousAssign2::readParams()
{
	XmippMetadataProgram::readParams();
    fnVol = getParam("--ref");
    maxShift = getDoubleParam("--max_shift");
    maxScale = getDoubleParam("--max_scale");
    maxAngularChange = getDoubleParam("--max_angular_change");
    maxResol = getDoubleParam("--max_resolution");
    Ts = getDoubleParam("--sampling");
    Rmax = getIntParam("--Rmax");
    pad = getIntParam("--padding");
    optimizeGrayValues = checkParam("--optimizeGray");
    optimizeShift = checkParam("--optimizeShift");
    optimizeScale = checkParam("--optimizeScale");
    optimizeAngles = checkParam("--optimizeAngles");
    originalImageLabel = getParam("--applyTo");
}

// Show ====================================================================
void ProgAngularContinuousAssign2::show()
{
    if (!verbose)
        return;
	XmippMetadataProgram::show();
    std::cout
    << "Reference volume:    " << fnVol              << std::endl
    << "Max. Shift:          " << maxShift           << std::endl
    << "Max. Scale:          " << maxScale           << std::endl
    << "Max. Angular Change: " << maxAngularChange   << std::endl
    << "Max. Resolution:     " << maxResol           << std::endl
    << "Sampling:            " << Ts                 << std::endl
    << "Max. Radius:         " << Rmax               << std::endl
    << "Padding factor:      " << pad                << std::endl
    << "Optimize gray:       " << optimizeGrayValues << std::endl
    << "Optimize shifts:     " << optimizeShift      << std::endl
    << "Optimize scale:      " << optimizeScale      << std::endl
    << "Optimize angles:     " << optimizeAngles     << std::endl
    << "Apply to:            " << originalImageLabel << std::endl
    ;
}

// usage ===================================================================
void ProgAngularContinuousAssign2::defineParams()
{
    addUsageLine("Make a continuous angular assignment");
	defaultComments["-i"].clear();
	defaultComments["-i"].addComment("Metadata with initial alignment");
	defaultComments["-o"].clear();
	defaultComments["-o"].addComment("Metadata with output alignment");
    XmippMetadataProgram::defineParams();
    addParamsLine("   --ref <volume>              : Reference volume");
    addParamsLine("  [--max_shift <s=-1>]         : Maximum shift allowed in pixels");
    addParamsLine("  [--max_scale <s=0.02>]       : Maximum scale change");
    addParamsLine("  [--max_angular_change <a=5>] : Maximum angular change allowed");
    addParamsLine("  [--max_resolution <f=4>]     : Maximum resolution (A)");
    addParamsLine("  [--sampling <Ts=1>]          : Sampling rate (A/pixel)");
    addParamsLine("  [--Rmax <R=-1>]              : Maximum radius (px). -1=Half of volume size");
    addParamsLine("  [--padding <p=2>]            : Padding factor");
    addParamsLine("  [--optimizeGray]             : Optimize gray values");
    addParamsLine("  [--optimizeShift]            : Optimize shift");
    addParamsLine("  [--optimizeScale]            : Optimize scale");
    addParamsLine("  [--optimizeAngles]           : Optimize angles");
    addParamsLine("  [--applyTo <label=image>]    : Which is the source of images to apply the final transformation");
    addExampleLine("A typical use is:",false);
    addExampleLine("xmipp_angular_continuous_assign2 -i anglesFromDiscreteAssignment.xmd --ref reference.vol -o assigned_angles.xmd");
}

// Produce side information ================================================
void ProgAngularContinuousAssign2::preProcess()
{
    // Read the reference volume
    V.read(fnVol);
    V().setXmippOrigin();
    Xdim=XSIZE(V());

    Ip().initZeros(Xdim,Xdim);
    E().initZeros(Xdim,Xdim);
    Ifilteredp().initZeros(Xdim,Xdim);

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

    // Transformation matrix
    A.initIdentity(3);
}

double tranformImage(FourierProjector &projector, double rot, double tilt, double psi,
		Projection &P, const MultidimArray<double> &I, const MultidimArray<double> &Ifiltered,
		MultidimArray<double> &Ip, MultidimArray<double> &Ifilteredp, MultidimArray<double> &E,
		const MultidimArray<int> &mask2D, double iMask2Dsum, double a, double b, Matrix2D<double> &A, int degree)
{
    projectVolume(projector, P, (int)XSIZE(I), (int)XSIZE(I),  rot, tilt, psi);

    double cost=0;
	applyGeometry(degree,Ip,I,A,IS_NOT_INV,DONT_WRAP,0.);
	applyGeometry(degree,Ifilteredp,Ifiltered,A,IS_NOT_INV,DONT_WRAP,0.);
	const MultidimArray<double> &mP=P();
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I)
	{
		if (DIRECT_MULTIDIM_ELEM(mask2D,n))
		{
			DIRECT_MULTIDIM_ELEM(Ip,n)=a*DIRECT_MULTIDIM_ELEM(Ip,n)+b;
			DIRECT_MULTIDIM_ELEM(Ifilteredp,n)=a*DIRECT_MULTIDIM_ELEM(Ifilteredp,n)+b;
			DIRECT_MULTIDIM_ELEM(E,n)=DIRECT_MULTIDIM_ELEM(mP,n)-DIRECT_MULTIDIM_ELEM(Ifilteredp,n);
			//DIRECT_MULTIDIM_ELEM(E,n)=DIRECT_MULTIDIM_ELEM(mP,n)-DIRECT_MULTIDIM_ELEM(Ip,n);
			cost+=fabs(DIRECT_MULTIDIM_ELEM(E,n));
		}
		else
		{
			DIRECT_MULTIDIM_ELEM(Ip,n)=0;
			DIRECT_MULTIDIM_ELEM(Ifilteredp,n)=0;
		}
	}
	cost*=iMask2Dsum;
	return cost;
}


double L1cost(double *x, void *_prm)
{
	double a=x[1];
	double b=x[2];
	double deltax=x[3];
	double deltay=x[4];
	double scalex=x[5];
	double scaley=x[6];
	double deltaRot=x[7];
	double deltaTilt=x[8];
	double deltaPsi=x[9];
	ProgAngularContinuousAssign2 *prm=(ProgAngularContinuousAssign2 *)_prm;
	if (deltax*deltax+deltay*deltay>prm->maxShift*prm->maxShift)
		return 1e38;
	if (fabs(scalex)>prm->maxScale || fabs(scaley)>prm->maxScale)
		return 1e38;
	if (fabs(deltaRot)>prm->maxAngularChange || fabs(deltaTilt)>prm->maxAngularChange || fabs(deltaPsi)>prm->maxAngularChange)
		return 1e38;
	MAT_ELEM(prm->A,0,0)=1+scalex;
	MAT_ELEM(prm->A,1,1)=1+scaley;
	MAT_ELEM(prm->A,0,2)=prm->old_shiftX+deltax;
	MAT_ELEM(prm->A,1,2)=prm->old_shiftY+deltay;
	return tranformImage(*prm->projector,prm->old_rot+deltaRot, prm->old_tilt+deltaTilt, prm->old_psi+deltaPsi,
			prm->P, prm->I(), prm->Ifiltered(), prm->Ip(), prm->Ifilteredp(), prm->E(),
			prm->mask2D, prm->iMask2Dsum, a, b, prm->A, LINEAR);
}

// Predict =================================================================
//#define DEBUG
void ProgAngularContinuousAssign2::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    rowOut=rowIn;

    // Read input image and initial parameters
//  ApplyGeoParams geoParams;
//	geoParams.only_apply_shifts=false;
//	geoParams.wrap=DONT_WRAP;

	rowIn.getValue(MDL_ANGLE_ROT,old_rot);
	rowIn.getValue(MDL_ANGLE_TILT,old_tilt);
	rowIn.getValue(MDL_ANGLE_PSI,old_psi);
	rowIn.getValue(MDL_SHIFT_X,old_shiftX);
	rowIn.getValue(MDL_SHIFT_Y,old_shiftY);
	double old_scaleX=0, old_scaleY=0, old_grayA=1, old_grayB=0, old_contShiftX=0, old_contShiftY=0;
	if (rowIn.containsLabel(MDL_CONTINUOUS_GRAY_A))
	{
		rowIn.getValue(MDL_CONTINUOUS_GRAY_A,old_grayA);
		rowIn.getValue(MDL_CONTINUOUS_GRAY_B,old_grayB);
		rowIn.getValue(MDL_CONTINUOUS_SCALE_X,old_scaleX);
		rowIn.getValue(MDL_CONTINUOUS_SCALE_Y,old_scaleY);
		rowIn.getValue(MDL_CONTINUOUS_X,old_contShiftX);
		rowIn.getValue(MDL_CONTINUOUS_Y,old_contShiftY);
	}

	I.read(fnImg);
	I().setXmippOrigin();

    Ifiltered()=I();
    filter.applyMaskSpace(Ifiltered());

    Matrix1D<double> p(9), steps(9);
    p(0)=old_grayA; // a in I'=a*I+b
    p(1)=old_grayB; // b in I'=a*I+b
    if (optimizeGrayValues && !rowIn.containsLabel(MDL_CONTINUOUS_GRAY_A))
	{
        double Iavg, Istd, Pavg, Pstd;
        I().computeAvgStdev(Iavg,Istd);
        P().computeAvgStdev(Pavg,Pstd);
		p(0)=Pstd/Istd;
		p(1)=Pavg-Iavg;
	}
    p(2)=old_contShiftX;
    p(3)=old_contShiftY;
    p(4)=old_scaleX;
    p(5)=old_scaleY;

    // Optimize
	double cost=-1;
	if (old_contShiftX*old_contShiftX+old_contShiftY*old_contShiftY>maxShift*maxShift ||
        fabs(old_scaleX)>maxScale || fabs(old_scaleY)>maxScale)
    	rowOut.setValue(MDL_ENABLED,-1);
	else
	{
		cost=1e38;
		int iter;
		steps.initZeros();
		if (optimizeGrayValues)
			steps(0)=steps(1)=1.;
		if (optimizeShift)
			steps(2)=steps(3)=1.;
		if (optimizeScale)
			steps(4)=steps(5)=1.;
		if (optimizeAngles)
			steps(6)=steps(7)=steps(8)=1.;
		powellOptimizer(p, 1, 9, &L1cost, this, 0.01, cost, iter, steps, false);
		if (cost>1e30)
			rowOut.setValue(MDL_ENABLED,-1);

		// Apply
		FileName fnOrig;
		rowIn.getValue(MDL::str2Label(originalImageLabel),fnOrig);
		I.read(fnImg);
		if (XSIZE(Ip())!=XSIZE(I()))
		{
			scaleToSize(BSPLINE3,Ip(),I(),XSIZE(Ip()),YSIZE(Ip()));
			I()=Ip();
		}
		A(0,2)=p(2)+old_shiftX;
		A(1,2)=p(3)+old_shiftY;
		A(0,0)=1+p(4);
		A(1,1)=1+p(5);

		bool flip;
		rowIn.getValue(MDL_FLIP,flip);
		if (flip)
		{
			// M*A*M
			A(0,1)*=-1;
			A(1,0)*=-1;
			A(0,2)*=-1;
		}
		applyGeometry(BSPLINE3,Ip(),I(),A,IS_NOT_INV,DONT_WRAP);

		Ip.write(fnImgOut);
	}
    rowOut.setValue(MDL_IMAGE_ORIGINAL, fnImg);
    rowOut.setValue(MDL_IMAGE, fnImgOut);
    rowOut.setValue(MDL_ANGLE_ROT,  old_rot+p(6));
    rowOut.setValue(MDL_ANGLE_TILT, old_tilt+p(7));
    rowOut.setValue(MDL_ANGLE_PSI,  old_psi+p(8));
    rowOut.setValue(MDL_SHIFT_X,    0.);
    rowOut.setValue(MDL_SHIFT_Y,    0.);
    rowOut.setValue(MDL_COST,       cost);
    rowOut.setValue(MDL_CONTINUOUS_GRAY_A,p(0));
    rowOut.setValue(MDL_CONTINUOUS_GRAY_B,p(1));
    rowOut.setValue(MDL_CONTINUOUS_SCALE_X,p(4));
    rowOut.setValue(MDL_CONTINUOUS_SCALE_Y,p(5));
    rowOut.setValue(MDL_CONTINUOUS_X,old_shiftX+p(2));
    rowOut.setValue(MDL_CONTINUOUS_Y,old_shiftY+p(3));

#ifdef DEBUG
    Image<double> save;
    save()=P();
    save.write("PPPprojection.xmp");
    save()=I();
    save.write("PPPexperimental.xmp");
    Ip.write("PPPexperimentalp.xmp");
    Ifiltered.write("PPPexperimentalFiltered.xmp");
    Ifilteredp.write("PPPexperimentalFilteredp.xmp");
    E.write("PPPresidual.xmp");
    std::cout << A << std::endl;
    std::cout << fnImgOut << " rewritten\n";
    std::cout << "Press any key" << std::endl;
    char c; std::cin >> c;
#endif
}
#undef DEBUG

void ProgAngularContinuousAssign2::postProcess()
{
	double minCost=1e38;
	MetaData &ptrMdOut=*getOutputMd();
	FOR_ALL_OBJECTS_IN_METADATA(ptrMdOut)
	{
		double cost;
		ptrMdOut.getValue(MDL_COST,cost,__iter.objId);
		if (cost<minCost)
			minCost=cost;
	}
	FOR_ALL_OBJECTS_IN_METADATA(ptrMdOut)
	{
		double cost;
		ptrMdOut.getValue(MDL_COST,cost,__iter.objId);
		ptrMdOut.setValue(MDL_WEIGHT_CONTINUOUS2,minCost/cost,__iter.objId);
	}
	ptrMdOut.write(fn_out.replaceExtension("xmd"));
}
