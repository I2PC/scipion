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

#include "reconstruct_random.h"

// Define params
void ProgRecRandom::defineParams()
{
    //usage
    addUsageLine("Generate 3D reconstructions from projections using random orientations");
    //params
    addParamsLine("   -i <md_file>                : Metadata file with input projections");
    addParamsLine("  [--oroot <volume_file=\"rec_random\">]  : Filename for output rootname");
    addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections");
    addParamsLine("  [--randomIter <N=7>]         : Number of iterations with random assignment");
    addParamsLine("  [--greedyIter <N=3>]         : Number of iterations with greedy assignment");
    addParamsLine("  [--rejection <p=25>]         : Percentage of images to reject for reconstruction");
    addParamsLine("  [--initial <file=\"\">]      : Initial volume if available");
}

// Read arguments ==========================================================
void ProgRecRandom::readParams()
{
    fnIn = getParam("-i");
    fnRoot = getParam("--oroot");
    fnSym = getParam("--sym");
    NiterRandom = getIntParam("--randomIter");
    NiterGreedy = getIntParam("--greedyIter");
    rejection = getDoubleParam("--rejection");
    fnInit = getParam("--initial");
}

// Show ====================================================================
void ProgRecRandom::show()
{
    if (verbose > 0)
    {
        std::cout << "Input metadata              : "  << fnIn        << std::endl;
        std::cout << "Output rootname             : "  << fnRoot      << std::endl;
        std::cout << "Number of random iterations : "  << NiterRandom << std::endl;
        std::cout << "Number of greedy iterations : "  << NiterGreedy << std::endl;
        std::cout << "Rejection percentage        : "  << rejection   << std::endl;
        if (fnSym != "")
            std::cout << "Symmetry for projections    : "  << fnSym << std::endl;
        if (fnInit !="")
            std::cout << "Initial volume              : "  << fnInit << std::endl;
    }
}

// Main routine ------------------------------------------------------------
void ProgRecRandom::run()
{
    show();
    produceSideinfo();

    bool finish=false;
    iter=1;
    do
    {
    	// Generate projections from the volume
    	generateProjections();

    	// Align the input images to the projections
    	size_t nImg=0;
    	double avgCorr=0, avgImprovement=0;
    	mdReconstruction.clear();
    	init_progress_bar(mdIn.size());
    	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
    	{
    		double corr, improvement;
    		alignSingleImage(nImg++,__iter.objId, corr, improvement);
    		avgCorr+=corr;
    		avgImprovement+=improvement;
    		progress_bar(nImg);
    	}
    	progress_bar(mdIn.size());
    	avgCorr/=mdIn.size();
    	avgImprovement/=mdIn.size();
    	std::cout << "Iter " << iter << " avg.correlation=" << avgCorr << " avg.improvement=" << avgImprovement << std::endl;

    	// Remove too good and too bad images
    	filterByCorrelation();
    	//mdReconstruction.write(fnAngles);

    	// Reconstruct
    	reconstructCurrent();

    	finish=(iter==(NiterRandom+NiterGreedy));
    	iter++;
    } while (!finish);
}

//#define DEBUG
void ProgRecRandom::alignSingleImage(size_t nImg, size_t id, double &newCorr, double &improvementFraction)
{
	mCurrentImage.aliasImageInStack(inputImages(),nImg);
	mCurrentImage.setXmippOrigin();

	MultidimArray<double> allCorrs;
	allCorrs.resizeNoCopy(mdGallery.size());

	size_t nGallery=0;
	size_t improvementCount=0;
	double oldCorr;
	mdIn.getValue(MDL_MAXCC,oldCorr,id);
#ifdef DEBUG
	FileName fnImg;
	mdIn.getValue(MDL_IMAGE,fnImg,id);
	std::cout << "Image: " << fnImg << " oldCorr=" << oldCorr << std::endl;
#endif
	std::vector< Matrix2D<double> > allM;
	FOR_ALL_OBJECTS_IN_METADATA(mdGallery)
	{
		mCurrentImageAligned=mCurrentImage;
		mGalleryProjection.aliasImageInStack(gallery(),nGallery);
		mGalleryProjection.setXmippOrigin();
		double corr=alignImagesConsideringMirrors(mGalleryProjection,mCurrentImageAligned,M,DONT_WRAP);
#ifdef DEBUG
		mdGallery.getValue(MDL_MAXCC,corr,__iter.objId);
#endif
		A1D_ELEM(allCorrs,nGallery)=corr;
		allM.push_back(M);
		++nGallery;
		if (corr>oldCorr)
		{
			++improvementCount;
#ifdef DEBUG
		mdGallery.getValue(MDL_IMAGE,fnImg,__iter.objId);
		std::cout << "   Matching Gallery: " << fnImg << " " << corr << std::endl;
#endif
		}
	}
	improvementFraction=(double)improvementCount/NSIZE(gallery());
#ifdef DEBUG
		mdGallery.write("PPPgallery.xmd");
		std::cout << "   oldcorr=" << oldCorr << std::endl;
#endif

	MultidimArray<double> scaledCorrs;
	allCorrs.cumlativeDensityFunction(scaledCorrs);
	nGallery=0;
	double bestCorr, bestAngleRot, bestAngleTilt, bestAnglePsi, bestShiftX, bestShiftY, bestWeight;
	bool bestFlip;
	bestCorr=-1;
	double correctionFactor=1;
	/*
	if (improvementCount>1)
		correctionFactor=1.0/improvementCount;
	else if (iter<=NiterRandom)
		correctionFactor=1.0/(0.98*XSIZE(scaledCorrs));
	else
		correctionFactor=1;*/
	FOR_ALL_OBJECTS_IN_METADATA(mdGallery)
	{
		bool getThis=(A1D_ELEM(allCorrs,nGallery)>oldCorr) || (improvementCount==0 && A1D_ELEM(scaledCorrs,nGallery)>=0.98);
		if (getThis)
		{
			bool flip;
			double scale, psi, shiftX, shiftY, angleRot, angleTilt, anglePsi;
			transformationMatrix2Parameters2D(allM[nGallery],flip,scale,shiftX,shiftY,anglePsi);
			anglePsi*=-1;
			double weight=A1D_ELEM(scaledCorrs,nGallery)*correctionFactor;
			double corr=A1D_ELEM(allCorrs,nGallery);
			mdGallery.getValue(MDL_ANGLE_ROT,angleRot,__iter.objId);
			mdGallery.getValue(MDL_ANGLE_TILT,angleTilt,__iter.objId);
#ifdef DEBUG
			mdGallery.getValue(MDL_IMAGE,fnImg,__iter.objId);
			std::cout << "   Getting Gallery: " << fnImg << " corr=" << corr << " cdf=" << weight << " rot=" << angleRot
					  << " tilt=" << angleTilt << std::endl;
#endif

			if (iter<=NiterRandom)
			{
				size_t recId=mdReconstruction.addObject();
				FileName fnImg;
				mdIn.getValue(MDL_IMAGE,fnImg,id);
				mdReconstruction.setValue(MDL_IMAGE,fnImg,recId);
				mdReconstruction.setValue(MDL_ENABLED,1,recId);
				mdReconstruction.setValue(MDL_MAXCC,corr,recId);
				mdReconstruction.setValue(MDL_ANGLE_ROT,angleRot,recId);
				mdReconstruction.setValue(MDL_ANGLE_TILT,angleTilt,recId);
				mdReconstruction.setValue(MDL_ANGLE_PSI,anglePsi,recId);
				mdReconstruction.setValue(MDL_SHIFT_X,shiftX,recId);
				mdReconstruction.setValue(MDL_SHIFT_Y,shiftY,recId);
				mdReconstruction.setValue(MDL_FLIP,flip,recId);
				mdReconstruction.setValue(MDL_WEIGHT,weight,recId);
			}
			if (corr>bestCorr)
			{
				bestCorr=corr;
				bestAngleRot=angleRot;
				bestAngleTilt=angleTilt;
				bestAnglePsi=anglePsi;
				bestShiftX=shiftX;
				bestShiftY=shiftY;
				bestFlip=flip;
				bestWeight=weight;
			}
		}

		++nGallery;
	}

	if (iter>NiterRandom)
	{
		size_t recId=mdReconstruction.addObject();
		FileName fnImg;
		mdIn.getValue(MDL_IMAGE,fnImg,id);
		mdReconstruction.setValue(MDL_IMAGE,fnImg,recId);
		mdReconstruction.setValue(MDL_ENABLED,1,recId);
		mdReconstruction.setValue(MDL_MAXCC,bestCorr,recId);
		mdReconstruction.setValue(MDL_ANGLE_ROT,bestAngleRot,recId);
		mdReconstruction.setValue(MDL_ANGLE_TILT,bestAngleTilt,recId);
		mdReconstruction.setValue(MDL_ANGLE_PSI,bestAnglePsi,recId);
		mdReconstruction.setValue(MDL_SHIFT_X,bestShiftX,recId);
		mdReconstruction.setValue(MDL_SHIFT_Y,bestShiftY,recId);
		mdReconstruction.setValue(MDL_FLIP,bestFlip,recId);
		mdReconstruction.setValue(MDL_WEIGHT,bestWeight,recId);
	}
	mdIn.setValue(MDL_MAXCC,bestCorr,id);
	newCorr=bestCorr;
#ifdef DEBUG
	mdReconstruction.write("PPPreconstruction.xmd");
	std::cout << "Press any key" << std::endl;
	char c;
	std::cin >> c;
#endif
}
#undef DEBUG

void ProgRecRandom::filterByCorrelation()
{
	MetaData mdAux;
	mdAux=mdReconstruction;
	mdAux.removeDisabled();

	std::vector<double> correlations;
	mdAux.getColumnValues(MDL_MAXCC,correlations);
	std::sort(correlations.begin(),correlations.end());
	size_t skip=(size_t)floor(correlations.size()*(rejection/100.0));
	double minCorr=correlations[skip];
	double maxCorr=correlations[correlations.size()-skip/4-1];

	std::vector<double> weights;
	mdAux.getColumnValues(MDL_WEIGHT,weights);
	std::sort(weights.begin(),weights.end());
	double minWeight=weights[skip];
	double maxWeight=weights[weights.size()-skip/4-1];

	FOR_ALL_OBJECTS_IN_METADATA(mdAux)
	{
		double cc, weight;
		mdAux.getValue(MDL_MAXCC,cc,__iter.objId);
		if (cc<minCorr) // COSS || cc>maxCorr)
			mdAux.setValue(MDL_ENABLED,-1,__iter.objId);
		mdAux.getValue(MDL_WEIGHT,weight,__iter.objId);
		if (weight<minWeight) // COSS || weight>maxWeight)
			mdAux.setValue(MDL_ENABLED,-1,__iter.objId);
	}
	mdAux.removeDisabled();
	mdAux.write(fnAngles);
}

void ProgRecRandom::reconstructCurrent()
{
	String args=formatString("-i %s -o %s --sym %s --weight -v 0",fnAngles.c_str(),fnVolume.c_str(),fnSym.c_str());
	String cmd=(String)"xmipp_reconstruct_fourier "+args;
	system(cmd.c_str());

	args=formatString("-i %s --mask circular %d -v 0",fnVolume.c_str(),-XSIZE(inputImages())/2);
	cmd=(String)"xmipp_transform_mask "+args;
	system(cmd.c_str());

	if (iter<=NiterRandom)
	{
		args=formatString("-i %s --select below 0 --substitute value 0 -v 0",fnVolume.c_str());
		cmd=(String)"xmipp_transform_threshold "+args;
		system(cmd.c_str());
	}
}

void ProgRecRandom::generateProjections()
{
	String args=formatString("-i %s -o %s --sampling_rate 5 --sym %s --compute_neighbors --angular_distance -1 --experimental_images %s -v 0",
			fnVolume.c_str(),fnGallery.c_str(),fnSym.c_str(),fnAngles.c_str());
	String cmd=(String)"xmipp_angular_project_library "+args;
	system(cmd.c_str());
	mdGallery.read(fnGalleryMetaData);
	gallery.read(fnGallery);
}

void ProgRecRandom::produceSideinfo()
{
	mdIn.read(fnIn);
	mdIn.removeDisabled();

	mdIn.fillConstant(MDL_MAXCC,"0.0");
	mdIn.fillRandom(MDL_ANGLE_ROT,"uniform",0,360);
	mdIn.fillRandom(MDL_ANGLE_TILT,"uniform",0,180);
	mdIn.fillRandom(MDL_ANGLE_PSI,"uniform",0,360);
	mdIn.fillConstant(MDL_SHIFT_X,"0.0");
	mdIn.fillConstant(MDL_SHIFT_Y,"0.0");
	mdIn.fillConstant(MDL_WEIGHT,"1.0");

	fnAngles=fnRoot+".xmd";
	fnVolume=fnRoot+".vol";
	fnGallery=fnRoot+"_gallery.stk";
	fnGalleryMetaData=fnRoot+"_gallery.doc";

	size_t xdim, ydim, zdim, ndim;
	getImageSize(mdIn,xdim, ydim, zdim, ndim);
	inputImages().resizeNoCopy(mdIn.size(), 1, ydim, xdim);

	FileName fnImg;
	Image<double> I;
	size_t n=0;
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		mdIn.getValue(MDL_IMAGE,fnImg,__iter.objId);
		I.read(fnImg);
		mCurrentImage.aliasImageInStack(inputImages(),n++);
		memcpy(MULTIDIM_ARRAY(mCurrentImage),MULTIDIM_ARRAY(I()),MULTIDIM_SIZE(mCurrentImage)*sizeof(double));
	}

	mdIn.write(fnAngles);
	if (fnInit=="")
		reconstructCurrent();
	else
		system(formatString("cp %s %s",fnInit.c_str(),fnVolume.c_str()).c_str());
}
