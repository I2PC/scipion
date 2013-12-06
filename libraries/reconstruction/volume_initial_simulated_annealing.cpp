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

#include "volume_initial_simulated_annealing.h"

// Define params
void ProgVolumeInitialSimulatedAnnealing::defineParams()
{
    //usage
    addUsageLine("Generate 3D reconstructions from projections using random orientations");
    //params
    addParamsLine("   -i <md_file>                : Metadata file with input projections");
    addParamsLine("  [--oroot <volume_file=\"rec_random\">]  : Filename for output rootname");
    addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections");
    addParamsLine("  [--randomIter <N=10>]        : Number of iterations with random assignment");
    addParamsLine("  [--greedyIter <N=0>]         : Number of iterations with greedy assignment");
    addParamsLine("  [--rejection <p=25>]         : Percentage of images to reject for reconstruction");
    addParamsLine("  [--initial <file=\"\">]      : Initial volume if available");
    addParamsLine("  [--keepIntermediateVolumes]  : Keep the volume of each iteration");
    addParamsLine("  [--dontApplyPositive]        : Do not apply positive constraint in the random iterations");
    addParamsLine("  [--T0 <T=0.1>]               : Initial temperature for simulated annealing");
    addParamsLine("  [--thr <n=1>]                : Number of threads");
}

// Read arguments ==========================================================
void ProgVolumeInitialSimulatedAnnealing::readParams()
{
    fnIn = getParam("-i");
    fnRoot = getParam("--oroot");
    fnSym = getParam("--sym");
    T0 = getDoubleParam("--T0");
    NiterRandom = getIntParam("--randomIter");
    NiterGreedy = getIntParam("--greedyIter");
    rejection = getDoubleParam("--rejection");
    fnInit = getParam("--initial");
    Nthr = getIntParam("--thr");
    positiveConstraint = !checkParam("--dontApplyPositive");
    keepIntermediateVolumes = checkParam("--keepIntermediateVolumes");
}

// Show ====================================================================
void ProgVolumeInitialSimulatedAnnealing::show()
{
    if (verbose > 0)
    {
        std::cout << "Input metadata              : "  << fnIn        << std::endl;
        std::cout << "Output rootname             : "  << fnRoot      << std::endl;
        std::cout << "T0                          : "  << T0          << std::endl;
        std::cout << "Number of random iterations : "  << NiterRandom << std::endl;
        std::cout << "Number of greedy iterations : "  << NiterGreedy << std::endl;
        std::cout << "Rejection percentage        : "  << rejection   << std::endl;
        std::cout << "Number of threads           : "  << Nthr        << std::endl;
        std::cout << "Apply positive constraint   : "  << positiveConstraint << std::endl;
        std::cout << "Keep intermediate volumes   : "  << keepIntermediateVolumes << std::endl;
        if (fnSym != "")
            std::cout << "Symmetry for projections    : "  << fnSym << std::endl;
        if (fnInit !="")
            std::cout << "Initial volume              : "  << fnInit << std::endl;
    }
}

// Alignment of a single image ============================================
//#define DEBUG
void alignSingleImage(size_t nImg, ProgVolumeInitialSimulatedAnnealing &prm, MetaData &mdReconstruction, double &newCorr, double &improvementFraction)
{
	MultidimArray<double> mGalleryProjection, mCurrentImage, mCurrentImageAligned;
	Matrix2D<double> M;

	mCurrentImage.aliasImageInStack(prm.inputImages(),nImg);
	mCurrentImage.setXmippOrigin();

	MultidimArray<double> allCorrs;
	allCorrs.resizeNoCopy(prm.mdGallery.size());

	size_t improvementCount=0;
	double oldCorr=prm.mdInp[nImg].maxcc;
#ifdef DEBUG
	FileName fnImg;
	prm.mdIn.getValue(MDL_IMAGE,fnImg,id);
	std::cout << "Image: " << fnImg << " oldCorr=" << oldCorr << std::endl;
#endif
	std::vector< Matrix2D<double> > allM;
	for (size_t nGallery=0; nGallery<XSIZE(allCorrs); ++nGallery)
	{
		mCurrentImageAligned=mCurrentImage;
		mGalleryProjection.aliasImageInStack(prm.gallery(),nGallery);
		mGalleryProjection.setXmippOrigin();
		double corr=alignImagesConsideringMirrors(mGalleryProjection,mCurrentImageAligned,M,DONT_WRAP);
#ifdef DEBUG
		mdGallery.getValue(MDL_MAXCC,corr,__iter.objId);
#endif
		A1D_ELEM(allCorrs,nGallery)=corr;
		allM.push_back(M);
		if (corr>oldCorr)
		{
			++improvementCount;
#ifdef DEBUG
		prm.mdGallery.getValue(MDL_IMAGE,fnImg,__iter.objId);
		std::cout << "   Matching Gallery: " << fnImg << " " << corr << std::endl;
#endif
		}
	}
	improvementFraction=(double)improvementCount/NSIZE(prm.gallery());
#ifdef DEBUG
		mdGallery.write("PPPgallery.xmd");
		std::cout << "   oldcorr=" << oldCorr << std::endl;
#endif

	MultidimArray<double> scaledCorrs;
	allCorrs.cumlativeDensityFunction(scaledCorrs);
	double bestCorr, bestAngleRot, bestAngleTilt, bestAnglePsi, bestShiftX, bestShiftY, bestWeight;
	bool bestFlip;
	bestCorr=-1;
	double correctionFactor=1;
	double currentT=prm.T0*(1.0-((double)prm.iter)/prm.NiterRandom);
	double icurrentT=1.0/currentT;
	for (size_t nGallery=0; nGallery<XSIZE(allCorrs); ++nGallery)
	{
		bool getThis=(A1D_ELEM(allCorrs,nGallery)>oldCorr) || (improvementCount==0 && A1D_ELEM(scaledCorrs,nGallery)>=0.98);
		if (A1D_ELEM(allCorrs,nGallery)<oldCorr && !getThis && prm.iter<=prm.NiterRandom)
		{
			// Simulated annealing
			double diff=oldCorr-A1D_ELEM(allCorrs,nGallery);
			double p=rnd_unif();
			getThis=(p<exp(-diff*icurrentT));
		}
		if (getThis)
		{
			bool flip;
			double scale, shiftX, shiftY, anglePsi;
			transformationMatrix2Parameters2D(allM[nGallery],flip,scale,shiftX,shiftY,anglePsi);
			anglePsi*=-1;
			double weight=A1D_ELEM(scaledCorrs,nGallery)*correctionFactor;
			double corr=A1D_ELEM(allCorrs,nGallery);
			double angleRot=prm.mdGallery[nGallery].rot;
			double angleTilt=prm.mdGallery[nGallery].tilt;
#ifdef DEBUG
			fnImg=prm.mdGallery[nGallery].fnImg;
			std::cout << "   Getting Gallery: " << fnImg << " corr=" << corr << " cdf=" << weight << " rot=" << angleRot
					  << " tilt=" << angleTilt << std::endl;
#endif

			if (prm.iter<=prm.NiterRandom)
			{
				size_t recId=mdReconstruction.addObject();
				FileName fnImg;
				fnImg=prm.mdInp[nImg].fnImg;
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

	if (prm.iter>prm.NiterRandom)
	{
		size_t recId=mdReconstruction.addObject();
		FileName fnImg;
		fnImg=prm.mdInp[nImg].fnImg;
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
	newCorr=bestCorr;
#ifdef DEBUG
	mdReconstruction.write("PPPreconstruction.xmd");
	std::cout << "Press any key" << std::endl;
	char c;
	std::cin >> c;
#endif
}
#undef DEBUG

// Subset alignment =======================================================
void threadAlignSubset(ThreadArgument &thArg)
{
	ProgVolumeInitialSimulatedAnnealing &prm=*((ProgVolumeInitialSimulatedAnnealing *)thArg.workClass);
	ThreadVolumeInitialAlignment results=prm.threadResults[thArg.thread_id];

	results.sumCorr=results.sumImprovement=0.0;
	results.mdReconstruction.clear();
	int nMax=(int)prm.mdInp.size();
	for (int nImg=0; nImg<nMax; ++nImg)
	{
		if ((nImg+1)%prm.Nthr==thArg.thread_id)
		{
			double corr, improvement;
			alignSingleImage(nImg, prm, results.mdReconstruction, corr, improvement);
			results.sumCorr+=corr;
			results.sumImprovement+=improvement;
			prm.mdInp[nImg].maxcc=corr;
		}

		if (thArg.thread_id==0)
			progress_bar(nImg+1);
	}

	// Update all MAXCC
	prm.mutexMaxCC.lock();

	prm.mdReconstruction.unionAll(results.mdReconstruction);
	prm.sumCorr+=results.sumCorr;
	prm.sumImprovement+=results.sumImprovement;

	prm.mutexMaxCC.unlock();
}

// Main routine ------------------------------------------------------------
void ProgVolumeInitialSimulatedAnnealing::run()
{
    show();
    produceSideinfo();

    bool finish=false;
    iter=1;
    ThreadManager thMgr(Nthr,this);
    do
    {
    	// Generate projections from the volume
    	generateProjections();

    	// Align the input images to the projections
    	mdReconstruction.clear();
    	sumCorr=sumImprovement=0;
    	init_progress_bar(mdIn.size());
    	thMgr.run(threadAlignSubset);
    	progress_bar(mdIn.size());
    	size_t nImg=0;
    	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
    		mdIn.setValue(MDL_MAXCC,mdInp[nImg++].maxcc,__iter.objId);
    	std::cout << "Iter " << iter << " avg.correlation=" << sumCorr/mdIn.size()
    			  << " avg.improvement=" << sumImprovement/mdIn.size() << std::endl;

    	// Remove too good and too bad images
    	filterByCorrelation();

    	// Reconstruct
    	reconstructCurrent();

    	finish=(iter==(NiterRandom+NiterGreedy));
    	iter++;
    } while (!finish);
    deleteFile(fnRoot+"_gallery_sampling.xmd");
    deleteFile(fnRoot+"_gallery.stk");
    deleteFile(fnRoot+"_gallery.doc");
}

void ProgVolumeInitialSimulatedAnnealing::filterByCorrelation()
{
	MetaData mdAux;
	mdAux=mdReconstruction;
	mdAux.removeDisabled();

	std::vector<double> correlations;
	mdAux.getColumnValues(MDL_MAXCC,correlations);
	std::sort(correlations.begin(),correlations.end());
	size_t skip=(size_t)floor(correlations.size()*(rejection/100.0));
	double minCorr=correlations[skip];
	//double maxCorr=correlations[correlations.size()-skip/4-1];

	std::vector<double> weights;
	mdAux.getColumnValues(MDL_WEIGHT,weights);
	std::sort(weights.begin(),weights.end());
	double minWeight=weights[skip];
	//double maxWeight=weights[weights.size()-skip/4-1];

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

void ProgVolumeInitialSimulatedAnnealing::reconstructCurrent()
{
	String args=formatString("-i %s -o %s --sym %s --weight --thr %d -v 0",fnAngles.c_str(),fnVolume.c_str(),fnSym.c_str(),Nthr);
	String cmd=(String)"xmipp_reconstruct_fourier "+args;
	if (system(cmd.c_str())==-1)
		REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");

	args=formatString("-i %s --mask circular %d -v 0",fnVolume.c_str(),-XSIZE(inputImages())/2);
	cmd=(String)"xmipp_transform_mask "+args;
	if (system(cmd.c_str())==-1)
		REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");

	if (iter<=NiterRandom && positiveConstraint)
	{
		args=formatString("-i %s --select below 0 --substitute value 0 -v 0",fnVolume.c_str());
		cmd=(String)"xmipp_transform_threshold "+args;
		if (system(cmd.c_str())==-1)
			REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
	}

	if (keepIntermediateVolumes)
	{
		cmd=formatString("cp %s %s_iter%02d.vol",fnVolume.c_str(),fnRoot.c_str(),iter);
		if (system(cmd.c_str())==-1)
			REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
	}
}

void ProgVolumeInitialSimulatedAnnealing::generateProjections()
{
	String args=formatString("-i %s -o %s --sampling_rate 5 --sym %s --compute_neighbors --angular_distance -1 --experimental_images %s -v 0",
			fnVolume.c_str(),fnGallery.c_str(),fnSym.c_str(),fnAngles.c_str());
	String cmd=(String)"xmipp_angular_project_library "+args;
	if (system(cmd.c_str())==-1)
		REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
	MetaData mdAux(fnGalleryMetaData);
	mdGallery.clear();
	FOR_ALL_OBJECTS_IN_METADATA(mdAux)
	{
		GalleryImage I;
		mdAux.getValue(MDL_IMAGE,I.fnImg,__iter.objId);
		mdAux.getValue(MDL_ANGLE_ROT,I.rot,__iter.objId);
		mdAux.getValue(MDL_ANGLE_TILT,I.tilt,__iter.objId);
		mdGallery.push_back(I);
	}
	gallery.read(fnGallery);
}

void ProgVolumeInitialSimulatedAnnealing::produceSideinfo()
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

	Image<double> I;
	size_t n=0;
	MultidimArray<double> mCurrentImage;
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		InputImage Ip;
		mdIn.getValue(MDL_IMAGE,Ip.fnImg,__iter.objId);
		Ip.maxcc=0.0;
		mdInp.push_back(Ip);
		I.read(Ip.fnImg);
		mCurrentImage.aliasImageInStack(inputImages(),n++);
		memcpy(MULTIDIM_ARRAY(mCurrentImage),MULTIDIM_ARRAY(I()),MULTIDIM_SIZE(mCurrentImage)*sizeof(double));
	}

	mdIn.write(fnAngles);
	iter=0;
	if (fnInit=="")
		reconstructCurrent();
	else
	{
		if (system(formatString("cp %s %s",fnInit.c_str(),fnVolume.c_str()).c_str())==-1)
			REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
	}

	threadResults=new ThreadVolumeInitialAlignment[Nthr];
}
