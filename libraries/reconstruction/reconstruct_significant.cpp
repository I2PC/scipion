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

#include "reconstruct_significant.h"
#include <algorithm>

// Define params
void ProgReconstructSignificant::defineParams()
{
    //usage
    addUsageLine("Generate 3D reconstructions from projections using random orientations");
    //params
    addParamsLine("   -i <md_file>                : Metadata file with input projections");
    addParamsLine("  --initvolumes <md_file>      : Set of initial volumes");
    addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
    addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections");
    addParamsLine("  [--iter <N=10>]              : Number of iterations");
    addParamsLine("  [--alpha0 <N=0.05>]          : Initial significance");
    addParamsLine("  [--alphaF <N=0.005>]         : Final significance");
    addParamsLine("  [--keepIntermediateVolumes]  : Keep the volume of each iteration");
    addParamsLine("  [--thr <n=1>]                : Number of threads");
    addParamsLine("  [--angularSampling <a=5>]    : Angular sampling in degrees for generating the projection gallery");
    addParamsLine("  [--maxShift <s=-1>]          : Maximum shift allowed (+-this amount)");
    addParamsLine("  [--minTilt <t=0>]            : Minimum tilt angle");
    addParamsLine("  [--maxTilt <t=90>]           : Maximum tilt angle");
    addParamsLine("  [--useImed]                  : Use Imed for weighting");
    addParamsLine("  [--angDistance <a=10>]       : Angular distance");
}

// Read arguments ==========================================================
void ProgReconstructSignificant::readParams()
{
    fnIn = getParam("-i");
    fnDir = getParam("--odir");
    fnSym = getParam("--sym");
    fnInit = getParam("--initvolumes");
    alpha0 = getDoubleParam("--alpha0");
    alphaF = getDoubleParam("--alphaF");
    Niter = getIntParam("--iter");
    Nthr = getIntParam("--thr");
    keepIntermediateVolumes = checkParam("--keepIntermediateVolumes");
    angularSampling=getDoubleParam("--angularSampling");
    maxShift=getDoubleParam("--maxShift");
    tilt0=getDoubleParam("--minTilt");
    tiltF=getDoubleParam("--maxTilt");
    useImed=checkParam("--useImed");
    angDistance=getDoubleParam("--angDistance");
}

// Show ====================================================================
void ProgReconstructSignificant::show()
{
    if (verbose > 0)
    {
        std::cout << "Input metadata              : "  << fnIn        << std::endl;
        std::cout << "Output directory            : "  << fnDir       << std::endl;
        std::cout << "Initial significance        : "  << alpha0      << std::endl;
        std::cout << "Final significance          : "  << alphaF      << std::endl;
        std::cout << "Number of iterations        : "  << Niter       << std::endl;
        std::cout << "Number of threads           : "  << Nthr        << std::endl;
        std::cout << "Keep intermediate volumes   : "  << keepIntermediateVolumes << std::endl;
        std::cout << "Angular sampling            : "  << angularSampling << std::endl;
        std::cout << "Maximum shift               : "  << maxShift << std::endl;
        std::cout << "Minimum tilt                : "  << tilt0 << std::endl;
        std::cout << "Maximum tilt                : "  << tiltF << std::endl;
        std::cout << "Use Imed                    : "  << useImed << std::endl;
        std::cout << "Angular distance            : "  << angDistance << std::endl;
        if (fnSym != "")
            std::cout << "Symmetry for projections    : "  << fnSym << std::endl;
        if (fnInit !="")
            std::cout << "Initial volume              : "  << fnInit << std::endl;
    }
}

// Image alignment ========================================================
//#define DEBUG
void progReconstructSignificantThreadAlign(ThreadArgument &thArg)
{
    AlignmentAux aux;
    CorrelationAux aux2;
    RotationalCorrelationAux aux3;

    ProgReconstructSignificant &prm=*((ProgReconstructSignificant *)thArg.workClass);
	MultidimArray<double> mGalleryProjection, mCurrentImage, mCurrentImageAligned;
	Matrix2D<double> M;
	std::vector< Matrix2D<double> > allM;

	size_t Nimgs=ZSIZE(prm.cc);
	size_t Nvols=YSIZE(prm.cc);
	size_t Ndirs=XSIZE(prm.cc);
	MultidimArray<double> imgcc(Nvols*Ndirs), imgimed(Nvols*Ndirs);
	MultidimArray<double> cdfcc, cdfimed;
	double one_alpha=1-prm.currentAlpha;

	for (size_t nVolume=0; nVolume<Nvols; ++nVolume)
	{
		prm.mdReconstructionPartial[nVolume*prm.Nthr+thArg.thread_id].clear();
		prm.mdProjectionMatching[nVolume*prm.Nthr+thArg.thread_id].clear();
	}

	FileName fnImg;
	for (size_t nImg=0; nImg<Nimgs; ++nImg)
	{
		if (((int)nImg+1)%prm.Nthr==thArg.thread_id)
		{
#ifdef DEBUG
			std::cout << "Processing: " << prm.mdInp[nImg] << std::endl;
#endif
			mCurrentImage.aliasImageInStack(prm.inputImages(),nImg);
			mCurrentImage.setXmippOrigin();
			fnImg=prm.mdInp[nImg];
			allM.clear();

			double bestCorr=-2, bestRot, bestTilt, bestImed=1e38, worstImed=-1e38;
			Matrix2D<double> bestM;
			int bestVolume=-1;

			// Compute all correlations
	    	for (size_t nVolume=0; nVolume<Nvols; ++nVolume)
		    	for (size_t nDir=0; nDir<Ndirs; ++nDir)
				{
					mCurrentImageAligned=mCurrentImage;
					mGalleryProjection.aliasImageInStack(prm.gallery[nVolume](),nDir);
					mGalleryProjection.setXmippOrigin();
					double corr=alignImagesConsideringMirrors(mGalleryProjection,mCurrentImageAligned,M,aux,aux2,aux3,DONT_WRAP);
					// double corr=alignImagesConsideringMirrors(mGalleryProjection,mCurrentImageAligned,M,DONT_WRAP);
					M=M.inv();
					double imed=imedDistance(mGalleryProjection, mCurrentImageAligned);

//					if (corr>0.99)
//					{
//					std::cout << prm.mdGallery[nVolume][nDir].fnImg << " corr= " << corr << " imed= " << imed << std::endl;
//					std::cout << "Matrix=" << M << std::endl;
//					Image<double> save;
//					save()=mGalleryProjection;
//					save.write("PPPgallery.xmp");
//					save()=mCurrentImage;
//					save.write("PPPcurrentImage.xmp");
//					save()=mCurrentImageAligned;
//					save.write("PPPcurrentImageAligned.xmp");
//					char c; std::cin >> c;
//					}

					DIRECT_A3D_ELEM(prm.cc,nImg,nVolume,nDir)=corr;
					size_t idx=nVolume*Ndirs+nDir;
					DIRECT_A1D_ELEM(imgcc,idx)=corr;
					DIRECT_A1D_ELEM(imgimed,idx)=imed;
					allM.push_back(M);

					if (corr>bestCorr)
					{
						bestM=M;
						bestCorr=corr;
						bestVolume=(int)nVolume;
						bestRot=prm.mdGallery[nVolume][nDir].rot;
						bestTilt=prm.mdGallery[nVolume][nDir].tilt;
					}

					if (imed<bestImed)
						bestImed=imed;
					else if (imed>worstImed)
						worstImed=imed;
				}

	    	// Keep the best assignment for the projection matching
			MetaData &mdProjectionMatching=prm.mdProjectionMatching[bestVolume*prm.Nthr+thArg.thread_id];
			double scale, shiftX, shiftY, anglePsi;
			bool flip;
			transformationMatrix2Parameters2D(bestM,flip,scale,shiftX,shiftY,anglePsi);

			if (prm.maxShift<0 || (prm.maxShift>0 && fabs(shiftX)<prm.maxShift && fabs(shiftY)<prm.maxShift))
			{
				size_t recId=mdProjectionMatching.addObject();
				mdProjectionMatching.setValue(MDL_IMAGE,fnImg,recId);
				mdProjectionMatching.setValue(MDL_ENABLED,1,recId);
				mdProjectionMatching.setValue(MDL_MAXCC,bestCorr,recId);
				mdProjectionMatching.setValue(MDL_ANGLE_ROT,bestRot,recId);
				mdProjectionMatching.setValue(MDL_ANGLE_TILT,bestTilt,recId);
				mdProjectionMatching.setValue(MDL_ANGLE_PSI,anglePsi,recId);
				mdProjectionMatching.setValue(MDL_SHIFT_X,-shiftX,recId);
				mdProjectionMatching.setValue(MDL_SHIFT_Y,-shiftY,recId);
				mdProjectionMatching.setValue(MDL_FLIP,flip,recId);
			}

			// Compute lower limit of correlation
			double rl=bestCorr*one_alpha;
			double z=0.5*log((1+rl)/(1-rl));
			double zl=z-2.96*sqrt(1.0/MULTIDIM_SIZE(mCurrentImage));
			double ccl=tanh(zl);

	    	// Get the best images
			imgcc.cumlativeDensityFunction(cdfcc);
			imgimed.cumlativeDensityFunction(cdfimed);
			for (size_t nVolume=0; nVolume<Nvols; ++nVolume)
			{
				MetaData &mdPartial=prm.mdReconstructionPartial[nVolume*prm.Nthr+thArg.thread_id];
				for (size_t nDir=0; nDir<Ndirs; ++nDir)
				{
					size_t idx=nVolume*Ndirs+nDir;
					double cdfccthis=DIRECT_A1D_ELEM(cdfcc,idx);
					double cdfimedthis=DIRECT_A1D_ELEM(cdfimed,idx);
					double cc=DIRECT_A1D_ELEM(imgcc,idx);
//					if (cc>ccl)
//						std::cout << "    " << prm.mdGallery[nVolume][nDir].fnImg << " ***cc=" << cc << " cdfcc=" << cdfccthis << " cdfimedthis=" << cdfimedthis << std::endl;
//					if (cdfccthis>=one_alpha)
//						std::cout << "    " << prm.mdGallery[nVolume][nDir].fnImg << " cc=" << cc << " ***cdfcc=" << cdfccthis << " cdfimedthis=" << cdfimedthis << std::endl;
//					if (cdfimedthis<=prm.currentAlpha)
//						std::cout << "    " << prm.mdGallery[nVolume][nDir].fnImg << " cc=" << cc << " cdfcc=" << cdfccthis << " ***cdfimedthis=" << cdfimedthis << std::endl;
					bool condition=!prm.useImed || (prm.useImed && cdfimedthis<=prm.currentAlpha);
					if (cdfccthis>=one_alpha && cc>ccl && condition)
					{
						double imed=DIRECT_A1D_ELEM(imgimed,idx);
						transformationMatrix2Parameters2D(allM[nVolume*Ndirs+nDir],flip,scale,shiftX,shiftY,anglePsi);
						if (prm.maxShift>0)
							if (fabs(shiftX)>prm.maxShift || fabs(shiftY)>prm.maxShift)
								continue;
						if (flip)
							shiftX*=-1;

						double weight=cdfccthis*(cc/bestCorr);
						if (prm.useImed)
							weight*=(1-cdfimedthis)*(bestImed/imed);
						DIRECT_A3D_ELEM(prm.weight,nImg,nVolume,nDir)=weight;
						double angleRot=prm.mdGallery[nVolume][nDir].rot;
						double angleTilt=prm.mdGallery[nVolume][nDir].tilt;
			#ifdef DEBUG
						std::cout << "   Getting Gallery: " << prm.mdGallery[nVolume][nDir].fnImg
								  << " corr=" << cc << " imed=" << imed << " weight=" << weight << " rot=" << angleRot
								  << " tilt=" << angleTilt << std::endl
						          << "Matrix=" << allM[nVolume*Ndirs+nDir] << std::endl
						          << "shiftX=" << shiftX << " shiftY=" << shiftY << std::endl;
			#endif

						size_t recId=mdPartial.addObject();
						mdPartial.setValue(MDL_IMAGE,fnImg,recId);
						mdPartial.setValue(MDL_ENABLED,1,recId);
						mdPartial.setValue(MDL_MAXCC,cc,recId);
						mdPartial.setValue(MDL_COST,imed,recId);
						mdPartial.setValue(MDL_ANGLE_ROT,angleRot,recId);
						mdPartial.setValue(MDL_ANGLE_TILT,angleTilt,recId);
						mdPartial.setValue(MDL_ANGLE_PSI,anglePsi,recId);
						mdPartial.setValue(MDL_SHIFT_X,-shiftX,recId);
						mdPartial.setValue(MDL_SHIFT_Y,-shiftY,recId);
						mdPartial.setValue(MDL_FLIP,flip,recId);
						mdPartial.setValue(MDL_IMAGE_IDX,(size_t)nImg,recId);
						mdPartial.setValue(MDL_REF,(int)nDir,recId);
						mdPartial.setValue(MDL_REF3D,(int)nVolume,recId);
						mdPartial.setValue(MDL_WEIGHT,weight,recId);
					}
				}
			}
#ifdef DEBUG
			std::cout << "Press any key" << std::endl;
			char c; std::cin >> c;
#endif
		}

		if (thArg.thread_id==0)
			progress_bar(nImg+1);
	}
}
#undef DEBUG

// Main routine ------------------------------------------------------------
void ProgReconstructSignificant::run()
{
    show();
    produceSideinfo();

    ThreadManager thMgr(Nthr,this);
    currentAlpha=alpha0;
    double deltaAlpha;
    if (Niter>1)
    	deltaAlpha=(alphaF-alpha0)/(Niter-1);
    else
    	deltaAlpha=0;

    MetaData mdReconstruction, mdPM, mdAux;
    MDRow auxRow;
	size_t Nimgs=mdInp.size();
	Image<double> save;
	MultidimArray<double> ccdir, cdfccdir;
	std::vector<double> ccdirNeighbourhood;
	ccdirNeighbourhood.resize(10*Nimgs);
	bool emptyVolumes=false;
	Matrix1D<double> dir1, dir2;
    for (iter=1; iter<=Niter; iter++)
    {
    	// Generate projections from the different volumes
    	generateProjections();

    	size_t Nvols=mdGallery.size();
    	size_t Ndirs=mdGallery[0].size();
    	cc.initZeros(Nimgs,Nvols,Ndirs);
    	weight=cc;

    	// Align the input images to the projections
    	std::cout << "Current significance: " << 1-currentAlpha << std::endl;
    	std::cerr << "Aligning images ...\n";
    	init_progress_bar(mdIn.size());
    	thMgr.run(progReconstructSignificantThreadAlign);
    	progress_bar(mdIn.size());

    	// Reweight according to each direction
		for (size_t nVolume=0; nVolume<Nvols; ++nVolume)
		{
			for (size_t nDir=0; nDir<Ndirs; ++nDir)
			{
				// Look for the best correlation for this direction
				double bestCorr=-2;
				for (size_t nImg=0; nImg<Nimgs; ++nImg)
				{
					double ccimg=DIRECT_A3D_ELEM(cc,nImg,nVolume,nDir);
					ccdirNeighbourhood.push_back(ccimg);
					if (ccimg>bestCorr)
						bestCorr=ccimg;
				}

				// Look in the neighbourhood
				GalleryImage &g1= mdGallery[nVolume][nDir];
				for (size_t nDir2=0; nDir2<Ndirs; ++nDir2)
				{
					if (nDir!=nDir2)
					{
						GalleryImage &g2= mdGallery[nVolume][nDir2];
						double ang=Euler_distanceBetweenAngleSets(g1.rot,g1.tilt,0.0,g2.rot,g2.tilt,0.0,true);
						if (ang<angDistance)
						{
							for (size_t nImg=0; nImg<Nimgs; ++nImg)
							{
								double ccimg=DIRECT_A3D_ELEM(cc,nImg,nVolume,nDir);
								ccdirNeighbourhood.push_back(ccimg);
								if (ccimg>bestCorr)
									bestCorr=ccimg;
							}
						}
					}
				}

				ccdir=ccdirNeighbourhood;
				ccdir.cumlativeDensityFunction(cdfccdir);

				// Reweight all images for this direction
				double iBestCorr=1.0/bestCorr;
				for (size_t nImg=0; nImg<Nimgs; ++nImg)
				{
					double ccimg=DIRECT_A3D_ELEM(cc,nImg,nVolume,nDir);
					double cdfThis=DIRECT_A1D_ELEM(cdfccdir,nImg);
					DIRECT_A3D_ELEM(weight,nImg,nVolume,nDir)*=ccimg*iBestCorr*cdfThis;
				}
			}
		}

    	// Write the corresponding angular metadata
		for (size_t nVolume=0; nVolume<mdInit.size(); ++nVolume)
		{
			mdReconstruction.clear();
			mdPM.clear();
			for (int thr=0; thr<Nthr; thr++)
			{
				// Readjust weights with direction weights
				MetaData &mdThr=mdReconstructionPartial[nVolume*Nthr+thr];
				FOR_ALL_OBJECTS_IN_METADATA(mdThr)
				{
					size_t nImg;
					int nVol, nDir;
					mdThr.getValue(MDL_IMAGE_IDX,nImg,__iter.objId);
					mdThr.getValue(MDL_REF,nDir,__iter.objId);
					mdThr.getValue(MDL_REF3D,nVol,__iter.objId);
					mdThr.setValue(MDL_WEIGHT,DIRECT_A3D_ELEM(weight,nImg,nVol,nDir),__iter.objId);
				}

				// Union all
				mdReconstruction.unionAll(mdThr);
				mdPM.unionAll(mdProjectionMatching[nVolume*Nthr+thr]);
			}

			// Remove zero-weight images
			mdAux.clear();
			FOR_ALL_OBJECTS_IN_METADATA(mdReconstruction)
			{
				double thisWeight;
				mdReconstruction.getValue(MDL_WEIGHT,thisWeight,__iter.objId);
				if (thisWeight>0)
				{
					mdReconstruction.getRow(auxRow,__iter.objId);
					mdAux.addRow(auxRow);
				}
			}

			String fnAngles=formatString("%s/angles_iter%02d_%02d.xmd",fnDir.c_str(),iter,nVolume);
			if (mdAux.size()>0)
				mdAux.write(fnAngles);
			else
			{
				std::cout << formatString("%s/angles_iter%02d_%02d.xmd is empty. Not written.",fnDir.c_str(),iter,nVolume) << std::endl;
				emptyVolumes=true;
			}
			if (mdPM.size()>0)
			{
				String fnImages=formatString("%s/images_iter%02d_%02d.xmd",fnDir.c_str(),iter,nVolume);
				mdPM.write(fnImages);

				// Remove from mdPM those images that do not participate in angles
				String fnAux=formatString("%s/aux_iter%02d_%02d.xmd",fnDir.c_str(),iter,nVolume);
				// Take only the image name from angles.xmd
				String cmd=(String)"xmipp_metadata_utilities -i "+fnAngles+" --operate keep_column image -o "+fnAux;
				if (system(cmd.c_str())!=0)
					REPORT_ERROR(ERR_MD,(String)"Cannot execute "+cmd);
				// Remove duplicated images
				cmd=(String)"xmipp_metadata_utilities -i "+fnAux+" --operate remove_duplicates image";
				if (system(cmd.c_str())!=0)
					REPORT_ERROR(ERR_MD,(String)"Cannot execute "+cmd);
				// Intersect with images.xmd
				String fnImagesSignificant=formatString("%s/images_significant_iter%02d_%02d.xmd",fnDir.c_str(),iter,nVolume);
				cmd=(String)"xmipp_metadata_utilities -i "+fnImages+" --set intersection "+fnAux+" image image -o "+fnImagesSignificant;
				if (system(cmd.c_str())!=0)
					REPORT_ERROR(ERR_MD,(String)"Cannot execute "+cmd);
		    	deleteFile(fnAux);
			}
			else
				std::cout << formatString("%s/images_iter%02d_%02d.xmd empty. Not written.",fnDir.c_str(),iter,nVolume) << std::endl;
	    	deleteFile(formatString("%s/gallery_iter%02d_%02d_sampling.xmd",fnDir.c_str(),iter,nVolume));
	    	deleteFile(formatString("%s/gallery_iter%02d_%02d.doc",fnDir.c_str(),iter,nVolume));
	    	deleteFile(formatString("%s/gallery_iter%02d_%02d.stk",fnDir.c_str(),iter,nVolume));
		}

		if (verbose>=2)
		{
			save()=cc;
			save.write(formatString("%s/cc_iter%02d.vol",fnDir.c_str(),iter));
			save()=weight;
			save.write(formatString("%s/weight_iter%02d.vol",fnDir.c_str(),iter));
		}

    	// Reconstruct
    	reconstructCurrent();

    	currentAlpha+=deltaAlpha;

    	if (emptyVolumes)
    		break;
    }
}

void ProgReconstructSignificant::reconstructCurrent()
{
	std::cerr << "Reconstructing volumes ..." << std::endl;
	for (size_t nVolume=0; nVolume<mdInit.size(); ++nVolume)
	{
		FileName fnAngles=formatString("%s/angles_iter%02d_%02d.xmd",fnDir.c_str(),iter,nVolume);
		if (!fnAngles.exists())
			continue;
		FileName fnVolume=formatString("%s/volume_iter%02d_%02d.vol",fnDir.c_str(),iter,nVolume);
		String args=formatString("-i %s -o %s --sym %s --weight --thr %d -v 0",fnAngles.c_str(),fnVolume.c_str(),fnSym.c_str(),Nthr);
		String cmd=(String)"xmipp_reconstruct_fourier "+args;
		if (system(cmd.c_str())==-1)
			REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");

		args=formatString("-i %s --mask circular %d -v 0",fnVolume.c_str(),-XSIZE(inputImages())/2);
		cmd=(String)"xmipp_transform_mask "+args;
		if (system(cmd.c_str())==-1)
			REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
	}
}

void ProgReconstructSignificant::generateProjections()
{
	int Nvol=(int)mdInit.size();
	FileName fnVol, fnGallery, fnAngles, fnGalleryMetaData;
	std::vector<GalleryImage> galleryNames;
	mdGallery.clear();
	for (int n=0; n<Nvol; n++)
	{
		fnVol=formatString("%s/volume_iter%02d_%02d.vol",fnDir.c_str(),iter-1,n);
		fnGallery=formatString("%s/gallery_iter%02d_%02d.stk",fnDir.c_str(),iter,n);
		fnAngles=formatString("%s/angles_iter%02d_%02d.xmd",fnDir.c_str(),iter-1,n);
		fnGalleryMetaData=formatString("%s/gallery_iter%02d_%02d.doc",fnDir.c_str(),iter,n);
		String args=formatString("-i %s -o %s --sampling_rate %f --sym %s --compute_neighbors --angular_distance -1 --experimental_images %s --min_tilt_angle %f --max_tilt_angle %f -v 0",
				fnVol.c_str(),fnGallery.c_str(),angularSampling,fnSym.c_str(),fnAngles.c_str(),tilt0,tiltF);

		String cmd=(String)"xmipp_angular_project_library "+args;
		if (system(cmd.c_str())==-1)
			REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");

		mdGallery.push_back(galleryNames);
		MetaData mdAux(fnGalleryMetaData);
		galleryNames.clear();
		FOR_ALL_OBJECTS_IN_METADATA(mdAux)
		{
			GalleryImage I;
			mdAux.getValue(MDL_IMAGE,I.fnImg,__iter.objId);
			mdAux.getValue(MDL_ANGLE_ROT,I.rot,__iter.objId);
			mdAux.getValue(MDL_ANGLE_TILT,I.tilt,__iter.objId);
			mdGallery[n].push_back(I);
		}
		gallery[n].read(fnGallery);
	}
}

void ProgReconstructSignificant::produceSideinfo()
{
	mdIn.read(fnIn);
	mdIn.removeDisabled();

	mdIn.fillConstant(MDL_MAXCC,"0.0");
	mdIn.fillConstant(MDL_ANGLE_ROT,"0.0");
	mdIn.fillConstant(MDL_ANGLE_TILT,"0.0");
	mdIn.fillConstant(MDL_ANGLE_PSI,"0.0");
	mdIn.fillConstant(MDL_SHIFT_X,"0.0");
	mdIn.fillConstant(MDL_SHIFT_Y,"0.0");
	mdIn.fillConstant(MDL_WEIGHT,"1.0");

	// Read all input images in memory
	size_t xdim, ydim, zdim, ndim;
	getImageSize(mdIn,xdim, ydim, zdim, ndim);
	inputImages().resizeNoCopy(mdIn.size(), 1, ydim, xdim);

	Image<double> I;
	size_t n=0;
	MultidimArray<double> mCurrentImage;
	FileName fnImg;
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		mdIn.getValue(MDL_IMAGE,fnImg,__iter.objId);
		mdInp.push_back(fnImg);
		I.read(fnImg);
		mCurrentImage.aliasImageInStack(inputImages(),n++);
		memcpy(MULTIDIM_ARRAY(mCurrentImage),MULTIDIM_ARRAY(I()),MULTIDIM_SIZE(mCurrentImage)*sizeof(double));
	}

	// Copy all input values as iteration 0 volumes
	mdInit.read(fnInit);
	FileName fnVol, fnAngles;
	Image<double> V, galleryDummy;
	int idx=0;
	FOR_ALL_OBJECTS_IN_METADATA(mdInit)
	{
		mdInit.getValue(MDL_IMAGE,fnVol,__iter.objId);
		V.read(fnVol);
		fnVol=formatString("%s/volume_iter00_%02d.vol",fnDir.c_str(),idx);
		V.write(fnVol);
		mdInit.setValue(MDL_IMAGE,fnVol,__iter.objId);
		fnAngles=formatString("%s/angles_iter00_%02d.xmd",fnDir.c_str(),idx);
		mdIn.write(fnAngles);
		gallery.push_back(galleryDummy);

		// Create a partial metadata for each thread and each volume
		for (int thr=0; thr<Nthr; ++thr)
		{
			mdReconstructionPartial.push_back(MetaData());
			mdProjectionMatching.push_back(MetaData());
		}
		idx++;
	}

	iter=0;
}
