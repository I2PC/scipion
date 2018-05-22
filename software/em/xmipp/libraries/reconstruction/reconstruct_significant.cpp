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
ProgReconstructSignificant::ProgReconstructSignificant()
{
	rank=0;
	Nprocessors=1;
	randomize_random_generator();
	deltaAlpha2=0;
}

void ProgReconstructSignificant::defineParams()
{
    //usage
    addUsageLine("Generate 3D reconstructions from projections using random orientations");
    //params
    addParamsLine("   -i <md_file>                : Metadata file with input projections");
    addParamsLine("  [--numberOfVolumes <N=1>]    : Number of volumes to reconstruct");
    addParamsLine("  [--initvolumes <md_file=\"\">] : Set of initial volumes. If none is given, a single volume");
    addParamsLine("                               : is reconstructed with random angles assigned to the images");
    addParamsLine("  [--initgallery <md_file=\"\">]   : Gallery of projections from a single volume");
    addParamsLine("  [--odir <outputDir=\".\">]   : Output directory");
    addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections");
    addParamsLine("  [--iter <N=10>]              : Number of iterations");
    addParamsLine("  [--alpha0 <N=0.05>]          : Initial significance");
    addParamsLine("  [--alphaF <N=0.005>]         : Final significance");
    addParamsLine("  [--keepIntermediateVolumes]  : Keep the volume of each iteration");
    addParamsLine("  [--angularSampling <a=5>]    : Angular sampling in degrees for generating the projection gallery");
    addParamsLine("  [--maxShift <s=-1>]          : Maximum shift allowed (+-this amount)");
    addParamsLine("  [--minTilt <t=0>]            : Minimum tilt angle");
    addParamsLine("  [--maxTilt <t=90>]           : Maximum tilt angle");
    addParamsLine("  [--useImed]                  : Use Imed for weighting");
    addParamsLine("  [--strictDirection]          : Images not significant for a direction are also discarded");
    addParamsLine("  [--angDistance <a=10>]       : Angular distance");
    addParamsLine("  [--dontApplyFisher]          : Do not select directions using Fisher");
    addParamsLine("  [--dontReconstruct]          : Do not reconstruct");
    addParamsLine("  [--useForValidation <numOrientationsPerParticle=10>] : Use the program for validation. This number defines the number of possible orientations per particle");
    addParamsLine("  [--dontCheckMirrors]         : Don't check mirrors in the alignment process");

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
    keepIntermediateVolumes = checkParam("--keepIntermediateVolumes");
    angularSampling=getDoubleParam("--angularSampling");
    maxShift=getDoubleParam("--maxShift");
    tilt0=getDoubleParam("--minTilt");
    tiltF=getDoubleParam("--maxTilt");
    useImed=checkParam("--useImed");
    strict=checkParam("--strictDirection");
    angDistance=getDoubleParam("--angDistance");
    Nvolumes=getIntParam("--numberOfVolumes");
    applyFisher=checkParam("--dontApplyFisher");
    fnFirstGallery=getParam("--initgallery");
    doReconstruct=!checkParam("--dontReconstruct");
    useForValidation=checkParam("--useForValidation");
    numOrientationsPerParticle = getIntParam("--useForValidation");
    dontCheckMirrors = checkParam("--dontCheckMirrors");

    if (!doReconstruct)
    {
    	if (rank==0)
    		std::cout << "Setting iterations to 1 because of --dontReconstruct\n";
    	Niter=1; // Make sure that there is a single iteration
    }
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
        std::cout << "Keep intermediate volumes   : "  << keepIntermediateVolumes << std::endl;
        std::cout << "Angular sampling            : "  << angularSampling << std::endl;
        std::cout << "Maximum shift               : "  << maxShift << std::endl;
        std::cout << "Minimum tilt                : "  << tilt0 << std::endl;
        std::cout << "Maximum tilt                : "  << tiltF << std::endl;
        std::cout << "Use Imed                    : "  << useImed << std::endl;
        std::cout << "Angular distance            : "  << angDistance << std::endl;
        std::cout << "Apply Fisher                : "  << applyFisher << std::endl;
        std::cout << "Reconstruct                 : "  << doReconstruct << std::endl;
        std::cout << "useForValidation            : "  << useForValidation << std::endl;
        std::cout << "dontCheckMirrors            : "  << dontCheckMirrors << std::endl;


        if (fnSym != "")
            std::cout << "Symmetry for projections    : "  << fnSym << std::endl;
        if (fnFirstGallery=="")
        {
			if (fnInit !="")
				std::cout << "Initial volume              : "  << fnInit << std::endl;
			else
				std::cout << "Number of volumes           : "  << Nvolumes << std::endl;
        }
        else
            std::cout <<     "Gallery                     : "  << fnFirstGallery << std::endl;

        if (useForValidation)
        	std::cout <<  " numOrientationsPerParticle : " << numOrientationsPerParticle << std::endl;
    }
}

// Image alignment ========================================================
//#define DEBUG
void ProgReconstructSignificant::alignImagesToGallery()
{
    AlignmentAux aux;
    CorrelationAux aux2;
    RotationalCorrelationAux aux3;

	Matrix2D<double> M;
	std::vector< Matrix2D<double> > allM;

	size_t Nvols=YSIZE(cc);
	size_t Ndirs=XSIZE(cc);

	// Clear the previous assignment
	for (size_t nvol=0; nvol<Nvols; ++nvol)
	{
		mdReconstructionPartial[nvol].clear();
		mdReconstructionProjectionMatching[nvol].clear();
	}

	MultidimArray<double> imgcc(Nvols*Ndirs), imgimed(Nvols*Ndirs);
	MultidimArray<double> cdfcc, cdfimed;
	double one_alpha=1-currentAlpha-deltaAlpha2;

	FileName fnImg;
	size_t nImg=0;
	Image<double> I;
	MultidimArray<double> mCurrentImageAligned, mGalleryProjection;
	if (rank==0)
	{
		std::cout << "Current significance: " << one_alpha << std::endl;
		std::cerr << "Aligning images ...\n";
		init_progress_bar(mdIn.size());
	}

	MultidimArray<double> ccVol(Nvols);
	MDRow row;
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		if ((nImg+1)%Nprocessors==rank)
		{
			mdIn.getValue(MDL_IMAGE,fnImg,__iter.objId);
			mdIn.getRow(row,__iter.objId);
#ifdef DEBUG
			std::cout << "Processing: " << fnImg << std::endl;
#endif
			I.read(fnImg);
			MultidimArray<double> &mCurrentImage=I();
			mCurrentImage.setXmippOrigin();
			allM.clear();

			double bestCorr=-2, bestRot, bestTilt, bestImed=1e38, worstImed=-1e38;
			Matrix2D<double> bestM;
			int bestVolume=-1;

			// Compute all correlations
	    	for (size_t nVolume=0; nVolume<Nvols; ++nVolume)
	    	{
	    		AlignmentTransforms *transforms=galleryTransforms[nVolume];
		    	for (size_t nDir=0; nDir<Ndirs; ++nDir)
				{
					mCurrentImageAligned=mCurrentImage;
					mGalleryProjection.aliasImageInStack(gallery[nVolume](),nDir);
					mGalleryProjection.setXmippOrigin();
					double corr;
					if (! dontCheckMirrors)
						corr=alignImagesConsideringMirrors(mGalleryProjection,transforms[nDir],
								mCurrentImageAligned,M,aux,aux2,aux3,DONT_WRAP);
					else
						corr = alignImages(mGalleryProjection, mCurrentImageAligned,
						                   M, DONT_WRAP);

//					double corr=alignImagesConsideringMirrors(mGalleryProjection,
//							mCurrentImageAligned,M,aux,aux2,aux3,DONT_WRAP);
					M=M.inv();
					double scale, shiftX, shiftY, anglePsi;
					bool flip;
					transformationMatrix2Parameters2D(M,flip,scale,shiftX,shiftY,anglePsi);

					double imed=imedDistance(mGalleryProjection, mCurrentImageAligned);
					if (maxShift>0 && (fabs(shiftX)>maxShift || fabs(shiftY)>maxShift))
					{
						corr/=3;
						imed*=3;
					}

//					//if (corr>0.99)
//					//{
//					//std::cout << prm.mdGallery[nVolume][nDir].fnImg << " corr= " << corr << " imed= " << imed << std::endl;
//					//std::cout << "Matrix=" << M << std::endl;
//					//Image<double> save;
//					//save()=mGalleryProjection;
//					//save.write("PPPgallery.xmp");
//					//save()=mCurrentImage;
//					//save.write("PPPcurrentImage.xmp");
//					//save()=mCurrentImageAligned;
//					//save.write("PPPcurrentImageAligned.xmp");
//					//char c; std::cin >> c;
//					//}

					DIRECT_A3D_ELEM(cc,nImg,nVolume,nDir)=corr;
					// For the paper plot: std::cout << corr << " " << imed << std::endl;
					size_t idx=nVolume*Ndirs+nDir;
					DIRECT_A1D_ELEM(imgcc,idx)=corr;
					DIRECT_A1D_ELEM(imgimed,idx)=imed;
					allM.push_back(M);

					if (corr>bestCorr)
					{
						bestM=M;
						bestCorr=corr;
						bestVolume=(int)nVolume;
						bestRot=mdGallery[nVolume][nDir].rot;
						bestTilt=mdGallery[nVolume][nDir].tilt;
						// std::cout << "nDir=" << nDir << " bestCorr=" << bestCorr << " imed=" << imed << " (bestImed=" << bestImed << ") M=" << M << std::endl;
                    }

					if (imed<bestImed)
						bestImed=imed;
					else if (imed>worstImed)
						worstImed=imed;
				}
	    	}

	    	// Keep the best assignment for the projection matching
	    	// Each process keeps a list of the images for each volume
			MetaData &mdProjectionMatching=mdReconstructionProjectionMatching[bestVolume];
			double scale, shiftX, shiftY, anglePsi;
			bool flip;
			transformationMatrix2Parameters2D(bestM,flip,scale,shiftX,shiftY,anglePsi);
			if (useForValidation && dontCheckMirrors)
				flip = false;

			if (maxShift<0 || (maxShift>0 && fabs(shiftX)<maxShift && fabs(shiftY)<maxShift))
			{
				size_t recId=mdProjectionMatching.addRow(row);
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

			// Compute the cumulative distributions
			imgcc.cumlativeDensityFunction(cdfcc);
			imgimed.cumlativeDensityFunction(cdfimed);

			// Check if force 1 volume
			size_t firstVolume=0;
			size_t lastVolume=Nvols-1;

	    	// Get the best images
			for (size_t nVolume=firstVolume; nVolume<=lastVolume; ++nVolume)
			{
				MetaData &mdPartial=mdReconstructionPartial[nVolume];
				for (size_t nDir=0; nDir<Ndirs; ++nDir)
				{
					size_t idx=nVolume*Ndirs+nDir;
					double cdfccthis=DIRECT_A1D_ELEM(cdfcc,idx);
					double cdfimedthis=DIRECT_A1D_ELEM(cdfimed,idx);
					double cc=DIRECT_A1D_ELEM(imgcc,idx);
					// bool condition=!useImed || (useImed && cdfimedthis<=currentAlpha);
//					if (cc>ccl)
//						std::cout << "Image " << nImg << " " << fnImg << " qualifies by Fisher to dir=" << nDir << std::endl;
//					if (!(cdfccthis>=one_alpha) && cc>ccl)
//						std::cout << "Image " << nImg << " " << fnImg << " does not qualify by correlation percentile to " << nDir << " -> " << cdfccthis << " " << one_alpha << std::endl;
//					if (!condition && cc>ccl)
//						std::cout << "Image " << nImg << " " << fnImg << " does not qualify by imed percentile to " << nDir << " -> " << cdfimedthis << " " << currentAlpha<< std::endl;
					bool condition=true;
					condition=condition && ((applyFisher && cc>=ccl) || !applyFisher);
					condition=condition && cdfccthis>=one_alpha;
					if (condition)
					{
						// std::cout << fnImg << " is selected for dir=" << nDir << std::endl;
						double imed=DIRECT_A1D_ELEM(imgimed,idx);
						transformationMatrix2Parameters2D(allM[nVolume*Ndirs+nDir],flip,scale,shiftX,shiftY,anglePsi);
						if (useForValidation && dontCheckMirrors)
							flip = false;

						if (maxShift>0)
							if (fabs(shiftX)>maxShift || fabs(shiftY)>maxShift)
								continue;
						if (flip)
							shiftX*=-1;

						double thisWeight=cdfccthis*(cc/bestCorr);
						// COSS: To promote sparsity in the volume assignment: sum_i(cc_i^p)/sum_i(cc_i)*cc_i^p/cc_i
						if (useImed)
							thisWeight*=(1-cdfimedthis)*(bestImed/imed);
						DIRECT_A3D_ELEM(weight,nImg,nVolume,nDir)=thisWeight;
						double angleRot=mdGallery[nVolume][nDir].rot;
						double angleTilt=mdGallery[nVolume][nDir].tilt;
			#ifdef DEBUG
						std::cout << "   Getting Gallery: " << mdGallery[nVolume][nDir].fnImg
								  << " corr=" << cc << " imed=" << imed << " bestImed=" << bestImed << " weight=" << thisWeight << " rot=" << angleRot
								  << " tilt=" << angleTilt << std::endl
								  << " weight by corr=" << cdfccthis*(cc/bestCorr) << " weight by imed=" << (1-cdfimedthis)*(bestImed/imed) << std::endl
						          << "Matrix=" << allM[nVolume*Ndirs+nDir] << std::endl
						          << "shiftX=" << shiftX << " shiftY=" << shiftY << std::endl;
			#endif

						size_t recId=mdPartial.addRow(row);
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
						mdPartial.setValue(MDL_WEIGHT,thisWeight,recId);
						mdPartial.setValue(MDL_WEIGHT_SIGNIFICANT,thisWeight,recId);
					}
				}
			}
#ifdef DEBUG
			std::cout << "Press any key" << std::endl;
			char c; std::cin >> c;
#endif
		}

		if (rank==0)
			progress_bar(nImg+1);
		nImg++;
	}
	if (rank==0)
		progress_bar(mdIn.size());
}
#undef DEBUG

// Main routine ------------------------------------------------------------
void ProgReconstructSignificant::run()
{
	if (rank==0)
		show();
    produceSideinfo();

    /*currentAlpha=alpha0;
    double deltaAlpha;
    if (Niter>1)
    	deltaAlpha=(alphaF-alpha0)/(Niter-1);
    else
    	deltaAlpha=0;
*/
    MetaData mdAux;
	size_t Nimgs=mdIn.size();
	Image<double> save;
	MultidimArray<double> ccdir, cdfccdir;
	std::vector<double> ccdirNeighbourhood;
	ccdirNeighbourhood.resize(10*Nimgs);
	bool emptyVolumes=false;
	Matrix1D<double> dir1, dir2;
    for (iter=1; iter<=Niter; iter++)
    {
    	if (rank==0)
    		std::cout << "\nIteration " << iter << std::endl;
    	// Generate projections from the different volumes
    	generateProjections();
    	if (useForValidation)
    		numberOfProjections();

	    currentAlpha=alpha0;
	    double deltaAlpha;
	    if (Niter>1)
	    	deltaAlpha=(alphaF-alpha0)/(Niter-1);
	    else
	    	deltaAlpha=0;

    	size_t Nvols=mdGallery.size();
    	size_t Ndirs=mdGallery[0].size();
    	cc.initZeros(Nimgs,Nvols,Ndirs);
    	weight=cc;
    	double oneAlpha=1-currentAlpha-deltaAlpha2;

    	// Align the input images to the projections
    	alignImagesToGallery();
    	synchronize();
    	gatherAlignment();

    	// Reweight according to each direction
    	if (rank==0)
    	{
    		std::cerr << "Readjusting weights ..." << std::endl;
    		init_progress_bar(Nvols);
			for (size_t nVolume=0; nVolume<Nvols; ++nVolume)
			{
				for (size_t nDir=0; nDir<Ndirs; ++nDir)
				{
					// Look for the best correlation for this direction
					double bestCorr=-2;
					ccdirNeighbourhood.clear();
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
						if ((cdfThis>=oneAlpha || !strict) && DIRECT_A3D_ELEM(weight,nImg,nVolume,nDir)>0)
						{
							// std::cout << "Neighborhood " << nDir << " accepts image " << nImg << " with weight " << DIRECT_A3D_ELEM(weight,nImg,nVolume,nDir) << "*" << ccimg << "*" << iBestCorr << "*" << cdfThis << "=" << DIRECT_A3D_ELEM(weight,nImg,nVolume,nDir)*ccimg*iBestCorr*cdfThis << std::endl;
							DIRECT_A3D_ELEM(weight,nImg,nVolume,nDir)*=ccimg*iBestCorr*cdfThis;
						}
						else
							DIRECT_A3D_ELEM(weight,nImg,nVolume,nDir)=0;
					}
				}
				progress_bar(nVolume);
			}
			progress_bar(Nvols);

			// Write the corresponding angular metadata
			for (size_t nVolume=0; nVolume<(size_t)Nvolumes; ++nVolume)
			{
				MetaData &mdReconstruction=mdReconstructionPartial[nVolume];
				// Readjust weights with direction weights
				FOR_ALL_OBJECTS_IN_METADATA(mdReconstruction)
				{
					size_t nImg;
					int nVol, nDir;
					mdReconstruction.getValue(MDL_IMAGE_IDX,nImg,__iter.objId);
					mdReconstruction.getValue(MDL_REF,nDir,__iter.objId);
					mdReconstruction.getValue(MDL_REF3D,nVol,__iter.objId);
					mdReconstruction.setValue(MDL_WEIGHT,DIRECT_A3D_ELEM(weight,nImg,nVol,nDir),__iter.objId);
					mdReconstruction.setValue(MDL_WEIGHT_SIGNIFICANT,DIRECT_A3D_ELEM(weight,nImg,nVol,nDir),__iter.objId);
					//if (DIRECT_A3D_ELEM(weight,nImg,nVol,nDir)==0)
					//	std::cout << "Direction " << nDir << " does not accept img " << nImg << std::endl;
				}

				// Remove zero-weight images
				mdAux.clear();
				if (mdReconstruction.size()>0)
					mdAux.importObjects(mdReconstruction,MDValueGT(MDL_WEIGHT,0.0));

				String fnAngles=formatString("%s/angles_iter%03d_%02d.xmd",fnDir.c_str(),iter,nVolume);
				if (mdAux.size()>0)
					mdAux.write(fnAngles);
				else
				{
					std::cout << formatString("%s/angles_iter%03d_%02d.xmd is empty. Not written.",fnDir.c_str(),iter,nVolume) << std::endl;
					emptyVolumes=true;
				}

				MetaData &mdPM=mdReconstructionProjectionMatching[nVolume];
				if (mdPM.size()>0 && fileExists(fnAngles.c_str()))
				{
					String fnImages=formatString("%s/images_iter%03d_%02d.xmd",fnDir.c_str(),iter,nVolume);
					mdPM.write(fnImages);

					// Remove from mdPM those images that do not participate in angles
					String fnAux=formatString("%s/aux_iter%03d_%02d.xmd",fnDir.c_str(),iter,nVolume);
					// Take only the image name from angles.xmd
					String cmd=(String)"xmipp_metadata_utilities -i "+fnAngles+" --operate keep_column image -o "+fnAux;
					if (system(cmd.c_str())!=0)
						REPORT_ERROR(ERR_MD,(String)"Cannot execute "+cmd);
					// Remove duplicated images
					cmd=(String)"xmipp_metadata_utilities -i "+fnAux+" --operate remove_duplicates image";
					if (system(cmd.c_str())!=0)
						REPORT_ERROR(ERR_MD,(String)"Cannot execute "+cmd);
					// Intersect with images.xmd
					String fnImagesSignificant=formatString("%s/images_significant_iter%03d_%02d.xmd",fnDir.c_str(),iter,nVolume);
					cmd=(String)"xmipp_metadata_utilities -i "+fnImages+" --set intersection "+fnAux+" image image -o "+fnImagesSignificant;
					if (system(cmd.c_str())!=0)
						REPORT_ERROR(ERR_MD,(String)"Cannot execute "+cmd);
					deleteFile(fnAux);
				}
				else
					std::cout << formatString("%s/images_iter%03d_%02d.xmd empty. Not written.",fnDir.c_str(),iter,nVolume) << std::endl;
				deleteFile(formatString("%s/gallery_iter%03d_%02d_sampling.xmd",fnDir.c_str(),iter,nVolume));
				deleteFile(formatString("%s/gallery_iter%03d_%02d.doc",fnDir.c_str(),iter,nVolume));
				deleteFile(formatString("%s/gallery_iter%03d_%02d.stk",fnDir.c_str(),iter,nVolume));
				if (iter>=1 && !keepIntermediateVolumes)
				{
					deleteFile(formatString("%s/volume_iter%03d_%02d.vol",fnDir.c_str(),iter-1,nVolume));
					deleteFile(formatString("%s/images_iter%03d_%02d.xmd",fnDir.c_str(),iter-1,nVolume));
					deleteFile(formatString("%s/angles_iter%03d_%02d.xmd",fnDir.c_str(),iter-1,nVolume));
					deleteFile(formatString("%s/images_significant_iter%03d_%02d.xmd",fnDir.c_str(),iter-1,nVolume));
				}
			}

			if (verbose>=2)
			{
				save()=cc;
				save.write(formatString("%s/cc_iter%03d.vol",fnDir.c_str(),iter));
				save()=weight;
				save.write(formatString("%s/weight_iter%03d.vol",fnDir.c_str(),iter));
			}
    	}
    	synchronize();

    	// Reconstruct
    	if (doReconstruct)
    		reconstructCurrent();
    	else
    		break;

    	currentAlpha+=deltaAlpha;

    	if (emptyVolumes)
    		break;
    }
}

void ProgReconstructSignificant::reconstructCurrent()
{
	if (rank==0)
		std::cerr << "Reconstructing volumes ..." << std::endl;
	MetaData MD;
	for (size_t nVolume=0; nVolume<(size_t)Nvolumes; ++nVolume)
	{
		if ((nVolume+1)%Nprocessors!=rank)
			continue;

		FileName fnAngles=formatString("%s/angles_iter%03d_%02d.xmd",fnDir.c_str(),iter,nVolume);
		if (!fnAngles.exists())
			continue;
		MD.read(fnAngles);
		std::cout << "Volume " << nVolume << ": number of images=" << MD.size() << std::endl;
		FileName fnVolume=formatString("%s/volume_iter%03d_%02d.vol",fnDir.c_str(),iter,nVolume);
		String args=formatString("-i %s -o %s --sym %s --weight -v 0",fnAngles.c_str(),fnVolume.c_str(),fnSym.c_str());
		String cmd=(String)"xmipp_reconstruct_fourier "+args;
		std::cout << cmd << std::endl;
		if (system(cmd.c_str())==-1)
			REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");

		if (fnSym!="c1")
		{
			args=formatString("-i %s --sym %s -v 0",fnVolume.c_str(),fnSym.c_str());
			cmd=(String)"xmipp_transform_symmetrize "+args;
			std::cout << cmd << std::endl;
			if (system(cmd.c_str())==-1)
				REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
		}

		args=formatString("-i %s --mask circular %d -v 0",fnVolume.c_str(),-Xdim/2);
		cmd=(String)"xmipp_transform_mask "+args;
		std::cout << cmd << std::endl;
		if (system(cmd.c_str())==-1)
			REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
	}
}

void ProgReconstructSignificant::generateProjections()
{
	FileName fnGallery, fnGalleryMetaData;
	if (iter>1 || fnFirstGallery=="")
	{
		FileName fnVol, fnAngles;
		// Generate projections
		for (int n=0; n<Nvolumes; n++)
		{
			if ((n+1)%Nprocessors!=rank)
				continue;
			fnVol=formatString("%s/volume_iter%03d_%02d.vol",fnDir.c_str(),iter-1,n);
			fnGallery=formatString("%s/gallery_iter%03d_%02d.stk",fnDir.c_str(),iter,n);
			fnAngles=formatString("%s/angles_iter%03d_%02d.xmd",fnDir.c_str(),iter-1,n);
			fnGalleryMetaData=formatString("%s/gallery_iter%03d_%02d.doc",fnDir.c_str(),iter,n);
			String args=formatString("-i %s -o %s --sampling_rate %f --sym %s --compute_neighbors --angular_distance -1 --experimental_images %s --min_tilt_angle %f --max_tilt_angle %f -v 0",
					fnVol.c_str(),fnGallery.c_str(),angularSampling,fnSym.c_str(),fnAngles.c_str(),tilt0,tiltF);

			String cmd=(String)"xmipp_angular_project_library "+args;
			if (system(cmd.c_str())==-1)
				REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
		}
		synchronize();
	}

	// Read projection galleries
	std::vector<GalleryImage> galleryNames;
	mdGallery.clear();

	CorrelationAux aux;
	AlignmentAux aux2;
	MultidimArray<double> mGalleryProjection;
	for (int n=0; n<Nvolumes; n++)
	{
		mdGallery.push_back(galleryNames);
		if (iter>1 || fnFirstGallery=="")
		{
			fnGalleryMetaData=formatString("%s/gallery_iter%03d_%02d.doc",fnDir.c_str(),iter,n);
			fnGallery=formatString("%s/gallery_iter%03d_%02d.stk",fnDir.c_str(),iter,n);
		}
		else
		{
			fnGalleryMetaData=fnFirstGallery;
			fnGallery=fnFirstGallery.replaceExtension("stk");
		}
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

		// Calculate transforms of this gallery
		size_t kmax=NSIZE(gallery[n]());
		if (galleryTransforms[n]==NULL)
		{
			delete galleryTransforms[n];
			galleryTransforms[n]=new AlignmentTransforms[kmax];
		}
		AlignmentTransforms *transforms=galleryTransforms[n];
		for (size_t k=0; k<kmax; ++k)
		{
			mGalleryProjection.aliasImageInStack(gallery[n](),k);
			mGalleryProjection.setXmippOrigin();
			aux.transformer1.FourierTransform((MultidimArray<double> &)mGalleryProjection, transforms[k].FFTI, true);
		    normalizedPolarFourierTransform(mGalleryProjection, transforms[k].polarFourierI, false,
		                                    XSIZE(mGalleryProjection) / 5, XSIZE(mGalleryProjection) / 2, aux2.plans, 1);
		}
	}
}

void ProgReconstructSignificant::numberOfProjections()
{

	double number_of_projections = (double)mdGallery[0].size();
/*
	double angDist[] =        {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0};
	double numWithOutSymm[] = {165016, 41258, 18338, 10318, 6600, 4586, 3367, 2586, 2042, 1652, 1367, 1148, 977, 843, 732, 643, 569, 510, 460, 412, 340, 288, 245, 211, 184, 161, 146};
	SymList SL;
	SL.readSymmetryFile(fnSym);
	size_t minIndex = 0;
	double tmp = 1e3;
	double numProjects = 0;
	for (int idx = 0; idx < 27; idx++)
	{
		if ( std::abs(angularSampling-angDist[idx]) < tmp)
		{
			minIndex = idx;
			numProjects = (numWithOutSymm[idx]/((2*(SL.symsNo()+1))));
			tmp = std::abs(angularSampling-angDist[idx]);
		}
	}
	angularSampling = angDist[minIndex];
*/

	alpha0 = numOrientationsPerParticle/number_of_projections;
	alphaF = alpha0;
	deltaAlpha2 = 1/(2*number_of_projections);

    if (rank==0)
    {
      std::cout << " Changing angular sampling to " << angularSampling << std::endl;
      std::cout << " Number of projections taking into account angular sampling and symmetry " << number_of_projections << std::endl;
      std::cout << " changing alpha0 to " << alpha0 << std::endl;
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

	size_t Ydim,Zdim,Ndim;
	getImageSize(mdIn,Xdim,Ydim,Zdim,Ndim);

	// Adjust alpha
	if ( (fnSym!="c1") && !useForValidation )
	{
		SymList SL;
		SL.readSymmetryFile(fnSym);
		alpha0*=SL.true_symNo;
		alphaF*=SL.true_symNo;
		if (alpha0>1)
			REPORT_ERROR(ERR_ARG_INCORRECT,"Alpha values are too large: reduce the error such that the error times the symmetry number is smaller than 1");
	}
	// If there is not any input volume, create a random one
	if (fnFirstGallery=="")
	{
		bool deleteInit=false;
		if (fnInit=="")
		{
			deleteInit=true;
			MetaData mdAux;
			for (int n=0; n<Nvolumes; ++n)
			{
				fnInit=fnDir+"/volume_random.xmd";
				if (rank==0)
				{
					MetaData mdRandom;
					mdRandom=mdIn;
					FOR_ALL_OBJECTS_IN_METADATA(mdRandom)
					{
						mdRandom.setValue(MDL_ANGLE_ROT,rnd_unif(0,360),__iter.objId);
						mdRandom.setValue(MDL_ANGLE_TILT,rnd_unif(0,180),__iter.objId);
						mdRandom.setValue(MDL_ANGLE_PSI,rnd_unif(0,360),__iter.objId);
						mdRandom.setValue(MDL_SHIFT_X,0.0,__iter.objId);
						mdRandom.setValue(MDL_SHIFT_Y,0.0,__iter.objId);
					}
					FileName fnAngles=fnDir+formatString("/angles_random_%02d.xmd",n);
					FileName fnVolume=fnDir+formatString("/volume_random_%02d.vol",n);
					mdRandom.write(fnAngles);
					String args=formatString("-i %s -o %s --sym %s -v 0",fnAngles.c_str(),fnVolume.c_str(),fnSym.c_str());
					String cmd=(String)"xmipp_reconstruct_fourier "+args;
					if (system(cmd.c_str())==-1)
						REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
					if (!keepIntermediateVolumes)
						deleteFile(fnAngles);

					// Symmetrize with many different possibilities to have a spherical volume
					args=formatString("-i %s --sym i1 --spline 1 -v 0",fnVolume.c_str());
					cmd=(String)"xmipp_transform_symmetrize "+args;
					if (system(cmd.c_str())==-1)
						REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");

					args=formatString("-i %s --sym i3 --spline 1 -v 0",fnVolume.c_str());
					if (system(cmd.c_str())==-1)
						REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");

					args=formatString("-i %s --sym i2 --spline 1 -v 0",fnVolume.c_str());
					if (system(cmd.c_str())==-1)
						REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
					deleteFile(fnAngles);
					mdAux.setValue(MDL_IMAGE,fnVolume,mdAux.addObject());
				}
			}
			if (rank==0)
				mdAux.write(fnInit);
			synchronize();
		}

		// Copy all input values as iteration 0 volumes
		MetaData mdInit;
		mdInit.read(fnInit);
		FileName fnVol;
		Image<double> V;
		int idx=0;
		FOR_ALL_OBJECTS_IN_METADATA(mdInit)
		{
			if (rank==0)
			{
				mdInit.getValue(MDL_IMAGE,fnVol,__iter.objId);
				V.read(fnVol);
				fnVol=formatString("%s/volume_iter000_%02d.vol",fnDir.c_str(),idx);
				V.write(fnVol);
				mdInit.setValue(MDL_IMAGE,fnVol,__iter.objId);
			}
			idx++;
		}
		Nvolumes=(int)mdInit.size();

		iter=0;
		synchronize();
		if (deleteInit && rank==0)
			deleteFile(fnInit);
	}
	else
		Nvolumes=1;
	synchronize();

	// Copy all input values as iteration 0 volumes
	FileName fnAngles;
	Image<double> galleryDummy;
	MetaData mdPartial, mdProjMatch;
	for (int idx=0; idx<Nvolumes; ++idx)
	{
		if (rank==0)
		{
			fnAngles=formatString("%s/angles_iter000_%02d.xmd",fnDir.c_str(),idx);
			mdIn.write(fnAngles);
		}
		gallery.push_back(galleryDummy);
		galleryTransforms.push_back(NULL);
		mdReconstructionPartial.push_back(mdPartial);
		mdReconstructionProjectionMatching.push_back(mdProjMatch);
	}

	iter=0;
	synchronize();
}
