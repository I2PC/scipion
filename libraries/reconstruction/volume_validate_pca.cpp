/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
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


#include "volume_validate_pca.h"

// Define params
void ProgVolumeValidationPCA::defineParams()
{
    //usage
    addUsageLine("Validate an obtained volume from a set of class averages");
    //params
    addParamsLine("   -i <md_file>                : Metadata file with input classes");
    addParamsLine("   -o <md_file>                : Metadata file with output information");
    addParamsLine("  [ --vol <file=\"\">]         : Input volume");
    addParamsLine("  [--sym <symfile=c1>]         : Enforce symmetry in projections");
    addParamsLine("  [--numVols <N=5>]            : Number of intermediate volumes to generate");
    addParamsLine("  [--numClasses <N=8>]         : Number of classes to generate the intermediate volumes");
}


// Read arguments ==========================================================
void ProgVolumeValidationPCA::readParams()
{
	fnClasses = getParam("-i");
    fnOut = getParam("-o");
    fnSym = getParam("--sym");
    NVols = getIntParam("--numVols");
    NClasses = getIntParam("--numClasses");
}

// Show ====================================================================
void ProgVolumeValidationPCA::show()
{
    if (verbose > 0)
    {
        std::cout << "Input classes metadata      				 : "  << fnClasses  << std::endl;
        std::cout << "Output metadata             				 : "  << fnOut      << std::endl;
        std::cout << "Number of intermediate volumes to generate : "  << NVols      << std::endl;
        std::cout << "Number of classes to be used               : "  << NClasses   << std::endl;
        if (fnSym != "")
            std::cout << "Symmetry for projections    : "  << fnSym << std::endl;
    }
}

void ProgVolumeValidationPCA::produceSideinfo()
{
	mdClasses.read(fnClasses);
	mdClasses.removeDisabled();
	getImageSize(mdClasses,xdim, ydim, zdim, ndim);

	/*
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
		system(formatString("cp %s %s",fnInit.c_str(),fnVolume.c_str()).c_str());

	threadResults=new ThreadVolumeInitialAlignment[Nthr];
	*/
}

void ProgVolumeValidationPCA::reconstructCurrent()
{
	//number of threads in the reconstruction
	int Nthr = 1;

	String args=formatString("-i %s -o %s --operate random_subset %d --mode overwrite",fnClasses.c_str(),fnAngles.c_str(),NClasses);
	std::cout << args << std::endl;
	String cmd=(String)"xmipp_metadata_utilities "+args;
	system(cmd.c_str());

	args=formatString("-i %s -o %s --sym %s --weight --thr %d -v 0",fnAngles.c_str(),fnVol.c_str(),fnSym.c_str(),Nthr);
	cmd=(String)"xmipp_reconstruct_fourier "+args;
	system(cmd.c_str());

	args=formatString("-i %s --mask circular %d -v 0",fnVol.c_str(),-xdim/2);
	cmd=(String)"xmipp_transform_mask "+args;
	system(cmd.c_str());

	args=formatString("-i %s --select below 0 --substitute value 0 -v 0",fnVol.c_str());
	cmd=(String)"xmipp_transform_threshold "+args;
	system(cmd.c_str());

}

void ProgVolumeValidationPCA::run()
{
    show();
    produceSideinfo();
    int index=0;

    std::stringstream ss;

    for( int index = 0; index< NVols; index++)
    {
        ss << index;
        fnAngles=fnClasses.removeAllExtensions()+ss.str();
        fnAngles+=".xmd";
        ss.str(std::string());

        ss << index;
        fnVol=fnClasses.removeAllExtensions()+ss.str();
        fnVol+=".vol";

        reconstructCurrent();
        std::cout << "File : " << fnVol << std::endl;
        ss.str(std::string());
    }
}

/*
// Main routine ------------------------------------------------------------
void ProgVolumeInitialH::run()
{
    show();
    produceSideinfo();

    boolinish=false;
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

void ProgVolumeInitialH::filterByCorrelation()
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

void ProgVolumeInitialH::reconstructCurrent()
{
	String args=formatString("-i %s -o %s --sym %s --weight --thr %d -v 0",fnAngles.c_str(),fnVolume.c_str(),fnSym.c_str(),Nthr);
	String cmd=(String)"xmipp_reconstruct_fourier "+args;
	system(cmd.c_str());

	args=formatString("-i %s --mask circular %d -v 0",fnVolume.c_str(),-XSIZE(inputImages())/2);
	cmd=(String)"xmipp_transform_mask "+args;
	system(cmd.c_str());

	if (iter<=NiterRandom && positiveConstraint)
	{
		args=formatString("-i %s --select below 0 --substitute value 0 -v 0",fnVolume.c_str());
		cmd=(String)"xmipp_transform_threshold "+args;
		system(cmd.c_str());
	}

	if (keepIntermediateVolumes)
	{
		cmd=formatString("cp %s %s_iter%02d.vol",fnVolume.c_str(),fnRoot.c_str(),iter);
		system(cmd.c_str());
	}
}

void ProgVolumeInitialH::generateProjections()
{
	String args=formatString("-i %s -o %s --sampling_rate 5 --sym %s --compute_neighbors --angular_distance -1 --experimental_images %s -v 0",
			fnVolume.c_str(),fnGallery.c_str(),fnSym.c_str(),fnAngles.c_str());
	String cmd=(String)"xmipp_angular_project_library "+args;
	system(cmd.c_str());
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

void ProgVolumeInitialH::produceSideinfo()
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
		system(formatString("cp %s %s",fnInit.c_str(),fnVolume.c_str()).c_str());

	threadResults=new ThreadVolumeInitialAlignment[Nthr];
}

*/
