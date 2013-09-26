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
    addParamsLine("  [--sym <symfile=c1>]              : Enforce symmetry in projections");
}

// Read arguments ==========================================================
void ProgRecRandom::readParams()
{
    fnIn = getParam("-i");
    fnRoot = getParam("--oroot");
    fnSym = getParam("--sym");
}

// Show ====================================================================
void ProgRecRandom::show()
{
    if (verbose > 0)
    {
        std::cout << " Input metadata            : "  << fnIn   << std::endl;
        std::cout << " Output rootname           : "  << fnRoot << std::endl;
        if (fnSym != "")
            std::cout << " Symmetry file for projections : "  << fnSym << std::endl;
    }
}

// Main routine ------------------------------------------------------------
void ProgRecRandom::run()
{
    show();
    produceSideinfo();

    bool finish=false;
    int iter=1;
    do
    {
    	// Generate projections from the volume
    	generateProjections();

    	// Align the input images to the projections
    	size_t nImg=0;
    	double avgCorr=0, avgImprovement=0;
    	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
    	{
    		double corr, improvement;
    		alignSingleImage(nImg++,__iter.objId, corr, improvement);
    		avgCorr+=corr;
    		avgImprovement+=improvement;
    	}
    	avgCorr/=mdIn.size();
    	avgImprovement/=mdIn.size();
    	std::cout << "Iter " << iter << " avg.correlation=" << avgCorr << " avg.improvement=" << avgImprovement << std::endl;

    	// Reconstruct
    	reconstructCurrent();

    	finish=(iter==10);
    	iter++;
    } while (!finish);
}

//#define DEBUG
//#define DEBUG_MORE
void ProgRecRandom::alignSingleImage(size_t nImg, size_t id, double &newCorr, double &improvementFraction)
{
	double bestCorr, angleRot, angleTilt, anglePsi, shiftX, shiftY;
	mdIn.getValue(MDL_MAXCC,bestCorr,id);
	mdIn.getValue(MDL_ANGLE_ROT,angleRot,id);
	mdIn.getValue(MDL_ANGLE_TILT,angleTilt,id);
	mdIn.getValue(MDL_ANGLE_PSI,anglePsi,id);
	mdIn.getValue(MDL_SHIFT_X,shiftX,id);
	mdIn.getValue(MDL_SHIFT_Y,shiftY,id);
	mCurrentImage.aliasImageInStack(inputImages(),nImg);
	mCurrentImage.setXmippOrigin();

	size_t nGallery=0;
	size_t improvementCount=0;
	double oldCorr=bestCorr;
	double bestRandom=-1.0;
	FOR_ALL_OBJECTS_IN_METADATA(mdGallery)
	{
		mCurrentImageAligned=mCurrentImage;
		mGalleryProjection.aliasImageInStack(gallery(),nGallery);
		mGalleryProjection.setXmippOrigin();
		double corr=alignImages(mGalleryProjection,mCurrentImageAligned,M,DONT_WRAP);
#ifdef DEBUG_MORE
	std::cout << "Image " << nImg << " corr=" << corr << std::endl;
	Image<double> save;
	save()=mGalleryProjection;
	save.write("PPPgalleryProjection.xmp");
	save()=mCurrentImageAligned;
	save.write("PPPimageAligned.xmp");
	std::cout << "Press any key\n";
	char c; std::cin >> c;
#endif
		if (corr>oldCorr)
		{
			++improvementCount;
			double random=abs(rnd_gaus(0.0,corr));
			if (random>bestRandom)
			{
				bool flip;
				double scale, psi;
				transformationMatrix2Parameters2D(M,flip,scale,shiftX,shiftY,anglePsi);
				anglePsi*=-1;
				bestRandom=random;
				bestCorr=corr;
				mdGallery.getValue(MDL_ANGLE_ROT,angleRot,__iter.objId);
				mdGallery.getValue(MDL_ANGLE_TILT,angleTilt,__iter.objId);
			}
		}
		++nGallery;
	}

	mdIn.setValue(MDL_MAXCC,bestCorr,id);
	mdIn.setValue(MDL_ANGLE_ROT,angleRot,id);
	mdIn.setValue(MDL_ANGLE_TILT,angleTilt,id);
	mdIn.setValue(MDL_ANGLE_PSI,anglePsi,id);
	mdIn.setValue(MDL_SHIFT_X,shiftX,id);
	mdIn.setValue(MDL_SHIFT_Y,shiftY,id);
	newCorr=bestCorr;
	improvementFraction=(double)improvementCount/NSIZE(gallery());
#ifdef DEBUG
	std::cout << "Image " << nImg << " oldCorr=" << oldCorr << " newCorr=" << bestCorr << " improvement=" << improvementCount << std::endl;
	std::cout << "Press any key\n";
	char c; std::cin >> c;
#endif
}
#undef DEBUG

void ProgRecRandom::reconstructCurrent()
{
	mdIn.write(fnAngles);
	String args=formatString("-i %s -o %s --sym %s -v 0",fnAngles.c_str(),fnVolume.c_str(),fnSym.c_str());
	String cmd=(String)"xmipp_reconstruct_fourier "+args;
	system(cmd.c_str());
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

	fnAngles=fnRoot+".xmd";
	fnVolume=fnRoot+".vol";
	fnGallery=fnRoot+"_gallery.stk";
	fnGalleryMetaData=fnRoot+"_gallery.doc";

	reconstructCurrent();

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
	inputImages.write("PPPinputimages.stk");
}
