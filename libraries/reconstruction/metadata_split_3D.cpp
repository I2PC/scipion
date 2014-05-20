/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2014)
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

#include "metadata_split_3D.h"
#include <data/geometry.h>
#include <data/mask.h>
#include <data/filters.h>
#include <data/basic_pca.h>

// Read arguments ==========================================================
void ProgMetadataSplit3D::readParams()
{
	fn_in    = getParam("-i");
    fn_vol   = getParam("--vol");
    fn_sym   = getParam("--sym");
    fn_oroot = getParam("--oroot");
    angularSampling = getDoubleParam("--angSampling");
    maxDist = getDoubleParam("--maxDist");
}

// Show ====================================================================
void ProgMetadataSplit3D::show()
{
    if (!verbose)
        return;
    std::cout
    << "Input images:     " << fn_in << std::endl
    << "Reference volume: " << fn_vol << std::endl
    << "Output rootname:  " << fn_oroot << std::endl
    << "Symmetry:         " << fn_sym << std::endl
    << "Angular sampling: " << angularSampling << std::endl
    << "Maximum distance: " << maxDist << std::endl
    ;
}

// usage ===================================================================
void ProgMetadataSplit3D::defineParams()
{
    addUsageLine("Separate projections according to a volume");
    addParamsLine("   -i <metadata>               : Metadata with the images to separate. Make sure they have an angular assignment");
    addParamsLine("   --vol <volume>              : Reference volume");
    addParamsLine("  [--oroot <rootname=split>]   : Rootname for the output files");
    addParamsLine("  [--sym <symmetry_file=c1>]   : Symmetry file if any");
    addParamsLine("                               :+The definition of the symmetry is described at [[transform_symmetrize_v3][transform_symmetrize]]");
    addParamsLine("  [--angSampling <a=5>]        : Angular sampling in degrees");
    addParamsLine("  [--maxDist <a=10>]           : Maximum angular distance in degrees");
    addExampleLine("xmipp_metadata_split_3D -i projections.sel --vol volume.vol --oroot split");
}

void getNeighbours(MetaData &mdIn, const Matrix1D<double> &projectionDir, MetaData &mdNeighbours, double maxDist)
{
	Matrix1D<double> projectionDir2;
	FileName fnImg;
	mdNeighbours.clear();
	MetaData mdAux;
	size_t refno;
	double cc;
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		double rot, tilt;
		mdIn.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
		mdIn.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
		mdIn.getValue(MDL_IMAGE_IDX,refno,__iter.objId);
		mdIn.getValue(MDL_MAXCC,cc,__iter.objId);
		Euler_direction(rot,tilt,0,projectionDir2);

		double angle=acos(dotProduct(projectionDir,projectionDir2));
		if (angle<maxDist)
		{
			mdIn.getValue(MDL_IMAGE,fnImg,__iter.objId);
			size_t id=mdAux.addObject();
			mdAux.setValue(MDL_IMAGE,fnImg,id);
			mdAux.setValue(MDL_IMAGE_IDX,refno,id);
			mdAux.setValue(MDL_MAXCC,cc,id);
		}
	}
	std::cout << "mdAux\n" << mdAux << std::endl;
	mdNeighbours.removeDuplicates(mdAux,MDL_IMAGE_IDX);
	std::cout << "mdNeighbours\n" << mdNeighbours << std::endl;
}

#define DEBUG
void analyzeNeighbours(MetaData &mdNeighbours, const FileName &fnRef, const MultidimArray<int> &imask)
{
    MultidimArray<float> v;
    v.resizeNoCopy((int)imask.sum());

	Image<double> Iref, I;
	Iref.read(fnRef);
	Iref().setXmippOrigin();
#ifdef DEBUG
	Image<double> Iavg;
	Iavg().initZeros(Iref());
	std::cout << "Reference image: " << fnRef << std::endl;
#endif
	PCAMahalanobisAnalyzer analyzer;
	FileName fnImg;
	Matrix2D<double> M;
	AlignmentAux aux;
	CorrelationAux aux2;
	RotationalCorrelationAux aux3;
	FOR_ALL_OBJECTS_IN_METADATA(mdNeighbours)
	{
		mdNeighbours.getValue(MDL_IMAGE,fnImg,__iter.objId);
		I.read(fnImg);
		I().setXmippOrigin();
		alignImagesConsideringMirrors(Iref(),I(),M,aux,aux2,aux3,WRAP,&imask);
		const MultidimArray<double> &mI=I();
		int idx=0;
		FOR_ALL_ELEMENTS_IN_ARRAY2D(imask)
		if (A2D_ELEM(imask,i,j))
		{
			A1D_ELEM(v,idx)=(float)A2D_ELEM(mI,i,j);
			++idx;
		}
		analyzer.addVector(v);

#ifdef DEBUG
		Iavg()+=I();
		std::cout << "   " << fnImg << std::endl;
#endif
	}
#ifdef DEBUG
	Iavg()*=1.0/(double)mdNeighbours.size();
	Iavg.write("PPPavg.xmp");
#endif

	// Construct PCA space
	analyzer.subtractAvg();
	analyzer.learnPCABasis(1, 10);

#ifdef DEBUG
	Iavg().initZeros();
	int idx2=0;
	FOR_ALL_ELEMENTS_IN_ARRAY2D(imask)
	if (A2D_ELEM(imask,i,j))
	{
		Iavg(i,j)=A1D_ELEM(analyzer.PCAbasis[0],idx2);
		++idx2;
	}
	Iavg.write("PPPpca1.xmp");
#endif

	// Threshold PCA
	double mu, sigma;
	MultidimArray<double> &pca=analyzer.PCAbasis[0];
	pca.computeAvgStdev(mu,sigma);
	double sigma2=2*sigma;
	double maskArea=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pca)
	{
		double p=DIRECT_MULTIDIM_ELEM(pca,n)-mu;
		if (p>sigma2)
		{
			DIRECT_MULTIDIM_ELEM(pca,n)=1;
			maskArea+=1;
		}
		else if (p<-sigma2)
		{
			DIRECT_MULTIDIM_ELEM(pca,n)=-1;
			maskArea+=1;
		}
		else
			DIRECT_MULTIDIM_ELEM(pca,n)=0;
	}

#ifdef DEBUG
	Iavg().initZeros();
	idx2=0;
	FOR_ALL_ELEMENTS_IN_ARRAY2D(imask)
	if (A2D_ELEM(imask,i,j))
	{
		Iavg(i,j)=A1D_ELEM(analyzer.PCAbasis[0],idx2);
		++idx2;
	}
	Iavg.write("PPPpca1_thresholded.xmp");
#endif

	// Project onto PCA
    Matrix2D<double> proj;
	analyzer.projectOnPCABasis(proj);
	int idx=0;
	FOR_ALL_OBJECTS_IN_METADATA(mdNeighbours)
	{
		mdNeighbours.setValue(MDL_COST,MAT_ELEM(proj,0,idx)/maskArea,__iter.objId);
		++idx;
	}
#ifdef DEBUG
	mdNeighbours.write("PPPclassification.xmd");
	std::cout << "Press any key\n";
	char c;
	std::cin >> c;
#endif
}
#undef DEBUG

void analyzeNeighbours2(MetaData &mdNeighbours, const FileName &fnRef, const MultidimArray<int> &imask)
{
	std::vector<double> cc;
	mdNeighbours.getColumnValues(MDL_MAXCC,cc);
	std::sort(cc.begin(),cc.end());

	double ccMedian=cc[cc.size()/2];
	FOR_ALL_OBJECTS_IN_METADATA(mdNeighbours)
	{
		double cci;
		mdNeighbours.getValue(MDL_MAXCC,cci,__iter.objId);
		if (cci>ccMedian)
			mdNeighbours.setValue(MDL_COST,1.0,__iter.objId);
		else
			mdNeighbours.setValue(MDL_COST,-1.0,__iter.objId);
	}
}

// usage ===================================================================
void ProgMetadataSplit3D::run()
{
	// Generate projections
	std::cerr << "Generating projections ..." << std::endl;
	String cmd=formatString("xmipp_angular_project_library -i %s -o %s_gallery.stk --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s -v 0",
			fn_vol.c_str(),fn_oroot.c_str(),angularSampling,fn_sym.c_str(),fn_in.c_str());
	if (system(cmd.c_str())!=0)
		REPORT_ERROR(ERR_UNCLASSIFIED,"Error when generating projections");

	// Read reference and input metadatas
	mdIn.read(fn_in);
	mdIn.removeDisabled();
	mdRef.read(fn_oroot+"_gallery.doc");

	// Create mask
	size_t Ndim, Zdim, Ydim, Xdim;
	getImageSize(mdIn,Xdim,Ydim,Zdim,Ndim);
	Mask mask;
	mask.type = BINARY_CIRCULAR_MASK;
	mask.mode = INNER_MASK;
	mask.R1 = Xdim/2;
	mask.generate_mask(Ydim,Xdim);

	// Get the maximum reference number
	size_t maxRef=0, refno;
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		mdIn.getValue(MDL_IMAGE_IDX,refno,__iter.objId);
		if (refno>maxRef)
			maxRef=refno;
	}

	// Calculate coocurrence matrix
	coocurrence.initZeros(maxRef+1,maxRef+1);
	Matrix1D<double> projectionDir;
	std::vector<FileName> fnNeighbours;
	FileName fnImg;
	maxDist=DEG2RAD(maxDist);
    const MultidimArray<int> &imask=mask.get_binary_mask();
    MetaData mdNeighbours;
    std::vector<size_t> refs;
    std::vector<double> pcaProjection;
    int i=0;
    std::cerr << "Classifying projections ...\n";
    init_progress_bar(mdRef.size());
	FOR_ALL_OBJECTS_IN_METADATA(mdRef)
	{
		// Get the projection direction of this image
		double rot, tilt;
		mdRef.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
		mdRef.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
		Euler_direction(rot,tilt,0,projectionDir);

		// Get all images in the input metadata that are close to this one
		getNeighbours(mdIn,projectionDir,mdNeighbours,maxDist);

		// Compute the projection on the first PCA component
		if (mdNeighbours.size()>0)
		{
			mdRef.getValue(MDL_IMAGE,fnImg,__iter.objId);
			//analyzeNeighbours(mdNeighbours,fnImg,imask);
			analyzeNeighbours2(mdNeighbours,fnImg,imask);
			mdNeighbours.getColumnValues(MDL_IMAGE_IDX,refs);
			mdNeighbours.getColumnValues(MDL_COST,pcaProjection);

			size_t imax=refs.size();
			for (size_t i=0; i<imax; ++i)
				for (size_t j=0; j<imax; ++j)
					if ((pcaProjection[i]<0 && pcaProjection[j]<0) || (pcaProjection[i]>0 && pcaProjection[j]>0))
						MAT_ELEM(coocurrence,refs[i],refs[j])+=
								0.5*(fabs(pcaProjection[i])+fabs(pcaProjection[j]));
			std::cout << "Maximum=" << coocurrence.computeMax() << std::endl;
		}
		++i;
		progress_bar(i);
	}
	progress_bar(mdRef.size());
	coocurrence.write("PPPcoocurrence.txt");
}

