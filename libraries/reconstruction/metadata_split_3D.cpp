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
	mdNeighbours.removeDuplicates(mdAux,MDL_IMAGE_IDX);
	std::vector<MDLabel> groupBy;
	groupBy.push_back(MDL_IMAGE_IDX);
	groupBy.push_back(MDL_IMAGE);
	if (mdAux.size()>0)
		mdNeighbours.aggregateGroupBy(mdAux,AGGR_MAX,groupBy,MDL_MAXCC,MDL_MAXCC);
}

void analyzeNeighbours(MetaData &mdNeighbours, const FileName &fnRef)
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
	if (!mdIn.containsLabel(MDL_IMAGE_IDX))
		REPORT_ERROR(ERR_MD_MISSINGLABEL,"Input metadata with images does not contain an imageIndex column");
	mdIn.removeDisabled();
	mdRef.read(fn_oroot+"_gallery.doc");
	deleteFile(fn_oroot+"_gallery.doc");
	deleteFile(fn_oroot+"_gallery.stk");
	deleteFile(fn_oroot+"_gallery_sampling.xmd");

	// Get the maximum reference number
	size_t maxRef=0, refno;
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		mdIn.getValue(MDL_IMAGE_IDX,refno,__iter.objId);
		if (refno>maxRef)
			maxRef=refno;
	}

	// Calculate coocurrence matrix
	correlatesWell.initZeros(maxRef+1);
	Matrix1D<double> projectionDir;
	std::vector<FileName> fnNeighbours;
	FileName fnImg;
	maxDist=DEG2RAD(maxDist);
    MetaData mdNeighbours;
    std::vector<size_t> refs;
    std::vector<double> upperHalf;
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

		// Check if it is upper or lower half
		if (mdNeighbours.size()>0)
		{
			mdRef.getValue(MDL_IMAGE,fnImg,__iter.objId);
			analyzeNeighbours(mdNeighbours,fnImg);
			mdNeighbours.getColumnValues(MDL_IMAGE_IDX,refs);
			mdNeighbours.getColumnValues(MDL_COST,upperHalf);

			size_t imax=refs.size();
			for (size_t i=0; i<imax; ++i)
				VEC_ELEM(correlatesWell,refs[i])+=upperHalf[i];
		}
		++i;
		progress_bar(i);
	}
	progress_bar(mdRef.size());

	// Split in two metadatas
	MetaData mdUpper, mdLower;
	MDRow row;
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		mdIn.getValue(MDL_IMAGE_IDX,refno,__iter.objId);
		mdIn.getRow(row,__iter.objId);
		row.setValue(MDL_COST,(double)VEC_ELEM(correlatesWell,refno));

		if (VEC_ELEM(correlatesWell,refno)>0)
			mdUpper.addRow(row);
		else if (VEC_ELEM(correlatesWell,refno)<0)
			mdLower.addRow(row);
	}
	mdUpper.write(fn_oroot+"_upper.xmd");
	mdLower.write(fn_oroot+"_lower.xmd");
}
