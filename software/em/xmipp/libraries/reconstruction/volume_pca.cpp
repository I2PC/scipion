/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2015)
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

#include <data/args.h>
#include <data/filters.h>
#include "volume_pca.h"

// Read arguments ==========================================================
void ProgVolumePCA::readParams()
{
    fnVols = getParam("-i");
    fnVolsOut = getParam("-o");
    NPCA = getIntParam("--Npca");
    fnBasis = getParam("--saveBasis");
    fnAvgVol = getParam("--avgVolume");
    fnOutStack = getParam("--opca");
    if (checkParam("--generatePCAVolumes"))
    	getListParam("--generatePCAVolumes",listOfPercentiles);
    if (checkParam("--mask"))
    	mask.readParams(this);
}

// Show ====================================================================
void ProgVolumePCA::show() const
{
	if (verbose==0)
		return;
    std::cout << "Input volumes:  " << fnVols     << std::endl
    		  << "Output metadata:" << fnVolsOut  << std::endl
    		  << "Number of PCAs: " << NPCA       << std::endl
    		  << "Basis:          " << fnBasis    << std::endl
    		  << "Avg. volume:    " << fnAvgVol   << std::endl
    		  << "Output PCA vols:" << fnOutStack << std::endl;
    std::cout << "Percentiles:    ";
    for (size_t i=0; i<listOfPercentiles.size(); i++)
    	std::cout << listOfPercentiles[i] << " ";
    std::cout << std::endl;
    if (mask.type!=NO_MASK)
    	mask.show();
}

// usage ===================================================================
void ProgVolumePCA::defineParams()
{
	addUsageLine("Compute the PCA of a set of volumes, within a mask (optional).");
    addParamsLine("   -i <volumes>             : Metadata with aligned volumes");
    addParamsLine("  [-o <volumes=\"\">]       : Output metadata with PCA projections");
    addParamsLine("  [--Npca <N=1>]            : Number of PCA bases");
    addParamsLine("  [--saveBasis <stack=\"\">]: Save the bases as a stack of volumes");
    addParamsLine("  [--generatePCAVolumes <...>]: List of percentiles (typically, \"10 90\"), to generate volumes along the 1st PCA basis");
    addParamsLine("  [--avgVolume <volume=\"\">] : Volume on which to add the PCA basis");
    addParamsLine("  [--opca <stack=\"\">]     : Stack of generated volumes");
    mask.defineParams(this,INT_MASK);
}

// Produce side information ================================================
void ProgVolumePCA::produce_side_info()
{
	mdVols.read(fnVols);
	mdVols.removeDisabled();

	size_t Xdim, Ydim, Zdim, Ndim;
	getImageSize(mdVols, Xdim, Ydim, Zdim, Ndim);
	if (mask.type!=NO_MASK)
		mask.generate_mask(Zdim,Ydim,Xdim);
	else
	{
		mask.imask.resizeNoCopy(Zdim,Ydim,Xdim);
		mask.imask.initConstant(1);
	}
}

void ProgVolumePCA::run()
{
    show();
    produce_side_info();

    const MultidimArray<int> &imask=mask.imask;
    size_t Nvoxels=imask.sum();
    MultidimArray<float> v;
    v.initZeros(Nvoxels);

    // Add all volumes to the analyzer
    FileName fnVol;
    FOR_ALL_OBJECTS_IN_METADATA(mdVols)
    {
    	mdVols.getValue(MDL_IMAGE,fnVol,__iter.objId);
    	V.read(fnVol);

    	// Construct vector
    	const MultidimArray<double> &mV=V();
    	size_t idx=0;
    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mV)
    	{
    		if (DIRECT_MULTIDIM_ELEM(imask,n))
    			DIRECT_MULTIDIM_ELEM(v,idx++)=DIRECT_MULTIDIM_ELEM(mV,n);
    	}

    	analyzer.addVector(v);
    }

    // Construct PCA basis
    analyzer.subtractAvg();
    analyzer.learnPCABasis(NPCA,100);

    // Project onto the PCA basis
    Matrix2D<double> proj;
    analyzer.projectOnPCABasis(proj);
    std::vector<double> dimredProj;
    dimredProj.resize(NPCA);
    int i=0;
    FOR_ALL_OBJECTS_IN_METADATA(mdVols)
    {
        memcpy(&dimredProj[0],&MAT_ELEM(proj,i,0),NPCA*sizeof(double));
        mdVols.setValue(MDL_DIMRED,dimredProj,__iter.objId);
        i++;
    }
    if (fnVolsOut!="")
    	mdVols.write(fnVolsOut);
    else
    	mdVols.write(fnVols);

    // Save the basis
	const MultidimArray<double> &mV=V();
	for (int i=NPCA-1; i>=0; --i)
	{
	    V().initZeros();
    	size_t idx=0;
    	const MultidimArray<double> &mPCA=analyzer.PCAbasis[i];
    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mV)
    	{
    		if (DIRECT_MULTIDIM_ELEM(imask,n))
    			DIRECT_MULTIDIM_ELEM(mV,n)=DIRECT_MULTIDIM_ELEM(mPCA,idx++);
    	}
    	if (fnBasis!="")
    		V.write(fnBasis,i+1,true,WRITE_OVERWRITE);
    }

	// Generate the PCA volumes
	if (listOfPercentiles.size()>0 && fnOutStack!="" && fnAvgVol!="")
	{
		Image<double> Vavg;
		if (fnAvgVol!="")
			Vavg.read(fnAvgVol);
		else
			Vavg().initZeros(V());

		Matrix1D<double> p;
		proj.toVector(p);
		Matrix1D<double> psorted=p.sort();

		Image<double> Vpca;
		Vpca()=Vavg();
		createEmptyFile(fnOutStack,(int)XSIZE(Vavg()),(int)YSIZE(Vavg()),(int)ZSIZE(Vavg()),listOfPercentiles.size());
		std::cout << "listOfPercentiles.size()=" << listOfPercentiles.size() << std::endl;
		for (size_t i=0; i<listOfPercentiles.size(); i++)
		{
			int idx=(int)round(textToFloat(listOfPercentiles[i].c_str())/100.0*VEC_XSIZE(p));
			std::cout << "Percentile " << listOfPercentiles[i] << " -> idx=" << idx << " p(idx)=" << psorted(idx) << std::endl;
			Vpca()+=psorted(idx)*V();
			Vpca.write(fnOutStack,i+1,true,WRITE_REPLACE);
		}
	}
}
