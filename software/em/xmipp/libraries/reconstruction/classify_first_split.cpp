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

#include "classify_first_split.h"
#include <data/symmetries.h>

// Read arguments ==========================================================
void ProgClassifyFirstSplit::readParams()
{
    fnClasses = getParam("-i");
    fnRoot = getParam("--oroot");
    Nrec = getIntParam("--Nrec");
    Nsamples = getIntParam("--Nsamples");
    fnSym = getParam("--sym");
    alpha = getDoubleParam("--alpha");
    mask.allowed_data_types = INT_MASK;
    if ((externalMask=checkParam("--mask")))
        mask.readParams(this);
}

// Show ====================================================================
void ProgClassifyFirstSplit::show()
{
    if (!verbose)
        return;
    std::cout
    << "Input classes:       " << fnClasses          << std::endl
    << "Output root:         " << fnRoot             << std::endl
    << "N. reconstructions:  " << Nrec               << std::endl
    << "N. samples:          " << Nsamples           << std::endl
    << "Symmetry:            " << fnSym              << std::endl
	<< "Alpha:               " << alpha              << std::endl
    ;
}

// usage ===================================================================
void ProgClassifyFirstSplit::defineParams()
{
    addUsageLine("Produce a first volume split from a set of directional classes");
    addParamsLine("   -i <metadata>               : Metadata with the list of directional classes with angles");
    addParamsLine("  [--oroot <fnroot=split>]     : Rootname for the output");
    addParamsLine("  [--Nrec <n=100>]             : Number of reconstructions");
    addParamsLine("  [--Nsamples <n=8>]           : Number of images in each reconstruction");
    addParamsLine("  [--sym <sym=c1>]             : Symmetry");
    addParamsLine("  [--alpha <a=0.05>]           : Alpha for the generation of the two separated volumes");
    mask.defineParams(this,INT_MASK);
}

void ProgClassifyFirstSplit::run()
{
    show();

    MetaData md, mdRec;
    md.read(fnClasses);

    // Generate the mean
    std::cerr << "Reconstructing average" << std::endl;
    String command=formatString("xmipp_reconstruct_fourier -i %s -o %s_avg.vol --max_resolution 0.25 -v 0",fnClasses.c_str(),fnRoot.c_str());
    int retval=system(command.c_str());
    Image<double> Vavg;
    Vavg.read(fnRoot+"_avg.vol");
    Vavg().setXmippOrigin();

    SymList SL;
    SL.readSymmetryFile(fnSym);
    int Nsym=SL.symsNo()+1;
    Matrix2D<double> E, L, R;

    FileName fnSubset=fnRoot+"_subset.xmd";
    FileName fnSubsetVol=fnRoot+"_subset.vol";
    command=formatString("xmipp_reconstruct_fourier -i %s -o %s --max_resolution 0.25 -v 0",fnSubset.c_str(),fnSubsetVol.c_str());

    Image<double> V;
    Nvols = 0;
    std::cerr << "Generating reconstructions from random subsets ...\n";
    init_progress_bar(Nrec);
    MultidimArray<double> zn(Nrec);
    pca.maxzn=2; // Skip outliers
    for (int n=0; n<Nrec; n++)
    {
    	// Generate random subset and randomize angles according to symmetry
    	mdRec.selectRandomSubset(md,Nsamples);
    	if (Nsym>1)
    	{
    		double rot, tilt, psi, rotp, tiltp, psip;
    		FOR_ALL_OBJECTS_IN_METADATA(mdRec)
			{
    			int idx=round(rnd_unif(0,Nsym));
    			if (idx>0)
    			{
					mdRec.getValue(MDL_ANGLE_ROT, rot, __iter.objId);
					mdRec.getValue(MDL_ANGLE_TILT, tilt, __iter.objId);
					mdRec.getValue(MDL_ANGLE_PSI, psi, __iter.objId);

    	            SL.getMatrices(idx - 1, L, R, false);
					Euler_apply_transf(L, R, rot, tilt, psi, rotp, tiltp, psip);

					mdRec.setValue(MDL_ANGLE_ROT, rotp, __iter.objId);
					mdRec.setValue(MDL_ANGLE_TILT, tiltp, __iter.objId);
					mdRec.setValue(MDL_ANGLE_PSI, psip, __iter.objId);
    			}
			}
    	}
    	mdRec.write(fnSubset);

    	// Perform reconstruction
    	retval=system(command.c_str());
    	V.read(fnSubsetVol);
    	V().setXmippOrigin();
//    	V.write(formatString("%s_%05d.vol",fnRoot.c_str(),n));
    	V()-=Vavg();
//    	V.write(formatString("%s_%05d_diff.vol",fnRoot.c_str(),n));
    	updateWithNewVolume(V());
    	zn(n)=pca.getCurrentProjection();
//    	std::cout << "zn=" << zn(n) << std::endl;
//    	char c; std::cin >> c;

    	progress_bar(n);
    }
    progress_bar(Nrec);
    deleteFile(fnSubset);
    deleteFile(fnSubsetVol);

    // Save average and first principal component
    Image<double> V1, V2, Vdiff;
    vectorToVolume(pca.c1,V1());
    V1.write(fnRoot+"_pc1.vol");
    vectorToVolume(pca.ysum,V1());
    V1()/=pca.N;
    V1()+=Vavg();
    V2()=V1();

    // Analyze now the projections
    MultidimArray<double> znSorted;
    zn.sort(znSorted);
    double z1=znSorted(int(alpha/2*Nrec));
    double z2=znSorted(int((1-alpha/2)*Nrec));
    std::cout << "z1=" << z1 << " z2=" << z2 << std::endl;

    v=pca.c1;
    v*=z1;
    vectorToVolume(v,Vdiff());
    V1()+=Vdiff();
    V1.write(fnRoot+"_v1.vol");

    v=pca.c1;
    v*=z2;
    vectorToVolume(v,Vdiff());
    V2()+=Vdiff();
    V2.write(fnRoot+"_v2.vol");
}

void ProgClassifyFirstSplit::volumeToVector(const MultidimArray<double> &V, MultidimArray<double> &v)
{
	v.resizeNoCopy(maskSize);
	const MultidimArray<int> &mmask = mask.get_binary_mask();
	size_t idx=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mmask)
	if (DIRECT_MULTIDIM_ELEM(mmask,n))
		DIRECT_MULTIDIM_ELEM(v,idx++)=DIRECT_MULTIDIM_ELEM(V,n);
}

void ProgClassifyFirstSplit::vectorToVolume(const MultidimArray<double> &v, MultidimArray<double> &V)
{
	const MultidimArray<int> &mmask = mask.get_binary_mask();
	V.initZeros(mmask);
	size_t idx=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mmask)
	if (DIRECT_MULTIDIM_ELEM(mmask,n))
		DIRECT_MULTIDIM_ELEM(V,n)=DIRECT_MULTIDIM_ELEM(v,idx++);
}

void ProgClassifyFirstSplit::updateWithNewVolume(const MultidimArray<double> &V)
{
	if (Nvols==0)
	{
	    if (!externalMask)
	    {
	    	mask.type = BINARY_CIRCULAR_MASK;
	    	mask.mode = INNER_MASK;
	    	mask.R1 = XSIZE(V)/2;
	    }
	    mask.generate_mask(V);

	    // Resize some internal variables
	    maskSize=(size_t) mask.get_binary_mask().sum();
	}

	volumeToVector(V,v);
	pca.addVector(v);
	Nvols++;
}

