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

#include "pdb_restore_with_dictionary.h"

void ProgRestoreWithPDBDictionary::defineParams()
{
    // Usage
    addUsageLine("This program takes a set of PDB files at low and high resolution and constructs a dictionary for them.");

    // Parameters
    addParamsLine(" -i <volume>      : Volume to restore");
    addParamsLine(" [-o <volume=\"\">] : Restored volume");
    addParamsLine(" [--root <root=dictionary>] : Rootname of the dictionary");
	ProgPDBDictionary::defineParams();
}

void ProgRestoreWithPDBDictionary::readParams()
{
	fnIn=getParam("-i");
	fnOut=getParam("-o");
	if (fnOut=="")
		fnOut=fnIn;
	fnRoot=getParam("--root");
	ProgPDBDictionary::readParams();
}

/* Show -------------------------------------------------------------------- */
void ProgRestoreWithPDBDictionary::show()
{
    if (verbose)
    {
		std::cout
		<< "Input volume:        " << fnIn              << std::endl
		<< "Output volume:       " << fnOut             << std::endl
		;
		ProgPDBDictionary::show();
    }
}

//#define DEBUG
void ProgRestoreWithPDBDictionary::run()
{
    show();
    loadDictionaries();
    constructRotationGroup3D();

    Image<double> V, Vhigh;
    MultidimArray<double> weightHigh;
    V.read(fnIn);
	const MultidimArray<double> &mV=V();
	double min, max, mean, std=0;
	mV.computeStats(mean,std,min,max);

    MultidimArray<double> patchLow, patchLowNormalized, canonicalPatch, patchHigh;
    patchSize_2=patchSize/2;
    if (mode==0)
    	patchLow.resizeNoCopy(patchSize,patchSize,patchSize);
    else
    	patchLow.resizeNoCopy(patchSize,patchSize);
    patchLow.setXmippOrigin();
	size_t Npatches=0, NcandidatePatches=0;

	std::vector< size_t > selectedPatchesIdx;
	std::vector<double> weight;
	Matrix1D<double> alpha, canonicalSignature;
	init_progress_bar(ZSIZE(mV));
	Vhigh().initZeros(mV);
	weightHigh.initZeros(mV);
	MultidimArray<double> &mVhigh=Vhigh();
	for (int k=patchSize_2; k<(int)ZSIZE(mV)-patchSize_2; ++k)
	{
		for (int i=patchSize_2; i<(int)ZSIZE(mV)-patchSize_2; ++i)
			for (int j=patchSize_2; j<(int)ZSIZE(mV)-patchSize_2; ++j)
			{
				 ++Npatches;
            	 extractPatch(mV,patchLow,k,i,j);

				double minPatchLow, maxPatchLow, meanPatchLow, stdPatchLow=0;
				patchLow.computeStats(meanPatchLow,stdPatchLow,minPatchLow,maxPatchLow);
				double R2=0;
				patchHigh.initZeros(patchLow);

				if (stdPatchLow > stdThreshold*std)
				{
					++NcandidatePatches;
					patchLowNormalized=patchLow;
					double norm=sqrt(patchLow.sum2());
					patchLowNormalized*=1.0/norm;
					size_t idxTransf=0;
					if (mode==0)
						idxTransf=canonicalOrientation3D(patchLowNormalized,canonicalPatch,canonicalSignature);
					else
						idxTransf=canonicalOrientation2D(patchLowNormalized,canonicalPatch,canonicalSignature);
					selectDictionaryPatches(canonicalPatch, canonicalSignature, selectedPatchesIdx, weight);
					if (selectedPatchesIdx.size()>0)
					{
						R2=approximatePatch(canonicalPatch,selectedPatchesIdx,weight,alpha);
					    reconstructPatch(idxTransf,selectedPatchesIdx,alpha,patchHigh);
					    patchHigh*=norm;
					    R2*=norm;
//						Image<double> save;
//						save()=patchLow;
//						save.write("PPPlow.vol");
//						save()=patchHigh;
//						save.write("PPPhigh.vol");
//						std::cout << "R2=" << R2 << " alpha=" << alpha << std::endl;
//						std::cout << "Press any key" << std::endl;
//						char c; std::cin >> c;
					}
				}
				// Insert patchHigh in Vhigh
				insertPatch(mVhigh,weightHigh,patchHigh,k,i,j,R2);
			}
		progress_bar(k);
	}
	progress_bar(ZSIZE(mV));

	// Correct by the Vhigh weights
//	Image<double> save;
//	save()=weightHigh;
//	save.write("PPPweightHigh.vol");
//	save()=mVhigh;
//	save.write("PPPVHigh.vol");
	FOR_ALL_ELEMENTS_IN_ARRAY3D(mVhigh)
	if (A3D_ELEM(weightHigh,k,i,j)>0)
		A3D_ELEM(mVhigh,k,i,j)/=A3D_ELEM(weightHigh,k,i,j);

	Vhigh.write(fnOut);
}

