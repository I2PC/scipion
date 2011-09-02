/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2011)
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

#include "mpi_classify_CL2D_core_analysis.h"
#include <classification/analyze_cluster.h>
#include <reconstruction/image_sort_by_statistics.h>

// Empty constructor =======================================================
ProgClassifyCL2DCore::ProgClassifyCL2DCore(int argc, char **argv)
{
    node=new MpiNode(argc,argv);
    if (!node->isMaster())
        verbose=0;
    taskDistributor=NULL;
    maxLevel=-1;
}

// MPI destructor
ProgClassifyCL2DCore::~ProgClassifyCL2DCore()
{
    delete taskDistributor;
    delete node;
}

// Read arguments ==========================================================
void ProgClassifyCL2DCore::readParams()
{
    fnRoot = getParam("-i");
    thGoodClass = getDoubleParam("--goodClass");
    thZscore = getDoubleParam("--thZscore");
    thPCAZscore = getDoubleParam("--thPCAZscore");
}

// Show ====================================================================
void ProgClassifyCL2DCore::show()
{
    if (!verbose)
        return;
    std::cout
    << "CL2D rootname:        " << fnRoot << std::endl
    << "Threshold good class: " << thGoodClass << std::endl
    << "Threshold Zscore:     " << thZscore << std::endl
    << "Threshold PCA Zscore: " << thPCAZscore << std::endl
    ;
}

// usage ===================================================================
void ProgClassifyCL2DCore::defineParams()
{
    addUsageLine("Compute the core of a CL2D clustering");
    addParamsLine("    -i <rootname>              : Rootname of the CL2D");
    addParamsLine("  [--goodClass <th=50>]        : A class is a good class if at least 50% (by default) of ");
    addParamsLine("                               : its images have been together in the whole hierarchy of classes.");
    addParamsLine("  [--thZscore <th=3>]          : Threshold for being considered an outlier in a class");
    addParamsLine("  [--thPCAZscore <th=3>]       : Threshold for being considered an outlier in a class after PCA projection");
    addExampleLine("mpirun -np 4 `which xmipp_mpi_classify_CL2D_core_analysis` -i CL2D/results");
}

// Produce side info =====================================================
void ProgClassifyCL2DCore::produceSideInfo()
{
	// Get maximum CL2D level
	maxLevel=0;
	FileName fnLevel, fnLevelCore;
	do
		fnLevel=formatString("%s_level_%02d_classes.xmd",fnRoot.c_str(),maxLevel++);
	while (fnLevel.exists());
	maxLevel-=2;
	if (maxLevel==-1)
		REPORT_ERROR(ERR_ARG_MISSING,"Cannot find any CL2D analysis in the directory given");

	// Read all the blocks available in all MetaData
	StringVector blocksAux;
	CL2DBlock block;
	for (int level=0; level<=maxLevel; level++)
	{
		fnLevel=formatString("%s_level_%02d_classes.xmd",fnRoot.c_str(),maxLevel);
		getBlocksInMetaDataFile(fnLevel,blocksAux);
		block.level=level;
		block.fnLevel=fnLevel;
		fnLevelCore=fnLevel.insertBeforeExtension("_core");
		if (fnLevelCore.exists() && node->rank==0)
			unlink(fnLevelCore.c_str());
		for (int i=0; i<blocksAux.size(); i++)
		{
			if (blocksAux[i].find("class_")!=std::string::npos)
			{
				block.block=blocksAux[i];
				blocks.push_back(block);
			}
		}
	}

	// Create a task file distributor for all blocks
	size_t Nblocks=blocks.size();
    taskDistributor=new FileTaskDistributor(Nblocks,1,node);
}

// Outliers ===============================================================
void ProgClassifyCL2DCore::removeOutliers()
{
	size_t first, last;
	ProgAnalyzeCluster analyzeCluster;
	analyzeCluster.verbose=0;
	analyzeCluster.NPCA=2;
	analyzeCluster.Niter=10;
	analyzeCluster.distThreshold=thPCAZscore;
    analyzeCluster.dontMask=false;

    ProgSortByStatistics sortJunk;
	sortJunk.verbose=0;
    sortJunk.multivariate=true;
    sortJunk.addToInput=true;
    sortJunk.cutoff=thZscore;

    MetaData MD;
    while (taskDistributor->getTasks(first, last))
        for (size_t idx=first; idx<=last; ++idx)
        {
        	// Remove outliers in the PCA projection
        	analyzeCluster.fnSel=blocks[idx].block+"@"+blocks[idx].fnLevel;
        	analyzeCluster.fnOut=blocks[idx].block+"@"+blocks[idx].fnLevel.insertBeforeExtension("_core");
        	analyzeCluster.run();

        	// Remove outliers in the
            sortJunk.fn = analyzeCluster.fnOut;
            sortJunk.run();

            // Remove outliers from file
            MD.read(sortJunk.fn);
            MD.removeDisabled();
            MD.write(sortJunk.fn,MD_APPEND);
        }
}

// Run ====================================================================
void ProgClassifyCL2DCore::run()
{
    show();
    produceSideInfo();
    removeOutliers();
}
