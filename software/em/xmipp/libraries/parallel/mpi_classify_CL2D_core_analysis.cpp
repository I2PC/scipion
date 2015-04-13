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

// Show block ==============================================================
std::ostream & operator << (std::ostream &out, const CL2DBlock &block)
{
    out << block.block << "@" << block.fnLevel << " -> level " << block.level << std::endl;
    return out;
}

// Empty constructor =======================================================
ProgClassifyCL2DCore::ProgClassifyCL2DCore(int argc, char **argv)
{
    node=new MpiNode(argc,argv);
    if (!node->isMaster())
        verbose=0;
    taskDistributor=NULL;
    maxLevel=-1;
    tolerance=0;
    thPCAZscore=3;
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
    fnRoot = getParam("--root");
    fnODir = getParam("--dir");
    if (checkParam("--computeCore"))
    {
        thPCAZscore = getDoubleParam("--computeCore",0);
        NPCA = getIntParam("--computeCore",1);
        action=COMPUTE_CORE;
    }
    else if (checkParam("--computeStableCore"))
    {
        tolerance = getIntParam("--computeStableCore",0);
        action=COMPUTE_STABLE_CORE;
    }
}

// Show ====================================================================
void ProgClassifyCL2DCore::show()
{
    if (!verbose)
        return;
    std::cout << "CL2D rootname:        " << fnRoot << std::endl
    << "CL2D output dir:      " << fnODir << std::endl;
    if (action==COMPUTE_STABLE_CORE)
        std::cout << "Tolerance:                " << tolerance << std::endl;
    else
        std::cout << "Threshold PCA Zscore:     " << thPCAZscore << std::endl
    			  << "Number of PCA dimensions: " << NPCA << std::endl;
}

// usage ===================================================================
void ProgClassifyCL2DCore::defineParams()
{
    addUsageLine("Compute the core of a CL2D clustering");
    addParamsLine("    --root <rootname>          : Rootname of the CL2D");
    addParamsLine("    --dir <dir>                : Output directory of the CL2D");
    addParamsLine("    --computeCore <thPCAZscore=3> <NPCA=2>: The class cores are computed by thresholding the Zscore of");
    addParamsLine("                               : the class images and their projections onto a nD PCA space (by default, n=2)");
    addParamsLine("or  --computeStableCore <tolerance=1> : The stable core is formed by all the images in the class that have been");
    addParamsLine("                               : in the same class (except in tolerance levels) in the whole hierarchy.");
    addParamsLine("                               : If tolerance=0, then the stable core is formed by the set of images that");
    addParamsLine("                               : have been always together. For this reason, the stable core can be computed only");
    addParamsLine("                               : for the level such that level>tolerance");
    addExampleLine("mpirun -np 4 `which xmipp_mpi_classify_CL2D_core_analysis` -i 2D/CL2D/run_001/results --computeCore");
    addExampleLine("mpirun -np 4 `which xmipp_mpi_classify_CL2D_core_analysis` -i 2D/CL2D/run_001/results --computeStableCore");
}

// Produce side info =====================================================
void ProgClassifyCL2DCore::produceSideInfo()
{
    // Get maximum CL2D level
    maxLevel=0;
    FileName fnLevel;
    do
        fnLevel=formatString("%s/level_%02d/%s_classes.xmd",fnODir.c_str(),maxLevel++,fnRoot.c_str());
    while (fnLevel.exists());
    maxLevel-=2;
    if (maxLevel==-1)
        REPORT_ERROR(ERR_ARG_MISSING,"Cannot find any CL2D analysis in the directory given");

    // Read all the blocks available in all MetaData
    StringVector blocksAux;
    CL2DBlock block;
    for (int level=0; level<=maxLevel; level++)
    {
        fnLevel=formatString("%s/level_%02d/%s_classes.xmd",fnODir.c_str(),level,fnRoot.c_str());
        getBlocksInMetaDataFile(fnLevel,blocksAux);
        block.level=level;
        block.fnLevel=fnLevel;
        block.fnLevelCore=fnLevel.insertBeforeExtension("_core");
        for (size_t i=0; i<blocksAux.size(); i++)
        {
            if (blocksAux[i].find("class")!=std::string::npos && 
                blocksAux[i].find("images")!=std::string::npos)
            {
                block.block=blocksAux[i];
                blocks.push_back(block);
            }
        }
    }

    // Create a task file distributor for all blocks
    size_t Nblocks=blocks.size();
    taskDistributor=new MpiTaskDistributor(Nblocks,1,node);

    // Get image dimensions
    if (Nblocks>0)
    {
        size_t Zdim, Ndim;
        getImageSizeFromFilename(blocks[0].block+"@"+blocks[0].fnLevel,Xdim,Ydim,Zdim,Ndim);
    }
}

// Outliers ===============================================================
void ProgClassifyCL2DCore::computeCores()
{
    if (verbose && node->rank==0)
        std::cerr << "Computing cores ...\n";
    ProgAnalyzeCluster analyzeCluster;
    analyzeCluster.verbose=0;
    analyzeCluster.NPCA=NPCA;
    analyzeCluster.Niter=10;
    analyzeCluster.distThreshold=thPCAZscore;
    analyzeCluster.dontMask=false;

    MetaData MD;
    size_t first, last;
    size_t Nblocks=blocks.size();
    if (verbose && node->rank==0)
        init_progress_bar(Nblocks);
    while (taskDistributor->getTasks(first, last))
        for (size_t idx=first; idx<=last; ++idx)
        {
            // Remove outliers in the PCA projection
            analyzeCluster.SFin.clear();
            analyzeCluster.fnSel=blocks[idx].block+"@"+blocks[idx].fnLevel;
            analyzeCluster.fnOut=blocks[idx].fnLevel.insertBeforeExtension((String)"_core_"+blocks[idx].block);
            analyzeCluster.run();

            // Remove outliers from file
            MD.read(analyzeCluster.fnOut);
            MD.removeDisabled();
            MD.write(analyzeCluster.fnOut,MD_APPEND);

            if (verbose && node->rank==0)
                progress_bar(idx);
        }
    taskDistributor->wait();
    if (verbose && node->rank==0)
        progress_bar(Nblocks);

    // Gather all results
    gatherResults(0,"core");
}

// Outliers ===============================================================
void ProgClassifyCL2DCore::computeStableCores()
{
    if (verbose && node->rank==0)
        std::cerr << "Computing stable cores ...\n";
    MetaData thisClass, anotherClass, commonImages, thisClassCore;
    MDRow row;
    size_t first, last;
    Matrix2D<unsigned char> coocurrence;
    Matrix1D<unsigned char> maximalCoocurrence;
    int Nblocks=blocks.size();
    taskDistributor->reset();
    std::vector<size_t> commonIdx;
    std::map<String,size_t> thisClassOrder;
    String fnImg;
    while (taskDistributor->getTasks(first, last))
        for (size_t idx=first; idx<=last; ++idx)
        {
            // Read block
            CL2DBlock &thisBlock=blocks[idx];
            if (thisBlock.level<=tolerance)
                continue;
            if (!existsBlockInMetaDataFile(thisBlock.fnLevelCore, thisBlock.block))
                continue;
            thisClass.read(thisBlock.block+"@"+thisBlock.fnLevelCore);
            thisClassCore.clear();

            // Add MDL_ORDER
            if (thisClass.size()>0)
            {
                size_t order=0;
                thisClassOrder.clear();
                FOR_ALL_OBJECTS_IN_METADATA(thisClass)
                {
                    thisClass.getValue(MDL_IMAGE,fnImg,__iter.objId);
                    thisClassOrder[fnImg]=order++;
                }

                // Calculate coocurrence within all blocks whose level is inferior to this
                size_t NthisClass=thisClass.size();
                if (NthisClass>0)
                {
                    try {
                       coocurrence.initZeros(NthisClass,NthisClass);
                    } catch (XmippError e)
                    {
                       std::cerr << e << std::endl;
                       std::cerr << "There is a memory allocation error. Most likely there are too many images in this class ("
                                 << NthisClass << " images). Consider increasing the number of initial and final classes\n";
                       REPORT_ERROR(ERR_MEM_NOTENOUGH,"While computing stable class");
                    }
                    for (int n=0; n<Nblocks; n++)
                    {
                        CL2DBlock &anotherBlock=blocks[n];
                        if (anotherBlock.level>=thisBlock.level)
                            break;
                        if (!existsBlockInMetaDataFile(anotherBlock.fnLevelCore, anotherBlock.block))
                            continue;
                        anotherClass.read(anotherBlock.block+"@"+anotherBlock.fnLevelCore);
                        anotherClass.intersection(thisClass,MDL_IMAGE);
                        commonImages.join1(anotherClass, thisClass, MDL_IMAGE,LEFT);
                        commonIdx.resize(commonImages.size());
                        size_t idx=0;
                        FOR_ALL_OBJECTS_IN_METADATA(commonImages)
                        {
                            commonImages.getValue(MDL_IMAGE,fnImg,__iter.objId);
                            commonIdx[idx++]=thisClassOrder[fnImg];
                        }
                        size_t Ncommon=commonIdx.size();
                        for (size_t i=0; i<Ncommon; i++)
                        {
                            size_t idx_i=commonIdx[i];
                            for (size_t j=i+1; j<Ncommon; j++)
                            {
                                size_t idx_j=commonIdx[j];
                                MAT_ELEM(coocurrence,idx_i,idx_j)+=1;
                            }
                        }
                    }
                }

                // Take only those elements whose coocurrence is maximal
                maximalCoocurrence.initZeros(NthisClass);
                int aimedCoocurrence=thisBlock.level-tolerance;
                FOR_ALL_ELEMENTS_IN_MATRIX2D(coocurrence)
                if (MAT_ELEM(coocurrence,i,j)==aimedCoocurrence)
                    VEC_ELEM(maximalCoocurrence,i)=VEC_ELEM(maximalCoocurrence,j)=1;

                // Now compute core
                FOR_ALL_OBJECTS_IN_METADATA(thisClass)
                {
                    thisClass.getValue(MDL_IMAGE,fnImg,__iter.objId);
                    size_t idx=thisClassOrder[fnImg];
                    if (VEC_ELEM(maximalCoocurrence,idx))
                    {
                        thisClass.getRow(row,__iter.objId);
                        thisClassCore.addRow(row);
                    }
                }
            }
            thisClassCore.write(thisBlock.fnLevel.insertBeforeExtension((String)"_stable_core_"+thisBlock.block),MD_APPEND);
        }
    taskDistributor->wait();

    // Gather all results
    gatherResults(tolerance+1,"stable_core");
}

void ProgClassifyCL2DCore::gatherResults(int firstLevel, const String &suffix)
{
    node->barrierWait();
    if (node->rank==0)
    {
        FileName fnBlock, fnClass, fnSummary;
        Image<double> classAverage;
        // Compute class averages
        MetaData classes, MD, MDoriginal;
        int Nblocks=blocks.size();
        for (int level=firstLevel; level<=maxLevel; level++)
        {
            classes.clear();
            fnSummary=formatString("%s/level_%02d/%s_classes_%s",fnODir.c_str(),level,fnRoot.c_str(),suffix.c_str());
            for (int idx=0; idx<Nblocks; idx++)
            {
                if (blocks[idx].level!=level)
                    continue;
                fnBlock=fnSummary+"_"+blocks[idx].block+".xmd";
                if (fileExists(fnBlock))
                {
                    MD.read(fnBlock);
                    MDoriginal.read(formatString("%s@%s/level_%02d/%s_classes.xmd",blocks[idx].block.c_str(),fnODir.c_str(),level,fnRoot.c_str()));
                    int classNo=textToInteger(blocks[idx].block.substr(6,6));
                    size_t classSize=MD.size();
                    fnClass.compose(classNo,fnSummary,"stk");
                    if (classSize>0)
                        getAverageApplyGeo(MD, classAverage());
                    else
                        classAverage().initZeros(Ydim,Xdim);
                    classAverage.write(fnClass);

                    size_t id=classes.addObject();
                    classes.setValue(MDL_REF,classNo,id);
                    classes.setValue(MDL_IMAGE,fnClass,id);
                    classes.setValue(MDL_CLASS_COUNT,classSize,id);
                    classes.setValue(MDL_MODELFRAC,((double)classSize)/MDoriginal.size(),id);
                }
            }
            classes.write((String)"classes@"+fnSummary+".xmd");
        }

        // Write the rest of blocks
        for (int idx=0; idx<Nblocks; idx++)
        {
            fnSummary=formatString("%s/level_%02d/%s_classes_%s",fnODir.c_str(),blocks[idx].level,fnRoot.c_str(),suffix.c_str());
            fnBlock=fnSummary+"_"+blocks[idx].block+".xmd";
            if (fileExists(fnBlock))
            {
                MD.read(fnBlock);
                unlink(fnBlock.c_str());
                MD.write(blocks[idx].block+"@"+fnSummary+".xmd",MD_APPEND);
            }
        }
    }
    node->barrierWait();
}

// Run ====================================================================
void ProgClassifyCL2DCore::run()
{
    show();
    produceSideInfo();
    if (action==COMPUTE_CORE)
        computeCores();
    else
        computeStableCores();
}
