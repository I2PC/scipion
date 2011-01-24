/***************************************************************************
 *
 * Authors:    Carlos Oscar           coss@cnb.csic.es (2010)
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

#include "analyze_cluster.h"
#include <data/args.h>
#include <data/filters.h>
#include <data/mask.h>

// Read arguments ==========================================================
void ProgAnalyzeCluster::readParams()
{
    fnSel = getParam("-i");
    fnOut = getParam("-o");
    fnRef = getParam("--ref");
    //FIXME
    //DROP block ROB
    //if (checkParam("--block"))
    //    block = getParam("--block");
    if (checkParam("--produceAligned"))
        fnOutAligned = getParam("--produceAligned");
    if (checkParam("--basis"))
        fnOutBasis = getParam("--basis");
    quiet = checkParam("--quiet");
    NPCA  = getIntParam("--NPCA");
    Niter  = getIntParam("--iter");
    distThreshold = getDoubleParam("--maxDist");
    dontMask = checkParam("--dontMask");
}

// Show ====================================================================
void ProgAnalyzeCluster::show()
{
    if (!quiet)
        std::cerr
        << "Input metadata file:    " << fnSel         << std::endl
        //<< "Block:                  " << block         << std::endl
        << "Reference:              " << fnRef         << std::endl
        << "Output metadata:        " << fnOut         << std::endl
        << "Output aligned stack:   " << fnOutAligned  << std::endl
        << "Output basis stack:     " << fnOutBasis    << std::endl
        << "PCA dimension:          " << NPCA          << std::endl
        << "Iterations:             " << Niter         << std::endl
        << "Maximum distance:       " << distThreshold << std::endl
        << "Don't mask:             " << dontMask      << std::endl
        ;
}

// usage ===================================================================
void ProgAnalyzeCluster::defineParams()
{
	addUsageLine("Score the images in a cluster according to their PCA projection");
    addParamsLine("   -i <metadatafile>             : metadata file  with images assigned to the cluster");
    addParamsLine("   -o <metadatafile>             : output metadata");
    addParamsLine("   --ref <image>                 : class representative");
    //addParamsLine("  [--block <blockName>]          : block name within the input metadata");
    addParamsLine("  [--produceAligned <stackName>] : write the aligned images");
    addParamsLine("  [--basis <stackName>]          : write the average and basis of the PCA in a stack");
    addParamsLine("  [--NPCA <dim=2>]               : PCA dimension");
    addParamsLine("  [--iter <N=10>]                : Number of iterations");
    addParamsLine("  [--maxDist <d=3>]              : Maximum distance");
    addParamsLine("                                 : Set to -1 if you don't want to filter images");
    addParamsLine("  [--dontMask]                   : Don't use a circular mask");
    addParamsLine("  [--quiet]                      : Don't show anything on screen");
}

// Produce side info  ======================================================
//#define DEBUG
void ProgAnalyzeCluster::produceSideInfo()
{
    bool align = fnOutAligned!="";
    basis = fnOutBasis!="";

    // Read input selfile and reference
    SFin.read(fnSel);
    if (SFin.size()==0)
        return;
    if (SFin.containsLabel(MDL_ENABLED))
    	SFin.removeObjects(MDValueEQ(MDL_ENABLED, -1));

    // Image holding current reference
    Image<double> Iref;
    Iref.read(fnRef);
    Iref().setXmippOrigin();

    // Prepare mask
    mask.resize(Iref());
    mask.setXmippOrigin();
    if (dontMask)
        mask.initConstant(1);
    else
        BinaryCircularMask(mask,XSIZE(Iref())/2, INNER_MASK);
    int Npixels=(int)mask.sum();

    // Read all images in the class and subtract the mean
    // once aligned
    Image<double> Iaux;
    FileName auxFn, fnOutIdx;
    Matrix2D<double> M;
    int idxStk=0;
    pcaAnalyzer.reserve(SFin.size());
    if (align)
    {
        Ialigned.reserve(SFin.size());
        if (exists(fnOutAligned))
        	unlink(fnOutAligned.c_str());
    }
    FOR_ALL_OBJECTS_IN_METADATA(SFin)
    {
        SFin.getValue( MDL_IMAGE, auxFn );
        Iaux.read( auxFn );
        Iaux().setXmippOrigin();

        // Choose between this image and its mirror
        MultidimArray<double> I, Imirror;
        I=Iaux();
        Imirror=I;
        Imirror.selfReverseX();
        Imirror.setXmippOrigin();

#ifdef DEBUG

        Iref.write("PPPref.xmp");
        Iaux->write("PPPclass.xmp");
#endif

        alignImages(Iref(),I,M);
        alignImages(Iref(),Imirror,M);
        double corr=correlation_index(Iref(),I,&mask);
        double corrMirror=correlation_index(Iref(),Imirror,&mask);

        if (corr>corrMirror)
            Iaux()=I;
        else
            Iaux()=Imirror;

        // Produce aligned
        SFout.addObject();
        if (align)
        {
            fnOutIdx.compose(idxStk,fnOutAligned);
            Iaux.write(fnOutAligned,idxStk,true,WRITE_APPEND);
            SFout.setValue(MDL_IMAGE,fnOutIdx);
            SFout.setValue(MDL_IMAGE_ORIGINAL,auxFn);
            idxStk++;
        }
        else
            SFout.setValue(MDL_IMAGE,auxFn);

        MultidimArray<float> v;
        v.initZeros(Npixels);
        int idx=0;
        FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
        if (A2D_ELEM(mask,i,j))
            A1D_ELEM(v,idx++)=IMGPIXEL(Iaux,i,j);
        pcaAnalyzer.addVector(v);
        if (basis)
            Ialigned.push_back(v); // The vector is duplicated because
        // the pcaAnalyzer normalizes the input vectors
        // and then, they cannot be reused

#ifdef DEBUG

        std::cout << "Correlation=" << corr << " mirror=" << corrMirror
        << std::endl;
        Iaux.write("PPPclassAligned.xmp");
        std::cout << "Press any key\n";
        char c;
        std::cin >> c;
#endif

    }
}
#undef DEBUG

// Run  ====================================================================
void ProgAnalyzeCluster::run()
{
    show();
    produceSideInfo();

    pcaAnalyzer.evaluateZScore(NPCA, Niter);

    // Output
    MultidimArray<double> IalignedAvg, Istddev;
    int N=SFin.size();
    if (N>0 && basis)
    {
        IalignedAvg.initZeros(XSIZE(Ialigned[0]));
        Istddev.initZeros(XSIZE(Ialigned[0]));
    }
    double Ngood=0;
    for (int n=0; n<N; n++)
    {
        int trueIdx=pcaAnalyzer.getSorted(n);
        double zscore=pcaAnalyzer.getSortedZscore(n);
        SFout.goToObject(trueIdx+1);
        SFout.setValue(MDL_ZSCORE, zscore);
        if (zscore<distThreshold || distThreshold<0)
        {
            SFout.setValue(MDL_ENABLED,1);
            if (basis)
            {
                const MultidimArray<float> &Ialigned_trueIdx=Ialigned[trueIdx];
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(IalignedAvg)
                {
                	double pixval=DIRECT_MULTIDIM_ELEM(Ialigned_trueIdx,n);
                	DIRECT_MULTIDIM_ELEM(IalignedAvg,n)+=pixval;
                	DIRECT_MULTIDIM_ELEM(Istddev,n)+=pixval*pixval;
                }
                Ngood++;
            }
        }
        else
            SFout.setValue(MDL_ENABLED,-1);
    }
    SFout.write(fnOut);
    if (basis && Ngood>0)
    {
    	if (exists(fnOutBasis))
    		unlink(fnOutBasis.c_str());
    	double iNgood=1.0/Ngood;
    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(IalignedAvg)
    	{
    		DIRECT_MULTIDIM_ELEM(IalignedAvg,n)*=iNgood;
    		DIRECT_MULTIDIM_ELEM(Istddev,n)*=iNgood;
    		DIRECT_MULTIDIM_ELEM(Istddev,n)=sqrt(DIRECT_MULTIDIM_ELEM(Istddev,n)-
    			DIRECT_MULTIDIM_ELEM(IalignedAvg,n)*DIRECT_MULTIDIM_ELEM(IalignedAvg,n));
    	}

        Image<double> save;
        save().initZeros(mask);
        int idx=0;
        FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
        {
            if (A2D_ELEM(mask,i,j))
                IMGPIXEL(save,i,j)=A1D_ELEM(IalignedAvg,idx++);
        }
        save.write(fnOutBasis,0,true,WRITE_APPEND);
        idx=0;
        double avgStd=Istddev.computeAvg();
        FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
        {
            if (A2D_ELEM(mask,i,j))
                IMGPIXEL(save,i,j)=A1D_ELEM(Istddev,idx++);
            else
            	IMGPIXEL(save,i,j)=avgStd;
        }
        save.write(fnOutBasis,1,true,WRITE_APPEND);

        for (int ii=0; ii<pcaAnalyzer.PCAbasis.size(); ii++)
        {
            save().initZeros(mask);
            idx=0;
            FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
            if (A2D_ELEM(mask,i,j))
                IMGPIXEL(save,i,j)=A1D_ELEM(pcaAnalyzer.PCAbasis[ii],idx++);
            save.write(fnOutBasis,ii+2,true,WRITE_APPEND);
        }
    }
}
