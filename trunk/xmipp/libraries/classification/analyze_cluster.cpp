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
void Prog_analyze_cluster_prm::read(int argc, char **argv)
{
    fnSel = getParameter(argc, argv, "-i");
    fnRef = getParameter(argc, argv, "-ref");
    oext  = getParameter(argc, argv, "-oext", "");
    align = checkParameter(argc, argv, "-produceAligned");
    NPCA  = textToInteger(getParameter(argc, argv, "-NPCA", "2"));
    Niter  = textToInteger(getParameter(argc, argv, "-iter", "10"));
    distThreshold = textToFloat(getParameter(argc, argv, "-maxDist", "2"));
    dontMask = checkParameter(argc, argv, "-dontMask");
}

// Show ====================================================================
void Prog_analyze_cluster_prm::show()
{
    std::cerr << "Input metadata file:    " << fnSel         << std::endl
    << "Reference:              " << fnRef         << std::endl
    << "Produce aligned:        " << align         << std::endl
    << "Output extension:       " << oext          << std::endl
    << "PCA dimension:          " << NPCA          << std::endl
    << "Iterations:             " << Niter         << std::endl
    << "Maximum distance:       " << distThreshold << std::endl
    << "Don't mask:       " << dontMask      << std::endl
    ;
}

// usage ===================================================================
void Prog_analyze_cluster_prm::usage()
{
    std::cerr << "Usage:  " << std::endl
    << "   -i <metadatafile>  : metadata file  with images assigned to the cluster\n"
    << "   -ref <image>       : class representative\n"
    << "  [-produceAligned]   : write the aligned images\n"
    << "  [-oext <ext=''>]    : in case you want to produce aligned images\n"
    << "                        use this flag to change the output extension\n"
    << "                        or the input images will be modified\n"
    << "  [-NPCA <dim=2>]     : PCA dimension\n"
    << "  [-iter <N=10>]      : Number of iterations\n"
    << "  [-maxDist <d=2>]    : Maximum distance\n"
    << "  [-dontMask]         : Don't use a circular mask\n"
    ;
}

// Produce side info  ======================================================
//#define DEBUG
void Prog_analyze_cluster_prm::produceSideInfo()
{
    // Read input selfile and reference
    SFin.read(fnSel, NULL);
    SFin.removeObjects(MDValueEqual(MDL_ENABLED, -1));

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

    // Read all images in the class and substract the mean
    // once aligned
    FOR_ALL_OBJECTS_IN_METADATA(SFin)
    {
        Image<double> Iaux;
        FileName auxFn;
        SFin.getValue( MDL_IMAGE, auxFn );
        Iaux.read( auxFn );
        Iaux().setXmippOrigin();
        SFout.addObject();
        SFout.setValue(MDL_IMAGE,auxFn);

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

        Matrix2D<double> M;
        alignImages(Iref(),I,M);
        alignImages(Iref(),Imirror,M);
        double corr=correlation_index(Iref(),I,&mask);
        double corrMirror=correlation_index(Iref(),Imirror,&mask);

        if (corr>corrMirror)
            Iaux()=I;
        else
            Iaux()=Imirror;

        // Produce aligned
        if (align)
        {
            if (oext=="")
                Iaux.write();
            else
            {
                FileName fnRoot=auxFn.without_extension();
                Iaux.write(fnRoot+"."+oext);
                SFout.setValue(MDL_IMAGE,fnRoot+"."+oext);
            }
        }

        MultidimArray<float> v;
        v.initZeros(Npixels);
        int idx=0;
        FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
        if (mask(i,j))
            v(idx++)=Iaux(i,j);
        pcaAnalyzer.addVector(v);
        Ialigned.push_back(v);

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
void Prog_analyze_cluster_prm::run()
{
	pcaAnalyzer.evaluateZScore(NPCA, Niter);

    // Output
    FileName fnRoot=fnSel.without_extension();
    MetaData SFout_good, SFout_bad;
    MultidimArray<double> IalignedAvg;
    IalignedAvg.initZeros(XSIZE(Ialigned[0]));
    double Ngood=0;
    int N=SFin.size();
    for (int n=0; n<N; n++)
    {
    	int trueIdx=pcaAnalyzer.getSorted(n);
    	double zscore=pcaAnalyzer.getSortedZscore(n);
    	FileName fnImg;
    	SFout.getValue(MDL_IMAGE,fnImg,trueIdx+1);
        if (zscore<distThreshold)
        {
            SFout_good.addObject();
            SFout_good.setValue( MDL_IMAGE, fnImg);
            SFout_good.setValue( MDL_ZSCORE, zscore);
            MultidimArray<float>& Iaux=Ialigned[trueIdx];
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(IalignedAvg)
				DIRECT_A1D_ELEM(IalignedAvg,i)+=DIRECT_A1D_ELEM(Iaux,i);
            Ngood++;
        }
        else
        {
            SFout_bad.addObject();
            SFout_bad.setValue( MDL_IMAGE, fnImg);
            SFout_bad.setValue( MDL_ZSCORE, zscore);
        }
    }
    SFout_good.write(fnRoot+"_pca.sel");
    SFout_bad.write(fnRoot+"_outliers.sel");
    if (Ngood>0)
    {
        IalignedAvg/=Ngood;
        int idx=0;
        Image<double> save;
        save().initZeros(mask);
        FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
        if (mask(i,j))
            save(i,j)=IalignedAvg(idx++);
            idx++;
        save.write(fnRoot+"_pca.xmp");
    }

    for (int ii=0; ii<NPCA; ii++)
    {
        Image<double> save;
        save().initZeros(mask);
        int idx=0;
        FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
        if (mask(i,j))
            save(i,j)=pcaAnalyzer.PCAbasis[ii](idx++);
        save.write(fnRoot+"_pcabasis_"+integerToString(ii,2)+".xmp");
    }
}
