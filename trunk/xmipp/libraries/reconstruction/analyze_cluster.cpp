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
              << "  [-oext <ext="">]    : in case you want to produce aligned images\n"
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
    MetaData SF;
    SF.read(fnSel,NULL);
    SF.removeObjects( MDL_ENABLED, -1 );
    //remove descarted images
    Iref.read(fnRef);
    Iref().setXmippOrigin();
    
    // Prepare mask
    mask.resize(Iref());
    mask.setXmippOrigin();
    if (dontMask)
        mask.initConstant(1);
    else
        BinaryCircularMask(mask,XSIZE(Iref())/2, INNER_MASK);
    Npixels=(int)mask.sum();

    // Read all images in the class and substract the mean
    // once aligned
    Matrix2D<double> Iavg;
    Iavg.initZeros(Iref());
    int currentIdx=-1;
    do
    {
        ImageXmipp Iaux;
        FileName auxFn;
        SF.getValue( MDL_IMAGE, auxFn );
        Iaux.read( auxFn );
        Iaux().setXmippOrigin();
        classfile.push_back(Iaux.name());
        currentIdx++;
        
        // Choose between this image and its mirror
        Matrix2D<double> I, Imirror;
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
            if (oext=="") Iaux.write();
            else
            {
                FileName fnRoot=Iaux.name().without_extension();
                Iaux.write(fnRoot+"."+oext);
                classfile[currentIdx]=Iaux.name();
            }
        }

        // Mask and store
        FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
            if (!mask(i,j)) Iaux(i,j)=0;
        Iavg+=Iaux();

        Matrix1D<float> *v=new Matrix1D<float>(Npixels);
        Matrix1D<float> *v2=new Matrix1D<float>(Npixels);
        int idx=0;
        FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
            if (mask(i,j))
            {
                (*v2)(idx)=(*v)(idx)=I(i,j);
                idx++;
            }
        Iclassorig.push_back(v2);
        Iclass.push_back(v);
        
        #ifdef DEBUG
            std::cout << "Correlation=" << corr << " mirror=" << corrMirror
                      << std::endl;
            Iaux.write("PPPclassAligned.xmp");
            std::cout << "Press any key\n";
            char c; std::cin >> c;
        #endif
    }
    while (SF.nextObject()!= MetaData::NO_MORE_OBJECTS);
    Iavg/=Iclass.size();

    // Compute the difference to the mean
    for (int ii=0; ii<Iclass.size(); ii++)
    {
        int idx=0;
        FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
            if (mask(i,j))
            {
                (*(Iclass[ii]))(idx)-=Iavg(i,j);
                idx++;
            }
        (*(Iclass[ii])).statisticsAdjust(0,1);
        #ifdef DEBUG
            ImageXmipp save;
            save().resize(mask);
            idx=0;
            FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
                if (mask(i,j))
                {
                    save(i,j)=(*(Iclass[ii]))(idx);
                    idx++;
                }
            save.write("PPPdiff.xmp");
            std::cout << "Press any key\n";
            char c; std::cin >> c;
        #endif
    }
    
    // Take the first differences for the PCA basis
    NPCA=XMIPP_MIN(NPCA,Iclass.size());
    for (int ii=0; ii<NPCA; ii++)
    {
        Matrix1D<double> *IPCA=new Matrix1D<double>(Npixels);
        typeCast(*Iclass[ii],*IPCA);
        PCAbasis.push_back(IPCA);
    }
}
#undef DEBUG

// projectOnPCABasis  ======================================================
void Prog_analyze_cluster_prm::projectOnPCABasis(Matrix2D<double> &CtY)
{
    int Nimgs=Iclass.size();
    CtY.initZeros(NPCA,Nimgs);
    for (int ii=0; ii<Nimgs; ii++)
    {
        const Matrix1D<float> &Iii=*(Iclass[ii]);
        for (int jj=0; jj<NPCA; jj++)
        {
            const Matrix1D<double> &Ijj=*(PCAbasis[jj]);
            CtY(jj,ii)=0;
            FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(Iii)
                DIRECT_MAT_ELEM(CtY,jj,ii)+=
                    DIRECT_VEC_ELEM(Iii,i)*DIRECT_VEC_ELEM(Ijj,i);
        }
    }
}

// learnPCABasis  ==========================================================
/* See Roweis, "EM algorithms for PCA and SPCA",
    Neural Information Processing Systems 10 (NIPS'97) pp.626-632 */
//#define DEBUG
void Prog_analyze_cluster_prm::learnPCABasis()
{
    int Nimgs=Iclass.size();
    for (int n=0; n<Niter; n++)
    {
        // E-step ..........................................................
        // Compute C^t*C
        Matrix2D<double> CtC(NPCA,NPCA);
        for (int ii=0; ii<NPCA; ii++)
        {
            const Matrix1D<double> &Iii=*(PCAbasis[ii]);
            for (int jj=ii; jj<NPCA; jj++)
            {
                const Matrix1D<double> &Ijj=*(PCAbasis[jj]);
                CtC(ii,jj)=0;
                FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(Ijj)
                    DIRECT_MAT_ELEM(CtC,ii,jj)+=
                        DIRECT_VEC_ELEM(Iii,i)*DIRECT_VEC_ELEM(Ijj,i);
                if (ii!=jj)
                    CtC(jj,ii)=CtC(ii,jj);
            }
        }
        
        // Compute C^t*Y
        Matrix2D<double> CtY, X;
        projectOnPCABasis(CtY);
        X=CtC.inv()*CtY;
        
        // M-step ..........................................................
        Matrix2D<double> XtXXtinv=X.transpose()*(X*X.transpose()).inv();
        for (int ii=0; ii<NPCA; ii++)
        {
            Matrix1D<double> &Ipca=*(PCAbasis[ii]);
            #ifdef DEBUG
                ImageXmipp save;
                save().resize(mask);
                int idx=0;
                FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
                    if (mask(i,j))
                    {
                        save(i,j)=Ipca(idx);
                        idx++;
                    }
                save.write("PPPoldbasis.xmp");
            #endif
            Ipca.initZeros();
            for (int jj=0; jj<Nimgs; jj++)
            {
                const Matrix1D<float> &I=*(Iclass[jj]);
                double val=XtXXtinv(jj,ii);
                FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX1D(I)
                    DIRECT_VEC_ELEM(Ipca,i)+=DIRECT_VEC_ELEM(I,i)*val;
            }
            #ifdef DEBUG
                idx=0;
                FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
                    if (mask(i,j))
                    {
                        save(i,j)=Ipca(idx);
                        idx++;
                    }
                save.write("PPPnewbasis.xmp");
                std::cout << "New basis: " << i << ". Press any key" << std::endl;
                char c; std::cin >> c;
            #endif
        }
    }
}
#undef DEBUG

// projectOnPCABasis  ======================================================
void Prog_analyze_cluster_prm::evaluateDistance(Matrix2D<double> &proj)
{
    int Nimgs=Iclass.size();
    
    // Estimate covariance matrix
    Matrix2D<double> cov(NPCA,NPCA);
    for (int ii=0; ii<Nimgs; ii++)
    {
        for (int i=0; i<NPCA; i++)
            for (int j=0; j<NPCA; j++)
                cov(i,j)+=proj(i,ii)*proj(j,ii);
    }
    cov/=Nimgs;
    
    Matrix2D<double> covinv=cov.inv();
    distance.initZeros(Nimgs);
    for (int ii=0; ii<Nimgs; ii++)
    {
        Matrix1D<double> x,aux;
        x.resize(NPCA);
        for (int i=0; i<NPCA; i++)
            x(i)=proj(i,ii);
        aux=x.transpose()*covinv*x;
        distance(ii)=ABS(aux(0));
    }
}

// Run  ====================================================================
void Prog_analyze_cluster_prm::run()
{
    learnPCABasis();
    Matrix2D<double> proj;
    projectOnPCABasis(proj);
    evaluateDistance(proj);

    // Output    
    FileName fnRoot=fnSel.without_extension();
    Matrix1D<int> idx=distance.indexSort();
    MetaData SFout_good, SFout_bad;
    Matrix1D<float> Iavg(Npixels);
    double N=0;
    std::ofstream fhOut;
    fhOut.open((fnRoot+"_score.txt").c_str());
    if (!fhOut)
        REPORT_ERROR(1,(std::string)"Cannot open "+fnRoot+"_score.txt for output");
    FOR_ALL_ELEMENTS_IN_MATRIX1D(idx) {
        fhOut << classfile[idx(i)-1] << "    " << distance(idx(i)-1) << std::endl;
        if (distance(idx(i)-1)<distThreshold)
        {
        	SFout_good.addObject();
        	SFout_good.setValue( MDL_IMAGE, classfile[idx(i)-1]);
            SFout_good.setValue( MDL_ENABLED, 1);
            Iavg+=*(Iclassorig[idx(i)-1]);
            N++;
        }
        else
        {
        	SFout_bad.addObject();
        	SFout_bad.setValue( MDL_IMAGE, classfile[idx(i)-1]);
        	SFout_bad.setValue( MDL_ENABLED, 1);
        }
    }
    fhOut.close();
    if (N>0)
    {
        SFout_good.write(fnRoot+"_pca.sel");
        SFout_bad.write(fnRoot+"_outliers.sel");
        Iavg/=N;
        int idx=0;
        ImageXmipp save;
        save().initZeros(mask);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
            if (mask(i,j))
            {
                save(i,j)=Iavg(idx);
                idx++;
            }
        save.write(fnRoot+"_pca.xmp");
    }
    
    for (int ii=0; ii<NPCA; ii++)
    {
        ImageXmipp save;
        save().initZeros(mask);
        int idx=0;
        FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
            if (mask(i,j))
            {
                save(i,j)=(*(PCAbasis[ii]))(idx);
                idx++;
            }
        save.write(fnRoot+"_pcabasis_"+integerToString(ii,2)+".xmp");
    }
}
