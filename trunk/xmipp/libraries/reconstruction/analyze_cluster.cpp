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
}

// Show ====================================================================
void Prog_analyze_cluster_prm::show()
{
    std::cerr << "Input selfile:    " << fnSel         << std::endl
              << "Reference:        " << fnRef         << std::endl
              << "Produce aligned:  " << align         << std::endl
              << "Output extension: " << oext          << std::endl
              << "PCA dimension:    " << NPCA          << std::endl
              << "Iterations:       " << Niter         << std::endl
              << "Maximum distance: " << distThreshold << std::endl
    ;
}

// usage ===================================================================
void Prog_analyze_cluster_prm::usage()
{
    std::cerr << "Usage:  " << std::endl
              << "   -i <selfile>       : selfile with images assigned to the cluster\n"
              << "   -ref <image>       : class representative\n"
              << "  [-produceAligned]   : write the aligned images\n"
              << "  [-oext <ext="">]    : in case you want to produce aligned images\n"
              << "                        use this flag to change the output extension\n"
              << "                        or the input images will be modified\n"
              << "  [-NPCA <dim=2>]     : PCA dimension\n"
              << "  [-iter <N=10>]      : Number of iterations\n"
              << "  [-maxDist <d=2>]    : Maximum distance\n"
    ;
}

// Produce side info  ======================================================
//#define DEBUG
void Prog_analyze_cluster_prm::produceSideInfo()
{
    // Read input selfile and reference
    SelFile SF;
    SF.read(fnSel);
    Iref.read(fnRef);
    Iref().setXmippOrigin();
    
    // Prepare mask
    mask.resize(Iref());
    mask.setXmippOrigin();
    BinaryCircularMask(mask,XSIZE(Iref())/2, INNER_MASK);

    // Read all images in the class and substract the mean
    // once aligned
    Matrix2D<double> Iavg;
    Iavg.initZeros(Iref());
    int currentIdx=-1;
    while (!SF.eof())
    {
        ImageXmipp *Iaux=new ImageXmipp;
        Iaux->read(SF.NextImg());
        (*Iaux)().setXmippOrigin();
        classfile.push_back(Iaux->name());
        currentIdx++;
        
        // Choose between this image and its mirror
        Matrix2D<double> I, Imirror;
        I=(*Iaux)();
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
            (*Iaux)()=I;
        else
            (*Iaux)()=Imirror;

        // Produce aligned
        if (align)
        {
            if (oext=="") Iaux->write();
            else
            {
                FileName fnRoot=Iaux->name().without_extension();
                Iaux->write(fnRoot+"."+oext);
                classfile[currentIdx]=Iaux->name();
            }
        }

        // Mask and store
        FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
            if (!mask(i,j)) (*Iaux)(i,j)=0;
        Iavg+=(*Iaux)();

        ImageXmipp *Iaux2=new ImageXmipp;
        *Iaux2=*Iaux;
        Iclassorig.push_back(Iaux2);
        Iclass.push_back(Iaux);
        
        #ifdef DEBUG
            std::cout << "Correlation=" << corr << " mirror=" << corrMirror
                      << std::endl;
            Iaux->write("PPPclassAligned.xmp");
            std::cout << "Press any key\n";
            char c; std::cin >> c;
        #endif
    }
    Iavg/=Iclass.size();

    // Compute the difference to the mean
    for (int i=0; i<Iclass.size(); i++)
    {
        (*(Iclass[i]))()-=Iavg;
        (*(Iclass[i]))().statisticsAdjust(0,1);
        #ifdef DEBUG
            Iclass[i]->write("PPPdiff.xmp");
            std::cout << "Press any key\n";
            char c; std::cin >> c;
        #endif
    }
    
    // Take the first differences for the PCA basis
    NPCA=XMIPP_MIN(NPCA,Iclass.size());
    for (int i=0; i<NPCA; i++)
    {
        ImageXmipp *IPCA=new ImageXmipp;
        *IPCA=*Iclass[i];
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
        const Matrix2D<double> &Iii=(*(Iclass[ii]))();
        for (int jj=0; jj<NPCA; jj++)
        {
            const Matrix2D<double> &Ijj=(*(PCAbasis[jj]))();
            CtY(jj,ii)=0;
            FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
                if (mask(i,j))
                    DIRECT_MAT_ELEM(CtY,jj,ii)+=
                        MAT_ELEM(Iii,i,j)*MAT_ELEM(Ijj,i,j);
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
            const Matrix2D<double> &Iii=(*(PCAbasis[ii]))();
            for (int jj=ii; jj<NPCA; jj++)
            {
                const Matrix2D<double> &Ijj=(*(PCAbasis[jj]))();
                CtC(ii,jj)=0;
                FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
                    if (mask(i,j))
                        DIRECT_MAT_ELEM(CtC,ii,jj)+=
                            MAT_ELEM(Iii,i,j)*MAT_ELEM(Ijj,i,j);
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
        for (int i=0; i<NPCA; i++)
        {
            Matrix2D<double> &Ipca=(*(PCAbasis[i]))();
            #ifdef DEBUG
                ImageXmipp save;
                save()=Ipca; save.write("PPPoldbasis.xmp");
            #endif
            Ipca.initZeros();
            for (int j=0; j<Nimgs; j++)
            {
                const Matrix2D<double> &I=(*(Iclass[j]))();
                Ipca+=I*XtXXtinv(j,i);
            }
            #ifdef DEBUG
                save()=Ipca; save.write("PPPnewbasis.xmp");
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
    SelFile SFout_good, SFout_bad;
    ImageXmipp Iavg;
    Iavg().initZeros(Iref());
    double N=0;
    std::ofstream fhOut;
    fhOut.open((fnRoot+"_score.txt").c_str());
    if (!fhOut)
        REPORT_ERROR(1,(std::string)"Cannot open "+fnRoot+"_score.txt for output");
    FOR_ALL_ELEMENTS_IN_MATRIX1D(idx) {
        fhOut << classfile[idx(i)-1] << "    " << distance(idx(i)-1) << std::endl;
        if (distance(idx(i)-1)<distThreshold)
        {
            SFout_good.insert(classfile[idx(i)-1]);
            Iavg()+=(*(Iclassorig[idx(i)-1]))();
            N++;
        }
        else
            SFout_bad.insert(classfile[idx(i)-1]);
    }
    fhOut.close();
    if (N>0)
    {
        SFout_good.write(fnRoot+"_pca.sel");
        SFout_bad.write(fnRoot+"_outliers.sel");
        Iavg()/=N;
        Iavg.write(fnRoot+"_pca.xmp");
    }
    
    for (int i=0; i<NPCA; i++)
        PCAbasis[i]->write(fnRoot+"_pcabasis_"+integerToString(i,2)+".xmp");
}
