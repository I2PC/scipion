/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2012)
 *             Wang Xia                          wangxia@lsec.cc.ac.cn
 *             Guoliang Xu                       xuguo@lsec.cc.ac.cn
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

#include "mpi_classify_FTTRI.h"
#include <data/metadata_extension.h>
#include <data/numerical_tools.h>
#include <data/xmipp_fftw.h>
#include <data/polar.h>
#include <data/histogram.h>
#include <data/filters.h>

//#define FTTRI_ALREADY_COMPUTED

// Empty constructor =======================================================
ProgClassifyFTTRI::ProgClassifyFTTRI(int argc, char **argv)
{
    node=new MpiNode(argc,argv);
    if (!node->isMaster())
        verbose=0;
    taskDistributor=NULL;
}

// MPI destructor
ProgClassifyFTTRI::~ProgClassifyFTTRI()
{
    delete taskDistributor;
    delete node;
}

// Read arguments ==========================================================
void ProgClassifyFTTRI::readParams()
{
    fnIn = getParam("-i");
    fnRoot = getParam("--oroot");
    pad = getDoubleParam("--padding");
    fmax = getDoubleParam("--maxfreq");
    zoom = getDoubleParam("--zoom");
    sigma1 = getDoubleParam("--sigma1");
    sigma2 = getDoubleParam("--sigma2");
    nref = getIntParam("--nref");
    nMinImages = getIntParam("--nmin");
    Niter = getIntParam("--iter");
    doPhase  = checkParam("--doPhase");
}

// Show ====================================================================
void ProgClassifyFTTRI::show()
{
    if (!verbose)
        return;
    std::cout
    << "Input file:  " << fnIn       << std::endl
    << "Output root: " << fnRoot     << std::endl
    << "Padding:     " << pad        << std::endl
    << "MaxFreq:     " << fmax       << std::endl
    << "Zoom:        " << zoom       << std::endl
    << "Sigma1:      " << sigma1     << std::endl
    << "Sigma2:      " << sigma2     << std::endl
    << "N.Classes:   " << nref       << std::endl
    << "N.Minimum:   " << nMinImages << std::endl
    << "Iterations:  " << Niter      << std::endl
    << "Do phase:    " << doPhase    << std::endl
    ;
}

// usage ===================================================================
void ProgClassifyFTTRI::defineParams()
{
    addUsageLine("Classify in 2D using Fourier Transform based Translational and Rotational Invariants");
    addParamsLine("    -i <infile>                : Metadata or stack with input images");
    addParamsLine("    --oroot <rootname>         : Rootname for output files");
    addParamsLine("    --nref <n>                 : Desired number of classes");
    addParamsLine("   [--padding <p=4>]           : Padding factor (it can be a non-integer factor");
    addParamsLine("   [--maxfreq <f=0.25>]        : Maximum frequency for the spectrum classification");
    addParamsLine("                               : Supply -1 for automatic estimation");
    addParamsLine("   [--zoom <f=1>]              : Perform polar transformation with this zoom factor (>=1) at low frequencies");
    addParamsLine("                               : Log-polar corresponds to an approximate factor of 2.8");
    addParamsLine("   [--nmin <n=5>]              : Minimum number of images in a class to be considered as such");
    addParamsLine("   [--iter <n=10>]             : Number of iterations for FRTTI classfication.");
    addParamsLine("                               : At each iteration, those classes whose size");
    addParamsLine("                               : is smaller than nmin are removed from the classification");
    addParamsLine("   [--sigma1+ <s=0.707>]       : First weight in the FTTRI, see paper for documentation");
    addParamsLine("   [--sigma2+ <s=1.5>]         : Second weight in the FTTRI, see paper for documentation");
    addParamsLine("   [--doPhase]                 : Do also an amplitude and phase classification");
    addExampleLine("mpirun -np 4 `which xmipp_mpi_classify_FTTRI` -i images.xmd --oroot class --nref 64");
}

// Produce side info =====================================================
void ProgClassifyFTTRI::produceSideInfo()
{
    mdIn.read(fnIn);

    // Get input size
    size_t zdim, ydim, xdim, ndim;
    getImageSize(mdIn, xdim, ydim, zdim, ndim);
    if (node->isMaster())
    {
        if (zdim!=1)
            REPORT_ERROR(ERR_MULTIDIM_DIM,"This program is only intended for images, not volumes");
        if (xdim!=ydim)
            REPORT_ERROR(ERR_MULTIDIM_DIM,"This program is only intended for squared images");
    }
    padXdim=(int)(pad*xdim);
    Rmax=(size_t)floor(fmax*padXdim);

    size_t blockSize=mdIn.size()/(5*node->size);
    if (blockSize==0)
        blockSize=1;
    mdIn.findObjects(imgsId);
    taskDistributor = new FileTaskDistributor(mdIn.size(), blockSize, node);
    fnFTTRI=fnRoot+"_FTTRI.mrcs";
    FTTRIXdim=(int)((Rmax+1)*0.35);
    FTTRIYdim=(int)((Rmax+1)*0.55);
    if (node->isMaster())
    {
        // Create output mask
        Image<int> mask;
        mask().initZeros(xdim,xdim);
        mask().setXmippOrigin();
        double R2=0.25*xdim*xdim;
        MultidimArray<int> &mMask=mask();
        FOR_ALL_ELEMENTS_IN_ARRAY2D(mMask)
        if (i*i+j*j<R2)
            A2D_ELEM(mMask,i,j)=1;
        mask.write(fnRoot+"_mask.mrc");

        // Create output FTTRI
#ifndef FTTRI_ALREADY_COMPUTED

        if (fileExists(fnFTTRI))
            fnFTTRI.deleteFile();
        createEmptyFile(fnFTTRI, 1+LAST_XMIPP_INDEX(FTTRIXdim), FTTRIYdim, 1, mdIn.size(), true, WRITE_OVERWRITE);
#endif

    }
    node->barrierWait();
}

void ProgClassifyFTTRI::produceFTTRI()
{
    if (verbose && node->rank==0)
        std::cerr << "Computing FTTRI ...\n";

    // Read mask
    Image<int> mask;
    mask.read(fnRoot+"_mask.mrc");
    MultidimArray<int> &mMask=mask();
    mMask.setXmippOrigin();

    // Compute weights
    Matrix1D<double> R, weight1, weight2;

    // Process all images
    size_t Nimgs=mdIn.size();
    if (verbose && node->rank==0)
        init_progress_bar(Nimgs);
    size_t first, last;
    FileName fnImg;
    MultidimArray<double> padI(padXdim,padXdim), magFTI, filteredMagFTI, polarFilteredMagFTI, magFTpolarFilteredMagFTI;
    padI.setXmippOrigin();
    Image<double> I, avgMagFTI, centralMagFTpolarFilteredMagFTI;
    avgMagFTI().initZeros(padI);
    MultidimArray< std::complex<double> > FTI, FTpolarFilteredMagFTI;
    FourierTransformer transformer1, transformer2;
    while (taskDistributor->getTasks(first, last))
        for (size_t idx=first; idx<=last; ++idx)
        {
            if (verbose && node->rank==0)
                progress_bar(idx);
            mdIn.getValue(MDL_IMAGE,fnImg,imgsId[idx]);
            I.read(fnImg);
            I().setXmippOrigin();
            MultidimArray<double> &mI=I();

            // Mask input images
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mMask)
            if (!DIRECT_MULTIDIM_ELEM(mMask,n))
                DIRECT_MULTIDIM_ELEM(mI,n)=0.0;

            // Pad input images
            mI.window(padI,FIRST_XMIPP_INDEX(padXdim),FIRST_XMIPP_INDEX(padXdim),
                      LAST_XMIPP_INDEX(padXdim),LAST_XMIPP_INDEX(padXdim));

            // Compute full Fourier transform and its magnitude
            transformer1.completeFourierTransform(padI,FTI);
            FFT_magnitude(FTI,magFTI);
            CenterFFT(magFTI,true);
            magFTI.window(filteredMagFTI,
                          FIRST_XMIPP_INDEX(Rmax),
                          FIRST_XMIPP_INDEX(Rmax),
                          LAST_XMIPP_INDEX(Rmax),
                          LAST_XMIPP_INDEX(Rmax));

            // Now express it in polar coordinates and apply weight
            image_convertCartesianToPolar_ZoomAtCenter(filteredMagFTI, polarFilteredMagFTI, R, zoom, 3, Rmax/2, Rmax, 0, PI, Rmax);
            if (VEC_XSIZE(weight1)==0)
            {
                weight1.resizeNoCopy(R);
                weight2.resizeNoCopy(R);
                FOR_ALL_ELEMENTS_IN_MATRIX1D(R)
                {
                    VEC_ELEM(weight1,i)=pow(VEC_ELEM(R,i),sigma1);
                    VEC_ELEM(weight2,i)=pow((Rmax-VEC_ELEM(R,i)),sigma2);
                }
            }
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(polarFilteredMagFTI)
            A2D_ELEM(polarFilteredMagFTI,i,j)*=VEC_ELEM(weight1,j);

            // And take new Fourier transform
            transformer2.completeFourierTransform(polarFilteredMagFTI,FTpolarFilteredMagFTI);
            FFT_magnitude(FTpolarFilteredMagFTI,magFTpolarFilteredMagFTI);
            CenterFFT(magFTpolarFilteredMagFTI,true);
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(magFTpolarFilteredMagFTI)
            A2D_ELEM(magFTpolarFilteredMagFTI,i,j)*=VEC_ELEM(weight2,j);
            magFTpolarFilteredMagFTI.setXmippOrigin();
            magFTpolarFilteredMagFTI.window(centralMagFTpolarFilteredMagFTI(),
                                            FIRST_XMIPP_INDEX(FTTRIYdim),0,
                                            LAST_XMIPP_INDEX(FTTRIYdim),LAST_XMIPP_INDEX(FTTRIXdim));
            centralMagFTpolarFilteredMagFTI().rangeAdjust(1,255);
            centralMagFTpolarFilteredMagFTI().selfLog10();
            centralMagFTpolarFilteredMagFTI.write(fnFTTRI,idx+1,true,WRITE_REPLACE);
        }
    if (verbose && node->rank==0)
        progress_bar(Nimgs);
}

void ProgClassifyFTTRI::estimateEpsilonInitialRange()
{
    if (node->isMaster())
    {
        dMin=1e3;
        dMax=-1;
        MultidimArray<int> perm;
        int N=mdIn.size();
        randomPermutation(N,perm);
        N=XMIPP_MIN(N,50);
        std::vector< MultidimArray<double> > fttri;
        Image<double> aux;
        FileName fn;
        for(int i=0; i<N; i++)
        {
            fn.compose(A1D_ELEM(perm,i)+1,fnFTTRI);
            aux.read(fn);
            fttri.push_back(aux());
        }

        int idx=0;
        for(int i=0; i<N-1; i++)
        {
            const MultidimArray<double> &fttri_i=fttri[i];
            for(int j=i+1; j<N; j++, idx++)
            {
                const MultidimArray<double> &fttri_j=fttri[j];
                double d=fttri_distance(fttri_i,fttri_j);
                dMin=std::min(d,dMin);
                dMax=std::max(d,dMax);
            }
        }
        if (verbose)
            std::cout << "Initial epsilon range: [" << dMin << "," << dMax << "]\n";
    }
}

double ProgClassifyFTTRI::fttri_distance(const MultidimArray<double> &fttri_i,
        const MultidimArray<double> &fttri_j)
{
    double retval=0;
    const size_t unroll=4;
    size_t nmax=unroll*(MULTIDIM_SIZE(fttri_i)/unroll);
    double* ptrI=NULL;
    double* ptrJ=NULL;
    size_t n;
    for (n=0, ptrI=MULTIDIM_ARRAY(fttri_i),ptrJ=MULTIDIM_ARRAY(fttri_j);
         n<nmax; n+=unroll, ptrI+=unroll, ptrJ+=unroll)
    {
        double diff0=*(ptrI)-*(ptrJ);
        retval+=diff0*diff0;
        double diff1=*(ptrI+1)-*(ptrJ+1);
        retval+=diff1*diff1;
        double diff2=*(ptrI+1)-*(ptrJ+1);
        retval+=diff2*diff2;
        double diff3=*(ptrI+1)-*(ptrJ+1);
        retval+=diff3*diff3;
    }
    for (n=nmax, ptrI=MULTIDIM_ARRAY(fttri_i)+nmax, ptrJ=MULTIDIM_ARRAY(fttri_j)+nmax;
         n<MULTIDIM_SIZE(fttri_i); ++n, ++ptrI, ++ptrJ)
    {
        double diff0=*(ptrI)-*(ptrJ);
        retval+=diff0*diff0;
    }
    return retval/(MULTIDIM_SIZE(fttri_i));
}

// Epsilon classifcation ==================================================
void ProgClassifyFTTRI::skipRandomNumberOfUnassignedClasses(
    size_t &currentPointer, size_t remaining)
{
    // Now the master selects the first class at random
    if (node->isMaster())
    {
        // Move current pointer to next not assigned image
        while (VEC_ELEM(notAssigned,currentPointer)==0)
            currentPointer=(currentPointer+1)%VEC_XSIZE(notAssigned);

        // Now skip some non empty
        size_t skip=(size_t) (remaining*rnd_unif());
        while (skip>0)
        {
            currentPointer=(currentPointer+1)%VEC_XSIZE(notAssigned);
            // Adjust to next not assigned image
            while (VEC_ELEM(notAssigned,currentPointer)==0)
                currentPointer=(currentPointer+1)%VEC_XSIZE(notAssigned);
            --skip;
        }
    }
    MPI_Bcast(&currentPointer, 1, MPI_LONG, 0, MPI_COMM_WORLD);
}

void receiveListFromRank(std::vector<size_t> &listToKeepIt, int rank)
{
    MPI_Status status;
    Matrix1D<size_t> aux;
    size_t length;
    MPI_Recv(&length, 1, MPI_LONG, rank, 0, MPI_COMM_WORLD, &status);
    aux.resize(length);
    MPI_Recv(&(VEC_ELEM(aux,0)), length, MPI_LONG, rank, 0, MPI_COMM_WORLD, &status);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(aux)
    listToKeepIt.push_back(VEC_ELEM(aux,i));
}

void sendListToRank(std::vector<size_t> &listToSend, int rank)
{
    Matrix1D<size_t> aux;
    aux.resizeNoCopy(listToSend.size());
    FOR_ALL_ELEMENTS_IN_MATRIX1D(aux)
    VEC_ELEM(aux,i)=listToSend[i];
    MPI_Send(&(VEC_XSIZE(aux)), 1, MPI_LONG, rank, 0, MPI_COMM_WORLD);
    MPI_Send(&(VEC_ELEM(aux,0)), VEC_XSIZE(aux), MPI_LONG, rank, 0, MPI_COMM_WORLD);
}

void ProgClassifyFTTRI::epsilonClassification(double epsilon)
{
    epsilonClasses.clear();
    notAssigned=notAssigned0;
    size_t remaining=(size_t)notAssigned.sum();
    size_t currentPointer=0;
    EpsilonClass newClass;
    FileName fnSeed, fnCandidate;
    Image<double> fttriSeed, fftriCandidate;
    std::vector<size_t> inNewClass;

    while (remaining>0)
    {
        // Select new class at random
        skipRandomNumberOfUnassignedClasses(currentPointer,remaining);
        fnSeed.compose(currentPointer+1,fnFTTRI);
        fttriSeed.read(fnSeed);
        MultidimArray<double> &mfftriSeed=fttriSeed();
        if (node->isMaster())
            inNewClass.push_back(currentPointer);
        VEC_ELEM(notAssigned,currentPointer)=0;

        // Check if any of the unassigned images belongs to this class
        FOR_ALL_ELEMENTS_IN_MATRIX1D(notAssigned)
        if (VEC_ELEM(notAssigned,i)==1 && (i+1)%node->size==node->rank)
        {
            fnCandidate.compose(i+1,fnFTTRI);
            fftriCandidate.read(fnCandidate);
            double distance=fttri_distance(mfftriSeed,fftriCandidate());
            if (distance<=epsilon)
                inNewClass.push_back(i);
        }

        // Synchronize lists
        if (node->isMaster())
        {
            // Receive all new elements
            for (size_t rank=1; rank<node->size; rank++)
                receiveListFromRank(inNewClass,rank);
            // And redistribute the new list
            for (size_t rank=1; rank<node->size; rank++)
                sendListToRank(inNewClass,rank);
        }
        else
        {
            // Send my own new elements
            sendListToRank(inNewClass,0);
            // Receive the whole list of new elements
            receiveListFromRank(inNewClass,0);
        }

        // Update internal structures
        size_t imax=inNewClass.size();
        for (size_t i=0; i<imax; ++i)
            VEC_ELEM(notAssigned,inNewClass[i])=0;
        remaining=(size_t)notAssigned.sum();
        if (node->isMaster())
        {
            newClass.memberIdx=inNewClass;
            epsilonClasses.push_back(newClass);
        }
        inNewClass.clear();
    }
}

size_t ProgClassifyFTTRI::wrapperFitness(double epsilon)
{
    MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Workers are waiting at searchOptimalEpsilon
    epsilonClassification(epsilon);
    size_t retval=epsilonClasses.size();
    double objective=fabs(retval-nref);
    if (objective<bestObjective)
    {
        bestObjective=objective;
        bestEpsilon=epsilon;
        bestEpsilonClasses=epsilonClasses;
    }
    std::cout << "Trying " << epsilon << " -> " << retval << std::endl;
    return retval;
}

void ProgClassifyFTTRI::searchOptimalEpsilon()
{
    if (node->isMaster())
    {
        std::cerr << "Computing epsilon classification ..." << std::endl;
        bestEpsilon=-1;
        bestObjective=1e38;
        double epsilonLeft=dMin;
        double epsilonRight=dMax;
        size_t Nleft=wrapperFitness(epsilonLeft);
        size_t Nright=wrapperFitness(epsilonRight);
        if (Nleft!=nref && Nright!=nref)
        {
        	size_t Ndivisions=0;
            do
            {
                // Adjust margins
                if (nref>Nleft)
                {
                    epsilonRight=epsilonLeft;
                    Nright=Nleft;
                    epsilonLeft*=0.9;
                    Nleft=wrapperFitness(epsilonLeft);
                    if (Nleft==0)
                        break;
                }
                else if (nref<Nright)
                {
                    epsilonLeft=epsilonRight;
                    Nleft=Nright;
                    epsilonRight*=1.1;
                    Nright=wrapperFitness(epsilonRight);
                    if (Nright==0)
                        break;
                }
                else
                {
                    // Interval is correct, look in the middle
                    double epsilonCentral=0.5*(epsilonLeft+epsilonRight);
                    size_t Ncentral=wrapperFitness(epsilonCentral);
                    Ndivisions++;
                    if (Ncentral==nref || Ndivisions==10)
                        break;
                    if (Ncentral>nref)
                    {
                        Nleft=Ncentral;
                        epsilonLeft=epsilonCentral;
                    }
                    else
                    {
                        Nright=Ncentral;
                        epsilonRight=epsilonCentral;
                    }
                }
            }
            while (1);
        }
        double finalize=-1e6;
        MPI_Bcast(&finalize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    else
    {
        // Wait for orders from the master
        while (1)
        {
            double epsilon;
            MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            if (epsilon>-5e5)
                epsilonClassification(epsilon);
            else
                break;
        }
    }
}

// Remove small classes
struct SDescendingClusterSort
{
    bool operator()(EpsilonClass const& rpStart, EpsilonClass const& rpEnd)
    {
        return rpStart.memberIdx.size() > rpEnd.memberIdx.size();
    }
};

void ProgClassifyFTTRI::removeSmallClasses()
{
    if (node->isMaster())
    {
        std::sort(bestEpsilonClasses.begin(), bestEpsilonClasses.end(), SDescendingClusterSort());

        std::vector< EpsilonClass > goodEpsilonClasses;
        size_t imax=bestEpsilonClasses.size();
        size_t imagesRemoved=0;
        size_t classesRemoved=0;
        imax=XMIPP_MIN(imax,nref);
        for (size_t i=0; i<imax; i++)
        {
            const std::vector<size_t> class_i=bestEpsilonClasses[i].memberIdx;
            size_t nmax=class_i.size();
            if (nmax>=nMinImages)
                goodEpsilonClasses.push_back(bestEpsilonClasses[i]);
            else
            {
                imagesRemoved+=nmax;
                classesRemoved++;
                for (size_t n=0; n<nmax; ++n)
                    VEC_ELEM(notAssigned0,class_i[n])=0;
            }
        }
        if (imagesRemoved>0)
            std::cout << "Removal of " << imagesRemoved << " images from "
            << classesRemoved << " classes for being in too small classes\n";
        bestEpsilonClasses=goodEpsilonClasses;
    }
    MPI_Bcast(MATRIX1D_ARRAY(notAssigned0), VEC_XSIZE(notAssigned0), MPI_CHAR, 0, MPI_COMM_WORLD);
}

// Split classes ==========================================================
int ProgClassifyFTTRI::findFarthest(const MultidimArray<double> &seed,
                                    const EpsilonClass &class_i, bool FTTRI)
{
    Image<double> candidate;
    FileName fnCandidate;
    int nmax=class_i.memberIdx.size();
    double maxDistance;
    if (FTTRI)
        maxDistance=-1;
    else
        maxDistance=1;
    int nMaxDistance=-1;
    const std::vector<size_t> &class_i_members=class_i.memberIdx;
    double d;
    Matrix2D<double> M;
    AlignmentAux aux;
    CorrelationAux aux2;
    RotationalCorrelationAux aux3;
    for (int n=0; n<nmax; n++)
    {
        if (FTTRI)
            fnCandidate.compose(class_i_members[n]+1,fnFTTRI);
        else
            mdIn.getValue(MDL_IMAGE,fnCandidate,imgsId[class_i_members[n]]);
        candidate.read(fnCandidate);
        candidate().setXmippOrigin();
        if (FTTRI)
            d=fttri_distance(seed,candidate());
        else
            d=alignImages(seed,candidate(),M,WRAP,aux,aux2,aux3);
        if ((d>maxDistance && FTTRI) || (d<maxDistance && !FTTRI))
        {
            maxDistance=d;
            nMaxDistance=n;
        }
    }
    return nMaxDistance;
}

void ProgClassifyFTTRI::splitLargeClasses(bool FTTRI)
{
    if (node->isMaster())
    {
        std::cerr << "Splitting large classes ..." << std::endl;
        Image<double> seed1, seed2, candidate;
        MultidimArray<double> candidateCopy;
        FileName fnSeed, fnCandidate;
        EpsilonClass c1, c2;
        std::sort(bestEpsilonClasses.begin(), bestEpsilonClasses.end(), SDescendingClusterSort());
        Matrix2D<double> M;
        AlignmentAux aux;
        CorrelationAux aux2;
        RotationalCorrelationAux aux3;
        while (bestEpsilonClasses.size()!=nref)
        {
            EpsilonClass &class_0=bestEpsilonClasses[0];
            std::vector<size_t> &class_0_members=class_0.memberIdx;
            if (class_0_members.size()<nMinImages)
                break;

            // Find the image that is farthest from the center
            if (FTTRI)
                fnSeed.compose(class_0_members[0]+1,fnFTTRI);
            else
                mdIn.getValue(MDL_IMAGE,fnSeed,imgsId[class_0_members[0]]);
            seed1.read(fnSeed);
            seed1().setXmippOrigin();
            int n1=findFarthest(seed1(),class_0,FTTRI);
            if (FTTRI)
                fnSeed.compose(class_0_members[n1]+1,fnFTTRI);
            else
                mdIn.getValue(MDL_IMAGE,fnSeed,imgsId[class_0_members[n1]]);
            seed1.read(fnSeed);
            seed1().setXmippOrigin();

            // Now find the one that is farthest from i1
            int n2=findFarthest(seed1(),class_0,FTTRI);
            if (FTTRI)
                fnSeed.compose(class_0_members[n2]+1,fnFTTRI);
            else
                mdIn.getValue(MDL_IMAGE,fnSeed,imgsId[class_0_members[n2]]);
            seed2.read(fnSeed);
            seed2().setXmippOrigin();

            // Now split
            c1.memberIdx.clear();
            c2.memberIdx.clear();
            c1.memberIdx.push_back(class_0_members[n1]);
            c2.memberIdx.push_back(class_0_members[n2]);
            int nmax=class_0_members.size();
            const MultidimArray<double> &mSeed1=seed1();
            const MultidimArray<double> &mSeed2=seed2();
            for (int n=0; n<nmax; ++n)
            {
                if (n==n1 || n==n2)
                    continue;
                size_t trueIdx=class_0_members[n];
                if (FTTRI)
                    fnCandidate.compose(trueIdx+1,fnFTTRI);
                else
                    mdIn.getValue(MDL_IMAGE,fnCandidate,imgsId[trueIdx]);
                candidate.read(fnCandidate);
                if (FTTRI)
                {
                    const MultidimArray<double> &mCandidate=candidate();
                    double d1=fttri_distance(mSeed1,mCandidate);
                    double d2=fttri_distance(mSeed2,mCandidate);
                    if (d1<d2)
                        c1.memberIdx.push_back(trueIdx);
                    else
                        c2.memberIdx.push_back(trueIdx);
                }
                else
                {
                    candidate().setXmippOrigin();
                    candidateCopy=candidate();
                    double d1=alignImages(mSeed1,candidate(),M,WRAP,aux,aux2,aux3);
                    double d2=alignImages(mSeed2,candidateCopy,M,WRAP,aux,aux2,aux3);
                    if (d1>d2)
                        c1.memberIdx.push_back(trueIdx);
                    else
                        c2.memberIdx.push_back(trueIdx);
                }
            }
            bestEpsilonClasses.erase(bestEpsilonClasses.begin());
            if (c1.memberIdx.size()>nMinImages)
                bestEpsilonClasses.push_back(c1);
            if (c2.memberIdx.size()>nMinImages)
                bestEpsilonClasses.push_back(c2);
            std::sort(bestEpsilonClasses.begin(), bestEpsilonClasses.end(), SDescendingClusterSort());
        }
        // Broadcast the best classes
        for (size_t i=0; i<nref; ++i)
        {
            int classSize=bestEpsilonClasses[i].memberIdx.size();
            MPI_Bcast(&classSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&(bestEpsilonClasses[i].memberIdx[0]), classSize, MPI_LONG, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        // Receive the best classes
        bestEpsilonClasses.clear();
        size_t classSize;
        size_t buffer[50000];
        EpsilonClass C;
        for (size_t i=0; i<nref; ++i)
        {
            MPI_Bcast(&classSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(buffer, classSize, MPI_LONG, 0, MPI_COMM_WORLD);
            C.memberIdx.clear();
            for (size_t n=0; n<classSize; ++n)
                C.memberIdx.push_back(buffer[n]);
            bestEpsilonClasses.push_back(C);
        }
    }
}

// Compute centroids =======================================================
void ProgClassifyFTTRI::computeClassCentroids(bool FTTRI)
{
    FileName fnCentroids, fnCandidate;
    if (FTTRI)
        fnCentroids=fnRoot+"_FTTRI_centroids.mrcs";
    else
        fnCentroids=fnRoot+"_image_centroids.mrcs";
    if (node->isMaster())
    {
        std::cerr << "Computing class centroids ..." << std::endl;
        if (FTTRI)
            createEmptyFile(fnCentroids, 1+LAST_XMIPP_INDEX(FTTRIXdim), FTTRIYdim, 1, nref, true, WRITE_OVERWRITE);
        else
        {
            size_t Xdim, Ydim, Zdim, Ndim;
            getImageSize(mdIn, Xdim, Ydim, Zdim, Ndim);
            createEmptyFile(fnCentroids,Xdim,Ydim,1,nref,true,WRITE_OVERWRITE);
        }
    }
    node->barrierWait();

    Image<double> centroid, candidate;
    if (node->isMaster())
        init_progress_bar(nref);
    Image<int> mask;
    if (FTTRI)
    {
        centroid().resizeNoCopy(FTTRIYdim,1+LAST_XMIPP_INDEX(FTTRIXdim));
        classEpsilon.resizeNoCopy(nref);
        classEpsilon.initConstant(-1);
    }
    else
    {
        mask.read(fnRoot+"_mask.mrc");
        centroid().resizeNoCopy(mask());
    }
    MultidimArray<double> &mCentroid=centroid();
    MultidimArray<int> &mMask=mask();
    MultidimArray<double> intraclassDistance, sortedDistance;
    MetaData MDclass;
    for (size_t i=0; i<nref; i++)
        if (((i+1)%node->size)==node->rank)
        {
            const std::vector<size_t> &class_i=bestEpsilonClasses[i].memberIdx;
            size_t nmax=class_i.size();
            mCentroid.initZeros();

            if (nmax!=0)
            {
                if (!FTTRI)
                    MDclass.clear();
                for (size_t n=0; n<nmax; n++)
                {
                    size_t trueIdx=class_i[n];
                    if (FTTRI)
                    {
                        fnCandidate.compose(trueIdx+1,fnFTTRI);
                        candidate.read(fnCandidate);
                        mCentroid+=candidate();
                    }
                    else
                    {
                        mdIn.getValue(MDL_IMAGE,fnCandidate,imgsId[trueIdx]);
                        MDclass.setValue(MDL_IMAGE,fnCandidate,MDclass.addObject());
                    }
                }
                if (FTTRI)
                {
                    mCentroid/=(double)nmax;

                    // Compute now class epsilon
                    intraclassDistance.resizeNoCopy(nmax);
                    for (size_t n=0; n<nmax; n++)
                    {
                        size_t trueIdx=class_i[n];
                        fnCandidate.compose(trueIdx+1,fnFTTRI);
                        candidate.read(fnCandidate);
                        A1D_ELEM(intraclassDistance,n)=fttri_distance(mCentroid,candidate());
                    }
                    intraclassDistance.sort(sortedDistance);
                    int idxLimit=std::min((size_t)floor(nmax*0.8),nmax-1);
                    double limit=A1D_ELEM(sortedDistance,idxLimit);
                    VEC_ELEM(classEpsilon,i)=A1D_ELEM(intraclassDistance,nmax-1);

                    // Robust estimate of the class centroid
                    mCentroid.initZeros();
                    double nactual=0;
                    for (size_t n=0; n<nmax; n++)
                    {
                        if (A1D_ELEM(intraclassDistance,n)>limit)
                            continue;
                        size_t trueIdx=class_i[n];
                        fnCandidate.compose(trueIdx+1,fnFTTRI);
                        candidate.read(fnCandidate);
                        mCentroid+=candidate();
                        nactual++;
                    }
                    mCentroid/=nactual;
                }
                else
                {
                    alignSetOfImages(MDclass,mCentroid,5,false);
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mMask)
                    if (!DIRECT_MULTIDIM_ELEM(mMask,n))
                        DIRECT_MULTIDIM_ELEM(mCentroid,n)=0.0;
                }
            }

            // Write new centroid
            centroid.write(fnCentroids,i+1,true,WRITE_REPLACE);
            if (node->isMaster())
                progress_bar(i);
        }
    if (node->isMaster())
        progress_bar(nref);
    if (FTTRI)
    {
        MPI_Allreduce(MPI_IN_PLACE, &(VEC_ELEM(classEpsilon,0)), nref, MPI_DOUBLE,
                      MPI_MAX, MPI_COMM_WORLD);
        epsilonMax=classEpsilon.computeMax();
        if (node->isMaster())
            std::cerr << "Maximum epsilon: " << epsilonMax << std::endl;
    }
}

// Compute class neighbours ================================================
void ProgClassifyFTTRI::computeClassNeighbours(bool FTTRI)
{
    if (node->isMaster())
        std::cerr << "Computing class neighbours ..." << std::endl;
    MultidimArray<double> *ptrCentroids=NULL;
    if (FTTRI)
    {
        fttriCentroids.read(fnRoot+"_FTTRI_centroids.mrcs");
        ptrCentroids=&fttriCentroids();
    }
    else
    {
        imageCentroids.read(fnRoot+"_image_centroids.mrcs");
        ptrCentroids=&imageCentroids();
    }
    const MultidimArray<double> &mCentroids=*ptrCentroids;
    for (size_t i1=0; i1<nref; ++i1)
        bestEpsilonClasses[i1].neighbours.clear();

    MultidimArray<double> centroid_i1, centroid_i2, centroid_i2p;
    if (FTTRI)
    {
        double limit=2*epsilonMax;
        for (size_t i1=0; i1<nref-1; ++i1)
        {
            centroid_i1.aliasImageInStack(mCentroids,i1);
            for (size_t i2=i1+1; i2<nref; ++i2)
            {
                centroid_i2.aliasImageInStack(mCentroids,i2);
                double d=fttri_distance(centroid_i1,centroid_i2);
                if (d<limit)
                {
                    bestEpsilonClasses[i1].neighbours.push_back(i2);
                    bestEpsilonClasses[i2].neighbours.push_back(i1);
                }
            }
        }
    }
    else
    {
        MultidimArray<double> corr(nref,nref);
        Matrix2D<double> M;
        AlignmentAux aux;
        CorrelationAux aux2;
        RotationalCorrelationAux aux3;
        corr.initConstant(-1);
        int count=0;
        if (node->isMaster())
            init_progress_bar(nref);
        for (size_t i1=0; i1<nref-1; ++i1)
        {
            centroid_i1.aliasImageInStack(mCentroids,i1);
            centroid_i1.setXmippOrigin();
            for (size_t i2=i1+1; i2<nref; ++i2, ++count)
            {
                if (((count+1)%node->size)==node->rank)
                {
                    centroid_i2p.aliasImageInStack(mCentroids,i2);
                    centroid_i2p.setXmippOrigin();
                    centroid_i2=centroid_i2p;
                    A2D_ELEM(corr,i2,i1)=A2D_ELEM(corr,i1,i2)=alignImages(centroid_i1,centroid_i2,M,WRAP,aux,aux2,aux3);
                }
            }
            if (node->isMaster())
                progress_bar(i1);
        }
        MPI_Allreduce(MPI_IN_PLACE, MULTIDIM_ARRAY(corr), nref*nref, MPI_DOUBLE,
                      MPI_MAX, MPI_COMM_WORLD);
        if (node->isMaster())
            progress_bar(nref);

        // For each class compute new neighbours
        MultidimArray<double> corrRow;
        MultidimArray<int> idx;
        const int Kneighbours=2;
        for (size_t i1=0; i1<nref; ++i1)
        {
            corr.getRow(i1,corrRow);
            corrRow.indexSort(idx);
            std::vector<int> &neighbours_i=bestEpsilonClasses[i1].neighbours;
            for (int k = 0; k < Kneighbours; k++)
                neighbours_i.push_back(idx(nref - k - 1) - 1);
        }
    }
}

size_t ProgClassifyFTTRI::reassignImagesToClasses(bool FTTRI)
{
    Matrix1D<short int> newAssignment;
    newAssignment.resizeNoCopy(mdIn.size());
    newAssignment.initConstant(-1);

    if (node->isMaster())
    {
        std::cerr << "Reassigning images to classes ..." << std::endl;
        init_progress_bar(nref);
    }

    MultidimArray<double> *ptrCentroids=NULL;
    if (FTTRI)
        ptrCentroids=&fttriCentroids();
    else
        ptrCentroids=&imageCentroids();
    const MultidimArray<double> &mCentroids=*ptrCentroids;

    Image<double> candidate;
    FileName fnCandidate;
    MultidimArray<double> own_class, neighbour, candidateCopy;
    size_t changes=0;
    AlignmentAux aux;
    CorrelationAux aux2;
    RotationalCorrelationAux aux3;
    Matrix2D<double> M;
    for (size_t i=0; i<nref; i++)
        if (((i+1)%node->size)==node->rank)
        {
            std::vector<size_t> &class_i=bestEpsilonClasses[i].memberIdx;
            std::vector<int> &neighbours_i=bestEpsilonClasses[i].neighbours;
            int nmax=class_i.size();
            int neighMax=neighbours_i.size();
            own_class.aliasImageInStack(mCentroids,i);
            own_class.setXmippOrigin();
            for (int n=0; n<nmax; ++n)
            {
                // Get image n
                int trueIdx=class_i[n];
                if (FTTRI)
                    fnCandidate.compose(trueIdx+1,fnFTTRI);
                else
                    mdIn.getValue(MDL_IMAGE,fnCandidate,imgsId[trueIdx]);
                candidate.read(fnCandidate);
                candidate().setXmippOrigin();

                // Compare to its own class
                const MultidimArray<double> &mCandidate=candidate();
                double bestD;
                if (FTTRI)
                    bestD=fttri_distance(own_class,mCandidate);
                else
                {
                    candidateCopy=candidate();
                    bestD=alignImages(own_class,candidateCopy,M,WRAP,aux,aux2,aux3);
                }
                VEC_ELEM(newAssignment,trueIdx)=i;

                // Now compare to the rest of neighbours
                bool alreadyChanged=false;
                for (int neigh=0; neigh<neighMax; neigh++)
                {
                    int neighbourIdx=neighbours_i[neigh];
                    neighbour.aliasImageInStack(mCentroids,neighbourIdx);
                    neighbour.setXmippOrigin();
                    double d;
                    if (FTTRI)
                        d=fttri_distance(neighbour,mCandidate);
                    else
                    {
                        candidateCopy=candidate();
                        d=alignImages(neighbour,candidateCopy,M,WRAP,aux,aux2,aux3);
                    }
                    if ((d<bestD && FTTRI) || (d>bestD && !FTTRI))
                    {
                        bestD=d;
                        VEC_ELEM(newAssignment,trueIdx)=neighbourIdx;
                        if (!alreadyChanged)
                        {
                            ++changes;
                            alreadyChanged=true;
                        }
                    }
                }
            }
            if (node->isMaster())
                progress_bar(i);
        }

    // Share assignments
    MPI_Allreduce(MPI_IN_PLACE, &(VEC_ELEM(newAssignment,0)), VEC_XSIZE(newAssignment), MPI_SHORT,
                  MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &changes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (node->isMaster())
    {
        progress_bar(nref);
        std::cerr << "Number of images changed: " << changes << std::endl;
    }

    // Create new best classes
    std::vector<EpsilonClass> newBestEpsilonClasses;
    newBestEpsilonClasses.resize(nref);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(newAssignment)
    if (VEC_ELEM(newAssignment,i)>=0)
        newBestEpsilonClasses[VEC_ELEM(newAssignment,i)].memberIdx.push_back(i);
    bestEpsilonClasses=newBestEpsilonClasses;
    return changes;
}

// Write results ==========================================================
void ProgClassifyFTTRI::writeResults(bool FTTRI)
{
    int imax=bestEpsilonClasses.size();
    MetaData MDclass, MDsummary;
    FileName fnImg, fnCandidate;
    FileName fnClasses=fnRoot+"_classes.xmd";
    if (fnClasses.exists())
        fnClasses.deleteFile();

    // Sort each class
    Image<double> candidate;
    MultidimArray<double> centroid;
    MultidimArray<double> distance;
    MultidimArray<int> idx;
    Matrix2D<double> M;
    AlignmentAux aux;
    CorrelationAux aux2;
    RotationalCorrelationAux aux3;
    for (int i=0; i<imax; i++)
    {
        MDclass.clear();
        const std::vector<size_t> &class_i=bestEpsilonClasses[i].memberIdx;
        size_t nmax=class_i.size();

        // Write this class in summary
        size_t id=MDsummary.addObject();
        MDsummary.setValue(MDL_REF, i+1, id);
        MDsummary.setValue(MDL_CLASS_COUNT,nmax,id);
        MDsummary.setValue(MDL_CLASSIFICATION_INTRACLASS_DISTANCE,classEpsilon(i),id);

        // Measure the distance of each class and sort
        if (FTTRI)
        {
            centroid.aliasImageInStack(fttriCentroids(),i);
            distance.resizeNoCopy(nmax);
            for (size_t n=0; n<nmax; n++)
            {
                fnCandidate.compose(class_i[n]+1,fnFTTRI);
                candidate.read(fnCandidate);
                A1D_ELEM(distance,n)=fttri_distance(centroid,candidate());
            }
        }
        else
        {
            centroid.aliasImageInStack(imageCentroids(),i);
            centroid.setXmippOrigin();
            distance.resizeNoCopy(nmax);
            for (size_t n=0; n<nmax; n++)
            {
                int trueIdx=class_i[n];
                mdIn.getValue(MDL_IMAGE,fnCandidate,imgsId[trueIdx]);
                candidate.read(fnCandidate);
                candidate().setXmippOrigin();
                A1D_ELEM(distance,n)=1-alignImages(centroid,candidate(),M,WRAP,aux,aux2,aux3);
            }
        }
        distance.indexSort(idx);

        // Write the class
        MDRow row;
        for (size_t n=0; n<nmax; n++)
        {
            int trueIdx=A1D_ELEM(idx,n)-1;
            size_t trueId=imgsId[class_i[trueIdx]];
            mdIn.getRow(row,trueId);
            size_t idClass=MDclass.addRow(row);
            MDclass.setValue(MDL_COST,A1D_ELEM(distance,trueIdx),idClass);
            MDclass.setValue(MDL_REF,i+1,idClass);
            mdIn.setValue(MDL_REF,i+1,trueId);
        }
        MDclass.write(formatString("class%06d_images@%s",i+1,fnClasses.c_str()),MD_APPEND);
    }
    MDsummary.write(formatString("classes@%s",fnClasses.c_str()),MD_APPEND);
    mdIn.write(fnRoot+"_images.xmd");
}

// Align images within classes ============================================
void ProgClassifyFTTRI::alignImagesWithinClasses()
{
    FileName fnCentroids, fnCandidate;
    fnCentroids=fnRoot+"_image_centroids.mrcs";
    if (node->isMaster())
    {
        std::cerr << "Aligning images within class ..." << std::endl;
        size_t Xdim, Ydim, Zdim, Ndim;
        getImageSize(mdIn, Xdim, Ydim, Zdim, Ndim);
        createEmptyFile(fnCentroids,Xdim,Ydim,1,nref,true,WRITE_OVERWRITE);
    }
    node->barrierWait();

    Image<double> centroid, candidate;
    if (node->isMaster())
        init_progress_bar(nref);
    Image<int> mask;
    mask.read(fnRoot+"_mask.mrc");
    centroid().resizeNoCopy(mask());

    MultidimArray<double> &mCentroid=centroid();
    MultidimArray<int> &mMask=mask();
    MetaData MDclass, MDaux;
    Matrix2D<double> M;
    MDRow row;
    for (size_t i=0; i<nref; i++)
        if (((i+1)%node->size)==node->rank)
        {
            const std::vector<size_t> &class_i=bestEpsilonClasses[i].memberIdx;
            size_t nmax=class_i.size();
            mCentroid.initZeros();

            if (nmax!=0)
            {
                MDclass.clear();
                for (size_t n=0; n<nmax; n++)
                {
                    size_t trueIdx=class_i[n];
                    size_t trueId=imgsId[trueIdx];
                    mdIn.getRow(row,trueId);
                    MDclass.addRow(row);
                }

                // Create centroid
                alignSetOfImages(MDclass,mCentroid,5,false);
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mMask)
                if (!DIRECT_MULTIDIM_ELEM(mMask,n))
                    DIRECT_MULTIDIM_ELEM(mCentroid,n)=0.0;

                // Align images within class to centroid
                FOR_ALL_OBJECTS_IN_METADATA(MDclass)
                {
                	MDclass.getValue(MDL_IMAGE,fnCandidate,__iter.objId);
                    candidate.read(fnCandidate);
                    candidate().setXmippOrigin();
                    double corr=alignImages(mCentroid,candidate(),M);
                    bool flip;
                    double scale, shiftx, shifty, psi;
                    transformationMatrix2Parameters2D(M, flip, scale, shiftx, shifty, psi);
                    MDclass.setValue(MDL_SHIFT_X, shiftx, __iter.objId);
                    MDclass.setValue(MDL_SHIFT_Y, shifty, __iter.objId);
                    MDclass.setValue(MDL_ANGLE_PSI, psi, __iter.objId);
                    MDclass.setValue(MDL_MAXCC, corr, __iter.objId);
                }
            }

            // Write new centroid
            centroid.write(fnCentroids,i+1,true,WRITE_REPLACE);
            MDaux.sort(MDclass,MDL_MAXCC,false);
            MDaux.write(formatString("%s_class%06d.xmd",fnRoot.c_str(),i+1));
            if (node->isMaster())
                progress_bar(i);
        }
    node->barrierWait();
    if (node->isMaster())
    {
        progress_bar(nref);
        FileName fnClasses=fnRoot+"_classes.xmd", fnClass;
        for (size_t i=0; i<nref; i++)
        {
        	fnClass=formatString("%s_class%06d.xmd",fnRoot.c_str(),i+1);
        	MDclass.read(fnClass);
        	MDclass.write(formatString("class%06d_images@%s",i+1,fnClasses.c_str()),MD_APPEND);
        	fnClass.deleteFile();
        }

        MetaData MDsummary;
        FileName classesBlock=(String)"classes@"+fnClasses;
        MDsummary.read(classesBlock);
        int nref=1;
        FileName fnRef;
        FOR_ALL_OBJECTS_IN_METADATA(MDsummary)
        {
        	fnRef.compose(nref++,fnCentroids);
        	MDsummary.setValue(MDL_IMAGE,fnRef,__iter.objId);
        }
        MDsummary.write(classesBlock,MD_APPEND);
    }
}

// Run ====================================================================
void ProgClassifyFTTRI::run()
{
    show();
    produceSideInfo();
#ifndef FTTRI_ALREADY_COMPUTED

    produceFTTRI();
#endif

    estimateEpsilonInitialRange();
    notAssigned0.resizeNoCopy(imgsId.size());
    notAssigned0.initConstant(1);

    searchOptimalEpsilon();
    if (node->isMaster())
        std::cout << "Final epsilon: " << bestEpsilon << "\n";
    removeSmallClasses();

    // Amplitude classification
    if (node->isMaster())
        std::cout << std::endl << "Amplitude classification ..." << std::endl;
    for (int n=0; n<Niter; n++)
    {
        if (node->isMaster())
            std::cout << std::endl << "Iteration " << n << std::endl;
        splitLargeClasses(true);
        computeClassCentroids(true);
        computeClassNeighbours(true);
        reassignImagesToClasses(true);
        removeSmallClasses();
    }

    computeClassCentroids(true);
    if (node->isMaster())
        writeResults(true);

    // Amplitude and phase classification
    if (doPhase)
    {
		if (node->isMaster())
			std::cout << std::endl << "Amplitude and phase classification ..." << std::endl;
		for (int n=0; n<Niter; n++)
		{
			if (node->isMaster())
				std::cout << std::endl << "Iteration " << n << std::endl;
			splitLargeClasses(false);
			computeClassCentroids(false);
			computeClassNeighbours(false);
			reassignImagesToClasses(false);
			removeSmallClasses();
		}
		if (node->isMaster())
			writeResults(false);
    }

    // Align images within classes
    node->barrierWait();
    alignImagesWithinClasses();

}
