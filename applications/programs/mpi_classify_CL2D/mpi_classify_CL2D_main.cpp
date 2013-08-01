/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2009)
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

#include "mpi_classify_CL2D.h"
#include <data/filters.h>
#include <data/mask.h>
#include <data/polar.h>
#include <data/xmipp_image_generic.h>

// Pointer to parameters
ProgClassifyCL2D *prm = NULL;
FILE * _logCL2D = NULL;

//#define DEBUG_WITH_LOG
#ifdef DEBUG_WITH_LOG
#define CREATE_LOG() _logCL2D = fopen(formatString("nodo%02d.log", node->rank).c_str(), "w+")
#define LOG(msg) do{fprintf(_logCL2D, "%s\t%s\n", getCurrentTimeString(), msg); fflush(_logCL2D); }while(0)
#define CLOSE_LOG() fclose(_logCL2D)
#else
#define CREATE_LOG() ;
#define LOG(msg) ;
#define CLOSE_LOG() ;
#endif

/* CL2D Assigned basics ------------------------------------------------ */
std::ostream & operator <<(std::ostream &out, const CL2DAssignment& assigned)
{
    out << "(" << assigned.objId << ", " << assigned.corr << ")";
    return out;
}

CL2DAssignment::CL2DAssignment()
{
    flip = false;
    likelihood = corr = psi = shiftx = shifty = 0;
    objId = BAD_OBJID;
}

void CL2DAssignment::readAlignment(const Matrix2D<double> &M)
{
    double scale;
    transformationMatrix2Parameters2D(M, flip, scale, shiftx, shifty, psi);
}

void CL2DAssignment::copyAlignment(const CL2DAssignment &alignment)
{
    corr = alignment.corr;
    likelihood = alignment.likelihood;
    flip = alignment.flip;
    psi = alignment.psi;
    shiftx = alignment.shiftx;
    shifty = alignment.shifty;
}

/* CL2DClass basics ---------------------------------------------------- */
CL2DClass::CL2DClass()
{
    plans = NULL;
    P.initZeros(prm->Ydim, prm->Xdim);
    P.setXmippOrigin();
    Pupdate = P;
}

CL2DClass::CL2DClass(const CL2DClass &other)
{
    plans = NULL;

    CL2DAssignment assignment;
    assignment.corr = 1;
    Pupdate = other.P;
    nextListImg.push_back(assignment);
    transferUpdate();

    currentListImg = other.currentListImg;
    histClass = other.histClass;
    histNonClass = other.histNonClass;
    neighboursIdx = other.neighboursIdx;
}

CL2DClass::~CL2DClass()
{
    delete plans;
}

void CL2DClass::updateProjection(const MultidimArray<double> &I,
                                 const CL2DAssignment &assigned,
                                 bool force)
{
    if ((assigned.corr > 0 && assigned.objId != BAD_OBJID) || force)
    {
        Pupdate += I;
        nextListImg.push_back(assigned);
    }
}

//#define DEBUG
void CL2DClass::transferUpdate()
{
    if (nextListImg.size() > 0)
    {
        // Take from Pupdate
        double iNq = 1.0 / nextListImg.size();
        Pupdate *= iNq;
        if (prm->normalizeImages)
        	Pupdate.statisticsAdjust(0, 1);
        P = Pupdate;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(P)
        if (!DIRECT_MULTIDIM_ELEM(prm->mask,n))
        	DIRECT_MULTIDIM_ELEM(P,n) = 0;
        Pupdate.initZeros(P);

        // Make sure the image is centered
        centerImage(P,corrAux,rotAux);
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(P)
        if (!DIRECT_A2D_ELEM(prm->mask,i,j))
            DIRECT_A2D_ELEM(P,i,j) = 0;

        // Compute the polar Fourier transform of the full image
        normalizedPolarFourierTransform(P, polarFourierP, false, XSIZE(P) / 5,
                                        XSIZE(P) / 2-2, plans, 1);
        size_t finalSize = 2 * polarFourierP.getSampleNoOuterRing() - 1;
        if (XSIZE(rotationalCorr) != finalSize)
            rotationalCorr.resize(finalSize);
        rotAux.local_transformer.setReal(rotationalCorr);

        // Take the list of images
        currentListImg = nextListImg;
        nextListImg.clear();

        // Update the histogram of corrs of elements inside and outside the class
        if (!prm->classicalMultiref)
        {
            MultidimArray<double> classCorr, nonClassCorr;

            classCorr.resizeNoCopy(currentListImg.size());
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(classCorr)
            DIRECT_A1D_ELEM(classCorr,i) = currentListImg[i].corr;

            nonClassCorr.resizeNoCopy(nextNonClassCorr.size());
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(nonClassCorr)
            DIRECT_A1D_ELEM(nonClassCorr,i) = nextNonClassCorr[i];
            nextNonClassCorr.clear();

            double minC=0., maxC=0.;
            classCorr.computeDoubleMinMax(minC, maxC);
            double minN=0., maxN=0.;
            nonClassCorr.computeDoubleMinMax(minN, maxN);
            double c0 = XMIPP_MIN(minC,minN);
            double cF = XMIPP_MAX(maxC,maxN);
            compute_hist(classCorr, histClass, c0, cF, 200);
            compute_hist(nonClassCorr, histNonClass, c0, cF, 200);
            histClass += 1; // Laplace correction
            histNonClass += 1;
            histClass /= histClass.sum();
            histNonClass /= histNonClass.sum();
        }

#ifdef DEBUG
        histClass.write("PPPclass.txt");
        histNonClass.write("PPPnonClass.txt");
        std::cout << "Histograms have been written. Press any key\n";
        char c;
        std::cin >> c;
#endif

    }
    else
    {
        currentListImg.clear();
        P.initZeros();
    }
}
#undef DEBUG

//#define DEBUG
//#define DEBUG_MORE
void CL2DClass::fitBasic(MultidimArray<double> &I, CL2DAssignment &result,
                         bool reverse)
{
    if (reverse)
    {
        I.selfReverseX();
        I.setXmippOrigin();
    }

    Matrix2D<double> ARS, ASR, R(3, 3);
    ARS.initIdentity(3);
    ASR = ARS;
    MultidimArray<double> IauxSR = I, IauxRS = I;
    Polar<std::complex<double> > polarFourierI;
#ifdef DEBUG_MORE
    Image<double> save2;
    save2()=P;
    save2.write("PPPI1.xmp");
#endif

	// Align the image with the node
	for (int i = 0; i < 3; i++) {
		double shiftX, shiftY, bestRot;

        // Shift then rotate
        bestShift(P, IauxSR, shiftX, shiftY, corrAux);
        MAT_ELEM(ASR,0,2) += shiftX;
        MAT_ELEM(ASR,1,2) += shiftY;
        applyGeometry(LINEAR, IauxSR, I, ASR, IS_NOT_INV, WRAP);
#ifdef DEBUG_MORE
        save2()=IauxSR;
        save2.write("PPPIauxSR_afterShift.xmp");
        std::cout << "ASR\n" << ASR << std::endl;
#endif

        normalizedPolarFourierTransform(IauxSR, polarFourierI, true,
                                        XSIZE(P) / 5, XSIZE(P) / 2-2, plans, 1);

        bestRot = best_rotation(polarFourierP, polarFourierI, rotAux);
        rotation2DMatrix(bestRot, R);
        SPEED_UP_tempsDouble;
        M3x3_BY_M3x3(ASR,R,ASR);
        applyGeometry(LINEAR, IauxSR, I, ASR, IS_NOT_INV, WRAP);
#ifdef DEBUG_MORE
        save2()=IauxSR;
        save2.write("PPPIauxSR_afterShiftAndRotation.xmp");
        std::cout << "ASR\n" << ASR << std::endl;
#endif

        // Rotate then shift
        normalizedPolarFourierTransform(IauxRS, polarFourierI, true,
                                        XSIZE(P) / 5, XSIZE(P) / 2-2, plans, 1);
        bestRot = best_rotation(polarFourierP, polarFourierI, rotAux);
        rotation2DMatrix(bestRot, R);
        M3x3_BY_M3x3(ARS,R,ARS);
        applyGeometry(LINEAR, IauxRS, I, ARS, IS_NOT_INV, WRAP);
#ifdef DEBUG_MORE
        save2()=IauxRS;
        save2.write("PPPIauxRS_afterRotation.xmp");
        std::cout << "ARS\n" << ARS << std::endl;
#endif

        bestShift(P, IauxRS, shiftX, shiftY, corrAux);
        MAT_ELEM(ARS,0,2) += shiftX;
        MAT_ELEM(ARS,1,2) += shiftY;
        applyGeometry(LINEAR, IauxRS, I, ARS, IS_NOT_INV, WRAP);
#ifdef DEBUG_MORE
        save2()=IauxRS;
        save2.write("PPPIauxRS_afterRotationAndShift.xmp");
        std::cout << "ARS\n" << ARS << std::endl;

        char c;
        std::cout << "Press any key\n";
        std::cin >> c;
#endif
    }

    // Compute the correntropy
    double corrRS=0.0, corrSR=0.0;
    const MultidimArray<int> &imask = prm->mask;
    if (prm->useCorrelation)
    {
        corrRS = corrSR = 0;
        long N = 0;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(P)
        {
            if (DIRECT_MULTIDIM_ELEM(imask,n))
            {
                double Pval = DIRECT_MULTIDIM_ELEM(P, n);
                corrRS += Pval * DIRECT_MULTIDIM_ELEM(IauxRS, n);
                corrSR += Pval * DIRECT_MULTIDIM_ELEM(IauxSR, n);
                ++N;
            }
        }
        corrRS /= N;
        corrSR /= N;
    }
    else
    {
        corrRS = fastCorrentropy(P, IauxRS, prm->sigma,
                                 prm->gaussianInterpolator, imask);
        corrSR = fastCorrentropy(P, IauxSR, prm->sigma,
                                 prm->gaussianInterpolator, imask);
    }
#ifdef DEBUG_MORE
       std::cout << "corrRS=" << corrRS << std::endl;
       std::cout << "corrSR=" << corrSR << std::endl;
#endif


    // Prepare result
    if (reverse)
    {
        // COSS: Check that this is correct
        MAT_ELEM(ARS,0,0) *= -1;
        MAT_ELEM(ARS,1,0) *= -1;
        MAT_ELEM(ASR,0,0) *= -1;
        MAT_ELEM(ASR,1,0) *= -1;
    }

    CL2DAssignment candidateRS, candidateSR;
    candidateRS.readAlignment(ARS);
    candidateSR.readAlignment(ASR);
    if (candidateRS.shiftx * candidateRS.shiftx +
    	candidateRS.shifty * candidateRS.shifty > prm->maxShift2)
        candidateRS.corr = 0;
    else
        candidateRS.corr = corrRS;
    if (candidateSR.shiftx * candidateSR.shiftx +
    	candidateSR.shifty * candidateSR.shifty > prm->maxShift2)
        candidateSR.corr = 0;
    else
        candidateSR.corr = corrSR;

    if (candidateRS.corr >= candidateSR.corr)
    {
        I = IauxRS;
        result.copyAlignment(candidateRS);
    }
    else if (corrSR > corrRS)
    {
        I = IauxSR;
        result.copyAlignment(candidateSR);
    }
    else
        result.corr = 0;

#ifdef DEBUG

    Image<double> save;
    save()=P;
    save.write("PPPI1.xmp");
    save()=I;
    save.write("PPPI2.xmp");
    save()=P-I;
    save.write("PPPdiff.xmp");
    std::cout << "sigma=" << prm->sigma << " corr=" << result.corr
    << ". Press" << std::endl;
    char c;
    std::cin >> c;
#endif
}
#undef DEBUG
#undef DEBUG_MORE

void CL2DClass::fit(MultidimArray<double> &I, CL2DAssignment &result)
{
    if (currentListImg.size() == 0)
        return;

    // Try this image
    MultidimArray<double> Idirect = I;
    CL2DAssignment resultDirect;
    fitBasic(Idirect, resultDirect);

    // Try its mirror
    MultidimArray<double> Imirror = I;
    CL2DAssignment resultMirror;
    fitBasic(Imirror, resultMirror, true);

    if (resultMirror.corr > resultDirect.corr)
    {
        I = Imirror;
        result.copyAlignment(resultMirror);
    }
    else
    {
        I = Idirect;
        result.copyAlignment(resultDirect);
    }

    result.likelihood = 0;
    if (!prm->classicalMultiref)
    {
        // Find the likelihood
        int idx;
        histClass.val2index(result.corr, idx);
        if (idx < 0)
            return;
        double likelihoodClass = 0;
        for (int i = 0; i <= idx; i++)
            likelihoodClass += DIRECT_A1D_ELEM(histClass,i);
        histNonClass.val2index(result.corr, idx);
        if (idx < 0)
            return;
        double likelihoodNonClass = 0;
        for (int i = 0; i <= idx; i++)
            likelihoodNonClass += DIRECT_A1D_ELEM(histNonClass,i);
        result.likelihood = likelihoodClass * likelihoodNonClass;
    }
}

/* Look for K neighbours in a list ----------------------------------------- */
//#define DEBUG
void CL2DClass::lookForNeighbours(const std::vector<CL2DClass *> listP, int K)
{
    int Q = listP.size();
    neighboursIdx.clear();
    if (K == Q)
    {
        // As many neighbours as codes
        for (int q = 0; q < Q; q++)
            neighboursIdx.push_back(q);
    }
    else
    {
        MultidimArray<double> distanceCode;
        distanceCode.initZeros(Q);
        CL2DAssignment assignment;
        MultidimArray<double> I;
        for (int q = 0; q < Q; q++)
        {
            if (listP[q] == this)
                distanceCode(q) = 1;
            else
            {
                I = listP[q]->P;
                fit(I, assignment);
                A1D_ELEM(distanceCode,q) = assignment.corr;
            }
        }

        MultidimArray<int> idx;
        distanceCode.indexSort(idx);
        for (int k = 0; k < K; k++)
            neighboursIdx.push_back(idx(Q - k - 1) - 1);
    }
#ifdef DEBUG
    Image<double> save;
    save()=P;
    save.write("PPPthis.xmp");
    for (int k=0; k<K; k++)
    {
        save()=listP[neighboursIdx[k]]->P;
        save.write((std::string)"PPPneigh"+integerToString(k,1));
    }
    std::cout << "Neighbours saved. Press any key\n";
    char c;
    std::cin >> c;
#endif
}
#undef DEBUG

/* Share assignments and classes -------------------------------------- */
void CL2D::shareAssignments(bool shareAssignment, bool shareUpdates,
                            bool shareNonCorr)
{
    if (shareAssignment)
    {
        // Put assignment in common
        std::vector<int> nodeRef;
        if (SF->containsLabel(MDL_REF))
            SF->getColumnValues(MDL_REF, nodeRef);
        else
            nodeRef.resize(SF->size(), -1);

        MPI_Allreduce(MPI_IN_PLACE, &(nodeRef[0]), nodeRef.size(), MPI_INT,
                      MPI_MAX, MPI_COMM_WORLD);
        SF->setColumnValues(MDL_REF, nodeRef);
#ifdef DEBUG_WITH_LOG
    	FileName fnImg;
        FOR_ALL_OBJECTS_IN_METADATA(*SF)
        {
        	int classRef;
        	SF->getValue(MDL_IMAGE,fnImg,__iter.objId);
        	SF->getValue(MDL_REF,classRef,__iter.objId);
        	LOG(formatString("Image %s: class: %d",fnImg.c_str(),classRef).c_str());
        }
#endif
    }

    // Share code updates
    if (shareUpdates)
    {
        std::vector<CL2DAssignment> auxList;
        std::vector<double> auxList2;
        int Q = P.size();
        for (int q = 0; q < Q; q++)
        {
            MPI_Allreduce(MPI_IN_PLACE, MULTIDIM_ARRAY(P[q]->Pupdate),
                          MULTIDIM_SIZE(P[q]->Pupdate), MPI_DOUBLE, MPI_SUM,
                          MPI_COMM_WORLD);

            // Share nextClassCorr and nextNonClassCorr
            std::vector<double> receivedNonClassCorr;
            std::vector<CL2DAssignment> receivedNextListImage;
            int listSize;
            for (size_t rank = 0; rank < prm->node->size; rank++)
            {
                if (rank == prm->node->rank)
                {
                    listSize = P[q]->nextListImg.size();
                    MPI_Bcast(&listSize, 1, MPI_INT, rank, MPI_COMM_WORLD);
                    MPI_Bcast(&(P[q]->nextListImg[0]),
                              P[q]->nextListImg.size() * sizeof(CL2DAssignment),
                              MPI_CHAR, rank, MPI_COMM_WORLD);
                    if (shareNonCorr)
                    {
                        listSize = P[q]->nextNonClassCorr.size();
                        MPI_Bcast(&listSize, 1, MPI_INT, rank, MPI_COMM_WORLD);
                        MPI_Bcast(&(P[q]->nextNonClassCorr[0]),
                                  P[q]->nextNonClassCorr.size(), MPI_DOUBLE,
                                  rank, MPI_COMM_WORLD);
                    }
                }
                else
                {
                    MPI_Bcast(&listSize, 1, MPI_INT, rank, MPI_COMM_WORLD);
                    auxList.resize(listSize);
                    MPI_Bcast(&(auxList[0]), listSize * sizeof(CL2DAssignment),
                              MPI_CHAR, rank, MPI_COMM_WORLD);
                    receivedNextListImage.reserve(
                        receivedNextListImage.size() + listSize);
                    for (int j = 0; j < listSize; j++)
                        receivedNextListImage.push_back(auxList[j]);
                    if (shareNonCorr)
                    {
                        MPI_Bcast(&listSize, 1, MPI_INT, rank, MPI_COMM_WORLD);
                        auxList2.resize(listSize);
                        MPI_Bcast(&(auxList2[0]), listSize, MPI_DOUBLE, rank,
                                  MPI_COMM_WORLD);
                        receivedNonClassCorr.reserve(
                            receivedNonClassCorr.size() + listSize);
                        for (int j = 0; j < listSize; j++)
                            receivedNonClassCorr.push_back(auxList2[j]);
                    }
                }
            }
            // Copy the received elements
            listSize = receivedNextListImage.size();
            P[q]->nextListImg.reserve(P[q]->nextListImg.size() + listSize);
            for (int j = 0; j < listSize; j++)
                P[q]->nextListImg.push_back(receivedNextListImage[j]);
            listSize = receivedNonClassCorr.size();
            P[q]->nextNonClassCorr.reserve(
                P[q]->nextNonClassCorr.size() + listSize);
            for (int j = 0; j < listSize; j++)
                P[q]->nextNonClassCorr.push_back(receivedNonClassCorr[j]);
        }

        transferUpdates();
    }
}

void CL2D::shareSplitAssignments(Matrix1D<int> &assignment, CL2DClass *node1,
                                 CL2DClass *node2) const
{
    // Share assignment
    int imax = VEC_XSIZE(assignment);
    MPI_Allreduce(MPI_IN_PLACE, MATRIX1D_ARRAY(assignment), imax, MPI_INT,
                  MPI_MAX, MPI_COMM_WORLD);

    // Share code updates
    std::vector<CL2DAssignment> auxList;
    std::vector<double> auxList2;
    CL2DClass *node = NULL;
    for (int q = 0; q < 2; q++)
    {
        if (q == 0)
            node = node1;
        else
            node = node2;
        MPI_Allreduce(MPI_IN_PLACE, MULTIDIM_ARRAY(node->Pupdate),
                      MULTIDIM_SIZE(node->Pupdate), MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);

        // Share nextClassCorr and nextNonClassCorr
        std::vector<double> receivedNonClassCorr;
        std::vector<CL2DAssignment> receivedNextListImage;
        int listSize;
        for (size_t rank = 0; rank < prm->node->size; rank++)
        {
            if (rank == prm->node->rank)
            {
                // Transmit node 1 next list of images
                listSize = node->nextListImg.size();
                MPI_Bcast(&listSize, 1, MPI_INT, rank, MPI_COMM_WORLD);
                MPI_Bcast(&(node->nextListImg[0]),
                          node->nextListImg.size() * sizeof(CL2DAssignment),
                          MPI_CHAR, rank, MPI_COMM_WORLD);
                // Transmit node 1 non class corr
                listSize = node->nextNonClassCorr.size();
                MPI_Bcast(&listSize, 1, MPI_INT, rank, MPI_COMM_WORLD);
                MPI_Bcast(&(node->nextNonClassCorr[0]),
                          node->nextNonClassCorr.size(), MPI_DOUBLE, rank,
                          MPI_COMM_WORLD);
            }
            else
            {
                // Receive node 1 next list of images
                MPI_Bcast(&listSize, 1, MPI_INT, rank, MPI_COMM_WORLD);
                auxList.resize(listSize);
                MPI_Bcast(&(auxList[0]), listSize * sizeof(CL2DAssignment),
                          MPI_CHAR, rank, MPI_COMM_WORLD);
                receivedNextListImage.reserve(
                    receivedNextListImage.size() + listSize);
                for (int j = 0; j < listSize; j++)
                    receivedNextListImage.push_back(auxList[j]);
                // Receive node 1 non class corr
                MPI_Bcast(&listSize, 1, MPI_INT, rank, MPI_COMM_WORLD);
                auxList2.resize(listSize);
                MPI_Bcast(&(auxList2[0]), listSize, MPI_DOUBLE, rank,
                          MPI_COMM_WORLD);
                receivedNonClassCorr.reserve(
                    receivedNonClassCorr.size() + listSize);
                for (int j = 0; j < listSize; j++)
                    receivedNonClassCorr.push_back(auxList2[j]);
            }
        }
        // Copy the received elements
        listSize = receivedNextListImage.size();
        node->nextListImg.reserve(node->nextListImg.size() + listSize);
        for (int j = 0; j < listSize; j++)
            node->nextListImg.push_back(receivedNextListImage[j]);
        listSize = receivedNonClassCorr.size();
        node->nextNonClassCorr.reserve(node->nextNonClassCorr.size() + listSize);
        for (int j = 0; j < listSize; j++)
            node->nextNonClassCorr.push_back(receivedNonClassCorr[j]);
    }

    node1->transferUpdate();
    node2->transferUpdate();
}

/* Destructor --------------------------------------------------------- */
CL2D::~CL2D()
{
	int qmax=P.size();
	for (int q=0; q<qmax; q++)
		delete P[q];
}

/* Read image --------------------------------------------------------- */
void CL2D::readImage(Image<double> &I, size_t objId, bool applyGeo) const
{
    if (applyGeo)
        I.readApplyGeo(*SF, objId);
    else
    {
        FileName fnImg;
        SF->getValue(MDL_IMAGE, fnImg, objId);
        I.read(fnImg);
    }
    I().setXmippOrigin();
    if (prm->normalizeImages)
    	I().statisticsAdjust(0, 1);
}

/* CL2D initialization ------------------------------------------------ */
//#define DEBUG
void CL2D::initialize(MetaData &_SF,
                      std::vector<MultidimArray<double> > &_codes0)
{
    if (prm->node->rank == 0)
        std::cout << "Initializing ...\n";

    SF = &_SF;
    Nimgs = SF->size();

    // Start with _Ncodes0 codevectors
    CL2DAssignment assignment;
    assignment.corr = 1;
    bool initialCodesGiven = _codes0.size() > 0;
    for (int q = 0; q < prm->Ncodes0; q++)
    {
        P.push_back(new CL2DClass());
        if (initialCodesGiven)
        {
            P[q]->Pupdate = _codes0[q];
            P[q]->nextListImg.push_back(assignment);
            P[q]->transferUpdate();
        }
    }

    // Estimate sigma and if no previous classes have been given,
    // assign randomly
    prm->sigma = 0;
    if (prm->node->rank == 0)
        init_progress_bar(Nimgs);
    Image<double> I;
    MultidimArray<double> Iaux, Ibest;
    bool oldUseCorrelation = prm->useCorrelation;
    prm->useCorrelation = true; // Since we cannot make the assignment before calculating sigma
    CL2DAssignment bestAssignment;
    size_t first, last;
    SF->fillConstant(MDL_REF,"-1");
    while (prm->taskDistributor->getTasks(first, last))
    {
        for (size_t idx = first; idx <= last; ++idx)
        {
            int q = -1;
            if (!initialCodesGiven)
                q = idx % (prm->Ncodes0);
            size_t objId = prm->objId[idx];
            readImage(I, objId, true);

            // Measure the variance of the signal outside a circular mask
            double avg, stddev;
            I().computeAvgStdev_within_binary_mask(prm->mask, avg, stddev);
            prm->sigma += stddev;

            // Put it randomly in one of the classes
            bestAssignment.objId = assignment.objId = objId;
            if (!initialCodesGiven)
            {
                bestAssignment.corr = 1;
                P[q]->updateProjection(I(), bestAssignment);
            }
            else
            {
                bestAssignment.corr = 0;
                q = -1;
                for (int qp = 0; qp < prm->Ncodes0; qp++)
                {
                    Iaux = I();
                    P[qp]->fit(Iaux, assignment);
                    if (assignment.corr > bestAssignment.corr)
                    {
                        bestAssignment = assignment;
                        Ibest = Iaux;
                        q = qp;
                    }
                }
                if (q != -1)
                    P[q]->updateProjection(Ibest, bestAssignment);
            }
            SF->setValue(MDL_REF, q + 1, objId);
            if (idx % 100 == 0 && prm->node->rank == 0)
                progress_bar(idx);
        }
    }
    if (prm->node->rank == 0)
        progress_bar(Nimgs);
    prm->useCorrelation = oldUseCorrelation;

    // Put sigma in common
    MPI_Allreduce(MPI_IN_PLACE, &(prm->sigma), 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    prm->sigma *= sqrt(2.0) / Nimgs;

    // Share all assignments
    shareAssignments(true, true, false);

    // Now compute the histogram of corr values
    if (!prm->classicalMultiref)
    {
        if (prm->node->rank == 0)
        {
            std::cout << "Computing histogram of correlation values\n";
            init_progress_bar(Nimgs);
        }

        CL2DAssignment inClass, outClass;
        prm->taskDistributor->reset();
        while (prm->taskDistributor->getTasks(first, last))
        {
            for (size_t idx = first; idx <= last; ++idx)
            {
                size_t objId = prm->objId[idx];
                readImage(I, objId, false);

                int q;
                SF->getValue(MDL_REF, q, objId);
                if (q == -1)
                    continue;
                q -= 1;

                outClass.objId = inClass.objId = objId;
                if (q != -1)
                {
                    Iaux = I();
                    P[q]->fit(Iaux, inClass);
                    P[q]->updateProjection(Iaux, inClass);
                    if (prm->Ncodes0 > 1)
                        for (int qp = 0; qp < prm->Ncodes0; qp++)
                        {
                            if (qp == q)
                                continue;
                            Iaux = I();
                            P[qp]->fit(Iaux, outClass);
                            P[qp]->updateNonProjection(outClass.corr);
                        }
                }
                if (prm->node->rank == 0 && idx % 100 == 0)
                    progress_bar(idx);
            }
        }
        if (prm->node->rank == 0)
            progress_bar(Nimgs);

        shareAssignments(false, true, true);
    }
}
#undef DEBUG

/* CL2D write --------------------------------------------------------- */
void CL2D::write(const FileName &fnODir, const FileName &fnRoot, int level) const
{
    int Q = P.size();
    MetaData SFout;
    Image<double> I;
    FileName fnResultsDir=formatString("%s/level_%02d",fnODir.c_str(),level);
    FileName fnOut = formatString("%s/%s_classes.stk",fnResultsDir.c_str(),fnRoot.c_str()), fnClass;
    fnOut.deleteFile();
    for (int q = 0; q < Q; q++)
    {
        fnClass.compose(q + 1, fnOut);
        I() = P[q]->P;
        I.write(fnClass, q, true, WRITE_APPEND);
        size_t id = SFout.addObject();
        SFout.setValue(MDL_REF, q + 1, id);
        SFout.setValue(MDL_IMAGE, fnClass, id);
        SFout.setValue(MDL_CLASS_COUNT,P[q]->currentListImg.size(), id);
    }
    FileName fnSFout = formatString("%s/%s_classes.xmd",fnResultsDir.c_str(),fnRoot.c_str());
    SFout.write("classes@"+fnSFout, MD_APPEND);

    // Make the selfiles of each class
    FileName fnImg;
    MDRow row;
    for (int q = 0; q < Q; q++)
    {
        MetaData SFq;
        std::vector<CL2DAssignment> &currentListImg = P[q]->currentListImg;
        int imax = currentListImg.size();
        for (int i = 0; i < imax; i++)
        {
            const CL2DAssignment &assignment = currentListImg[i];
            SF->getRow(row,assignment.objId);
            row.setValue(MDL_FLIP, assignment.flip);
            row.setValue(MDL_SHIFT_X, assignment.shiftx);
            row.setValue(MDL_SHIFT_Y, assignment.shifty);
            row.setValue(MDL_ANGLE_PSI, assignment.psi);
            SFq.addRow(row);
        }
        MetaData SFq_sorted;
        SFq_sorted.sort(SFq, MDL_IMAGE);
        SFq_sorted.write(formatString("class%06d_images@%s",q+1,fnSFout.c_str()),MD_APPEND);
    }
}

void CL2D::lookNode(MultidimArray<double> &I, int oldnode, int &newnode,
                    CL2DAssignment &bestAssignment)
{
    int Q = P.size();
    int bestq = -1;
    MultidimArray<double> bestImg, Iaux;
    Matrix1D<double> corrList;
    corrList.resizeNoCopy(Q);
    CL2DAssignment assignment;
    bestAssignment.likelihood = bestAssignment.corr = 0;
    size_t objId = bestAssignment.objId;
    for (int q = 0; q < Q; q++)
    {
        // Check if q is neighbour of the oldnode
        bool proceed = false;
        if (oldnode >= 0)
        {
            int imax = P[oldnode]->neighboursIdx.size();
            for (int i = 0; i < imax; i++)
            {
                if (P[oldnode]->neighboursIdx[i] == q)
                {
                    proceed = true;
                    break;
                }
            }
            if (!proceed)
            {
                double threshold = 3.0 * P[oldnode]->currentListImg.size();
                threshold = XMIPP_MAX(threshold,1000);
                threshold = (double) (XMIPP_MIN(threshold,Nimgs)) / Nimgs;
                proceed = (rnd_unif(0, 1) < threshold);
            }
        }
        else
            proceed = true;

		if (proceed) {
			// Try this image
			Iaux = I;
			P[q]->fit(Iaux, assignment);
			VEC_ELEM(corrList,q) = assignment.corr;
			if ((!prm->classicalMultiref && assignment.likelihood > bestAssignment.likelihood) ||
				(prm->classicalMultiref && assignment.corr > bestAssignment.corr) ||
				 prm->classifyAllImages) {
				bestq = q;
				bestImg = Iaux;
				bestAssignment = assignment;
			}
		}
	}

    I = bestImg;
    newnode = bestq;
    bestAssignment.objId = objId;

    // Assign it to the new node and remove it from the rest
    // of nodes if it was among the best
    if (newnode != -1)
    {
        P[newnode]->updateProjection(I, bestAssignment);
        if (!prm->classicalMultiref)
            for (int q = 0; q < Q; q++)
                if (q != newnode && corrList(q) > 0)
                    P[q]->updateNonProjection(corrList(q));
    }
}

void CL2D::transferUpdates()
{
    int Q = P.size();
    for (int q = 0; q < Q; q++)
        P[q]->transferUpdate();
}

/* Run CL2D ------------------------------------------------------------------ */
//#define DEBUG
void CL2D::run(const FileName &fnODir, const FileName &fnOut, int level)
{
	size_t Q = P.size();

    if (prm->node->rank == 0)
        std::cout << "Quantizing with " << Q << " codes...\n";

    std::vector<int> oldAssignment, newAssignment;

    int iter = 1;
    bool goOn = true;
    MetaData MDChanges;
    Image<double> I;
    int progressStep = XMIPP_MAX(1,Nimgs/60);
    CL2DAssignment assignment;
    FileName fnResultsDir=formatString("%s/level_%02d",fnODir.c_str(),level);
    fnResultsDir.makePath(0755);
    while (goOn)
    {
        if (prm->node->rank == 0)
        {
            std::cout << "Iteration " << iter << " ...\n";
            std::cerr << "Iteration " << iter << " ...\n";
            init_progress_bar(Nimgs);
        }

        size_t K = std::min((size_t)prm->Nneighbours+1,Q);
        if (K == 0)
            K = Q;
        for (size_t q = 0; q < Q; q++)
            P[q]->lookForNeighbours(P, K);

        int node;
        double corrSum = 0;
        SF->getColumnValues(MDL_REF, oldAssignment);
        int *ptrOld = &(oldAssignment[0]);
        for (size_t n = 0; n < Nimgs; ++n, ++ptrOld)
            *ptrOld -= 1;
        SF->fillConstant(MDL_REF, "-1");
        prm->taskDistributor->reset();
        size_t first, last;
        while (prm->taskDistributor->getTasks(first, last))
        {
            for (size_t idx = first; idx <= last; ++idx)
            {
                size_t objId = prm->objId[idx];
                readImage(I, objId, false);
                LOG(((String)"Processing image: "+I.name()).c_str());

                assignment.objId = objId;
                lookNode(I(), oldAssignment[idx], node, assignment);
                SF->setValue(MDL_REF, node + 1, objId);
                corrSum += assignment.corr;
                if (prm->node->rank == 0 && idx % progressStep == 0)
                    progress_bar(idx);
            }
        }

        // Gather all pieces computed by nodes
        MPI_Allreduce(MPI_IN_PLACE, &corrSum, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        shareAssignments(true, true, true);

        // Some report
        size_t idMdChanges;
        if (prm->node->rank == 0)
        {
            progress_bar(Nimgs);
            double avgSimilarity = corrSum / Nimgs;
            std::cout << "\nAverage correlation with input vectors="
            << avgSimilarity << std::endl;
            idMdChanges = MDChanges.addObject();
            MDChanges.setValue(MDL_ITER, iter, idMdChanges);
            MDChanges.setValue(MDL_CL2D_SIMILARITY, avgSimilarity, idMdChanges);
        }

        // Count changes
        SF->getColumnValues(MDL_REF, newAssignment);
        int *ptrNew = &(newAssignment[0]);
        for (size_t n = 0; n < Nimgs; ++n, ++ptrNew)
            *ptrNew -= 1;
        size_t Nchanges = 0;
        if (iter > 1)
        {
            int *ptrNew = &(newAssignment[0]);
            ptrOld = &(oldAssignment[0]);
            for (size_t n = 0; n < Nimgs; ++n, ++ptrNew, ++ptrOld)
                if (*ptrNew != *ptrOld)
                    ++Nchanges;
        }
        if (prm->node->rank == 0)
        {
            FileName fnClasses=formatString("%s/%s_classes.xmd",fnResultsDir.c_str(),fnOut.c_str());
        	if (iter==1)
        	{
        		MetaData dummy;
        		dummy.write(formatString("classes@%s",fnClasses.c_str()));
        	}

            std::cout << "Number of assignment changes=" << Nchanges
            << std::endl;
            MDChanges.setValue(MDL_CL2D_CHANGES, (int)Nchanges, idMdChanges);
            MDChanges.write(formatString("info@%s",fnClasses.c_str()),MD_APPEND);
        }

        // Check if there are empty nodes
        if (Q > 1)
        {
            bool smallNodes;
            do
            {
                smallNodes = false;
                int largestNode = -1, sizeLargestNode = -1, smallNode = -1;
                size_t sizeSmallestNode = Nimgs + 1;
                for (size_t q = 0; q < Q; q++)
                {
                    if (P[q]->currentListImg.size() < sizeSmallestNode)
                    {
                        smallNode = q;
                        sizeSmallestNode = P[q]->currentListImg.size();
                    }
                    if ((int) (P[q]->currentListImg.size()) > sizeLargestNode)
                    {
                        sizeLargestNode = P[q]->currentListImg.size();
                        largestNode = q;
                    }
                }
                if (sizeLargestNode==0)
                	REPORT_ERROR(ERR_UNCLASSIFIED,"All classes are of size 0: normally this happens when images are too noisy");
                if (largestNode == -1 || smallNode == -1)
                    break;
                if (sizeSmallestNode < prm->PminSize * Nimgs / Q * 0.01 && sizeSmallestNode<0.25*sizeLargestNode)
                {
                    if (prm->node->rank == 0 && prm->verbose)
                        std::cout << "Splitting node " << largestNode
                        << " by overwriting " << smallNode << std::endl;
                    smallNodes = true;

                    // Clear the old assignment of the images in the small node
                    std::vector<CL2DAssignment> &currentListImgSmall =
                        P[smallNode]->currentListImg;
                    size_t iimax = currentListImgSmall.size();
                    for (size_t ii = 0; ii < iimax; ii++)
                        SF->setValue(MDL_REF, -1, currentListImgSmall[ii].objId);
                    delete P[smallNode];

                    // Clear the old assignment of the images in the large node
                    std::vector<CL2DAssignment> &currentListImgLargest =
                        P[largestNode]->currentListImg;
                    iimax = currentListImgLargest.size();
                    for (size_t ii = 0; ii < iimax; ii++)
                        SF->setValue(MDL_REF, -1,
                                     currentListImgLargest[ii].objId);

                    // Now split the largest node
                    CL2DClass *node1 = new CL2DClass();
                    CL2DClass *node2 = new CL2DClass();
                    std::vector<size_t> splitAssignment;
                    splitNode(P[largestNode], node1, node2, splitAssignment);
                    delete P[largestNode];

                    // Keep the results of the split
                    P[largestNode] = node1;
                    P[smallNode] = node2;
                    iimax = splitAssignment.size();
                    for (size_t ii = 0; ii < iimax; ii += 2)
                    {
                        if (splitAssignment[ii + 1] == 1)
                            SF->setValue(MDL_REF, largestNode + 1,
                                         splitAssignment[ii]);
                        else
                            SF->setValue(MDL_REF, smallNode + 1,
                                         splitAssignment[ii]);
                    }
                }
            }
            while (smallNodes);
        }

        if (prm->node->rank == 0)
            write(fnODir,fnOut,level);

        if ((iter > 1 && Nchanges < 0.005 * Nimgs && Q > 1) || iter >= prm->Niter)
            goOn = false;
        iter++;
    }

    std::sort(P.begin(), P.end(), SDescendingClusterSort());
}

/* Clean ------------------------------------------------------------------- */
int CL2D::cleanEmptyNodes()
{
    int retval = 0;
    std::vector<CL2DClass *>::iterator ptr = P.begin();
    while (ptr != P.end())
        if ((*ptr)->currentListImg.size() == 0)
        {
            ptr = P.erase(ptr);
            retval++;
        }
        else
            ptr++;
    return retval;
}

/* Split ------------------------------------------------------------------- */
//#define DEBUG
void CL2D::splitNode(CL2DClass *node, CL2DClass *&node1, CL2DClass *&node2,
                     std::vector<size_t> &splitAssignment) const
{
	const CL2DClass* firstNode=node;
    std::vector<CL2DClass *> toDelete;
    Matrix1D<int> newAssignment, oldAssignment, firstSplitAssignment;
    Image<double> I;
    MultidimArray<double> Iaux1, Iaux2, corrList;
    MultidimArray<int> idx;
    CL2DAssignment assignment, assignment1, assignment2;
    CL2DClass *firstSplitNode1 = NULL;
    CL2DClass *firstSplitNode2 = NULL;
    size_t minAllowedSize = (size_t)(prm->PminSize * 0.01 * node->currentListImg.size());

    bool oldclassicalMultiref = prm->classicalMultiref;
    prm->classicalMultiref = false;
    bool finish;
    bool success = true;
    do
    {
#ifdef DEBUG
	std::cout << "Splitting node " << node << "(" << node->currentListImg.size() << ") into " << node1 << " and " << node2 << std::endl;
#endif

        finish = true;
        node2->neighboursIdx = node1->neighboursIdx = node->neighboursIdx;
        node2->P = node1->P = node->P;

        size_t imax = node->currentListImg.size();
        if (imax < minAllowedSize)
        {
            toDelete.push_back(node1);
            toDelete.push_back(node2);
#ifdef DEBUG
            std::cout << "Pushing to delete " << node1 << std::endl;
            std::cout << "Pushing to delete " << node2 << std::endl;
#endif
            success = false;
            break;
        }

        // Compute the corr histogram
        if (prm->node->rank == 0 && prm->verbose >= 2)
            std::cerr << "Calculating corr distribution at split ..."
            << std::endl;
        corrList.initZeros(imax);
        for (size_t i = 0; i < imax; i++)
        {
            if ((i + 1) % (prm->node->size) == prm->node->rank)
            {
                readImage(I, node->currentListImg[i].objId, false);
                node->fit(I(), assignment);
                A1D_ELEM(corrList,i) = assignment.corr;
            }
            if (prm->node->rank == 0 && i % 25 == 0 && prm->verbose >= 2)
                progress_bar(i);
        }
        if (prm->node->rank == 0 && prm->verbose >= 2)
            progress_bar(imax);
        MPI_Allreduce(MPI_IN_PLACE, MULTIDIM_ARRAY(corrList), imax, MPI_DOUBLE,
                      MPI_MAX, MPI_COMM_WORLD);
        newAssignment.initZeros(imax);

        // Compute threshold
        corrList.indexSort(idx);
        double corrThreshold = corrList(idx(XSIZE(idx)/2)-1);
#ifdef DEBUG
        corrList.printStats();
        std::cout << " corrThreshold= " << corrThreshold << std::endl;
#endif
        if (corrThreshold == 0)
        {
            if (firstSplitNode1 != NULL)
            {
                toDelete.push_back(node1);
                toDelete.push_back(node2);
#ifdef DEBUG
            std::cout << "Pushing to delete " << node1 << std::endl;
            std::cout << "Pushing to delete " << node2 << std::endl;
#endif
                success = false;
                break;
            }
            else
            {
                // Split at random
                for (size_t i = 0; i < imax; i++)
                {
                    assignment.objId = node->currentListImg[i].objId;
                    readImage(I, assignment.objId, false);
                    node->fit(I(), assignment);
                    if ((i + 1) % 2 == 0)
                    {
                        node1->updateProjection(I(), assignment,true);
                        VEC_ELEM(newAssignment,i) = 1;
                        node2->updateNonProjection(assignment.corr, true);
                    }
                    else
                    {
                        node2->updateProjection(I(), assignment,true);
                        VEC_ELEM(newAssignment,i) = 2;
                        node1->updateNonProjection(assignment.corr, true);
                    }
                }
                node1->transferUpdate();
                node2->transferUpdate();
#ifdef DEBUG
            std::cout << "Splitting at random " << node1->currentListImg.size() << " " << node2->currentListImg.size() << std::endl;
#endif
                success = true;
                break;
            }
        }
        else
        {
            // Split according to corr
            if (prm->node->rank == 0 && prm->verbose >= 2)
                std::cerr << "Splitting by corr threshold ..." << std::endl;
            for (size_t i = 0; i < imax; i++)
            {
                if ((i + 1) % (prm->node->size) == prm->node->rank)
                {
                    assignment.objId = node->currentListImg[i].objId;
                    readImage(I, assignment.objId, false);
                    node->fit(I(), assignment);
                    if (assignment.corr < corrThreshold)
                    {
                        node1->updateProjection(I(), assignment);
                        VEC_ELEM(newAssignment,i) = 1;
                        node2->updateNonProjection(assignment.corr);
                    }
                    else
                    {
                        node2->updateProjection(I(), assignment);
                        VEC_ELEM(newAssignment,i) = 2;
                        node1->updateNonProjection(assignment.corr);
                    }
                }
                if (prm->node->rank == 0 && i % 25 == 0 && prm->verbose >= 2)
                    progress_bar(i);
            }
            if (prm->node->rank == 0 && prm->verbose >= 2)
                progress_bar(imax);
            shareSplitAssignments(newAssignment, node1, node2);
        }

#ifdef DEBUG
	std::cout << "After first split Node1: " << node1->currentListImg.size() << " Node2 " << node2->currentListImg.size() << std::endl;
#endif

	// Backup the first split in case it fails
        if (firstSplitNode1 == NULL)
        {
#ifdef DEBUG
	std::cout << "Creating backup\n";
#endif
            firstSplitAssignment = newAssignment;
            firstSplitNode1 = new CL2DClass(*node1);
            firstSplitNode2 = new CL2DClass(*node2);
        }

        // Split iterations
        for (int it = 0; it < prm->Niter; it++)
        {
            if (prm->node->rank == 0 && prm->verbose >= 2)
            {
                init_progress_bar(imax);
            }

            oldAssignment = newAssignment;
            newAssignment.initZeros();
            for (size_t i = 0; i < imax; i++)
            {
                if ((i + 1) % (prm->node->size) == prm->node->rank)
                {
                    // Read image
                    assignment1.objId = assignment2.objId = assignment.objId
                                                            = node->currentListImg[i].objId;
                    readImage(I, assignment.objId, false);

                    Iaux1 = I();
                    node1->fit(Iaux1, assignment1);
                    Iaux2 = I();
                    node2->fit(Iaux2, assignment2);

                    if (assignment1.likelihood > assignment2.likelihood)
                    {
                        node1->updateProjection(Iaux1, assignment1);
                        VEC_ELEM(newAssignment,i) = 1;
                        node2->updateNonProjection(assignment2.corr);
                    }
                    else if (assignment2.likelihood > assignment1.likelihood)
                    {
                        node2->updateProjection(Iaux2, assignment2);
                        VEC_ELEM(newAssignment,i) = 2;
                        node1->updateNonProjection(assignment1.corr);
                    }
                }
                if (prm->node->rank == 0 && i % 25 == 0 && prm->verbose >= 2)
                    progress_bar(i);
            }
            if (prm->node->rank == 0 && prm->verbose >= 2)
                progress_bar(imax);
            shareSplitAssignments(newAssignment, node1, node2);

#ifdef DEBUG
	std::cout << "After refinement iteration: " << it << " Node1: " << node1->currentListImg.size() << " Node2 " << node2->currentListImg.size() << std::endl;
#endif

	        int Nchanges = 0;
            FOR_ALL_ELEMENTS_IN_MATRIX1D(newAssignment)
            if (newAssignment(i) != oldAssignment(i))
                Nchanges++;
            if (prm->node->rank == 0 && prm->verbose >= 2)
                std::cout << "Number of assignment split changes=" << Nchanges
                << std::endl;

            // Check if one of the nodes is too small
            if (node1->currentListImg.size() < minAllowedSize
                || node2->currentListImg.size() < minAllowedSize
                || Nchanges < 0.005 * imax)
                break;
        }

        if (node1->currentListImg.size() < minAllowedSize && node2->currentListImg.size() < minAllowedSize)
        {
            if (prm->node->rank == 0 && prm->verbose >= 2)
                std::cout << "Removing both nodes, they are too small "
                << node1->currentListImg.size() << " "
                << minAllowedSize << "...\n";
            if (node1 != node)
                delete node1;
            if (node2 != node)
                delete node2;
            success = false;
            finish = true;
        }
        else if (node1->currentListImg.size() < minAllowedSize)
        {
            if (prm->node->rank == 0 && prm->verbose >= 2)
                std::cout << "Removing node1, it's too small "
                << node1->currentListImg.size() << " "
                << minAllowedSize << "...\n";
            if (node1 != node)
                delete node1;
            node1 = new CL2DClass();
            toDelete.push_back(node2);
#ifdef DEBUG
            std::cout << "Pushing to delete " << node2 << std::endl;
#endif
            node = node2;
            node2 = new CL2DClass();
            finish = false;
        }
        else if (node2->currentListImg.size() < minAllowedSize)
        {
            if (prm->node->rank == 0 && prm->verbose >= 2)
                std::cout << "Removing node2, it's too small "
                << node2->currentListImg.size() << " "
                << minAllowedSize << "...\n";
            if (node2 != node)
                delete node2;
            node2 = new CL2DClass();
            toDelete.push_back(node1);
#ifdef DEBUG
            std::cout << "Pushing to delete " << node1 << std::endl;
#endif
            node = node1;
            node1 = new CL2DClass();
            finish = false;
        }
    }
    while (!finish);
    for (size_t i = 0; i < toDelete.size(); i++)
        if (toDelete[i] != firstNode)
        {
#ifdef DEBUG
            std::cout << "deleting from list " << toDelete[i] << std::endl;
#endif
            delete toDelete[i];
        }

    if (success)
    {
        for (size_t i = 0; i < node1->currentListImg.size(); i++)
        {
            splitAssignment.push_back(node1->currentListImg[i].objId);
            splitAssignment.push_back(1);
        }
        for (size_t i = 0; i < node2->currentListImg.size(); i++)
        {
            splitAssignment.push_back(node2->currentListImg[i].objId);
            splitAssignment.push_back(2);
        }
        delete firstSplitNode1;
        delete firstSplitNode2;
    }
    else
    {
        node1 = firstSplitNode1;
        node2 = firstSplitNode2;
        splitAssignment.reserve(VEC_XSIZE(firstSplitAssignment));
        FOR_ALL_ELEMENTS_IN_MATRIX1D(firstSplitAssignment)
        splitAssignment.push_back(VEC_ELEM(firstSplitAssignment,i));
    }
    prm->classicalMultiref = oldclassicalMultiref;
}
#undef DEBUG

void CL2D::splitFirstNode()
{
    std::sort(P.begin(), P.end(), SDescendingClusterSort());
    int Q = P.size();
    P.push_back(new CL2DClass());
    P.push_back(new CL2DClass());
    std::vector<size_t> splitAssignment;
    splitNode(P[0], P[Q], P[Q + 1], splitAssignment);
    delete P[0];
    P[0] = NULL;
    P.erase(P.begin());
}

/* MPI constructor --------------------------------------------------------- */
ProgClassifyCL2D::ProgClassifyCL2D(int argc, char** argv)
{
    node = new MpiNode(argc, argv);
    if (!node->isMaster())
        verbose = 0;
    taskDistributor = NULL;
}

/* Destructor -------------------------------------------------------------- */
ProgClassifyCL2D::~ProgClassifyCL2D()
{
    delete node;
    delete taskDistributor;
}

/* VQPrm I/O --------------------------------------------------------------- */
void ProgClassifyCL2D::readParams()
{
    fnSel = getParam("-i");
    fnOut = getParam("--oroot");
    fnODir = getParam("--odir");
    fnCodes0 = getParam("--ref0");
    Niter = getIntParam("--iter");
    Nneighbours = getIntParam("--neigh");
    Ncodes0 = getIntParam("--nref0");
    Ncodes = getIntParam("--nref");
    PminSize = getDoubleParam("--minsize");
    String aux;
    aux = getParam("--distance");
    useCorrelation = aux == "correlation";
    classicalMultiref = checkParam("--classicalMultiref");
    maxShift = getDoubleParam("--maxShift");
	classifyAllImages = checkParam("--classifyAllImages");
	normalizeImages = !checkParam("--dontNormalizeImages");
}

void ProgClassifyCL2D::show() const {
	if (!verbose)
		return;
	std::cout << "Input images:            " << fnSel << std::endl
			<< "Output root:             " << fnOut << std::endl
			<< "Output dir:              " << fnODir << std::endl
			<< "Iterations:              " << Niter << std::endl
			<< "CodesSel0:               " << fnCodes0 << std::endl
			<< "Codes0:                  " << Ncodes0 << std::endl
			<< "Codes:                   " << Ncodes << std::endl
			<< "Neighbours:              " << Nneighbours << std::endl
			<< "Minimum node size:       " << PminSize << std::endl
			<< "Use Correlation:         " << useCorrelation << std::endl
			<< "Classical Multiref:      " << classicalMultiref << std::endl
			<< "Maximum shift:           " << maxShift << std::endl
			<< "Classify all images:     " << classifyAllImages << std::endl
			<< "Normalize images:        " << normalizeImages << std::endl
	;
}

void ProgClassifyCL2D::defineParams()
{
    addUsageLine("Divide a selfile into the desired number of classes. ");
    addUsageLine("+Vector quantization with correntropy and a probabilistic criterion is used for creating the subdivisions.");
    addUsageLine("+Correlation and the standard maximum correlation criterion can also be used and normally produce good results.");
    addUsageLine("+Correntropy and the probabilistic clustering criterion are recommended for images with very low SNR or cases in which the correlation have clear difficulties to converge.");
    addUsageLine("+");
    addUsageLine("+The algorithm is fully described in [[http://www.ncbi.nlm.nih.gov/pubmed/20362059][this article]].");
    addUsageLine("+");
    addUsageLine("+An interesting convergence criterion is the number of images changing classes between iterations. If a low percentage of the image change class, then the clustering is rather stable and clear.");
    addUsageLine("+If many images change class, it is likely that there is not enough SNR to determine so many classes. It is recommended to reduce the number of classes");
    addSeeAlsoLine("mpi_image_sort");
    addParamsLine("    -i <selfile>             : Selfile with the input images");
    addParamsLine("   [--odir <dir=\".\">]      : Output directory");
    addParamsLine("   [--oroot <root=class>]    : Output rootname, by default, class");
    addParamsLine("   [--iter <N=20>]           : Number of iterations");
    addParamsLine("   [--nref0 <N=2>]           : Initial number of code vectors");
    addParamsLine("or  --ref0 <selfile=\"\">    : Selfile with initial code vectors");
    addParamsLine("   [--nref <N=16>]           : Final number of code vectors");
    addParamsLine("   [--neigh+ <N=4>]          : Number of neighbour code vectors");
    addParamsLine("                             : Set -1 for all");
    addParamsLine("   [--minsize+ <N=20>]       : Percentage minimum node size");
    addParamsLine("                             : Let's say that we have 1000 images and we are currently estimating 2 classes.");
    addParamsLine("                             : On average each class should have about 500 images. ");
    addParamsLine("                             : If one of the classes drops below 20% (by default) of 500,");
    addParamsLine("                             : then that class is removed and the remaining class is splitted in two. ");
    addParamsLine("                             : The same reasoning is applied when the number of classes is 4, but 20% ");
    addParamsLine("                             : is calculated now on 250 (=1000/4).");
    addParamsLine("   [--distance <type=correntropy>]       : Distance type");
    addParamsLine("            where <type>");
    addParamsLine("                       correntropy correlation: See CL2D paper for the definition of correntropy");
    addParamsLine("   [--classicalMultiref]     : Instead of enhanced clustering");
    addParamsLine("   [--maxShift <d=10>]       : Maximum allowed shift");
	addParamsLine("   [--classifyAllImages]     : By default, some images may not be classified. Use this option to classify them all.");
	addParamsLine("   [--dontNormalizeImages]   : By default, input images are normalized to have 0 mean and standard deviation 1");
    addExampleLine("mpirun -np 3 `which xmipp_mpi_classify_CL2D` -i images.stk --nref 256 --oroot class --odir CL2Dresults --iter 10");
}

void ProgClassifyCL2D::produceSideInfo()
{
    maxShift2 = maxShift * maxShift;

    gaussianInterpolator.initialize(6, 60000, false);

    // Get image dimensions
    SF.read(fnSel);
    SF.removeDisabled();

    size_t Zdim, Ndim;
    getImageSize(SF, Xdim, Ydim, Zdim, Ndim);

    // Prepare the Task distributor
    SF.findObjects(objId);
    size_t Nimgs = objId.size();
    taskDistributor = new FileTaskDistributor(Nimgs,
                      XMIPP_MAX(1,Nimgs/(5*node->size)), node);

    // Prepare mask for evaluating the noise outside
    mask.resize(prm->Ydim, prm->Xdim);
    mask.setXmippOrigin();
    BinaryCircularMask(mask, prm->Xdim / 2, INNER_MASK);

    // Read input code vectors if available
    std::vector<MultidimArray<double> > codes0;
    if (fnCodes0 != "")
    {
        Image<double> I;
        MetaData SFCodes(fnCodes0);

        size_t Xdim0, Ydim0, Zdim0, Ndim0;
        getImageSize(SFCodes, Xdim0, Ydim0, Zdim0, Ndim0);
        if (Xdim0!=Xdim || Ydim0!=Ydim)
        	REPORT_ERROR(ERR_MULTIDIM_SIZE,"Input reference and images are not of the same size");

        FOR_ALL_OBJECTS_IN_METADATA(SFCodes)
        {
            I.readApplyGeo(SFCodes, __iter.objId);
            I().setXmippOrigin();
            codes0.push_back(I());
        }
        Ncodes0 = codes0.size();
    }
    vq.initialize(SF, codes0);
}

void ProgClassifyCL2D::run()
{
	CREATE_LOG();
    show();
    produceSideInfo();

    // Run all iterations
    int level = 0;
    int Q = vq.P.size();
    bool originalClassicalMultiref=classicalMultiref;
    if (Q==1)
    	classicalMultiref=true;
    vq.run(fnODir, fnOut, level);

    while (Q < Ncodes)
    {
        if (node->rank == 0)
            std::cout << "Spliting nodes ...\n";

        int Nclean = vq.cleanEmptyNodes();
        int Nsplits = XMIPP_MIN(Q,Ncodes-Q) + Nclean;

        for (int i = 0; i < Nsplits; i++)
        {
            vq.splitFirstNode();
            if (node->rank == 0)
                std::cout << "Currently there are " << vq.P.size() << " nodes"
                << std::endl;
        }

        Q = vq.P.size();
        level++;
        classicalMultiref=originalClassicalMultiref || Q==1;
        vq.run(fnODir, fnOut, level);
    }
    if (node->rank == 0)
    {
        std::sort(vq.P.begin(), vq.P.end(), SDescendingClusterSort());
        Q = vq.P.size();
        MetaData SFq, SFclassified, SFaux, SFaux2;
        for (int q = 0; q < Q; q++)
        {
            SFq.read(formatString("class%06d_images@%s/level_%02d/%s_classes.xmd",q+1,fnODir.c_str(),level,fnOut.c_str()));
            SFq.fillConstant(MDL_REF, integerToString(q + 1));
            SFq.fillConstant(MDL_ENABLED, "1");
            SFclassified.unionAll(SFq);
        }
        SFaux = SF;
        SFaux.subtraction(SFclassified, MDL_IMAGE);
        SFaux.fillConstant(MDL_ENABLED, "-1");
        SFaux2.join(SFclassified, SF, MDL_IMAGE, LEFT);
        SFclassified.clear();
        SFaux2.unionAll(SFaux);
        SFaux.clear();
        SFaux.sort(SFaux2, MDL_IMAGE);
        SFaux.write(fnODir+"/images.xmd");
    }
    CLOSE_LOG();
}

/* Main -------------------------------------------------------------------- */
int main(int argc, char** argv)
{
    ProgClassifyCL2D progprm(argc, argv);
    progprm.read(argc, argv);
    prm = &progprm;
    int e = progprm.tryRun();
    return e;
}
