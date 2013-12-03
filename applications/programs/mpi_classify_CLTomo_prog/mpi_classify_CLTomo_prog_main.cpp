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

#include "mpi_classify_CLTomo.h"
#include <data/mask.h>
#include <data/polar.h>
#include <data/xmipp_image_generic.h>
#include <reconstruction/symmetrize.h>

// Pointer to parameters
ProgClassifyCL3D *prm = NULL;
FILE * _logCL3D = NULL;

//#define DEBUG_WITH_LOG
#ifdef DEBUG_WITH_LOG
#define CREATE_LOG() _logCL3D = fopen(formatString("nodo%02d.log", node->rank).c_str(), "w+")
#define LOG(msg) do{fprintf(_logCL3D, "%s\t%s\n", getCurrentTimeString(), msg.c_str()); fflush(_logCL3D); }while(0)
#define CLOSE_LOG() fclose(_logCL3D)
#else
#define CREATE_LOG() ;
#define LOG(msg) ;
#define CLOSE_LOG() ;
#endif

/* CL3D Assigned basics ------------------------------------------------ */
std::ostream & operator <<(std::ostream &out, const CL3DAssignment& assigned)
{
    out << "(" << assigned.objId << " Angles=" << assigned.rot << "," << assigned.tilt << "," << assigned.psi
    << " Shifts=" << assigned.shiftx << "," << assigned.shifty << "," << assigned.shiftz << " -> " << assigned.score << ")";
    return out;
}

CL3DAssignment::CL3DAssignment()
{
    score = psi = tilt = rot = shiftx = shifty = shiftz = 0;
    objId = BAD_OBJID;
}

void CL3DAssignment::readAlignment(const Matrix2D<double> &M)
{
    double scale;
    bool flip;
    transformationMatrix2Parameters3D(M, flip, scale, shiftx, shifty, shiftz, rot, tilt, psi);
    rot=realWRAP(rot,-180,180);
    tilt=realWRAP(tilt,-180,180);
    psi=realWRAP(psi,-180,180);
}

void CL3DAssignment::copyAlignment(const CL3DAssignment &alignment)
{
    score = alignment.score;
    psi = alignment.psi;
    tilt = alignment.tilt;
    rot = alignment.rot;
    shiftx = alignment.shiftx;
    shifty = alignment.shifty;
    shiftz = alignment.shiftz;
}

/* CL3DClass basics ---------------------------------------------------- */
CL3DClass::CL3DClass()
{
    Paux.initZeros(prm->Zdim, prm->Ydim, prm->Xdim);
    Paux.setXmippOrigin();
    transformer.setReal(Paux);
    //PupdateReal=Paux;
    Pupdate.initZeros(transformer.fFourier);
    PupdateMask.initZeros(Pupdate);
    pyIfourierMaskFRM=NULL;
    //weightSum=0;
}

CL3DClass::CL3DClass(const CL3DClass &other)
{
    CL3DAssignment assignment;
    assignment.score = -1e38;
    updateProjection((MultidimArray<double> &)other.P,assignment);
    transferUpdate();

    currentListImg = other.currentListImg;
    neighboursIdx = other.neighboursIdx;
}

//#define DEBUG
void CL3DClass::updateProjection(MultidimArray<double> &I,
                                 const CL3DAssignment &assigned,
                                 bool force)
{
    if ((fabs(assigned.score) <=1 && assigned.objId != BAD_OBJID && assigned.score>0) || force)
    {
    	constructFourierMask(I);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Ifourier)
        if (DIRECT_MULTIDIM_ELEM(IfourierMask,n))
        {
            double *ptrIfourier=(double*)&DIRECT_MULTIDIM_ELEM(Ifourier,n);
            double *ptrPupdate=(double*)&DIRECT_MULTIDIM_ELEM(Pupdate,n);
            *(ptrPupdate)+=(*(ptrIfourier))*assigned.score;
            *(ptrPupdate+1)+=(*(ptrIfourier+1))*assigned.score;
            DIRECT_MULTIDIM_ELEM(PupdateMask,n)+=assigned.score;
        }
        //PupdateReal+=I*assigned.score;
        //weightSum+=assigned.score;
        nextListImg.push_back(assigned);

#ifdef DEBUG
        Image<double> save;
        save()=I;
        save.write("PPPupdateI.vol");
        save()=PupdateReal;
        save.write("PPPupdateUpdateReal.vol");
        save()=P;
        if (MULTIDIM_SIZE(P)>0)
        	save.write("PPPupdateP.vol");

        // Take from Pupdate
        MultidimArray< std::complex<double> > PfourierAux;
        PfourierAux=Pupdate;
        save().initZeros(I);
        double *ptrPupdate=(double*)&DIRECT_MULTIDIM_ELEM(PfourierAux,0);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(PupdateMask)
        {
            double maskVal=DIRECT_MULTIDIM_ELEM(PupdateMask,n);
            if (maskVal>0)
            {
                double iMask=1./maskVal;
                *(ptrPupdate)*=iMask;
                *(ptrPupdate+1)*=iMask;
            }
            ptrPupdate+=2;
        }
        transformer.inverseFourierTransform(PfourierAux,save());
        save.write("PPPupdateUpdateFourier.vol");
        std::cout << "Updating. Press any key\n";
        char c;
        std::cin >> c;
#endif
    }
}
#undef DEBUG

//#define DEBUG
void CL3DClass::transferUpdate()
{
    if (nextListImg.size() > 0)
    {
        // Take from Pupdate
        double *ptrPupdate=(double*)&DIRECT_MULTIDIM_ELEM(Pupdate,0);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(PupdateMask)
        {
            double maskVal=DIRECT_MULTIDIM_ELEM(PupdateMask,n);
            if (maskVal>0 && DIRECT_MULTIDIM_ELEM(prm->maxFreqMask,n))
            {
                double iMask=1./maskVal;
                *(ptrPupdate)*=iMask;
                *(ptrPupdate+1)*=iMask;
            }
            else
            {
            	*(ptrPupdate)=*(ptrPupdate+1)=0.;
            }
            ptrPupdate+=2;
        }
        Pfourier=Pupdate;

        transformer.inverseFourierTransform(Pupdate,Paux);
        //P=PupdateReal/weightSum;

        // Compact support in real space
#ifdef DEBUG
        Image<double> save;
        save()=Paux;
        save.write("PPPtransfer0.vol");
#endif

        double mean;
        detectBackground(Paux,bgMask,0.01,mean);
        if (bgMask.sum()<0.9*(XSIZE(Paux)-4)*(YSIZE(Paux)-4)*(ZSIZE(Paux)-4))
        {
            FOR_ALL_ELEMENTS_IN_ARRAY3D(Paux)
            if (A3D_ELEM(bgMask,k,i,j))
                A3D_ELEM(Paux,k,i,j) = 0;
            else if (k<=STARTINGZ(Paux)+2  || i<=STARTINGY(Paux)+2  || j<=STARTINGX(Paux)+2 ||
        		     k>=FINISHINGZ(Paux)-2 || i>=FINISHINGY(Paux)-2 || j>=FINISHINGX(Paux)-2)
        	    A3D_ELEM(Paux,k,i,j) = 0;
        }

#ifdef DEBUG
        save()=Paux;
        save.write("PPPtransfer1.vol");
#endif
        // Compact support in wavelet space
        forceDWTSparsity(Paux,1.0-prm->DWTsparsity);

#ifdef DEBUG
        save()=Paux;
        save.write("PPPtransfer2.vol");
#endif
        // Symmetrize
        symmetrizeVolume(prm->SL,Paux,P,WRAP);
#ifdef DEBUG
        save()=P;
        save.write("PPPtransfer3.vol");
#endif

        // Normalize and remove outside sphere
        P.statisticsAdjust(0, 1);
        const MultidimArray<int> &mask=prm->mask.get_binary_mask();
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(P)
        if (!DIRECT_MULTIDIM_ELEM(mask,n))
            DIRECT_MULTIDIM_ELEM(P,n) = 0;
        Pupdate.initZeros();
        //PupdateReal.initZeros();
        PupdateMask.initZeros();
        //weightSum=0;

#ifdef DEBUG
        save()=P;
        save.write("PPPtransfer4.vol");
#endif
        /* COSS: STILL TO PUT IN 3D
        // Make sure the image is centered
        centerImage(P,corrAux,rotAux);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(P)
        if (!DIRECT_MULTIDIM_ELEM(mask,n))
            DIRECT_MULTIDIM_ELEM(P,n) = 0;
        */

        // Take the list of images
        currentListImg = nextListImg;
        nextListImg.clear();
    }
    else
    {
        currentListImg.clear();
        P.initZeros();
    }
#ifdef DEBUG
        std::cout << "Transfer finished. Press any key\n";
        char c;
        std::cin >> c;
#endif
}
#undef DEBUG

void CL3DClass::constructFourierMask(MultidimArray<double> &I)
{
    transformer.FourierTransform(I,Ifourier,false);

    FFT_magnitude(Ifourier,IfourierMag);
    IfourierMag.resize(1,1,1,MULTIDIM_SIZE(IfourierMag));
    IfourierMag.sort(IfourierMagSorted);
    double minMagnitude=A1D_ELEM(IfourierMagSorted,(int)(prm->sparsity*XSIZE(IfourierMag)));

    IfourierMask.initZeros(Ifourier);
    for (size_t n=1; n<MULTIDIM_SIZE(Ifourier); ++n)
        if (DIRECT_MULTIDIM_ELEM(IfourierMag,n)>minMagnitude)
        	DIRECT_MULTIDIM_ELEM(IfourierMask,n)=1;
}

void CL3DClass::constructFourierMaskFRM()
{
	int xdim=(int)(prm->Xdim);
	int xdim_2=(int)xdim/2;
	IfourierMaskFRM.initZeros(prm->Xdim,prm->Xdim,prm->Xdim);
	IfourierMaskFRM.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY3D(IfourierMask)
	if (A3D_ELEM(IfourierMask,k,i,j) && A3D_ELEM(prm->maxFreqMask,k,i,j))
	{
		int kp=k, ip=i;
		if (k>xdim_2)
			kp-=xdim;
		if (i>xdim_2)
			ip-=xdim;
		A3D_ELEM(IfourierMaskFRM,j,ip,kp)=A3D_ELEM(IfourierMaskFRM,-j,-ip,-kp)=1;
	}

    PyObject *pyMask=convertToNumpy(IfourierMaskFRM);
	PyObject *arglist = Py_BuildValue("(Oi)", pyMask,0);
	if (pyIfourierMaskFRM!=NULL)
		Py_DECREF(pyIfourierMaskFRM);
	pyIfourierMaskFRM = PyObject_CallObject(prm->wedgeClass, arglist);
	/*if (PyErr_Occurred()!=NULL)
		PyErr_Print();*/
	Py_DECREF(arglist);
	Py_DECREF(pyMask);
}

//#define DEBUG
void CL3DClass::fitBasic(MultidimArray<double> &I, CL3DAssignment &result)
{
#ifdef DEBUG
    Image<double> save2;
    save2()=I;
    save2.write("PPPfitBasicI0.xmp");
    save2()=P;
    save2.write("PPPfitBasicI1.xmp");
#endif

    if (currentListImg.size() == 0)
        return;

    Matrix2D<double> A;
    double frmScore;
    constructFourierMask(I);
    constructFourierMaskFRM();
    if (!prm->dontAlign)
		alignVolumesFRM(prm->frmFunc, P, I, pyIfourierMaskFRM, result.rot, result.tilt, result.psi, result.shiftx, result.shifty, result.shiftz,
				frmScore,A,prm->maxShift, prm->maxFreq);
    if (fabs(result.shiftx)>prm->maxShiftX || fabs(result.shifty)>prm->maxShiftY || fabs(result.shiftz)>prm->maxShiftZ ||
    	fabs(result.rot)>prm->maxRot || fabs(result.tilt)>prm->maxTilt || fabs(result.psi)>prm->maxPsi)
    	result.score=-1e38;
    else
    {
		applyGeometry(LINEAR, Iaux, I, A, IS_NOT_INV, DONT_WRAP);
		apply_binary_mask(prm->mask.get_binary_mask(),Iaux,I,0.0);

		constructFourierMask(I);
		result.score=frmScore;
    }

#ifdef DEBUG

    Image<double> save;
    save()=I;
    save.write("PPPfitBasicI2.xmp");
    save()=P-I;
    save.write("PPPfitBasicdiff.xmp");
    std::cout << "final score=" << result.score << ". Press" << std::endl;
    char c;
    std::cin >> c;
#endif
}
#undef DEBUG

/* Look for K neighbours in a list ----------------------------------------- */
//#define DEBUG
void CL3DClass::lookForNeighbours(const std::vector<CL3DClass *> listP, int K)
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
        CL3DAssignment assignment;
        MultidimArray<double> I;
        for (int q = 0; q < Q; q++)
        {
            if (listP[q] == this)
                distanceCode(q) = 1;
            else
            {
                I = listP[q]->P;
                fitBasic(I, assignment);
                A1D_ELEM(distanceCode,q) = 1-assignment.score;
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
void CL3D::shareAssignments(bool shareAssignment, bool shareUpdates)
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
    }

    // Share code updates
    if (shareUpdates)
    {
        std::vector<CL3DAssignment> auxList;
        std::vector<double> auxList2;
        int Q = P.size();
        for (int q = 0; q < Q; q++)
        {
            MPI_Allreduce(MPI_IN_PLACE, MULTIDIM_ARRAY(P[q]->Pupdate),
                          2*MULTIDIM_SIZE(P[q]->Pupdate), MPI_DOUBLE, MPI_SUM,
                          MPI_COMM_WORLD);

            // Share nextClassCorr and nextNonClassCorr
            std::vector<CL3DAssignment> receivedNextListImage;
            int listSize;
            for (size_t rank = 0; rank < prm->node->size; rank++)
            {
                if (rank == prm->node->rank)
                {
                    listSize = P[q]->nextListImg.size();
                    MPI_Bcast(&listSize, 1, MPI_INT, rank, MPI_COMM_WORLD);
                    MPI_Bcast(&(P[q]->nextListImg[0]),
                              P[q]->nextListImg.size() * sizeof(CL3DAssignment),
                              MPI_CHAR, rank, MPI_COMM_WORLD);
                }
                else
                {
                    MPI_Bcast(&listSize, 1, MPI_INT, rank, MPI_COMM_WORLD);
                    auxList.resize(listSize);
                    MPI_Bcast(&(auxList[0]), listSize * sizeof(CL3DAssignment),
                              MPI_CHAR, rank, MPI_COMM_WORLD);
                    receivedNextListImage.reserve(
                        receivedNextListImage.size() + listSize);
                    for (int j = 0; j < listSize; j++)
                        receivedNextListImage.push_back(auxList[j]);
                }
            }
            // Copy the received elements
            listSize = receivedNextListImage.size();
            P[q]->nextListImg.reserve(P[q]->nextListImg.size() + listSize);
            for (int j = 0; j < listSize; j++)
                P[q]->nextListImg.push_back(receivedNextListImage[j]);
        }

        transferUpdates();
    }
}

void CL3D::shareSplitAssignments(Matrix1D<int> &assignment, CL3DClass *node1,
                                 CL3DClass *node2) const
{
    // Share assignment
    int imax = VEC_XSIZE(assignment);
    MPI_Allreduce(MPI_IN_PLACE, MATRIX1D_ARRAY(assignment), imax, MPI_INT,
                  MPI_MAX, MPI_COMM_WORLD);

    // Share code updates
    std::vector<CL3DAssignment> auxList;
    std::vector<double> auxList2;
    CL3DClass *node = NULL;
    for (int q = 0; q < 2; q++)
    {
        if (q == 0)
            node = node1;
        else
            node = node2;
        MPI_Allreduce(MPI_IN_PLACE, MULTIDIM_ARRAY(node->Pupdate),
                      2*MULTIDIM_SIZE(node->Pupdate), MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);

        // Share nextClassCorr and nextNonClassCorr
        std::vector<CL3DAssignment> receivedNextListImage;
        int listSize;
        for (size_t rank = 0; rank < prm->node->size; rank++)
        {
            if (rank == prm->node->rank)
            {
                // Transmit node 1 next list of images
                listSize = node->nextListImg.size();
                MPI_Bcast(&listSize, 1, MPI_INT, rank, MPI_COMM_WORLD);
                MPI_Bcast(&(node->nextListImg[0]),
                          node->nextListImg.size() * sizeof(CL3DAssignment),
                          MPI_CHAR, rank, MPI_COMM_WORLD);
            }
            else
            {
                // Receive node 1 next list of images
                MPI_Bcast(&listSize, 1, MPI_INT, rank, MPI_COMM_WORLD);
                auxList.resize(listSize);
                MPI_Bcast(&(auxList[0]), listSize * sizeof(CL3DAssignment),
                          MPI_CHAR, rank, MPI_COMM_WORLD);
                receivedNextListImage.reserve(
                    receivedNextListImage.size() + listSize);
                for (int j = 0; j < listSize; j++)
                    receivedNextListImage.push_back(auxList[j]);
            }
        }
        // Copy the received elements
        listSize = receivedNextListImage.size();
        node->nextListImg.reserve(node->nextListImg.size() + listSize);
        for (int j = 0; j < listSize; j++)
            node->nextListImg.push_back(receivedNextListImage[j]);
    }

    node1->transferUpdate();
    node2->transferUpdate();
}

/* Read image --------------------------------------------------------- */
void CL3D::readImage(Image<double> &I, size_t objId, bool applyGeo) const
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
    I().statisticsAdjust(0, 1);
}

/* CL3D initialization ------------------------------------------------ */
//#define DEBUG
void CL3D::initialize(MetaData &_SF,
                      std::vector<MultidimArray<double> > &_codes0)
{
    if (prm->node->rank == 0)
        std::cout << "Initializing ...\n";

    SF = &_SF;
    Nimgs = SF->size();

    // Start with _Ncodes0 codevectors
    CL3DAssignment assignment;
    assignment.score=1;
    bool initialCodesGiven = _codes0.size() > 0;
    for (int q = 0; q < prm->Ncodes0; q++)
    {
        P.push_back(new CL3DClass());
        if (initialCodesGiven)
        {
        	P[q]->updateProjection(_codes0[q],assignment,true);
            P[q]->transferUpdate();
        }
    }

    // If no previous classes have been given, assign randomly
    if (prm->node->rank == 0)
        init_progress_bar(Nimgs);
    Image<double> I;
    MultidimArray<double> Iaux, Ibest;
    CL3DAssignment bestAssignment;
    size_t first, last;
    while (prm->taskDistributor->getTasks(first, last))
    {
        for (size_t idx = first; idx <= last; ++idx)
        {
            int q = -1;
            if (!initialCodesGiven)
                q = idx % (prm->Ncodes0);
            size_t objId = prm->objId[idx];
            readImage(I, objId, true);

            // Put it randomly in one of the classes
            bestAssignment.objId = assignment.objId = objId;
            if (!initialCodesGiven)
            {
                bestAssignment.score = 1;
                if (prm->randomizeStartingOrientation)
                {
					// Randomize the orientation of the volume, to avoid aligning all the missing wedges
					double rot=rnd_unif(0,360);
					double tilt=rnd_unif(0,360);
					double psi=rnd_unif(0,360);
					Matrix2D<double> E;
					Euler_angles2matrix(rot,tilt,psi,E,true);
					selfApplyGeometry(LINEAR,I(),E,IS_NOT_INV,DONT_WRAP,0.0);
                }
                P[q]->updateProjection(I(), bestAssignment);
            }
            else
            {
                bestAssignment.score = -1e38;
                q = -1;
                for (int qp = 0; qp < prm->Ncodes0; qp++)
                {
                    Iaux = I();
                    P[qp]->fitBasic(Iaux, assignment);
                    if (assignment.score > bestAssignment.score)
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

    // Share all assignments
    shareAssignments(true, true);
}
#undef DEBUG

/* CL3D write --------------------------------------------------------- */
void CL3D::write(const FileName &fnRoot, int level) const
{
    int Q = P.size();
    MetaData SFout;
    Image<double> I;
    FileName fnOut = formatString("%s_classes_level_%02d.stk",fnRoot.c_str(),level), fnClass;
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
    FileName fnSFout = formatString("%s_classes_level_%02d.xmd",fnRoot.c_str(),level);
    SFout.write(formatString("classes@%s", fnSFout.c_str()), MD_APPEND);

    // Make the selfiles of each class
    FileName fnImg;
    MDRow row;
    for (int q = 0; q < Q; q++)
    {
        MetaData SFq;
        std::vector<CL3DAssignment> &currentListImg = P[q]->currentListImg;
        int imax = currentListImg.size();
        for (int i = 0; i < imax; i++)
        {
            const CL3DAssignment &assignment = currentListImg[i];
            SF->getRow(row,assignment.objId);
            row.setValue(MDL_SHIFT_X, assignment.shiftx);
            row.setValue(MDL_SHIFT_Y, assignment.shifty);
            row.setValue(MDL_SHIFT_Z, assignment.shiftz);
            row.setValue(MDL_ANGLE_ROT, assignment.rot);
            row.setValue(MDL_ANGLE_TILT, assignment.tilt);
            row.setValue(MDL_ANGLE_PSI, assignment.psi);
            SFq.addRow(row);
        }
        MetaData SFq_sorted;
        SFq_sorted.sort(SFq, MDL_IMAGE);
        SFq_sorted.write(formatString("class%06d_images@%s",q+1, fnSFout.c_str()), MD_APPEND);
    }
}

void CL3D::lookNode(MultidimArray<double> &I, int oldnode, int &newnode,
                    CL3DAssignment &bestAssignment)
{
    int Q = P.size();
    int bestq = -1;
    MultidimArray<double> bestImg, Iaux;
    Matrix1D<double> corrList;
    corrList.resizeNoCopy(Q);
    CL3DAssignment assignment;
    bestAssignment.score = -1e38;
    size_t objId = bestAssignment.objId;
    for (int q = 0; q < Q; q++)
    {
        // Check if q is neighbour of the oldnode
        bool proceed = false;
        if (oldnode >= 0)
        {
            int imax = P[oldnode]->neighboursIdx.size();
            for (int i = 0; i < imax; i++)
                if (P[oldnode]->neighboursIdx[i] == q)
                {
                    proceed = true;
                    break;
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

        if (proceed)
        {
            // Try this image
            Iaux = I;
            P[q]->fitBasic(Iaux, assignment);

            VEC_ELEM(corrList,q) = assignment.score;
            if (assignment.score > bestAssignment.score ||
                prm->classifyAllImages)
            {
                bestq = q;
                bestImg = Iaux;
                bestAssignment = assignment;
            }
        }
    }

    I = bestImg;
    newnode = bestq;
    bestAssignment.objId = objId;

    // Assign it to the new node
    if (newnode != -1)
    	P[newnode]->updateProjection(I, bestAssignment);
}

void CL3D::transferUpdates()
{
    int Q = P.size();
    for (int q = 0; q < Q; q++)
        P[q]->transferUpdate();
}

/* Run CL3D ------------------------------------------------------------------ */
//#define DEBUG
void CL3D::run(const FileName &fnOut, int level)
{
    int Q = P.size();

    if (prm->node->rank == 0)
        std::cout << "Quantizing with " << Q << " codes...\n";

    std::vector<int> oldAssignment, newAssignment;

    int iter = 1;
    bool goOn = true;
    MetaData MDChanges;
    Image<double> I;
    int progressStep = XMIPP_MAX(1,Nimgs/60);
    CL3DAssignment assignment;
    while (goOn)
    {
        if (prm->node->rank == 0)
        {
            std::cout << "Iteration " << iter << " ...\n";
            std::cerr << "Iteration " << iter << " ...\n";
            init_progress_bar(Nimgs);
        }
        LOG(formatString("Iteration %d",iter));

        int K = XMIPP_MIN(prm->Nneighbours+1,Q);
        if (K == 0)
            K = Q;
        for (int q = 0; q < Q; q++)
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

                assignment.objId = objId;
                lookNode(I(), oldAssignment[idx], node, assignment);
                LOG(formatString("Analyzing %s oldAssignment=%d newAssignment=%d",I.name().c_str(),oldAssignment[idx], node));
                SF->setValue(MDL_REF, node + 1, objId);
                if (assignment.score>0)
                    corrSum += assignment.score;
                if (prm->node->rank == 0 && idx % progressStep == 0)
                    progress_bar(idx);
            }
        }
        FileName fnAux;
        for (int q=0; q<Q; q++)
        {
            for (size_t n=0; n<P[q]->currentListImg.size(); n++)
            {
                SF->getValue(MDL_IMAGE,fnAux,P[q]->currentListImg[n].objId);
                LOG(formatString("In node %d (%d): %s",q,n,fnAux.c_str()));
            }
        }

        // Gather all pieces computed by nodes
        MPI_Allreduce(MPI_IN_PLACE, &corrSum, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        shareAssignments(true, true);

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
        int Nchanges = 0;
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
            std::cout << "Number of assignment changes=" << Nchanges
            << std::endl;
            MDChanges.setValue(MDL_CL2D_CHANGES, Nchanges, idMdChanges);
            MDChanges.write(
                formatString("info@%s_classes_level_%02d.xmd",fnOut.c_str(), level));
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
                for (int q = 0; q < Q; q++)
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
                        std::cout << "Splitting node " << largestNode << " (" << sizeLargestNode
                        << ") by overwriting " << smallNode << " (" << sizeSmallestNode << " )" << std::endl;
                    smallNodes = true;

                    // Clear the old assignment of the images in the small node
                    std::vector<CL3DAssignment> &currentListImgSmall =
                        P[smallNode]->currentListImg;
                    size_t iimax = currentListImgSmall.size();
                    for (size_t ii = 0; ii < iimax; ii++)
                        SF->setValue(MDL_REF, -1, currentListImgSmall[ii].objId);
                    delete P[smallNode];

                    // Clear the old assignment of the images in the large node
                    std::vector<CL3DAssignment> &currentListImgLargest =
                        P[largestNode]->currentListImg;
                    iimax = currentListImgLargest.size();
                    for (size_t ii = 0; ii < iimax; ii++)
                        SF->setValue(MDL_REF, -1,
                                     currentListImgLargest[ii].objId);

                    // Now split the largest node
                    CL3DClass *node1 = new CL3DClass();
                    CL3DClass *node2 = new CL3DClass();
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
            write(fnOut,level);

        if ((iter > 1 && Nchanges < 0.005 * Nimgs && Q > 1) || iter >= prm->Niter)
            goOn = false;
        iter++;
    }

    std::sort(P.begin(), P.end(), SDescendingClusterSort());
}

/* Clean ------------------------------------------------------------------- */
int CL3D::cleanEmptyNodes()
{
    int retval = 0;
    std::vector<CL3DClass *>::iterator ptr = P.begin();
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
void CL3D::splitNode(CL3DClass *node, CL3DClass *&node1, CL3DClass *&node2,
                     std::vector<size_t> &splitAssignment,
                     bool iterate) const
{
    LOG(formatString("Splitting node 0: listsize=%d (volsize=%d)",node->currentListImg.size(),XSIZE(node->P)));
    std::vector<CL3DClass *> toDelete;
    Matrix1D<int> newAssignment, oldAssignment, firstSplitAssignment;
    Image<double> I;
    MultidimArray<double> Iaux1, Iaux2, corrList;
    MultidimArray<int> idx;
    CL3DAssignment assignment, assignment1, assignment2;
    CL3DClass *firstSplitNode1 = NULL;
    CL3DClass *firstSplitNode2 = NULL;
    size_t minAllowedSize = (size_t)(prm->PminSize * 0.01 * node->currentListImg.size());

    bool finish;
    bool success = true;
    do
    {
        finish = true;
        node2->neighboursIdx = node1->neighboursIdx = node->neighboursIdx;
        node1->P = node->P;
        node2->P = node->P;

        size_t imax = node->currentListImg.size();
        if (imax < minAllowedSize)
        {
            toDelete.push_back(node1);
            toDelete.push_back(node2);
            success = false;
            break;
        }

        // Compute the score histogram
        if (prm->node->rank == 0 && prm->verbose >= 2)
            std::cerr << "Calculating score distribution at split ..."
            << std::endl;
        corrList.initZeros(imax);
        for (size_t i = 0; i < imax; i++)
        {
            if ((i + 1) % (prm->node->size) == prm->node->rank)
            {
                readImage(I, node->currentListImg[i].objId, false);
                node->fitBasic(I(), assignment);
                A1D_ELEM(corrList,i) = assignment.score;
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
        LOG(formatString("Splitting node corrThreshold=%f",corrThreshold));
        if (corrThreshold == 0)
        {
            if (firstSplitNode1 != NULL)
            {
                toDelete.push_back(node1);
                toDelete.push_back(node2);
                success = false;
                break;
            }
            else
            {
                LOG(((String)"Splitting at random"));
                // Split at random
                for (size_t i = 0; i < imax; i++)
                {
                    assignment.objId = node->currentListImg[i].objId;
                    readImage(I, assignment.objId, false);
                    node->fitBasic(I(), assignment);
                    if ((i + 1) % 2 == 0)
                    {
                        node1->updateProjection(I(), assignment);
                        VEC_ELEM(newAssignment,i) = 1;
                    }
                    else
                    {
                        node2->updateProjection(I(), assignment);
                        VEC_ELEM(newAssignment,i) = 2;
                    }
                }
                success = true;
                break;
            }
        }

        // Split according to score
        if (prm->node->rank == 0 && prm->verbose >= 2)
            std::cerr << "Splitting by score threshold ..." << std::endl;
        LOG(((String)"Splitting by threshold"));
        for (size_t i = 0; i < imax; i++)
        {
            if ((i + 1) % (prm->node->size) == prm->node->rank)
            {
                assignment.objId = node->currentListImg[i].objId;
                readImage(I, assignment.objId, false);
                node->fitBasic(I(), assignment);
                if (assignment.score > corrThreshold)
                {
                    node1->updateProjection(I(), assignment);
                    VEC_ELEM(newAssignment,i) = 1;
                }
                else
                {
                    node2->updateProjection(I(), assignment);
                    VEC_ELEM(newAssignment,i) = 2;
                }
            }
            if (prm->node->rank == 0 && i % 25 == 0 && prm->verbose >= 2)
                progress_bar(i);
        }
        if (prm->node->rank == 0 && prm->verbose >= 2)
            progress_bar(imax);
        shareSplitAssignments(newAssignment, node1, node2);

        // Backup the first split in case it fails
        if (firstSplitNode1 == NULL)
        {
            firstSplitAssignment = newAssignment;
            firstSplitNode1 = new CL3DClass(*node1);
            firstSplitNode2 = new CL3DClass(*node2);
        }

        LOG(formatString("Splitting node 1: listsize1=%d (volsize1=%d)",node1->currentListImg.size(),XSIZE(node1->P)));
        LOG(formatString("Splitting node 1: listsize2=%d (volsize2=%d)",node2->currentListImg.size(),XSIZE(node2->P)));

        // Split iterations
        if (iterate)
        {
			for (int it = 0; it < prm->Niter; it++)
			{
				if (prm->node->rank == 0 && prm->verbose >= 2)
				{
					std::cerr << "Split iteration " << it << std::endl;
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
						node1->fitBasic(Iaux1, assignment1);
						Iaux2 = I();
						node2->fitBasic(Iaux2, assignment2);

						if (assignment1.score > assignment2.score && assignment1.score>0)
						{
							node1->updateProjection(Iaux1, assignment1);
							VEC_ELEM(newAssignment,i) = 1;
						}
						else if (assignment2.score > assignment1.score && assignment2.score>0)
						{
							node2->updateProjection(Iaux2, assignment2);
							VEC_ELEM(newAssignment,i) = 2;
						}
					}
					if (prm->node->rank == 0 && i % 25 == 0 && prm->verbose >= 2)
						progress_bar(i);
				}
				if (prm->node->rank == 0 && prm->verbose >= 2)
					progress_bar(imax);
				shareSplitAssignments(newAssignment, node1, node2);

				int Nchanges = 0;
				FOR_ALL_ELEMENTS_IN_MATRIX1D(newAssignment)
				if (newAssignment(i) != oldAssignment(i))
					Nchanges++;
				if (prm->node->rank == 0 && prm->verbose >= 2)
					std::cout << "Number of assignment split changes=" << Nchanges
					<< std::endl;

				// Check if one of the nodes is too small
				LOG(formatString("Splitting node 1.5: list1.size=%d list2.size=%d",node1->currentListImg.size(),node2->currentListImg.size()));
				if (node1->currentListImg.size() < minAllowedSize
					|| node2->currentListImg.size() < minAllowedSize
					|| Nchanges < 0.005 * imax)
					break;
			}
        }

        if (node1->currentListImg.size() < minAllowedSize)
        {
            if (prm->node->rank == 0 && prm->verbose >= 2)
                std::cout << "Removing node1, it's too small "
                << node1->currentListImg.size() << " "
                << minAllowedSize << "...\n";
            if (node1 != node)
                delete node1;
            node1 = new CL3DClass();
            toDelete.push_back(node2);
            node = node2;
            node2 = new CL3DClass();
            finish = false;
            LOG(formatString("Splitting node 2: Deleting small node 1: listsize1=%d (volsize1=%d)",node1->currentListImg.size(),XSIZE(node1->P)));
        }
        else if (node2->currentListImg.size() < minAllowedSize)
        {
            if (prm->node->rank == 0 && prm->verbose >= 2)
                std::cout << "Removing node2, it's too small "
                << node2->currentListImg.size() << " "
                << minAllowedSize << "...\n";
            if (node2 != node)
                delete node2;
            node2 = new CL3DClass();
            toDelete.push_back(node1);
            node = node1;
            node1 = new CL3DClass();
            finish = false;
            LOG(formatString("Splitting node 2: Deleting small node 2: listsize2=%d (volsize2=%d)",node2->currentListImg.size(),XSIZE(node2->P)));
        }
    }
    while (!finish);
    for (size_t i = 0; i < toDelete.size(); i++)
        if (toDelete[i] != node)
            delete toDelete[i];

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
}
#undef DEBUG

void CL3D::splitFirstNode()
{
    std::sort(P.begin(), P.end(), SDescendingClusterSort());
    int Q = P.size();
    P.push_back(new CL3DClass());
    P.push_back(new CL3DClass());
    std::vector<size_t> splitAssignment;
    splitNode(P[0], P[Q], P[Q + 1], splitAssignment, P.size()>1);
    delete P[0];
    P[0] = NULL;
    P.erase(P.begin());
}

/* MPI constructor --------------------------------------------------------- */
ProgClassifyCL3D::ProgClassifyCL3D(int argc, char** argv)
{
    node = new MpiNode(argc, argv);
    if (!node->isMaster())
        verbose = 0;
    taskDistributor = NULL;
    mask.allowed_data_types=INT_MASK;
}

/* Destructor -------------------------------------------------------------- */
ProgClassifyCL3D::~ProgClassifyCL3D()
{
    delete node;
    delete taskDistributor;
}

/* VQPrm I/O --------------------------------------------------------------- */
void ProgClassifyCL3D::readParams()
{
    fnSel = getParam("-i");
    fnOut = getParam("--oroot");
    fnCodes0 = getParam("--ref0");
    Niter = getIntParam("--iter");
    Nneighbours = getIntParam("--neigh");
    Ncodes0 = getIntParam("--nref0");
    Ncodes = getIntParam("--nref");
    PminSize = getDoubleParam("--minsize");
    sparsity = getDoubleParam("--sparsity");
    DWTsparsity = getDoubleParam("--DWTsparsity");
    maxShiftZ = getDoubleParam("--maxShiftZ");
    maxShiftY = getDoubleParam("--maxShiftY");
    maxShiftX = getDoubleParam("--maxShiftX");
    maxShift=std::max(maxShiftX,maxShiftY);
    maxShift=std::max(maxShift,maxShiftZ);
    maxFreq=getDoubleParam("--maxFreq");
    maxRot = getDoubleParam("--maxRot");
    maxTilt = getDoubleParam("--maxTilt");
    maxPsi = getDoubleParam("--maxPsi");
    classifyAllImages = checkParam("--classifyAllImages");
    randomizeStartingOrientation = checkParam("--randomizeStartingOrientation");
    fnSym = getParam("--sym");
    if (checkParam("--mask"))
        mask.readParams(this);
    dontAlign = checkParam("--dontAlign");
}

void ProgClassifyCL3D::show() const
{
    if (!verbose)
        return;
    std::cout << "Input images:            " << fnSel << std::endl
    << "Output images:           " << fnOut << std::endl
    << "Iterations:              " << Niter << std::endl;
    if (fnCodes0!="")
        std::cout << "CodesSel0:               " << fnCodes0 << std::endl;
    else
        std::cout << "Codes0:                  " << Ncodes0 << std::endl;
    std::cout << "Codes:                   " << Ncodes << std::endl
    << "Neighbours:              " << Nneighbours << std::endl
    << "Minimum node size:       " << PminSize << std::endl
    << "Sparsity:                " << sparsity << std::endl
    << "DWT Sparsity:            " << DWTsparsity << std::endl
    << "Maximum shift Z:         " << maxShiftZ << std::endl
    << "Maximum shift Y:         " << maxShiftY << std::endl
    << "Maximum shift X:         " << maxShiftX << std::endl
    << "Maximum rot:             " << maxRot << std::endl
    << "Maximum tilt:            " << maxTilt << std::endl
    << "Maximum psi:             " << maxPsi << std::endl
    << "Classify all images:     " << classifyAllImages << std::endl
    << "Symmetry:                " << fnSym << std::endl
    << "Don't align:             " << dontAlign << std::endl
    ;
    mask.show();
}

void ProgClassifyCL3D::defineParams()
{
    addUsageLine("Divide a selfile of volumes into the desired number of classes. ");
    addUsageLine("+Vector quantization with correntropy and a probabilistic criterion is used for creating the subdivisions.");
    addUsageLine("+Correlation and the standard maximum correlation criterion can also be used and normally produce good results.");
    addUsageLine("+Correntropy and the probabilistic clustering criterion are recommended for images with very low SNR or cases in which the correlation have clear difficulties to converge.");
    addUsageLine("+");
    addUsageLine("+The algorithm is fully described in [[http://www.ncbi.nlm.nih.gov/pubmed/20362059][this article]].");
    addUsageLine("+");
    addUsageLine("+An interesting convergence criterion is the number of images changing classes between iterations. If a low percentage of the image change class, then the clustering is rather stable and clear.");
    addUsageLine("+If many images change class, it is likely that there is not enough SNR to determine so many classes. It is recommended to reduce the number of classes");
    addParamsLine("    -i <selfile>             : Selfile with the input images");
    addParamsLine("   [--oroot <root=class>]    : Output rootname, by default, class");
    addParamsLine("   [--iter <N=20>]           : Number of iterations");
    addParamsLine("   [--nref0 <N=2>]           : Initial number of code vectors");
    addParamsLine("or  --ref0 <selfile=\"\">    : Selfile with initial code vectors");
    addParamsLine("   [--nref <N=16>]           : Final number of code vectors");
    addParamsLine("   [--neigh+ <N=4>]          : Number of neighbour code vectors");
    addParamsLine("                             : Set -1 for all");
    addParamsLine("   [--minsize+ <N=20>]       : Percentage minimum node size");
    addParamsLine("   [--sparsity+ <f=0.975>]   : Percentage of Fourier coefficients to drop");
    addParamsLine("   [--DWTsparsity+ <f=0.99>] : Percentage of wavelet coefficients to drop");
    addParamsLine("   [--maxShiftX <d=10>]      : Maximum allowed shift in X");
    addParamsLine("   [--maxShiftY <d=10>]      : Maximum allowed shift in Y");
    addParamsLine("   [--maxShiftZ <d=10>]      : Maximum allowed shift in Z");
    addParamsLine("   [--maxRot <d=360>]        : Maximum allowed rotational angle in absolute value");
    addParamsLine("   [--maxTilt <d=360>]       : Maximum allowed tilt angle in absolute value");
    addParamsLine("   [--maxPsi <d=360>]        : Maximum allowed in-plane angle in absolute value");
    addParamsLine("   [--classifyAllImages]     : By default, some images may not be classified. Use this option to classify them all.");
    addParamsLine("   [--sym <s=c1>]            : Symmetry of the classes to be reconstructed");
    addParamsLine("   [--maxFreq <w=0.2>]       : Maximum frequency to be reconstructed");
    addParamsLine("   [--randomizeStartingOrientation] : Use this option to avoid aligning all missing wedges");
    addParamsLine("   [--dontAlign]             : Do not align volumes, only classify");
    Mask::defineParams(this,INT_MASK,NULL,NULL,true);
    addExampleLine("mpirun -np 3 `which xmipp_mpi_classify_CL3D` -i images.stk --nref 256 --oroot class --iter 10");
}

void ProgClassifyCL3D::produceSideInfo()
{
    gaussianInterpolator.initialize(6, 60000, false);

    // Get image dimensions
    SF.read(fnSel);
    size_t Ndim;
    getImageSize(SF, Xdim, Ydim, Zdim, Ndim);

    // MaxFreq mask
    MultidimArray<double> V;
    MultidimArray< std::complex<double> > VF;
    V.initZeros(Zdim,Ydim,Xdim);
    FourierTransformer transformer;
    transformer.FourierTransform(V,VF);
    maxFreqMask.initZeros(VF);
    Matrix1D<int> idx(3);
    Matrix1D<double> freq(3);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(maxFreqMask)
    {
    	XX(idx)=j;
    	YY(idx)=i;
    	ZZ(idx)=k;
    	FFT_idx2digfreq(V,idx,freq);
    	if (freq.module()<=maxFreq)
    		A3D_ELEM(maxFreqMask,k,i,j)=1;
    }

    // Prepare symmetry list
    SL.readSymmetryFile(fnSym);

    // Prepare the Task distributor
    SF.findObjects(objId);
    size_t Nimgs = objId.size();
    taskDistributor = new FileTaskDistributor(Nimgs,
                      XMIPP_MAX(1,Nimgs/(5*node->size)), node);

    // Prepare mask for evaluating the noise outside
    MultidimArray<int> sphericalMask;
    sphericalMask.resize(Zdim, Ydim, Xdim);
    sphericalMask.setXmippOrigin();
    BinaryCircularMask(sphericalMask, Xdim / 2, INNER_MASK);

    mask.generate_mask(Zdim, Ydim, Xdim);
    MultidimArray<int> &mMask=mask.get_binary_mask();
    mMask.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY3D(sphericalMask)
    if (A3D_ELEM(sphericalMask,k,i,j)==0)
        A3D_ELEM(mMask,k,i,j)=0;

    // Read input code vectors if available
    std::vector<MultidimArray<double> > codes0;
    if (fnCodes0 != "")
    {
        Image<double> I;
        MetaData SFCodes(fnCodes0);

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

void ProgClassifyCL3D::run()
{
    String xmippPython;
	initializeXmippPython(xmippPython);
    frmFunc = getPointerToPythonFRMFunction();
    wedgeClass = getPointerToPythonGeneralWedgeClass();

    CREATE_LOG();
    show();
    produceSideInfo();

    // Run all iterations
    int level = 0;
    int Q = vq.P.size();
    vq.run(fnOut, level);

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
        vq.run(fnOut, level);
    }
    if (node->rank == 0)
    {
        std::sort(vq.P.begin(), vq.P.end(), SDescendingClusterSort());
        Q = vq.P.size();
        MetaData SFq, SFclassified, SFaux, SFaux2;
        for (int q = 0; q < Q; q++)
        {
            SFq.read(formatString("class%06d_images@%s_classes_level_%02d.xmd", q + 1,
                                  fnOut.c_str(), level));
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
        SFaux.write(fnOut + "_images.xmd");
    }
    CLOSE_LOG();
}

/* Main -------------------------------------------------------------------- */
int main(int argc, char** argv)
{
    ProgClassifyCL3D progprm(argc, argv);
    progprm.read(argc, argv);
    prm = &progprm;
    return progprm.tryRun();
}
