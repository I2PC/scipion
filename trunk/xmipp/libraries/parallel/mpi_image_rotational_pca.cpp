/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano      coss@cnb.csic.es (2002)
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

// Translated from MATLAB code by Yoel Shkolnisky

#include "mpi_image_rotational_pca.h"
#include <data/mask.h>
#include <data/metadata_extension.h>
#include <mpi.h>

// Empty constructor =======================================================
ProgImageRotationalPCA::ProgImageRotationalPCA(int argc, char **argv)
{
    node=new MpiNode(argc,argv);
    if (!node->isMaster())
        verbose=0;
}

// MPI destructor
ProgImageRotationalPCA::~ProgImageRotationalPCA()
{
	delete fileMutex;
    delete node;
    for (int n=0; n<HbufferMax; n++)
    	delete [] Hbuffer[n];
}

// Read arguments ==========================================================
void ProgImageRotationalPCA::readParams()
{
    fnIn = getParam("-i");
    fnRoot = getParam("--oroot");
    Neigen = getIntParam("--eigenvectors");
    Nits = getIntParam("--iterations");
    max_shift_change = getDoubleParam("--max_shift_change");
    psi_step = getDoubleParam("--psi_step");
    shift_step = getDoubleParam("--shift_step");
}

// Show ====================================================================
void ProgImageRotationalPCA::show()
{
    if (!verbose)
        return;
    std::cout
    << "Input:               " << fnIn << std::endl
    << "Output root:         " << fnRoot << std::endl
    << "Number eigenvectors: " << Neigen << std::endl
    << "Number iterations:   " << Nits << std::endl
    << "Psi step:            " << psi_step << std::endl
    << "Max shift change:    " << max_shift_change << " step: " << shift_step << std::endl
    ;
}

// usage ===================================================================
void ProgImageRotationalPCA::defineParams()
{
    addUsageLine("Makes a rotational invariant representation of the image collection");
    addParamsLine("    -i <selfile>               : Selfile with experimental images");
    addParamsLine("   --oroot <rootname>          : Rootname for output");
    addParamsLine("  [--eigenvectors <N=200>]     : Number of eigenvectors");
    addParamsLine("  [--iterations <N=2>]         : Number of iterations");
    addParamsLine("  [--max_shift_change <r=0>]   : Maximum change allowed in shift");
    addParamsLine("  [--psi_step <ang=1>]         : Step in psi in degrees");
    addParamsLine("  [--shift_step <r=1>]         : Step in shift in pixels");
    addExampleLine("Typical use:",false);
    addExampleLine("xmipp_mpi_image_rotational_pca -i images.stk --oroot images_eigen");
}

// Produce side info =====================================================
void ProgImageRotationalPCA::produceSideInfo()
{
	time_config();

    MD.read(fnIn);
    Nimg=MD.size();
    int Ydim, Zdim;
    size_t Ndim;
    ImgSize(MD, Xdim, Ydim, Zdim, Ndim);
    Nangles=floor(360.0/psi_step);
    Nshifts=(2*max_shift_change+1)/shift_step;
    Nshifts*=Nshifts;

    // Construct mask
    mask.resizeNoCopy(Xdim,Xdim);
    mask.setXmippOrigin();
    double R2=0.25*Xdim*Xdim;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
    	A2D_ELEM(mask,i,j)=(i*i+j*j<R2);
    mask.initConstant(1); // COSS Temporarily the mask is full
    Npixels=(int)mask.sum();

    W.resizeNoCopy(Npixels, Neigen+2);

    if (node->isMaster())
    {
        // F (#images*#shifts*#angles) x (#eigenvectors+2)*(its+1)
        createEmptyFileWithGivenLength(fnRoot+"_matrixF.raw",Nimg*Nshifts*Nangles*(Neigen+2)*(Nits+1)*sizeof(double));
        // H (#images*#shifts*#angles) x (#eigenvectors+2)
        createEmptyFileWithGivenLength(fnRoot+"_matrixH.raw",Nimg*Nshifts*Nangles*(Neigen+2)*sizeof(double));

        // Initialize with random numbers between -1 and 1
        FOR_ALL_ELEMENTS_IN_MATRIX2D(W)
        MAT_ELEM(W,i,j)=rnd_unif(-1.0,1.0);

        // Send to workers
        MPI_Bcast(&MAT_ELEM(W,0,0),MAT_XSIZE(W)*MAT_YSIZE(W),MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    else
        // Receive W
        MPI_Bcast(&MAT_ELEM(W,0,0),MAT_XSIZE(W)*MAT_YSIZE(W),MPI_DOUBLE,0,MPI_COMM_WORLD);
    F.mapToFile(fnRoot+"_matrixF.raw",(Neigen+2)*(Nits+1),Nimg*Nshifts*Nangles);
    H.mapToFile(fnRoot+"_matrixH.raw",Nimg*Nshifts*Nangles,Neigen+2);
    Hblock.resizeNoCopy(Nangles*Nshifts,Neigen+2);

    // Prepare buffer
    fileMutex = new MpiFileMutex(node);
    for (int n=0; n<HbufferMax; n++)
    	Hbuffer.push_back(new double[MAT_XSIZE(Hblock)*MAT_YSIZE(Hblock)]);
}

// Buffer =================================================================
void ProgImageRotationalPCA::writeToHBuffer(double *dest)
{
	int n=HbufferDestination.size();
    memcpy(Hbuffer[n],&MAT_ELEM(Hblock,0,0),MAT_XSIZE(Hblock)*MAT_YSIZE(Hblock)*sizeof(double));
    HbufferDestination.push_back(dest);
    if (n==(HbufferMax-1))
    	flushHBuffer();
}

void ProgImageRotationalPCA::flushHBuffer()
{
	int nmax=HbufferDestination.size();
	fileMutex->lock();
	for (int n=0; n<nmax; ++n)
		memcpy(HbufferDestination[n],Hbuffer[n],MAT_XSIZE(Hblock)*MAT_YSIZE(Hblock)*sizeof(double));
	fileMutex->unlock();
	HbufferDestination.clear();
}

// Apply T ================================================================
void ProgImageRotationalPCA::applyT()
{
	TimeStamp t0;
	annotate_time(&t0);
    W.initZeros(Npixels,MAT_XSIZE(H));
    Wnode.initZeros(Npixels,MAT_XSIZE(H));

    FileName fnImg;
    size_t idx=0;
    const int unroll=8;
    const int jmax=(MAT_XSIZE(Wnode)/unroll)*unroll;

    FOR_ALL_OBJECTS_IN_METADATA(MD)
    {
        // Process image only if it is in my rank
        if ((++idx)%(node->size)!=node->rank)
            continue;

        // Read image
        MD.getValue(MDL_IMAGE,fnImg,__iter.objId);
        I.read(fnImg);
        const MultidimArray<double> &mI=I();

        // Locate the corresponding index in Matrix H
        // and copy a block in memory to speed up calculations
        size_t Hidx=(idx-1)*Nangles*Nshifts;
        memcpy(&MAT_ELEM(Hblock,0,0),&MAT_ELEM(H,Hidx,0),MAT_XSIZE(Hblock)*MAT_YSIZE(Hblock)*sizeof(double));

        // For each rotation and shift
        int block_idx=0;
        for (double psi=0; psi<360; psi+=psi_step)
        {
            rotation2DMatrix(psi,A,true);
            for (double y=-max_shift_change; y<=max_shift_change; y+=shift_step)
            {
                MAT_ELEM(A,1,2)=y;
                for (double x=-max_shift_change; x<=max_shift_change; x+=shift_step, ++block_idx)
                {
                    MAT_ELEM(A,0,2)=x;

                    // Rotate and shift image
                    applyGeometry(1,Iaux,mI,A,IS_INV,true);

                    // Update Wnode
                    for (int i=0; i<MAT_YSIZE(Wnode); ++i)
                    {
                        double pixval=DIRECT_MULTIDIM_ELEM(Iaux,i);
                        double *ptrWnode=&MAT_ELEM(Wnode,i,0);
                        unsigned char *ptrMask=&DIRECT_MULTIDIM_ELEM(mask,0);
                        double *ptrHblock=&MAT_ELEM(Hblock,block_idx,0);
                        for (int j=0; j<jmax; j+=unroll, ptrHblock+=unroll, ptrMask+=unroll)
                        {
                            if (*ptrMask    ) (*ptrWnode++) +=pixval*(*ptrHblock  );
                            if (*(ptrMask+1)) (*ptrWnode++) +=pixval*(*(ptrHblock+1));
                            if (*(ptrMask+2)) (*ptrWnode++) +=pixval*(*(ptrHblock+2));
                            if (*(ptrMask+3)) (*ptrWnode++) +=pixval*(*(ptrHblock+3));
                            if (*(ptrMask+4)) (*ptrWnode++) +=pixval*(*(ptrHblock+4));
                            if (*(ptrMask+5)) (*ptrWnode++) +=pixval*(*(ptrHblock+5));
                            if (*(ptrMask+6)) (*ptrWnode++) +=pixval*(*(ptrHblock+6));
                            if (*(ptrMask+7)) (*ptrWnode++) +=pixval*(*(ptrHblock+7));
                        }
                        for (int j=jmax; j<MAT_XSIZE(Wnode); ++j, ptrHblock+=1, ptrMask+=1)
                            if (*ptrMask) (*ptrWnode++) +=pixval*(*ptrHblock);
                    }
                }
            }
        }
    }
    MPI_Allreduce(MATRIX2D_ARRAY(Wnode), MATRIX2D_ARRAY(W), MAT_XSIZE(W)*MAT_YSIZE(W),
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    std::cout << "T:"; print_elapsed_time(t0);
}

// Apply T ================================================================
void ProgImageRotationalPCA::applyTt()
{
	TimeStamp t0;
	annotate_time(&t0);

	// Compute W transpose to accelerate memory access
    Wtranspose.resizeNoCopy(MAT_XSIZE(W),MAT_YSIZE(W));
    FOR_ALL_ELEMENTS_IN_MATRIX2D(Wtranspose)
    MAT_ELEM(Wtranspose,i,j) = MAT_ELEM(W,j,i);

    const size_t unroll=8;
    const size_t nmax=(MAT_XSIZE(Wtranspose)/unroll)*unroll;

    FileName fnImg;
    size_t idx=0;
    FOR_ALL_OBJECTS_IN_METADATA(MD)
    {
        // Process image only if it is in my rank
        if ((++idx)%(node->size)!=node->rank)
            continue;

        // Read image
        MD.getValue(MDL_IMAGE,fnImg,__iter.objId);
        I.read(fnImg);
        const MultidimArray<double> &mI=I();

        // For each rotation and shift
        int block_idx=0;
        for (double psi=0; psi<360; psi+=psi_step)
        {
            rotation2DMatrix(psi,A,true);
            for (double y=-max_shift_change; y<=max_shift_change; y+=shift_step)
            {
                MAT_ELEM(A,1,2)=y;
                for (double x=-max_shift_change; x<=max_shift_change; x+=shift_step, ++block_idx)
                {
                    MAT_ELEM(A,0,2)=x;

                    // Rotate and shift image
                    applyGeometry(1,Iaux,mI,A,IS_INV,true);

                    // Update Hblock
                    for (int j=0; j<MAT_XSIZE(Hblock); j++)
                    {
                        double dotproduct=0;
                        const double *ptrIaux=MULTIDIM_ARRAY(Iaux);
                        unsigned char *ptrMask=&DIRECT_MULTIDIM_ELEM(mask,0);
                        const double *ptrWtranspose=&MAT_ELEM(Wtranspose,j,0);
                        for (size_t n=0; n<nmax; n+=unroll, ptrIaux+=unroll, ptrMask+=unroll)
                        {
                            if (*(ptrMask  )) dotproduct+=(*(ptrIaux  ))*(*ptrWtranspose++);
                            if (*(ptrMask+1)) dotproduct+=(*(ptrIaux+1))*(*ptrWtranspose++);
                            if (*(ptrMask+2)) dotproduct+=(*(ptrIaux+2))*(*ptrWtranspose++);
                            if (*(ptrMask+3)) dotproduct+=(*(ptrIaux+3))*(*ptrWtranspose++);
                            if (*(ptrMask+4)) dotproduct+=(*(ptrIaux+4))*(*ptrWtranspose++);
                            if (*(ptrMask+5)) dotproduct+=(*(ptrIaux+5))*(*ptrWtranspose++);
                            if (*(ptrMask+6)) dotproduct+=(*(ptrIaux+6))*(*ptrWtranspose++);
                            if (*(ptrMask+7)) dotproduct+=(*(ptrIaux+7))*(*ptrWtranspose++);
                        }
                        for (size_t n=nmax; n<MAT_XSIZE(Wtranspose); ++n, ++ptrMask, ++ptrIaux)
                            if (*ptrMask) dotproduct+=(*ptrIaux)*(*ptrWtranspose++);
                        MAT_ELEM(Hblock,block_idx,j)=dotproduct;
                    }
                }
            }
        }

        // Locate the corresponding index in Matrix H
        // and copy block to disk
        size_t Hidx=(idx-1)*Nangles*Nshifts;
        writeToHBuffer(&MAT_ELEM(H,Hidx,0));
    }
    flushHBuffer();
    std::cout << "Tt:"; print_elapsed_time(t0);
}

// QR =====================================================================
int ProgImageRotationalPCA::QR()
{
	TimeStamp t0;
	annotate_time(&t0);
    size_t jQ=0;
    Matrix1D<double> qj1, qj2;
    int iBlockMax=MAT_XSIZE(F)/4;

    for (int j1=0; j1<MAT_YSIZE(F); j1++)
    {
        F.getRow(j1,qj1);
        // Project twice in the already established subspace
        // One projection should be enough but Gram-Schmidt suffers
        // from numerical problems
        for (int it=0; it<2; it++)
        {
            for (int j2=0; j2<jQ; j2++)
            {
                F.getRow(j2,qj2);

                // Compute dot product
                double s12=qj1.dotProduct(qj2);

                // Subtract the part of qj2 from qj1
                double *ptr1=&VEC_ELEM(qj1,0);
                const double *ptr2=&VEC_ELEM(qj2,0);
                for (int i=0; i<iBlockMax; i++)
                {
                    (*ptr1++)-=s12*(*ptr2++);
                    (*ptr1++)-=s12*(*ptr2++);
                    (*ptr1++)-=s12*(*ptr2++);
                    (*ptr1++)-=s12*(*ptr2++);
                }
                for (int i=iBlockMax*4; i<MAT_XSIZE(F); ++i)
                    (*ptr1++)-=s12*(*ptr2++);
            }
        }

        // Keep qj1 in Q if it has enough norm
        double Rii=qj1.module();
        if (Rii>1e-14)
        {
            // Make qj1 to be unitary and store in Q
            qj1/=Rii;
            F.setRow(jQ++,qj1);
        }
    }
    std::cout << "QR:"; print_elapsed_time(t0);
    return jQ;
}

// Copy H to F ============================================================
void ProgImageRotationalPCA::copyHtoF(int block)
{
	if (node->isMaster())
	{
		size_t Hidx=block*MAT_XSIZE(H);
		FOR_ALL_ELEMENTS_IN_MATRIX2D(H)
			MAT_ELEM(F,Hidx+j,i)=MAT_ELEM(H,i,j);
	}
	node->barrierWait();
}

// Run ====================================================================
void ProgImageRotationalPCA::run()
{
    show();
    produceSideInfo();

    // Compute matrix F:
    // Set H pointing to the first block of F
    applyTt(); // H=Tt(W)
    copyHtoF(0);
    for (int it=0; it<Nits; it++)
    {
        applyT(); // W=T(H)
        applyTt(); // H=Tt(W)
        copyHtoF(it+1);
    }
    H.clear();

    // QR decomposition of matrix F
    int qrDim;
    if (node->isMaster())
    {
        qrDim=QR();
        MPI_Bcast(&qrDim,1,MPI_INT,0,MPI_COMM_WORLD);
    }
    else
        MPI_Bcast(&qrDim,1,MPI_INT,0,MPI_COMM_WORLD);

    // Load the first qrDim columns of F in matrix H
	if (node->isMaster())
	{
        createEmptyFileWithGivenLength(fnRoot+"_matrixH.raw",Nimg*Nshifts*Nangles*qrDim*sizeof(double));
        H.mapToFile(fnRoot+"_matrixH.raw",MAT_XSIZE(F),qrDim);
		FOR_ALL_ELEMENTS_IN_MATRIX2D(H)
        MAT_ELEM(H,i,j)=MAT_ELEM(F,j,i);
		node->barrierWait();
	}
	else
	{
		node->barrierWait();
        H.mapToFile(fnRoot+"_matrixH.raw",MAT_XSIZE(F),qrDim);
	}

	// Apply T
    Hblock.resizeNoCopy(Nangles*Nshifts,qrDim);
    applyT();

    // Apply SVD and extract the basis
    if (node->isMaster())
    {
    	TimeStamp t0;
    	annotate_time(&t0);
        // SVD of W
        Matrix2D<double> U,V;
        Matrix1D<double> S;
        svdcmp(W,U,S,V);
        std::cout << "SVD:"; print_elapsed_time(t0);

        // Keep the first Neigen images from U
    	annotate_time(&t0);
        Image<double> I;
        I().resizeNoCopy(Xdim,Xdim);
        const MultidimArray<double> &mI=I();
        FileName fnImg;
        MetaData MD;
        for (int eig=0; eig<Neigen; eig++)
        {
        	int Un=0;
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mI)
            	if (DIRECT_MULTIDIM_ELEM(mask,n))
            		DIRECT_MULTIDIM_ELEM(mI,n)=MAT_ELEM(U,Un++,eig);
            fnImg.compose(eig+1,fnRoot,"stk");
            I.write(fnImg);
            size_t id=MD.addObject();
            MD.setValue(MDL_IMAGE,fnImg,id);
            MD.setValue(MDL_WEIGHT,VEC_ELEM(S,eig),id);
        }
        MD.write(fnRoot+".xmd");
        std::cout << "Output:"; print_elapsed_time(t0);
    }
}
