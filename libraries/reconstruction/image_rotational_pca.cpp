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
#include "image_rotational_pca.h"
#include <data/mask.h>
#include <data/metadata_extension.h>

// Empty constructor =======================================================
ProgImageRotationalPCA::ProgImageRotationalPCA()
{
  rank = 0;
  verbose = 1;
  fileMutex = NULL;
  threadMutex = NULL;
  taskDistributor = NULL;
  thMgr = NULL;
}

// MPI destructor
ProgImageRotationalPCA::~ProgImageRotationalPCA()
{
  clearHbuffer();
  delete fileMutex;
  delete threadMutex;
  delete taskDistributor;
  delete thMgr;
}

void
ProgImageRotationalPCA::clearHbuffer()
{
  int nmax = Hbuffer.size();
  for (int n = 0; n < nmax; n++)
    delete[] Hbuffer[n];
  Hbuffer.clear();
}

// Read arguments ==========================================================
void
ProgImageRotationalPCA::readParams()
{
  fnIn = getParam("-i");
  fnRoot = getParam("--oroot");
  Neigen = getIntParam("--eigenvectors");
  Nits = getIntParam("--iterations");
  max_shift_change = getDoubleParam("--max_shift_change");
  psi_step = getDoubleParam("--psi_step");
  shift_step = getDoubleParam("--shift_step");
  maxNimgs = getIntParam("--maxImages");
  Nthreads = getIntParam("--thr");
}

// Show ====================================================================
void
ProgImageRotationalPCA::show()
{
  if (!verbose)
    return;
  std::cout << "Input:               " << fnIn << std::endl
      << "Output root:         " << fnRoot << std::endl
      << "Number eigenvectors: " << Neigen << std::endl
      << "Number iterations:   " << Nits << std::endl << "Psi step:            "
      << psi_step << std::endl << "Max shift change:    " << max_shift_change
      << " step: " << shift_step << std::endl << "Max images:          "
      << maxNimgs << std::endl << "Number of threads:   " << Nthreads
      << std::endl;
}

// usage ===================================================================
void
ProgImageRotationalPCA::defineParams()
{
  addUsageLine(
      "Makes a rotational invariant representation of the image collection");
  addParamsLine(
      "    -i <selfile>               : Selfile with experimental images");
  addParamsLine("   --oroot <rootname>          : Rootname for output");
  addParamsLine("  [--eigenvectors <N=200>]     : Number of eigenvectors");
  addParamsLine("  [--iterations <N=2>]         : Number of iterations");
  addParamsLine(
      "  [--max_shift_change <r=0>]   : Maximum change allowed in shift");
  addParamsLine("  [--psi_step <ang=1>]         : Step in psi in degrees");
  addParamsLine("  [--shift_step <r=1>]         : Step in shift in pixels");
  addParamsLine("  [--maxImages <N=-1>]         : Maximum number of images");
  addParamsLine("  [--thr <N=1>]                : Number of threads");
  addExampleLine("Typical use (4 nodes with 4 processors):", false);
  addExampleLine(
      "mpirun -np 4 `which xmipp_mpi_image_rotational_pca` -i images.stk --oroot images_eigen --thr 4");
}

void ProgImageRotationalPCA::selectPartFromMd(MetaData &MDin)
{
    MetaData MDaux;
    MDaux.randomize(MDin);
    MDin.selectPart(MDaux, 0, maxNimgs);
}
void ProgImageRotationalPCA::comunicateMatrix(Matrix2D<double> &W)
{
}

void ProgImageRotationalPCA::createMutexes(size_t Nimgs)
{
  fileMutex = new Mutex();
  threadMutex = new Mutex();
  taskDistributor = new ThreadTaskDistributor(Nimgs, XMIPP_MAX(1,Nimgs/5));
}

// Produce side info =====================================================
void ProgImageRotationalPCA::produceSideInfo()
{
  time_config();
  MetaData MDin(fnIn);

  if (maxNimgs > 0)
    selectPartFromMd(MDin);

    Nimg = MDin.size();
    size_t Ydim, Zdim, Ndim;
    getImageSize(MDin, Xdim, Ydim, Zdim, Ndim);
    Nangles = (int)floor(360.0 / psi_step);
    Nshifts = (int)((2 * max_shift_change + 1) / shift_step);
    Nshifts *= Nshifts;

    // Construct mask
    mask.resizeNoCopy(Xdim, Xdim);
    mask.setXmippOrigin();
    double R2 = 0.25 * Xdim * Xdim;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
      A2D_ELEM(mask,i,j)=(i*i+j*j<R2);
    Npixels = (int) mask.sum();

    // Thread Manager
    thMgr = new ThreadManager(Nthreads, this);
    Image<double> dummy;
    Matrix2D<double> dummyMatrix, dummyHblock, dummyW;
    dummyW.resizeNoCopy(Npixels, Neigen + 2);
    dummyHblock.resizeNoCopy(2 * Nangles * Nshifts, Neigen + 2);
    for (int n = 0; n < Nthreads; ++n)
    {
      Wnode.push_back(dummyW);
      Hblock.push_back(dummyHblock);
      A.push_back(dummyMatrix);
      I.push_back(dummy);
      Iaux.push_back(dummy());
      MD.push_back(MDin);
    }

    Matrix2D<double> &W = Wnode[0];
    FileName fnMatrixF(fnRoot + "_matrixF.raw");
    FileName fnMatrixH(fnRoot + "_matrixH.raw");

    if (IS_MASTER)
    {
      // F (#images*#shifts*#angles) x (#eigenvectors+2)*(its+1)
      fnMatrixF.createEmptyFileWithGivenLength(
          Nimg * 2 * Nangles * Nshifts * (Neigen + 2) * (Nits + 1)
              * sizeof(double));
      // H (#images*#shifts*#angles) x (#eigenvectors+2)
      fnMatrixH.createEmptyFileWithGivenLength(
          Nimg * 2 * Nangles * Nshifts * (Neigen + 2) * sizeof(double));

      // Initialize with random numbers between -1 and 1
      FOR_ALL_ELEMENTS_IN_MATRIX2D(W)
        MAT_ELEM(W,i,j)=rnd_unif(-1.0,1.0);
    }

    comunicateMatrix(W);

    F.mapToFile(fnMatrixF, (Neigen + 2) * (Nits + 1),
        Nimg * 2 * Nangles * Nshifts);
    H.mapToFile(fnMatrixH, Nimg * 2 * Nangles * Nshifts, Neigen + 2);

    // Prepare buffer
    for (int n = 0; n < HbufferMax; n++)
      Hbuffer.push_back(
          new double[MAT_XSIZE(dummyHblock) * MAT_YSIZE(dummyHblock)]);

    // Construct a FileTaskDistributor
    MDin.findObjects(objId);
    size_t Nimgs = objId.size();
    createMutexes(Nimgs);
}

// Buffer =================================================================
void ProgImageRotationalPCA::writeToHBuffer(int idx, double *dest)
{
  threadMutex->lock();
  int n=HbufferDestination.size();
  const Matrix2D<double> &Hblock_idx=Hblock[idx];
  // COSS std::cout << "Copying block " << idx << " into " << n << "(" << Hbuffer.size() << ")" << std::endl;
  memcpy(Hbuffer[n],&MAT_ELEM(Hblock_idx,0,0),MAT_XSIZE(Hblock_idx)*MAT_YSIZE(Hblock_idx)*sizeof(double));
  HbufferDestination.push_back(dest);
  if (n==(HbufferMax-1))
  flushHBuffer();
  threadMutex->unlock();
}

void ProgImageRotationalPCA::flushHBuffer()
{
  int nmax=HbufferDestination.size();
  const Matrix2D<double> &Hblock_0=Hblock[0];
  fileMutex->lock();
  for (int n=0; n<nmax; ++n)
  memcpy(HbufferDestination[n],Hbuffer[n],MAT_XSIZE(Hblock_0)*MAT_YSIZE(Hblock_0)*sizeof(double));
  fileMutex->unlock();
  HbufferDestination.clear();
}

// Apply T ================================================================
void threadApplyT(ThreadArgument &thArg)
{
  ProgImageRotationalPCA *self=(ProgImageRotationalPCA *) thArg.workClass;
  //MpiNode *node=self->node;
  int rank = self->rank;
  ThreadTaskDistributor *taskDistributor=self->taskDistributor;
  std::vector<size_t> &objId=self->objId;
  MetaData &MD=self->MD[thArg.thread_id];

  Image<double> &I=self->I[thArg.thread_id];
  MultidimArray<double> &Iaux=self->Iaux[thArg.thread_id];
  Matrix2D<double> &Wnode=self->Wnode[thArg.thread_id];
  Matrix2D<double> &Hblock=self->Hblock[thArg.thread_id];
  Matrix2D<double> &A=self->A[thArg.thread_id];
  Matrix2D<double> &H=self->H;
  MultidimArray< unsigned char > &mask=self->mask;
  Wnode.initZeros(self->Npixels,MAT_XSIZE(H));

  const int unroll=8;
  const int jmax=(MAT_XSIZE(Wnode)/unroll)*unroll;

  size_t first, last;
  if (IS_MASTER && thArg.thread_id==0)
  {
    std::cerr << "Applying T ...\n";
    init_progress_bar(objId.size());
  }
  while (taskDistributor->getTasks(first, last))
  {
    for (size_t idx=first; idx<=last; ++idx)
    {
      // Read image
      I.readApplyGeo(MD,objId[idx]);
      MultidimArray<double> &mI=I();

      // Locate the corresponding index in Matrix H
      // and copy a block in memory to speed up calculations
      size_t Hidx=idx*2*self->Nangles*self->Nshifts;
      memcpy(&MAT_ELEM(Hblock,0,0),&MAT_ELEM(H,Hidx,0),MAT_XSIZE(Hblock)*MAT_YSIZE(Hblock)*sizeof(double));

      // For each rotation, shift and mirror
      int block_idx=0;
      for (int mirror=0; mirror<2; ++mirror)
      {
        if (mirror)
        {
          mI.selfReverseX();
          mI.setXmippOrigin();
        }
        for (double psi=0; psi<360; psi+=self->psi_step)
        {
          rotation2DMatrix(psi,A,true);
          for (double y=-self->max_shift_change; y<=self->max_shift_change; y+=self->shift_step)
          {
            MAT_ELEM(A,1,2)=y;
            for (double x=-self->max_shift_change; x<=self->max_shift_change; x+=self->shift_step, ++block_idx)
            {
              MAT_ELEM(A,0,2)=x;

              // Rotate and shift image
              applyGeometry(1,Iaux,mI,A,IS_INV,true);

              // Update Wnode
              int i=0;
              FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iaux)
              {
                if (DIRECT_MULTIDIM_ELEM(mask,n)==0)
                continue;
                double pixval=DIRECT_MULTIDIM_ELEM(Iaux,n);
                double *ptrWnode=&MAT_ELEM(Wnode,i,0);
                double *ptrHblock=&MAT_ELEM(Hblock,block_idx,0);
                for (int j=0; j<jmax; j+=unroll, ptrHblock+=unroll, ptrWnode+=unroll)
                {
                  (*(ptrWnode )) +=pixval*(*ptrHblock );
                  (*(ptrWnode+1)) +=pixval*(*(ptrHblock+1));
                  (*(ptrWnode+2)) +=pixval*(*(ptrHblock+2));
                  (*(ptrWnode+3)) +=pixval*(*(ptrHblock+3));
                  (*(ptrWnode+4)) +=pixval*(*(ptrHblock+4));
                  (*(ptrWnode+5)) +=pixval*(*(ptrHblock+5));
                  (*(ptrWnode+6)) +=pixval*(*(ptrHblock+6));
                  (*(ptrWnode+7)) +=pixval*(*(ptrHblock+7));
                }
                for (size_t j=jmax; j<MAT_XSIZE(Wnode); ++j, ptrHblock+=1, ptrWnode+=1)
                (*(ptrWnode )) +=pixval*(*ptrHblock );
                ++i;
              }
            }
          }
        }
      }
    }
    if (IS_MASTER && thArg.thread_id==0)
    progress_bar(last);
  }
}

void ProgImageRotationalPCA::allReduceApplyT(Matrix2D<double> &Wnode_0)
{
}

void ProgImageRotationalPCA::applyT()
{
  Matrix2D<double> &Wnode_0=Wnode[0];
  Wnode_0.initZeros(Npixels,MAT_XSIZE(H));
  taskDistributor->reset();
  thMgr->run(threadApplyT);

  // Gather all Wnodes from all threads
  for (int n=1; n<Nthreads; ++n)
    Wnode_0 += Wnode[n];
  allReduceApplyT(Wnode_0);
  if (IS_MASTER)
    progress_bar(objId.size());
}

// Apply T ================================================================
void threadApplyTt(ThreadArgument &thArg)
{
  ProgImageRotationalPCA *self=(ProgImageRotationalPCA *) thArg.workClass;
  //MpiNode *node=self->node;
  int rank = self->rank;
  ThreadTaskDistributor *taskDistributor=self->taskDistributor;
  std::vector<size_t> &objId=self->objId;
  MetaData &MD=self->MD[thArg.thread_id];

  Image<double> &I=self->I[thArg.thread_id];
  MultidimArray<double> &Iaux=self->Iaux[thArg.thread_id];
  Matrix2D<double> &A=self->A[thArg.thread_id];
  Matrix2D<double> &Hblock=self->Hblock[thArg.thread_id];
  Matrix2D<double> &H=self->H;
  Matrix2D<double> &Wtranspose=self->Wtranspose;
  MultidimArray< unsigned char > &mask=self->mask;

  const size_t unroll=8;
  const size_t nmax=(MAT_XSIZE(Wtranspose)/unroll)*unroll;

  size_t first, last;
  if (IS_MASTER && thArg.thread_id==0)
  {
    std::cerr << "Applying Tt ...\n";
    init_progress_bar(objId.size());
  }
  while (taskDistributor->getTasks(first, last))
  {
    for (size_t idx=first; idx<=last; ++idx)
    {
      // Read image
      I.readApplyGeo(MD,objId[idx]);
      MultidimArray<double> &mI=I();

      // For each rotation and shift
      int block_idx=0;
      for (int mirror=0; mirror<2; ++mirror)
      {
        if (mirror)
        {
          mI.selfReverseX();
          mI.setXmippOrigin();
        }
        for (double psi=0; psi<360; psi+=self->psi_step)
        {
          rotation2DMatrix(psi,A,true);
          for (double y=-self->max_shift_change; y<=self->max_shift_change; y+=self->shift_step)
          {
            MAT_ELEM(A,1,2)=y;
            for (double x=-self->max_shift_change; x<=self->max_shift_change; x+=self->shift_step, ++block_idx)
            {
              MAT_ELEM(A,0,2)=x;

              // Rotate and shift image
              applyGeometry(1,Iaux,mI,A,IS_INV,true);

              // Update Hblock
              for (size_t j=0; j<MAT_XSIZE(Hblock); j++)
              {
                double dotproduct=0;
                const double *ptrIaux=MULTIDIM_ARRAY(Iaux);
                unsigned char *ptrMask=&DIRECT_MULTIDIM_ELEM(mask,0);
                const double *ptrWtranspose=&MAT_ELEM(Wtranspose,j,0);
                for (size_t n=0; n<nmax; n+=unroll, ptrIaux+=unroll, ptrMask+=unroll)
                {
                  if (*(ptrMask ))
                  dotproduct+=(*(ptrIaux ))*(*ptrWtranspose++);
                  if (*(ptrMask+1))
                  dotproduct+=(*(ptrIaux+1))*(*ptrWtranspose++);
                  if (*(ptrMask+2))
                  dotproduct+=(*(ptrIaux+2))*(*ptrWtranspose++);
                  if (*(ptrMask+3))
                  dotproduct+=(*(ptrIaux+3))*(*ptrWtranspose++);
                  if (*(ptrMask+4))
                  dotproduct+=(*(ptrIaux+4))*(*ptrWtranspose++);
                  if (*(ptrMask+5))
                  dotproduct+=(*(ptrIaux+5))*(*ptrWtranspose++);
                  if (*(ptrMask+6))
                  dotproduct+=(*(ptrIaux+6))*(*ptrWtranspose++);
                  if (*(ptrMask+7))
                  dotproduct+=(*(ptrIaux+7))*(*ptrWtranspose++);
                }
                for (size_t n=nmax; n<MAT_XSIZE(Wtranspose); ++n, ++ptrMask, ++ptrIaux)
                if (*ptrMask)
                dotproduct+=(*ptrIaux)*(*ptrWtranspose++);
                MAT_ELEM(Hblock,block_idx,j)=dotproduct;
              }
            }
          }
        }
      }

      // Locate the corresponding index in Matrix H
      // and copy block to disk
      size_t Hidx=idx*2*self->Nangles*self->Nshifts;
      // COSS std::cout << "idx=" << idx << " Hidx=" << Hidx << std::endl;
      // COSS std::cout << "H shape " << MAT_YSIZE(H) << " x " << MAT_XSIZE(H) << std::endl;
      self->writeToHBuffer(thArg.thread_id,&MAT_ELEM(H,Hidx,0));
    }
    if (IS_MASTER && thArg.thread_id==0)
      progress_bar(last);
  }
}

void ProgImageRotationalPCA::applyTt()
{
  // Compute W transpose to accelerate memory access
  Matrix2D<double> &W=Wnode[0];
  Wtranspose.resizeNoCopy(MAT_XSIZE(W),MAT_YSIZE(W));
  FOR_ALL_ELEMENTS_IN_MATRIX2D(Wtranspose)
  MAT_ELEM(Wtranspose,i,j) = MAT_ELEM(W,j,i);

  taskDistributor->reset();
  thMgr->run(threadApplyTt);
  flushHBuffer();
  if (IS_MASTER)
    progress_bar(objId.size());
}

// QR =====================================================================
int ProgImageRotationalPCA::QR()
{
  size_t jQ=0;
  Matrix1D<double> qj1, qj2;
  int iBlockMax=MAT_XSIZE(F)/4;

  for (size_t j1=0; j1<MAT_YSIZE(F); j1++)
  {
    F.getRow(j1,qj1);
    // Project twice in the already established subspace
    // One projection should be enough but Gram-Schmidt suffers
    // from numerical problems
    for (int it=0; it<2; it++)
    {
      for (size_t j2=0; j2<jQ; j2++)
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
        for (size_t i=iBlockMax*4; i<MAT_XSIZE(F); ++i)
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
  return jQ;
}

// Copy H to F ============================================================
void ProgImageRotationalPCA::copyHtoF(int block)
{
    size_t Hidx=block*MAT_XSIZE(H);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(H)
    MAT_ELEM(F,Hidx+j,i)=MAT_ELEM(H,i,j);
}

void ProgImageRotationalPCA::comunicateQrDim(int &qrDim)
{
    std::cerr << "Performing QR decomposition ..." << std::endl;
    qrDim = QR();
}

void ProgImageRotationalPCA::mapMatrix(int qrDim)
{
    FileName(fnRoot+"_matrixH.raw").createEmptyFileWithGivenLength(Nimg*2*Nangles*Nshifts*qrDim*sizeof(double));
    H.mapToFile(fnRoot+"_matrixH.raw",MAT_XSIZE(F),qrDim);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(H)
      MAT_ELEM(H,i,j)=MAT_ELEM(F,j,i);
}

void ProgImageRotationalPCA::applySVD()
{
  // Apply SVD and extract the basis
    std::cerr << "Performing SVD decomposition ..." << std::endl;
    // SVD of W
    Matrix2D<double> U,V;
    Matrix1D<double> S;
    svdcmp(Wnode[0],U,S,V);

    // Keep the first Neigen images from U
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
}


// Run ====================================================================
void ProgImageRotationalPCA::run()
{
  show();
  produceSideInfo();

  // Compute matrix F:
  // Set H pointing to the first block of F
  applyTt();// H=Tt(W)
  copyHtoF(0);
  for (int it=0; it<Nits; it++)
  {
    applyT(); // W=T(H)
    applyTt();// H=Tt(W)
    copyHtoF(it+1);
  }
  H.clear();
  clearHbuffer();

  // QR decomposition of matrix F
  int qrDim;

  comunicateQrDim(qrDim);

  if (qrDim==0)
  REPORT_ERROR(ERR_VALUE_INCORRECT,"No subspace have been found");

  mapMatrix(qrDim);
  F.clear();

  // Apply T
  for (int n=0; n<Nthreads; ++n)
  Hblock[n].resizeNoCopy(2*Nangles*Nshifts,qrDim);
  applyT();

  // Free memory
  H.clear();
  Hblock.clear();

  applySVD();

  // Clean files
  if (IS_MASTER)
  {
    unlink((fnRoot+"_matrixH.raw").c_str());
    unlink((fnRoot+"_matrixF.raw").c_str());
  }
}
