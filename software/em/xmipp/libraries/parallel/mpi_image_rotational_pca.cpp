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

// Empty constructor =======================================================
MpiProgImageRotationalPCA::MpiProgImageRotationalPCA(int argc, char **argv)
{
    node = new MpiNode(argc,argv);
    rank = node->rank;
    if (!IS_MASTER)
        verbose = 0;
}

MpiProgImageRotationalPCA::~MpiProgImageRotationalPCA()
{
  delete node;
}

void MpiProgImageRotationalPCA::selectPartFromMd(MetaData &MDin)
{
  if (IS_MASTER)
  {
    ProgImageRotationalPCA::selectPartFromMd(MDin);
    MDin.write(fnRoot + "_temp.xmd");
    node->barrierWait();
    node->barrierWait();
    unlink((fnRoot + "_temp.xmd").c_str());
  }
  else
  {
    node->barrierWait();
    MDin.read(fnRoot + "_temp.xmd");
    node->barrierWait();
  }
}

void MpiProgImageRotationalPCA::comunicateMatrix(Matrix2D<double> &W)
{
    MPI_Bcast(&MAT_ELEM(W,0,0),MAT_XSIZE(W)*MAT_YSIZE(W),MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void MpiProgImageRotationalPCA::createMutexes(size_t Nimgs)
{
  fileMutex = new MpiFileMutex(node);
  threadMutex = new Mutex();
  taskDistributor = new MpiTaskDistributor(Nimgs, XMIPP_MAX(1,Nimgs/(5*node->size)), node);
}

/** Last part of function applyT */
void  MpiProgImageRotationalPCA::allReduceApplyT(Matrix2D<double> &Wnode_0)
{
    MPI_Allreduce(MPI_IN_PLACE, MATRIX2D_ARRAY(Wnode_0), MAT_XSIZE(Wnode_0)*MAT_YSIZE(Wnode_0),
        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void MpiProgImageRotationalPCA::comunicateQrDim(int &qrDim)
{
  if (IS_MASTER)
    ProgImageRotationalPCA::comunicateQrDim(qrDim);
  node->barrierWait();
  MPI_Bcast(&qrDim,1,MPI_INT,0,MPI_COMM_WORLD);
}

void MpiProgImageRotationalPCA::mapMatrix(int qrDim)
{
// Load the first qrDim columns of F in matrix H
  if (IS_MASTER)
  {
    ProgImageRotationalPCA::mapMatrix(qrDim);
    node->barrierWait();
  }
  else
  {
    node->barrierWait();
    H.mapToFile(fnRoot+"_matrixH.raw",MAT_XSIZE(F),qrDim);
  }
}

void MpiProgImageRotationalPCA::applySVD()
{
  if (IS_MASTER)
    ProgImageRotationalPCA::applySVD();
  node->barrierWait();
}

// Copy H to F ============================================================
void MpiProgImageRotationalPCA::copyHtoF(int block)
{
  if (IS_MASTER)
    ProgImageRotationalPCA::copyHtoF(block);
  node->barrierWait();
}

