/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2002)
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
#ifndef _MPI_PROG_IMAGE_ROTATIONAL_PCA
#define _MPI_PROG_IMAGE_ROTATIONAL_PCA

#include <parallel/xmipp_mpi.h>
#include <data/metadata.h>
#include <classification/pca.h>
#include <reconstruction/image_rotational_pca.h>

/**@defgroup RotationalPCA Rotational invariant PCA
   @ingroup ReconsLibrary */
//@{
/** Rotational invariant PCA parameters. */
class MpiProgImageRotationalPCA: public ProgImageRotationalPCA
{
public:
    // Mpi node
    MpiNode *node;
    // H buffer

    /// Empty constructor
    MpiProgImageRotationalPCA(int argc, char **argv);

    /// Destructor
    ~MpiProgImageRotationalPCA();

    /** Read input images */
    virtual void selectPartFromMd(MetaData &MDin);

    /** Comunicate matrix, only meanful for MPI */
    virtual void comunicateMatrix(Matrix2D<double> &W);

    /** Create mutexes and distributor */
    virtual void createMutexes(size_t Nimgs);

    /** Last part of function applyT */
    virtual void allReduceApplyT(Matrix2D<double> &Wnode_0);

    /** Comunicate int param */
    virtual void comunicateQrDim(int &qrDim);

    /** Map matrix */
    virtual void mapMatrix(int qrDim);

    /** Apply SVD */
    virtual void applySVD();

    /** Copy H to F. */
    virtual void copyHtoF(int block);
};
//@}
#endif
