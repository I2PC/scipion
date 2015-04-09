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
#ifndef _PROG_IMAGE_ROTATIONAL_PCA
#define _PROG_IMAGE_ROTATIONAL_PCA

#include <data/metadata.h>
#include <data/xmipp_program.h>
#include <data/xmipp_threads.h>
#include <classification/pca.h>

#define IS_MASTER (rank == 0)

/**@defgroup RotationalPCA Rotational invariant PCA
   @ingroup ReconsLibrary */
//@{
/** Rotational invariant PCA parameters. */
class ProgImageRotationalPCA: public XmippProgram
{
public:
	/** Input selfile */
	FileName fnIn;
	/** Output root */
	FileName fnRoot;
	/** Number of eigenvectors */
	int Neigen;
	/** Number of iterations */
	int Nits;
    /** Psi step */
    double psi_step;
    /** Maximum shift change */
    double max_shift_change;
    /** Shift step */
    double shift_step;
    /** Maximum number of images */
    int maxNimgs;
    /** Number of threads */
    int Nthreads;
    /** Rank, used later for MPI */
    int rank;

public:
    // Input metadata
    std::vector<MetaData> MD;
    // Number of images
    size_t Nimg;
    // Number of angles
    int Nangles;
    // Number of shifts
    int Nshifts;
    // Image size
    size_t Xdim;
    // Number of pixels
    int Npixels;
    // Mpi node
    //MpiNode *node;
    // H buffer
    std::vector<double *> Hbuffer;
    // H buffer destination addresses
    std::vector<double *> HbufferDestination;
    // Buffer maximum length
    static const int HbufferMax=20;
    // Mpi file lock
    Mutex *fileMutex;
    // Thread mutex
    Mutex *threadMutex;
    // SVD matrix
    Matrix2D<double> H;
    // SVD matrix
    Matrix2D<double> F;
public:
    // Input image
    std::vector< Image<double> > I;
    // Rotated and shifted image
    std::vector< MultidimArray<double> > Iaux;
    // Geometric transformation
    std::vector< Matrix2D<double> > A;
    // H block
    std::vector< Matrix2D<double> > Hblock;
    // W node
    std::vector< Matrix2D<double> > Wnode;
    // W transpose
    Matrix2D<double> Wtranspose;
    // Mask
    MultidimArray< unsigned char > mask;
    // FileTaskDistributor
    ThreadTaskDistributor *taskDistributor;
    // Thread Manager
    ThreadManager *thMgr;
    // Vector of object ids
    std::vector<size_t> objId;
public:
    /// Empty constructor
    ProgImageRotationalPCA();

    /// Destructor
    virtual ~ProgImageRotationalPCA();

    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Usage
    void defineParams();

    /// Produce side info
    void produceSideInfo();

    /** Write to H buffer */
    void writeToHBuffer(int idx, double *dest);

    /** Flush buffer */
    void flushHBuffer();

    /** Clear H buffer */
    void clearHbuffer();

    /** Apply T.
     * W=T(H).
     */
    void applyT();

    /** Apply Tt.
     * H=Tt(W).
     */
    void applyTt();

    /** QR decomposition.
     * In fact, only Q is computed. It returns the number of columns
     * of Q different from 0.
     */
    int QR();

    /** Run. */
    void run();

    /********************** Following functions should be overwritten in MPI version ******/
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
