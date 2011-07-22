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
#include <classification/pca.h>
#include <parallel/xmipp_mpi.h>

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
public:
    // Input metadata
    MetaData MD;
    // Number of images
    size_t Nimg;
    // Number of angles
    int Nangles;
    // Number of shifts
    int Nshifts;
    // Image size
    int Xdim;
    // Mpi node
    MpiNode *node;
    // SVD matrix
    Matrix2D<double> H;
public:
    // Input image
    Image<double> I;
    // Rotated and shifted image
    MultidimArray<double> Iaux;
    // Geometric transformation
    Matrix2D<double> A;
    // H block
    Matrix2D<double> Hblock;
    // W
    Matrix2D<double> W;
    // W node
    Matrix2D<double> Wnode;
    // W transpose
    Matrix2D<double> Wtranspose;
public:
    /// Empty constructor
    ProgImageRotationalPCA(int argc, char **argv);

    /// Destructor
    ~ProgImageRotationalPCA();

    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Usage
    void defineParams();

    /// Produce side info
    void produceSideInfo();

    /** Apply T.
     * W=T(H).
     */
    void applyT();

    /** Apply Tt.
     * H=Tt(W).
     */
    void applyTt();

    /** Run. */
    void run();
};
//@}
#endif
