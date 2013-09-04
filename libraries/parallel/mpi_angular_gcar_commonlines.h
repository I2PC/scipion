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

#include <parallel/xmipp_mpi.h>
#include <data/metadata.h>
#include <classification/pca.h>


/**@defgroup RotationalPCA Rotational invariant PCA
   @ingroup ReconsLibrary */
//@{
/** Rotational invariant PCA parameters. */
class ProgAngularGCARCommonLines: public XmippProgram
{
public:
	/** Input selfile */
	FileName fnIn;
	/** Output root */
	FileName fnRoot;
    /** Psi step */
    double psi_step;
    /** Maximum shift change */
    double max_shift_change;
    /** Shift step */
    double shift_step;
public:
    // Mpi node
    MpiNode *node;
    // SVD matrix
    Matrix2D<double> W;
public:
    /// Empty constructor
    ProgAngularGCARCommonLines(int argc, char **argv);

    /// Destructor
    ~ProgAngularGCARCommonLines();

    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Usage
    void defineParams();

    /// Produce side info
    void produceSideInfo();

    /** Run. */
    void run();
};
//@}
#endif
