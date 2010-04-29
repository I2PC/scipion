/***************************************************************************
 *
 * Authors:    Carlos Oscar           coss@cnb.csic.es (2010)
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
#ifndef _PROG_ANALYZE_CLUSTER
#define _PROG_ANALYZE_CLUSTER

#include <data/funcs.h>
#include <data/metadata.h>

/**@defgroup AnalyzeClusterProgram analyze cluster
   @ingroup ReconsLibraryPrograms */
//@{
/** Analyze cluster parameters. */
class Prog_analyze_cluster_prm
{
public:
    /** Filename selection file containing the images */
    FileName fnSel;

    /** Filename reference image */
    FileName fnRef;

    /**  Filename output extension */
    FileName oext;

    /** Boolean to align input images */
    bool align;

    /** PCA dimension */
    int NPCA;

    /** Number of iterations */
    int Niter;

    /** Distance Threshold */
    double distThreshold;

    /** Don't mask*/
    bool dontMask;
public:
    // SelFile images
    std::vector< FileName > classfile;

    // Image holding current reference
    ImageXmipp Iref;

    // Mask of the background
    Matrix2D<int> mask;

    // Set of images assigned to the class
    std::vector< Matrix1D<float> * > Iclass;

    // Set of images assigned to the class
    std::vector< Matrix1D<float> * > Iclassorig;

    // Set of basis functions
    std::vector< Matrix1D<double> * > PCAbasis;

    // Set of distances
    Matrix1D<double> distance;

    // Number of pixels in the mask
    int Npixels;
public:
    /// Read argument from command line
    void read(int argc, char **argv);

    /// Show
    void show();

    /// Usage
    void usage();

    /// Produce side info
    void produceSideInfo();

    /// Learn basis
    void learnPCABasis();

    /// Project on basis
    void projectOnPCABasis(Matrix2D<double> &CtY);

    /// Evaluate distance
    void evaluateDistance(Matrix2D<double> &proj);

    /// Main routine
    void run();
};
//@}
#endif
