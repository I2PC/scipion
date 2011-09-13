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

#include <data/xmipp_funcs.h>
#include <data/metadata.h>
#include <data/metadata_extension.h>
#include <data/basic_pca.h>
#include <data/xmipp_program.h>

/**@defgroup AnalyzeClusterProgram Analyze cluster
   @ingroup ClassificationLibrary */
//@{
/** Analyze cluster parameters. */
class ProgAnalyzeCluster: public XmippProgram
{
public:
    /** Filename selection file containing the images */
    FileName fnSel;

    /**  Output metadata */
    FileName fnOut;

    /**  Output basis stack */
    FileName fnOutBasis;

    /** PCA dimension */
    int NPCA;

    /** Number of iterations */
    int Niter;
    
    /** Distance Threshold */
    double distThreshold;
    
    /** Don't mask*/
    bool dontMask;
public:
    // Produce PCA basis
    bool basis;

    // SelFiles of images
    MetaData SFin, SFout;

    // Mask of the background
    MultidimArray<int> mask;

    // PCA Analyzer
    PCAMahalanobisAnalyzer pcaAnalyzer;

public:
    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Define parameters
    void defineParams();

    /// Produce side info
    void produceSideInfo();

    /// Produce basis
    void produceBasis(MultidimArray<double> &basis);

    /// Main routine
    void run();
};
//@}
#endif
