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
#ifndef _PROG_CLASSIFY_CL2D_CORE_ANALYSIS
#define _PROG_CLASSIFY_CL2D_CORE_ANALYSIS

#include <parallel/xmipp_mpi.h>
#include <data/metadata.h>

/**@defgroup ClassifyCL2DCore Core analysis for CL2D
   @ingroup ReconsLibrary */
//@{
/** CL2D block description
 */
class CL2DBlock {
public:
	String   block;
	FileName fnLevel;
	FileName fnLevelCore;
	int      level;
};

/// Show CL2Dblock
std::ostream & operator << (std::ostream &out, const CL2DBlock &block);

enum ClassifyCL2DCoreAction
{
	COMPUTE_CORE,
	COMPUTE_STABLE_CORE
};

/** Core analysis parameters. */
class ProgClassifyCL2DCore: public XmippProgram
{
public:
	/** CL2D rootname */
	FileName fnRoot;
	/** CL2D output dir */
	FileName fnODir;
	/** Number of PCA dimensions */
	int NPCA;
	/** Threshold PCA Zscore */
	double thPCAZscore;
	/** Tolerance: How many levels before are allowed two images not to coincide */
	int tolerance;
	/** Action */
	ClassifyCL2DCoreAction action;
public:
    // Mpi node
    MpiNode *node;
    // FileTaskDistributor
    MpiTaskDistributor *taskDistributor;
	// CL2D blocks
	std::vector<CL2DBlock> blocks;
    // MaxLevel
    int maxLevel;
    // Projection size
    size_t Ydim, Xdim;
public:
    /// Empty constructor
    ProgClassifyCL2DCore(int argc, char **argv);

    /// Destructor
    ~ProgClassifyCL2DCore();

    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Usage
    void defineParams();

    /// Produce side info
    void produceSideInfo();

    /// Remove outliers
    void computeCores();

    /// Compute cores
    void computeStableCores();

    /// Gather results
    void gatherResults(int firstLevel, const String &suffix);

    /** Run. */
    void run();
};
//@}
#endif
