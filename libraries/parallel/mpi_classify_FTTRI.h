/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2012)
 *             Wang Xia                wangxia@lsec.cc.ac.cn
 *             Guoliang Xu             xuguo@lsec.cc.ac.cn
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
#ifndef _PROG_CLASSIFY_FTTRI
#define _PROG_CLASSIFY_FTTRI

#include <parallel/xmipp_mpi.h>
#include <data/metadata.h>


/**@defgroup ClassifyFTTRI Classification with Fourier Transform based Translational and Rotational Invariants
   @ingroup ReconsLibrary */
//@{

/** Epsilon Class. */
class EpsilonClass
{
public:
	std::vector<size_t> memberIdx;
	std::vector<int> neighbours;
};

/** Core analysis parameters. */
class ProgClassifyFTTRI: public XmippProgram
{
public:
	/// Input selfile
	FileName fnIn;

	/// Output rootname
	FileName fnRoot;

	/// Padding factor
	double pad;

	/// Maximum frequency (normalized to 0.5)
	double fmax;

	/// Maximum frequency in pixels
	size_t Rmax;

	/// Zoom factor for polar conversion
	double zoom;

	/// First weight
	double sigma1;

	/// Second weight
	double sigma2;

	/// Desired number of classes
	size_t nref;

	/// Minimum number of images in a class
	size_t nMinImages;

	/// Number of iterations
	int Niter;

	/// Do phase optimization
	bool doPhase;
public:
	// Name of the invariants;
	FileName fnFTTRI;
	// Final size of the invariant
	int FTTRIXdim;
	int FTTRIYdim;
    // Mpi node
    MpiNode *node;
    // FileTaskDistributor
    FileTaskDistributor *taskDistributor;
    // MetaData In
    MetaData mdIn;
    // Object Ids
    std::vector<size_t> imgsId;
    // Padded size
    int padXdim;
    // Set of classes
    std::vector< EpsilonClass > epsilonClasses;
    // Initial range for epsilon
    double dMin, dMax;
    // Copy of bestObjective
    double bestObjective;
    // Copy of bestEpsilon
    double bestEpsilon;
    // Copy of bestClass
    std::vector< EpsilonClass > bestEpsilonClasses;
    // Epsilon of each class
	Matrix1D<double> classEpsilon;
	// Epsilon max is the maximum of all class epsilons
	double epsilonMax;
    // Assigned to class
    Matrix1D<unsigned char> notAssigned;
    // Assigned to class init
    Matrix1D<unsigned char> notAssigned0;
    // All FTTRI centroids are loaded in memory
	Image<double> fttriCentroids;
	// All Image centroids are loaded in memory
	Image<double> imageCentroids;
public:
    /// Empty constructor
    ProgClassifyFTTRI(int argc, char **argv);

    /// Destructor
    ~ProgClassifyFTTRI();

    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Usage
    void defineParams();

    /// Produce side info
    void produceSideInfo();

    /// Produce invariants
    void produceFTTRI();

    /// Estimate first epsilon range
    void estimateEpsilonInitialRange();

    /// Distance between two invariants
    double fttri_distance(const MultidimArray<double> &fttri_i,
                          const MultidimArray<double> &fttri_j);

    /// Skip random number of unassigned classes
    void skipRandomNumberOfUnassignedClasses(
    		size_t &currentPointer, size_t remaining);

    /// Function to optimize
    size_t wrapperFitness(double epsilon);

    /// Epsilon classification
    void epsilonClassification(double epsilon);

    /// Search for optimal epsilon
    void searchOptimalEpsilon();

    /// Remove small classes
    void removeSmallClasses();

    /// Find farthest image
    int findFarthestFTTRI(const MultidimArray<double> &fttriSeed,
    		const EpsilonClass &class_i);

    int findFarthest(const MultidimArray<double> &seed,
            const EpsilonClass &class_i, bool FTTRI);

    /// Split large classes
    void splitLargeClasses(bool FTTRI);

    /// Split large classes
    void splitLargeFTTRIClasses();

    /// Compute centroid neighbours
    void computeClassNeighbours(bool FTTRI);


    /// Reassign images to image classes
    size_t reassignImagesToClasses(bool FTTRI);

    /** Write results */
    void writeResults(bool FTTRI);

    /** Align images within classes */
    void alignImagesWithinClasses();

    /** compute class centroids */
    void computeClassCentroids(bool FTTRI);

    /** Run. */
    void run();
};
//@}
#endif
