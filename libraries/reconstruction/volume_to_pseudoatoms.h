/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
/* This file contains functions related to the Radon Transform */

#ifndef _CONVERT_VOL2PSEUDO_HH
#define _CONVERT_VOL2PSEUDO_HH

#include <data/xmipp_image.h>
#include <data/mask.h>
#include <data/xmipp_threads.h>
#include <data/xmipp_program.h>
#include <vector>

/**@defgroup ConvertVol2Pseudo ConvertVol2Pseudo
   @ingroup ReconsLibrary */
//@{
/// Pseudoatom class
class PseudoAtom
{
public:
    /// Location
    Matrix1D<double> location;
    
    /// Intensity
    double           intensity;
    
    /// Empty constructor
    PseudoAtom();

    /// Show pseudo atom
    friend std::ostream& operator << (std::ostream &o, const PseudoAtom &f);
};

/// Comparison between pseudo atoms
bool operator <(const PseudoAtom &a, const PseudoAtom &b);

// Forward declaration
class ProgVolumeToPseudoatoms;

// Thread parameters
struct Prog_Convert_Vol2Pseudo_ThreadParams
{
    int myThreadID;
    ProgVolumeToPseudoatoms *parent;
    int Nintensity;
    int Nmovement;
};

class ProgVolumeToPseudoatoms: public XmippProgram
{
public:
    /// Volume to convert
    FileName fnVol;
    
    /// Output volume
    FileName fnOut;
    
    // Mask
    Mask mask_prm;

    // Use mask
    bool useMask;

    /// Sigma
    double sigma;
    
    /// Maximum error (as a percentage)
    double targetError;

    /// Stop criterion
    double stop;
    
    /// Initial seeds
    int initialSeeds;
    
    /// Grow seeds
    double growSeeds;
    
    /// Allow gaussians to move
    bool allowMovement;
    
    /// Allow gaussians to vary intensity
    bool allowIntensity;

    /** Intensity fraction.
        In case intensity is not allowed to change, this fraction
        is multiplied by the intensity range and all atoms will have
        this intensity value. */
    double intensityFraction;

    /// Column for the intensity (if any)
    std::string intensityColumn;
    
    /// Mindistance
    double minDistance;

    /// Penalization
    double penalty;
    
    /// Number of threads
    int numThreads;
    
    /// Sampling rate
    double sampling;
    
    /// N closest atoms for the distance histogram
    size_t Nclosest;

    /// Don't scale the atom weights at the end
    bool dontScale;

    /// Binarize
    bool binarize;

    /// Threshold for the binarization
    double threshold;
public:
    /// Read parameters from command line
    void readParams();
    
    /// show parameters
    void show() const;
    
    /// define parameters
    void defineParams();
    
    /// Prepare side info
    void produceSideInfo();
    
    /// Run
    void run();

    /// Place seeds
    void placeSeeds(int Nseeds);
    
    /// Remove seeds
    void removeSeeds(int Nseeds);

    /// Remove too close seeds
    void removeTooCloseSeeds();

    /// Compute average of a volume
    double computeAverage(int k, int i, int j, MultidimArray<double> &V);

    /// Draw a Gaussian on a volume
    void drawGaussian(double k, double i, double j, MultidimArray<double> &V,
        double intensity);

    /// Draw approximation
    void drawApproximation();

    /// Extract region around a Gaussian
    void extractRegion(int idxGaussian, MultidimArray<double> &region,
        bool extended=false) const;

    /// Insert region
    void insertRegion(const MultidimArray<double> &region);

    /// Evaluate region
    double evaluateRegion(const MultidimArray<double> &region) const;

    /// Optimize current atoms
    void optimizeCurrentAtoms();
    
    /// Optimize current atoms (thread)
    static void* optimizeCurrentAtomsThread(void * threadArgs);
    
    /// Write results
    void writeResults();
public:
    // Input volume
    Image<double> Vin;
    
    // Current approximation volume
    Image<double> Vcurrent;
    
    // Energy of the difference
    double energyDiff;
    
    // Maximum percentage diffence
    double percentageDiff;
    
    // Original energy
    double energyOriginal;
    
    // List of atoms
    std::vector< PseudoAtom > atoms;
    
    // Maximum radius
    double sigma3;
    
    // Percentil 1
    double percentil1;

    // Range
    double range;

    // Small atom intensity
    double smallAtom;

    // Gaussian table
    MultidimArray<double> gaussianTable;
    
    // Barrier
    barrier_t barrier;
    
#define KILLTHREAD -1
#define WORKTHREAD  0
    // Thread operation code
    int threadOpCode;
    
    // Pointer to thread arguments
    Prog_Convert_Vol2Pseudo_ThreadParams *threadArgs;
    
    // Pointer to threads
    pthread_t *threadIds;
};
//@}
#endif
