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

#include <data/volume.h>
#include <data/mask.h>
#include <vector>

/**@defgroup ConvertVol2Pseudo ConvertVol2Pseudo
   @ingroup ReconsLibraryPrograms */
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

class Prog_Convert_Vol2Pseudo
{
public:
    /// Volume to convert
    FileName fnVol;
    
    /// Output volume
    FileName fnOut;
    
    // Mask
    Mask_Params mask_prm;

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

    /// Column for the intensity (if any)
    std::string intensityColumn;
    
    /// Mindistance
    double minDistance;
public:
    /// Read parameters from command line
    void read(int argc, char **argv);
    
    /// show parameters
    void show() const;
    
    /// show usage
    void usage() const;
    
    /// Prepare side info
    void produceSideInfo();
    
    /// Run
    void run();

    /// Place seeds
    void placeSeeds(int Nseeds);
    
    /// Remove seeds
    void removeSeeds(int Nseeds);

    /// Compute average of a volume
    double computeAverage(int k, int i, int j, Matrix3D<double> &V);

    /// Draw a Gaussian on a volume
    void drawGaussian(double k, double i, double j, Matrix3D<double> &V,
        double intensity);

    /// Draw approximation
    void drawApproximation();

    /// Extract region around a Gaussian
    void extractRegion(int idxGaussian, Matrix3D<double> &region,
        bool extended=false) const;

    /// Insert region
    void insertRegion(const Matrix3D<double> &region);

    /// Evaluate region
    double evaluateRegion(const Matrix3D<double> &region) const;

    /// Optimize current atoms
    void optimizeCurrentAtoms();
    
    /// Write results
    void writeResults();
public:
    // Input volume
    VolumeXmipp Vin;
    
    // Current approximation volume
    VolumeXmipp Vcurrent;
    
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
    Matrix1D<double> gaussianTable;
};
//@}
#endif
