/***************************************************************************
 *
 * Authors:    Slavica Jonic            slavica.jonic@impmc.upmc.fr
 *             Carlos Oscar             coss@cnb.uam.es
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#ifndef _PROG_NMA_ALIGNMENT
#define _PROG_NMA_ALIGNMENT

#include <vector>
#include <data/xmipp_program.h>
#include <data/metadata.h>
#include <data/xmipp_image.h>
#include "volume_from_pdb.h"
#include "../../external/condor/ObjectiveFunction.h"
#include "../../external/condor/Vector.h"

/**@defgroup NMAAlignment Alignment with Normal modes
   @ingroup ReconsLibrary */
//@{
/** NMA Alignment Parameters. */
class ProgNmaAlignment: public XmippMetadataProgram
{
public:
    /** MPI version */
    bool MPIversion;
    
    /** Resume computations */
    bool resume;

    /// Reference atomic or pseudo-atomic structure in PDB format
    FileName fnPDB;

    /// Output directory
    FileName fnOutDir;

    /// File with a list of mode filenames
    FileName fnModeList;

    /// Positive scaling factor to scale the initial trust region radius
    double trustradius_scale;
    
    /// Pixel size in Angstroms
    double sampling_rate;
    
    /// Mask file for 2D masking of the projections of the deformed volume
    FileName fnmask;

    /// Center the PDB structure
    bool do_centerPDB;

    /// Real-space projection matching instead of the default combined global (wavelet-space) and local (Fourier central slice) projection matching
    bool projMatch;

    /// Angular sampling step for computing the reference projections (this is the minimum step for real-space projection matching while the minimum step does not exist in the default method)
    double discrAngStep;

    /// Sigma of Gaussian weigthing in Fourier space (parameter of central-slice method)
    double gaussian_DFT_sigma;

    /// Sigma of Gaussian weigthing in real space for spline interpolation in Fourier space (parameter of central-slice method)
    double gaussian_Real_sigma;

    /// Zero-frequency weight (parameter of central-slice method)
    double weight_zero_freq;

    /// Low-pass filter the volume from PDB 
    bool do_FilterPDBVol;

    /// Low-pass cut-off frequency
    double cutoff_LPfilter;

    /// Use pseudo-atoms instead of atoms
    bool useFixedGaussian;

    /// Gaussian standard deviation for pseudo-atoms
    double sigmaGaussian;
 
public:

    // Random generator seed
    int rangen;
    
    // All estimated parameters (with the cost)
    Matrix1D<double> parameters;

    // Trial parameters
    Matrix1D<double> trial;
    
    // Best trial parameters
    Matrix1D<double> trial_best;
    
    // Best fitness
    Matrix1D<double> fitness_min;
    
    // Number of modes
    int numberOfModes;

    // Size of the images in the selfile
    int imgSize;
    
    // Current image being considered
    FileName currentImgName;
    
    // Current stage of optimization
    int currentStage;
    
    // Solution of the first stage
    Matrix1D<double> bestStage1;
    
    // Temporary document file
    MetaData DF_out;

    // Template for temporal filename generation
    char nameTemplate[256];

    // Volume from PDB
    ProgPdbConverter* progVolumeFromPDB;

public:
    /// Empty constructor
    ProgNmaAlignment();

    /// Destructor
    ~ProgNmaAlignment();

    /// Define params
    void defineParams();

    /// Read arguments from command line
    void readParams();

    /// Show
    void show();

   /** Create deformed PDB */
    FileName createDeformedPDB(int pyramidLevel) const;

    /** Perform a complete search with the given image and reference
        volume at the given level of pyramid. Return the values
    in the last five positions of trial. */
    void performCompleteSearch(const FileName &fnRandom,
        int pyramidLevel) const;

    /** Perform a continuous search with the given image and reference
        volume at the given pyramid level. Return the values
    in the last five positions of trial. */
    double performContinuousAssignment(const FileName &fnRandom, int pyramidLevel) const;

    /** Computes the fitness of a set of trial parameters */
    double computeFitness(Matrix1D<double> &trial) const;

    /** Update the best fitness and the corresponding best trial*/
    void updateBestFit(double fitness, int dim);

    /** Create the processing working files.
     * The working files are:
     * nmaTodo.xmd for images to process (nmaTodo = mdIn - nmaDone)
     * nmaDone.xmd image already processed (could exists from a previous run)
     */
    virtual void createWorkFiles();

    /** Produce side info.
        An exception is thrown if any of the files is not found*/
    virtual void preProcess();
    /** Assign NMA and Alignment parameters to an image */
    virtual void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

    /** Write the final parameters. */
    virtual void finishProcessing();

    /** Write the parameters found for one image */
    virtual void writeImageParameters(const FileName &fnImg);
};

class ObjFunc_nma_alignment: public UnconstrainedObjectiveFunction
{
  public:
    ObjFunc_nma_alignment(int _t, int _n=0);
    ~ObjFunc_nma_alignment(){};
    double eval(Vector v, int *nerror=NULL);
};

//@}
#endif
