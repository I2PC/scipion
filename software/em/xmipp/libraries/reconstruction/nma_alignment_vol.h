/***************************************************************************
 *
 * Authors:    Carlos Oscar             coss@cnb.csic.es
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
#include <condor/ObjectiveFunction.h>
#include <condor/Vector.h>

#include "data/xmipp_program.h"
#include "data/metadata.h"
#include "data/xmipp_image.h"
#include "volume_from_pdb.h"

/**@defgroup NMAAlignmentVol Alignment of volumes with Normal modes
   @ingroup ReconsLibrary */
//@{
/** NMA Alignment Parameters. */
class ProgNmaAlignmentVol: public XmippMetadataProgram
{
public:
    /** MPI version */
    bool MPIversion;
    
    /** Resume computations */
    bool resume;

    /// Reference atomic or pseudo-atomic structure in PDB format
    FileName fnPDB;

    /// Output PDB
    FileName fnOutPDB;

    /// Output directory
    FileName fnOutDir;

    /// File with a list of mode filenames
    FileName fnModeList;

    /// Pixel size in Angstroms
    double sampling_rate;
    
    /// Mask file for 2D masking of the projections of the deformed volume
    FileName fnmask;

    /// Center the PDB structure
    bool do_centerPDB;

    /// Low-pass filter the volume from PDB 
    bool do_FilterPDBVol;

    /// Low-pass cut-off frequency
    double cutoff_LPfilter;

    /// Use pseudo-atoms instead of atoms
    bool useFixedGaussian;

    /// Gaussian standard deviation for pseudo-atoms
    double sigmaGaussian;
 
    /// Align volumes
    bool alignVolumes;

    /// Trust radius scale
    double trustradius_scale;
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
    double fitness_min;
    
    // Number of modes
    int numberOfModes;

    // Size of the volumes in the selfile
    int imgSize;
    
    // Current volume being considered
    FileName currentVolName;
    
    // Template for temporal filename generation
    char nameTemplate[256];

    // Volume from PDB
    ProgPdbConverter* progVolumeFromPDB;

    // Volume that is being fitted
    Image<double> V, Vdeformed;

    // Mask
    MultidimArray<int> mask;

public:
    /// Empty constructor
    ProgNmaAlignmentVol();

    /// Destructor
    ~ProgNmaAlignmentVol();

    /// Define params
    void defineParams();

    /// Read arguments from command line
    void readParams();

    /// Show
    void show();

   /** Create deformed PDB */
    FileName createDeformedPDB() const;

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
    /** Assign NMA and Alignment parameters to a volume */
    virtual void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

    /** Write the final parameters. */
    virtual void finishProcessing();

    /** Write the parameters found for one image */
    virtual void writeVolumeParameters(const FileName &fnImg);
};

class ObjFunc_nma_alignment_vol: public UnconstrainedObjectiveFunction
{
  public:
    ObjFunc_nma_alignment_vol(int _t, int _n=0);
    ~ObjFunc_nma_alignment_vol(){};
    double eval(Vector v, int *nerror=NULL);
};

//@}
#endif
