/***************************************************************************
 *
 * Authors:    Slavica Jonic            slavica.jonic@impmc.jussieu.fr
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
#include <data/metadata.h>
#include <data/progs.h>
#include "../../external/condor/ObjectiveFunction.h"
#include "../../external/condor/Vector.h"

/**@defgroup NMAAlignment Alignment with Normal modes
   @ingroup ReconsLibraryPrograms */
//@{
/** NMA Alignment Parameters. */
class Prog_nma_alignment_prm: public Prog_parameters
{
public:
    /** MPI version */
    bool MPIversion;
    
    /// PDB file
    FileName fnPDB;
    
    /// File zith a list of modes
    FileName fnModeList;

    /// Scaling factor to scale deformation amplitude
    double scale_defamp;
    
    /// Sampling rate
    double sampling_rate;
    
    ///File name of the symmetry file (or point group)
    FileName symmetry;

    /// Mask file
    FileName fnmask;

    /// Use fixed Gaussian instead of scattering factors
    bool useFixedGaussian;

    /// Fixed Gaussian standard deviation
    double sigmaGaussian;

    /// Center PDB
    bool do_centerPDB;

    /// Output file
    FileName fnOut;    
public:

    // Random generator seed
    int rangen;
    
    // List of modes
    std::vector<std::string> modeList;
    
    // All estimated parameters (with the cost)
    Matrix1D<double> parameters;

    // Trial parameters
    Matrix1D<double> trial;
    
    // Best trial parameters
    Matrix1D<double> trial_best;
    
    // Best fitness
    Matrix1D<double> fitness_min;
    
    // Size of the images in the selfile
    int imgSize;
    
    // Current image being considered
    const ImageXmipp *currentImg;
    
    // Current stage of optimization
    int currentStage;
    
    // Solution of the first stage
    Matrix1D<double> bestStage1;
    
    // List of all parameters
    std::vector< Matrix1D<double> > listAssignments;
    
    // List of all image names
    std::vector< FileName > img_names;   
    
    // Temporary document file
    MetaData DF_out;
      
public:
    /// Empty constructor
    Prog_nma_alignment_prm();

    /// Read argument from command line
    void read(int argc, char **argv);

    /// Show
    void show();

    /// Usage
    void usage();

    /** Produce side info.
        An exception is thrown if any of the files is not found*/
    void produce_side_info(int rank=0);
    
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
    void update_bestfit(double fitness, int dim) const;

    /** Assign NMA and Alignment parameters to an image */
    void assignParameters(ImageXmipp &img);

    /** Finish processing.
        Close all output files. */
    void finish_processing();
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
