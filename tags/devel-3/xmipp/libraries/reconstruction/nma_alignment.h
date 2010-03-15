/***************************************************************************
 *
 * Authors:    Slavica Jonic            slavica.jonic@epfl.ch (2004)
 *             Carlos Oscar             coss@cnb.csic.es
 *
 * Biomedical Imaging Group, EPFL.
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

#ifndef _PROG_NMA_ALIGNMENT
#define _PROG_NMA_ALIGNMENT

#include <vector>
#include <data/progs.h>

/**@defgroup NMAAlignment Alignment with Normal modes
   @ingroup ReconsLibraryPrograms */
//@{
/** NMA Alignment Parameters. */
class Prog_nma_alignment_prm: public Prog_parameters
{
public:
    /// PDB file
    FileName fnPDB;
    
    /// File zith a list of modes
    FileName fnModeList;

    /// Don't modify header
    bool dont_modify_header;
    
    /// Sampling rate
    double sampling_rate;

    /// PDB file
    FileName fnOut;    
public:
    /// List of modes
    std::vector<std::string> modeList;
    
    /// Size of the images in the selfile
    int imgSize;
    
    /// Current image being considered
    const ImageXmipp *currentImg;
    
    /// Current stage of optimization
    int currentStage;
    
    /// Solution of the first stage
    Matrix1D<double> bestStage1;
    
    /// List of all parameters
    std::vector< Matrix1D<double> > listAssignments;
    
    /// List of all image names
    std::vector< FileName > img_names;
      
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
    void produce_side_info();
    
    /** Create deformed PDB */
    FileName createDeformedPDB(const Matrix1D<double> &trial,
        int pyramidLevel) const;

    /** Perform a complete search with the given image and reference
        volume at the given level of pyramid. Return the values
	in the last five positions of trial. */
    void performCompleteSearch(const FileName &fnRandom,
        int pyramidLevel,  int pyramidLevelCont, Matrix1D<double> &trial) const;

    /** Perform a continuous search with the given image and reference
        volume at the given pyramid level. Return the values
	in the last five positions of trial. */
    double performContinuousAssignment(const FileName &fnRandom,
        int pyramidLevel, int pyramidLevelDiscrete,
	Matrix1D<double> &trial) const;

    /** Computes the fitness of a set of trial parameters */
    double computeFitness(Matrix1D<double> &trial) const;

    /** Assign NMA and Alignment parameters to an image */
    Matrix1D<double> assignParameters(const ImageXmipp &img);

    /** Finish processing.
        Close all output files. */
    void finish_processing();
};
//@}
#endif
