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
#ifndef _PROG_FILTER_PROJECTIONS_HH
#define _PROG_FILTER_PROJECTIONS_HH

#include <data/docfile.h>
#include <data/funcs.h>
#include <data/volume.h>
#include <data/mask.h>
#include <vector>

/**@defgroup Filter_projections Select which projections proceed to reconstruction
   @ingroup ReconsLibraryPrograms */
//@{
/* Filter projection Program Parameters ----------------------------------- */
/** Parameter class for the project program */
class Prog_Filter_Projections_Parameters
{
public:
    /// Filename of the input docfile
    FileName fn_in;

    /// Rootname of the output docfile and selfile
    FileName fn_out;

    /// Filename of the docfile with the score
    FileName fn_score;

    /// Percentil of the score
    double percentil_score;

    /// Filename of the docfile with the cost
    FileName fn_cost;

    /// Percentil of the cost
    double percentil_cost;

    /// FileName of the first docfile for movements
    FileName fn_movement0;
    
    /// Maximum change in angles
    double angleLimit;
    
    /// Maximum change in shifts
    double shiftLimit;

    /// FileName of the volume for fitting
    FileName fn_vol;
    
    /// First fitting radius
    double r1;
    
    /// Second fitting radius
    double r2;
    
    /// Percentil of the normalization
    double percentil_normalization;

    /// Number of threads
    int Nthreads;
public:
    // Docfile with the images
    DocFile DF_in;

    // Docfile with the score
    DocFile DF_score;

    // Docfile with the cost
    DocFile DF_cost;

    // Docfile with the movements 0
    DocFile DF_movement0;

    // Volume for the fitting
    VolumeXmipp V;

    // Mask for the normalization
    Matrix2D<double> dhardMask;

    // Mask for the normalization
    Matrix2D<int> ihardMask;
    
    // Mask for the reconstruction
    Matrix2D<double> softMask;

    // Vector of which images are valid and which are not
    std::vector<bool> valid;
    
    // Vector of correlations
    std::vector<double> correlations;
public:
    /** Read from a command line.
        An exception might be thrown by any of the internal conversions,
        this would mean that there is an error in the command line and you
        might show a usage message. */
    void read(int argc, char **argv);

    /** Usage message.
        This function shows the way of introducing these parameters. */
    void usage();

    /** Show parameters. */
    void show();

    /** Produce side information. */
    void produce_side_info();

    /** Apply filters */
    void run();
};
//@}
#endif
