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

#ifndef _PROG_COMMONLINES_HH
#  define _PROG_COMMONLINES_HH

#include <data/xmipp_funcs.h>
#include <data/metadata.h>
#include <data/multidim_array.h>
#include <data/numerical_tools.h>
#include <data/xmipp_program.h>
#include <iostream>

/**@defgroup CommonLinesProgram Common Lines (find common lines between projections)
   @ingroup ReconsLibrary */
//@{
/* Parameters -------------------------------------------------------------- */
/// Commonline
class CommonLine
{
public:
    /// Angle of the best common line in image i
    double angi;
    /// Angle of the best common line in image j
    double angj;
    /// Distance between both common lines
    double distanceij;
public:
    /// Empty constructor
    CommonLine();
};

/// CommonLine Parameters
class ProgCommonLine: public XmippProgram
{
public:
    /// input file
    FileName        fn_sel;
    /// output file
    FileName        fn_out;
    /// Output style
    String          outputStyle;
    /// Scale output measure
    bool scaleDistance;
    /// Low pass filter
    double          lpf;
    /// High pass filter
    double          hpf;
    /// Angular sampling
    double          stepAng;
    /// Qualify
    bool            qualify;

    /// Memory limit
    double          mem;
    /// Number of threads
    int             Nthr;
    /// Number of processors
    int             Nmpi;
public:
    /** Read parameters from command line. */
    void readParams();

    /** Define parameters */
    void defineParams();

    /** Produce Side Information */
    void produceSideInfo();

    /** Process block */
    void processBlock(int i, int j);

    /** Get and prepare block */
    void getAndPrepareBlock(int i,
        std::vector< MultidimArray<std::complex<double> > > &blockRTFs,
        std::vector<MultidimArray<double> > &blockRTs);

    /** Show parameters */
    void show();

    /** Qualify common lines */
    void qualifyCommonLines();

    /** Write results */
    void writeResults();
    
    /** Run */
    void run(int rank=0);
public:
    // Block size
    int Nblock;

    // Input selfile
    MetaData SF;

    // Number of images
    int Nimg;

    // Xdim size of the images
    int Xdim;

    // Common line matrix
    std::vector<CommonLine> CLmatrix;

    // Common line matrix
    std::vector<double> qualification;
};
//@}

#endif
