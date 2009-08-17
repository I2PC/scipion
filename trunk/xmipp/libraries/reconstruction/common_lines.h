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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#ifndef _PROG_COMMONLINES_HH
#  define _PROG_COMMONLINES_HH

#include <data/funcs.h>
#include <data/selfile.h>
#include <iostream>

/**@defgroup CommonLinesProgram Common Lines (find common lines between projections)
   @ingroup ReconsLibraryPrograms */
//@{
/* Parameters -------------------------------------------------------------- */
typedef enum {CORRENTROPY=0, CORRELATION=1} DistanceType;

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
class CommonLine_Parameters
{
public:
    /// input file
    FileName        fn_sel;
    /// output file
    FileName        fn_out;
    /// Distance measure
    DistanceType    distance;
    /// Low pass filter
    double          lpf;
    /// High pass filter
    double          hpf;
    /// Angular sampling
    double          stepAng;

    /// Memory limit
    double          mem;
    /// Number of threads
    int             Nthr;
    /// Number of processors
    int             Nmpi;
public:
    /** Read parameters from command line. */
    void read(int argc, char **argv);

    /** Usage */
    void usage();

    /** Produce Side Information */
    void produceSideInfo();

    /** Process block */
    void processBlock(int i, int j);

    /** Get and prepare block */
    void getAndPrepareBlock(int i,
        std::vector< Matrix2D<double> > &blockImgs);

    /** Show parameters */
    friend std::ostream & operator << (std::ostream &out,
        const CommonLine_Parameters &prm);

    /** Write results */
    void writeResults();
public:
    // Block size
    int Nblock;

    // Input selfile
    SelFile SF;

    // Number of images
    int Nimg;

    // Sigma for the correntropy
    double sigma;

    // Common line matrix
    std::vector<CommonLine> CLmatrix;
};

/** Main program */
void ROUT_commonlines(CommonLine_Parameters &prm, int rank=0);
//@}

#endif
