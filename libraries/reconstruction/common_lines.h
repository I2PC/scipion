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
    /// Index of the maximum
    /// jmax=-5 -> line j has to be shifted 5 pixels to the left  to match line i
    /// jmax= 5 -> line j has to be shifted 5 pixels to the right to match line i
    int jmax;
    /// Percentile (good common lines have very high percentiles)
    double percentile;
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
    /// Outlier fraction
    double outlierFraction;
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
    /// MPI Rank
    int rank;
public:
    /** Empty constructor */
    ProgCommonLine(): rank(0) {};

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

    /** Solve for shifts */
    void solveForShifts();

    /** Write results */
    void writeResults();
    
    /** Run */
    void run();
public:
    // Block size
    int Nblock;

    // Input selfile
    MetaData SF;

    // Number of images
    int Nimg;

    // Xdim size of the images
    size_t Xdim;

    // Common line matrix
    std::vector<CommonLine> CLmatrix;

    // Estimated shifts
    Matrix1D<double> shift;
};

#define POW2(x) (x*x)


/**
* Class to hold the common lines info between
* two projection images k1 and k2
* It will store the indexes of commonline in k1 and k2
* and the vector direction
*/
class CommonLineInfo
{
public:
    size_t nRays; //Number of rays in Fourier space
    int k[2]; //image projection indexes;
    int idx[2]; //commonlines indexes in k1 and k2;
    DVector vector; //direction vector

    CommonLineInfo(int n, int k1=0, int k2=0)
    {
        nRays = n;
        setImages(k1, k2);
        idx[0] = idx[1] = 0;
        vector.initZeros(3);
    }

    int getImage(int i)
    {
        return k[i];
    }

    int getIndex(int i)
    {
        return idx[i];
    }

    void setImages(int k1, int k2)
    {
        k[0] = k1;
        k[1] = k2;
    }

    void setIndex(int i, double theta)
    {
        idx[i] = ROUND(theta/TWOPI * nRays) % nRays;
    }
}
;//class CommonLineInfo

void randomQuaternions(int k, DMatrix &qArray);

void saveMatrix(const char *fn, DMatrix &array);

void quaternionToMatrix(const DVector &q, DMatrix &rotMatrix);

void quaternionCommonLines(const DMatrix &quaternions, CommonLineInfo &clInfo);

void commonlineMatrixCheat(const DMatrix &quaternions, size_t nRays,
                           DMatrix &clMatrix, DMatrix &clCorr);

void anglesRotationMatrix(const DMatrix &clMatrix, size_t nRays, int i, int j, DMatrix &U);

#define SMALL_TRIANGLE -101
/** Negative output means error
 * -101 Triangle too small
 */
int tripletRotationMatrix(const DMatrix &clMatrix, size_t nRays, int k1, int k2, int k3, DMatrix &R);

void computeSyncMatrix(const DMatrix &clMatrix, size_t nRays, DMatrix &sMatrix, DMatrix * pQuaternions=NULL);

void rotationsFromSyncMatrix(const DMatrix &sMatrix, DMatrix * pQuaternions = NULL);



//@}

#endif
