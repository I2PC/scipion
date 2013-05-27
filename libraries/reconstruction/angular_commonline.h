/***************************************************************************
 *
 * Authors:     Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es)
 *
 * Universidad San Pablo C.E.U.
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

#ifndef _PROG_ANGULAR_COMMON_LINES
#define _PROG_ANGULAR_COMMON_LINES

#include <data/xmipp_funcs.h>
#include <data/metadata.h>
#include <data/numerical_tools.h>
#include <data/matrix2d.h>
#include <vector>
#include <data/symmetries.h>

/**@defgroup AngularCommonLines angular_common_lines
    (angular assignment through common lines)
   @ingroup ReconsLibrary */
//@{

// Forward declaration
class Prog_Angular_CommonLine;

// Class for solving the optimization problem
class EulerSolver: public DESolver {
public:
    bool show;
    int NToSolve, Ndim, Nimg;
    double roti, tilti, psii;
    double rotj, tiltj, psij;
    Matrix1D<double> normali, normalj, commonline, commonlinei, commonlinej;
    const Prog_Angular_CommonLine *parent;
    const Matrix1D<int> * alreadyOptimized;
    const Matrix1D<double> * currentSolution;
    const MultidimArray<int> * imgIdx;
    
    Matrix1D<double> imgAvgCorrelation;
    Matrix1D<double> imgMinCorrelation;
    MultidimArray<double> correlationMatrix;
public:
    // Constructor
    EulerSolver(int dim, int pop,
        const Matrix1D<int> &newAlreadyOptimized,
        const Matrix1D<double> &newCurrentSolution,
        const MultidimArray<int> &newImgIdx,
        const Prog_Angular_CommonLine *newParent);

    // Energy function for the solver
    double EnergyFunction(double trial[],bool &bAtSolution);

    // Distance between two common lines
    double similarityBetweenTwoLines(int imgi, int imgj);
    
    // Set Debug
    void setShow(bool newShow);
};

/** Class for running the program */
class Prog_Angular_CommonLine {
public:
    /** Filename of the selfile with the images */
    FileName fnSel;
    
    /** Output filename */
    FileName fnOut;
    
    /** Symmetry file */
    FileName fnSym;
    
    /** Number of generations */
    int NGen;
    
    /** Number of images in a group */
    int NGroup;

    /** Try initial solution */
    bool tryInitial;

    /** Read parameters from command line */
    void read(int argc, const char **argv);

    /** Show parameters */
    void show() const;
    
    /** Usage */
    void usage() const;

    /** Produce side info */
    void produceSideInfo();
    
    /** Optimize a single group*/
    double optimizeGroup(const Matrix1D<int> &imgIdx,
        Matrix1D<double> &solution, bool show=false);
    
    /** Optimize all images */
    void optimize(Matrix1D<double> &solution);
    
    /** Optimize what we have done so far */
    double realignCurrentSolution();
    
    /** Try a solution */
    double trySolution(const Matrix1D<double> &solution);
    
    /** Compute clusters */
    double computeClusters(const MultidimArray<double> &correlationMatrix,
        std::vector< std::vector<int> > &clusters, 
        MultidimArray<double> &worseCorrelationMatrix,
        MultidimArray<double> &bestCorrelationMatrix, bool show=false) const;

    /** Split the alignment in two clusters and find an image to remove */
    std::vector<int> removeViaClusters(
        const MultidimArray<double> &correlationMatrix);

    /** Run */
    void run();
public:
    // Selfile with images
    MetaData SF;
    
    // Optimizer
    EulerSolver *solver;
    
    // Input images
    std::vector< MultidimArray<double> > img;
    
    // Initial solution
    Matrix1D<double> initialSolution;
    
    // Already optimized images
    Matrix1D<int> alreadyOptimized;
    
    // Current optimized solution, angles of the already
    // optimized images
    Matrix1D<double> currentSolution;

    // Correlation of the images individually
    Matrix1D<double> currentImgAvgCorrelation;
    
    // Worse Correlation of the images individually
    Matrix1D<double> currentImgMinCorrelation;
    
    // Correlation between any pair of images
    MultidimArray<double> currentCorrelationMatrix;
    
    // Best line in image i matching with image j
    Matrix2D<int> bestLine;

    // Best correlation of a line in image i matching with image j
    Matrix2D<double> bestLineCorrelation;

    // Symmetry list
    SymList SL;

    // Radon transform of the images
    std::vector< std::vector< MultidimArray<double> > > radon;

    // Left matrices for the symmetry transformations
    std::vector< Matrix2D<double> > L;

    // Right matrices for the symmetry transformations
    std::vector< Matrix2D<double> > R;
};

//@}
#endif
