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
/*****************************************************************************/
/* Variable and prototype definitions for the Numerical Core                 */
/*****************************************************************************/
#ifndef _NUMERICAL_TOOLS_HH
#define _NUMERICAL_TOOLS_HH

#include "numerical_recipes.h"
#include "matrix1d.h"
#include "matrix2d.h"
#include "multidim_array.h"

template<typename T> class Matrix1D;
template<typename T> class Matrix2D;

/// @defgroup NumericalTools Numerical Tools
/// @ingroup DataLibrary

/** Generate random permutation
 * @ingroup NumericalTools
 *
 * Generate a random permutation of the numbers between 0 and N-1
 */
void randomPermutation(int N, MultidimArray<int>& result);

/** Optimize using Powell's method.
  * @ingroup NumericalTools
  *
  * See Numerical Recipes Chapter 10.
  *
  * Problem: Minimize f(x) starting at point p. n is the dimension of x.
  * If changes are smaller than ftol then stop. The number of iterations is
  * returned in iter, fret contains the function value at the minimum and p
  * contains the minimum.
  *
  * Watch out that the evaluating function f must consider that x starts at
  * index 1, at least, and goes until n. i0 is used to allow optimization in
  * large spaces where only one part is going to be optimized. Thus, if in a
  * vector of dimension 20 you want to optimize the first 3 components then
  * i0=1, n=3; if you want to optimize it all, i0=1, n=20; and if you want to
  * optimize the last five components i0=15, n=5.
  *
  * The steps define the allowed steps on each variable. When a variable has
  * large ranges this step should be larger to allow bigger movements. The steps
  * should have only the dimension of the optimized part (3,20 and 5 in the
  * examples above).
  *
  * The option show forces the routine to show the convergence path
  *
  * If your function needs extra parameters you can pass them through
  * the void pointer prm. If you don't need this feature set it to NULL.
  *
  * Example of use:
  *
  * @code
  * MultidimArray<double> x(8), steps(8);
  * double fitness;
  * int iter;
  * steps.initConstant(1);
  * x.initZeros();
  * powellOptimizer(x,1,8,&wrapperFitness,NULL,0.01,fitness,iter,steps,true);
  * @endcode
  *
  */
void powellOptimizer(Matrix1D< double >& p,
                     int i0, int n,
                     double(*f)(double* , void *),
                     void *prm,
                     double ftol,
                     double& fret,
                     int& iter,
                     const Matrix1D< double >& steps,
                     bool show = false);

/** Gaussian interpolator
 * @ingroup NumericalTools
 *
 * This class helps to perform a quick evaluation of the N(0,1) Gaussian.
 * 1/sqrt(2*PI)*exp(-x^2)
 *
 * @code
 *  GaussianInterpolator GI;
 *  GI.initialize(6,60000);
 *  double g=GI.getValue(1.96);
 * @endcode
 */
class GaussianInterpolator {
    MultidimArray<double> v;
    double xstep;
    double xmax;
    double ixstep;
public:
    /** Constructor.
        xmax is the maximum value considered by the interpolator.
        N is the number of samples between 0 and xmax.
        If normalize is set to true, then the factor 1/sqrt(2*PI)
        is introduced. */
    void initialize(double xmax, int N, bool normalize=true);

    /** Value */
    inline double getValue(double x) const {
        x=fabs(x);
        if (x>xmax) return 0;
        else
        {
            int iaux=(int)round(x*ixstep);
            return DIRECT_A1D_ELEM(v,iaux);
        }
    }
};

/** Solve equation system, nonnegative solution
 * @ingroup NumericalTools
 *
 * The equation system is defined by Ax=b, it is solved for x. x is forced to be
 * nonnegative. It is designed to cope with large equation systems. This
 * function is borrowed from LAPACK nnls.
 *
 * The norm of the vector Ax-b is returned.
 */
double solveNonNegative(const Matrix2D< double >& A, const Matrix1D< double >& b,
                    Matrix1D< double >& result);

/** Solve equation system, symmetric positive-definite matrix
 * @ingroup NumericalTools
 *
 * The equation system is defined by Ax=b, it is solved for x. This method can
 * only be applied if A is positive-definite matrix and symmetric. It applies a
 * Cholesky factorization and backsubstitution (see Numerical Recipes).
 */
void solveViaCholesky(const Matrix2D< double >& A,
                        const Matrix1D< double >& b,
                        Matrix1D< double >& result);

/** Evaluate quadratic form
 * @ingroup NumericalTools
 *
 * Given x, c and H this function returns the value of the quadratic form
 * val=c^t*x+0.5*x^t*H^t*H*x and the gradient of the quadratic form at x
 * grad=c+H*x.
 *
 * Exceptions are thrown if the vectors and matrices do not have consistent
 * dimensions.
 */
void evaluateQuadratic(const Matrix1D< double >& x, const Matrix1D< double >& c,
                    const Matrix2D< double >& H, double& val,
                    Matrix1D< double >& grad);

/** Solves Quadratic programming subproblem
 * @ingroup NumericalTools
 *
 * @code
 * min 0.5*x'Cx + d'x   subject to:  A*x <= b
 *  x                                Aeq*x=beq
 *                                   bl<=x<=bu
 * @endcode
 */
void quadraticProgramming(const Matrix2D< double >& C, const Matrix1D< double >& d,
              const Matrix2D< double >& A, const Matrix1D< double >& b,
              const Matrix2D< double >& Aeq, const Matrix1D< double >& beq,
              Matrix1D< double >& bl, Matrix1D< double >& bu,
              Matrix1D< double >& x);


/** Solves the least square problem
 * @ingroup NumericalTools
 *
 * @code
 * min 0.5*(Norm(C*x-d))   subject to:  A*x <= b
 * x                                    Aeq*x=beq
 *                                      bl<=x<=bu
 * @endcode
 */
void leastSquare(const Matrix2D< double >& C, const Matrix1D< double >& d,
            const Matrix2D< double >& A, const Matrix1D< double >& b,
            const Matrix2D< double >& Aeq, const Matrix1D< double >& beq,
            Matrix1D< double >& bl, Matrix1D< double >& bu,
            Matrix1D< double >& x);

/** Solves the regularized least squares problem
 * @ingroup NumericalTools
 *
 * @code
 * min Norm(A*x-d) + lambda * norm (G*x)
 * @endcode
 *
 * Give an empty G matrix (NULL matrix) if G is the identity matrix
 * If AtA is not a NULL matrix, then the product AtA is not computed.
 */
void regularizedLeastSquare(const Matrix2D< double >& A,
    const Matrix1D< double >& d, double lambda,
    const Matrix2D< double >& G, Matrix1D< double >& x);

/** Solves a linear equation system by LU decomposition.
 * @ingroup NumericalTools
 */
template<typename T>
void solve(const Matrix2D<T>& A, const Matrix1D<T>& b, Matrix1D<T>& result)
{
    if (MAT_XSIZE(A) == 0)
        REPORT_ERROR(ERR_MATRIX_EMPTY, "Solve: Matrix is empty");

    if (MAT_XSIZE(A) != MAT_YSIZE(A))
        REPORT_ERROR(ERR_MATRIX_SIZE, "Solve: Matrix is not squared");

    if (MAT_XSIZE(A) != b.size())
        REPORT_ERROR(ERR_MATRIX_SIZE, "Solve: Different sizes of Matrix and Vector");

    if (b.isRow())
        REPORT_ERROR(ERR_MATRIX_DIM, "Solve: Not correct vector shape");

    // Perform LU decomposition and then solve
    Matrix1D< int > indx;
    T d;
    Matrix2D<T> LU;
    ludcmp(A, LU, indx, d);
    result = b;
    lubksb(LU, indx, result);
}

/** Solves a linear equation system by SVD decomposition.
 * @ingroup NumericalTools
 */
template<typename T>
void solveBySVD(const Matrix2D< T >& A, const Matrix1D< T >& b,
                  Matrix1D< double >& result, double tolerance)
{
    if (A.Xdim() == 0)
        REPORT_ERROR(ERR_MATRIX_EMPTY, "Solve: Matrix is empty");

    if (A.Xdim() != A.Ydim())
        REPORT_ERROR(ERR_MATRIX_SIZE, "Solve: Matrix is not squared");

    if (A.Xdim() != b.size())
        REPORT_ERROR(ERR_MATRIX_SIZE, "Solve: Different sizes of Matrix and Vector");

    if (b.isRow())
        REPORT_ERROR(ERR_MATRIX_DIM, "Solve: Not correct vector shape");

    // First perform de single value decomposition
    // Xmipp interface that calls to svdcmp of numerical recipes
    Matrix2D< double > u, v;
    Matrix1D< double > w;
    svdcmp(A, u, w, v);

    // Here is checked if eigenvalues of the svd decomposition are acceptable
    // If a value is lower than tolerance, the it's zeroed, as this increases
    // the precision of the routine.
    for (int i = 0; i < w.size(); i++)
    	if (w(i) < tolerance)
    		w(i) = 0;

    // Set size of matrices
    result.resize(b.size());

    // Xmipp interface that calls to svdksb of numerical recipes
    Matrix1D< double > bd;
    typeCast(b, bd);
    svbksb(u, w, v, bd, result);
}

/** Solves a linear equation system by Gaussian elimination.
 * @ingroup NumericalTools
 */
template<typename T>
void solve(const Matrix2D<T>& A, const Matrix2D<T>& b, Matrix2D<T>& result)
{
    if (A.Xdim() == 0)
        REPORT_ERROR(ERR_MATRIX_EMPTY, "Solve: Matrix is empty");

    if (A.Xdim() != A.Ydim())
        REPORT_ERROR(ERR_MATRIX_SIZE, "Solve: Matrix is not squared");

    if (A.Ydim() != b.Ydim())
        REPORT_ERROR(ERR_MATRIX_SIZE, "Solve: Different sizes of A and b");

    // Solve
    result = b;
    Matrix2D<T> Aux = A;
    gaussj(Aux.adaptForNumericalRecipes2(), Aux.Ydim(),
           result.adaptForNumericalRecipes2(), b.Xdim());
}

// Differential Evolution Solver Class
// Based on algorithms developed by Dr. Rainer Storn & Kenneth Price
// Written By: Lester E. Godwin
//             PushCorp, Inc.
//             Dallas, Texas
//             972-840-0208 x102
//             godwin@pushcorp.com
// Created: 6/8/98
// Last Modified: 6/8/98
// Revision: 1.0

#define stBest1Exp   0
#define stRand1Exp   1
#define stRandToBest1Exp 2
#define stBest2Exp   3
#define stRand2Exp   4
#define stBest1Bin   5
#define stRand1Bin   6
#define stRandToBest1Bin 7
#define stBest2Bin   8
#define stRand2Bin   9

class DESolver;

typedef void(DESolver::*StrategyFunction)(int);

/** Differential evolution class.
 *  For further information visit: http://www.icsi.berkeley.edu/~storn/code.html
  * Example of use:
  * @code
  * // Declaration of the optimizer
  * class AlignmentSolver : public DESolver {
  * public:
  *    AlignmentSolver(int dim,int pop) : DESolver(dim,pop), count(0) {;}
  *    double EnergyFunction(double trial[],bool &bAtSolution) {
  *       double result=wrapper_fitness_individual(trial);
  *       if (count++ % (5*nPop) == 0)
  *           std::cout << "Evaluations= " << count/nPop
  *                     << " energy= " << Energy()
  * 	            << " axisRot= " << Solution()[0]
  * 	            << " axisTilt= " << Solution()[1]
  * 	            << " axisShiftX= " << Solution()[2]
  * 	            << " axisShiftY= " << Solution()[3]
  * 	            << " axisShiftZ= " << Solution()[4]
  * 		    << std::endl;
  *       bAtSolution=false;
  *       return(result);
  *    }
  * private:
  *    int count;
  * };
  *
  * // Optimizer use:
  *        AlignmentSolver solver(length,length*Npop);
  *        solver.Setup(MULTIDIM_ARRAY(min_allowed),
  *                     MULTIDIM_ARRAY(max_allowed), stBest2Bin, 0.5, 0.8);
  *        solver.Solve(NGenerations);
  *        double current_energy=solver.Energy();
  *        double* bestSolution=solver.Solution();
  * @endcode
  */
class DESolver
{
public:
    /// Empty constructor
    DESolver(int dim, int popSize);

    /// Destructor
    virtual ~DESolver(void);

    /// Setup() must be called before solve to set min, max, strategy etc.
    void Setup(double min[],
               double max[],
               int deStrategy,
               double diffScale,
               double crossoverProb);

    /** Solve() returns true if EnergyFunction() returns true.
        Otherwise it runs maxGenerations generations and returns false.
    */
    virtual bool Solve(int maxGenerations);

    /** EnergyFunction must be overridden for problem to solve
        testSolution[] is nDim array for a candidate solution
        setting bAtSolution = true indicates solution is found
        and Solve() immediately returns true.
    */
    virtual double EnergyFunction(double testSolution[], bool& bAtSolution) = 0;

    /// Return dimension
    int Dimension() const
    {
        return (nDim);
    }

    /// Return population
    int Population() const
    {
        return (nPop);
    }

    /// Call these functions after Solve() to get results.
    double Energy() const
    {
        return (bestEnergy);
    }

    /// Return best solution
    double* Solution(void)
    {
        return (bestSolution);
    }

    /// Return the number of generations
    int Generations() const
    {
        return (generations);
    }

protected:
    void SelectSamples(int candidate,
                       int* r1,
                       int* r2 = 0,
                       int* r3 = 0,
                       int* r4 = 0,
                       int* r5 = 0);

    int nDim;
    int nPop;
    int generations;

    int strategy;
    StrategyFunction calcTrialSolution;
    double scale;
    double probability;

    double trialEnergy;
    double bestEnergy;

    double* trialSolution;
    double* bestSolution;
    double* popEnergy;
    double* population;
private:
    void Best1Exp(int candidate);
    void Rand1Exp(int candidate);
    void RandToBest1Exp(int candidate);
    void Best2Exp(int candidate);
    void Rand2Exp(int candidate);
    void Best1Bin(int candidate);
    void Rand1Bin(int candidate);
    void RandToBest1Bin(int candidate);
    void Best2Bin(int candidate);
    void Rand2Bin(int candidate);
};

/** Check the randomness of a sequence.
 *  The function returns a z-score of randomness. The highest the Z-score in absolute value,
 *  the less random the sequence is. The sequence is supposed to be formed by two symbols:
 *  + and -, 0 and 1, A and B, ...
 */
double checkRandomness(const std::string &sequence);
#endif
