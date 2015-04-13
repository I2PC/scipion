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
#include "numerical_tools.h"

/* Random permutation ------------------------------------------------------ */
void randomPermutation(int N, MultidimArray<int>& result)
{
    MultidimArray<double> aux;
    aux.resize(N);
    aux.initRandom(0,1);

    aux.indexSort(result);
    result-=1;
}

/* Powell's optimizer ------------------------------------------------------ */
void powellOptimizer(Matrix1D<double> &p, int i0, int n,
                     double(*f)(double *x, void *), void * prm,
                     double ftol, double &fret,
                     int &iter, const Matrix1D<double> &steps, bool show)
{
    double *xi = NULL;

    // Adapt indexes of p
    double *pptr = p.adaptForNumericalRecipes();
    double *auxpptr = pptr + (i0 - 1);

    // Form direction matrix
    ask_Tvector(xi, 1, n*n);
    for (int i = 1, ptr = 1; i <= n; i++)
        for (int j = 1; j <= n; j++, ptr++)
            xi[ptr] = (i == j) ? steps(i - 1) : 0;

    // Optimize
    xi -= n; // This is because NR works with matrices starting at [1,1]
    powell(auxpptr, xi, n, ftol, iter, fret, f, prm, show);
    xi += n;

    // Exit
    free_Tvector(xi, 1, n*n);
}

/* Gaussian interpolator -------------------------------------------------- */
void GaussianInterpolator::initialize(double _xmax, int N, bool normalize)
{
    xmax=_xmax;
    xstep=xmax/N;
    ixstep=1.0/xstep;
    v.initZeros(N);
    double inorm=1.0/sqrt(2*PI);
    FOR_ALL_ELEMENTS_IN_ARRAY1D(v)
    {
        double x=i*xstep;
        v(i)=exp(-x*x/2);
        if (normalize)
            v(i)*=inorm;
    }
}

/* Solve Cx=d, nonnegative x */
double solveNonNegative(const Matrix2D<double> &C, const Matrix1D<double> &d,
                        Matrix1D<double> &result)
{
    if (C.Xdim() == 0)
        REPORT_ERROR(ERR_MATRIX_EMPTY, "Solve_nonneg: Matrix is empty");
    if (C.Ydim() != d.size())
        REPORT_ERROR(ERR_MATRIX_SIZE, "Solve_nonneg: Different sizes of Matrix and Vector");
    if (d.isRow())
        REPORT_ERROR(ERR_MATRIX_DIM, "Solve_nonneg: Not correct vector shape");

    Matrix2D<double> Ct;
    Ct=C.transpose();

    result.initZeros(Ct.Ydim());
    double rnorm;

    // Watch out that matrix Ct is transformed.
    int success = nnls(MATRIX2D_ARRAY(Ct), Ct.Xdim(), Ct.Ydim(),
                       MATRIX1D_ARRAY(d),
                       MATRIX1D_ARRAY(result),
                       &rnorm, NULL, NULL, NULL);
    if (success == 1)
        std::cerr << "Warning, too many iterations in nnls\n";
    else if (success == 2)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "Solve_nonneg: Not enough memory");
    return rnorm;
}

/* Solve Ax=b, A definite positive and symmetric --------------------------- */
void solveViaCholesky(const Matrix2D<double> &A, const Matrix1D<double> &b,
                      Matrix1D<double> &result)
{
    Matrix2D<double> Ap = A;
    Matrix1D<double> p(A.Xdim());
    result.resize(A.Xdim());
    choldc(Ap.adaptForNumericalRecipes2(), A.Xdim(),
           p.adaptForNumericalRecipes());
    cholsl(Ap.adaptForNumericalRecipes2(), A.Xdim(),
           p.adaptForNumericalRecipes(), b.adaptForNumericalRecipes(),
           result.adaptForNumericalRecipes());
}


/* Quadratic form ---------------------------------------------------------- */
void evaluateQuadratic(const Matrix1D<double> &x, const Matrix1D<double> &c,
                       const Matrix2D<double> &H, double &val, Matrix1D<double> &grad)
{
    if (x.size() != c.size())
        REPORT_ERROR(ERR_MATRIX_SIZE, "Eval_quadratic: Not compatible sizes in x and c");
    if (H.Xdim() != x.size())
        REPORT_ERROR(ERR_MATRIX_SIZE, "Eval_quadratic: Not compatible sizes in x and H");

    // H*x, store in grad
    grad.initZeros(x.size());
    for (size_t i = 0; i < H.Ydim(); i++)
        for (size_t j = 0; j < x.size(); j++)
            grad(i) += H(i, j) * x(j);

    // Now, compute c^t*x+1/2*x^t*H*x
    // Add c to the gradient
    double quad = 0;
    val = 0;
    for (size_t j = 0; j < x.size(); j++)
    {
        quad += grad(j) * grad(j); // quad=x^t*H^t*H*x
        val += c(j) * x(j);  // val=c^t*x

        grad(j) += c(j);     // grad+=c
    }
    val += 0.5 * quad;
}

/* Quadprog and Lsqlin ----------------------------------------------------- */
/* Structure to pass the objective function and constraints to cfsqp*/
typedef struct
{
    Matrix2D<double> C;
    Matrix2D<double> D;
    Matrix2D<double> A;
    Matrix2D<double> B;
}
CDAB;

/*/////////////////////////////////////////////////////////////////////////
      Internal functions used by the quadraticProgramming function
/////////////////////////////////////////////////////////////////////////*/
/* To calculate the value of the objective function */
void quadraticProgramming_obj32(int nparam, int j, double* x, double* fj, void* cd)
{
    CDAB* in = (CDAB *)cd;
    Matrix2D<double> X(nparam,1);
    for (int i=0; i<nparam; ++i)
        X(i,0)=x[i];
    Matrix2D<double> result;
    result = 0.5 * X.transpose() * in->C * X + in->D.transpose() * X;

    *fj = result(0, 0);
}

/* To calculate the value of the jth constraint */
void quadraticProgramming_cntr32(int nparam, int j, double* x, double* gj, void* cd)
{
    CDAB* in = (CDAB *)cd;
    *gj = 0;
    for (int k = 0; k < nparam; k++)
        *gj += in->A(j - 1, k) * x[k];
    *gj -= in->B(j - 1, 0);
}

/* To calculate the value of the derivative of objective function */
void quadraticProgramming_grob32(int nparam, int j, double* x, double* gradfj, void(*mydummy)(int, int, double*, double*, void*), void *cd)
{
    CDAB* in = (CDAB *)cd;
    Matrix2D<double> X(1,nparam);
    for (int i=0; i<nparam; ++i)
        X(0,i)=x[i];

    Matrix2D<double> gradient;
    gradient = in->C * X + in->D;
    for (int k = 0; k < nparam; k++)
        gradfj[k] = gradient(k, 0);
}

/* To calculate the value of the derivative of jth constraint */
void quadraticProgramming_grcn32(int nparam, int j, double *x, double *gradgj, void(*mydummy)(int, int, double*, double*, void*), void *cd)
{
    CDAB* in = (CDAB *)cd;
    for (int k = 0; k < nparam; k++)
        gradgj[k] = in->A(j - 1, k);
}

/**************************************************************************

   Solves Quadratic programming subproblem.

  min 0.5*x'Cx + d'x   subject to:  A*x <= b
   x                                Aeq*x=beq
                                bl<=x<=bu

**************************************************************************/
void quadraticProgramming(const Matrix2D<double> &C, const Matrix1D<double> &d,
                          const Matrix2D<double> &A,   const Matrix1D<double> &b,
                          const Matrix2D<double> &Aeq, const Matrix1D<double> &beq,
                          Matrix1D<double> &bl,        Matrix1D<double> &bu,
                          Matrix1D<double> &x)
{
    CDAB prm;
    prm.C = C;
    prm.D.fromVector(d);
    prm.A.initZeros(A.Ydim() + Aeq.Ydim(), A.Xdim());
    prm.B.initZeros(prm.A.Ydim(), 1);


    // Copy Inequalities
    for (size_t i = 0; i < A.Ydim(); i++)
    {
        for (size_t j = 0; j < A.Xdim(); j++)
            prm.A(i, j) = A(i, j);
        prm.B(i, 0) = b(i);
    }

    // Copy Equalities
    for (size_t i = 0; i < Aeq.Ydim(); i++)
    {
        for (size_t j = 0; j < Aeq.Xdim(); j++)
            prm.A(i + A.Ydim(), j) = Aeq(i, j);
        prm.B(i + A.Ydim(), 0) = beq(i);
    }

    double bigbnd = 1e30;
    // Bounds
    if (bl.size() == 0)
    {
        bl.resize(C.Xdim());
        bl.initConstant(-bigbnd);
    }
    if (bu.size() == 0)
    {
        bu.resize(C.Xdim());
        bu.initConstant(bigbnd);
    }

    // Define intermediate variables
    int    mode = 100;  // CFSQP mode
    int    iprint = 0;  // Debugging
    int    miter = 1000;  // Maximum number of iterations
    double eps = 1e-4; // Epsilon
    double epsneq = 1e-4; // Epsilon for equalities
    double udelta = 0.e0; // Finite difference approximation
    // of the gradients. Not used in this function
    int    nparam = C.Xdim(); // Number of variables
    int    nf = 1;          // Number of objective functions
    int    neqn = Aeq.Ydim();        // Number of nonlinear equations
    int    nineqn = A.Ydim();      // Number of nonlinear inequations
    int    nineq = A.Ydim();  // Number of linear inequations
    int    neq = Aeq.Ydim();  // Number of linear equations
    int    inform;
    int    ncsrl = 0, ncsrn = 0, nfsr = 0, mesh_pts[] = {0};

    if (x.size() == 0)
        x.initZeros(nparam);
    Matrix1D<double> f(nf), g(nineq + neq), lambda(nineq + neq + nf + nparam);

    // Call the minimization routine
    cfsqp(nparam, nf, nfsr, nineqn, nineq, neqn, neq, ncsrl, ncsrn, mesh_pts,
          mode, iprint, miter, &inform, bigbnd, eps, epsneq, udelta,
          MATRIX1D_ARRAY(bl), MATRIX1D_ARRAY(bu),
          MATRIX1D_ARRAY(x),
          MATRIX1D_ARRAY(f), MATRIX1D_ARRAY(g),
          MATRIX1D_ARRAY(lambda),
          //  quadprg_obj32,quadprog_cntr32,quadprog_grob32,quadprog_grcn32,
          quadraticProgramming_obj32, quadraticProgramming_cntr32, grobfd, grcnfd,
          (void*)&prm);

#ifdef DEBUG

    if (inform == 0)
        std::cout << "SUCCESSFUL RETURN. \n";
    if (inform == 1 || inform == 2)
        std::cout << "\nINITIAL GUESS INFEASIBLE.\n";
    if (inform == 3)
        printf("\n MAXIMUM NUMBER OF ITERATIONS REACHED.\n");
    if (inform > 3)
        printf("\ninform=%d\n", inform);
#endif
}

/**************************************************************************

   Solves the least square problem

  min 0.5*(Norm(C*x-d))   subject to:  A*x <= b
   x                                   Aeq*x=beq
                                       bl<=x<=bu
**************************************************************************/
void leastSquare(const Matrix2D<double> &C, const Matrix1D<double> &d,
                 const Matrix2D<double> &A,   const Matrix1D<double> &b,
                 const Matrix2D<double> &Aeq, const Matrix1D<double> &beq,
                 Matrix1D<double> &bl,        Matrix1D<double> &bu,
                 Matrix1D<double> &x)
{
    // Convert d to Matrix2D for multiplication
    Matrix2D<double> P;
    P.fromVector(d);
    P = -2.0 * (P.transpose() * C);
    P = P.transpose();

    //Convert back to vector for passing it to quadraticProgramming
    Matrix1D<double> newd;
    P.toVector(newd);

    quadraticProgramming(C.transpose()*C, newd, A, b, Aeq, beq, bl, bu, x);
}

/* Regularized least squares ----------------------------------------------- */
void regularizedLeastSquare(const Matrix2D< double >& A,
                            const Matrix1D< double >& d, double lambda,
                            const Matrix2D< double >& G, Matrix1D< double >& x)
{
    int Nx=A.Xdim(); // Number of variables

    Matrix2D<double> X(Nx,Nx); // X=(A^t * A +lambda *G^t G)
    // Compute A^t*A
    FOR_ALL_ELEMENTS_IN_MATRIX2D(X)
    // Compute the dot product of the i-th and j-th columns of A
    for (size_t k=0; k<A.Ydim(); k++)
        X(i,j) += A(k,i) * A(k,j);

    // Compute lambda*G^t*G
    if (G.Xdim()==0)
        for (int i=0; i<Nx; i++)
            X(i,i)+=lambda;
    else
        FOR_ALL_ELEMENTS_IN_MATRIX2D(X)
        // Compute the dot product of the i-th and j-th columns of G
        for (size_t k=0; k<G.Ydim(); k++)
            X(i,j) += G(k,i) * G(k,j);

    // Compute A^t*d
    Matrix1D<double> Atd(Nx);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(Atd)
    // Compute the dot product of the i-th column of A and d
    for (size_t k=0; k<A.Ydim(); k++)
        Atd(i) += A(k,i) * d(k);

    // Compute the inverse of X
    Matrix2D<double> Xinv;
    X.inv(Xinv);

    // Now multiply Xinv * A^t * d
    x.initZeros(Nx);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(Xinv)
    x(i) += Xinv(i,j) * Atd(j);
}


////////DE_solver

#define Element(a,b,c)  a[b*nDim+c]
#define RowVector(a,b)  (&a[b*nDim])
#define CopyVector(a,b) memcpy((a),(b),nDim*sizeof(double))

DESolver::DESolver(int dim, int popSize) :
nDim(dim), nPop(popSize),
generations(0), strategy(stRand1Exp),
scale(0.7), probability(0.5), bestEnergy(0.0),
trialSolution(0), bestSolution(0),
popEnergy(0), population(0)
{
    trialSolution = new double[nDim];
    bestSolution  = new double[nDim];
    popEnergy   = new double[nPop];
    population   = new double[nPop * nDim];

    randomize_random_generator();
}

DESolver::~DESolver(void)
{
    if (trialSolution)
        delete [] trialSolution;
    if (bestSolution)
        delete [] bestSolution;
    if (popEnergy)
        delete [] popEnergy;
    if (population)
        delete [] population;

    trialSolution = bestSolution = popEnergy = population = NULL;
    return;
}

void DESolver::Setup(double min[], double max[],
                     int deStrategy, double diffScale, double crossoverProb)
{
    int i,j;

    strategy = deStrategy;
    scale  = diffScale;
    probability     = crossoverProb;
    bestEnergy      = 1.0E20;

    for (i = 0; i < nPop; i++)
    {
        for (j = 0; j < nDim; j++)
        {
            population[i*nDim+j] = rnd_unif(min[j], max[j]);
            //Element(population, i, j) = rnd_unif(min[j], max[j]);
        }

        popEnergy[i] = 1.0E20;
    }

    for (i = 0; i < nDim; i++)
        bestSolution[i] = (min[i] + max[i]) / 2.0;

    switch (strategy)
    {
    case stBest1Exp:
        calcTrialSolution = &DESolver::Best1Exp;
        break;

    case stRand1Exp:
        calcTrialSolution = &DESolver::Rand1Exp;
        break;

    case stRandToBest1Exp:
        calcTrialSolution = &DESolver::RandToBest1Exp;
        break;

    case stBest2Exp:
        calcTrialSolution = &DESolver::Best2Exp;
        break;

    case stRand2Exp:
        calcTrialSolution = &DESolver::Rand2Exp;
        break;

    case stBest1Bin:
        calcTrialSolution = &DESolver::Best1Bin;
        break;

    case stRand1Bin:
        calcTrialSolution = &DESolver::Rand1Bin;
        break;

    case stRandToBest1Bin:
        calcTrialSolution = &DESolver::RandToBest1Bin;
        break;

    case stBest2Bin:
        calcTrialSolution = &DESolver::Best2Bin;
        break;

    case stRand2Bin:
        calcTrialSolution = &DESolver::Rand2Bin;
        break;
    }

    return;
}

bool DESolver::Solve(int maxGenerations)
{
    bool bAtSolution = false;
    int generation;

    for (generation = 0;(generation < maxGenerations) && !bAtSolution;generation++)
        for (int candidate = 0; candidate < nPop; candidate++)
        {
            (this->*calcTrialSolution)(candidate);
            trialEnergy = EnergyFunction(trialSolution, bAtSolution);

            if (trialEnergy < popEnergy[candidate])
            {
                // New low for this candidate
                popEnergy[candidate] = trialEnergy;
                CopyVector(RowVector(population, candidate), trialSolution);

                // Check if all-time low
                if (trialEnergy < bestEnergy)
                {
                    bestEnergy = trialEnergy;
                    CopyVector(bestSolution, trialSolution);
                }
            }
        }

    generations = generation;
    return(bAtSolution);
}

void DESolver::Best1Exp(int candidate)
{
    int r1, r2;
    int n;

    SelectSamples(candidate, &r1, &r2);
    n = (int)rnd_unif(0.0, (double)nDim);

    CopyVector(trialSolution, RowVector(population, candidate));
    for (int i = 0; (rnd_unif(0.0, 1.0) < probability) && (i < nDim); i++)
    {
        trialSolution[n] = bestSolution[n]
                           + scale * (Element(population, r1, n)
                                      - Element(population, r2, n));
        n = (n + 1) % nDim;
    }

    return;
}

void DESolver::Rand1Exp(int candidate)
{
    int r1, r2, r3;
    int n;

    SelectSamples(candidate, &r1, &r2, &r3);
    n = (int)rnd_unif(0.0, (double)nDim);

    CopyVector(trialSolution, RowVector(population, candidate));
    for (int i = 0; (rnd_unif(0.0, 1.0) < probability) && (i < nDim); i++)
    {
        trialSolution[n] = Element(population, r1, n)
                           + scale * (Element(population, r2, n)
                                      - Element(population, r3, n));
        n = (n + 1) % nDim;
    }

    return;
}

void DESolver::RandToBest1Exp(int candidate)
{
    int r1, r2;
    int n;

    SelectSamples(candidate, &r1, &r2);
    n = (int)rnd_unif(0.0, (double)nDim);

    CopyVector(trialSolution, RowVector(population, candidate));
    for (int i = 0; (rnd_unif(0.0, 1.0) < probability) && (i < nDim); i++)
    {
        trialSolution[n] += scale * (bestSolution[n] - trialSolution[n])
                            + scale * (Element(population, r1, n)
                                       - Element(population, r2, n));
        n = (n + 1) % nDim;
    }

    return;
}

void DESolver::Best2Exp(int candidate)
{
    int r1, r2, r3, r4;
    int n;

    SelectSamples(candidate, &r1, &r2, &r3, &r4);
    n = (int)rnd_unif(0.0, (double)nDim);

    CopyVector(trialSolution, RowVector(population, candidate));
    for (int i = 0; (rnd_unif(0.0, 1.0) < probability) && (i < nDim); i++)
    {
        trialSolution[n] = bestSolution[n] +
                           scale * (Element(population, r1, n)
                                    + Element(population, r2, n)
                                    - Element(population, r3, n)
                                    - Element(population, r4, n));
        n = (n + 1) % nDim;
    }

    return;
}

void DESolver::Rand2Exp(int candidate)
{
    int r1, r2, r3, r4, r5;
    int n;

    SelectSamples(candidate, &r1, &r2, &r3, &r4, &r5);
    n = (int)rnd_unif(0.0, (double)nDim);

    CopyVector(trialSolution, RowVector(population, candidate));
    for (int i = 0; (rnd_unif(0.0, 1.0) < probability) && (i < nDim); i++)
    {
        trialSolution[n] = Element(population, r1, n)
                           + scale * (Element(population, r2, n)
                                      + Element(population, r3, n)
                                      - Element(population, r4, n)
                                      - Element(population, r5, n));
        n = (n + 1) % nDim;
    }

    return;
}

void DESolver::Best1Bin(int candidate)
{
    int r1, r2;
    int n;

    SelectSamples(candidate, &r1, &r2);
    n = (int)rnd_unif(0.0, (double)nDim);

    CopyVector(trialSolution, RowVector(population, candidate));
    for (int i = 0; i < nDim; i++)
    {
        if ((rnd_unif(0.0, 1.0) < probability) || (i == (nDim - 1)))
            trialSolution[n] = bestSolution[n]
                               + scale * (Element(population, r1, n)
                                          - Element(population, r2, n));
        n = (n + 1) % nDim;
    }

    return;
}

void DESolver::Rand1Bin(int candidate)
{
    int r1, r2, r3;
    int n;

    SelectSamples(candidate, &r1, &r2, &r3);
    n = (int)rnd_unif(0.0, (double)nDim);

    CopyVector(trialSolution, RowVector(population, candidate));
    for (int i = 0; i < nDim; i++)
    {
        if ((rnd_unif(0.0, 1.0) < probability) || (i  == (nDim - 1)))
            trialSolution[n] = Element(population, r1, n)
                               + scale * (Element(population, r2, n)
                                          - Element(population, r3, n));
        n = (n + 1) % nDim;
    }

    return;
}

void DESolver::RandToBest1Bin(int candidate)
{
    int r1, r2;
    int n;

    SelectSamples(candidate, &r1, &r2);
    n = (int)rnd_unif(0.0, (double)nDim);

    CopyVector(trialSolution, RowVector(population, candidate));
    for (int i = 0; i < nDim; i++)
    {
        if ((rnd_unif(0.0, 1.0) < probability) || (i  == (nDim - 1)))
            trialSolution[n] += scale * (bestSolution[n] - trialSolution[n])
                                + scale * (Element(population, r1, n)
                                           - Element(population, r2, n));
        n = (n + 1) % nDim;
    }

    return;
}

void DESolver::Best2Bin(int candidate)
{
    int r1, r2, r3, r4;
    int n;

    SelectSamples(candidate, &r1, &r2, &r3, &r4);
    n = (int)rnd_unif(0.0, (double)nDim);

    CopyVector(trialSolution, RowVector(population, candidate));
    for (int i = 0; i < nDim; i++)
    {
        if ((rnd_unif(0.0, 1.0) < probability) || (i  == (nDim - 1)))
            trialSolution[n] = bestSolution[n]
                               + scale * (Element(population, r1, n)
                                          + Element(population, r2, n)
                                          - Element(population, r3, n)
                                          - Element(population, r4, n));
        n = (n + 1) % nDim;
    }

    return;
}

void DESolver::Rand2Bin(int candidate)
{
    int r1, r2, r3, r4, r5;
    int n;

    SelectSamples(candidate, &r1, &r2, &r3, &r4, &r5);
    n = (int)rnd_unif(0.0, (double)nDim);

    CopyVector(trialSolution, RowVector(population, candidate));
    for (int i = 0; i < nDim; i++)
    {
        if ((rnd_unif(0.0, 1.0) < probability) || (i  == (nDim - 1)))
            trialSolution[n] = Element(population, r1, n)
                               + scale * (Element(population, r2, n)
                                          + Element(population, r3, n)
                                          - Element(population, r4, n)
                                          - Element(population, r5, n));
        n = (n + 1) % nDim;
    }

    return;
}

void DESolver::SelectSamples(int candidate, int *r1, int *r2,
                             int *r3, int *r4, int *r5)
{
    if (r1)
    {
        do
        {
            *r1 = (int)rnd_unif(0.0, (double)nPop);
        }
        while (*r1 == candidate);
    }

    if (r2)
    {
        do
        {
            *r2 = (int)rnd_unif(0.0, (double)nPop);
        }
        while ((*r2 == candidate) || (*r2 == *r1));
    }

    if (r3)
    {
        do
        {
            *r3 = (int)rnd_unif(0.0, (double)nPop);
        }
        while ((*r3 == candidate) || (*r3 == *r2) || (*r3 == *r1));
    }

    if (r4)
    {
        do
        {
            *r4 = (int)rnd_unif(0.0, (double)nPop);
        }
        while ((*r4 == candidate) || (*r4 == *r3) || (*r4 == *r2) || (*r4 == *r1));
    }

    if (r5)
    {
        do
        {
            *r5 = (int)rnd_unif(0.0, (double)nPop);
        }
        while ((*r5 == candidate) || (*r5 == *r4) || (*r5 == *r3)
               || (*r5 == *r2) || (*r5 == *r1));
    }

    return;
}

/* Check Randomness ------------------------------------------------------- */
/* See http://home.ubalt.edu/ntsbarsh/Business-stat/opre504.htm for the formulas */
/* This is called runs test */
double checkRandomness(const std::string &sequence)
{
    int imax=sequence.size();
    if (imax<=1)
        return 0;
    double n0=1, n1=0, R=1;
    int current=0;
    for (int i=1; i<imax; ++i)
    {
        if (sequence[i]!=sequence[i-1])
        {
        	current=(current+1)%2;
        	R++;
            if (current==0)
                n0++;
            else
                n1++;
        }
        else
        {
            if (current==0)
                n0++;
            else
                n1++;
        }
    }
    double m=1+2*n0*n1/(n0+n1);
    double s=sqrt(2*n0*n1*(2*n0*n1-n0-n1)/((n0+n1)*(n0+n1)*(n0+n1-1)));
    double z=(R-m)/s;
    return ABS(z);
}
