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

/* ------------------------------------------------------------------------- */
/* MATRICES                                                                  */
/* ------------------------------------------------------------------------- */
#include "matrix2d.h"

/* Interface to numerical recipes: svbksb ---------------------------------- */
void svbksb(Matrix2D<double> &u, Matrix1D<double> &w, Matrix2D<double> &v,
            Matrix1D<double> &b, Matrix1D<double> &x)
{
    // Call to the numerical recipes routine. Results will be stored in X
    svbksb(u.adaptForNumericalRecipes2(),
           w.adaptForNumericalRecipes(),
           v.adaptForNumericalRecipes2(),
           u.rowNumber(), u.colNumber(),
           b.adaptForNumericalRecipes(),
           x.adaptForNumericalRecipes());
}

/* Solve Cx=d, nonnegative x */
double solveNonNegative(const Matrix2D<double> &C, const Matrix1D<double> &d,
                    Matrix1D<double> &result)
{
    if (XSIZE(C) == 0)
        REPORT_ERROR(1108, "Solve_nonneg: Matrix is empty");
    if (YSIZE(C) != XSIZE(d))
        REPORT_ERROR(1102, "Solve_nonneg: Different sizes of Matrix and Vector");
    if (d.isRow())
        REPORT_ERROR(1107, "Solve_nonneg: Not correct vector shape");

    Matrix2D<double> Ct(XSIZE(C), YSIZE(C)); // Ct=C^t
    FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(Ct)
        DIRECT_MAT_ELEM(Ct, i, j) = DIRECT_MAT_ELEM(C, j, i);

    result.initZeros(YSIZE(Ct));
    double rnorm;

    // Watch out that matrix Ct is transformed.
    int success = nnls(MULTIDIM_ARRAY(Ct), XSIZE(Ct), YSIZE(Ct),
                       MULTIDIM_ARRAY(d),
                       MULTIDIM_ARRAY(result),
                       &rnorm, NULL, NULL, NULL);
    if (success == 1)
        std::cerr << "Warning, too many iterations in nnls\n";
    else if (success == 2)
        REPORT_ERROR(1, "Solve_nonneg: Not enough memory");
    return rnorm;
}

/* Solve Ax=b, A definite positive and symmetric --------------------------- */
void solveViaCholesky(const Matrix2D<double> &A, const Matrix1D<double> &b,
                        Matrix1D<double> &result)
{
    Matrix2D<double> Ap = A;
    Matrix1D<double> p(XSIZE(A));
    result.resize(XSIZE(A));
    choldc(Ap.adaptForNumericalRecipes2(), XSIZE(A),
           p.adaptForNumericalRecipes());
    cholsl(Ap.adaptForNumericalRecipes2(), XSIZE(A),
           p.adaptForNumericalRecipes(), b.adaptForNumericalRecipes(),
           result.adaptForNumericalRecipes());
}

// Special case for complex numbers
template <>
void applyGeometryBSpline(Matrix2D< std::complex<double> > &M2,
                        const Matrix2D<double> &A, const Matrix2D< std::complex<double> > &M1,
                        int Splinedegree, bool inv, bool wrap, std::complex<double> outside)
{
    Matrix2D<double> re, im, rotre, rotim;
    double outre, outim;
    re.resize(YSIZE(M1), XSIZE(M1));
    im.resize(YSIZE(M1), XSIZE(M1));
    outre = outside.real();
    outim = outside.imag();
    Complex2RealImag(MULTIDIM_ARRAY(M1),
                     MULTIDIM_ARRAY(re), MULTIDIM_ARRAY(im),
                     MULTIDIM_SIZE(M1));
    applyGeometryBSpline(rotre, A, re, Splinedegree, inv, wrap, outre);
    applyGeometryBSpline(rotim, A, im, Splinedegree, inv, wrap, outim);
    M2.resize(M1);
    RealImag2Complex(MULTIDIM_ARRAY(rotre), MULTIDIM_ARRAY(rotim),
                     MULTIDIM_ARRAY(M2), MULTIDIM_SIZE(re));
}


/* Is diagonal ------------------------------------------------------------- */
template <>
bool Matrix2D< std::complex<double> >::isDiagonal() const
{
    if (XSIZE(*this) != YSIZE(*this))
        return false;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(*this)
        if (i != j && abs(DIRECT_MAT_ELEM(*this, i, j)) > XMIPP_EQUAL_ACCURACY)
            return false;
    return true;
}

/* Is Scalar --------------------------------------------------------------- */
template <>
bool Matrix2D< std::complex<double> >::isScalar() const
{
    if (!isDiagonal())
        return false;
    for (int i = 1; i < YSIZE(*this); i++)
        if (abs(DIRECT_MAT_ELEM(*this, i, i) - DIRECT_MAT_ELEM(*this, 0, 0)) >
            XMIPP_EQUAL_ACCURACY)
            return false;
    return true;
}

/* Rotation 2D ------------------------------------------------------------- */
Matrix2D<double> rotation2DMatrix(double ang)
{
    Matrix2D<double> result(3, 3);
    double cosine, sine;

    ang = DEG2RAD(ang);
    cosine = cos(ang);
    sine = sin(ang);

    DIRECT_MAT_ELEM(result, 0, 0) = cosine;
    DIRECT_MAT_ELEM(result, 0, 1) = -sine;
    DIRECT_MAT_ELEM(result, 0, 2) = 0;

    DIRECT_MAT_ELEM(result, 1, 0) = sine;
    DIRECT_MAT_ELEM(result, 1, 1) = cosine;
    DIRECT_MAT_ELEM(result, 1, 2) = 0;

    DIRECT_MAT_ELEM(result, 2, 0) = 0;
    DIRECT_MAT_ELEM(result, 2, 1) = 0;
    DIRECT_MAT_ELEM(result, 2, 2) = 1;

    return result;
}

/* Translation 2D ---------------------------------------------------------- */
Matrix2D<double> translation2DMatrix(Matrix1D<double> v)
{
    if (XSIZE(v) != 2)
        REPORT_ERROR(1002, "Translation2D_matrix: vector is not in R2");

    Matrix2D<double> result(3, 3);

    result.initIdentity();
    DIRECT_MAT_ELEM(result, 0, 2) = XX(v);
    DIRECT_MAT_ELEM(result, 1, 2) = YY(v);

    return result;
}

/* Rotation 3D around the system axes -------------------------------------- */
Matrix2D<double> rotation3DMatrix(double ang, char axis)
{
    Matrix2D<double> result(4, 4);
    double cosine, sine;

    ang = DEG2RAD(ang);
    cosine = cos(ang);
    sine = sin(ang);

    result.initZeros();
    DIRECT_MAT_ELEM(result, 3, 3) = 1;
    switch (axis)
    {
    case 'Z':
        DIRECT_MAT_ELEM(result, 0, 0) = cosine;
        DIRECT_MAT_ELEM(result, 0, 1) = -sine;
        DIRECT_MAT_ELEM(result, 1, 0) = sine;
        DIRECT_MAT_ELEM(result, 1, 1) = cosine;
        DIRECT_MAT_ELEM(result, 2, 2) = 1;
        break;
    case 'Y':
        DIRECT_MAT_ELEM(result, 0, 0) = cosine;
        DIRECT_MAT_ELEM(result, 0, 2) = -sine;
        DIRECT_MAT_ELEM(result, 2, 0) = sine;
        DIRECT_MAT_ELEM(result, 2, 2) = cosine;
        DIRECT_MAT_ELEM(result, 1, 1) = 1;
        break;
    case 'X':
        DIRECT_MAT_ELEM(result, 1, 1) = cosine;
        DIRECT_MAT_ELEM(result, 1, 2) = -sine;
        DIRECT_MAT_ELEM(result, 2, 1) = sine;
        DIRECT_MAT_ELEM(result, 2, 2) = cosine;
        DIRECT_MAT_ELEM(result, 0, 0) = 1;
        break;
    default:
        REPORT_ERROR(1105, "rotation3DMatrix: Unknown axis");
    }
    return result;
}

/* Align a vector with Z axis */
Matrix2D<double> alignWithZ(const Matrix1D<double> &axis)
{
    Matrix1D<double>  Axis;
    Matrix2D<double>  A(4, 4);

    if (XSIZE(axis) != 3)
        REPORT_ERROR(1002, "alignWithZ: Axis is not in R3");

    // Copy axis and compute length of the projection on YZ plane
    Axis = axis;
    Axis.selfNormalize();
    double proj_mod = sqrt(YY(Axis) * YY(Axis) + ZZ(Axis) * ZZ(Axis));

    A(3, 3) = 1;
    if (proj_mod > XMIPP_EQUAL_ACCURACY)
    { // proj_mod!=0
        // Build Matrix A, which makes the turning axis coincident with Z
        A(0, 0) = proj_mod;
        A(0, 1) = -XX(Axis) * YY(Axis) / proj_mod;
        A(0, 2) = -XX(Axis) * ZZ(Axis) / proj_mod;
        A(1, 0) = 0;
        A(1, 1) = ZZ(Axis) / proj_mod;
        A(1, 2) = -YY(Axis) / proj_mod;
        A(2, 0) = XX(Axis);
        A(2, 1) = YY(Axis);
        A(2, 2) = ZZ(Axis);
    }
    else
    {
        // I know that the Axis is the X axis
        A(0, 0) = 0;
        A(0, 1) = 0;
        A(0, 2) = -1;
        A(1, 0) = 0;
        A(1, 1) = 1;
        A(1, 2) = 0;
        A(2, 0) = 1;
        A(2, 1) = 0;
        A(2, 2) = 0;
    }
    return A;
}

/* Rotation 3D around any axis -------------------------------------------- */
Matrix2D<double> rotation3DMatrix(double ang, const Matrix1D<double> &axis)
{
    // Compute a matrix which makes the turning axis coincident with Z
    // And turn around this axis
    Matrix2D<double> A = alignWithZ(axis);
    return A.transpose() * rotation3DMatrix(ang, 'Z') * A;
}

/* Translation 3D ---------------------------------------------------------- */
Matrix2D<double> translation3DMatrix(const Matrix1D<double> &v)
{
    if (XSIZE(v) != 3)
        REPORT_ERROR(1002, "Translation3D_matrix: vector is not in R3");

    Matrix2D<double> result(4, 4);

    result.initIdentity();
    DIRECT_MAT_ELEM(result, 0, 3) = XX(v);
    DIRECT_MAT_ELEM(result, 1, 3) = YY(v);
    DIRECT_MAT_ELEM(result, 2, 3) = ZZ(v);

    return result;
}

/* Scale 3D ---------------------------------------------------------------- */
Matrix2D<double> scale3DMatrix(const Matrix1D<double> &sc)
{
    if (XSIZE(sc) != 3)
        REPORT_ERROR(1002, "Scale3D_matrix: vector is not in R3");

    Matrix2D<double> result(4, 4);

    result.initIdentity();
    DIRECT_MAT_ELEM(result, 0, 0) = XX(sc);
    DIRECT_MAT_ELEM(result, 1, 1) = YY(sc);
    DIRECT_MAT_ELEM(result, 2, 2) = ZZ(sc);

    return result;
}

// Show a complex matrix ---------------------------------------------------
std::ostream& operator<<(std::ostream& ostrm,
    const Matrix2D< std::complex<double> >& v)
{
    if (XSIZE(v) == 0 || YSIZE(v) == 0)
        ostrm << "NULL matrix\n";
    else
    {
        for (int i = STARTINGY(v); i <= FINISHINGY(v); i++)
        {
            for (int j = STARTINGX(v); j <= FINISHINGX(v); j++)
                ostrm << MAT_ELEM(v, i, j) << ' ';
            ostrm << std::endl;
        }
    }

    return ostrm;
}

/* Quadratic form ---------------------------------------------------------- */
void evaluateQuadratic(const Matrix1D<double> &x, const Matrix1D<double> &c,
                    const Matrix2D<double> &H, double &val, Matrix1D<double> &grad)
{
    if (XSIZE(x) != XSIZE(c))
        REPORT_ERROR(1102, "Eval_quadratic: Not compatible sizes in x and c");
    if (XSIZE(H) != XSIZE(x))
        REPORT_ERROR(1102, "Eval_quadratic: Not compatible sizes in x and H");

    // H*x, store in grad
    grad.initZeros(XSIZE(x));
    for (int i = 0; i < YSIZE(H); i++)
        for (int j = 0; j < XSIZE(x); j++)
            DIRECT_VEC_ELEM(grad, i) += DIRECT_MAT_ELEM(H, i, j) *
                                        DIRECT_VEC_ELEM(x, j);

    // Now, compute c^t*x+1/2*x^t*H*x
    // Add c to the gradient
    double quad = 0;
    val = 0;
    for (int j = 0; j < XSIZE(x); j++)
    {
        quad += DIRECT_VEC_ELEM(grad, j) * DIRECT_VEC_ELEM(grad, j); // quad=x^t*H^t*H*x
        val += DIRECT_VEC_ELEM(c, j) * DIRECT_VEC_ELEM(x, j);  // val=c^t*x

        DIRECT_VEC_ELEM(grad, j) += DIRECT_VEC_ELEM(c, j);     // grad+=c
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
    prm.A.initZeros(YSIZE(A) + YSIZE(Aeq), XSIZE(A));
    prm.B.initZeros(YSIZE(prm.A), 1);


    // Copy Inequalities
    for (int i = 0; i < YSIZE(A); i++)
    {
        for (int j = 0; j < XSIZE(A); j++)
            prm.A(i, j) = A(i, j);
        prm.B(i, 0) = b(i);
    }

    // Copy Equalities
    for (int i = 0; i < YSIZE(Aeq); i++)
    {
        for (int j = 0; j < XSIZE(Aeq); j++)
            prm.A(i + YSIZE(A), j) = Aeq(i, j);
        prm.B(i + YSIZE(A), 0) = beq(i);
    }

    double bigbnd = 1e30;
    // Bounds
    if (XSIZE(bl) == 0)
    {
        bl.resize(XSIZE(C));
        bl.initConstant(-bigbnd);
    }
    if (XSIZE(bu) == 0)
    {
        bu.resize(XSIZE(C));
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
    int    nparam = XSIZE(C); // Number of variables
    int    nf = 1;          // Number of objective functions
    int    neqn = YSIZE(Aeq);        // Number of nonlinear equations
    int    nineqn = YSIZE(A);      // Number of nonlinear inequations
    int    nineq = YSIZE(A);  // Number of linear inequations
    int    neq = YSIZE(Aeq);  // Number of linear equations
    int    inform;
    int    ncsrl = 0, ncsrn = 0, nfsr = 0, mesh_pts[] = {0};

    if (XSIZE(x) == 0)
        x.initZeros(nparam);
    Matrix1D<double> f(nf), g(nineq + neq), lambda(nineq + neq + nf + nparam);

    // Call the minimization routine
    cfsqp(nparam, nf, nfsr, nineqn, nineq, neqn, neq, ncsrl, ncsrn, mesh_pts,
          mode, iprint, miter, &inform, bigbnd, eps, epsneq, udelta,
          MULTIDIM_ARRAY(bl), MULTIDIM_ARRAY(bu),
          MULTIDIM_ARRAY(x),
          MULTIDIM_ARRAY(f), MULTIDIM_ARRAY(g),
          MULTIDIM_ARRAY(lambda),
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
    P = -2 * P.transpose() * C;
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
    int Nd=YSIZE(A); // Number of data samples
    int Nx=XSIZE(A); // Number of variables

    Matrix2D<double> X(Nx,Nx); // X=(A^t * A +lambda *G^t G)
    // Compute A^t*A
    FOR_ALL_ELEMENTS_IN_MATRIX2D(X)
        // Compute the dot product of the i-th and j-th columns of A
        for (int k=0; k<YSIZE(A); k++)
            DIRECT_MAT_ELEM(X,i,j)+=
                DIRECT_MAT_ELEM(A,k,i)*DIRECT_MAT_ELEM(A,k,j);

    // Compute lambda*G^t*G
    if (XSIZE(G)==0)
        for (int i=0; i<Nx; i++)
            DIRECT_MAT_ELEM(X,i,i)+=lambda;
    else
        FOR_ALL_ELEMENTS_IN_MATRIX2D(X)
            // Compute the dot product of the i-th and j-th columns of G
            for (int k=0; k<YSIZE(G); k++)
                DIRECT_MAT_ELEM(X,i,j)+=
                    DIRECT_MAT_ELEM(G,k,i)*DIRECT_MAT_ELEM(G,k,j);

    // Compute A^t*d
    Matrix1D<double> Atd(Nx);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(Atd)
        // Compute the dot product of the i-th column of A and d
        for (int k=0; k<YSIZE(A); k++)
            DIRECT_VEC_ELEM(Atd,i)+=
                DIRECT_MAT_ELEM(A,k,i)*DIRECT_VEC_ELEM(d,k);

    // Compute the inverse of X
    Matrix2D<double> Xinv;
    X.inv(Xinv);

    // Now multiply Xinv * A^t * d
    x.initZeros(Nx);
    FOR_ALL_ELEMENTS_IN_MATRIX2D(Xinv)
        DIRECT_VEC_ELEM(x,i)+=DIRECT_MAT_ELEM(Xinv,i,j)*
            DIRECT_VEC_ELEM(Atd,j);
}
