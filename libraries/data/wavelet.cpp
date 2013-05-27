/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Antonio Jose Rodriguez Sanchez (ajr@cnb.csic.es)
 *              Arun Kulshreshth        (arun_2000_iitd@yahoo.com)
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
/* WAVELETS                                                                  */
/* ------------------------------------------------------------------------- */

#include "wavelet.h"
#include "args.h"
#include "numerical_tools.h"
#include "histogram.h"
#include "mask.h"

#include "xmipp_image.h"
#include <external/bilib/headers/wavelet.h>
#include "xmipp_fftw.h"


/* Wavelet ----------------------------------------------------------------- */

void Bilib_DWT(const MultidimArray<double> &input,
               MultidimArray<double> &result, int iterations, int isign)
{
    if (iterations < 1)
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, "Bilib_DWT: iterations must be >=1");
    int size_multiple = (int)pow(2.0, (double) iterations);
    if (XSIZE(input) % size_multiple != 0)
        REPORT_ERROR(ERR_MULTIDIM_SIZE,
                     (std::string)"Bilib_DWT: Xsize must be a multiple of " +
                     integerToString(size_multiple));
    if (YSIZE(input) > 1 && YSIZE(input) % size_multiple != 0)
        REPORT_ERROR(ERR_MULTIDIM_SIZE,
                     (std::string)"Bilib_DWT 2D: Ysize must be a multiple of " +
                     integerToString(size_multiple));
    if (ZSIZE(input) > 1 && ZSIZE(input) % size_multiple != 0)
        REPORT_ERROR(ERR_MULTIDIM_SIZE,
                     (std::string)"Bilib_DWT 3D: Zsize must be a multiple of " +
                     integerToString(size_multiple));

    result.initZeros(input);
    TWaveletStruct TW;
    if (isign == 1)
        TW.Operation = "Analysis";
    else if (isign == -1)
        TW.Operation = "Synthesis";
    else
        REPORT_ERROR(ERR_VALUE_INCORRECT, (std::string)"waveletTransform: unrecognized isign");
    TW.Filter = "Orthonormal Spline";
    TW.BoundaryConditions = "Mirror";
    TW.Order = "1";
    TW.Alpha = 0;

    if (isign == 1)
    {
        // First iteration
        TW.Input = MULTIDIM_ARRAY(input);
        TW.NxOutput = TW.NxInput = XSIZE(input);
        TW.NyOutput = TW.NyInput = YSIZE(input);
        TW.NzOutput = TW.NzInput = ZSIZE(input);
        TW.Output = MULTIDIM_ARRAY(result);
        Wavelet(&TW);

        // Subsequent iterations
        for (int i = 1; i < iterations; i++)
        {
            // Size of the Lowest subband
            int xsize = XMIPP_MAX(1, XSIZE(input) / (int)pow(2.0, (double)i));
            int ysize = XMIPP_MAX(1, YSIZE(input) / (int)pow(2.0, (double)i));
            int zsize = XMIPP_MAX(1, ZSIZE(input) / (int)pow(2.0, (double)i));

            // Pick the Lowest subband
            MultidimArray<double> input_aux, result_aux;
            input_aux.resize(zsize, ysize, xsize);
            result_aux.resize(zsize, ysize, xsize);
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(input_aux)
                DIRECT_A3D_ELEM(input_aux, k, i, j) = DIRECT_A3D_ELEM(result, k, i, j);

            // DWT
            TW.Input = MULTIDIM_ARRAY(input_aux);
            TW.NxOutput = TW.NxInput = XSIZE(input_aux);
            TW.NyOutput = TW.NyInput = YSIZE(input_aux);
            TW.NzOutput = TW.NzInput = ZSIZE(input_aux);
            TW.Output = MULTIDIM_ARRAY(result_aux);
            Wavelet(&TW);

            // Return the subband to the output
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(input_aux)
                DIRECT_A3D_ELEM(result, k, i, j) = DIRECT_A3D_ELEM(result_aux, k, i, j);
        }
    }
    else if (isign == -1)
    {
        // Subsequent iterations
        for (int i = iterations - 1; i >= 1; i--)
        {
            // Size of the Lowest subband
            int xsize = XMIPP_MAX(1, XSIZE(input) / (int)pow(2.0, (double)i));
            int ysize = XMIPP_MAX(1, YSIZE(input) / (int)pow(2.0, (double)i));
            int zsize = XMIPP_MAX(1, ZSIZE(input) / (int)pow(2.0, (double)i));

            // Pick the Lowest subband
            MultidimArray<double> input_aux, result_aux;
            input_aux.resize(zsize, ysize, xsize);
            result_aux.resize(zsize, ysize, xsize);
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(input_aux)
                DIRECT_A3D_ELEM(input_aux, k, i, j) = DIRECT_A3D_ELEM(input, k, i, j);

            // DWT
            TW.Input = MULTIDIM_ARRAY(input_aux);
            TW.NxOutput = TW.NxInput = XSIZE(input_aux);
            TW.NyOutput = TW.NyInput = YSIZE(input_aux);
            TW.NzOutput = TW.NzInput = ZSIZE(input_aux);
            TW.Output = MULTIDIM_ARRAY(result_aux);
            Wavelet(&TW);

            // Return the subband to the output
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(input_aux)
                DIRECT_A3D_ELEM(input, k, i, j) = DIRECT_A3D_ELEM(result_aux, k, i, j);
        }

        // First iteration
        TW.Input = MULTIDIM_ARRAY(input);
        TW.NxOutput = TW.NxInput = XSIZE(input);
        TW.NyOutput = TW.NyInput = YSIZE(input);
        TW.NzOutput = TW.NzInput = ZSIZE(input);
        TW.Output = MULTIDIM_ARRAY(result);
        Wavelet(&TW);
    }
}

// Set the DWT type --------------------------------------------------------
void set_DWT_type(int DWT_type)
{
    pwtset(DWT_type);
}

void IDWT(const MultidimArray<double> &v, MultidimArray<double> &result)
{
    DWT(v, result, -1);
}

// Lowpass DWT -------------------------------------------------------------
void DWT_lowpass2D(const MultidimArray<double> &v, MultidimArray<double> &result)
{
    MultidimArray<double> dwt, aux;
    result.initZeros(YSIZE(v), XSIZE(v) / 2);
    DWT(v, dwt);
    int Nx = Get_Max_Scale(XSIZE(v));
    for (int s = 0; s < Nx; s++)
    {
        // Perform the inverse DWT transform of the low pass
        dwt.resize(XSIZE(dwt) / 2, YSIZE(dwt) / 2);
        IDWT(dwt, aux);
        // Copy the result to the 01 quadrant of the result
        int x1, y1, x2, y2, x, y, i, j;
        SelectDWTBlock(s, v, "01", x1, x2, y1, y2);
        for (y = y1, i = 0; y <= y2; y++, i++)
            for (x = x1, j = 0; x <= x2; x++, j++)
                A2D_ELEM(result, y, x) = A2D_ELEM(aux, i, j);
    }
}

// Select block ------------------------------------------------------------
// Quadrant .---------------------------------------------------------------
std::string Quadrant2D(int q)
{
    switch (q)
    {
    case 0:
        return "00";
        break;
    case 1:
        return "01";
        break;
    case 2:
        return "10";
        break;
    case 3:
        return "11";
        break;
    default:
    	REPORT_ERROR(ERR_ARG_INCORRECT,"Unknown quadrant");
    }
}

std::string Quadrant3D(int q)
{
    switch (q)
    {
    case 0:
        return "000";
        break;
    case 1:
        return "001";
        break;
    case 2:
        return "010";
        break;
    case 3:
        return "011";
        break;
    case 4:
        return "100";
        break;
    case 5:
        return "101";
        break;
    case 6:
        return "110";
        break;
    case 7:
        return "111";
        break;
    default:
    	REPORT_ERROR(ERR_ARG_INCORRECT,"Unknown quadrant");
    }
}

// Provide block -----------------------------------------------------------
#define DWT_Scale(i,smax) ((int)((i==0)?smax-1:(ABS((CEIL(log10((double)(i+1))/log10(2.0))-smax)))))
#define DWT_Quadrant1D(i,s,smax) ((s!=smax-1)?'1':((i==0)?'0':'1'))
#define DWT_QuadrantnD(i,s,sp,smax) \
    ((s!=sp)?'0':DWT_Quadrant1D(i,s,smax))
#define DWT_icoef1D(i,s,smax) ()

void Get_Scale_Quadrant(int size_x, int x,
                        int &scale, std::string &quadrant)
{
    double Nx = Get_Max_Scale(size_x);
    quadrant = "x";
    scale = DWT_Scale(x, Nx);
    quadrant[0] = DWT_Quadrant1D(x, scale, Nx);
}

void Get_Scale_Quadrant(int size_x, int size_y, int x, int y,
                        int &scale, std::string &quadrant)
{
    double Nx = Get_Max_Scale(size_x);
    double Ny = Get_Max_Scale(size_y);
    quadrant = "xy";
    double scalex = DWT_Scale(x, Nx);
    double scaley = DWT_Scale(y, Ny);
    scale = (int)(XMIPP_MIN(scalex, scaley));
    quadrant[1] = DWT_QuadrantnD(y, scaley, scale, Ny);
    quadrant[0] = DWT_QuadrantnD(x, scalex, scale, Nx);
}

void Get_Scale_Quadrant(int size_x, int size_y, int size_z,
                        int x, int y, int z,
                        int &scale, std::string &quadrant)
{
    double Nx = Get_Max_Scale(size_x);
    double Ny = Get_Max_Scale(size_y);
    double Nz = Get_Max_Scale(size_z);
    quadrant = "xyz";
    double scalex = DWT_Scale(x, Nx);
    double scaley = DWT_Scale(y, Ny);
    double scalez = DWT_Scale(z, Nz);
    scale = (int)(XMIPP_MIN(scalez, XMIPP_MIN(scalex, scaley)));
    quadrant[2] = DWT_QuadrantnD(z, scalez, scale, Nz);
    quadrant[1] = DWT_QuadrantnD(y, scaley, scale, Ny);
    quadrant[0] = DWT_QuadrantnD(x, scalex, scale, Nx);
}

// Clean quadrant ----------------------------------------------------------
void clean_quadrant2D(MultidimArray<double> &I, int scale, const std::string &quadrant)
{
    Matrix1D<int> corner1(2), corner2(2);
    Matrix1D<double> r(2);
    SelectDWTBlock(scale, I, quadrant, XX(corner1), XX(corner2), YY(corner1), YY(corner2));
    FOR_ALL_ELEMENTS_IN_ARRAY2D_BETWEEN(corner1, corner2)
        I(r) = 0;
}

void clean_quadrant3D(MultidimArray<double> &I, int scale, const std::string &quadrant)
{
    int x1, y1, z1, x2, y2, z2;
    SelectDWTBlock(scale, I, quadrant, x1, x2, y1, y2, z1, z2);
    Matrix1D<int> corner1(3), corner2(3);
    Matrix1D<double> r(3);
    SelectDWTBlock(scale, I, quadrant, XX(corner1), XX(corner2),
                   YY(corner1), YY(corner2), ZZ(corner1), ZZ(corner2));
    FOR_ALL_ELEMENTS_IN_ARRAY3D_BETWEEN(corner1, corner2) I(r) = 0;
}

void soft_thresholding(MultidimArray<double> &I, double th)
{
    FOR_ALL_ELEMENTS_IN_ARRAY3D(I)
        if (ABS(A3D_ELEM(I, k, i, j)) > th)
        if (A3D_ELEM(I, k, i, j) > 0)
            A3D_ELEM(I, k, i, j) -= th;
        else
            A3D_ELEM(I, k, i, j) += th;
    else
        A3D_ELEM(I, k, i, j) = 0;
}

// Adaptive soft thresholding ----------------------------------------------
void adaptive_soft_thresholding_block2D(MultidimArray<double> &I, int scale,
                                        const std::string &quadrant, double sigma)
{
    // Compute block variance
    Matrix1D<int> corner1(2), corner2(2);
    Matrix1D<double> r(2);
    SelectDWTBlock(scale, I, quadrant,
                   XX(corner1), XX(corner2), YY(corner1), YY(corner2));
    double dummy, avg, stddev;
    I.computeStats(avg, stddev, dummy, dummy, corner1, corner2);

    // Now denoise
    double th = sigma * sigma / stddev;
    FOR_ALL_ELEMENTS_IN_ARRAY2D_BETWEEN(corner1, corner2)
    {
        if (ABS(I(r)) > th)
            if (I(r) > 0)
                I(r) -= th;
            else
                I(r) += th;
        else
            I(r) = 0;
    }
}

double compute_noise_power(MultidimArray<double> &I)
{
    // Compute histogram of the absolute values of the DWT coefficients
    // at scale=0
    Histogram1D hist;
    double avg, stddev, min_val, max_val;
    I.computeStats(avg, stddev, min_val, max_val);
    hist.init(0, XMIPP_MAX(ABS(min_val), ABS(max_val)), 100);

    Matrix1D<int> corner1(2), corner2(2);
    Matrix1D<double> r(2);
    SelectDWTBlock(0, I, "01",
                   XX(corner1), XX(corner2), YY(corner1), YY(corner2));
    FOR_ALL_ELEMENTS_IN_ARRAY2D_BETWEEN(corner1, corner2)
    {
    	double value=fabs(I(r));
    	INSERT_VALUE(hist,value);
    }

    SelectDWTBlock(0, I, "10",
                   XX(corner1), XX(corner2), YY(corner1), YY(corner2));
    FOR_ALL_ELEMENTS_IN_ARRAY2D_BETWEEN(corner1, corner2)
    {
    	double value=fabs(I(r));
    	INSERT_VALUE(hist,value);
    }

    SelectDWTBlock(0, I, "11",
                   XX(corner1), XX(corner2), YY(corner1), YY(corner2));
    FOR_ALL_ELEMENTS_IN_ARRAY2D_BETWEEN(corner1, corner2)
    {
    	double value=fabs(I(r));
    	INSERT_VALUE(hist,value);
    }

    return hist.percentil(50) / 0.6745;
}

void adaptive_soft_thresholding2D(MultidimArray<double> &I, int scale)
{
    double sigma = compute_noise_power(I);
    for (int s = 0; s <= scale; s++)
    {
        adaptive_soft_thresholding_block2D(I, s, "01", sigma);
        adaptive_soft_thresholding_block2D(I, s, "10", sigma);
        adaptive_soft_thresholding_block2D(I, s, "11", sigma);
    }
}

// Keep central part -------------------------------------------------------
void DWT_keep_central_part(MultidimArray<double> &I, double R)
{
    Mask mask(INT_MASK);
    mask.type = BINARY_DWT_CIRCULAR_MASK;
    if (R == -1)
        mask.R1 = (double)XSIZE(I) / 2 + 1;
    else
        mask.R1 = R;
    mask.smin = 0;
    mask.smax = Get_Max_Scale(XSIZE(I));
    mask.quadrant = "xx";
    mask.resize(I);
    mask.generate_mask();
    mask.apply_mask(I, I);
}

// Bayesian Wiener filtering -----------------------------------------------
//DWT_Bijaoui_denoise_LL -- Bijaoui denoising at a perticular scale.
void DWT_Bijaoui_denoise_LL2D(MultidimArray<double> &WI, int scale,
                            const std::string &orientation,
                            double mu, double S, double N)
{
    Matrix1D<int> x0(2), xF(2), r(2);
    SelectDWTBlock(scale, WI, orientation, XX(x0), XX(xF), YY(x0), YY(xF));

    double SN = S + N;
    double S_N = S / SN;
    if (S < 1e-6 && N < 1e-6)
    {
        FOR_ALL_ELEMENTS_IN_ARRAY2D_BETWEEN(x0, xF)
		A2D_ELEM(WI,YY(r),XX(r)) = 0;
    }
    else
    {
    	double iSN=1.0/SN;
    	double iN=1.0/N;
        FOR_ALL_ELEMENTS_IN_ARRAY2D_BETWEEN(x0, xF)
        {
        	double &WI_r=A2D_ELEM(WI,YY(r),XX(r));
            double y = WI_r;
            double ymu = y - mu;
            double ymu2 = -0.5 * ymu * ymu;
            double expymu2SN = exp(ymu2 * iSN);
            double den = exp(ymu2 * iN) + expymu2SN;
            if (den > 1e-10)
                WI_r *= S_N * expymu2SN / den;
        }
    }
}

void DWT_Bijaoui_denoise_LL3D(MultidimArray<double> &WI, int scale,
                            const std::string &orientation,
                            double mu, double S, double N)
{
    Matrix1D<int> x0(3), xF(3), r(3);
    SelectDWTBlock(scale, WI, orientation, XX(x0), XX(xF), YY(x0), YY(xF),
                   ZZ(x0), ZZ(xF));

    double SN = S + N;
    double S_N = S / SN;
    if (S < 1e-6 && N < 1e-6)
    {
        FOR_ALL_ELEMENTS_IN_ARRAY3D_BETWEEN(x0, xF) WI(r) = 0;
    }
    else
    {
        FOR_ALL_ELEMENTS_IN_ARRAY3D_BETWEEN(x0, xF)
        {
            double y = WI(r);
            double ymu = y - mu;
            double ymu2 = -0.5 * ymu * ymu;
            double expymu2SN = exp(ymu2 / SN);
            double den = exp(ymu2 / N) + expymu2SN;
            if (den > 1e-10)
                WI(r) = S_N * expymu2SN / den * y;
        }
    }
}

//#define DEBUG
void bayesian_solve_eq_system(
    const Matrix1D<double> &power,
    const Matrix1D<double> &average,
    const Matrix1D<double> &Ncoefs,
    double SNR0,
    double SNRF,
    double powerI,
    double power_rest,
    bool white_noise,
    int tell,
    Matrix1D<double> &estimatedS)
{

    int scale_dim = power.size();

    Matrix2D<double> A;
    int extra_constraints = 0;
    if (white_noise)
        extra_constraints += (scale_dim - 1);
    A.initZeros(2*(scale_dim - 1) + 2*scale_dim + 2 + extra_constraints, 2*scale_dim);
    for (int i = 1;i < scale_dim;i++)
    {
        A(i - 1, i - 1) = 1;
        A(i - 1, i) = -1;
        A(i - 1 + scale_dim - 1, i - 1 + scale_dim) = 1;
        A(i - 1 + scale_dim - 1, i + scale_dim) = -1;
    }
    for (int i = 0;i < 2*scale_dim;i++)
        A(i + 2*(scale_dim - 1), i) = -1;

    // Constraints on the SNR
    Matrix1D<double> aux0coefs(scale_dim);
    for (int j = 0;j < scale_dim;j++)
        aux0coefs(j) = Ncoefs(j) * SNR0;
    Matrix1D<double> auxFcoefs(scale_dim);
    for (int j = 0;j < scale_dim;j++)
        auxFcoefs(j) = Ncoefs(j) * SNRF;

    //initializing the second last row of A
    for (int j = 0;j < scale_dim;j++)
        A(2*(scale_dim - 1) + 2*scale_dim, j) = (-1) * auxFcoefs(j);
    for (int j = scale_dim;j < 2*scale_dim;j++)
        A(2*(scale_dim - 1) + 2*scale_dim, j) = Ncoefs(j - scale_dim);

    //initializing the last row of A
    for (int j = 0;j < scale_dim;j++)
        A(2*(scale_dim - 1) + 2*scale_dim + 1, j) = aux0coefs(j);
    for (int j = scale_dim;j < 2*scale_dim;j++)
        A(2*(scale_dim - 1) + 2*scale_dim + 1, j) = (-1) * Ncoefs(j - scale_dim);

    // White noise constraints
    if (white_noise)
        for (int i = 0; i < scale_dim - 1; i++)
        {
            A(A.mdimy - (scale_dim - 1) + i, i)  = -1.01;
            A(A.mdimy - (scale_dim - 1) + i, i + 1) = 1;
        }

    //initialize the matrix b
    Matrix1D<double> b(MAT_YSIZE(A));

    // Initialize Aeq matrix
    Matrix2D<double> Aeq;
    Aeq.initZeros(1, 2*scale_dim);
    for (int j = 0;j < scale_dim;j++)
    {
        Aeq(0, j) = Ncoefs(j);
        Aeq(0, j + scale_dim) = Ncoefs(j);
    }

    //initialize beq matrix
    Matrix1D<double> beq;
    beq.initZeros(1);
    beq(0) = powerI - power_rest;

    //initialization of Matrix C (cost matrix)
    Matrix2D<double> C;
    C.initZeros(scale_dim, 2*scale_dim);
    for (int j = 0;j < scale_dim;j++)
    {
        C(j, j) = 1;
        C(j, j + scale_dim) = 1;
    }

    // initialise the estimatedS which will contain the solution vector
    estimatedS.initZeros(2*scale_dim);
#ifdef DEBUG
    //Writing the matrices to ASCII files for comparing with matlab
    C.write("./matrices/C.txt");
    power.write("./matrices/power.txt");
    A.write("./matrices/A.txt");
    b.write("./matrices/b.txt");
    Aeq.write("./matrices/Aeq.txt");
    beq.write("./matrices/beq.txt");

    std::cout << "Equation system Cx=d\n"
    << "C=\n" << C << std::endl
    << "d=" << (power / Ncoefs).transpose() << std::endl
    << "Constraints\n"
    << "Ax<=b\n"
    << "A=\n" << A << std::endl
    << "b=" << b.transpose() << std::endl
    << "Aeq x=beq\n"
    << "Aeq=\n" << Aeq << std::endl
    << "beq=" << beq.transpose() << std::endl;
#endif

    // Solve the system
    Matrix1D<double> bl, bu;
    leastSquare(C, power / Ncoefs, A, b, Aeq, beq, bl, bu, estimatedS);
    // COSS
    estimatedS /= 2;

#ifdef DEBUG

    std::cout << "scale_dim :: " << scale_dim << std::endl;
    std::cout << "--------estimatedS -------- \n";
    std::cout << estimatedS;
    std::cout << "--------------------------- \n";
    std::cout << "Inequality constraints agreement" << std::endl
    << (A*estimatedS).transpose() << std::endl;
    std::cout << "Equality constraints agreement" << std::endl
    << (Aeq*estimatedS).transpose() << std::endl;
    std::cout << "Goal function value: " << (C*estimatedS).transpose() << std::endl;
#endif
}
#undef DEBUG

//#define DEBUG
Matrix1D<double> bayesian_wiener_filtering2D(MultidimArray<double> &WI, int allowed_scale,
        double SNR0, double SNRF, bool white_noise, int tell, bool denoise)
{
    /*Calculate the power of the wavelet transformed image */
    double powerI = WI.sum2();

    /*Number of pixels and some constraints on SNR*/
    int Xdim = XSIZE(WI);
    int max_scale = ROUND(log(double(Xdim)) / log(2.0));

#ifdef DEBUG

    std::cout << "powerI= " << powerI << std::endl;
    double powerWI = WI.sum2();
    std::cout << "powerWI= " << powerWI << std::endl;
    std::cout << "Ydim = " << Ydim << "  Xdim = " << Xdim << "\n";
    std::cout << "Ncoef_total= " << Ncoef_total << std::endl;
    std::cout << "max_scale = " << max_scale << "\n";
#endif

    /*Calculate the power at each band*/
    //init the scale vector
    Matrix1D<int> scale(XMIPP_MIN(allowed_scale + 1, max_scale - 1));
    FOR_ALL_ELEMENTS_IN_MATRIX1D(scale) scale(i) = i;
    int scale_dim = scale.size();

    //define some vectors
    Matrix1D<double> power(scale_dim), average(scale_dim), Ncoefs(scale_dim);
    Matrix1D<int> x0(2), xF(2), r(2);
    std::vector<std::string> orientation;
    orientation.push_back("01");
    orientation.push_back("10");
    orientation.push_back("11");
    int orientationSize=orientation.size();
    int jmax=scale.size();
    for (int j = 0;j < jmax;j++)
    {
        for (size_t k = 0; k < orientation.size(); k++)
        {
            SelectDWTBlock(scale(j), WI, orientation[k],
                           XX(x0), XX(xF), YY(x0), YY(xF));
            FOR_ALL_ELEMENTS_IN_ARRAY2D_BETWEEN(x0, xF)
            {
            	double aux=DIRECT_A2D_ELEM(WI,YY(r),XX(r));
                VEC_ELEM(power,j) += aux*aux;
                VEC_ELEM(average,j) += aux;
            }
        }
        VEC_ELEM(Ncoefs,j) = (int)pow(2.0, 2 * (max_scale - VEC_ELEM(scale,j) - 1)) * orientationSize;
        VEC_ELEM(average,j) = VEC_ELEM(average,j) / VEC_ELEM(Ncoefs,j);
    }

    /*Evaluate the power of the unconsidered part of the image */
    double power_rest = 0.0;
    int Ncoefs_rest = 0;
    SelectDWTBlock(scale(scale_dim - 1), WI, "00", XX(x0), XX(xF), YY(x0), YY(xF));
    FOR_ALL_ELEMENTS_IN_ARRAY2D_BETWEEN(x0, xF)
    {
    	double aux=DIRECT_A2D_ELEM(WI,YY(r),XX(r));
    	power_rest += aux*aux;
    }
    Ncoefs_rest = (int)pow(2.0, 2 * (max_scale - 1 - scale(scale_dim - 1)));

    if (tell)
    {
        std::cout << "scale = " << std::endl << scale << std::endl;
        std::cout << "power= " << std::endl << power << "\n";
        std::cout << "average= " << std::endl << average << "\n";
        std::cout << "Ncoefs= " << std::endl << Ncoefs << "\n";
        std::cout << "power_rest= " << power_rest << "\n";
        std::cout << "Ncoefs_rest= " << Ncoefs_rest << "\n";
        std::cout << "powerI= " << powerI << std::endl;
        std::cout << "Total sum of powers = " << power.sum() + power_rest << std::endl;
        std::cout << "Difference powers = " << powerI - power.sum() - power_rest << std::endl;
    }

    /*Solve the Equation System*/
    Matrix1D<double> estimatedS;
    bayesian_solve_eq_system(power, average, Ncoefs,
                             SNR0, SNRF, powerI, power_rest, white_noise, tell, estimatedS);

    if (tell)
    {
        std::cout << "estimatedS =\n" << estimatedS << std::endl;
        double S = 0, N = 0;
        for (int i = 0; i < scale_dim; i++)
        {
            N += Ncoefs(i) * estimatedS(i);
            S += Ncoefs(i) * estimatedS(scale_dim + i);
        }
        std::cout << "SNR value=" << S / N << std::endl << std::endl;
    }

    /* Apply the Bijaoui denoising to all scales >= allowed_scale */
    if (denoise)
        bayesian_wiener_filtering2D(WI, allowed_scale, estimatedS);

    return estimatedS;
}
#undef DEBUG

void bayesian_wiener_filtering2D(MultidimArray<double> &WI,
                               int allowed_scale, Matrix1D<double> &estimatedS)
{
    std::vector<std::string> orientation;
    orientation.push_back("01");
    orientation.push_back("10");
    orientation.push_back("11");

    int max_scale = ROUND(log(double(XSIZE(WI))) / log(2.0));
    Matrix1D<int> scale(XMIPP_MIN(allowed_scale + 1, max_scale - 1));
    FOR_ALL_ELEMENTS_IN_MATRIX1D(scale) scale(i) = i;

    for (size_t i = 0;i < scale.size();i++)
    {
        double N = estimatedS(i);
        double S = estimatedS(i + scale.size());
        for (size_t k = 0; k < orientation.size(); k++)
            DWT_Bijaoui_denoise_LL2D(WI, scale(i), orientation[k], 0, S, N);
    }
}

//#define DEBUG
Matrix1D<double> bayesian_wiener_filtering3D(MultidimArray<double> &WI, int allowed_scale,
                                                  double SNR0, double SNRF, bool white_noise,
                                                  int tell, bool denoise)
{
    /*Calculate the power of the wavelet transformed image */
    double powerI = WI.sum2();

    /*Number of pixels and some constraints on SNR*/
    size_t Xdim=ZSIZE(WI);
    int max_scale = ROUND(log(double(Xdim)) / log(2.0));

#ifdef DEBUG

    std::cout << "powerI= " << powerI << std::endl;
    double powerWI = WI.sum2();
    std::cout << "powerWI= " << powerWI << std::endl;
    std::cout << "Zdim= " << Zdim << " Ydim = " << Ydim << "  Xdim = " << Xdim << "\n";
    std::cout << "Ncoef_total= " << Ncoef_total << std::endl;
    std::cout << "max_scale = " << max_scale << "\n";
#endif

    /*Calculate the power at each band*/
    //init the scale vector
    Matrix1D<int> scale(allowed_scale + 1);
    FOR_ALL_ELEMENTS_IN_MATRIX1D(scale) scale(i) = i;
    int scale_dim = scale.size();

    //define some vectors
    Matrix1D<double> power(scale_dim), average(scale_dim), Ncoefs(scale_dim);
    Matrix1D<int> x0(3), xF(3), r(3);
    std::vector<std::string> orientation;
    orientation.push_back("001");
    orientation.push_back("010");
    orientation.push_back("011");
    orientation.push_back("100");
    orientation.push_back("101");
    orientation.push_back("110");
    orientation.push_back("111");
    for (size_t j = 0;j < scale.size();j++)
    {
        for (size_t k = 0; k < orientation.size(); k++)
        {
            SelectDWTBlock(scale(j), WI, orientation[k],
                           XX(x0), XX(xF), YY(x0), YY(xF), ZZ(x0), ZZ(xF));
            FOR_ALL_ELEMENTS_IN_ARRAY3D_BETWEEN(x0, xF)
            {
                power(j) += WI(r) * WI(r);
                average(j) += WI(r);
            }
        }
        Ncoefs(j) = (int)pow(2.0, 3 * (max_scale - scale(j) - 1)) * orientation.size();
        average(j) = average(j) / Ncoefs(j);
    }

    /*Evaluate the power of the unconsidered part of the image */
    double power_rest = 0.0;
    int Ncoefs_rest = 0;
    SelectDWTBlock(scale(scale_dim - 1), WI, "000", XX(x0), XX(xF), YY(x0), YY(xF),
                   ZZ(x0), ZZ(xF));
    FOR_ALL_ELEMENTS_IN_ARRAY3D_BETWEEN(x0, xF)
    power_rest += WI(r) * WI(r);
    Ncoefs_rest = (int)pow(2.0, 3 * (max_scale - 1 - scale(scale_dim - 1)));

    if (tell)
    {
        std::cout << "scale = " << std::endl << scale << std::endl;
        std::cout << "power= " << std::endl << power << "\n";
        std::cout << "average= " << std::endl << average << "\n";
        std::cout << "Ncoefs= " << std::endl << Ncoefs << "\n";
        std::cout << "power_rest= " << power_rest << "\n";
        std::cout << "Ncoefs_rest= " << Ncoefs_rest << "\n";
        std::cout << "powerI= " << powerI << std::endl;
        std::cout << "Total sum of powers = " << power.sum() + power_rest << std::endl;
        std::cout << "Difference powers = " << powerI - power.sum() - power_rest << std::endl;
    }

    /*Solve the Equation System*/
    Matrix1D<double> estimatedS;
    bayesian_solve_eq_system(power, average, Ncoefs,
                             SNR0, SNRF, powerI, power_rest, white_noise, tell, estimatedS);
    if (tell)
    {
        std::cout << "estimatedS =\n" << estimatedS << std::endl;
        double S = 0, N = 0;
        for (int i = 0; i < scale_dim; i++)
        {
            N += Ncoefs(i) * estimatedS(i);
            S += Ncoefs(i) * estimatedS(scale_dim + i);
        }
        std::cout << "SNR value=" << S / N << std::endl << std::endl;
    }

    /* Apply the Bijaoui denoising to all scales >= allowed_scale */
    if (denoise)
        bayesian_wiener_filtering3D(WI, allowed_scale, estimatedS);
    return estimatedS;
}
#undef DEBUG

void bayesian_wiener_filtering3D(MultidimArray<double> &WI,
                               int allowed_scale, Matrix1D<double> &estimatedS)
{
    std::vector<std::string> orientation;
    orientation.push_back("001");
    orientation.push_back("010");
    orientation.push_back("011");
    orientation.push_back("100");
    orientation.push_back("101");
    orientation.push_back("110");
    orientation.push_back("111");

    int max_scale = ROUND(log(double(XSIZE(WI))) / log(2.0));
    MultidimArray<int> scale(XMIPP_MIN(allowed_scale + 1, max_scale - 1));
    FOR_ALL_ELEMENTS_IN_ARRAY1D(scale) scale(i) = i;

    for (size_t i = 0;i < XSIZE(scale);i++)
    {
        double N = estimatedS(i);
        double S = estimatedS(i + XSIZE(scale));
        for (size_t k = 0; k < orientation.size(); k++)
            DWT_Bijaoui_denoise_LL3D(WI, scale(i), orientation[k], 0, S, N);
    }
}

void phaseCongMono(MultidimArray< double >& I,
				   MultidimArray< double >& Ph,
				   MultidimArray< double >& Or,
				   MultidimArray< double >& Energy,
				   MultidimArray<double>& lowPass,
				   MultidimArray<double>& Radius,
				   MultidimArray< std::complex <double> >& H,
                                 int nScale,
                                 double minWaveLength,
                                 double mult,
                                 double sigmaOnf)
{

//    #define DEBUG
	double epsilon= .0001; // Used to prevent division by zero.
	//First we set the image origin in the image center
    I.setXmippOrigin();
	//Image size
    size_t NR, NC,NZ, NDim;
    I.getDimensions(NC,NR,NZ,NDim);

    if ( (NZ!=1) || (NDim!=1) )
        REPORT_ERROR(ERR_MULTIDIM_DIM,(std::string)"ZDim and NDim has to be equals to one");

    MultidimArray< std::complex <double> > fftIm;
    // Fourier Transformer
    FourierTransformer ftrans(FFTW_BACKWARD);//, ftransh(FFTW_BACKWARD), ftransf(FFTW_BACKWARD);
    ftrans.FourierTransform(I, fftIm, false);

    if ( (lowPass.xinit == 0) & (lowPass.yinit == 0) & (Radius.xinit == 0) & (Radius.yinit == 0)
    	  & (H.xinit == 0) & (H.yinit == 0) )
	{

    	H.resizeNoCopy(I);
        lowPass.resizeNoCopy(I);
        Radius.resizeNoCopy(I);

        double cutoff = .4;
        double n = 10;
        double wx,wy;

        for (int i=STARTINGY(I); i<=FINISHINGY(I); ++i)
        {
            FFT_IDX2DIGFREQ(i,YSIZE(I),wy);
            double wy2=wy*wy;
            for (int j=STARTINGX(I); j<=FINISHINGX(I); ++j)
            {
            	FFT_IDX2DIGFREQ(j,XSIZE(I),wx);

            	double wx2=wx*wx;
            	double radius=sqrt(wy2+wx2);
            	if (radius < 1e-10) radius=1;
            	A2D_ELEM(Radius,i,j)=radius;
            	double *ptr=(double*)&A2D_ELEM(H,i,j);
            	*ptr=wy/radius;
            	*(ptr+1)=wx/radius;
            	A2D_ELEM(lowPass,i,j)= 1.0/(1.0+std::pow(radius/cutoff,n));
            }

        }

	}

    MultidimArray<double> logGabor;
    logGabor.resizeNoCopy(I);
    //Bandpassed image in the frequency domain and in the spatial domain.
    //h is a Bandpassed monogenic filtering, real part of h contains
    //convolution result with h1, imaginary part
    //contains convolution result with h2.
    MultidimArray< std::complex <double> > fftImF, fullFftIm, f,Hc,h;

    fftImF.resizeNoCopy(I);
    Hc.resizeNoCopy(I);
    f.resizeNoCopy(I);
    h.resizeNoCopy(I);

    MultidimArray<double> An,AnTemp;
    MultidimArray<double> temp;
    MultidimArray<double> F;
    MultidimArray<double> h1;
    MultidimArray<double> h2;
    MultidimArray<double> tempMat;

    temp.resizeNoCopy(I);
    An.resizeNoCopy(I);
    F.resizeNoCopy(I);
    h1.resizeNoCopy(I);
    h2.resizeNoCopy(I);
    AnTemp.resizeNoCopy(I);
    ftrans.getCompleteFourier(fullFftIm);
    CenterFFT(fullFftIm,false);

    for (int num = 0; num < nScale; ++num)
    {
    	double waveLength = minWaveLength;

    	for(int numMult=0; numMult < num;numMult++)
    		waveLength*=mult;

    	double fo = 1.0/waveLength;
    	FOR_ALL_ELEMENTS_IN_ARRAY2D(fullFftIm)
    	{
    		double temp1 =(std::log(A2D_ELEM(Radius,i,j)/fo));
    		double temp2 = std::log(sigmaOnf);
    		A2D_ELEM(logGabor,i,j) = std::exp(-(temp1*temp1) /(2 * temp2*temp2))*A2D_ELEM(lowPass,i,j);
    		if (A2D_ELEM(Radius,i,j)<=1e-10) (A2D_ELEM(logGabor,i,j)=0);

    		double *ptr= (double*)&A2D_ELEM(fullFftIm,i,j);
    		temp1 = *ptr;
    		temp2 = *(ptr+1);
    		ptr= (double*)&A2D_ELEM(fftImF,i,j);
    		*ptr=temp1*A2D_ELEM(logGabor,i,j);
    		*(ptr+1)=temp2*A2D_ELEM(logGabor,i,j);
    	}

    	CenterFFT(fftImF,true);
    	Hc = fftImF*H;
    	ftrans.inverseFourierTransform(fftImF,f);
    	ftrans.inverseFourierTransform(Hc,h);

    	//ftransf.setFourier(fftImF);
    	//ftransf.inverseFourierTransform();
    	//ftransh.setFourier(Hc);
    	//ftransh.inverseFourierTransform();

    	f.getReal(temp);
    	F += temp;
    	AnTemp = temp*temp;

    	h.getReal(temp);
    	h1 += temp;
    	AnTemp += temp*temp;

    	h.getImag(temp);
    	h2 += temp;
    	AnTemp += temp*temp;

    	AnTemp.selfSQRT();
    	An += AnTemp;
    }

    Or.resizeNoCopy(I);
	Or.setXmippOrigin();

	Ph.resizeNoCopy(I);
	Ph.setXmippOrigin();

	Energy.resizeNoCopy(I);
	Energy.setXmippOrigin();

	FOR_ALL_ELEMENTS_IN_ARRAY2D(I)
  	{
   		double temph1 = A2D_ELEM(h1,i,j);
   		double temph2 = A2D_ELEM(h2,i,j);
   		double tempF =  A2D_ELEM(F,i,j);

   		A2D_ELEM(Or,i,j) = std::atan2(temph1,temph2);
   		A2D_ELEM(Ph,i,j)=  std::atan2(tempF, std::sqrt(temph1*temph1+
   						   temph2*temph2));
   		A2D_ELEM(Energy,i,j) = std::sqrt(tempF*tempF+temph1*temph1
   							 + temph2*temph2)+epsilon;
  	}

#ifdef DEBUG

    FileName fpName1    = "test1.txt";
    FileName fpName2    = "test2.txt";
    FileName fpName3    = "test3.txt";
    FileName fpName4    = "test4.txt";

    Or.write(fpName1);
    Ph.write(fpName2);
    Energy.write(fpName3);

    #endif
}
