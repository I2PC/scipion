/***************************************************************************
 *
 * Authors:    Slavica Jonic                slavica.jonic@a3.epfl.ch (2004)
 *             Carlos Oscar Sanchez Sorzano coss.eps@ceu.es
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

#include "angular_continuous_assign.h"

#include <bilib/configs.h>
#include <bilib/headers/linearalgebra.h>
#include <bilib/headers/error.h>
#include <bilib/headers/messagedisplay.h>
#include <bilib/headers/getputd.h>
#include <bilib/headers/kernel.h>
#include <bilib/headers/kerneldiff1.h>
#include <bilib/types/tsplinebasis.h>
#include <bilib/types/tboundaryconvention.h>
#include <bilib/headers/changebasis.h>
#include <bilib/headers/dft.h>

#include <data/xmipp_funcs.h>
#include <data/args.h>
#include <data/xmipp_fft.h>
#include <data/mask.h>

/* ------------------------------------------------------------------------- */
// Prototypes
int cstregistration(struct cstregistrationStruct *Data);

// Empty constructor =======================================================
ProgAngularContinuousAssign::ProgAngularContinuousAssign()
{
    produces_a_metadata = true;
}

// Read arguments ==========================================================
void ProgAngularContinuousAssign::readParams()
{
	XmippMetadataProgram::readParams();
    fn_ref = getParam("--ref");
    gaussian_DFT_sigma = getDoubleParam("--gaussian_Fourier");
    gaussian_Real_sigma = getDoubleParam("--gaussian_Real");
    weight_zero_freq = getDoubleParam("--zerofreq_weight");
    max_no_iter = getIntParam("--max_iter");
    max_shift = getDoubleParam("--max_shift");
    max_angular_change = getDoubleParam("--max_angular_change");
}

// Show ====================================================================
void ProgAngularContinuousAssign::show()
{
    if (!verbose)
        return;
	XmippMetadataProgram::show();
    std::cout << "Reference volume:    " << fn_ref              << std::endl
    << "Gaussian Fourier:    " << gaussian_DFT_sigma  << std::endl
    << "Gaussian Real:       " << gaussian_Real_sigma << std::endl
    << "Zero-frequency weight:"<< weight_zero_freq    << std::endl
    << "Max. Iter:           " << max_no_iter         << std::endl
    << "Max. Shift:          " << max_shift           << std::endl
    << "Max. Angular Change: " << max_angular_change  << std::endl
    ;
}

// usage ===================================================================
void ProgAngularContinuousAssign::defineParams()
{
    addUsageLine("Make a continuous angular assignment");
    addUsageLine("+This program assigns Euler angles to experimental projections by ");
    addUsageLine("+minimizing the difference between the given experimental image and ");
    addUsageLine("+the central slice of a reference volume in Fourier space. All ");
    addUsageLine("+interpolations are based on a B-spline model of the images and the ");
    addUsageLine("+volume. The translations are also optimized. Since an iterative ");
    addUsageLine("+optimization is performed, it must be initialized with some rough ");
    addUsageLine("+estimate of the angles. The output of xmipp_angular_predict can be ");
    addUsageLine("+used as initialization without any transformation.");
    addUsageLine("+The method is fully described at http://www.ncbi.nlm.nih.gov/pubmed/15885434");
	defaultComments["-i"].clear();
	defaultComments["-i"].addComment("Metadata with initial alignment");
	defaultComments["-o"].clear();
	defaultComments["-o"].addComment("Metadata with output alignment");
    XmippMetadataProgram::defineParams();
    addParamsLine("   --ref <volume>              : Reference volume");
    addParamsLine("  [--gaussian_Fourier <s=0.5>] : Weighting sigma in Fourier space");
    addParamsLine("                               :+Small values of this parameter concentrate ");
    addParamsLine("                               :+the optimization on low frequencies. This ");
    addParamsLine("                               :+value should be below 1.");
    addParamsLine("  [--gaussian_Real    <s=0.5>] : Weighting sigma in Real space");
    addParamsLine("                               :+Small values of this parameter concentrate ");
    addParamsLine("                               :+the optimization in the image center. This ");
    addParamsLine("                               :+value should be below 1.");
    addParamsLine("  [--zerofreq_weight  <s=0. >] : Zero-frequency weight");
    addParamsLine("  [--max_iter <max=60>]        : Maximum number of iterations");
    addParamsLine("                               :+Convergence might be reached before this number of iterations");
    addParamsLine("  [--max_shift <s=-1>]         : Maximum shift allowed");
    addParamsLine("                               :+Use this option to limit the maximum shift that the ");
    addParamsLine("                               :+optimization process can find. If the minimum is ");
    addParamsLine("                               :+found beyond this limit (note it is an absolute limit ");
    addParamsLine("                               :+on the shift and it is not relative to the initial shift), ");
    addParamsLine("                               :+then the solution given by the initial guess is kept ");
    addParamsLine("                               :+since the optimized one looks suspicious.");
    addParamsLine("  [--max_angular_change <a=-1>]: Maximum angular change allowed");
    addParamsLine("                               :+Use this option to limit the maximum angular ");
    addParamsLine("                               :+change with respect to the initial solution ");
    addParamsLine("                               :+that the optimization algorithm can find. If ");
    addParamsLine("                               :+the solution found is beyond this limit, then ");
    addParamsLine("                               :+the initial solution is returned instead of the ");
    addParamsLine("                               :+optimized one since this latter looks suspicious.");
    addExampleLine("A typical use is:",false);
    addExampleLine("xmipp_angular_continuous_assign -i anglesFromDiscreteAssignment.doc --ref reference.vol -o assigned_angles.xmd");
}

// Produce side information ================================================
void ProgAngularContinuousAssign::preProcess()
{
    // Read the reference volume
    Image<double> V;
    V.read(fn_ref);
    V().setXmippOrigin();

    // Prepare the masks in real space
    Mask mask_Real3D;
    mask_Real3D.type = mask_Real.type = GAUSSIAN_MASK;
    mask_Real3D.mode = mask_Real.mode = INNER_MASK;

    mask_Real3D.sigma = mask_Real.sigma = gaussian_Real_sigma * ((double)XSIZE(V()));

    mask_Real3D.generate_mask(V());
    mask_Real.generate_mask(YSIZE(V()), XSIZE(V()));

    double gs2 = 2. * PI * gaussian_Real_sigma * gaussian_Real_sigma * ((double)XSIZE(V()) * (double) XSIZE(V()));
    double gs3 = gs2 * sqrt(2. * PI) * gaussian_Real_sigma * ((double) XSIZE(V()));

    mask_Real3D.get_cont_mask() *= gs3;
    mask_Real.get_cont_mask() *= gs2;

    //For debugging
    /*Image<double> mask2D_test;
    mask2D_test.read("ident2D.spi"); //Identity image
    mask2D_test().setXmippOrigin();
    mask_Real.apply_mask(mask2D_test(), mask2D_test());
    mask2D_test.write("mask_2Dtest.xmp");*/

    // Prepare the weight function in Fourier space
    mask_Fourier.type = GAUSSIAN_MASK;
    mask_Fourier.mode = INNER_MASK;

    mask_Fourier.sigma = gaussian_DFT_sigma * ((double)XSIZE(V()));
    mask_Fourier.generate_mask(YSIZE(V()), XSIZE(V()));

    double gsf2 = 2. * PI * gaussian_DFT_sigma * gaussian_DFT_sigma * ((double)XSIZE(V()) * (double)XSIZE(V()));
    mask_Fourier.get_cont_mask() *= gsf2;
    mask_Fourier.get_cont_mask()(0, 0) *= weight_zero_freq;

    // Weight the input volume in real space
    mask_Real3D.apply_mask(V(), V());

    // Perform the DFT of the reference volume
    int Status;
    reDFTVolume = V();
    imDFTVolume.resize(V());
    CenterFFT(reDFTVolume, false);
    VolumeDftRealToRealImaginary(MULTIDIM_ARRAY(reDFTVolume),
                                 MULTIDIM_ARRAY(imDFTVolume), XSIZE(V()), YSIZE(V()), ZSIZE(V()),
                                 &Status);
    CenterFFT(reDFTVolume, true);
    CenterFFT(imDFTVolume, true);
}

// Predict =================================================================
void ProgAngularContinuousAssign::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    // Read the image and take its angles from the Metadata
    // if they are available. If not, take them from the header.
    // If not, set them to 0.
    Image<double> img;
    img.read(fnImg);
    img().setXmippOrigin();

    double old_rot, old_tilt, old_psi, old_shiftX, old_shiftY;
    rowIn.getValue(MDL_ANGLE_ROT,old_rot);
    rowIn.getValue(MDL_ANGLE_TILT,old_tilt);
    rowIn.getValue(MDL_ANGLE_PSI,old_psi);
    rowIn.getValue(MDL_SHIFT_X,old_shiftX);
    rowIn.getValue(MDL_SHIFT_Y,old_shiftY);

    Matrix1D<double> pose(5);
    pose(0) = old_rot;
    pose(1) = old_tilt;
    pose(2) = old_psi;
    pose(3) = -old_shiftX; // The convention of shifts is different
    pose(4) = -old_shiftY; // for Slavica

    mask_Real.apply_mask(img(), img());

    double cost = CSTSplineAssignment(reDFTVolume, imDFTVolume,
                                      img(), mask_Fourier.get_cont_mask(), pose, max_no_iter);

    Matrix2D<double> Eold, Enew;
    Euler_angles2matrix(old_rot,old_tilt,old_psi,Eold);
    Euler_angles2matrix(pose(0),pose(1),pose(2),Enew);
    double angular_change=Euler_distanceBetweenMatrices(Eold,Enew);
    double shift=sqrt(pose(3)*pose(3)+pose(4)*pose(4));
    double new_rot=old_rot;
    double new_tilt=old_tilt;
    double new_psi=old_psi;
    double new_shiftX=old_shiftX;
    double new_shiftY=old_shiftY;
    if (angular_change<max_angular_change || max_angular_change<0)
    {
    	new_rot    =  pose(0);
    	new_tilt   =  pose(1);
    	new_psi    =  pose(2);
    }
    else
        cost=-1;
    if (shift<max_shift || max_shift<0)
    {
    	new_shiftX = -pose(3);
    	new_shiftY = -pose(4);
    }
    else
        cost=-1;

    rowOut.setValue(MDL_ANGLE_ROT,  new_rot);
    rowOut.setValue(MDL_ANGLE_TILT, new_tilt);
    rowOut.setValue(MDL_ANGLE_PSI,  new_psi);
    rowOut.setValue(MDL_SHIFT_X,    new_shiftX);
    rowOut.setValue(MDL_SHIFT_Y,    new_shiftY);
    rowOut.setValue(MDL_COST,      cost);
}

/* ------------------------------------------------------------------------- */
/* Interface to Slavica's routines                                           */
/* ------------------------------------------------------------------------- */
double CSTSplineAssignment(
    MultidimArray<double> &ReDFTVolume,
    MultidimArray<double> &ImDFTVolume,
    MultidimArray<double> &image,
    MultidimArray<double> &weight,
    Matrix1D<double> &pose_parameters,
    int               max_no_iter
)
{
    // Build the parameter structure .........................................
    cstregistrationStruct Data;

    // Set the Volume input
    Data.ReDftVolume    = MULTIDIM_ARRAY(ReDFTVolume);
    Data.nx_ReDftVolume = XSIZE(ReDFTVolume);
    Data.ny_ReDftVolume = YSIZE(ReDFTVolume);
    Data.nz_ReDftVolume = ZSIZE(ReDFTVolume);

    Data.ImDftVolume    = MULTIDIM_ARRAY(ImDFTVolume);
    Data.nx_ImDftVolume = XSIZE(ImDFTVolume);
    Data.ny_ImDftVolume = YSIZE(ImDFTVolume);
    Data.nz_ImDftVolume = ZSIZE(ImDFTVolume);

    // Perform the DFT of the input image
    int Status;
    MultidimArray<double> realImg(image), imagImg;
    CenterFFT(realImg, false);
    imagImg.resize(image);
    VolumeDftRealToRealImaginary(MULTIDIM_ARRAY(realImg),
                                 MULTIDIM_ARRAY(imagImg), XSIZE(image), YSIZE(image), 1, &Status);
    CenterFFT(realImg, true);
    CenterFFT(imagImg, true);

    // Set the Image input
    Data.ReDftImage     = MULTIDIM_ARRAY(realImg);
    Data.nx_ReDftImage  = XSIZE(realImg);
    Data.ny_ReDftImage  = YSIZE(realImg);

    Data.ImDftImage     = MULTIDIM_ARRAY(imagImg);
    Data.nx_ImDftImage  = XSIZE(imagImg);
    Data.ny_ImDftImage  = YSIZE(imagImg);

    // Set the weight input
    Data.Weight         = MULTIDIM_ARRAY(weight);
    Data.nx_Weight      = XSIZE(weight);
    Data.ny_Weight      = YSIZE(weight);

    // Set the sampling rates
    Matrix1D<double> sampling_rate(3);
    sampling_rate.initConstant(1);
    Data.VoxelSize      = MATRIX1D_ARRAY(sampling_rate);
    Data.nx_VoxelSize   = 3;
    Data.PixelSize      = MATRIX1D_ARRAY(sampling_rate);
    Data.nx_PixelSize   = 2;

    // Set the initial pose parameters
    Data.Parameters     = MATRIX1D_ARRAY(pose_parameters);
    Data.nx_Parameters  = 5;

    // Set the final pose parameters.
    Matrix2D<double> output_pose(5, max_no_iter + 1);
    Data.OutputParameters = MATRIX2D_ARRAY(output_pose);
    Data.nx_OutputParameters = max_no_iter + 1;
    Data.ny_OutputParameters = 5;

    // Set performance parameters
    Matrix1D<double> Cost, TimePerIter, Failures;
    Cost.initZeros(max_no_iter + 1);
    TimePerIter.initZeros(max_no_iter + 1);
    Failures.initZeros(max_no_iter + 1);
    long             NumberIterPerformed, NumberSuccPerformed,
    NumberFailPerformed;
    Data.Cost             = MATRIX1D_ARRAY(Cost);
    Data.nx_Cost          = max_no_iter + 1;
    Data.TimePerIter      = MATRIX1D_ARRAY(TimePerIter);
    Data.nx_TimePerIter   = max_no_iter + 1;
    Data.NumberIterPerformed = &NumberIterPerformed;
    Data.NumberSuccPerformed = &NumberSuccPerformed;
    Data.NumberFailPerformed = &NumberFailPerformed;
    Data.Failures            = MATRIX1D_ARRAY(Failures);
    Data.nx_Failures         = max_no_iter + 1;

    // Set the parameters for the extracted central slice
    MultidimArray<double> dftProj(2, YSIZE(image), XSIZE(image));
    Data.dftProj          = MULTIDIM_ARRAY(dftProj);
    Data.nx_dftProj       = XSIZE(image);
    Data.ny_dftProj       = YSIZE(image);
    Data.nz_dftProj       = 2;

    // Set the optimizer parameters
    Data.ScaleLambda      = 2;
    Data.LambdaInitial    = 1000;
    Data.MaxNoIter        = max_no_iter;
    Data.MaxNoFailure     = (long)(0.3 * max_no_iter);
    Data.SatisfNoSuccess  = (long)(0.7 * max_no_iter);
    ;
    Data.MakeDesiredProj  = 0;
    Data.ToleranceAngle   = 0.0;
    Data.ToleranceShift   = 0.0;

    // Call Slavica's routine ...............................................
    Status = cstregistration(&Data);

    // Retrieve results .....................................................
    // Skip last iterations if they are failures
    double retval;
    if (Status != ERROR)
    {
        long last_iteration_performed = *(Data.NumberIterPerformed) + 1;
        while (Failures(last_iteration_performed - 1) > 0.0)
            last_iteration_performed--;

        // Get the pose parameters
        pose_parameters(0) = RAD2DEG(output_pose(0, last_iteration_performed - 1));
        pose_parameters(1) = RAD2DEG(output_pose(1, last_iteration_performed - 1));
        pose_parameters(2) = RAD2DEG(output_pose(2, last_iteration_performed - 1));
        pose_parameters(3) = output_pose(3, last_iteration_performed - 1);
        pose_parameters(4) = output_pose(4, last_iteration_performed - 1);

        // Get the cost
        retval = Cost(last_iteration_performed - 1);
    }
    else
    {
        retval = -1;
        std::cout << "There is a problem with one image, angles not assigned\n";
    }

    // Return
    return retval;
}

/* ------------------------------------------------------------------------- */
/* Boundary conditions                                                       */
/* ------------------------------------------------------------------------- */
long mirroring(long period, long k)
{
    long sqk, qk;

    if (period == 1L)
        qk = 0L;
    else
        qk = k - (2L * period - 2L) * (long) floor((double) k / (double)(2L * period - 2L));

    if ((qk < period) && (qk >= 0L))
        sqk = qk;
    else
        sqk = 2L * period - 2L - qk;

    return(sqk);
}

long periodization(long period, long k)
{
    long sqk, qk;

    if (period == 1L)
        qk = 0L;
    else if (k >= 0L)
        qk = k - period * (long) floor((double) k / (double) period);
    else
        qk = ABS(k + 1L) - period * (long) floor((double) ABS(k + 1L) / (double) period);

    if ((k >= 0L) || (period == 1L))
        sqk = qk;
    else
        sqk = period - 1L - qk;

    return(sqk);
}

#define PERIODIZATION(y, period, k) \
{ \
	long K=k; \
    long qk; \
    \
    if (period == 1L) \
        qk = 0L; \
    else if (K >= 0L) \
        qk = K - period * (long) floor((double) K / (double) period); \
    else \
    { \
    	long aux=ABS(K+1L); \
        qk = aux - period * (long) floor((double) aux / (double) period); \
    } \
    \
    if ((K >= 0L) || (period == 1L)) \
        y = qk; \
    else \
        y = period - 1L - qk; \
}

/* ------------------------------------------------------------------------- */
/* Gradient and Hessian at pixel                                             */
/* ------------------------------------------------------------------------- */
int gradhesscost_atpixel(
    double *Gradient,
    double *Hessian,
    double *cost,
    double Difference,
    double dp0,
    double dp1,
    double dp2,
    double dp3,
    double dp4,
    double Weight)
{
    long   m, l, CC;

    *cost += Weight * Difference * Difference;

    for (m = 0L; m < 5L; m++)
    {
        CC = m * 5L;
        if (m == 0L)
        {
            Gradient[0] += Weight * dp0 * Difference;
            for (l = 0L; l < 5L; l++)
            {
                if (l == 0L)
                {
                    Hessian[CC + l] += Weight * dp0 * dp0;
                }
                else if (l == 1L)
                {
                    Hessian[CC + l] += Weight * dp0 * dp1;
                }
                else if (l == 2L)
                {
                    Hessian[CC + l] += Weight * dp0 * dp2;
                }
                else if (l == 3L)
                {
                    Hessian[CC + l] += Weight * dp0 * dp3;
                }
                else
                {
                    Hessian[CC + l] += Weight * dp0 * dp4;
                }

            }
        }
        if (m == 1L)
        {
            Gradient[1] += Weight * dp1 * Difference;
            for (l = 1L; l < 5L; l++)
            {
                if (l == 1L)
                {
                    Hessian[CC + l] += Weight * dp1 * dp1;
                }
                else if (l == 2L)
                {
                    Hessian[CC + l] += Weight * dp1 * dp2;
                }
                else if (l == 3L)
                {
                    Hessian[CC + l] += Weight * dp1 * dp3;
                }
                else
                {
                    Hessian[CC + l] += Weight * dp1 * dp4;
                }

            }
        }
        if (m == 2L)
        {
            Gradient[2] += Weight * dp2 * Difference;
            for (l = 2L; l < 5L; l++)
            {
                if (l == 2L)
                {
                    Hessian[CC + l] += Weight * dp2 * dp2;
                }
                else if (l == 3L)
                {
                    Hessian[CC + l] += Weight * dp2 * dp3;
                }
                else
                {
                    Hessian[CC + l] += Weight * dp2 * dp4;
                }

            }
        }
        if (m == 3L)
        {
            Gradient[3] += Weight * dp3 * Difference;
            for (l = 3L; l < 5L; l++)
            {
                if (l == 3L)
                {
                    Hessian[CC + l] += Weight * dp3 * dp3;
                }
                else
                {
                    Hessian[CC + l] += Weight * dp3 * dp4;
                }

            }
        }
        if (m == 4L)
        {
            Gradient[4] += Weight * dp4 * Difference;
            Hessian[CC + 4L] += Weight * dp4 * dp4;
        }
    }

    return(!ERROR);
}/* End of gradhesscost_atpixel */

/* ------------------------------------------------------------------------- */
/* Optimizer                                                                 */
/* ------------------------------------------------------------------------- */
int return_gradhesscost(
              double *Gradient,
              double *Hessian,
              double *cost,
              double *Parameters,
              double *CoefRe,
              double *CoefIm,
              double *pntr_ReInp,
              double *pntr_ImInp,
              double *pntr_Weight,
              long   Nx,
              long   Ny,
              long   Nz,
              long   indx0,
              long   indy0,
              long   indz0,
              long   Mx,
              long   My,
              double *Left,
              double *Right,
              double *Q1,
              double *Q3,
              double scx1,
              double scy1,
              double S,
              long   SizeIm,
              int     DoDesProj,
              double *dftCompProj)
          {
              int      Status = !ERROR;
              long     indx, indy, l, l1, l2, m, m1, m2, n, n1, n2;
              long     i, j, CC1, CC, index;
              double   phi, theta, psi, x0, y0, Sinphi, Cosphi, Sinpsi, Cospsi, Sintheta, Costheta;
              double   *hlp, *Rz1, *Ry, *Rz2, *DRz1, *DRy, *DRz2;
              double   *DR0, *DR1, *DR2, *mn, *Arg;
              double   *Grad_re, *Grad_im, *Hessian_re, *Hessian_im, *Gradient_re, *Gradient_im;
              double   scx, scy, x, y, z;
              double   xminusl, yminusm, zminusn, re_Coeff, im_Coeff;
              double   re_columns, re_rows, im_columns, im_rows;
              double   re_slices, re_slices_d1, re_slices_d2, re_slices_d3, re_columns_d1, re_columns_d2, re_rows_d1;
              double   im_slices, im_slices_d1, im_slices_d2, im_slices_d3, im_columns_d1, im_columns_d2, im_rows_d1;
              double   a, b, c, SinC, CosC, redftproj, imdftproj;
              double   da, db, dc, da0[1], da1[1], da2[1], da3[1], da4[1];
              double   db0[1], db1[1], db2[1], db3[1], db4[1], dc0[1], dc1[1], dc2[1], dc3[1], dc4[1];
              double   prv_re, prv_im, Difference_re, Difference_im, Weight;
              double   dp0_re, dp0_im, dp1_re, dp1_im, dp2_re, dp2_im, dp3_re, dp3_im, dp4_re, dp4_im;
              double   cost_re, cost_im, *R, *Q, *Q2, *dQ0, *dQ1, *dQ2, *dQ2_0, *dQ2_1, *dQ2_2;
              double   sum_xx_re, sum_xx_im, sc_re, sc_im;
              double  *DP_0, *DP_1, *DP_2, *DP_3, *DP_4, *pntr_ReOut, *pntr_ImOut;
              double  *pntr_DP_re, *pntr_DP_im, *pntr_DP_0_re, *pntr_DP_0_im, *pntr_DP_1_re, *pntr_DP_1_im;
              double  *pntr_DP_2_re, *pntr_DP_2_im, *pntr_DP_3_re, *pntr_DP_3_im;
              double  *pntr_DP_4_re, *pntr_DP_4_im;

              phi = Parameters[0];
              theta = Parameters[1];
              psi = Parameters[2];
              x0 = Parameters[3];
              y0 = Parameters[4];

              Sinphi = sin(phi);
              Cosphi = cos(phi);
              Sinpsi = sin(psi);
              Cospsi = cos(psi);
              Sintheta = sin(theta);
              Costheta = cos(theta);

              Rz1 = (double *)malloc((size_t) 16L * sizeof(double));
              if (Rz1 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Rz1");
                  return(ERROR);
              }

              Ry = (double *)malloc((size_t) 16L * sizeof(double));
              if (Ry == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Ry");
                  free(Rz1);
                  return(ERROR);
              }

              Rz2 = (double *)malloc((size_t) 16L * sizeof(double));
              if (Rz2 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Rz2");
                  free(Rz1);
                  free(Ry);
                  return(ERROR);
              }

              if (GetIdentitySquareMatrix(Rz2, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by GetIdentitySquareMatrix");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  return(ERROR);
              }

              hlp = Rz2;
              *hlp++ = Cospsi;
              *hlp = Sinpsi;
              hlp += (std::ptrdiff_t)3L;
              *hlp++ = - Sinpsi;
              *hlp = Cospsi;

              if (GetIdentitySquareMatrix(Rz1, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by GetIdentitySquareMatrix");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  return(ERROR);
              }

              hlp = Rz1;
              *hlp++ = Cosphi;
              *hlp = Sinphi;
              hlp += (std::ptrdiff_t)3L;
              *hlp++ = - Sinphi;
              *hlp = Cosphi;

              if (GetIdentitySquareMatrix(Ry, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by GetIdentitySquareMatrix");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  return(ERROR);
              }

              hlp = Ry;
              *hlp = Costheta;
              hlp += (std::ptrdiff_t)2L;
              *hlp = - Sintheta;
              hlp += (std::ptrdiff_t)6L;
              *hlp = Sintheta;
              hlp += (std::ptrdiff_t)2L;
              *hlp = Costheta;

              R = (double *)malloc((size_t) 16L * sizeof(double));
              if (R == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for R");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  return(ERROR);
              }

              if (multiply_3Matrices(Rz2, Ry, Rz1, R, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  return(ERROR);
              }

              DRz1 = (double *)malloc((size_t) 16L * sizeof(double));
              if (DRz1 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DRz1");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  return(ERROR);
              }

              DRy = (double *)malloc((size_t) 16L * sizeof(double));
              if (DRy == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DRy");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  return(ERROR);
              }

              DRz2 = (double *)malloc((size_t) 16L * sizeof(double));
              if (DRz2 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DRz2");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  return(ERROR);
              }

              for (i = 0L; i < 16L; i++)
              {
                  DRz2[i] = 0.0;
                  DRz1[i] = 0.0;
                  DRy[i] = 0.0;
              }

              hlp = DRz2;
              *hlp++ = - Sinpsi;
              *hlp = Cospsi;
              hlp += (std::ptrdiff_t)3L;
              *hlp++ = - Cospsi;
              *hlp = - Sinpsi;

              hlp = DRz1;
              *hlp++ = - Sinphi;
              *hlp = Cosphi;
              hlp += (std::ptrdiff_t)3L;
              *hlp++ = - Cosphi;
              *hlp = - Sinphi;

              hlp = DRy;
              *hlp = - Sintheta;
              hlp += (std::ptrdiff_t)2L;
              *hlp = - Costheta;
              hlp += (std::ptrdiff_t)6L;
              *hlp = Costheta;
              hlp += (std::ptrdiff_t)2L;
              *hlp = - Sintheta;

              DR0 = (double *)malloc((size_t) 16L * sizeof(double));
              if (DR0 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DR0");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  return(ERROR);
              }
              DR1 = (double *)malloc((size_t) 16L * sizeof(double));
              if (DR1 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DR1");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  return(ERROR);
              }
              DR2 = (double *)malloc((size_t) 16L * sizeof(double));
              if (DR2 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DR2");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  return(ERROR);
              }

              if (multiply_3Matrices(Rz2, Ry, DRz1, DR0, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  return(ERROR);
              }
              if (multiply_3Matrices(Rz2, DRy, Rz1, DR1, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  return(ERROR);
              }
              if (multiply_3Matrices(DRz2, Ry, Rz1, DR2, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  return(ERROR);
              }

              mn = (double *)malloc((size_t) 4L * sizeof(double));
              if (mn == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for mn");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  return(ERROR);
              }

              Arg = (double *)malloc((size_t) 4L * sizeof(double));
              if (Arg == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Arg");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  return(ERROR);
              }

              Grad_re = (double *)malloc((size_t) 4L * sizeof(double));
              if (Grad_re == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Grad_re");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  return(ERROR);
              }
              Grad_im = (double *)malloc((size_t) 4L * sizeof(double));
              if (Grad_im == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Grad_im");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  return(ERROR);
              }

              Gradient_re = (double *)malloc((size_t) 5L * sizeof(double));
              if (Gradient_re == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Gradient_re");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  return(ERROR);

              }

              Gradient_im = (double *)malloc((size_t) 5L * sizeof(double));
              if (Gradient_im == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Gradient_im");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  return(ERROR);

              }

              Hessian_re = (double *)malloc((size_t) 25L * sizeof(double));
              if (Hessian_re == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Hessian_re");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  return(ERROR);

              }

              Hessian_im = (double *)malloc((size_t) 25L * sizeof(double));
              if (Hessian_im == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Hessian_im");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  return(ERROR);

              }

              AllocateVolumeDouble(&DP_0, Mx, My, 2L, &Status);
              if (Status == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DP_0");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  return(ERROR);
              }
              AllocateVolumeDouble(&DP_1, Mx, My, 2L, &Status);
              if (Status == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DP_1");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  return(ERROR);
              }
              AllocateVolumeDouble(&DP_2, Mx, My, 2L, &Status);
              if (Status == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DP_2");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  return(ERROR);
              }
              AllocateVolumeDouble(&DP_3, Mx, My, 2L, &Status);
              if (Status == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DP_3");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  return(ERROR);
              }
              AllocateVolumeDouble(&DP_4, Mx, My, 2L, &Status);
              if (Status == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for DP_4");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  return(ERROR);
              }

              pntr_ReOut = dftCompProj;
              pntr_ImOut = pntr_ReOut + (std::ptrdiff_t) SizeIm;

              pntr_DP_0_re = DP_0;
              pntr_DP_0_im = pntr_DP_0_re + (std::ptrdiff_t) SizeIm;

              pntr_DP_1_re = DP_1;
              pntr_DP_1_im = pntr_DP_1_re + (std::ptrdiff_t) SizeIm;

              pntr_DP_2_re = DP_2;
              pntr_DP_2_im = pntr_DP_2_re + (std::ptrdiff_t) SizeIm;

              pntr_DP_3_re = DP_3;
              pntr_DP_3_im = pntr_DP_3_re + (std::ptrdiff_t) SizeIm;

              pntr_DP_4_re = DP_4;
              pntr_DP_4_im = pntr_DP_4_re + (std::ptrdiff_t) SizeIm;



              Q = (double *)malloc((size_t) 16L * sizeof(double));
              if (Q == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Q");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  return(ERROR);
              }

              Q2 = (double *)malloc((size_t) 16L * sizeof(double));
              if (Q2 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for Q2");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  return(ERROR);
              }

              if (multiply_3Matrices(Left, R, Right, Q, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  return(ERROR);
              }
              if (MatrixTranspose(Q, Q2, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by MatrixTranspose");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  return(ERROR);
              }

              if (multiply_3Matrices(Q1, Q2, Q3, Q, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  return(ERROR);
              }


              dQ0 = (double *)malloc((size_t) 16L * sizeof(double));
              if (dQ0 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for dQ0");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  return(ERROR);
              }
              dQ1 = (double *)malloc((size_t) 16L * sizeof(double));
              if (dQ1 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for dQ1");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  free(dQ0);
                  return(ERROR);
              }
              dQ2 = (double *)malloc((size_t) 16L * sizeof(double));
              if (dQ2 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for dQ2");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  free(dQ0);
                  free(dQ1);
                  return(ERROR);
              }


              dQ2_0 = (double *)malloc((size_t) 16L * sizeof(double));
              if (dQ2_0 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for dQ2_0");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  free(dQ0);
                  free(dQ1);
                  free(dQ2);
                  return(ERROR);
              }
              dQ2_1 = (double *)malloc((size_t) 16L * sizeof(double));
              if (dQ2_1 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for dQ2_1");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  free(dQ0);
                  free(dQ1);
                  free(dQ2);
                  free(dQ2_0);
                  return(ERROR);
              }
              dQ2_2 = (double *)malloc((size_t) 16L * sizeof(double));
              if (dQ2_2 == (double *)NULL)
              {
                  WRITE_ERROR(return_gradhesscost, "ERROR - Not enough memory for dQ2_2");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  free(dQ0);
                  free(dQ1);
                  free(dQ2);
                  free(dQ2_0);
                  free(dQ2_1);
                  return(ERROR);
              }

              if (multiply_3Matrices(Left, DR0, Right, dQ0, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  free(dQ0);
                  free(dQ1);
                  free(dQ2);
                  free(dQ2_0);
                  free(dQ2_1);
                  free(dQ2_2);
                  return(ERROR);
              }
              if (MatrixTranspose(dQ0, dQ2_0, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by MatrixTranspose");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  free(dQ0);
                  free(dQ1);
                  free(dQ2);
                  free(dQ2_0);
                  free(dQ2_1);
                  free(dQ2_2);
                  return(ERROR);
              }
              if (multiply_3Matrices(Left, DR1, Right, dQ1, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  free(dQ0);
                  free(dQ1);
                  free(dQ2);
                  free(dQ2_0);
                  free(dQ2_1);
                  free(dQ2_2);
                  return(ERROR);
              }
              if (MatrixTranspose(dQ1, dQ2_1, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by MatrixTranspose");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  free(dQ0);
                  free(dQ1);
                  free(dQ2);
                  free(dQ2_0);
                  free(dQ2_1);
                  free(dQ2_2);
                  return(ERROR);
              }
              if (multiply_3Matrices(Left, DR2, Right, dQ2, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  free(dQ0);
                  free(dQ1);
                  free(dQ2);
                  free(dQ2_0);
                  free(dQ2_1);
                  free(dQ2_2);
                  return(ERROR);
              }
              if (MatrixTranspose(dQ2, dQ2_2, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by MatrixTranspose");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  free(dQ0);
                  free(dQ1);
                  free(dQ2);
                  free(dQ2_0);
                  free(dQ2_1);
                  free(dQ2_2);
                  return(ERROR);
              }


              if (multiply_3Matrices(Q1, dQ2_0, Q3, dQ0, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  free(dQ0);
                  free(dQ1);
                  free(dQ2);
                  free(dQ2_0);
                  free(dQ2_1);
                  free(dQ2_2);
                  return(ERROR);
              }

              if (multiply_3Matrices(Q1, dQ2_1, Q3, dQ1, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  free(dQ0);
                  free(dQ1);
                  free(dQ2);
                  free(dQ2_0);
                  free(dQ2_1);
                  free(dQ2_2);
                  return(ERROR);
              }

              if (multiply_3Matrices(Q1, dQ2_2, Q3, dQ2, 4L, 4L, 4L, 4L) == ERROR)
              {
                  WRITE_ERROR(return_gradhesscost, "Error returned by multiply_3Matrices");
                  free(Rz1);
                  free(Ry);
                  free(Rz2);
                  free(R);
                  free(DRz1);
                  free(DRy);
                  free(DRz2);
                  free(DR0);
                  free(DR1);
                  free(DR2);
                  free(mn);
                  free(Arg);
                  free(Grad_re);
                  free(Grad_im);
                  free(Gradient_re);
                  free(Gradient_im);
                  free(Hessian_re);
                  free(Hessian_im);
                  FreeVolumeDouble(&DP_0);
                  FreeVolumeDouble(&DP_1);
                  FreeVolumeDouble(&DP_2);
                  FreeVolumeDouble(&DP_3);
                  FreeVolumeDouble(&DP_4);
                  free(Q);
                  free(Q2);
                  free(dQ0);
                  free(dQ1);
                  free(dQ2);
                  free(dQ2_0);
                  free(dQ2_1);
                  free(dQ2_2);
                  return(ERROR);
              }


              sum_xx_re = 0.0;
              sum_xx_im = 0.0;

              scx = scx1 * x0;
              scy = scy1 * y0;

              double aux, aux2;
              for (indy = 0L; indy < My; indy++)
              {
                  for (indx = 0L; indx < Mx; indx++)
                  {

                      mn[0] = (double)(indx - Mx / 2L);
                      mn[1] = (double)(indy - My / 2L);
                      mn[2] = 0.0;
                      mn[3] = 1.0;

                      if (MatrixTimesVector(Q, mn, Arg, 4L, 4L) == ERROR)
                      {
                          WRITE_ERROR(return_gradhesscost, "Error returned by MatrixTimesVector");
                          free(Rz1);
                          free(Ry);
                          free(Rz2);
                          free(R);
                          free(DRz1);
                          free(DRy);
                          free(DRz2);
                          free(DR0);
                          free(DR1);
                          free(DR2);
                          free(mn);
                          free(Arg);
                          free(Grad_re);
                          free(Grad_im);
                          free(Gradient_re);
                          free(Gradient_im);
                          free(Hessian_re);
                          free(Hessian_im);
                          FreeVolumeDouble(&DP_0);
                          FreeVolumeDouble(&DP_1);
                          FreeVolumeDouble(&DP_2);
                          FreeVolumeDouble(&DP_3);
                          FreeVolumeDouble(&DP_4);
                          free(Q);
                          free(Q2);
                          free(dQ0);
                          free(dQ1);
                          free(dQ2);
                          free(dQ2_0);
                          free(dQ2_1);
                          free(dQ2_2);
                          return(ERROR);
                      }

                      x = Arg[0];
                      y = Arg[1];
                      z = Arg[2];

                      l1 = (long) ceil(x - 2.0);
                      l2 = l1 + 3L;
                      m1 = (long) ceil(y - 2.0);
                      m2 = m1 + 3L;
                      n1 = (long) ceil(z - 2.0);
                      n2 = n1 + 3L;

                      re_slices = 0.0;
                      re_slices_d1 = 0.0;
                      re_slices_d2 = 0.0;
                      re_slices_d3 = 0.0;
                      im_slices = 0.0;
                      im_slices_d1 = 0.0;
                      im_slices_d2 = 0.0;
                      im_slices_d3 = 0.0;
                      for (n = n1; n <= n2; n++)
                      {
                    	  long Laux;
                    	  PERIODIZATION(Laux,Nz, n - indz0);
                          CC1 = Nx * Ny * Laux;
                          re_columns = 0.0;
                          re_columns_d1 = 0.0;
                          re_columns_d2 = 0.0;
                          im_columns = 0.0;
                          im_columns_d1 = 0.0;
                          im_columns_d2 = 0.0;
                          for (m = m1; m <= m2; m++)
                          {
                        	  PERIODIZATION(Laux,Ny, m - indy0);
                              CC = CC1 + Nx * Laux;
                              re_rows = 0.0;
                              re_rows_d1 = 0.0;
                              im_rows = 0.0;
                              im_rows_d1 = 0.0;
                              for (l = l1; l <= l2; l++)
                              {
                                  xminusl = x - (double) l;
                            	  PERIODIZATION(Laux,Nx, l - indx0);
                            	  long Laux2=CC + Laux;
                                  re_Coeff = CoefRe[Laux2];
                                  im_Coeff = CoefIm[Laux2];
                                  BSPLINE03(aux,xminusl);
                                  BSPLINE03DIFF1(aux2,xminusl);
                                  re_rows += re_Coeff * aux;
                                  re_rows_d1 += re_Coeff * aux2;
                                  im_rows += im_Coeff * aux;
                                  im_rows_d1 += im_Coeff * aux2;
                              }
                              yminusm = y - (double) m;
                              BSPLINE03(aux,yminusm);
                              BSPLINE03DIFF1(aux2,yminusm);
                              re_columns +=  re_rows * aux;
                              re_columns_d1 +=  re_rows_d1 * aux;
                              re_columns_d2 +=  re_rows * aux2;
                              im_columns +=  im_rows * aux;
                              im_columns_d1 +=  im_rows_d1 * aux;
                              im_columns_d2 +=  im_rows * aux2;
                          }
                          zminusn = z - (double) n;
                          BSPLINE03(aux,zminusn);
                          BSPLINE03DIFF1(aux2,zminusn);
                          re_slices +=  re_columns * aux;
                          re_slices_d1 +=  re_columns_d1 * aux;
                          re_slices_d2 +=  re_columns_d2 * aux;
                          re_slices_d3 +=  re_columns * aux2;
                          im_slices +=  im_columns * aux;
                          im_slices_d1 +=  im_columns_d1 * aux;
                          im_slices_d2 +=  im_columns_d2 * aux;
                          im_slices_d3 +=  im_columns * aux2;
                      }

                      a = re_slices;
                      Grad_re[0] = re_slices_d1;
                      Grad_re[1] = re_slices_d2;
                      Grad_re[2] = re_slices_d3;
                      Grad_re[3] = 0.0;

                      b = im_slices;
                      Grad_im[0] = im_slices_d1;
                      Grad_im[1] = im_slices_d2;
                      Grad_im[2] = im_slices_d3;
                      Grad_im[3] = 0.0;

                      c = scx * mn[0] + scy * mn[1];

                      SinC = sin(c);
                      CosC = cos(c);

                      redftproj = (a * CosC + b * SinC) * S;
                      imdftproj = (b * CosC - a * SinC) * S;

                      index = indy * Mx + indx;

                      pntr_ReOut[index] = redftproj;
                      pntr_ImOut[index] = imdftproj;


                      if (!((indx == (Mx / 2L)) && (indy == (My / 2L))))
                      {
                          sum_xx_re += redftproj * redftproj;
                          sum_xx_im += imdftproj * imdftproj;
                      }

                      if (MatrixTimesVector(dQ0, mn, Arg, 4L, 4L) == ERROR)
                      {
                          WRITE_ERROR(return_gradhesscost, "Error returned by MatrixTimesVector");
                          free(Rz1);
                          free(Ry);
                          free(Rz2);
                          free(R);
                          free(DRz1);
                          free(DRy);
                          free(DRz2);
                          free(DR0);
                          free(DR1);
                          free(DR2);
                          free(mn);
                          free(Arg);
                          free(Grad_re);
                          free(Grad_im);
                          free(Gradient_re);
                          free(Gradient_im);
                          free(Hessian_re);
                          free(Hessian_im);
                          FreeVolumeDouble(&DP_0);
                          FreeVolumeDouble(&DP_1);
                          FreeVolumeDouble(&DP_2);
                          FreeVolumeDouble(&DP_3);
                          FreeVolumeDouble(&DP_4);
                          free(Q);
                          free(Q2);
                          free(dQ0);
                          free(dQ1);
                          free(dQ2);
                          free(dQ2_0);
                          free(dQ2_1);
                          free(dQ2_2);
                          return(ERROR);
                      }
                      if (VectorScalarProduct(Grad_re, Arg, da0, 4L) == ERROR)
                      {
                          WRITE_ERROR(return_gradhesscost, "Error returned by VectorScalarProduct");
                          free(Rz1);
                          free(Ry);
                          free(Rz2);
                          free(R);
                          free(DRz1);
                          free(DRy);
                          free(DRz2);
                          free(DR0);
                          free(DR1);
                          free(DR2);
                          free(mn);
                          free(Arg);
                          free(Grad_re);
                          free(Grad_im);
                          free(Gradient_re);
                          free(Gradient_im);
                          free(Hessian_re);
                          free(Hessian_im);
                          FreeVolumeDouble(&DP_0);
                          FreeVolumeDouble(&DP_1);
                          FreeVolumeDouble(&DP_2);
                          FreeVolumeDouble(&DP_3);
                          FreeVolumeDouble(&DP_4);
                          free(Q);
                          free(Q2);
                          free(dQ0);
                          free(dQ1);
                          free(dQ2);
                          free(dQ2_0);
                          free(dQ2_1);
                          free(dQ2_2);
                          return(ERROR);
                      }
                      if (VectorScalarProduct(Grad_im, Arg, db0, 4L) == ERROR)
                      {
                          WRITE_ERROR(return_gradhesscost, "Error returned by VectorScalarProduct");
                          free(Rz1);
                          free(Ry);
                          free(Rz2);
                          free(R);
                          free(DRz1);
                          free(DRy);
                          free(DRz2);
                          free(DR0);
                          free(DR1);
                          free(DR2);
                          free(mn);
                          free(Arg);
                          free(Grad_re);
                          free(Grad_im);
                          free(Gradient_re);
                          free(Gradient_im);
                          free(Hessian_re);
                          free(Hessian_im);
                          FreeVolumeDouble(&DP_0);
                          FreeVolumeDouble(&DP_1);
                          FreeVolumeDouble(&DP_2);
                          FreeVolumeDouble(&DP_3);
                          FreeVolumeDouble(&DP_4);
                          free(Q);
                          free(Q2);
                          free(dQ0);
                          free(dQ1);
                          free(dQ2);
                          free(dQ2_0);
                          free(dQ2_1);
                          free(dQ2_2);
                          return(ERROR);
                      }

                      if (MatrixTimesVector(dQ1, mn, Arg, 4L, 4L) == ERROR)
                      {
                          WRITE_ERROR(return_gradhesscost, "Error returned by MatrixTimesVector");
                          free(Rz1);
                          free(Ry);
                          free(Rz2);
                          free(R);
                          free(DRz1);
                          free(DRy);
                          free(DRz2);
                          free(DR0);
                          free(DR1);
                          free(DR2);
                          free(mn);
                          free(Arg);
                          free(Grad_re);
                          free(Grad_im);
                          free(Gradient_re);
                          free(Gradient_im);
                          free(Hessian_re);
                          free(Hessian_im);
                          FreeVolumeDouble(&DP_0);
                          FreeVolumeDouble(&DP_1);
                          FreeVolumeDouble(&DP_2);
                          FreeVolumeDouble(&DP_3);
                          FreeVolumeDouble(&DP_4);
                          free(Q);
                          free(Q2);
                          free(dQ0);
                          free(dQ1);
                          free(dQ2);
                          free(dQ2_0);
                          free(dQ2_1);
                          free(dQ2_2);
                          return(ERROR);
                      }
                      if (VectorScalarProduct(Grad_re, Arg, da1, 4L) == ERROR)
                      {
                          WRITE_ERROR(return_gradhesscost, "Error returned by VectorScalarProduct");
                          free(Rz1);
                          free(Ry);
                          free(Rz2);
                          free(R);
                          free(DRz1);
                          free(DRy);
                          free(DRz2);
                          free(DR0);
                          free(DR1);
                          free(DR2);
                          free(mn);
                          free(Arg);
                          free(Grad_re);
                          free(Grad_im);
                          free(Gradient_re);
                          free(Gradient_im);
                          free(Hessian_re);
                          free(Hessian_im);
                          FreeVolumeDouble(&DP_0);
                          FreeVolumeDouble(&DP_1);
                          FreeVolumeDouble(&DP_2);
                          FreeVolumeDouble(&DP_3);
                          FreeVolumeDouble(&DP_4);
                          free(Q);
                          free(Q2);
                          free(dQ0);
                          free(dQ1);
                          free(dQ2);
                          free(dQ2_0);
                          free(dQ2_1);
                          free(dQ2_2);
                          return(ERROR);
                      }
                      if (VectorScalarProduct(Grad_im, Arg, db1, 4L) == ERROR)
                      {
                          WRITE_ERROR(return_gradhesscost, "Error returned by VectorScalarProduct");
                          free(Rz1);
                          free(Ry);
                          free(Rz2);
                          free(R);
                          free(DRz1);
                          free(DRy);
                          free(DRz2);
                          free(DR0);
                          free(DR1);
                          free(DR2);
                          free(mn);
                          free(Arg);
                          free(Grad_re);
                          free(Grad_im);
                          free(Gradient_re);
                          free(Gradient_im);
                          free(Hessian_re);
                          free(Hessian_im);
                          FreeVolumeDouble(&DP_0);
                          FreeVolumeDouble(&DP_1);
                          FreeVolumeDouble(&DP_2);
                          FreeVolumeDouble(&DP_3);
                          FreeVolumeDouble(&DP_4);
                          free(Q);
                          free(Q2);
                          free(dQ0);
                          free(dQ1);
                          free(dQ2);
                          free(dQ2_0);
                          free(dQ2_1);
                          free(dQ2_2);
                          return(ERROR);
                      }

                      if (MatrixTimesVector(dQ2, mn, Arg, 4L, 4L) == ERROR)
                      {
                          WRITE_ERROR(return_gradhesscost, "Error returned by MatrixTimesVector");
                          free(Rz1);
                          free(Ry);
                          free(Rz2);
                          free(R);
                          free(DRz1);
                          free(DRy);
                          free(DRz2);
                          free(DR0);
                          free(DR1);
                          free(DR2);
                          free(mn);
                          free(Arg);
                          free(Grad_re);
                          free(Grad_im);
                          free(Gradient_re);
                          free(Gradient_im);
                          free(Hessian_re);
                          free(Hessian_im);
                          FreeVolumeDouble(&DP_0);
                          FreeVolumeDouble(&DP_1);
                          FreeVolumeDouble(&DP_2);
                          FreeVolumeDouble(&DP_3);
                          FreeVolumeDouble(&DP_4);
                          free(Q);
                          free(Q2);
                          free(dQ0);
                          free(dQ1);
                          free(dQ2);
                          free(dQ2_0);
                          free(dQ2_1);
                          free(dQ2_2);
                          return(ERROR);
                      }
                      if (VectorScalarProduct(Grad_re, Arg, da2, 4L) == ERROR)
                      {
                          WRITE_ERROR(return_gradhesscost, "Error returned by VectorScalarProduct");
                          free(Rz1);
                          free(Ry);
                          free(Rz2);
                          free(R);
                          free(DRz1);
                          free(DRy);
                          free(DRz2);
                          free(DR0);
                          free(DR1);
                          free(DR2);
                          free(mn);
                          free(Arg);
                          free(Grad_re);
                          free(Grad_im);
                          free(Gradient_re);
                          free(Gradient_im);
                          free(Hessian_re);
                          free(Hessian_im);
                          FreeVolumeDouble(&DP_0);
                          FreeVolumeDouble(&DP_1);
                          FreeVolumeDouble(&DP_2);
                          FreeVolumeDouble(&DP_3);
                          FreeVolumeDouble(&DP_4);
                          free(Q);
                          free(Q2);
                          free(dQ0);
                          free(dQ1);
                          free(dQ2);
                          free(dQ2_0);
                          free(dQ2_1);
                          free(dQ2_2);
                          return(ERROR);
                      }
                      if (VectorScalarProduct(Grad_im, Arg, db2, 4L) == ERROR)
                      {
                          WRITE_ERROR(return_gradhesscost, "Error returned by VectorScalarProduct");
                          free(Rz1);
                          free(Ry);
                          free(Rz2);
                          free(R);
                          free(DRz1);
                          free(DRy);
                          free(DRz2);
                          free(DR0);
                          free(DR1);
                          free(DR2);
                          free(mn);
                          free(Arg);
                          free(Grad_re);
                          free(Grad_im);
                          free(Gradient_re);
                          free(Gradient_im);
                          free(Hessian_re);
                          free(Hessian_im);
                          FreeVolumeDouble(&DP_0);
                          FreeVolumeDouble(&DP_1);
                          FreeVolumeDouble(&DP_2);
                          FreeVolumeDouble(&DP_3);
                          FreeVolumeDouble(&DP_4);
                          free(Q);
                          free(Q2);
                          free(dQ0);
                          free(dQ1);
                          free(dQ2);
                          free(dQ2_0);
                          free(dQ2_1);
                          free(dQ2_2);
                          return(ERROR);
                      }

                      da3[0] = 0.0;
                      db3[0] = 0.0;
                      da4[0] = 0.0;
                      db4[0] = 0.0;

                      dc0[0] = 0.0;
                      dc1[0] = 0.0;
                      dc2[0] = 0.0;

                      dc3[0] = scx1 * mn[0];
                      dc4[0] = scy1 * mn[1];

                      da = da0[0];
                      dc = dc0[0];
                      db = db0[0];
                      pntr_DP_re = pntr_DP_0_re;
                      pntr_DP_im = pntr_DP_0_im;
                      pntr_DP_re[index] =
                          (da * CosC - a * SinC * dc + db * SinC + b * CosC * dc) * S;
                      pntr_DP_im[index] =
                          (db * CosC - b * SinC * dc - da * SinC - a * CosC * dc) * S;


                      da = da1[0];
                      dc = dc1[0];
                      db = db1[0];
                      pntr_DP_re = pntr_DP_1_re;
                      pntr_DP_im = pntr_DP_1_im;
                      pntr_DP_re[index] =
                          (da * CosC - a * SinC * dc + db * SinC + b * CosC * dc) * S;
                      pntr_DP_im[index] =
                          ((db * CosC - b * SinC * dc - da * SinC - a * CosC * dc) * S);

                      da = da2[0];
                      dc = dc2[0];
                      db = db2[0];
                      pntr_DP_re = pntr_DP_2_re;
                      pntr_DP_im = pntr_DP_2_im;
                      pntr_DP_re[index] =
                          ((da * CosC - a * SinC * dc + db * SinC + b * CosC * dc) * S);
                      pntr_DP_im[index] =
                          ((db * CosC - b * SinC * dc - da * SinC - a * CosC * dc) * S);


                      da = da3[0];
                      dc = dc3[0];
                      db = db3[0];
                      pntr_DP_re = pntr_DP_3_re;
                      pntr_DP_im = pntr_DP_3_im;
                      pntr_DP_re[index] =
                          ((da * CosC - a * SinC * dc + db * SinC + b * CosC * dc) * S);
                      pntr_DP_im[index] =
                          ((db * CosC - b * SinC * dc - da * SinC - a * CosC * dc) * S);

                      da = da4[0];
                      dc = dc4[0];
                      db = db4[0];
                      pntr_DP_re = pntr_DP_4_re;
                      pntr_DP_im = pntr_DP_4_im;
                      pntr_DP_re[index] =
                          ((da * CosC - a * SinC * dc + db * SinC + b * CosC * dc) * S);
                      pntr_DP_im[index] =
                          ((db * CosC - b * SinC * dc - da * SinC - a * CosC * dc) * S);

                  }
              }

              if (!DoDesProj)
              {

                  cost_re = 0.0;
                  cost_im = 0.0;
                  for (i = 0L; i < 5L; i++)
                  {
                      Gradient_re[i] = 0.0;
                      Gradient_im[i] = 0.0;
                      for (j = 0L; j < 5L; j++)
                      {
                          Hessian_re[i * 5L + j] = 0.0;
                          Hessian_im[i * 5L + j] = 0.0;
                      }
                  }

                  sc_re = 1.0 / sqrt(sum_xx_re + sum_xx_im);
                  sc_im = sc_re;


                  for (indy = 0L; indy < My; indy++)
                  {
                      for (indx = 0L; indx < Mx; indx++)
                      {

                          index = indy * Mx + indx;

                          Weight = (double) pntr_Weight[index];

                          prv_re = (sc_re * (double) pntr_ReOut[index]);
                          Difference_re = (double)(prv_re - pntr_ReInp[index]);

                          prv_im = (sc_im * (double) pntr_ImOut[index]);
                          Difference_im = (double)(prv_im - pntr_ImInp[index]);

                          dp0_re = sc_re * (double) pntr_DP_0_re[index];
                          dp1_re = sc_re * (double) pntr_DP_1_re[index];
                          dp2_re = sc_re * (double) pntr_DP_2_re[index];
                          dp3_re = sc_re * (double) pntr_DP_3_re[index];
                          dp4_re = sc_re * (double) pntr_DP_4_re[index];

                          dp0_im = sc_im * (double) pntr_DP_0_im[index];
                          dp1_im = sc_im * (double) pntr_DP_1_im[index];
                          dp2_im = sc_im * (double) pntr_DP_2_im[index];
                          dp3_im = sc_im * (double) pntr_DP_3_im[index];
                          dp4_im = sc_im * (double) pntr_DP_4_im[index];

                          if (gradhesscost_atpixel(Gradient_re, Hessian_re, &cost_re, Difference_re, dp0_re, dp1_re, dp2_re, dp3_re, dp4_re, Weight) == ERROR)
                          {
                              WRITE_ERROR(return_gradhesscost, "Error returned by gradhesscost_atpixel");
                              free(Rz1);
                              free(Ry);
                              free(Rz2);
                              free(R);
                              free(DRz1);
                              free(DRy);
                              free(DRz2);
                              free(DR0);
                              free(DR1);
                              free(DR2);
                              free(mn);
                              free(Arg);
                              free(Grad_re);
                              free(Grad_im);
                              free(Gradient_re);
                              free(Gradient_im);
                              free(Hessian_re);
                              free(Hessian_im);
                              FreeVolumeDouble(&DP_0);
                              FreeVolumeDouble(&DP_1);
                              FreeVolumeDouble(&DP_2);
                              FreeVolumeDouble(&DP_3);
                              FreeVolumeDouble(&DP_4);
                              free(Q);
                              free(Q2);
                              free(dQ0);
                              free(dQ1);
                              free(dQ2);
                              free(dQ2_0);
                              free(dQ2_1);
                              free(dQ2_2);
                              return(ERROR);
                          }
                          if (gradhesscost_atpixel(Gradient_im, Hessian_im, &cost_im, Difference_im, dp0_im, dp1_im, dp2_im, dp3_im, dp4_im, Weight))
                          {
                              WRITE_ERROR(return_gradhesscost, "Error returned by gradhesscost_atpixel");
                              free(Rz1);
                              free(Ry);
                              free(Rz2);
                              free(R);
                              free(DRz1);
                              free(DRy);
                              free(DRz2);
                              free(DR0);
                              free(DR1);
                              free(DR2);
                              free(mn);
                              free(Arg);
                              free(Grad_re);
                              free(Grad_im);
                              free(Gradient_re);
                              free(Gradient_im);
                              free(Hessian_re);
                              free(Hessian_im);
                              FreeVolumeDouble(&DP_0);
                              FreeVolumeDouble(&DP_1);
                              FreeVolumeDouble(&DP_2);
                              FreeVolumeDouble(&DP_3);
                              FreeVolumeDouble(&DP_4);
                              free(Q);
                              free(Q2);
                              free(dQ0);
                              free(dQ1);
                              free(dQ2);
                              free(dQ2_0);
                              free(dQ2_1);
                              free(dQ2_2);
                              return(ERROR);
                          }

                      }
                  }


                  *cost = (cost_re + cost_im) / 2.0;

                  for (i = 0L; i < 5L; i++)
                  {
                      Gradient[i] = Gradient_re[i] + Gradient_im[i];
                      for (j = 0L; j < 5L; j++)
                          Hessian[i * 5L + j] = Hessian_re[i * 5L + j] + Hessian_im[i * 5L + j];
                  }


                  for (i = 0L; i < 4L; i++)
                      for (j = i + 1L; j < 5L; j++)
                          Hessian[j * 5L + i] = Hessian[i * 5L + j];

              }

              free(Rz1);
              free(Ry);
              free(Rz2);
              free(R);
              free(DRz1);
              free(DRy);
              free(DRz2);
              free(DR0);
              free(DR1);
              free(DR2);
              free(mn);
              free(Arg);
              free(Grad_re);
              free(Grad_im);
              free(Gradient_re);
              free(Gradient_im);
              free(Hessian_re);
              free(Hessian_im);
              FreeVolumeDouble(&DP_0);
              FreeVolumeDouble(&DP_1);
              FreeVolumeDouble(&DP_2);
              FreeVolumeDouble(&DP_3);
              FreeVolumeDouble(&DP_4);
              free(Q);
              free(Q2);
              free(dQ0);
              free(dQ1);
              free(dQ2);
              free(dQ2_0);
              free(dQ2_1);
              free(dQ2_2);

              return(!ERROR);
          }/* End of return_gradhesscost */

          /* ------------------------------------------------------------------------- */
          /* Optimizer                                                                 */
          /* ------------------------------------------------------------------------- */
          int levenberg_cst(
              double beta[],
              double *alpha,
              double *cost,
              double a[],
              double *CoefRe,
              double *CoefIm,
              double *pntr_ReInp,
              double *pntr_ImInp,
              double *pntr_Weight,
              long   Nx,
              long   Ny,
              long   Nz,
              long   indx0,
              long   indy0,
              long   indz0,
              long   Mx,
              long   My,
              double *Left,
              double *Right,
              double *Q1,
              double *Q3,
              double scx1,
              double scy1,
              double S,
              long   SizeIm,
              int    DoDesProj,
              double *dftCompProj,
              double OldCost,
              double *lambda,
              double LambdaScale,
              long   *iter,
              double tol_angle,
              double tol_shift,
              int    *IteratingStop)
          {
#ifndef DBL_EPSILON
#define DBL_EPSILON 2e-16
#endif
              double   *da, epsilon = DBL_EPSILON;
              double   *u,   *v, *w;
              double   *t;
              double   wmax, thresh;
              long   i, j, ma = 5L;

              double costchanged[1];

              da = (double *)malloc((size_t)ma * sizeof(double));
              if (da == (double *)NULL)
              {
                  WRITE_ERROR(levenberg_cst, "ERROR - Not enough memory for da in levenberg_cst");
                  return(ERROR);
              }
              u = (double *)malloc((size_t)(ma * ma) * sizeof(double));
              if (u == (double *)NULL)
              {
                  free(da);
                  WRITE_ERROR(levenberg_cst, "ERROR - Not enough memory for u in levenberg_cst");
                  return(ERROR);
              }
              v = (double *)malloc((size_t)(ma * ma) * sizeof(double));
              if (v == (double *)NULL)
              {
                  free(u);
                  free(da);
                  WRITE_ERROR(levenberg_cst, "ERROR - Not enough memory for v in levenberg_cst");
                  return(ERROR);
              }
              w = (double *)malloc((size_t)ma * sizeof(double));
              if (w == (double *)NULL)
              {
                  free(v);
                  free(u);
                  free(da);
                  WRITE_ERROR(levenberg_cst, "ERROR - Not enough memory for w in levenberg_cst");
                  return(ERROR);
              }


              t = u;
              for (i = 0L; (i < ma); t += (std::ptrdiff_t)(ma + 1L), i++)
              {
                  for (j = 0L; (j < ma); alpha++, j++)
                      *u++ = -*alpha;
                  *t *= 1.0 + *lambda;
              }
              u -= (std::ptrdiff_t)(ma * ma);
              alpha -= (std::ptrdiff_t)(ma * ma);

              int Status;
              if (SingularValueDecomposition(u, ma, ma, w, v, SVDMAXITER, &Status) == ERROR)
              {
                  free(w);
                  free(v);
                  free(u);
                  free(da);
                  WRITE_ERROR(levenberg_cst, "ERROR - Unable to perform svdcmp in levenberg_cst");
                  return(ERROR);
              }
              wmax = 0.0;
              t = w + (std::ptrdiff_t)ma;
              while (--t >= w)
                  if (*t > wmax)
                      wmax = *t;
              thresh = epsilon * wmax;
              w += (std::ptrdiff_t)ma;
              j = ma;
              while (--j >= 0L)
              {
                  if (*--w < thresh)
                  {
                      *w = 0.0;
                      for (i = 0; (i < ma); i++)
                      {
                          u[i * ma + j] = 0.0;
                          v[i * ma + j] = 0.0;
                      }
                  }
              }
              if (SingularValueBackSubstitution(u, w, v, ma, ma, beta, da, &Status) == ERROR)
              {
                  free(w);
                  free(v);
                  free(u);
                  free(da);
                  WRITE_ERROR(levenberg_cst, "ERROR - Unable to perform svbksb in levenberg_cst");
                  return(ERROR);
              }


              v = (double *)memcpy(v, a, (size_t)ma * sizeof(double));
              t = v + (std::ptrdiff_t)ma;
              a += (std::ptrdiff_t)ma;
              da += (std::ptrdiff_t)ma;
              while (--t >= v)
              {
                  da--;
                  *t = *--a;
                  *t += *da;
              }


              if (return_gradhesscost(w, u, costchanged, v,
                                      CoefRe, CoefIm, pntr_ReInp, pntr_ImInp, pntr_Weight,
                                      Nx, Ny, Nz, indx0, indy0, indz0, Mx, My,
                                      Left, Right, Q1, Q3, scx1, scy1, S, SizeIm,
                                      DoDesProj, dftCompProj) == ERROR)
              {
                  WRITE_ERROR(levenberg_cst, "Error returned by total_gradhesscost");
                  free(w);
                  free(v);
                  free(u);
                  free(da);
                  return(ERROR);
              }


              (*iter)++;
              if (costchanged[0] < OldCost)
              {
                  if ((fabs(a[0] - v[0]) < tol_angle) && (fabs(a[1] - v[1]) < tol_angle) && (fabs(a[2] - v[2]) < tol_angle))
                  {
                      if ((fabs(a[3] - v[3]) < tol_shift) && (fabs(a[4] - v[4]) < tol_shift))
                      {
                          *IteratingStop = 1;
                      }
                  }
              }

              if (*cost == -1.0)
              {
                  if (costchanged[0] < OldCost)
                  {
                      for (i = 0L; (i < ma); i++)
                      {
                          for (j = 0L; (j < ma); j++)
                              alpha[i * ma + j] = u[i * ma + j];
                          beta[i] = w[i];
                          a[i] = v[i];
                      }
                      *cost = costchanged[0];
                  }

                  free(w);
                  free(u);
                  free(da);
                  free(v);
                  return(!ERROR);
              }




              for (i = 0L; (i < ma); i++)
              {
                  for (j = 0L; (j < ma); j++)
                      alpha[i * ma + j] = u[i * ma + j];
                  beta[i] = w[i];
                  a[i] = v[i];
              }
              *cost = costchanged[0];


#ifndef DBL_MIN
#define DBL_MIN 1e-26
#endif
#ifndef DBL_MAX
#define DBL_MAX 1e+26
#endif


              if (costchanged[0] < OldCost)
              {
                  if (*lambda > DBL_MIN)
                  {
                      *lambda /= LambdaScale;
                  }
                  else
                  {
                      *IteratingStop = 1;
                  }
              }
              else
              {
                  if (*lambda < DBL_MAX)
                  {
                      *lambda *= LambdaScale;
                  }
                  else
                  {
                      *IteratingStop = 1;
                  }
              }
              free(w);
              free(v);
              free(u);
              free(da);
              return(!ERROR);
          } /* End of levenberg_cst */

          /* ------------------------------------------------------------------------- */
          /* Main routines                                                             */
          /* ------------------------------------------------------------------------- */
          /* Check the input parameters ---------------------------------------------- */
          int cstregistrationCheck(struct cstregistrationStruct *Data)
          {

              if (Data == (struct cstregistrationStruct *)NULL)
              {
                  WRITE_ERROR(cstregistrationCheck, "Invalid data pointer.");
                  return(ERROR);
              }

              if ((Data->nx_ReDftVolume <= 0L) || (Data->ny_ReDftVolume <= 0L)
                  || (Data->nz_ReDftVolume <= 0L))
              {
                  WRITE_ERROR(cstregistrationCheck, " Error - ReDftVolume is missing.");
                  return(ERROR);
              }

              if ((Data->nx_ImDftVolume <= 0L) || (Data->ny_ImDftVolume <= 0L)
                  || (Data->nz_ImDftVolume <= 0L))
              {
                  WRITE_ERROR(cstregistrationCheck, " Error - ImDftVolume is missing.");
                  return(ERROR);
              }

              if ((Data->nx_ReDftImage <= 0L) || (Data->ny_ReDftImage <= 0L))
              {
                  WRITE_ERROR(cstregistrationCheck, " Error - ReDftImage is missing.");
                  return(ERROR);
              }

              if ((Data->nx_ImDftImage <= 0L) || (Data->ny_ImDftImage <= 0L))
              {
                  WRITE_ERROR(cstregistrationCheck, " Error - ImDftImage is missing.");
                  return(ERROR);
              }

              if ((Data->nx_ReDftVolume != Data->nx_ImDftVolume)
                  || (Data->ny_ReDftVolume != Data->ny_ImDftVolume)
                  || (Data->nz_ReDftVolume != Data->nz_ImDftVolume))
              {
                  WRITE_ERROR(cstregistrationCheck,
                              " Error - ReDftVolume and ImDftVolume should have same dimensions.");
                  return(ERROR);
              }

              if ((Data->nx_ReDftImage != Data->nx_ImDftImage)
                  || (Data->ny_ReDftImage != Data->ny_ImDftImage))
              {
                  WRITE_ERROR(cstregistrationCheck,
                              " Error - ReDftImage and ImDftImage should have same dimensions.");
                  return(ERROR);
              }

              if ((Data->nx_Weight <= 0L) || (Data->ny_Weight <= 0L))
              {
                  WRITE_ERROR(cstregistrationCheck, " Error - Weight image is missing.");
                  return(ERROR);
              }

              if ((Data->nx_Weight != Data->nx_ImDftImage)
                  || (Data->ny_Weight != Data->ny_ImDftImage))
              {
                  WRITE_ERROR(cstregistrationCheck,
                              " Error - Weight should have same dimensions as ReDftImage and ImDftImage.");
                  return(ERROR);
              }

              if (Data->nx_VoxelSize != 3L)
              {
                  WRITE_ERROR(cstregistrationCheck,
                              " Error - VoxelSize is a 3-element vector.");
                  return(ERROR);
              }
              if (((double) *Data->VoxelSize <= 0.0) || ((double) *(Data->VoxelSize + (std::ptrdiff_t) 1L)
                      <= 0.0) || ((double) *(Data->VoxelSize + (std::ptrdiff_t) 2L) <= 0.0))
              {
                  WRITE_ERROR(cstregistrationCheck,
                              " Error - VoxelSize must have all positive elements.");
                  return(ERROR);
              }

              if (Data->nx_PixelSize != 2L)
              {
                  WRITE_ERROR(cstregistrationCheck,
                              " Error - PixelSize is a 2-element vector.");
                  return(ERROR);
              }
              if (((double) *Data->PixelSize <= 0.0) ||
                  ((double) *(Data->PixelSize + (std::ptrdiff_t) 1L) <= 0.0))
              {
                  WRITE_ERROR(cstregistrationCheck,
                              " Error - PixelSize must have positive elements.");
                  return(ERROR);
              }

              if (Data->nx_Parameters != 5L)
              {
                  WRITE_ERROR(cstregistrationCheck,
                              " Error - Rotation is a 5-element vector.");
                  return(ERROR);
              }
              if (Data->ScaleLambda <= 0.0)
              {
                  WRITE_ERROR(cstregistrationCheck,
                              " Error - ScaleLambda must be greater than 0.0.");
                  return(ERROR);
              }
              if (Data->LambdaInitial < 0.0)
              {
                  WRITE_ERROR(cstregistrationCheck,
                              " Error - LambdaInitial must be greater than or equal to 0.0.");
                  return(ERROR);
              }
              if (Data->MaxNoIter < 0L)
              {
                  WRITE_ERROR(cstregistrationCheck,
                              " Error - MaxNoIter must be greater or equal to 0.");
                  return(ERROR);
              }
              if (Data->MaxNoFailure < 0L)
              {
                  WRITE_ERROR(cstregistrationCheck,
                              " Error - MaxNoFailure must be greater or equal to 0.");
                  return(ERROR);
              }
              if (Data->SatisfNoSuccess < 0L)
              {
                  WRITE_ERROR(cstregistrationCheck,
                              " Error - SatisfNoSuccess must be greater or equal to 0.");
                  return(ERROR);
              }
              if (Data->SatisfNoSuccess + Data->MaxNoFailure > Data->MaxNoIter)
              {
                  WRITE_ERROR(cstregistrationCheck,
                              " Error - SatisfNoSuccess + MaxNoFailure must be smaller or equal to MaxNoIter.");
                  return(ERROR);
              }
              if ((Data->MakeDesiredProj != 0L) && (Data->MakeDesiredProj != 1L))
              {
                  WRITE_ERROR(cstregistrationCheck,
                              " Error - MakeDesiredProj should be 0 or 1.");
                  return(ERROR);
              }
              if (Data->ToleranceAngle < 0.0)
              {
                  WRITE_ERROR(cstregistrationCheck, " Error - ToleranceAngle should be positive.");
                  return(ERROR);
              }
              if (Data->ToleranceShift < 0.0)
              {
                  WRITE_ERROR(cstregistrationCheck, " Error - ToleranceShift should be positive.");
                  return(ERROR);
              }
              return(!ERROR);

          }

          /* Manage problem size ----------------------------------------------------- */
          int cstregistrationSize(struct cstregistrationStruct *Data)
          {

              long   MaxIter;

              if (Data == (struct cstregistrationStruct *)NULL)
              {
                  WRITE_ERROR(cstregistrationSize, "Invalid data pointer");
                  return(ERROR);
              }

              MaxIter = Data->MaxNoIter + 1L;

              Data->nx_dftProj = Data->nx_ReDftImage;
              Data->ny_dftProj = Data->ny_ReDftImage;
              Data->nz_dftProj = 2L;

              Data->nx_OutputParameters = MaxIter;
              Data->ny_OutputParameters = 5L;

              Data->nx_Cost = MaxIter;
              Data->nx_TimePerIter = MaxIter;
              Data->nx_Failures = MaxIter;

              return(!ERROR);

          }

          /* Registration ------------------------------------------------------------ */
          int  cstregistration(struct cstregistrationStruct *Data)
          {
              const long      OrderOfSpline = 3L;
              const double   epsilon = DBL_EPSILON;

              int      Status = !ERROR, DoDesProj, IteratingStop, FlagMaxIter;

              long   MaxIter, MaxIter1, MaxIter2, iter;
              long   Mx, My, Nx, Ny, Nz, indx0, indy0, indz0, indx, indy, index;
              long   SizeIm, MaxNumberOfFailures, SatisfNumberOfSuccesses, nSuccess, nFailure;
              double  *reDftVolume, *imDftVolume, *pntr_ReInp, *pntr_ImInp, *CoefRe, *CoefIm, *par;
              double  *pntr_par0, *pntr_par1, *pntr_par2, *pntr_par3, *pntr_par4, *pntr_cost;
              double  *pntr_time, *pntr_FailureIter, *pntr_RedftCompProj, *pntr_ImdftCompProj;
              double   LambdaScale, lambda, cost, OldCost, tol_angle, tol_shift;
              double   OneIterInSeconds, *Parameters, *Gradient, *Hessian;
              double   *Q1, *Q3, *As, *Ap, *hlp;
              double   vox_x, vox_y, vox_z, pix_x, pix_y, scx1, scy1, S;
              double   sum_xx_re, sum_xx_im, dftproj_inp_re, dftproj_inp_im, sc_re, sc_im;
              time_t   time1, time2, *tp1 = NULL, *tp2 = NULL;

              if (Data == (struct cstregistrationStruct *)NULL)
              {
                  WRITE_ERROR(cstregistration, "Invalid data pointer");
                  return(ERROR);
              }

              DoDesProj = (int) Data->MakeDesiredProj;

              reDftVolume = Data->ReDftVolume;
              imDftVolume = Data->ImDftVolume;

              Nx = Data->nx_ReDftVolume;
              Ny = Data->ny_ReDftVolume;
              Nz = Data->nz_ReDftVolume;

              indx0 = -Nx / 2L;
              indy0 = -Ny / 2L;
              indz0 = -Nz / 2L;

              pntr_ReInp = Data->ReDftImage;
              pntr_ImInp = Data->ImDftImage;

              Mx = Data->nx_ReDftImage;
              My = Data->ny_ReDftImage;

              SizeIm = Mx * My;

              par = Data->VoxelSize;
              vox_x = (double) *par++;
              vox_y = (double) *par++;
              vox_z = (double) *par;

              par = Data->PixelSize;
              pix_x = (double) *par++;
              pix_y = (double) *par;

              scx1 = 2.0 * PI  / ((double) Mx * pix_x);
              scy1 = 2.0 * PI  / ((double) My * pix_y);

              S = (double)(Nx * Ny * Nz) / (double)(Mx * My);

              LambdaScale = Data->ScaleLambda;
              lambda = Data->LambdaInitial;
              cost = 0.0;

              MaxIter = Data->MaxNoIter;
              MaxIter1 = MaxIter - 1L;
              MaxIter2 = MaxIter + 1L;

              MaxNumberOfFailures = Data->MaxNoFailure;
              SatisfNumberOfSuccesses = Data->SatisfNoSuccess;
              tol_angle = Data->ToleranceAngle;
              tol_shift = Data->ToleranceShift;

              Parameters = (double *)malloc((size_t) 5L * sizeof(double));
              if (Parameters == (double *)NULL)
              {
                  WRITE_ERROR(cstregistration, "ERROR - Not enough memory for Parameters");
                  return(ERROR);
              }

              par = Data->Parameters;
              Parameters[0] = (double)(*par++) * PI / 180.0;
              Parameters[1] = (double)(*par++) * PI / 180.0;
              Parameters[2] = (double)(*par++) * PI / 180.0;
              Parameters[3] = (double)(*par++);
              Parameters[4] = (double)(*par);

              AllocateVolumeDouble(&CoefRe, Nx, Ny, Nz, &Status);
              if (Status == ERROR)
              {
                  WRITE_ERROR(cstregistration, "ERROR - Not enough memory for CoefRe");
                  free(Parameters);
                  return(ERROR);
              }

              AllocateVolumeDouble(&CoefIm, Nx, Ny, Nz, &Status);
              if (Status == ERROR)
              {
                  WRITE_ERROR(cstregistration, "ERROR - Not enough memory for CoefIm");
                  free(Parameters);
                  FreeVolumeDouble(&CoefRe);
                  return(ERROR);
              }

              ChangeBasisVolume(reDftVolume, CoefRe, Nx, Ny, Nz, CardinalSpline,
                                BasicSpline, OrderOfSpline, Periodic, epsilon, &Status);
              if (Status == ERROR)
              {
                  WRITE_ERROR(cstregistration, "ERROR");
                  free(Parameters);
                  FreeVolumeDouble(&CoefRe);
                  FreeVolumeDouble(&CoefIm);
                  return(ERROR);
              }

              ChangeBasisVolume(imDftVolume, CoefIm, Nx, Ny, Nz, CardinalSpline,
                                BasicSpline, OrderOfSpline, Periodic, epsilon, &Status);
              if (Status == ERROR)
              {
                  WRITE_ERROR(cstregistration, "ERROR");
                  free(Parameters);
                  FreeVolumeDouble(&CoefRe);
                  FreeVolumeDouble(&CoefIm);
                  return(ERROR);
              }

              Gradient = (double *)malloc((size_t) 5L * sizeof(double));
              if (Gradient == (double *)NULL)
              {
                  WRITE_ERROR(cstregistration, "ERROR - Not enough memory for Gradient");
                  FreeVolumeDouble(&CoefRe);
                  FreeVolumeDouble(&CoefIm);
                  free(Parameters);
                  return(ERROR);
              }

              Hessian = (double *)malloc((size_t) 25L * sizeof(double));
              if (Hessian == (double *)NULL)
              {
                  WRITE_ERROR(cstregistration, "ERROR - Not enough memory for Hessian");
                  FreeVolumeDouble(&CoefRe);
                  FreeVolumeDouble(&CoefIm);
                  free(Parameters);
                  free(Gradient);
                  return(ERROR);
              }

              Q1 = (double *)malloc((size_t) 16L * sizeof(double));
              if (Q1 == (double *)NULL)
              {
                  WRITE_ERROR(cstregistration, "ERROR - Not enough memory for Q1");
                  FreeVolumeDouble(&CoefRe);
                  FreeVolumeDouble(&CoefIm);
                  free(Parameters);
                  free(Gradient);
                  free(Hessian);
                  return(ERROR);
              }

              if (GetIdentitySquareMatrix(Q1, 4L) == ERROR)
              {
                  WRITE_ERROR(cstregistration, "Error returned by GetIdentitySquareMatrix");
                  FreeVolumeDouble(&CoefRe);
                  FreeVolumeDouble(&CoefIm);
                  free(Parameters);
                  free(Gradient);
                  free(Hessian);
                  free(Q1);
                  return(ERROR);
              }

              hlp = Q1;
              *hlp = (double) Nx;
              hlp += (std::ptrdiff_t)5L;
              *hlp = (double) Ny;
              hlp += (std::ptrdiff_t)5L;
              *hlp = (double) Nz;

              Q3 = (double *)malloc((size_t) 16L * sizeof(double));
              if (Q3 == (double *)NULL)
              {
                  WRITE_ERROR(cstregistration, "ERROR - Not enough memory for Q3");
                  FreeVolumeDouble(&CoefRe);
                  FreeVolumeDouble(&CoefIm);
                  free(Parameters);
                  free(Gradient);
                  free(Hessian);
                  free(Q1);
                  return(ERROR);
              }

              if (GetIdentitySquareMatrix(Q3, 4L) == ERROR)
              {
                  WRITE_ERROR(cstregistration, "Error returned by GetIdentitySquareMatrix");
                  FreeVolumeDouble(&CoefRe);
                  FreeVolumeDouble(&CoefIm);
                  free(Parameters);
                  free(Gradient);
                  free(Hessian);
                  free(Q1);
                  free(Q3);
                  return(ERROR);
              }

              hlp = Q3;
              *hlp = 1.0 / (double) Mx;
              hlp += (std::ptrdiff_t)5L;
              *hlp = 1.0 / (double) My;


              As = (double *)malloc((size_t) 16L * sizeof(double));
              if (As == (double *)NULL)
              {
                  WRITE_ERROR(cstregistration, "ERROR - Not enough memory for As");
                  FreeVolumeDouble(&CoefRe);
                  FreeVolumeDouble(&CoefIm);
                  free(Parameters);
                  free(Gradient);
                  free(Hessian);
                  free(Q1);
                  free(Q3);
                  return(ERROR);
              }

              Ap = (double *)malloc((size_t) 16L * sizeof(double));
              if (Ap == (double *)NULL)
              {
                  WRITE_ERROR(cstregistration, "ERROR - Not enough memory for Ap");
                  FreeVolumeDouble(&CoefRe);
                  FreeVolumeDouble(&CoefIm);
                  free(Parameters);
                  free(Gradient);
                  free(Hessian);
                  free(Q1);
                  free(Q3);
                  free(As);
                  return(ERROR);
              }

              if (GetIdentitySquareMatrix(As, 4L) == ERROR)
              {
                  WRITE_ERROR(cstregistration, "Error returned by GetIdentitySquareMatrix");
                  FreeVolumeDouble(&CoefRe);
                  FreeVolumeDouble(&CoefIm);
                  free(Parameters);
                  free(Gradient);
                  free(Hessian);
                  free(Q1);
                  free(Q3);
                  free(As);
                  free(Ap);
                  return(ERROR);
              }

              hlp = As;
              *hlp = vox_x;
              hlp += (std::ptrdiff_t)5L;
              *hlp = vox_y;
              hlp += (std::ptrdiff_t)5L;
              *hlp = vox_z;


              if (GetIdentitySquareMatrix(Ap, 4L) == ERROR)
              {
                  WRITE_ERROR(cstregistration, "Error returned by GetIdentitySquareMatrix");
                  FreeVolumeDouble(&CoefRe);
                  FreeVolumeDouble(&CoefIm);
                  free(Parameters);
                  free(Gradient);
                  free(Hessian);
                  free(Q1);
                  free(Q3);
                  free(As);
                  free(Ap);
                  return(ERROR);
              }

              hlp = Ap;
              *hlp = 1.0 / pix_x;
              hlp += (std::ptrdiff_t)5L;
              *hlp = 1.0 / pix_y;

              // Normalize by the std. dev. the input dft image
              sum_xx_re = 0.0;
              sum_xx_im = 0.0;
              for (indy = 0L; indy < My; indy++)
              {
                  for (indx = 0L; indx < Mx; indx++)
                  {
                      index = indy * Mx + indx;
                      dftproj_inp_re = (double) pntr_ReInp[index];
                      dftproj_inp_im = (double) pntr_ImInp[index];
                      sum_xx_re += dftproj_inp_re * dftproj_inp_re;
                      sum_xx_im += dftproj_inp_im * dftproj_inp_im;
                  }
              }

              index = (My / 2L) * Mx + (Mx / 2L);
              dftproj_inp_re = (double) pntr_ReInp[index];
              dftproj_inp_im = (double) pntr_ImInp[index];

              sc_re = 1.0 / sqrt(sum_xx_re + sum_xx_im - dftproj_inp_re * dftproj_inp_re - dftproj_inp_im * dftproj_inp_im);
              sc_im = sc_re;

              for (indy = 0L; indy < My; indy++)
              {
                  for (indx = 0L; indx < Mx; indx++)
                  {
                      index = indy * Mx + indx;
                      dftproj_inp_re = (double) pntr_ReInp[index];
                      dftproj_inp_im = (double) pntr_ImInp[index];
                      pntr_ReInp[index] = (sc_re * dftproj_inp_re);
                      pntr_ImInp[index] = (sc_im * dftproj_inp_im);
                  }
              }

              pntr_RedftCompProj = Data->dftProj;
              pntr_ImdftCompProj = pntr_RedftCompProj + (std::ptrdiff_t) SizeIm;
              for (indy = 0L; indy < My; indy++)
              {
                  for (indx = 0L; indx < Mx; indx++)
                  {
                      index = indy * Mx + indx;
                      pntr_RedftCompProj[index] = 0.0;
                      pntr_ImdftCompProj[index] = 0.0;
                  }
              }

              pntr_par0 = Data->OutputParameters;
              pntr_par1 = pntr_par0 + (std::ptrdiff_t) MaxIter2;
              pntr_par2 = pntr_par1 + (std::ptrdiff_t) MaxIter2;
              pntr_par3 = pntr_par2 + (std::ptrdiff_t) MaxIter2;
              pntr_par4 = pntr_par3 + (std::ptrdiff_t) MaxIter2;
              pntr_cost = Data->Cost;
              pntr_time = Data->TimePerIter;
              pntr_FailureIter = Data->Failures;
              for (indx = 0L; indx < MaxIter2; indx++)
              {
                  *pntr_par0++ = 0.0;
                  *pntr_par1++ = 0.0;
                  *pntr_par2++ = 0.0;
                  *pntr_par3++ = 0.0;
                  *pntr_par4++ = 0.0;
                  *pntr_cost++ = 0.0;
                  *pntr_time++ = 0.0;
                  *pntr_FailureIter++ = 0.0;
              }

              time1 = time(tp1);

              if (return_gradhesscost(Gradient, Hessian, &cost, Parameters,
                                      CoefRe, CoefIm, pntr_ReInp, pntr_ImInp, Data->Weight,
                                      Nx, Ny, Nz, indx0, indy0, indz0, Mx, My,
                                      Ap, As, Q1, Q3, scx1, scy1, S, SizeIm,
                                      DoDesProj, Data->dftProj) == ERROR)
              {
                  WRITE_ERROR(cstregistration, "Error returned by return_gradhesscost");
                  FreeVolumeDouble(&CoefRe);
                  FreeVolumeDouble(&CoefIm);
                  free(Parameters);
                  free(Gradient);
                  free(Hessian);
                  free(Q1);
                  free(Q3);
                  free(As);
                  free(Ap);
                  return(ERROR);
              }

              time2 = time(tp2);
              OneIterInSeconds = difftime(time2, time1);


              if (DoDesProj)
              {
                  FreeVolumeDouble(&CoefRe);
                  FreeVolumeDouble(&CoefIm);
                  free(Parameters);
                  free(Gradient);
                  free(Hessian);
                  free(Q1);
                  free(Q3);
                  free(As);
                  free(Ap);
                  return(!ERROR);
              }

              pntr_par0 = Data->OutputParameters;
              pntr_par1 = pntr_par0 + (std::ptrdiff_t) MaxIter2;
              pntr_par2 = pntr_par1 + (std::ptrdiff_t) MaxIter2;
              pntr_par3 = pntr_par2 + (std::ptrdiff_t) MaxIter2;
              pntr_par4 = pntr_par3 + (std::ptrdiff_t) MaxIter2;

              pntr_cost = Data->Cost;
              pntr_time = Data->TimePerIter;

              pntr_FailureIter = Data->Failures;

              *pntr_par0++ = Parameters[0];
              *pntr_par1++ = Parameters[1];
              *pntr_par2++ = Parameters[2];
              *pntr_par3++ = Parameters[3];
              *pntr_par4++ = Parameters[4];
              *pntr_cost++ = cost;
              *pntr_time++ = OneIterInSeconds;
              pntr_FailureIter++;

              if ((MaxIter == 0L) && (!DoDesProj))
              {
                  FreeVolumeDouble(&CoefRe);
                  FreeVolumeDouble(&CoefIm);
                  free(Parameters);
                  free(Gradient);
                  free(Hessian);
                  free(Q1);
                  free(Q3);
                  free(As);
                  free(Ap);
                  return(!ERROR);
              }

              nSuccess = 0L;
              nFailure = 0L;
              OldCost   = cost;
              iter = - 1L;
              IteratingStop = 0;
              FlagMaxIter = (MaxIter != 1L);
              if (!FlagMaxIter)
                  cost = - 1.0;

              do
              {
                  time1 = time(tp1);

                  if (levenberg_cst(Gradient, Hessian, &cost, Parameters,
                                    CoefRe, CoefIm, pntr_ReInp, pntr_ImInp, Data->Weight,
                                    Nx, Ny, Nz, indx0, indy0, indz0, Mx, My,
                                    Ap, As, Q1, Q3, scx1, scy1, S, SizeIm,
                                    DoDesProj, Data->dftProj, OldCost, &lambda, LambdaScale,
                                    &iter, tol_angle, tol_shift, &IteratingStop) == ERROR)
                  {
                      WRITE_ERROR(cstregistration, "Error returned by levenberg_cst");
                      FreeVolumeDouble(&CoefRe);
                      FreeVolumeDouble(&CoefIm);
                      free(Parameters);
                      free(Gradient);
                      free(Hessian);
                      free(Q1);
                      free(Q3);
                      free(As);
                      free(Ap);
                      return(ERROR);
                  }

                  time2 = time(tp2);
                  OneIterInSeconds = difftime(time2, time1);

                  *pntr_par0++ = Parameters[0];
                  *pntr_par1++ = Parameters[1];
                  *pntr_par2++ = Parameters[2];
                  *pntr_par3++ = Parameters[3];
                  *pntr_par4++ = Parameters[4];
                  *pntr_cost++ = cost;
                  *pntr_time++ = OneIterInSeconds;

                  if (cost < OldCost)
                  {
                      OldCost = cost;
                      nSuccess++;
                      pntr_FailureIter++;
                      if (nSuccess >= SatisfNumberOfSuccesses)
                      {
                          break;
                      }
                      if (IteratingStop)
                      {
                          break;
                      }
                  }
                  else
                  {
                      nFailure++;
                      *pntr_FailureIter++ = (iter + 1L);
                  }

              }
              while ((nFailure <= MaxNumberOfFailures) && (iter < MaxIter1) && FlagMaxIter);

              *(Data->NumberIterPerformed) = (iter + 1L);
              *(Data->NumberSuccPerformed) = nSuccess;
              *(Data->NumberFailPerformed) = nFailure;

              FreeVolumeDouble(&CoefRe);
              FreeVolumeDouble(&CoefIm);
              free(Parameters);
              free(Gradient);
              free(Hessian);
              free(Q1);
              free(Q3);
              free(As);
              free(Ap);

              return(!ERROR);
          }
