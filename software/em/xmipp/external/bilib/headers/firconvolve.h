/**@defgroup FIRConvolve FIR Convolve
   @ingroup BilibLibrary */
//@{
/*--------------------------------------------------------------------------*/
/** FIR Convolve 1D.
 * The specified boundary convention applies to the input data only,
 *  not to the kernel. The boundary convention applied to the kernel
 *  is FiniteDataSupport. The input and the output have the same length.
 *  The origin for the kernel is given with respect to the leftmost sample [0].
 *
 *  success: return(!ERROR); failure: return(ERROR)
 *
 *  General structure is as follows:
 *  @code
 *   for (i = 0L; (i < SignalLength); i++) {
 *       Sum = 0.0;
 *       for (j = -Infinity; (j <= Infinity); j++)
 *            Sum += InputData[j] * Kernel[KernelOrigin + i - j];
 *       OutputData[i] = Sum;
 *   }
 *  @endcode
*/
extern int  FirConvolve
    (
        double InputData[],  /* data to process */
        double OutputData[],  /* result */
        long SignalLength,  /* length of the 1D data array */
        double Kernel[],   /* kernel */
        long KernelOrigin,  /* center of the kernel */
        long KernelLength,  /* length of the 1D kernel */
        enum TBoundaryConvention
        Convention   /* boundary convention */
    );

/*--------------------------------------------------------------------------*/
/** FIR Convolve Antisymmetric 1D.
  * The specified boundary convention applies to the input data only,
  * not to the kernel. The boundary convention applied to the kernel is
  * FiniteDataSupport. The input and the output have the same length.
  * The origin for the kernel is its leftmost sample; it corresponds to its
  * anti-symmetry axis. The value of the kernel at the origin is expected to be
  * HalfKernel[0] = 0.0. The full length of the symmetric kernel is
  * (2L * KernelHalfLength - 1L)
  *
  * success: return(!ERROR); failure: return(ERROR)
  *
  * General structure is as follows:
  * @code
  *  for (i = 0L; (i < SignalLength); i++) {
  *     Sum = 0.0;
  *     for (j = 1L; (j < KernelHalfLength); j++)
  *          Sum += (InputData[i + j] - InputData[i - j]) * HalfKernel[j];
  *     OutputData[i] = Sum;
  *  }
  * @endcode
*/
extern int  FirConvolveAntiSymmetric
    (
        double InputData[],  /* data to process */
        double OutputData[],  /* result */
        long SignalLength,  /* length of the 1D data array */
        double HalfKernel[],  /* causal part of the kernel */
        long KernelHalfLength, /* length of the causal part of the kernel */
        enum TBoundaryConvention
        Convention   /* boundary convention */
    );

/*--------------------------------------------------------------------------*/
/** FIR Convolve Antisymmetric 3D.
  * The specified boundary convention applies to the input data only,
  * not to the kernel. The boundary convention applied to the kernel is
  * FiniteDataSupport. VolumeSource is a (double)volume of size (Nx x Ny x Nz).
  * OutputData is a (double)volume of size (Nx x Ny x Nz). The origin for the
  * kernel is its leftmost sample; it corresponds to its anti-symmetry axis.
  * The full length of the amti-symmetric kernel is (2L * KernelHalfLength - 1L).
  * The 1D kernel is applied successively to each principal direction in a
  * separable fashion.
  *
  * success: return(!ERROR); failure: return(ERROR)
  *
  * General structure is as follows:
  * @code
  *  for (i = 0L; (i < SignalLength); i++) {
  *     Sum = InputData[i] * HalfKernel[0];
  *     for (j = 1L; (j < KernelHalfLength); j++)
  *         Sum += (InputData[i + j] - InputData[i - j]) * HalfKernel[j];
  *     OutputData[i] = Sum;
  *  }
  * @endcode
*/
extern int  FirConvolveAntiSymmetricVolume
    (
        double *VolumeSource,  /* data to process */
        double *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        double HalfKernel[],  /* causal part of the kernel */
        long KernelHalfLength, /* length of the causal part of the kernel */
        enum TBoundaryConvention
        Convention   /* boundary convention */
    );

/*--------------------------------------------------------------------------*/
/** FIR Convolve Symmetric 1D.
  * The specified boundary convention applies to the input data only, not
  * to the kernel. The boundary convention applied to the kernel is
  * FiniteDataSupport. The input and the output have the same length.
  * The origin for the kernel is its leftmost sample; it corresponds to its
  * symmetry axis. The full length of the symmetric kernel is
  * (2L * KernelHalfLength - 1L).
  *
  * success: return(!ERROR); failure: return(ERROR)
  *
  * General structure is as follows:
  * @code
  *  for (i = 0L; (i < SignalLength); i++) {
  *      Sum = InputData[i] * HalfKernel[0];
  *      for (j = 1L; (j < KernelHalfLength); j++)
  *          Sum += (InputData[i + j] + InputData[i - j]) * HalfKernel[j];
  *      OutputData[i] = Sum;
  *  }
  *  @endcode
*/

extern int  FirConvolveSymmetric
    (
        double InputData[],  /* data to process */
        double OutputData[],  /* result */
        long SignalLength,  /* length of the 1D data array */
        double HalfKernel[],  /* causal part of the kernel */
        long KernelHalfLength, /* length of the causal part of the kernel */
        enum TBoundaryConvention
        Convention   /* boundary convention */
    );

/*--------------------------------------------------------------------------*/
/** FIR Convolve Symmetric 3D.
 * The specified boundary convention applies to the input data only,
 *  not to the kernel. The boundary convention applied to the kernel is
 *  FiniteDataSupport. VolumeSource is a (double)volume of size (Nx x Ny x Nz).
 *  OutputData is a (double)volume of size (Nx x Ny x Nz). The origin for the
 *  kernel is its leftmost sample; it corresponds to its symmetry axis.
 *  The full length of the symmetric kernel is (2L * KernelHalfLength - 1L).
 *  The 1D kernel is applied successively to each principal direction in a
 *  separable fashion.
 *
 *  success: return(!ERROR); failure: return(ERROR)
 *
 *  General structure is as follows:
 *  @code
 *   for (i = 0L; (i < SignalLength); i++) {
 *     Sum = InputData[i] * HalfKernel[0];
 *     for (j = 1L; (j < KernelHalfLength); j++)
 *         Sum += (InputData[i + j] + InputData[i - j]) * HalfKernel[j];
 *     OutputData[i] = Sum;
 *   }
 *   @endcode
*/
extern int  FirConvolveSymmetricVolume
    (
        double *VolumeSource,  /* data to process */
        double *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        double HalfKernel[],  /* causal part of the kernel */
        long KernelHalfLength, /* length of the causal part of the kernel */
        enum TBoundaryConvention
        Convention   /* boundary convention */
    );
//@}

