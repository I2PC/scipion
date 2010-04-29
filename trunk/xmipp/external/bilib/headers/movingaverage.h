/**@defgroup BilibMA Moving average
   @ingroup BilibLibrary */
//@{
/*--------------------------------------------------------------------------*/
/** Moving average (1D).
  * The specified boundary convention applies to the input data only, not to
  * the kernel. The boundary convention applied to the kernel is
  * FiniteDataSupport.
  *
  * The kernel has a constant value and sums up to 1.0.
  * The input and the output have the same length.
  * The origin for the kernel is given with respect to the leftmost sample [0].
  *
  * success: return(!ERROR); failure: return(ERROR)
  *
  * General structure is as follows:
  * @code
  *  for (i = 0L; (i < SignalLength); i++) {
  *    Sum = 0.0;
  *    for (j = -Infinity; (j <= Infinity); j++)
  *       Sum += InputData[j] * Kernel[KernelOrigin + i - j];
  *    OutputData[i] = Sum;
  *  }
  * @endcode
*/
extern int  MovingAverage
    (
        double InputData[],  /* data to process */
        double OutputData[],  /* result */
        long SignalLength,  /* length of the 1D data array */
        long KernelOrigin,  /* center of the kernel */
        long KernelLength,  /* length of the 1D kernel */
        enum TBoundaryConvention
        Convention   /* boundary convention */
    );

/*--------------------------------------------------------------------------*/
/** Moving average (3D).
  * The specified boundary convention applies to the input data only,
  * not to the kernel. The boundary convention applied to the kernel is
  * FiniteDataSupport.
  * VolumeSource is a (float)volume of size (Nx x Ny x Nz).
  * OutputData is a (float)volume of size (Nx x Ny x Nz).
  * The origin for the kernel is given with respect to the leftmost sample [0].
  * The 1D kernel is applied successively to each principal direction in a
  * separable fashion.
  *
  * success: return(!ERROR); failure: return(ERROR).
  *
  * General structure is as follows:
  * @code
 *  for (i = 0L; (i < SignalLength); i++) {
 *     Sum = 0.0;
 *     for (j = i - KernelOrigin; (j < (KernelLength + i - KernelOrigin)); j++)
 *          Sum += InputData[j];
 *     OutputData[i] = Sum / KernelLength;
 *   }
 * @endcode
*/
extern int  MovingAverageVolume
    (
        float *VolumeSource,  /* data to process */
        float *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long KernelOrigin,  /* center of the kernel */
        long KernelLength,  /* length of the 1D kernel */
        enum TBoundaryConvention
        Convention   /* boundary convention */
    );
//@}
