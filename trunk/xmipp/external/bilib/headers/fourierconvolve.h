/**@defgroup FourierConvolve Fourier Convolve
   @ingroup BilibLibrary */
//@{
/*--------------------------------------------------------------------------*/
/** Fourier convolution with any kernel.
    Same conventions as for OneConvolveFourier.
    Preprocessing steps have to be carrried out elsewhere.

    success: return(!ERROR); failure: return(ERROR)
*/
extern int  ManyConvolveFourier
(
    double Data[],    /* data for in-place processing */
    double KernelDht[],  /* discrete Hartley transform of the kernel */
    double CaS[],    /* Hartley transform coefficients */
    double *ScratchBuffer,  /* scratch buffer */
    long SignalLength,  /* signal length */
    double Shift    /* additional translation */
);

/*--------------------------------------------------------------------------*/
/** Fourier convolution with a symmetric kernel.
    Same conventions as for OneConvolveFourierSymmetricKernel.
    Preprocessing steps have to be carrried out elsewhere

    success: return(!ERROR); failure: return(ERROR)
*/
extern int  ManyConvolveFourierSymmetricKernel
(
    double Data[],    /* data for in-place processing */
    double KernelDht[],  /* discrete Hartley transform of the kernel */
    double CaS[],    /* Hartley transform coefficients */
    double *ScratchBuffer,  /* scratch buffer */
    long SignalLength,  /* signal length */
    double Shift    /* additional translation */
);

/*--------------------------------------------------------------------------*/
/** Fourier convolution 1D.
    Fourier convolution with any kernel. The kernel has an infinite,
    periodic (SignalLength) support. The kernel origin (hot spot) is at
    index [0].

    The highest coordinate (SignalLength-2)/2 is at index
    [(SignalLength-2)/2] for SignalLength even. The highest coordinate
    (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd.

    The lowest coordinate -SignalLength/2 is at index [SignalLength/2]
    for SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at
    index [(SignalLength+1)/2] for SignalLength odd.

    The coordinate -1 is at index [SignalLength-1] for SignalLength even or odd.

    The signal has an infinite, periodic (SignalLength) support. The output
    has same length SignalLength than the input. The output has already
    been allocated. The result of the convolution overwrites the input.
    The inverse DHT of the kernel overwrites the kernel. If the kernel is
    finite support, don't forget to pad!

    success: return(!ERROR); failure: return(ERROR)

    The returned value is duplicated in Status
*/
extern int  OneConvolveFourier
(
    double Data[],    /* data for in-place processing */
    double Kernel[],   /* kernel for in-place processing */
    long SignalLength,  /* signal length */
    double Shift,    /* additional translation */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
/** Fourier convolution with a symmetrical kernel.
    The kernel has an infinite, periodic (SignalLength) support. The kernel
    origin (hot spot) is at index [0]. The highest coordinate
    (SignalLength-2)/2 is at index [(SignalLength-2)/2] for SignalLength even.
    For symmetry, kernel[(SignalLength-2)/2] = 0 for SignalLength even.

    The highest coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2]
    for SignalLength odd. No special symmetry requirement for SignalLength odd.

    The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
    SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at
    index [(SignalLength+1)/2] for SignalLength odd.

    The coordinate -1 is at index [SignalLength-1] for SignalLength even or odd.

    The signal has an infinite, periodic (SignalLength) support. The output has
    the same length SignalLength than the input. The output has already been
    allocated.

    The kernel is even-symmetric around 0 (k[x] = k[-x]). Only the coord. >= 0
    of the kernel need be given.

    The result of the convolution overwrites the input. The inverse DHT of
    the kernel overwrites the kernel. If the kernel is finite support, don't
    forget to pad!

    Observe the separate symmetries for SignalLength odd/even!
    For a symmetric kernel, DHT <=> DFT

    success: return(!ERROR); failure: return(ERROR)

    The returned value is duplicated in Status
*/
extern int  OneConvolveFourierSymmetricKernel
(
    double Data[],    /* data for in-place processing */
    double SymmetricKernel[], /* symmetric kernel for in-place processing */
    long SignalLength,  /* signal length */
    double Shift,    /* additional translation */
    int  *Status    /* error management */
);
//@}
