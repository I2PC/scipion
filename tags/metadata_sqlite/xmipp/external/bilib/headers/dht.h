/**@defgroup HartleyTransform Hartley Transform
   @ingroup BilibLibrary */
//@{
/*--------------------------------------------------------------------------*/
/** Discrete Hartley transform of a real signal.

    The input signal is given by Data (double). The computation is in-place
    (the output replaces the input).

    ScratchBuffer is a pre-allocated workspace of size SignalLength.
    The values returned in ScratchBuffer are meaningless.

    CaS is an input array of coefficients of size SignalLength
    (see GetCaS function). CaS is modified internally, but is restored
    when DHT returns.

    SignalLength is the length of the signal. No restriction on SignalLength,
    but best efficiency for radix 2, 3, 4, 5.

    Performs a DHT transform when (Forward == TRUE). Performs an inverse DHT
    transform when (Forward == FALSE).

    success: return(!ERROR); failure: return(ERROR);
*/
extern int  DiscreteHartleyTransform
    (
        double Data[],    /* signal */
        double *ScratchBuffer,  /* scratch buffer */
        double CaS[],    /* coefficients */
        long SignalLength,  /* signal length */
        int  Forward    /* direction */
    );
//@}

