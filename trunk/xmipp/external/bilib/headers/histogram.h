/**@defgroup BilibHistograms Histograms
   @ingroup BilibLibrary */
//@{
/*--------------------------------------------------------------------------*/
/** Build histogram.
    Computation of the frequencies of occurence of data values.
    VolumeSource is a (float)volume of size (Nx x Ny x Nz).
    Value[] is a (float)array of length previously determined by
    HistogramGetSize. The returned content of Value[] is sorted in strict
    ascending order. Frequency[] is a (double)array of length previously
    determined by HistogramGetSize.

    success: return(!ERROR); failure: return(ERROR)
*/
extern int  HistogramBuild
(
    float *VolumeSource,  /* data to process */
    long Nx,     /* width of the volume */
    long Ny,     /* height of the volume */
    long Nz,     /* depth of the volume */
    double Frequency[],  /* output vector of ordinates */
    float Value[],   /* output vector of abscissa */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
/** Equalize histogram.
  * Construction of the lookup table: Value[k] <-> EqualizedValue[k].
  * EqualizedValue[] satisfies:
  * @code
  *    sum(k in K(n0)) Frequency[k] ~=~ sum(k in K(n1)) Frequency[k]
  *    for all n0, n1 in [0L, NumberOfClasses - 1L]
  * @endcode
  * where ~=~ means "is about equal to", and where K(n) is a domain such that
  * @code
  *     EqualizedValue[k0] = (sum(k1 in K(n)) Frequency[k1] * Value[k1])
  *        / (sum(k2 in K(n)) Frequency[k2])
  *     for all k0 in K(n)
  * @endcode
  * under the constraint
  * @code
  *    DistinctElements(QuantizedValue[]) == NumberOfClasses
  * @endcode
  *
  * Frequency[] is a (double)array of length HistogramLength.
  * The content of Frequency[] must be strictly positive.
  * The content of Frequency[] must have unit sum.
  * Value[] is a (float)array of length HistogramLength.
  * The content of Value[] must be sorted in strictly ascending order.
  * EqualizedValue[] is a returned (float)array of length HistogramLength.
  *
  * On input, NumberOfClasses indicates the desired number of classes.
  * On output, NumberOfClasses returns the effective number of classes.
  * NumberOfClasses is no greater than (1.0 / max(Frequency[])), and never
  * increases.
  *
  * It may happen that the only solution that satisfies all constraints is
  * undesirable
  * e.g., Frequency[] = {0.9, 0.1};
  *       Value[] = {10.0F, 90.0F};
  *       NumberOfClasses = 2L; (desired)
  * results in
  *       QuantizedValues[] = {18.0F, 18.0F};
  *       NumberOfClasses = 1L; (actual)
  *
  * success: return(!ERROR); failure: return(ERROR)
*/
extern int  HistogramEqualize
(
    double Frequency[],  /* histogram ordinates */
    float Value[],   /* histogram abscissa */
    float EqualizedValue[], /* output vector of abscissa */
    long HistogramLength, /* length of the histogram */
    long *NumberOfClasses, /* number of classes, desired -> actual */
    double Tolerance   /* admissible relative error */
);

/*--------------------------------------------------------------------------*/
/** Get size of histogram.
    Determination of the number of differing data values in a volume.
    VolumeSource is a (float)volume of size (Nx x Ny x Nz).

    success: return(!ERROR); failure: return(ERROR) */
extern int  HistogramGetSize
(
    float *VolumeSource,  /* data to process */
    long Nx,     /* width of the volume */
    long Ny,     /* height of the volume */
    long Nz,     /* depth of the volume */
    long *HistogramLength, /* output length of the histogram */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
/** Histogram K means.
    Construction of the lookup table: Value[k] <-> QuantizedValue[k].
    Minimization of sum(k) Frequency[k] * (Value[k] - QuantizedValue[k])^2.
    under the constraint DistinctElements(QuantizedValue[]) == NumberOfClasses.

    Frequency[] is a (double)array of length HistogramLength. The content of
    Frequency[] must be strictly positive. The content of Frequency[] must have
    unit sum. Value[] is a (float)array of length HistogramLength.
    The content of Value[] must be sorted in strictly ascending order.
    QuantizedValue[] is a returned (float)array of length HistogramLength.

    On input, NumberOfClasses indicates the desired number of classes.
    On output, NumberOfClasses returns the effective number of classes.
    NumberOfClasses never increases.

    Important cases that go undetected (unfortunately):
    1. convergence to a non-global optimum
         e.g., Frequency[] = {0.25, 0.25, 0.25, 0.25};
               Value[] = {1.0F, 2.0F, 10.0F, 20.0F};
               NumberOfClasses = 3L;
         results in the local optimum
               QuantizedValues[] = {1.0F, 2.0F, 15.0F, 15.0F};
         instead of the true global optimum
               QuantizedValues[] = {1.5F, 1.5F, 10.0F, 20.0F};

    2. there may be more than one global optimum
         e.g., Frequency[] = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
               Value[] = {-1.0F, 0.0F, 1.0F};
               NumberOfClasses = 2L;
         results in the global optimum
               QuantizedValues[] = {-0.5F, -0.5F, 1.0F};
         the other global optimum is ignored
               QuantizedValues[] = {-1.0F, 0.5F, 0.5F};

    success: return(!ERROR); failure: return(ERROR) */
extern int  HistogramKMeans
(
    double Frequency[],  /* histogram ordinates */
    float Value[],   /* histogram abscissa */
    float QuantizedValue[], /* output vector of abscissa */
    long HistogramLength, /* length of the histogram */
    long *NumberOfClasses, /* number of classes, desired -> actual */     double Tolerance,   /* admissible relative error */
    int  *Status    /* error management */
);
//@}
