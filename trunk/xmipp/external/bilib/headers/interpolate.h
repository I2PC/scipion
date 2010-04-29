/**@defgroup SplineInterpolation Spline interpolation
   @ingroup BilibLibrary */
//@{
/*--------------------------------------------------------------------------*/
/** Interpolation 1D.
    Interpolation of the 1D array LineCoeff[] of length LineLength.
    Result = SUM(k): LineCoeff[k] * B-spline(Degree, Argument - k).
    Convention must be consistent with the B-spline coefficients stored
    in LineCoeff[]

    success: return(!ERROR); failure: return(ERROR) */
extern int  SplineInterpolateLine
(
    double LineCoeff[],  /* B-spline coefficients to interpolate */
    long LineLength,   /* length of the line */
    double Argument,   /* input abscissa */
    double *Result,   /* output ordinate */
    long Degree,    /* degree of the spline */
    enum TBoundaryConvention
    Convention,   /* boundary convention */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
/** Interpolation 2D.
    See coments on the 1D interpolation. */
extern int  SplineInterpolateImage
(
    float *ImageCoeff,  /* B-spline coefficients to interpolate */
    long Nx,     /* width of the image */
    long Ny,     /* height of the image */
    double Xargument,   /* input X abscissa */
    double Yargument,   /* input Y abscissa */
    double *Result,   /* output ordinate */
    long Degree,    /* degree of the spline */
    enum TBoundaryConvention
    Convention,   /* boundary convention */
    int  *Status    /* error management */
);

/*--------------------------------------------------------------------------*/
/** Interpolation 3D.
    See coments on the 1D interpolation. */
extern int  SplineInterpolateVolume
(
    float *VolumeCoeff,  /* B-spline coefficients to interpolate */
    long Nx,     /* width of the volume */
    long Ny,     /* height of the volume */
    long Nz,     /* depth of the volume */
    double Xargument,   /* input X abscissa */
    double Yargument,   /* input Y abscissa */
    double Zargument,   /* input Z abscissa */
    double *Result,   /* output ordinate */
    long Degree,    /* degree of the spline */
    enum TBoundaryConvention
    Convention,   /* boundary convention */
    int  *Status    /* error management */
);
//@}
