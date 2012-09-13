/**@defgroup BilibGradient Gradient
   @ingroup BilibLibrary */
//@{
/*--------------------------------------------------------------------------*/
/** Compute the gradient of 1D data.
    InputData and OutputGradient have size SignalLength.

    success: return(!ERROR); failure: return(ERROR)
*/
extern int  LinearGradient
    (
        double InputData[],  /* input 1D data */
        double OutputGradient[], /* output 1D gradient */
        long SignalLength,  /* signal length */
        long Degree,    /* degree of the spline space */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        double Tolerance,   /* admissible relative error */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Compute the gradient of 2D data.
    InputImage has size (Nx x Ny). OutputGradient is an array of two
    elements. Each element of OutputGradient is a (float *) pointer to an
    image. Each element has size (Nx x Ny). The first element is the gradient
    along x. The second element is the gradient along y.

    success: return(!ERROR); failure: return(ERROR)
*/
extern int  PlanarGradient
    (
        float *InputImage,  /* input 2D data */
        float *OutputGradient[], /* output 2D gradient array [x, y] */
        long Nx,     /* width of the image */
        long Ny,     /* height of the image */
        long Degree,    /* degree of the spline space */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        double Tolerance,   /* admissible relative error */
        int  *Status    /* error management */
    );

/*--------------------------------------------------------------------------*/
/** Compute the gradient of 3D data.
    InputVolume has size (Nx x Ny x Nz). OutputGradient is an array of
    three elements. Each element of OutputGradient is a (float *) pointer
    to a volume. Each element of OutputGradient has size (Nx x Ny x Nz).
    The first element is the gradient along x. The second element is the
    gradient along y. The third element is the gradient along z.

    success: return(!ERROR); failure: return(ERROR)
*/
extern int  VolumetricGradient
    (
        float *InputVolume,  /* input 3D data */
        float *OutputGradient[], /* output 3D gradient array [x, y, z] */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        long Degree,    /* degree of the spline space */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        double Tolerance,   /* admissible relative error */
        int  *Status    /* error management */
    );
//@}
