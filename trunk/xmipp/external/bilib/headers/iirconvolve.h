/**@name IIR Convolve */
//@{
/*--------------------------------------------------------------------------*/
/** Canonic progressive.
    Computes OutputData[i] = InputData[i] + SUM(k): OutputData[i + k + 1] *
       Kernel[k].
    The input and the output have the same length.
    The kernel and the initial output values have the same length.

    success: return(!ERROR); failure: return(ERROR) */
extern int  IirConvolveCanonicProgressive
    (
        double InputData[],  /* data to process */
        double OutputData[],  /* result */
        long SignalLength,  /* length of the 1D data array */
        double Kernel[],   /* kernel */
        double RightInit[],  /* progressive recursion initialization */
        long KernelLength  /* length of the 1D kernel */
    );

/*--------------------------------------------------------------------------*/
/** Canonic regressive.
    Computes OutputData[i] = InputData[i] + SUM(k): OutputData[i - k - 1]
       * Kernel[k]. The input and the output have the same length.
    The kernel and the initial output values have the same length.

    success: return(!ERROR); failure: return(ERROR) */
extern int  IirConvolveCanonicRegressive
    (
        double InputData[],  /* data to process */
        double OutputData[],  /* result */
        long SignalLength,  /* length of the 1D data array */
        double Kernel[],   /* kernel */
        double LeftInit[],   /* regressive recursion initialization */
        long KernelLength  /* length of the 1D kernel */
    );

/*--------------------------------------------------------------------------*/
/** The filter is defined by its poles (1D).
    The input and the output have the same length. In-place processing is
    allowed. No more than two poles are allowed for a finite support boundary.

    success: return(!ERROR); failure: return(ERROR) */
extern int  IirConvolvePoles
    (
        double InputData[],  /* data to process */
        double OutputData[],  /* result */
        long SignalLength,  /* length of the 1D data array */
        double RealPoles[],  /* array of real poles */
        long PoleNumber,   /* number of poles */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        double Tolerance   /* admissible relative error */
    );

/*--------------------------------------------------------------------------*/
/** The filter is defined by its poles (3D).
    No more than two poles are allowed for a finite support boundary.
    VolumeSource is a (double)volume of size (Nx x Ny x Nz).
    OutputData is a (double)volume of size (Nx x Ny x Nz).

    success: return(!ERROR); failure: return(ERROR) */
extern int  IirConvolvePolesVolume
    (
        double *VolumeSource,  /* data to process */
        double *VolumeDestination, /* result */
        long Nx,     /* width of the volume */
        long Ny,     /* height of the volume */
        long Nz,     /* depth of the volume */
        double RealPoles[],  /* array of real poles */
        long PoleNumber,   /* number of poles */
        enum TBoundaryConvention
        Convention,   /* boundary convention */
        double Tolerance,   /* admissible relative error */
        int  *Status    /* error management */
    );
//@}

