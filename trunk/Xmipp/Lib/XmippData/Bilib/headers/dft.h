/**@name Fourier Transform */
//@{

/*--------------------------------------------------------------------------*/
/** Amplitude/Phase --> Real Imaginary.
    Converts an (amplitude, phase) representation of a complex signal
    into a (real, imaginary) representation the input phase is in [rad].
    in-place processing. The input signal (amplitude Am2Re and phase Ph2Im) is
    replaced by the output signal (real Am2Re and imaginary Ph2Im)
    SignalLength is the signal length.

    success: return(!ERROR); failure: return(ERROR);
*/
extern int		AmplitudePhaseToRealImaginary
				(
					double	Am2Re[],			/* (amplitude -> real) */
					double	Ph2Im[],			/* (phase -> imaginary) */
					long	SignalLength		/* signal length */
				);

/*--------------------------------------------------------------------------*/
/** DFT of a complex signal.

    Computes the direct DFT of a complex signal given in
    (amplitude, phase) representation and returns an (amplitude, phase)
    representation. The input phase is in [rad]. The output phase is in
    [rad]; its domain is (-PI, PI). The origin is at index [0].

    The highest coordinate (SignalLength-2)/2 is at index
    [(SignalLength-2)/2] for SignalLength even. The highest coordinate
    (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd.
    
    The lowest coordinate -SignalLength/2 is at index [SignalLength/2]
    for SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at
    index [(SignalLength+1)/2] for SignalLength odd.
    
    The coordinate -1 is at index [SignalLength-1] for SignalLength even or
    odd.
    
    In-place processing. The input signal (amplitude Am2Am and phase Ph2Ph) is
    replaced by the output signal (amplitude Am2Am and phase Ph2Ph).
    
    (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each
    the values returned in (TmpRe, TmpIm) are meaningless.
    
    CaS is an input array of coefficients of size SignalLength
    (see GetCaS function).
    CaS is modified internally, but is restored when
    DftAmplitudePhaseToAmplitudePhase returns.

    SignalLength is the length of the signal.
    no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5.
    
    success: return(!ERROR); failure: return(ERROR);
*/
extern int		DftAmplitudePhaseToAmplitudePhase
				(
					double	Am2Am[],			/* amplitude -> amplitude */
					double	Ph2Ph[],			/* phase -> phase */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				);

/*--------------------------------------------------------------------------*/
/** DFT of a complex signal.

    Computes the direct DFT of a complex signal given in (amplitude, phase)
    representation and returns a (real, imaginary) representation. The input
    phase is in [rad]. The origin is at index [0].
    
    The highest coordinate (SignalLength-2)/2 is at index
    [(SignalLength-2)/2] for SignalLength even. The highest coordinate
    (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd.
    
    The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
    SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
    [(SignalLength+1)/2] for SignalLength odd.
    
    The coordinate -1 is at index [SignalLength-1] for SignalLength even or odd.
    
    In-place processing. The input signal (amplitude Am2Re and phase Ph2Im) is
    replaced by the output signal (real Am2Re and imaginary Ph2Im).
    
    (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each.
    The values returned in (TmpRe, TmpIm) are meaningless.
    
    CaS is an input array of coefficients of size SignalLength
    (see GetCaS function). CaS is modified internally, but is restored when
    DftAmplitudePhaseToRealImaginary returns.
    
    SignalLength is the length of the signal.
    No restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5.
    
    success: return(!ERROR); failure: return(ERROR);
*/
extern int		DftAmplitudePhaseToRealImaginary
				(
					double	Am2Re[],			/* amplitude -> real */
					double	Ph2Im[],			/* phase -> imaginary */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				);

/*--------------------------------------------------------------------------*/
/** DFT of a complex signal.

    Computes the direct DFT of a complex signal given in (real, imaginary)
    representation and returns an (amplitude, phase) representation. The
    output phase is in [rad]; its domain is (-PI, PI). The origin is at index
    [0].
    
    The highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2]
    for SignalLength even. The highest coordinate (SignalLength-1)/2 is at
    index [(SignalLength-1)/2] for SignalLength odd.
    
    The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
    SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
    [(SignalLength+1)/2] for SignalLength odd.
    
    The coordinate -1 is at index [SignalLength-1] for SignalLength even or
    odd.
    
    In-place processing. The input signal (real Re2Am and imaginary Im2Ph) is
    replaced by the output signal (amplitude Re2Am and phase Im2Ph).
    
    (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each.
    The values returned in (TmpRe, TmpIm) are meaningless.
    
    CaS is an input array of coefficients of size SignalLength
    (see GetCaS function). CaS is modified internally, but is restored
    when DftRealImaginaryToAmplitudePhase returns.
    
    SignalLength is the length of the signal. No restriction on SignalLength,
    but best efficiency for radix 2, 3, 4, 5.
    
    success: return(!ERROR); failure: return(ERROR);
*/
extern int		DftRealImaginaryToAmplitudePhase
				(
					double	Re2Am[],			/* real -> amplitude */
					double	Im2Ph[],			/* imaginary -> phase */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				);

/*--------------------------------------------------------------------------*/
/** DFT of a complex signal.

    Computes the direct DFT of a complex signal given in (real, imaginary)
    representation and returns a (real, imaginary) representation. The origin
    is at index [0].
    
    The highest coordinate (SignalLength-2)/2 is at index
    [(SignalLength-2)/2] for SignalLength even. The highest coordinate
    (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd.
    
    The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
    SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
    [(SignalLength+1)/2] for SignalLength odd.
    
    The coordinate -1 is at index [SignalLength-1] for SignalLength even or
    odd.
    
    In-place processing. The input signal (real Re2Re and imaginary Im2Im) is
    replaced by the output signal (real Re2Re and imaginary Im2Im).
    
    (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each.
    The values returned in (TmpRe, TmpIm) are meaningless.
    
    CaS is an input array of coefficients of size SignalLength
    (see GetCaS function). CaS is modified internally, but is restored
    when DftRealImaginaryToRealImaginary returns
    
    SignalLength is the length of the signal. No restriction on SignalLength,
    but best efficiency for radix 2, 3, 4, 5
    
    success: return(!ERROR); failure: return(ERROR);
*/
extern int		DftRealImaginaryToRealImaginary
				(
					double	Re2Re[],			/* real -> real */
					double	Im2Im[],			/* imaginary -> imaginary */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				);

/*--------------------------------------------------------------------------*/
/** DFT of a real signal.

    Computes the direct DFT of a real signal and returns an
    (amplitude, phase) representation. The output phase is in [rad];
    its domain is (-PI, PI). The origin is at index [0].
    
    The highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2]
    for SignalLength even. The highest coordinate (SignalLength-1)/2 is at
    index [(SignalLength-1)/2] for SignalLength odd.
    
    The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
    SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
    [(SignalLength+1)/2] for SignalLength odd.
    
    The coordinate -1 is at index [SignalLength-1] for SignalLength even or
    odd.
    
    In-place processing. Te input signal (real R2Am) is	replaced by the
    output signal (amplitude R2Am and phase PhOut).
    
    Tmp is a pre-allocated workspace of size SignalLength. 
    The values returned in Tmp are meaningless.
    
    CaS is an input array of coefficients of size SignalLength
    (see GetCaS function). CaS is modified internally, but is restored
    when DftRealToAmplitudePhase returns.
    
    SignalLength is the length of the signal. No restriction on SignalLength,
    but best efficiency for radix 2, 3, 4, 5
    
    success: return(!ERROR); failure: return(ERROR);
*/
extern int		DftRealToAmplitudePhase
				(
					double	R2Am[],				/* real -> amplitude */
					double	PhOut[],			/* output phase */
					double	*Tmp,				/* scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				);

/*--------------------------------------------------------------------------*/
/** DFT of a real signal.

    Computes the direct DFT of a real signal and returns a
    (real, imaginary) representation. The origin is at index [0].
    
    The highest coordinate (SignalLength-2)/2 is at index
    [(SignalLength-2)/2] for SignalLength even. The highest coordinate
    (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd.
    
    The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
    SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
    [(SignalLength+1)/2] for SignalLength odd.
    
    The coordinate -1 is at index [SignalLength-1] for SignalLength even or odd.
    
    In-place processing. The input signal (real R2Re) is replaced by the
    output signal (real R2Re and imaginary ImOut).
    
    Tmp is a pre-allocated workspace of size SignalLength. The values returned
    in Tmp are meaningless.
    
    CaS is an input array of coefficients of size SignalLength
    (see GetCaS function). CaS is modified internally, but is restored
    when DftRealToRealImaginary returns.
    
    SignalLength is the length of the signal. No restriction on SignalLength,
    but best efficiency for radix 2, 3, 4, 5
    
    success: return(!ERROR); failure: return(ERROR);
*/
extern int		DftRealToRealImaginary
				(
					double	R2Re[],				/* real -> real */
					double	ImOut[],			/* output imaginary */
					double	*Tmp,				/* scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				);

/*--------------------------------------------------------------------------*/
/** Compute the CaS array necessary for Fourier transforms.

    Computes an array of coefficients of size SignalLength. These
    coefficients are necessary for performing a Hartley transform. The
    Hartley transform is an alternate representation of Fourier for real
    signals. Hartley computations are more accurate (less roundoff errors)
    than Fourier. The same coefficients are used for direct and inverse
    transforms.
    
    success: return(!ERROR); failure: return(ERROR);
*/
extern int		GetCaS
				(
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				);

/*--------------------------------------------------------------------------*/
/** Inverse Fourier Transform of a complex signal.

    Computes the inverse DFT of a complex signal given in (amplitude, phase)
    representation and returns an (amplitude, phase) representation.
    The input phase is in [rad]. The output phase is in [rad];
    its domain is (-PI, PI). The origin is at index [0].
    
    The highest coordinate (SignalLength-2)/2 is at index
    [(SignalLength-2)/2] for SignalLength even. The highest coordinate
    (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd.
    
    The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
    SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at
    index [(SignalLength+1)/2] for SignalLength odd.
    
    The coordinate -1 is at index [SignalLength-1] for SignalLength even or odd.
    
    In-place processing. The input signal (amplitude Am2Am and phase Ph2Ph) is
    replaced by the output signal (amplitude Am2Am and phase Ph2Ph).
    
    (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each.
    The values returned in (TmpRe, TmpIm) are meaningless.
    
    CaS is an input array of coefficients of size SignalLength
    (see GetCaS function). CaS is modified internally, but is restored when
    DftAmplitudePhaseToAmplitudePhase returns.
    
    SignalLength is the length of the signal. No restriction on SignalLength,
    but best efficiency for radix 2, 3, 4, 5
    
    success: return(!ERROR); failure: return(ERROR);
*/
extern int		InvDftAmplitudePhaseToAmplitudePhase
				(
					double	Am2Am[],			/* amplitude -> amplitude */
					double	Ph2Ph[],			/* phase -> phase */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				);

/*--------------------------------------------------------------------------*/
/** Inverse Fourier Transform of a real signal.

    Computes the inverse DFT of a complex signal given in (amplitude, phase)
    representation and returns a real signal. The complex Fourier signal is
    symmetrized before the inverse transformation is applied. The input phase
    is in [rad]. The origin is at index [0].
    
    The highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2]
    for SignalLength even. The highest coordinate (SignalLength-1)/2 is at
    index [(SignalLength-1)/2] for SignalLength odd.
    
    The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
    SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
    [(SignalLength+1)/2] for SignalLength odd.
    
    The coordinate -1 is at index [SignalLength-1] for SignalLength even or odd.
    
    In-place processing. The input signal (amplitude Am2R and phase PhIn) is
    replaced by the output signal (real Am2R).
    
    Tmp is a pre-allocated workspace of size SignalLength. The values returned
    in Tmp are meaningless.
    
    CaS is an input array of coefficients of size SignalLength
    (see GetCaS function). CaS is modified internally, but is restored when
    InvDftAmplitudePhaseToReal returns.
    
    SignalLength is the length of the signal. No restriction on SignalLength,
    but best efficiency for radix 2, 3, 4, 5
    
    success: return(!ERROR); failure: return(ERROR);
*/
extern int		InvDftAmplitudePhaseToReal
				(
					double	Am2R[],				/* amplitude -> real */
					double	PhIn[],				/* input phase */
					double	*Tmp,				/* scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				);

/*--------------------------------------------------------------------------*/
/** Inverse Fourier Transform of a complex signal.

    Computes the inverse DFT of a complex signal given in (amplitude, phase)
    representation and returns a (real, imaginary) representation. The input
    phase is in [rad]. The origin is at index [0].
    
    The highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2]
    for SignalLength even. The highest coordinate (SignalLength-1)/2 is at
    index [(SignalLength-1)/2] for SignalLength odd.
    
    The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
    SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
    [(SignalLength+1)/2] for SignalLength odd.
    
    The coordinate -1 is at index [SignalLength-1] for SignalLength even or odd.

    In-place processing. The input signal (amplitude Am2Re and phase Ph2Im) is
    replaced by the output signal (real Am2Re and imaginary Ph2Im).
    
    (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each.
    Tthe values returned in (TmpRe, TmpIm) are meaningless.
    
    CaS is an input array of coefficients of size SignalLength
    (see GetCaS function). CaS is modified internally, but is restored when
    InvDftAmplitudePhaseToRealImaginary returns.
    
    SignalLength is the length of the signal.
    No restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5.
    
    success: return(!ERROR); failure: return(ERROR);
*/
extern int		InvDftAmplitudePhaseToRealImaginary
				(
					double	Am2Re[],			/* amplitude -> real */
					double	Ph2Im[],			/* phase -> imaginary */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				);

/*--------------------------------------------------------------------------*/
/** Inverse Fourier Transform for a complex signal.

    Computes the inverse DFT of a complex signal given in (real, imaginary)
    representation and returns an (amplitude, phase) representation.
    The output phase is in [rad]; its domain is (-PI, PI). The origin is
    at index [0].
    
    The highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2]
    for SignalLength even. The highest coordinate (SignalLength-1)/2 is at
    index [(SignalLength-1)/2] for SignalLength odd.
    
    The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
    SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
    [(SignalLength+1)/2] for SignalLength odd.
    
    The coordinate -1 is at index [SignalLength-1] for SignalLength even or odd.

    In-place processing. The input signal (real Re2Am and imaginary Im2Ph) is
    replaced by the output signal (amplitude Re2Am and phase Im2Ph).
    
    (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each.
    The values returned in (TmpRe, TmpIm) are meaningless.
    
    CaS is an input array of coefficients of size SignalLength
    (see GetCaS function). CaS is modified internally, but is restored when
    InvDftRealImaginaryToAmplitudePhase returns.
    
    SignalLength is the length of the signal. No restriction on SignalLength,
    but best efficiency for radix 2, 3, 4, 5
    
    success: return(!ERROR); failure: return(ERROR);
*/
extern int		InvDftRealImaginaryToAmplitudePhase
				(
					double	Re2Am[],			/* real -> amplitude */
					double	Im2Ph[],			/* imaginary -> phase */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				);

/*--------------------------------------------------------------------------*/
/** Inverse Fourier Transform of a real signal.

    Computes the inverse DFT of a complex signal given in (real, imaginary)
    representation and returns a real signal. The complex Fourier signal is
    symmetrized before the inverse transformation is applied. The origin is
    at index [0].
    
    The highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2]
    for SignalLength even. The highest coordinate (SignalLength-1)/2 is at
    index [(SignalLength-1)/2] for SignalLength odd.
    
    The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
    SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
    [(SignalLength+1)/2] for SignalLength odd.
    
    The coordinate -1 is at index [SignalLength-1] for SignalLength even or odd.
    
    In-place processing. The input signal (real Re2R and imaginary ImIn) is
    replaced by the output signal (real Re2R).
    
    Tmp is a pre-allocated workspace of size SignalLength. The values returned
    in Tmp are meaningless.
    
    CaS is an input array of coefficients of size SignalLength
    (see GetCaS function). CaS is modified internally, but is restored when
    InvDftRealImaginaryToReal returns.
    
    SignalLength is the length of the signal. No restriction on SignalLength,
    but best efficiency for radix 2, 3, 4, 5
    
    success: return(!ERROR); failure: return(ERROR);
*/
extern int		InvDftRealImaginaryToReal
				(
					double	Re2R[],				/* real -> real */
					double	ImIn[],				/* input imaginary */
					double	*Tmp,				/* scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				);

/*--------------------------------------------------------------------------*/
/** Inverse Fourier Transform of a complex signal.

    Computes the inverse DFT of a complex signal given in (real, imaginary)
    representation and returns a (real, imaginary) representation.
    The origin is at index [0].
    
    The highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2]
    for SignalLength even. The highest coordinate (SignalLength-1)/2 is at
    index [(SignalLength-1)/2] for SignalLength odd.
    
    The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
    SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
    [(SignalLength+1)/2] for SignalLength odd.
    
    The coordinate -1 is at index [SignalLength-1] for SignalLength even or odd.

    In-place processing. The input signal (real Re2Re and imaginary Im2Im) is
    replaced by the output signal (real Re2Re and imaginary Im2Im).
    
    (TmpRe, TmpIm) are pre-allocated workspaces of size SignalLength each.
    The values returned in (TmpRe, TmpIm) are meaningless.
    
    CaS is an input array of coefficients of size SignalLength
    (see GetCaS function). CaS is modified internally, but is restored when
    InvDftRealImaginaryToRealImaginary returns.
    
    SignalLength is the length of the signal. No restriction on SignalLength,
    but best efficiency for radix 2, 3, 4, 5
    
    success: return(!ERROR); failure: return(ERROR);
*/
extern int		InvDftRealImaginaryToRealImaginary
				(
					double	Re2Re[],			/* real -> real */
					double	Im2Im[],			/* imaginary -> imaginary */
					double	*TmpRe,				/* first scratch workspace */
					double	*TmpIm,				/* second scratch workspace */
					double	CaS[],				/* Hartley transform coefficients */
					long	SignalLength		/* signal length */
				);

/*--------------------------------------------------------------------------*/
/** Real/Imaginary --> Amplitude/Phase.

    Converts a (real, imaginary) representation of a complex signal
    into an (amplitude, phase) representation  The output phase is in [rad];
    its domain is (-PI, PI).
    
    In-place processing. The input signal (real Re2Am and imaginary Im2Ph) is
    replaced by the output signal (amplitude Re2Am and phase Im2Ph).
    
    SignalLength is the signal length.
    
    success: return(!ERROR); failure: return(ERROR);
*/
extern int		RealImaginaryToAmplitudePhase
				(
					double	Re2Am[],			/* real -> amplitude */
					double	Im2Ph[],			/* imaginary -> phase */
					long	SignalLength		/* signal length */
				);

/*--------------------------------------------------------------------------*/
/** Volume Amplitude/Phase --> Real/Imaginary.
   Converts an (amplitude, phase) representation of a complex signal
   into a (real, imaginary) representation.
   
   The input phase is in [rad]. In-place processing .
   The input signal (amplitude Am2Re and phase Ph2Im) is
   replaced by the output signal (real Am2Re and imaginary Ph2Im).
   
   Nx is the width of the volume. Ny is the height of the volume.
   Nz is the depth of the volume. 
   
   success: return(!ERROR); failure: return(ERROR); */
extern int		VolumeAmplitudePhaseToRealImaginary
				(
					float	Am2Re[],			/* (amplitude -> real) */
					float	Ph2Im[],			/* (phase -> imaginary) */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz					/* depth of the volume */
				);

/*--------------------------------------------------------------------------*/
/** Direct DFT of a complex signal.
   Computes the direct DFT of a complex signal given in (amplitude, phase)
   representation and returns an (amplitude, phase) representation.
   
   The input phase is in [rad]. The output phase is in [rad]; its domain
   is (-PI, PI). The origin is at index [0].
   
   The highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2]
   for SignalLength even. The highest coordinate (SignalLength-1)/2 is at
   index [(SignalLength-1)/2] for SignalLength odd.
   The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
   SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
   [(SignalLength+1)/2] for SignalLength odd.
   The coordinate -1 is at index [SignalLength-1] for SignalLength even or
   odd.
   
   No restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5.
   In the explanations above, SignalLength has to be replaced by Nx, Ny, Nz.
   
   Nx is the width of the volume. Ny is the height of the volume.
   Nz is the depth of the volume.
   
   In-place processing. The input signal (amplitude Am2Am and phase Ph2Ph) is
   replaced by the output signal (amplitude Am2Am and phase Ph2Ph).
   
   success: return(!ERROR); failure: return(ERROR);
   The returned value is duplicated in Status */
extern int		VolumeDftAmplitudePhaseToAmplitudePhase
				(
					float	*Am2Am,				/* amplitude -> amplitude */
					float	*Ph2Ph,				/* phase -> phase */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				);

/*--------------------------------------------------------------------------*/
/** Computes the direct DFT of a complex signal.
   Computes the direct DFT of a complex signal given in (amplitude, phase)
   representation and returns a (real, imaginary) representation.
   
   The input phase is in [rad]. The origin is at index [0].
   The highest coordinate (SignalLength-2)/2 is at index
   [(SignalLength-2)/2] for SignalLength even. The highest coordinate
   (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd.
   The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
   SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
   [(SignalLength+1)/2] for SignalLength odd. The coordinate -1 is at index
   [SignalLength-1] for SignalLength even or odd.
   
   No restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5.
   In the explanations above, SignalLength has to be replaced by Nx, Ny, Nz.
   Nx is the width of the volume. Ny is the height of the volume.
   Nz is the depth of the volume.
   
   In-place processing, the input signal (amplitude Am2Re and phase Ph2Im) is
   replaced by the output signal (real Am2Re and imaginary Ph2Im).
   
   success: return(!ERROR); failure: return(ERROR);
   The returned value is duplicated in Status */
extern int		VolumeDftAmplitudePhaseToRealImaginary
				(
					float	*Am2Re,				/* amplitude -> real */
					float	*Ph2Im,				/* phase -> imaginary */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				);

/*--------------------------------------------------------------------------*/
/** Computes the direct DFT of a complex signal.
   Computes the direct DFT of a complex signal given in (real, imaginary)
   representation and returns an (amplitude, phase) representation.
   
   The output phase is in [rad]; its domain is (-PI, PI).
   The origin is at index [0]. The highest coordinate (SignalLength-2)/2 is
   at index [(SignalLength-2)/2] for SignalLength even. The highest coordinate
   (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd.
   The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
   SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
   [(SignalLength+1)/2] for SignalLength odd. The coordinate -1 is at index
   [SignalLength-1] for SignalLength even or odd.
   
   No restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5.
   In the explanations above, SignalLength has to be replaced by Nx, Ny, Nz.
   Nx is the width of the volume. Ny is the height of the volume.
   Nz is the depth of the volume.
   
   In-place processing. The input signal (real Re2Am and imaginary Im2Ph) is
   replaced by the output signal (amplitude Re2Am and phase Im2Ph).
   
   Success: return(!ERROR); failure: return(ERROR);
   The returned value is duplicated in Status */
extern int		VolumeDftRealImaginaryToAmplitudePhase
				(
					float	*Re2Am,				/* real -> amplitude */
					float	*Im2Ph,				/* imaginary -> phase */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				);

/*--------------------------------------------------------------------------*/
/** Computes the direct DFT of a complex signal.
   Computes the direct DFT of a complex signal given in (real, imaginary)
   representation and returns a (real, imaginary) representation.
   
   The origin is at index [0]. The highest coordinate (SignalLength-2)/2 is
   at index [(SignalLength-2)/2] for SignalLength even. The highest coordinate
   (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd.
   The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
   SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
   [(SignalLength+1)/2] for SignalLength odd. The coordinate -1 is at index
   [SignalLength-1] for SignalLength even or odd. 
   
   No restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5.
   In the explanations above, SignalLength has to be replaced by Nx, Ny, Nz.
   Nx is the width of the volume. Ny is the height of the volume.
   Nz is the depth of the volume.
   
   In-place processing. The input signal (real Re2Re and imaginary Im2Im) is
   replaced by the output signal (real Re2Re and imaginary Im2Ph).
   
   Success: return(!ERROR); failure: return(ERROR);.
   The returned value is duplicated in Status */
extern int		VolumeDftRealImaginaryToRealImaginary
				(
					float	*Re2Re,				/* real -> real */
					float	*Im2Im,				/* imaginary -> imaginary */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				);

/*--------------------------------------------------------------------------*/
/** Computes the direct DFT of a real signal.
   Computes the direct DFT of a real signal and returns an (amplitude, phase)
   representation. The output phase is in [rad]; its domain is (-PI, PI).
   
   The origin is at index [0]. The highest coordinate (SignalLength-2)/2 is
   at index [(SignalLength-2)/2] for SignalLength even. The highest coordinate
   (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd.
   The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
   SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
   [(SignalLength+1)/2] for SignalLength odd. The coordinate -1 is at index
   [SignalLength-1] for SignalLength even or odd.
   
   No restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5.
   In the explanations above, SignalLength has to be replaced by Nx, Ny, Nz.
   Nx is the width of the volume. Ny is the height of the volume.
   Nz is the depth of the volume.
   
   In-place processing. The input signal Re2Am is replaced by the output
   signal (amplitude Re2Am and phase PhOut).
   
   success: return(!ERROR); failure: return(ERROR);
   The returned value is duplicated in Status */
extern int		VolumeDftRealToAmplitudePhase
				(
					float	*Re2Am,				/* real -> amplitude */
					float	*PhOut,				/* output phase */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				);

/*--------------------------------------------------------------------------*/
/** Computes the direct DFT of a real signal.
    Computes the direct DFT of a real signal and returns a (real, imaginary)
    representation.
    
    The origin is at index [0]. The highest coordinate (SignalLength-2)/2 is
    at index [(SignalLength-2)/2] for SignalLength even. The highest
    coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for
    SignalLength odd. The lowest coordinate -SignalLength/2 is at index
    [SignalLength/2] for SignalLength even. The lowest coordinate
    -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd.
    The coordinate -1 is at index [SignalLength-1] for SignalLength even or
    odd.
    
    No restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5.
    In the explanations above, SignalLength has to be replaced by Nx, Ny, Nz.
    Nx is the width of the volume. Ny is the height of the volume.
    Nz is the depth of the volume.
    
    In-place processing. The input signal Re2Re is
    replaced by the output signal (real Re2Re and imaginary ImOut).
    
    success: return(!ERROR); failure: return(ERROR);
    The returned value is duplicated in Status */
extern int		VolumeDftRealToRealImaginary
				(
					float	*Re2Re,				/* real -> real */
					float	*ImOut,				/* output imaginary */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				);

/*--------------------------------------------------------------------------*/
/** Computes the inverse DFT of a complex signal.
    Computes the inverse DFT of a complex signal given in (amplitude, phase)
    representation and returns an (amplitude, phase) representation.
    The input phase is in [rad]. The output phase is in [rad]; its domain is
    (-PI, PI).
    
    The origin is at index [0]. The highest coordinate (SignalLength-2)/2 is
    at index [(SignalLength-2)/2] for SignalLength even. The highest
    coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for
    SignalLength odd. The lowest coordinate -SignalLength/2 is at index
    [SignalLength/2] for SignalLength even. The lowest coordinate
    -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd.
    The coordinate -1 is at index [SignalLength-1] for SignalLength even or
    odd.
    
    No restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5.
    In the explanations above, SignalLength has to be replaced by Nx, Ny, Nz.
    Nx is the width of the volume. Ny is the height of the volume.
    Nz is the depth of the volume.
    
    In-place processing. The input signal (amplitude Am2Am and phase Ph2Ph) is
    replaced by the output signal (amplitude Am2Am and phase Ph2Ph).
    
    success: return(!ERROR); failure: return(ERROR);
    The returned value is duplicated in Status */
extern int		VolumeInvDftAmplitudePhaseToAmplitudePhase
				(
					float	*Am2Am,				/* amplitude -> amplitude */
					float	*Ph2Ph,				/* phase -> phase */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				);

/*--------------------------------------------------------------------------*/
/** Computes the inverse DFT of a complex signal.
    Computes the inverse DFT of a complex signal given in (amplitude, phase)
    representation and returns a real signal.
    The complex Fourier signal is symmetrized before the inverse
    transformation is applied. The input phase is in [rad].
    The origin is at index [0].
    
    The highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2]
    for SignalLength even. The highest coordinate (SignalLength-1)/2 is at
    index [(SignalLength-1)/2] for SignalLength odd. The lowest coordinate
    -SignalLength/2 is at index [SignalLength/2] for SignalLength even.
    The lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2]
    for SignalLength odd. The coordinate -1 is at index [SignalLength-1] for
    SignalLength even or odd.
    
    No restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5.
    In the explanations above, SignalLength has to be replaced by Nx, Ny, Nz.
    Nx is the width of the volume. Ny is the height of the volume.
    Nz is the depth of the volume.
    
    In-place processing. The input signal (amplitude Am2Re and phase PhIn) is
    replaced by the output signal (real Am2Re). PhIn is destroyed.
    
    success: return(!ERROR); failure: return(ERROR);
    The returned value is duplicated in Status */
extern int		VolumeInvDftAmplitudePhaseToReal
				(
					float	*Am2Re,				/* amplitude -> real */
					float	*PhIn,				/* input phase */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				);

/*--------------------------------------------------------------------------*/
/** Computes the inverse DFT of a complex signal.
    Computes the inverse DFT of a complex signal given in (amplitude, phase)
    representation and returns a (real, imaginary) representation.
    
    The input phase is in [rad]. The origin is at index [0].
    
    The highest coordinate (SignalLength-2)/2 is at index [(SignalLength-2)/2]
    for SignalLength even. The highest coordinate (SignalLength-1)/2 is at
    index [(SignalLength-1)/2] for SignalLength odd. The lowest coordinate
    -SignalLength/2 is at index [SignalLength/2] for SignalLength even.
    The lowest coordinate -(SignalLength-1)/2 is at index [(SignalLength+1)/2]
    for SignalLength odd. The coordinate -1 is at index [SignalLength-1] for
    SignalLength even or odd.
    
    No restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5.
    In the explanations above, SignalLength has to be replaced by Nx, Ny, Nz.
    Nx is the width of the volume. Ny is the height of the volume.
    Nz is the depth of the volume.
    
    In-place processing. The input signal (amplitude Am2Re and phase Ph2Im) is
    replaced by the output signal (real Am2Re and imaginary Ph2Im).
    
    success: return(!ERROR); failure: return(ERROR);
    The returned value is duplicated in Status */
extern int		VolumeInvDftAmplitudePhaseToRealImaginary
				(
					float	*Am2Re,				/* amplitude -> real */
					float	*Ph2Im,				/* phase -> imaginary */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				);

/*--------------------------------------------------------------------------*/
/** Computes the inverse DFT of a complex signal.
    Computes the inverse DFT of a complex signal given in (real, imaginary)
    representation and returns an (amplitude, phase) representation.
    The output phase is in [rad]; its domain is (-PI, PI).
    
    The origin is at index [0]. The highest coordinate (SignalLength-2)/2 is
    at index [(SignalLength-2)/2] for SignalLength even. The highest coordinate
    (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd.
    The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
    SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
    [(SignalLength+1)/2] for SignalLength odd. The coordinate -1 is at index
    [SignalLength-1] for SignalLength even or odd.
    
    No restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5.
    In the explanations above, SignalLength has to be replaced by Nx, Ny, Nz.
    Nx is the width of the volume. Ny is the height of the volume.
    Nz is the depth of the volume.
    
    In-place processing. The input signal (real Re2Am and imaginary Im2Ph) is
    replaced by the output signal (amplitude Re2Am and phase Im2Ph).
    
    success: return(!ERROR); failure: return(ERROR);
    The returned value is duplicated in Status */
extern int		VolumeInvDftRealImaginaryToAmplitudePhase
				(
					float	*Re2Am,				/* real -> amplitude */
					float	*Im2Ph,				/* imaginary -> phase */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				);

/*--------------------------------------------------------------------------*/
/** Computes the inverse DFT of a complex signal.
   Computes the inverse DFT of a complex signal given in (real, imaginary)
   representation and returns a real signal.
   The complex Fourier signal is symmetrized before the inverse transformation
   is applied.
   
   The origin is at index [0]. The highest coordinate (SignalLength-2)/2 is
   at index [(SignalLength-2)/2] for SignalLength even. The highest coordinate
   (SignalLength-1)/2 is at index [(SignalLength-1)/2] for SignalLength odd.
   The lowest coordinate -SignalLength/2 is at index [SignalLength/2] for
   SignalLength even. The lowest coordinate -(SignalLength-1)/2 is at index
   [(SignalLength+1)/2] for SignalLength odd. The coordinate -1 is at index
   [SignalLength-1] for SignalLength even or odd.
   
   No restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5.
   In the explanations above, SignalLength has to be replaced by Nx, Ny, Nz.
   Nx is the width of the volume. Ny is the height of the volume.
   Nz is the depth of the volume.
   
   In-place processing. The input signal (real Re2Re and imaginary ImIn) is
   replaced by the output signal (real Re2Re).
   
   success: return(!ERROR); failure: return(ERROR);
   The returned value is duplicated in Status */
extern int		VolumeInvDftRealImaginaryToReal
				(
					float	*Re2Re,				/* real -> real */
					float	*ImIn,				/* input imaginary */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				);

/*--------------------------------------------------------------------------*/
/** Computes the inverse DFT of a complex signal.
    Computes the inverse DFT of a complex signal given in (real, imaginary)
    representation and returns a (real, imaginary) representation.
    
    The origin is at index [0]. The highest coordinate (SignalLength-2)/2 is
    at index [(SignalLength-2)/2] for SignalLength even. The highest
    coordinate (SignalLength-1)/2 is at index [(SignalLength-1)/2] for
    SignalLength odd. The lowest coordinate -SignalLength/2 is at index
    [SignalLength/2] for SignalLength even. The lowest coordinate
    -(SignalLength-1)/2 is at index [(SignalLength+1)/2] for SignalLength odd.
    The coordinate -1 is at index [SignalLength-1] for SignalLength even or
    odd.
    
    No restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5.
    In the explanations above, SignalLength has to be replaced by Nx, Ny, Nz.
    Nx is the width of the volume. Ny is the height of the volume.
    Nz is the depth of the volume.
    
    In-place processing. The input signal (real Re2Re and imaginary Im2Im) is
    replaced by the output signal (real Re2Re and imaginary Im2Ph).
    
    success: return(!ERROR); failure: return(ERROR);
    The returned value is duplicated in Status */
extern int		VolumeInvDftRealImaginaryToRealImaginary
				(
					float	*Re2Re,				/* real -> real */
					float	*Im2Im,				/* imaginary -> imaginary */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					int		*Status				/* error management */
				);

/*--------------------------------------------------------------------------*/
/* Real/Imaginary --> Amplitude/Phase.
   Converts a (real, imaginary) representation of a complex signal
   into an (amplitude, phase) representation.
   The output phase is in [rad]; its domain is (-PI, PI).
   
   In-place processing. The input signal (real Re2Am and imaginary Im2Ph) is
   replaced by the output signal (amplitude Re2Am and phase Im2Ph).
   
   Nx is the width of the volume. Ny is the height of the volume.
   Nz is the depth of the volume.
   
   success: return(!ERROR); failure: return(ERROR);*/
extern int		VolumeRealImaginaryToAmplitudePhase
				(
					float	Re2Am[],			/* real -> amplitude */
					float	Im2Ph[],			/* imaginary -> phase */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz					/* depth of the volume */
				);
//@}