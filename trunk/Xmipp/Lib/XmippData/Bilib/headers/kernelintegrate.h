/**@name Kernel integrals */
//@{
/*--------------------------------------------------------------------------*/
/** Computes the integral(-Infinity, Argument) for a Blu interpolant function.
    success: return(!ERROR); failure: return(ERROR); */
extern int		BlipIntegrate
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Blu interpolant function
    of degree 0 (order 1) */
extern double	Blip00Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Blu interpolant function
    of degree 1 (order 2) */
extern double	Blip01Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Blu interpolant function
    of degree 3 (order 4) */
extern double	Blip03Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Computes the integral(-Infinity, Argument) for a Basic spline function
    of degree Degree.
    success: return(!ERROR); failure: return(ERROR); */
extern int		BsplineIntegrate
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Basic spline function
    of degree 0 (order 1) */
extern double	Bspline00Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Basic spline function
    of degree 1 (order 2) */
extern double	Bspline01Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Basic spline function
    of degree 2 (order 3) */
extern double	Bspline02Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Basic spline function
    of degree 3 (order 4) */
extern double	Bspline03Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Basic spline function
    of degree 4 (order 5) */
extern double	Bspline04Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Basic spline function
    of degree 5 (order 6) */
extern double	Bspline05Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Basic spline function
    of degree 6 (order 7) */
extern double	Bspline06Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Basic spline function
    of degree 7 (order 8) */
extern double	Bspline07Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Basic spline function
    of degree 8 (order 9) */
extern double	Bspline08Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Basic spline function
    of degree 9 (order 10) */
extern double	Bspline09Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Basic spline function
    of degree 10 (order 11) */
extern double	Bspline10Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Basic spline function
    of degree 11 (order 12) */
extern double	Bspline11Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Dodgson kernel
    (order 2) evaluated at Argument */
extern double	DodgsonIntegrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a cubic Keys optimal kernel
    (order 3) evaluated at Argument */
extern double	KeysOptimalIntegrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Computes the integral(-Infinity, Argument) for a Blu optimum function.
    success: return(!ERROR); failure: return(ERROR); */
extern int		OmomsIntegrate
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Blu interpolant function
    of degree 0 (order 1) */
extern double	Omoms00Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Blu interpolant function
    of degree 1 (order 2) */
extern double	Omoms01Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Blu interpolant function
    of degree 2 (order 3) */
extern double	Omoms02Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the integral(-Infinity, Argument) for a Blu interpolant function
    of degree 3 (order 4) */
extern double	Omoms03Integrate
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the approximate integral(-Infinity, Argument) for a sinc kernel
    evaluated at Argument. Maximum error is about 2.0E-7 */
extern double	SincIntegrate
				(
					double	Argument			/* input */
				);
//@}