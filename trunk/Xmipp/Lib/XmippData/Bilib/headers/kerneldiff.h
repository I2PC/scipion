/**@name Kernel derivatives */
//@{
/*--------------------------------------------------------------------------*/
/** Computes a derivative of order Derivative for a Blu interpolant function.
    success: return(!ERROR); failure: return(ERROR); */

extern int		BlipDiff
				(
					long	Degree,				/* degree */
					long	Derivative,			/* order of the derivative */
					double	Argument,			/* input */
					double	*Result				/* output */
				);

/*--------------------------------------------------------------------------*/
/** Computes a derivative of order Derivative for a Basic spline function
   of degree Degree.
   success: return(!ERROR); failure: return(ERROR); */
extern int		BsplineDiff
				(
					long	Degree,				/* degree */
					long	Derivative,			/* order of the derivative */
					double	Argument,			/* input */
					double	*Result				/* output */
				);

/*--------------------------------------------------------------------------*/
/** Computes a derivative of order Derivative for a Dodgson kernel
    (order 2) evaluated at Argument.
    success: return(!ERROR); failure: return(ERROR); */
extern int		DodgsonDiff
				(
					long	Derivative,			/* order of the derivative */
					double	Argument,			/* input */
					double	*Result				/* output */
				);

/*--------------------------------------------------------------------------*/
/** Computes a derivative of order Derivative for a cubic Keys optimal kernel
    (order 3) evaluated at Argument.
    success: return(!ERROR); failure: return(ERROR); */
extern int		KeysOptimalDiff
				(
					long	Derivative,			/* order of the derivative */
					double	Argument,			/* input */
					double	*Result				/* output */
				);

/*--------------------------------------------------------------------------*/
/** Computes a derivative of order Derivative for a Blu optimum function.
    success: return(!ERROR); failure: return(ERROR); */
extern int		OmomsDiff
				(
					long	Degree,				/* degree */
					long	Derivative,			/* order of the derivative */
					double	Argument,			/* input */
					double	*Result				/* output */
				);

/*--------------------------------------------------------------------------*/
/** Computes a derivative of order Derivative for a Positive kernel
    (order 1) evaluated at Argument.
    success: return(!ERROR); failure: return(ERROR); */
extern int		PositiveDiff
				(
					long	Derivative,			/* order of the derivative */
					double	Argument,			/* input */
					double	*Result				/* output */
				);

/*--------------------------------------------------------------------------*/
/** Computes a derivative of order Derivative for a sinc kernel evaluated
    at Argument. success: return(!ERROR); failure: return(ERROR); */
extern int		SincDiff
				(
					long	Derivative,			/* order of the derivative */
					double	Argument,			/* input */
					double	*Result				/* output */
				);
//@}

