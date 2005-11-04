#ifndef _BILIBKERNEL
   #define _BILIBKERNEL
/**@name Kernel definition */
//@{
/*--------------------------------------------------------------------------*/
/** Computes a Blu interpolant function.
    success: return(!ERROR); failure: return(ERROR); */
extern int		Blip
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blu interpolant function of degree 0 (order 1)
    evaluated at Argument. */
extern double	Blip00
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blu interpolant function of degree 1 (order 2)
    evaluated at Argument */
extern double	Blip01
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/* Returns the value of a Blu interpolant function of degree 3 (order 4)
   evaluated at Argument */
extern double	Blip03
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Computes a Basic spline function.
    success: return(!ERROR); failure: return(ERROR); */
extern int		Bspline
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blu interpolant function of degree 0 (order 1)
    evaluated at Argument */
extern double	Bspline00
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Basic spline function of degree 1 (order 2)
   evaluated at Argument */
extern double	Bspline01
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Basic spline function of degree 2 (order 3)
    evaluated at Argument */
extern double	Bspline02
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Basic spline function of degree 3 (order 4)
    evaluated at Argument */
extern double	Bspline03
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Basic spline function of degree 4 (order 5)
   evaluated at Argument */
extern double	Bspline04
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Basic spline function of degree 5 (order 6)
    evaluated at Argument */
extern double	Bspline05
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Basic spline function of degree 6 (order 7)
    evaluated at Argument */
extern double	Bspline06
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Basic spline function of degree 7 (order 8)
   evaluated at Argument */
extern double	Bspline07
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Basic spline function of degree 8 (order 9)
   evaluated at Argument */
extern double	Bspline08
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Basic spline function of degree 9 (order 10)
   evaluated at Argument */
extern double	Bspline09
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Basic spline function of degree 10 (order 11)
    evaluated at Argument */
extern double	Bspline10
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Basic spline function of degree 11 (order 12)
    evaluated at Argument */
extern double	Bspline11
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns 3 values for a Basic spline function of degree 2 (order 3).
   Evaluation is performed at {Argument - 1.0, Argument, Argument + 1.0}.
   Argument must be in [-0.5, 0.5]. Computational load: 3 indirections,
   4(double)assignments, 4 (double)additions, 2 (double)multiplications.
   
   success: return(!ERROR); failure: return(ERROR); */
extern int		BsplineArray02
				(
					double	Argument,			/* fractional input */
					double	*b2_minus1,			/* 1st returned coefficient */
					double	*b2_plus0,			/* 2nd returned coefficient */
					double	*b2_plus1			/* 3rd returned coefficient */
				);

/*--------------------------------------------------------------------------*/
/** Returns 4 values for a Basic spline function of degree 3 (order 4).
    Evaluation is performed at {Argument - 2.0, Argument - 1.0, Argument,
    Argument + 1.0}. Argument must be in [0.0, 1.0].
    Computational load: 7 indirections, 7(double)assignments,
    6 (double)additions, 8 (double)multiplications.
    
    success: return(!ERROR); failure: return(ERROR); */
extern int		BsplineArray03
				(
					double	Argument,			/* fractional input */
					double	*b3_minus2,			/* 1st returned coefficient */
					double	*b3_minus1,			/* 2nd returned coefficient */
					double	*b3_plus0,			/* 3rd returned coefficient */
					double	*b3_plus1			/* 4th returned coefficient */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Dodgson kernel evaluated at Argument (order 2) */
extern double	Dodgson
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns 3 values for a Dodgson kernel (order 2).
    Evaluation is performed at {Argument - 1.0, Argument, Argument + 1.0}.
    Argument must be in [-0.5, 0.5].
    Computational load: 5 indirections, 4(double)assignments,
    4 (double)additions, 3 (double)multiplications
    
    success: return(!ERROR); failure: return(ERROR); */
extern int		DodgsonArray
				(
					double	Argument,			/* fractional input */
					double	*d_minus1,			/* 1st returned coefficient */
					double	*d_plus0,			/* 2nd returned coefficient */
					double	*d_plus1			/* 3rd returned coefficient */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of the quartic German kernel (order 5) evaluated at
   Argument */
extern double	German04
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of the cubic Keys kernel evaluated at Argument */
extern double	Keys
				(
					double	Argument,			/* input */
					double	a					/* tuning parameter */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of the cubic Keys optimal kernel (order 3) evaluated at
   Argument */
extern double	KeysOptimal
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns 4 values for a Keys kernel).
    Evaluation is performed at {Argument - 2.0, Argument - 1.0, Argument,
    Argument + 1.0}. Argument must be in [0.0, 1.0].
    Computational load: 7 indirections, 6(double)assignments,
       6 (double)additions, 8 (double)multiplications.
    
    success: return(!ERROR); failure: return(ERROR); */
extern int		KeysOptimalArray
				(
					double	Argument,			/* fractional input */
					double	*k3_minus2,			/* 1st returned coefficient */
					double	*k3_minus1,			/* 2nd returned coefficient */
					double	*k3_plus0,			/* 3rd returned coefficient */
					double	*k3_plus1			/* 4th returned coefficient */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Meijering function of degree 5 (order ?)
    evaluated at Argument */
extern double	Meijering05
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Meijering function of degree 7 (order ?)
    evaluated at Argument */
extern double	Meijering07
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Computes a Blu optimum function.
    success: return(!ERROR); failure: return(ERROR); */
extern int		Omoms
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blu optimum function of degree 0 (order 1)
    evaluated at Argument */
extern double	Omoms00
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blu optimum function of degree 1 (order 2)
    evaluated at Argument */
extern double	Omoms01
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blu optimum function of degree 2 (order 3)
    evaluated at Argument */
extern double	Omoms02
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blu optimum function of degree 3 (order 4)
    evaluated at Argument */
extern double	Omoms03
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blu optimum function of degree 4 (order 5)
    evaluated at Argument */
extern double	Omoms04
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blu optimum function of degree 5 (order 6)
   evaluated at Argument */
extern double	Omoms05
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blu optimum function of degree 6 (order 7)
    evaluated at Argument */
extern double	Omoms06
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blu optimum function of degree 7 (order 8)
    evaluated at Argument */
extern double	Omoms07
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blu optimum function of degree 8 (order 9)
    evaluated at Argument */
extern double	Omoms08
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blu optimum function of degree 9 (order 10)
    evaluated at Argument */
extern double	Omoms09
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blu optimum function of degree 10 (order 11)
    evaluated at Argument */
extern double	Omoms10
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of a Blu optimum function of degree 11 (order 12)
    evaluated at Argument */
extern double	Omoms11
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns 4 values for an oMoms function of degree 3 (order 4).
    Evaluation is performed at {Argument - 2.0, Argument - 1.0, Argument,
    Argument + 1.0}. Argument must be in [0.0, 1.0]
    Computational load: 7 indirections, 6(double)assignments,
       11 (double)additions, 9 (double)multiplications
   
    success: return(!ERROR); failure: return(ERROR); */
extern int		OmomsArray03
				(
					double	Argument,			/* fractional input */
					double	*b3_minus2,			/* 1st returned coefficient */
					double	*b3_minus1,			/* 2nd returned coefficient */
					double	*b3_plus0,			/* 3rd returned coefficient */
					double	*b3_plus1			/* 4th returned coefficient */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of the Positive kernel of degree 3 evaluated at
    Argument */
extern double	Positive
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Computes a Schaum interpolating function.
    success: return(!ERROR); failure: return(ERROR); */
extern int		Schaum
				(
					long	Degree,				/* degree */
					double	Argument,			/* input */
					double	*Result				/* output */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of the Schaum kernel of degree 2 evaluated at Argument */
extern double	Schaum02
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of the Schaum kernel of degree 3 evaluated at Argument */
extern double	Schaum03
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of the Schaum kernel of degree 4 evaluated at Argument */
extern double	Schaum04
				(
					double	Argument			/* input */
				);

/*--------------------------------------------------------------------------*/
/** Returns the value of the sinc kernel evaluated at Argument */
extern double	Sinc
				(
					double	Argument			/* input */
				);
//@}
#endif
