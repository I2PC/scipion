/**@name Polynomial */
//@{
/*--------------------------------------------------------------------------*/
/** Differentiation of the polynomial p(x) = (a[0] + Sum a[k] x^k).
    The degree of the input polynomial is Degree.
    There are (Degree+1) input coefficients.
    The degree of the output polynomial is (Degree-1).
    There are Degree output coefficients.
    
    success: return(!ERROR); failure: return(ERROR) */
extern int		PolynomialDifferentiation
				(
					double	a[],				/* input polynomial coefficients */
					long	Degree,				/* degree of the input polynomial */
					double	b[]					/* resulting coefficients */
				);

/*--------------------------------------------------------------------------*/
/** Evaluates a polynomial: result = a[0] + Sum a[k] x^k.
    The degree of the polynomial is Degree.
    There are (Degree+1) coefficients.
    
    success: return(!ERROR); failure: return(ERROR) */
extern int		PolynomialEvaluation
				(
					double	x,					/* argument */
					double	a[],				/* polynomial coefficients */
					long	Degree,				/* degree of the polynomial */
					double	*Result				/* resulting value */
				);

/*--------------------------------------------------------------------------*/
/** Polynomial multuplication.
    Multiplication resulting in
    b[0] + Sum b[k] x^k = (a1[0] + Sum a1[k] x^k) (a2[0] + Sum a2[k] x^k).
    
    The degree of the 1st input polynomial is Degree1.
    The degree of the 2nd input polynomial is Degree2.
    There are (Degree1+1) input coefficients for the first multiplicand.
    There are (Degree2+1) input coefficients for the second multiplicand.
    The degree of the output polynomial is (Degree1+Degree2).
    There are (Degree1+Degree2+1) output coefficients.
    
    success: return(!ERROR); failure: return(ERROR) */
extern int		PolynomialMultiplication
				(
					double	a1[],				/* 1st input polynomial coefficients */
					long	Degree1,			/* degree of the 1st input polynomial */
					double	a2[],				/* 2nd input polynomial coefficients */
					long	Degree2,			/* degree of the 2nd input polynomial */
					double	b[]					/* resulting coefficients */
				);

/*--------------------------------------------------------------------------*/
/** Primitive of the polynomial p(t) = (a[0] + Sum a[k] t^k).
    The degree of the input polynomial is Degree.
    There are (Degree+1) input coefficients.
    The degree of the output polynomial is (Degree+1).
    There are (Degree+2) output coefficients.
    
    success: return(!ERROR); failure: return(ERROR) */
extern int		PolynomialPrimitive
				(
					double	a[],				/* input polynomial coefficients */
					long	Degree,				/* degree of the input polynomial */
					double	b[]					/* resulting coefficients */
				);

/*--------------------------------------------------------------------------*/
/** Find the real roots of the polynomial p(x) = (a[0] + Sum a[k] x^k).
    The degree of the polynomial is Degree.
    There are (Degree+1) input coefficients (a[]).
    The output array (RealRoot[]) must have size (Degree).
    Only the first (RealRootNumber) roots returned in RealRoot[] are valid.
    (RealRootNumber) -> -1 when the equation is indeterminate.
    The returned roots are sorted in ascendent order.
    
    success: return(!ERROR); failure: return(ERROR) */
extern int		PolynomialRealRoots
				(
					double	a[],				/* polynomial coefficients */
					long	Degree,				/* degree of the polynomial */
					double	RealRoot[],			/* resulting real roots */
					long	*RealRootNumber,	/* number of real roots */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				);

/*--------------------------------------------------------------------------*/
/** Returns the sign of x with Sign(0) = 0. */
extern double	Sign
				(
					double	x					/* argument */
				);

/*--------------------------------------------------------------------------*/
/** Computes the one-sided power function. xPlus(0, 0) = 1/2 */
extern double	xPlus
				(
					double	x,					/* argument */
					double	p					/* power */
				);
//@}