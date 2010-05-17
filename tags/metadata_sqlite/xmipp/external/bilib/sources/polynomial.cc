/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>

/*****************************************************************************
 *	Toolbox defines
 ****************************************************************************/
#include	"configs.h"
#include	"debug.h"
#include	"error.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"linearalgebra.h"
#include	"messagedisplay.h"
#include	"polynomial.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
#ifdef		DEBUG
#include	<float.h>
#include	<limits.h>
#include	<string.h>
#endif

/*****************************************************************************
 *	Declaration of static procedures
 ****************************************************************************/
/* None */

/*****************************************************************************
 *	Definition of static procedures
 ****************************************************************************/
/* None */

/*****************************************************************************
 *	Definition of extern procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
extern int		PolynomialDifferentiation
				(
					double	a[],				/* input polynomial coefficients */
					long	Degree,				/* degree of the input polynomial */
					double	b[]					/* resulting coefficients */
				)

/* differentiation of the polynomial p(x) = (a[0] + Sum a[k] x^k) */
/* the degree of the input polynomial is Degree */
/* there are (Degree+1) input coefficients */
/* the degree of the output polynomial is (Degree-1) */
/* there are Degree output coefficients */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin PolynomialDifferentiation */

	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(PolynomialDifferentiation, a, Status,
/**/	"Missing input polynomial coefficients")
/**/DEBUG_CHECK_NULL_POINTER(PolynomialDifferentiation, b, Status,
/**/	"Missing output polynomial coefficients")
/**/DEBUG_CHECK_RANGE_LONG(PolynomialDifferentiation, Degree, 0L, LONG_MAX, Status,
/**/	"Invalid polynomial degree (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(PolynomialDifferentiation, Status)
/**/DEBUG_WRITE_ENTERING(PolynomialDifferentiation,
/**/	"About to execute PolynomialDifferentiation")

	if (Degree == 0L) {
		*b = 0.0;
/**/	DEBUG_WRITE_LEAVING(PolynomialDifferentiation, "Done")
		return(Status);
	}
	for (i = 1L; (i <= Degree); i++) {
		*b++ = *++a * (double)i;
	}
/**/DEBUG_WRITE_LEAVING(PolynomialDifferentiation, "Done")
	return(Status);
} /* end PolynomialDifferentiation */

/*--------------------------------------------------------------------------*/
extern int		PolynomialEvaluation
				(
					double	x,					/* argument */
					double	a[],				/* polynomial coefficients */
					long	Degree,				/* degree of the polynomial */
					double	*Result				/* resulting value */
				)

/* evaluates a polynomial: result = a[0] + Sum a[k] x^k */
/* the degree of the polynomial is Degree */
/* there are (Degree+1) coefficients */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin PolynomialEvaluation */

	double	*p;
	double	Sum;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(PolynomialEvaluation, a, Status,
/**/	"Missing polynomial coefficients")
/**/DEBUG_CHECK_NULL_POINTER(PolynomialEvaluation, Result, Status,
/**/	"Missing Result")
/**/DEBUG_CHECK_RANGE_LONG(PolynomialEvaluation, Degree, 0L, LONG_MAX, Status,
/**/	"Invalid polynomial degree (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(PolynomialEvaluation, Status)
/**/DEBUG_WRITE_ENTERING(PolynomialEvaluation,
/**/	"About to evaluate a polynomial")

	p = a + (ptrdiff_t)Degree;
	Sum = *p--;
	while (a <= p) {
		Sum = x * Sum + *p--;
	}
	*Result = Sum;
/**/DEBUG_WRITE_LEAVING(PolynomialEvaluation, "Done")
	return(Status);
} /* end PolynomialEvaluation */

/*--------------------------------------------------------------------------*/
extern int		PolynomialMultiplication
				(
					double	a1[],				/* 1st input polynomial coefficients */
					long	Degree1,			/* degree of the 1st input polynomial */
					double	a2[],				/* 2nd input polynomial coefficients */
					long	Degree2,			/* degree of the 2nd input polynomial */
					double	b[]					/* resulting coefficients */
				)

/* multiplication resulting in */
/*   b[0] + Sum b[k] x^k = (a1[0] + Sum a1[k] x^k) (a2[0] + Sum a2[k] x^k) */
/* the degree of the 1st input polynomial is Degree1 */
/* the degree of the 2nd input polynomial is Degree2 */
/* there are (Degree1+1) input coefficients for the first multiplicand */
/* there are (Degree2+1) input coefficients for the second multiplicand */
/* the degree of the output polynomial is (Degree1+Degree2) */
/* there are (Degree1+Degree2+1) output coefficients */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin PolynomialMultiplication */

	double	*p, *q;
	double	Sum;
	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(PolynomialMultiplication, a1, Status,
/**/	"Missing first input polynomial coefficients")
/**/DEBUG_CHECK_NULL_POINTER(PolynomialMultiplication, a2, Status,
/**/	"Missing second input polynomial coefficients")
/**/DEBUG_CHECK_NULL_POINTER(PolynomialMultiplication, b, Status,
/**/	"Missing output polynomial coefficients")
/**/DEBUG_CHECK_RANGE_LONG(PolynomialMultiplication, Degree1, 0L, LONG_MAX, Status,
/**/	"Invalid first polynomial degree (should be positive)")
/**/DEBUG_CHECK_RANGE_LONG(PolynomialMultiplication, Degree2, 0L, LONG_MAX, Status,
/**/	"Invalid second polynomial degree (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(PolynomialMultiplication, Status)
/**/DEBUG_WRITE_ENTERING(PolynomialMultiplication,
/**/	"About to execute PolynomialMultiplication")

	for (i = 0L; (i <= (Degree1 + Degree2)); i++) {
		Sum = 0.0;
		p = a2;
		q = a1 + (ptrdiff_t)i;
		while (a1 <= q) {
			Sum += *q-- * *p++;
		}
		*b++ = Sum;
	}
/**/DEBUG_WRITE_LEAVING(PolynomialMultiplication, "Done")
	return(Status);
} /* end PolynomialMultiplication */

/*--------------------------------------------------------------------------*/
extern int		PolynomialPrimitive
				(
					double	a[],				/* input polynomial coefficients */
					long	Degree,				/* degree of the input polynomial */
					double	b[]					/* resulting coefficients */
				)

/* primitive of the polynomial p(t) = (a[0] + Sum a[k] t^k) */
/* the degree of the input polynomial is Degree */
/* there are (Degree+1) input coefficients */
/* the degree of the output polynomial is (Degree+1) */
/* there are (Degree+2) output coefficients */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin PolynomialPrimitive */

	long	i;
	int		Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(PolynomialPrimitive, a, Status,
/**/	"Missing input polynomial coefficients")
/**/DEBUG_CHECK_NULL_POINTER(PolynomialPrimitive, b, Status,
/**/	"Missing output polynomial coefficients")
/**/DEBUG_CHECK_RANGE_LONG(PolynomialPrimitive, Degree, 0L, LONG_MAX, Status,
/**/	"Invalid polynomial degree (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(PolynomialPrimitive, Status)
/**/DEBUG_WRITE_ENTERING(PolynomialPrimitive,
/**/	"About to execute PolynomialPrimitive")

	*b++ = 0.0;
	for (i = 0L; (i <= Degree); i++) {
		*b++ = *a++ / (double)(i + 1L);
	}
/**/DEBUG_WRITE_LEAVING(PolynomialPrimitive, "Done")
	return(Status);
} /* end PolynomialPrimitive */

/*--------------------------------------------------------------------------*/
extern int		PolynomialRealRoots
				(
					double	a[],				/* polynomial coefficients */
					long	Degree,				/* degree of the polynomial */
					double	RealRoot[],			/* resulting real roots */
					long	*RealRootNumber,	/* number of real roots */
					double	Tolerance,			/* admissible relative error */
					int		*Status				/* error management */
				)

/* find the real roots of the polynomial p(x) = (a[0] + Sum a[k] x^k) */
/* the degree of the polynomial is Degree */
/* there are (Degree+1) input coefficients (a[]) */
/* the output array (RealRoot[]) must have size (Degree) */
/* only the first (RealRootNumber) roots returned in RealRoot[] are valid */
/* (RealRootNumber) -> -1 when the equation is indeterminate */
/* the returned roots are sorted in ascendent order */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin PolynomialRealRoots */

	double	CubicResolvent[4], QuadraticResolvent[3];
	double	p, q, r, s, t;
	double	Det, Swap;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(PolynomialRealRoots, a, *Status,
/**/	"Missing polynomial coefficients")
/**/DEBUG_CHECK_NULL_POINTER(PolynomialRealRoots, RealRoot, *Status,
/**/	"Missing RealRoot")
/**/DEBUG_CHECK_NULL_POINTER(PolynomialRealRoots, RealRootNumber, *Status,
/**/	"Missing RealRootNumber")
/**/DEBUG_CHECK_RANGE_LONG(PolynomialRealRoots, Degree, 0L, LONG_MAX, *Status,
/**/	"Invalid polynomial degree (should be positive)")
/**/DEBUG_CHECK_RANGE_DOUBLE(PolynomialRealRoots, Tolerance, 0.0, DBL_MAX, *Status,
/**/	"Invalid Tolerance (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(PolynomialRealRoots, *Status)
/**/DEBUG_WRITE_ENTERING(PolynomialRealRoots,
/**/	"About to execute PolynomialRealRoots")

	switch (Degree) {
		case 0L:
			if (fabs(a[Degree]) <= Tolerance) {
				*RealRootNumber = -1L;
				WRITE_WARNING(PolynomialRealRoots, "Indeterminate equation")
				break;
			}
			*Status = ERROR;
			WRITE_ERROR(PolynomialRealRoots, "Absurd equation encountered (0 != 0)")
			break;
		case 1L:
			if (fabs(a[Degree]) <= Tolerance) {
				PolynomialRealRoots(a, Degree - 1L, RealRoot, RealRootNumber,
					Tolerance, Status);
				break;
			}
			RealRoot[0] = -a[0] / a[1];
			*RealRootNumber = 1L;
			break;
		case 2L:
			if (fabs(a[Degree]) <= Tolerance) {
				PolynomialRealRoots(a, Degree - 1L, RealRoot, RealRootNumber,
					Tolerance, Status);
				break;
			}
			Det = a[1] * a[1] - 4.0 * a[0] * a[2];
			if (Det < 0.0) {
				*RealRootNumber = 0L;
			}
			else if (Det == 0.0) {
				RealRoot[0] = -0.5 * a[1] / a[2];
				RealRoot[1] = RealRoot[0];
				*RealRootNumber = 2L;
			}
			else {
				Det = (a[1] < 0.0) ? (a[1] - sqrt(Det)) : (a[1] + sqrt(Det));
				RealRoot[0] = -0.5 * Det / a[2];
				RealRoot[1] = -2.0 * a[0] / Det;
				if (RealRoot[1] < RealRoot[0]) {
					Swap = RealRoot[0];
					RealRoot[0] = RealRoot[1];
					RealRoot[1] = Swap;
				}
				*RealRootNumber = 2L;
			}
			break;
		case 3L:
			if (fabs(a[Degree]) <= Tolerance) {
				PolynomialRealRoots(a, Degree - 1L, RealRoot, RealRootNumber,
					Tolerance, Status);
				break;
			}
			q = (3.0 * a[1] - a[2] * a[2] / a[3]) / (9.0 * a[3]);
			r = (a[1] * a[2] - 3.0 * a[0] * a[3] - (2.0 / 9.0) * a[2] * a[2] * a[2] / a[3])
				/ (6.0 * a[3] * a[3]);
			Det = q * q * q + r * r;
			if (Det < 0.0) {
				Det = 2.0 * sqrt(-q);
				r /= sqrt(-q * q * q);
				r = (1.0 / 3.0) * acos(r);
				q = (1.0 / 3.0) * a[2] / a[3];
				RealRoot[0] = Det * cos(r) - q;
				RealRoot[1] = Det * cos(r + (2.0 / 3.0) * PI) - q;
				RealRoot[2] = Det * cos(r + (4.0 / 3.0) * PI) - q;
				if (RealRoot[0] < RealRoot[1]) {
					if (RealRoot[2] < RealRoot[1]) {
						Swap = RealRoot[1];
						RealRoot[1] = RealRoot[2];
						RealRoot[2] = Swap;
						if (RealRoot[1] < RealRoot[0]) {
							Swap = RealRoot[0];
							RealRoot[0] = RealRoot[1];
							RealRoot[1] = Swap;
						}
					}
				}
				else {
					Swap = RealRoot[0];
					RealRoot[0] = RealRoot[1];
					RealRoot[1] = Swap;
					if (RealRoot[2] < RealRoot[1]) {
						Swap = RealRoot[1];
						RealRoot[1] = RealRoot[2];
						RealRoot[2] = Swap;
						if (RealRoot[1] < RealRoot[0]) {
							Swap = RealRoot[0];
							RealRoot[0] = RealRoot[1];
							RealRoot[1] = Swap;
						}
					}
				}
				*RealRootNumber = 3L;
			}
			else if (Det == 0.0) {
				p = 2.0 * pow(r, 1.0 / 3.0);
				q = (1.0 / 3.0) * a[2] / a[3];
				if (p < 0) {
					RealRoot[0] = p - q;
					RealRoot[1] = -0.5 * p - q;
					RealRoot[2] = RealRoot[1];
				}
				else {
					RealRoot[0] = -0.5 * p - q;
					RealRoot[1] = RealRoot[0];
					RealRoot[2] = p - q;
				}
				*RealRootNumber = 3L;
			}
			else {
				Det = sqrt(Det);
				q = ((r + Det) < 0.0) ? (-exp((1.0 / 3.0) * log(-r - Det)))
					: (exp((1.0 / 3.0) * log(r + Det)));
				r = (r < Det) ? (-exp((1.0 / 3.0) * log(Det - r)))
					: (exp((1.0 / 3.0) * log(r - Det)));
				Det = (-1.0 / 3.0) * a[2] / a[3];
				RealRoot[0] = q + r + Det;
				*RealRootNumber = 1L;
			}
			break;
		case 4L:
			if (fabs(a[Degree]) <= Tolerance) {
				PolynomialRealRoots(a, Degree - 1L, RealRoot, RealRootNumber,
					Tolerance, Status);
				break;
			}
			s = a[3] / a[4];
			p = a[2] / a[4] - (3.0 / 8.0) * s * s;
			s = 2.0 * a[4];
			t = a[3] * a[3];
			q = a[1] / a[4] + a[3] * (t - 4.0 * a[4] * a[2]) / (s * s * s);
			r = a[3] * (3.0 * t * a[3] + 16.0 * (s * s * a[1] - a[2] * a[3] * a[4]));
			s *= 4.0 * s;
			r = a[0] / a[4] - r / (s * s);
			s = 0.25 * a[3] / a[4];
			CubicResolvent[0] = 4.0 * p * r - q * q;
			CubicResolvent[1] = -4.0 * r;
			CubicResolvent[2] = -p;
			CubicResolvent[3] = 1.0;
			PolynomialRealRoots(CubicResolvent, 3L, RealRoot, RealRootNumber,
				Tolerance, Status);
			if (*Status == ERROR) {
				WRITE_ERROR(PolynomialRealRoots,
					"Unable to solve the cubic resolvent of a quartic equation")
				break;
			}
			switch (*RealRootNumber) {
				case 1L:
					t = RealRoot[0];
					break;
				case 3L:
					t = RealRoot[2];
					break;
				default:
					*Status = ERROR;
					WRITE_ERROR(PolynomialRealRoots,
						"Unexpected error (wrong number of cubic roots)")
/**/				DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
					return(*Status);
			}
			if (t < p) {
				*RealRootNumber = 0L;
/**/			DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
				return(*Status);
			}
			if ((t == p) || (fabs(q) <= Tolerance)) {
				QuadraticResolvent[0] = r;
				QuadraticResolvent[1] = p;
				QuadraticResolvent[2] = 1.0;
				PolynomialRealRoots(QuadraticResolvent, 2L, RealRoot, RealRootNumber,
					Tolerance, Status);
				if (*Status == ERROR) {
					WRITE_ERROR(PolynomialRealRoots,
						"Unable to solve the quadratic resolvent of a degenerate quartic equation")
					break;
				}
				switch (*RealRootNumber) {
					case 0L:
/**/					DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
						return(*Status);
					case 2L:
						if ((RealRoot[0] < 0.0) && (RealRoot[1] < 0.0)) {
							*RealRootNumber = 0L;
/**/						DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
							return(*Status);
						}
						if ((RealRoot[0] < 0.0) && (0.0 <= RealRoot[1])) {
							RealRoot[0] = -sqrt(RealRoot[1]) - s;
							RealRoot[1] = sqrt(RealRoot[1]) - s;
/**/						DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
							return(*Status);
						}
						*RealRootNumber = 4L;
						RealRoot[3] = sqrt(RealRoot[1]) - s;
						RealRoot[2] = sqrt(RealRoot[0]) - s;
						RealRoot[1] = -sqrt(RealRoot[0]) - s;
						RealRoot[0] = -RealRoot[3] - 2.0 * s;
/**/					DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
						return(*Status);
					default:
						*Status = ERROR;
						WRITE_ERROR(PolynomialRealRoots,
							"Unexpected error (wrong number of degenerate quadratic roots)")
/**/					DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
						return(*Status);
				}
			}
			r = sqrt(t - p);
			QuadraticResolvent[0] = 0.5 * (t - q / r);
			QuadraticResolvent[1] = r;
			QuadraticResolvent[2] = 1.0;
			PolynomialRealRoots(QuadraticResolvent, 2L, RealRoot, RealRootNumber,
				Tolerance, Status);
			if (*Status == ERROR) {
				WRITE_ERROR(PolynomialRealRoots,
					"Unable to solve the first quadratic resolvent of a quartic equation")
/**/			DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
				return(*Status);
			}
			switch (*RealRootNumber) {
				case 0L:
					QuadraticResolvent[0] = 0.5 * (t + q / r);
					QuadraticResolvent[1] = -r;
					PolynomialRealRoots(QuadraticResolvent, 2L, RealRoot, RealRootNumber,
						Tolerance, Status);
					if (*Status == ERROR) {
						WRITE_ERROR(PolynomialRealRoots,
							"Unable to solve the second quadratic resolvent of a quartic equation")
						break;
					}
					switch (*RealRootNumber) {
						case 0L:
/**/						DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
							return(*Status);
						case 2L:
							RealRoot[0] -= s;
							RealRoot[1] -= s;
/**/						DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
							return(*Status);
						default:
							*Status = ERROR;
							WRITE_ERROR(PolynomialRealRoots,
								"Unexpected error (wrong number of quadratic roots)")
							break;
					}
					break;
				case 2L:
					QuadraticResolvent[0] = 0.5 * (t + q / r);
					QuadraticResolvent[1] = -r;
					PolynomialRealRoots(QuadraticResolvent, 2L, RealRoot + (ptrdiff_t)2,
						RealRootNumber, Tolerance, Status);
					if (*Status == ERROR) {
						WRITE_ERROR(PolynomialRealRoots,
							"Unable to solve the second quadratic resolvent of a quartic equation")
						break;
					}
					switch (*RealRootNumber) {
						case 0L:
							RealRoot[0] -= s;
							RealRoot[1] -= s;
/**/						DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
							return(*Status);
						case 2L:
							*RealRootNumber = 4L;
							RealRoot[0] -= s;
							RealRoot[1] -= s;
							RealRoot[2] -= s;
							RealRoot[3] -= s;
							if (RealRoot[2] < RealRoot[0]) {
								r = RealRoot[0];
								RealRoot[0] = RealRoot[2];
								RealRoot[2] = r;
								r = RealRoot[1];
								RealRoot[1] = RealRoot[3];
								RealRoot[3] = r;
							}
							if (RealRoot[1] < RealRoot[2]) {
/**/							DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
								return(*Status);
							}
							else {
								r = RealRoot[1];
								RealRoot[1] = RealRoot[2];
								RealRoot[2] = r;
								if (RealRoot[2] < RealRoot[3]) {
/**/								DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
									return(*Status);
								}
								else {
									r = RealRoot[2];
									RealRoot[2] = RealRoot[3];
									RealRoot[3] = r;
/**/								DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
									return(*Status);
								}
							}
/**/						DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
							return(*Status);
						default:
							*Status = ERROR;
							WRITE_ERROR(PolynomialRealRoots,
								"Unexpected error (wrong number of quadratic roots)")
							break;
					}
					break;
				default:
					*Status = ERROR;
					WRITE_ERROR(PolynomialRealRoots,
						"Unexpected error (wrong number of quadratic roots)")
					break;
			}
			break;
		default:
			if (a[Degree] == 0.0) {
				PolynomialRealRoots(a, Degree - 1L, RealRoot, RealRootNumber,
					Tolerance, Status);
				break;
			}
			*Status = ERROR;
			WRITE_ERROR(PolynomialRealRoots, "Not yet implemented")
			break;
	}
/**/DEBUG_WRITE_LEAVING(PolynomialRealRoots, "Done")
	return(*Status);
} /* end PolynomialRealRoots */

/*--------------------------------------------------------------------------*/
extern double	Sign
				(
					double	x					/* argument */
				)

/* returns the sign of x with Sign(0) = 0 */

{ /* begin Sign */

/**/DEBUG_WRITE_ENTERING(Sign,
/**/	"About to execute Sign")

	if (x < 0.0) {
		x = -1.0;
	}
	else if (0.0 < x) {
		x = 1.0;
	}
/**/DEBUG_WRITE_LEAVING(Sign, "Done")
	return(x);
} /* end Sign */

/*--------------------------------------------------------------------------*/
extern double	xPlus
				(
					double	x,					/* argument */
					double	p					/* power */
				)

/* computes the one-sided power function with xPlus(0, 0) = 1/2 */

{ /* begin xPlus */

/**/DEBUG_WRITE_ENTERING(xPlus,
/**/	"About to execute xPlus")

	x = (x < 0.0) ? (0.0) : ((p == 0.0) ? ((x == 0.0) ? (0.5) : (1.0)) : (pow(x, p)));
/**/DEBUG_WRITE_LEAVING(xPlus, "Done")
	return(x);
} /* end xPlus */

