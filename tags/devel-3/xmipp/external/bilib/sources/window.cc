/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<math.h>
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
#include	"messagedisplay.h"
#include	"window.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
#ifdef		DEBUG
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
extern double	Bartlet
				(
					double	Argument,			/* input */
					long	WindowHalfLength	/* half size H of the window */
				)

/* returns the value of a Bartlet (triangular) window evaluated at Argument */
/* the width of the symmetric window is (2 * H) */

{ /* begin Bartlet */

#ifdef DEBUG
/**/int		Status = !ERROR;
#endif

/**/DEBUG_CHECK_RANGE_LONG(Bartlet, WindowHalfLength, 1L, LONG_MAX, Status,
/**/	"Invalid window half length (should be strictly positive)")
/**/DEBUG_WRITE_ENTERING(Bartlet,
/**/	"About to compute a Bartlet value")

#ifdef DEBUG
/**/if (Status == ERROR) {
/**/	WRITE_ERROR(Bartlet,
/**/		"Unable to perform Bartlet because of invalid window length")
/**/	DEBUG_WRITE_LEAVING(Bartlet, "Done")
/**/	return(HUGE_VAL);
/**/}
#endif
	Argument = fabs(Argument);
	if (Argument < WindowHalfLength) {
/**/	DEBUG_WRITE_LEAVING(Bartlet, "Done")
		return(1.0 - fabs(Argument) / (double)WindowHalfLength);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Bartlet, "Done")
		return(0.0);
	}
} /* end Bartlet */

/*--------------------------------------------------------------------------*/
extern double	Blackman
				(
					double	Argument,			/* input */
					long	WindowHalfLength	/* half size H of the window */
				)

/* returns the value of a Blackman window evaluated at Argument */
/* the width of the symmetric window is (2 * H) */

{ /* begin Blackman */

#ifdef DEBUG
/**/int		Status = !ERROR;
#endif

/**/DEBUG_CHECK_RANGE_LONG(Blackman, WindowHalfLength, 1L, LONG_MAX, Status,
/**/	"Invalid window half length (should be strictly positive)")
/**/DEBUG_WRITE_ENTERING(Blackman,
/**/	"About to compute a Blackman value")

#ifdef DEBUG
/**/if (Status == ERROR) {
/**/	WRITE_ERROR(Blackman,
/**/		"Unable to perform Blackman because of invalid window length")
/**/	DEBUG_WRITE_LEAVING(Blackman, "Done")
/**/	return(HUGE_VAL);
/**/}
#endif
	Argument = fabs(Argument);
	if (Argument < WindowHalfLength) {
/**/	DEBUG_WRITE_LEAVING(Blackman, "Done")
		return(3969.0 / 9304.0 + (1155.0 / 2326.0)
			* cos(Argument * PI / (double)WindowHalfLength) + (715.0 / 9304.0)
			* cos(Argument * (2.0 * PI) / (double)WindowHalfLength));
	}
	else if (Argument == WindowHalfLength) {
/**/	DEBUG_WRITE_LEAVING(Blackman, "Done")
		return(4.0 / 1163.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Blackman, "Done")
		return(0.0);
	}
} /* end Blackman */

/*--------------------------------------------------------------------------*/
extern double	Dirichlet
				(
					double	Argument,			/* input */
					long	WindowHalfLength	/* half size H of the window */
				)

/* returns the value of a Dirichlet (rectangular) window evaluated at Argument */
/* the width of the symmetric window is (2 * H) */

{ /* begin Dirichlet */

#ifdef DEBUG
/**/int		Status = !ERROR;
#endif

/**/DEBUG_CHECK_RANGE_LONG(Dirichlet, WindowHalfLength, 1L, LONG_MAX, Status,
/**/	"Invalid window half length (should be strictly positive)")
/**/DEBUG_WRITE_ENTERING(Dirichlet,
/**/	"About to compute a Dirichlet value")

#ifdef DEBUG
/**/if (Status == ERROR) {
/**/	WRITE_ERROR(Dirichlet,
/**/		"Unable to perform Dirichlet because of invalid window length")
/**/	DEBUG_WRITE_LEAVING(Dirichlet, "Done")
/**/	return(HUGE_VAL);
/**/}
#endif
	Argument = fabs(Argument);
	if (Argument < WindowHalfLength) {
/**/	DEBUG_WRITE_LEAVING(Dirichlet, "Done")
		return(1.0);
	}
	else if (Argument == WindowHalfLength) {
/**/	DEBUG_WRITE_LEAVING(Dirichlet, "Done")
		return(0.5);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Dirichlet, "Done")
		return(0.0);
	}
} /* end Dirichlet */

/*--------------------------------------------------------------------------*/
extern double	HammingClassic
				(
					double	Argument,			/* input */
					long	WindowHalfLength	/* half size H of the window */
				)

/* returns the value of a Hamming window evaluated at Argument */
/* classic Hamming weights are used */
/* the width of the symmetric window is (2 * H) */

{ /* begin HammingClassic */

#ifdef DEBUG
/**/int		Status = !ERROR;
#endif

/**/DEBUG_CHECK_RANGE_LONG(HammingClassic, WindowHalfLength, 1L, LONG_MAX, Status,
/**/	"Invalid window half length (should be strictly positive)")
/**/DEBUG_WRITE_ENTERING(HammingClassic,
/**/	"About to compute a Hamming value")

#ifdef DEBUG
/**/if (Status == ERROR) {
/**/	WRITE_ERROR(HammingClassic,
/**/		"Unable to perform HammingClassic because of invalid window length")
/**/	DEBUG_WRITE_LEAVING(HammingClassic, "Done")
/**/	return(HUGE_VAL);
/**/}
#endif
	if (Argument < WindowHalfLength) {
/**/	DEBUG_WRITE_LEAVING(HammingClassic, "Done")
		return(0.54 + 0.46 * cos(Argument * PI / (double)WindowHalfLength));
	}
	else if (Argument == WindowHalfLength) {
/**/	DEBUG_WRITE_LEAVING(HammingClassic, "Done")
		return(0.04);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(HammingClassic, "Done")
		return(0.0);
	}
} /* end HammingClassic */

/*--------------------------------------------------------------------------*/
extern double	HammingExact
				(
					double	Argument,			/* input */
					long	WindowHalfLength	/* half size H of the window */
				)

/* returns the value of a Hamming window evaluated at Argument */
/* optimal (non-classic) weights are used */
/* the width of the symmetric window is (2 * H) */

{ /* begin HammingExact */

#ifdef DEBUG
/**/int		Status = !ERROR;
#endif

/**/DEBUG_CHECK_RANGE_LONG(HammingExact, WindowHalfLength, 1L, LONG_MAX, Status,
/**/	"Invalid window half length (should be strictly positive)")
/**/DEBUG_WRITE_ENTERING(HammingExact,
/**/	"About to compute an exact Hamming value")

#ifdef DEBUG
/**/if (Status == ERROR) {
/**/	WRITE_ERROR(HammingExact,
/**/		"Unable to perform HammingExact because of invalid window length")
/**/	DEBUG_WRITE_LEAVING(HammingExact, "Done")
/**/	return(HUGE_VAL);
/**/}
#endif
	if (Argument < WindowHalfLength) {
/**/	DEBUG_WRITE_LEAVING(HammingExact, "Done")
		return(25.0 / 46.0 + (21.0 / 46.0) * cos(Argument * PI / (double)WindowHalfLength));
	}
	else if (Argument == WindowHalfLength) {
/**/	DEBUG_WRITE_LEAVING(HammingExact, "Done")
		return(1.0 / 23.0);
	}
	else {
/**/	DEBUG_WRITE_LEAVING(HammingExact, "Done")
		return(0.0);
	}
} /* end HammingExact */

/*--------------------------------------------------------------------------*/
extern double	Hanning
				(
					double	Argument,			/* input */
					long	WindowHalfLength	/* half size H of the window */
				)

/* returns the value of a Hanning window evaluated at Argument */
/* the width of the symmetric window is (2 * H) */

{ /* begin Hanning */

#ifdef DEBUG
/**/int		Status = !ERROR;
#endif

/**/DEBUG_CHECK_RANGE_LONG(Hanning, WindowHalfLength, 1L, LONG_MAX, Status,
/**/	"Invalid window half length (should be strictly positive)")
/**/DEBUG_WRITE_ENTERING(Hanning,
/**/	"About to compute a Hanning value")

#ifdef DEBUG
/**/if (Status == ERROR) {
/**/	WRITE_ERROR(Hanning,
/**/		"Unable to perform Hanning because of invalid window length")
/**/	DEBUG_WRITE_LEAVING(Hanning, "Done")
/**/	return(HUGE_VAL);
/**/}
#endif
	if (Argument < WindowHalfLength) {
/**/	DEBUG_WRITE_LEAVING(Hanning, "Done")
		return(1.0 / 2.0 + (1.0 / 2.0) * cos(Argument * PI / (double)WindowHalfLength));
	}
	else {
/**/	DEBUG_WRITE_LEAVING(Hanning, "Done")
		return(0.0);
	}
} /* end Hanning */

