/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>

/*****************************************************************************
 *	Toolbox defines
 ****************************************************************************/
#include	"configs.h"
#include	"debug.h"
#include	"error.h"
#include	"tboundaryconvention.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"../headers/convert.h"
#include	"fold.h"
#include	"getput.h"
#include	"interpolate.h"
#include	"kernel.h"
#include	"messagedisplay.h"

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
extern int		SplineInterpolateImage
				(
					float	*ImageCoeff,		/* B-spline coefficients to interpolate */
					long	Nx,					/* width of the image */
					long	Ny,					/* height of the image */
					double	Xargument,			/* input X abscissa */
					double	Yargument,			/* input Y abscissa */
					double	*Result,			/* output ordinate */
					long	Degree,				/* degree of the spline */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

{ /* begin SplineInterpolateImage */

	double	*wx = (double *)NULL, *wy = (double *)NULL;
	double	*ux, *uy;
	float	*pi;
	long	*v, *vx, *vy;
	long	*kx, *ky;
	double	Xfrac, Yfrac;
	double	w0, w1;
	double	qi, qj;
	long	i, j, k, l;
	long	x, y;
	long	x0, x1, y0, y1;
	const long
			HalfDegree = Degree / 2L, Width = Degree + 1L;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(SplineInterpolateImage, ImageCoeff, *Status,
/**/	"Missing ImageCoeff")
/**/DEBUG_CHECK_NULL_POINTER(SplineInterpolateImage, Result, *Status,
/**/	"Missing Result")
/**/DEBUG_CHECK_RANGE_LONG(SplineInterpolateImage, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid Nx (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SplineInterpolateImage, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid Ny (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SplineInterpolateImage, Degree, 0L, LONG_MAX, *Status,
/**/	"Invalid Degree (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(SplineInterpolateImage, *Status)
/**/DEBUG_WRITE_ENTERING(SplineInterpolateImage,
/**/	"About to execute SplineInterpolateImage")

	switch (Convention) {
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case FiniteDataSupport:
			if (2L <= Degree) {
				*Status = ERROR;
				WRITE_ERROR(SplineInterpolateImage, "Invalid boundary convention")
/**/			DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
				return(*Status);
			}
			break;
		case FiniteCoefficientSupport:
			Convention = FiniteDataSupport;
			break;
		case AntiMirrorOnBounds:
		default:
			*Status = ERROR;
			WRITE_ERROR(SplineInterpolateImage, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
			return(*Status);
	}
	switch (Degree) {
		case 0L:
		case 1L:
			wy = (double *)NULL;
			wx = (double *)NULL;
			ky = (long *)NULL;
			kx = (long *)NULL;
			break;
		default:
			AllocateLineDouble(&wy, Width, Status);
			if (*Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
				return(*Status);
			}
			AllocateLineDouble(&wx, Width, Status);
			if (*Status == ERROR) {
				FreeLineDouble(&wy);
/**/			DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
				return(*Status);
			}
			ky = (long *)malloc((size_t)(Width * (long)sizeof(long)));
			*Status = (ky == (long *)NULL);
			if (*Status == ERROR) {
				FreeLineDouble(&wx);
				FreeLineDouble(&wy);
/**/			DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
				return(*Status);
			}
			kx = (long *)malloc((size_t)(Width * (long)sizeof(long)));
			*Status = (kx == (long *)NULL);
			if (*Status == ERROR) {
				free(ky);
				FreeLineDouble(&wx);
				FreeLineDouble(&wy);
/**/			DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
				return(*Status);
			}
			break;
	}
	switch (Degree) {
		case 0L:
			j = ConvertDoubleToLong(Yargument);
			i = ConvertDoubleToLong(Xargument);
			if (fabs(Xargument - (double)i) == 0.5) {
				if (fabs(Yargument - (double)j) == 0.5) {
					SplineInterpolateImage(ImageCoeff, Nx, Ny, Xargument, Yargument,
						Result, 1L, Convention, Status);
				}
				else {
					SplineInterpolateImage(ImageCoeff, Nx, Ny, Xargument, (double)j,
						Result, 1L, Convention, Status);
				}
			}
			else {
				if (fabs(Yargument - (double)j) == 0.5) {
					SplineInterpolateImage(ImageCoeff, Nx, Ny, (double)i, Yargument,
						Result, 1L, Convention, Status);
				}
				else {
					if ((0L <= i) && (i < Nx) && (0L <= j) && (j < Ny)) {
						*Result = (double)ImageCoeff[Nx * j + i];
					}
					else {
						switch (Convention) {
							case FiniteDataSupport:
								*Result = 0.0;
								break;
							case MirrorOffBounds:
							case MirrorOnBounds:
							case Periodic:
								*Status = GetFoldedIndex(j, &y, Ny, Convention);
								if (*Status == ERROR) {
									break;
								}
								*Status = GetFoldedIndex(i, &x, Nx, Convention);
								if (*Status == ERROR) {
									break;
								}
								*Result = (double)ImageCoeff[Nx * y + x];
								break;
							default:
								*Status = ERROR;
								WRITE_ERROR(SplineInterpolateImage,
									"Invalid boundary convention")
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
								return(*Status);
						}
					}
				}
			}
/**/		DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
			return(*Status);
		case 1L:
			y = (long)floor(Yargument);
			x = (long)floor(Xargument);
			Yfrac = Yargument - (double)y;
			Xfrac = Xargument - (double)x;
			if ((0L <= x) && (x < (Nx - 1L)) && (0L <= y) && (y < (Ny - 1L))) {
				k = Nx * y + x;
				w0 = Xfrac * (double)ImageCoeff[k + 1L] + (1.0 - Xfrac)
					* (double)ImageCoeff[k];
				k += Nx;
				w1 = Xfrac * (double)ImageCoeff[k + 1L] + (1.0 - Xfrac)
					* (double)ImageCoeff[k];
				*Result = Yfrac * w1 + (1.0 - Yfrac) * w0;
			}
			else {
				switch (Convention) {
					case FiniteDataSupport:
						if ((x < -1L) || (Nx <= x) || (y < -1L) || (Ny <= y)) {
							*Result = 0.0;
							break;
						}
						if (x == -1L) {
							if (y == -1L) {
								*Result = Xfrac * Yfrac * (double)*ImageCoeff;
								break;
							}
							if (y == (Ny - 1L)) {
								*Result = Xfrac * (1.0 - Yfrac)
									* (double)ImageCoeff[Nx * y];
								break;
							}
							y *= Nx;
							w0 = (double)ImageCoeff[y];
							y += Nx;
							w1 = (double)ImageCoeff[y];
							*Result = Xfrac * (Yfrac * w1 + (1.0 - Yfrac) * w0);
							break;
						}
						if (x == (Nx - 1L)) {
							if (y == -1L) {
								*Result = (1.0 - Xfrac) * Yfrac * (double)ImageCoeff[x];
								break;
							}
							if (y == (Ny - 1L)) {
								*Result = (1.0 - Xfrac) * (1.0 - Yfrac)
									* (double)ImageCoeff[Nx * Ny - 1L];
								break;
							}
							k = Nx * y + x;
							w0 = (double)ImageCoeff[k];
							k += Nx;
							w1 = (double)ImageCoeff[k];
							*Result = (1.0 - Xfrac) * (Yfrac * w1 + (1.0 - Yfrac) * w0);
							break;
						}
						if (y == -1L) {
							*Result = Yfrac * (Xfrac * (double)ImageCoeff[x + 1L]
								+ (1.0 - Xfrac) * (double)ImageCoeff[x]);
							break;
						}
						if (y == (Ny - 1L)) {
							k = Nx * y + x;
							*Result = (1.0 - Yfrac) * (Xfrac * (double)ImageCoeff[k + 1L]
								+ (1.0 - Xfrac) * (double)ImageCoeff[k]);
							break;
						}
						*Status = ERROR;
						WRITE_ERROR(SplineInterpolateImage, "Unexpected internal condition")
/**/					DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
						return(*Status);
					case MirrorOffBounds:
					case MirrorOnBounds:
					case Periodic:
						*Status = GetFoldedIndex(y, &y0, Ny, Convention);
						if (*Status == ERROR) {
							break;
						}
						*Status = GetFoldedIndex(y + 1L, &y1, Ny, Convention);
						if (*Status == ERROR) {
							break;
						}
						*Status = GetFoldedIndex(x, &x0, Nx, Convention);
						if (*Status == ERROR) {
							break;
						}
						*Status = GetFoldedIndex(x + 1L, &x1, Nx, Convention);
						if (*Status == ERROR) {
							break;
						}
						y0 *= Nx;
						y1 *= Nx;
						w0 = Xfrac * (double)ImageCoeff[y0 + x1]
							+ (1.0 - Xfrac) * (double)ImageCoeff[y0 + x0];
						w1 = Xfrac * (double)ImageCoeff[y1 + x1]
							+ (1.0 - Xfrac) * (double)ImageCoeff[y1 + x0];
						*Result = Yfrac * w1 + (1.0 - Yfrac) * w0;
						break;
					default:
						*Status = ERROR;
						WRITE_ERROR(SplineInterpolateImage, "Invalid boundary convention")
/**/					DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
						return(*Status);
				}
			}
/**/		DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
			return(*Status);
		case 2L:
			y = ConvertDoubleToLong(Yargument);
			x = ConvertDoubleToLong(Xargument);
			BsplineArray02(Yargument - (double)y, wy, &(wy[1]), &(wy[2]));
			v = ky;
			*v++ = y + 1L;
			*v++ = y;
			*v = y - 1L;
			BsplineArray02(Xargument - (double)x, wx, &(wx[1]), &(wx[2]));
			v = kx;
			*v++ = x + 1L;
			*v++ = x;
			*v = x - 1L;
			if ((x < 1L) || ((Nx - 1L) <= x) || (y < 1L) || ((Ny - 1L) <= y)) {
				switch (Convention) {
					case FiniteDataSupport:
						for (k = 0L; (k < 3L); k++) {
							if (ky[k] < 0L) {
								wy[k] = 0.0;
								ky[k] = 0L;
							}
							else if (Ny <= ky[k]) {
								wy[k] = 0.0;
								ky[k] = 0L;
							}
							if (kx[k] < 0L) {
								wx[k] = 0.0;
								kx[k] = 0L;
							}
							else if (Nx <= kx[k]) {
								wx[k] = 0.0;
								kx[k] = 0L;
							}
						}
						break;
					case MirrorOffBounds:
					case MirrorOnBounds:
					case Periodic:
						for (k = 0L; (k < 3L); k++) {
							GetFoldedIndex(ky[k], &l, Ny, Convention);
							ky[k] = l;
							GetFoldedIndex(kx[k], &l, Nx, Convention);
							kx[k] = l;
						}
						break;
					default:
						*Status = ERROR;
						free(kx);
						free(ky);
						FreeLineDouble(&wx);
						FreeLineDouble(&wy);
						WRITE_ERROR(SplineInterpolateImage, "Invalid boundary convention")
/**/					DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
						return(*Status);
				}
			}
			uy = wy;
			vy = ky;
			ux = wx;
			vx = kx;
			pi = ImageCoeff + (ptrdiff_t)(Nx * *vy++);
			qi = *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux * (double)*(pi+ (ptrdiff_t)*vx);
			qj = *uy++ * qi;
			ux = wx;
			vx = kx;
			pi = ImageCoeff + (ptrdiff_t)(Nx * *vy++);
			qi = *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux * (double)*(pi+ (ptrdiff_t)*vx);
			qj += *uy++ * qi;
			ux = wx;
			vx = kx;
			pi = ImageCoeff + (ptrdiff_t)(Nx * *vy);
			qi = *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux * (double)*(pi+ (ptrdiff_t)*vx);
			qj += *uy * qi;
			*Result = qj;
			free(kx);
			free(ky);
			*Status = FreeLineDouble(&wx);
			if (*Status == ERROR) {
				FreeLineDouble(&wy);
/**/			DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
				return(*Status);
			}
			*Status = FreeLineDouble(&wy);
/**/		DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
			return(*Status);
		case 3L:
			y = ConvertDoubleToLong(Yargument);
			x = ConvertDoubleToLong(Xargument);
			Yargument -= (double)y;
			Xargument -= (double)x;
			if (Yargument < 0.0) {
				BsplineArray03(Yargument + 1.0, wy, &(wy[1]), &(wy[2]), &(wy[3]));
				v = ky;
				*v++ = ++y;
				*v++ = --y;
				*v++ = --y;
				*v = --y;
			}
			else {
				BsplineArray03(Yargument, wy, &(wy[1]), &(wy[2]), &(wy[3]));
				v = ky;
				y += 2L;
				*v++ = y;
				*v++ = --y;
				*v++ = --y;
				*v = --y;
			}
			if (Xargument < 0.0) {
				BsplineArray03(Xargument + 1.0, wx, &(wx[1]), &(wx[2]), &(wx[3]));
				v = kx;
				*v++ = ++x;
				*v++ = --x;
				*v++ = --x;
				*v = --x;
			}
			else {
				BsplineArray03(Xargument, wx, &(wx[1]), &(wx[2]), &(wx[3]));
				v = kx;
				x += 2L;
				*v++ = --x;
				*v++ = --x;
				*v++ = --x;
				*v = --x;
			}
			if ((kx[3] < 0L) || (Nx <= kx[0]) || (ky[3] < 0L) || (Ny <= ky[0])) {
				switch (Convention) {
					case FiniteDataSupport:
						for (k = 0L; (k < 4L); k++) {
							if (ky[k] < 0L) {
								wy[k] = 0.0;
								ky[k] = 0L;
							}
							else if (Ny <= ky[k]) {
								wy[k] = 0.0;
								ky[k] = 0L;
							}
							if (kx[k] < 0L) {
								wx[k] = 0.0;
								kx[k] = 0L;
							}
							else if (Nx <= kx[k]) {
								wx[k] = 0.0;
								kx[k] = 0L;
							}
						}
						break;
					case MirrorOffBounds:
					case MirrorOnBounds:
					case Periodic:
						for (k = 0L; (k < 4L); k++) {
							GetFoldedIndex(ky[k], &l, Ny, Convention);
							ky[k] = l;
							GetFoldedIndex(kx[k], &l, Nx, Convention);
							kx[k] = l;
						}
						break;
					default:
						*Status = ERROR;
						free(kx);
						free(ky);
						FreeLineDouble(&wx);
						FreeLineDouble(&wy);
						WRITE_ERROR(SplineInterpolateImage, "Invalid boundary convention")
/**/					DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
						return(*Status);
				}
			}
			uy = wy;
			vy = ky;
			ux = wx;
			vx = kx;
			pi = ImageCoeff + (ptrdiff_t)(Nx * *vy++);
			qi = *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux * (double)*(pi+ (ptrdiff_t)*vx);
			qj = *uy++ * qi;
			ux = wx;
			vx = kx;
			pi = ImageCoeff + (ptrdiff_t)(Nx * *vy++);
			qi = *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux * (double)*(pi+ (ptrdiff_t)*vx);
			qj += *uy++ * qi;
			ux = wx;
			vx = kx;
			pi = ImageCoeff + (ptrdiff_t)(Nx * *vy++);
			qi = *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux * (double)*(pi+ (ptrdiff_t)*vx);
			qj += *uy++ * qi;
			ux = wx;
			vx = kx;
			pi = ImageCoeff + (ptrdiff_t)(Nx * *vy);
			qi = *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
			qi += *ux * (double)*(pi+ (ptrdiff_t)*vx);
			qj += *uy * qi;
			*Result = qj;
			free(kx);
			free(ky);
			*Status = FreeLineDouble(&wx);
			if (*Status == ERROR) {
				FreeLineDouble(&wy);
/**/			DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
				return(*Status);
			}
			*Status = FreeLineDouble(&wy);
/**/		DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
			return(*Status);
		case 4L:
			y = ConvertDoubleToLong(Yargument) + 2L;
			x = ConvertDoubleToLong(Xargument) + 2L;
			Yargument -= (double)y;
			Xargument -= (double)x;
			for (k = 0L; (k < 5L); k++) {
				wy[k] = Bspline04(Yargument++);
				ky[k] = y--;
				wx[k] = Bspline04(Xargument++);
				ky[k] = x--;
			}
			break;
		case 5L:
			y = ConvertDoubleToLong(Yargument);
			x = ConvertDoubleToLong(Xargument);
			Yargument -= (double)y;
			Xargument -= (double)x;
			if (Yargument < 0.0) {
				Yargument -= 2.0;
				y += 2L;
			}
			else {
				Yargument -= 3.0;
				y += 3L;
			}
			if (Xargument < 0.0) {
				Xargument -= 2.0;
				x += 2L;
			}
			else {
				Xargument -= 3.0;
				x += 3L;
			}
			for (k = 0L; (k < 6L); k++) {
				wy[k] = Bspline05(Yargument++);
				ky[k] = y--;
				wx[k] = Bspline05(Xargument++);
				kx[k] = x--;
			}
			break;
		case 6L:
			y = ConvertDoubleToLong(Yargument) + 3L;
			x = ConvertDoubleToLong(Xargument) + 3L;
			Yargument -= (double)y;
			Xargument -= (double)x;
			for (k = 0L; (k < 7L); k++) {
				wy[k] = Bspline06(Yargument++);
				ky[k] = y--;
				wx[k] = Bspline06(Xargument++);
				ky[k] = x--;
			}
			break;
		case 7L:
			y = ConvertDoubleToLong(Yargument);
			x = ConvertDoubleToLong(Xargument);
			Yargument -= (double)y;
			Xargument -= (double)x;
			if (Yargument < 0.0) {
				Yargument -= 3.0;
				y += 3L;
			}
			else {
				Yargument -= 4.0;
				y += 4L;
			}
			if (Xargument < 0.0) {
				Xargument -= 3.0;
				x += 3L;
			}
			else {
				Xargument -= 4.0;
				x += 4L;
			}
			for (k = 0L; (k < 8L); k++) {
				wy[k] = Bspline07(Yargument++);
				ky[k] = y--;
				wx[k] = Bspline07(Xargument++);
				kx[k] = x--;
			}
			break;
		case 8L:
			y = ConvertDoubleToLong(Yargument) + 4L;
			x = ConvertDoubleToLong(Xargument) + 4L;
			Yargument -= (double)y;
			Xargument -= (double)x;
			for (k = 0L; (k < 9L); k++) {
				wy[k] = Bspline08(Yargument++);
				ky[k] = y--;
				wx[k] = Bspline08(Xargument++);
				ky[k] = x--;
			}
			break;
		case 9L:
			y = ConvertDoubleToLong(Yargument);
			x = ConvertDoubleToLong(Xargument);
			Yargument -= (double)y;
			Xargument -= (double)x;
			if (Yargument < 0.0) {
				Yargument -= 4.0;
				y += 4L;
			}
			else {
				Yargument -= 5.0;
				y += 5L;
			}
			if (Xargument < 0.0) {
				Xargument -= 4.0;
				x += 4L;
			}
			else {
				Xargument -= 5.0;
				x += 5L;
			}
			for (k = 0L; (k < 10L); k++) {
				wy[k] = Bspline09(Yargument++);
				ky[k] = y--;
				wx[k] = Bspline09(Xargument++);
				kx[k] = x--;
			}
			break;
		case 10L:
			y = ConvertDoubleToLong(Yargument) + 5L;
			x = ConvertDoubleToLong(Xargument) + 5L;
			Yargument -= (double)y;
			Xargument -= (double)x;
			for (k = 0L; (k < 11L); k++) {
				wy[k] = Bspline10(Yargument++);
				ky[k] = y--;
				wx[k] = Bspline10(Xargument++);
				ky[k] = x--;
			}
			break;
		case 11L:
			y = ConvertDoubleToLong(Yargument);
			x = ConvertDoubleToLong(Xargument);
			Yargument -= (double)y;
			Xargument -= (double)x;
			if (Yargument < 0.0) {
				Yargument -= 5.0;
				y += 5L;
			}
			else {
				Yargument -= 6.0;
				y += 6L;
			}
			if (Xargument < 0.0) {
				Xargument -= 5.0;
				x += 5L;
			}
			else {
				Xargument -= 6.0;
				x += 6L;
			}
			for (k = 0L; (k < 12L); k++) {
				wy[k] = Bspline11(Yargument++);
				ky[k] = y--;
				wx[k] = Bspline11(Xargument++);
				kx[k] = x--;
			}
			break;
		default:
			y = ConvertDoubleToLong(Yargument);
			x = ConvertDoubleToLong(Xargument);
			Yargument -= (double)y;
			Xargument -= (double)x;
			if (Degree & 1L) {
				if (Yargument < 0.0) {
					Yargument -= (double)HalfDegree;
					y += HalfDegree;
				}
				else {
					Yargument -= (double)HalfDegree + 1.0;
					y += HalfDegree + 1L;
				}
				if (Xargument < 0.0) {
					Xargument -= (double)HalfDegree;
					x += HalfDegree;
				}
				else {
					Xargument -= (double)HalfDegree + 1.0;
					x += HalfDegree + 1L;
				}
			}
			for (k = 0L; (k < Width); k++) {
				*Status = Bspline(Degree, Yargument++, &(wy[k]));
				if (*Status == ERROR) {
					free(kx);
					free(ky);
					FreeLineDouble(&wx);
					FreeLineDouble(&wy);
/**/				DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
					return(*Status);
				}
				ky[k] = y--;
				*Status = Bspline(Degree, Xargument++, &(wx[k]));
				if (*Status == ERROR) {
					free(kx);
					free(ky);
					FreeLineDouble(&wx);
					FreeLineDouble(&wy);
/**/				DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
					return(*Status);
				}
				kx[k] = x--;
			}
			break;
	}
	if ((kx[Degree] < 0L) || (Nx <= kx[0]) || (ky[Degree] < 0L) || (Ny <= ky[0])) {
		switch (Convention) {
			case FiniteDataSupport:
				for (k = 0L; (k < Width); k++) {
					if (ky[k] < 0L) {
						wy[k] = 0.0;
						ky[k] = 0L;
					}
					else if (Ny <= ky[k]) {
						wy[k] = 0.0;
						ky[k] = 0L;
					}
					if (kx[k] < 0L) {
						wx[k] = 0.0;
						kx[k] = 0L;
					}
					else if (Nx <= kx[k]) {
						wx[k] = 0.0;
						kx[k] = 0L;
					}
				}
				break;
			case MirrorOffBounds:
			case MirrorOnBounds:
			case Periodic:
				for (k = 0L; (k < Width); k++) {
					GetFoldedIndex(ky[k], &l, Ny, Convention);
					ky[k] = l;
					GetFoldedIndex(kx[k], &l, Nx, Convention);
					kx[k] = l;
				}
				break;
			default:
				*Status = ERROR;
				free(kx);
				free(ky);
				FreeLineDouble(&wx);
				FreeLineDouble(&wy);
				WRITE_ERROR(SplineInterpolateImage, "Invalid boundary convention")
/**/			DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
				return(*Status);
		}
	}
	uy = wy;
	vy = ky;
	qj = 0.0;
	for (j = 0L; (j < Width); j++) {
		ux = wx;
		vx = kx;
		pi = ImageCoeff + (ptrdiff_t)(Nx * *vy++);
		qi = 0.0;
		for (i = 0L; (i < Width); i++) {
			qi += *ux++ * (double)*(pi + (ptrdiff_t)*vx++);
		}
		qj += *uy++ * qi;
	}
	*Result = qj;
	free(kx);
	free(ky);
	*Status = FreeLineDouble(&wx);
	if (*Status == ERROR) {
		FreeLineDouble(&wy);
/**/	DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
		return(*Status);
	}
	*Status = FreeLineDouble(&wy);
/**/DEBUG_WRITE_LEAVING(SplineInterpolateImage, "Done")
	return(*Status);
} /* SplineInterpolateImage */

/*--------------------------------------------------------------------------*/
extern int		SplineInterpolateLine
				(
					double	LineCoeff[],		/* B-spline coefficients to interpolate */
					long	LineLength,			/* length of the line */
					double	Argument,			/* input abscissa */
					double	*Result,			/* output ordinate */
					long	Degree,				/* degree of the spline */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

/* interpolation of the 1D array LineCoeff[] of length LineLength */
/* Result = SUM(k): LineCoeff[k] * B-spline(Degree, Argument - k) */
/* Convention must be consistent with the B-spline coefficients stored in LineCoeff[] */
/* success: return(!ERROR); failure: return(ERROR) */

{ /* begin SplineInterpolateLine */

	double	*p, *x = (double *)NULL, *y = (double *)NULL;
	long	*n;
	double	Weight[12], Value[12];
	long	Index[12];
	double	u, v, w;
	long	i, j, k;

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(SplineInterpolateLine, LineCoeff, *Status,
/**/	"Missing LineCoeff")
/**/DEBUG_CHECK_NULL_POINTER(SplineInterpolateLine, Result, *Status,
/**/	"Missing Result")
/**/DEBUG_CHECK_RANGE_LONG(SplineInterpolateLine, LineLength, 1L, LONG_MAX, *Status,
/**/	"Invalid LineLength (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SplineInterpolateLine, Degree, 0L, LONG_MAX, *Status,
/**/	"Invalid Degree (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(SplineInterpolateLine, *Status)
/**/DEBUG_WRITE_ENTERING(SplineInterpolateLine,
/**/	"About to execute SplineInterpolateLine")

	switch (Convention) {
		case AntiMirrorOnBounds:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case FiniteDataSupport:
			if (2L <= Degree) {
				*Status = ERROR;
				WRITE_ERROR(SplineInterpolateLine, "Invalid boundary convention")
/**/			DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
				return(*Status);
			}
			break;
		case FiniteCoefficientSupport:
			Convention = FiniteDataSupport;
			break;
		default:
			*Status = ERROR;
			WRITE_ERROR(SplineInterpolateLine, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
			return(*Status);
	}
	switch (Degree) {
		case 0L:
			i = ConvertDoubleToLong(Argument);
			if (fabs(Argument - (double)i) == 0.5) {
				SplineInterpolateLine(LineCoeff, LineLength, Argument, Result, 1L,
					Convention, Status);
			}
			else {
				if ((0L <= i) && (i < LineLength)) {
					*Result = LineCoeff[i];
				}
				else {
					switch (Convention) {
						case AntiMirrorOnBounds:
						case FiniteDataSupport:
							*Status = GetFoldedValueDouble(LineCoeff, i, Result, LineLength,
								Convention);
							break;
						case MirrorOffBounds:
						case MirrorOnBounds:
						case Periodic:
							*Status = GetFoldedIndex(i, &j, LineLength, Convention);
							if (*Status == !ERROR) {
								*Result = (double)LineCoeff[j];
							}
							break;
						default:
							*Status = ERROR;
							WRITE_ERROR(SplineInterpolateLine,
								"Unexpected boundary convention")
							break;
					}
				}
			}
			break;
		case 1L:
			i = (long)floor(Argument);
			w = Argument - (double)i;
			if ((0L <= i) && (i < (LineLength - 1L))) {
				*Result = w * LineCoeff[i + 1L] + (1.0 - w) * LineCoeff[i];
			}
			else {
				switch (Convention) {
					case AntiMirrorOnBounds:
					case FiniteDataSupport:
						*Status = GetFoldedValueDouble(LineCoeff, i + 1L, &u, LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedValueDouble(LineCoeff, i, &v, LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Result = w * u + (1.0 - w) * v;
						break;
					case MirrorOffBounds:
					case MirrorOnBounds:
					case Periodic:
						*Status = GetFoldedIndex(i + 1L, &j, LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedIndex(i, &k, LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Result = w * LineCoeff[j] + (1.0 - w) * LineCoeff[k];
						break;
					default:
						*Status = ERROR;
						WRITE_ERROR(SplineInterpolateLine, "Unexpected boundary convention")
						break;
				}
			}
			break;
		case 2L:
			i = ConvertDoubleToLong(Argument);
			*Status = BsplineArray02(Argument - (double)i,
				&(Weight[0]), &(Weight[1]), &(Weight[2]));
			if (*Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
				return(*Status);
			}
			if ((1L <= i) && (i < (LineLength - 1L))) {
				*Result = Weight[0] * LineCoeff[i + 1L] + Weight[1] * LineCoeff[i]
					+ Weight[2] * LineCoeff[i - 1L];
			}
			else {
				switch (Convention) {
					case AntiMirrorOnBounds:
					case FiniteDataSupport:
						*Status = GetFoldedValueDouble(LineCoeff, i + 1L, &(Value[2]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedValueDouble(LineCoeff, i, &(Value[1]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedValueDouble(LineCoeff, i - 1L, &(Value[0]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Result = Weight[0] * Value[2] + Weight[1] * Value[1]
							+ Weight[2] * Value[0];
						break;
					case MirrorOffBounds:
					case MirrorOnBounds:
					case Periodic:
						*Status = GetFoldedIndex(i + 1L, &(Index[2]), LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedIndex(i, &(Index[1]), LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedIndex(i - 1L, &(Index[0]), LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Result = Weight[0] * LineCoeff[Index[2]] + Weight[1]
							* LineCoeff[Index[1]] + Weight[2] * LineCoeff[Index[0]];
						break;
					default:
						*Status = ERROR;
						WRITE_ERROR(SplineInterpolateLine, "Unexpected boundary convention")
						break;
				}
			}
			break;
		case 3L:
			i = (long)floor(Argument);
			*Status = BsplineArray03(Argument - (double)i,
				&(Weight[0]), &(Weight[1]), &(Weight[2]), &(Weight[3]));
			if (*Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
				return(*Status);
			}
			if ((1L <= i) && (i < (LineLength - 2L))) {
				*Result = Weight[0] * LineCoeff[i + 2L] + Weight[1] * LineCoeff[i + 1L]
					+ Weight[2] * LineCoeff[i] + Weight[3] * LineCoeff[i - 1L];
			}
			else {
				switch (Convention) {
					case AntiMirrorOnBounds:
					case FiniteDataSupport:
						*Status = GetFoldedValueDouble(LineCoeff, i + 2L, &(Value[3]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedValueDouble(LineCoeff, i + 1L, &(Value[2]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedValueDouble(LineCoeff, i, &(Value[1]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedValueDouble(LineCoeff, i - 1L, &(Value[0]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Result = Weight[0] * Value[3] + Weight[1] * Value[2]
							+ Weight[2] * Value[1] + Weight[3] * Value[0];
						break;
					case MirrorOffBounds:
					case MirrorOnBounds:
					case Periodic:
						*Status = GetFoldedIndex(i + 2L, &(Index[3]), LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedIndex(i + 1L, &(Index[2]), LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedIndex(i, &(Index[1]), LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedIndex(i - 1L, &(Index[0]), LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Result = Weight[0] * LineCoeff[Index[3]] + Weight[1]
							* LineCoeff[Index[2]] + Weight[2] * LineCoeff[Index[1]]
							+ Weight[3] * LineCoeff[Index[0]];
						break;
					default:
						*Status = ERROR;
						WRITE_ERROR(SplineInterpolateLine, "Unexpected boundary convention")
						break;
				}
			}
			break;
		case 4L:
			i = ConvertDoubleToLong(Argument);
			w = Argument - (double)i;
			Weight[0] = Bspline04(w - 2.0);
			Weight[1] = Bspline04(w - 1.0);
			Weight[2] = Bspline04(w);
			Weight[3] = Bspline04(w + 1.0);
			Weight[4] = Bspline04(w + 2.0);
			if ((2L <= i) && (i < (LineLength - 2L))) {
				*Result = Weight[0] * LineCoeff[i + 2L] + Weight[1] * LineCoeff[i + 1L]
					+ Weight[2] * LineCoeff[i] + Weight[3] * LineCoeff[i - 1L]
					+ Weight[4] * LineCoeff[i - 2L];
			}
			else {
				switch (Convention) {
					case AntiMirrorOnBounds:
					case FiniteDataSupport:
						*Status = GetFoldedValueDouble(LineCoeff, i + 2L, &(Value[4]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedValueDouble(LineCoeff, i + 1L, &(Value[3]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedValueDouble(LineCoeff, i, &(Value[2]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedValueDouble(LineCoeff, i - 1L, &(Value[1]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedValueDouble(LineCoeff, i - 2L, &(Value[0]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Result = Weight[0] * Value[4] + Weight[1] * Value[3]
							+ Weight[2] * Value[2] + Weight[3] * Value[1]
							+ Weight[4] * Value[0];
						break;
					case MirrorOffBounds:
					case MirrorOnBounds:
					case Periodic:
						*Status = GetFoldedIndex(i + 2L, &(Index[4]), LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedIndex(i + 1L, &(Index[3]), LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedIndex(i, &(Index[2]), LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedIndex(i - 1L, &(Index[1]), LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedIndex(i - 2L, &(Index[0]), LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Result = Weight[0] * LineCoeff[Index[4]] + Weight[1]
							* LineCoeff[Index[3]] + Weight[2] * LineCoeff[Index[2]]
							+ Weight[3] * LineCoeff[Index[1]] + Weight[4]
							* LineCoeff[Index[0]];
						break;
					default:
						*Status = ERROR;
						WRITE_ERROR(SplineInterpolateLine, "Unexpected boundary convention")
						break;
				}
			}
			break;
		case 5L:
			i = (long)floor(Argument);
			w = Argument - (double)i;
			Weight[0] = Bspline05(w - 3.0);
			Weight[1] = Bspline05(w - 2.0);
			Weight[2] = Bspline05(w - 1.0);
			Weight[3] = Bspline05(w);
			Weight[4] = Bspline05(w + 1.0);
			Weight[5] = Bspline05(w + 2.0);
			if ((2L <= i) && (i < (LineLength - 3L))) {
				*Result = Weight[0] * LineCoeff[i + 3L] + Weight[1] * LineCoeff[i + 2L]
					+ Weight[2] * LineCoeff[i + 1L] + Weight[3] * LineCoeff[i]
					+ Weight[4] * LineCoeff[i - 1L] + Weight[5] * LineCoeff[i - 2L];
			}
			else {
				switch (Convention) {
					case AntiMirrorOnBounds:
					case FiniteDataSupport:
						*Status = GetFoldedValueDouble(LineCoeff, i + 3L, &(Value[5]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedValueDouble(LineCoeff, i + 2L, &(Value[4]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedValueDouble(LineCoeff, i + 1L, &(Value[3]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedValueDouble(LineCoeff, i, &(Value[2]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedValueDouble(LineCoeff, i - 1L, &(Value[1]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedValueDouble(LineCoeff, i - 2L, &(Value[0]),
							LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Result = Weight[0] * Value[5] + Weight[1] * Value[4]
							+ Weight[2] * Value[3] + Weight[3] * Value[2]
							+ Weight[4] * Value[1] + Weight[5] * Value[0];
						break;
					case MirrorOffBounds:
					case MirrorOnBounds:
					case Periodic:
						*Status = GetFoldedIndex(i + 3L, &(Index[5]), LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedIndex(i + 2L, &(Index[4]), LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedIndex(i + 1L, &(Index[3]), LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedIndex(i, &(Index[2]), LineLength, Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedIndex(i - 1L, &(Index[1]), LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Status = GetFoldedIndex(i - 2L, &(Index[0]), LineLength,
							Convention);
						if (*Status == ERROR) {
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
						}
						*Result = Weight[0] * LineCoeff[Index[5]] + Weight[1]
							* LineCoeff[Index[4]] + Weight[2] * LineCoeff[Index[3]]
							+ Weight[3] * LineCoeff[Index[2]] + Weight[4]
							* LineCoeff[Index[1]] + Weight[5] * LineCoeff[Index[0]];
						break;
					default:
						*Status = ERROR;
						WRITE_ERROR(SplineInterpolateLine, "Unexpected boundary convention")
						break;
				}
			}
			break;
		case 6L:
			i = ConvertDoubleToLong(Argument);
			w = Argument - (double)i;
			Weight[0] = Bspline06(w - 3.0);
			Weight[1] = Bspline06(w - 2.0);
			Weight[2] = Bspline06(w - 1.0);
			Weight[3] = Bspline06(w);
			Weight[4] = Bspline06(w + 1.0);
			Weight[5] = Bspline06(w + 2.0);
			Weight[6] = Bspline06(w + 3.0);
			if ((3L <= i) && (i < (LineLength - 3L))) {
				*Result = Weight[0] * LineCoeff[i + 3L] + Weight[1] * LineCoeff[i + 2L]
					+ Weight[2] * LineCoeff[i + 1L] + Weight[3] * LineCoeff[i]
					+ Weight[4] * LineCoeff[i - 1L] + Weight[5] * LineCoeff[i - 2L]
					+ Weight[6] * LineCoeff[i - 3L];
			}
			else {
				switch (Convention) {
					case AntiMirrorOnBounds:
					case FiniteDataSupport:
						k = Degree;
						for (j = i + Degree / 2L; ((i - Degree / 2L) <= j); j--) {
							*Status = GetFoldedValueDouble(LineCoeff, j, &(Value[k--]),
								LineLength, Convention);
							if (*Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
						}
						*Result = Weight[0] * Value[6] + Weight[1] * Value[5]
							+ Weight[2] * Value[4] + Weight[3] * Value[3]
							+ Weight[4] * Value[2] + Weight[5] * Value[1]
							+ Weight[6] * Value[0];
						break;
					case MirrorOffBounds:
					case MirrorOnBounds:
					case Periodic:
						k = Degree;
						for (j = i + Degree / 2L; ((i - Degree / 2L) <= j); j--) {
							*Status = GetFoldedIndex(j, &(Index[k--]), LineLength,
								Convention);
							if (*Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
						}
						*Result = Weight[0] * LineCoeff[Index[6]] + Weight[1]
							* LineCoeff[Index[5]] + Weight[2] * LineCoeff[Index[4]]
							+ Weight[3] * LineCoeff[Index[3]] + Weight[4]
							* LineCoeff[Index[2]] + Weight[5] * LineCoeff[Index[1]]
							+ Weight[6] * LineCoeff[Index[0]];
						break;
					default:
						*Status = ERROR;
						WRITE_ERROR(SplineInterpolateLine, "Unexpected boundary convention")
						break;
				}
			}
			break;
		case 7L:
			i = (long)floor(Argument);
			w = Argument - (double)i;
			Weight[0] = Bspline07(w - 4.0);
			Weight[1] = Bspline07(w - 3.0);
			Weight[2] = Bspline07(w - 2.0);
			Weight[3] = Bspline07(w - 1.0);
			Weight[4] = Bspline07(w);
			Weight[5] = Bspline07(w + 1.0);
			Weight[6] = Bspline07(w + 2.0);
			Weight[7] = Bspline07(w + 3.0);
			if ((3L <= i) && (i < (LineLength - 4L))) {
				*Result = Weight[0] * LineCoeff[i + 4L] + Weight[1] * LineCoeff[i + 3L]
					+ Weight[2] * LineCoeff[i + 2L] + Weight[3] * LineCoeff[i + 1L]
					+ Weight[4] * LineCoeff[i] + Weight[5] * LineCoeff[i - 1L]
					+ Weight[6] * LineCoeff[i - 2L] + Weight[7] * LineCoeff[i - 3L];
			}
			else {
				switch (Convention) {
					case AntiMirrorOnBounds:
					case FiniteDataSupport:
						k = Degree;
						for (j = i + (Degree + 1L) / 2L; ((i - (Degree - 1L) / 2L) <= j);
							j--) {
							*Status = GetFoldedValueDouble(LineCoeff, j, &(Value[k--]),
								LineLength, Convention);
							if (*Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
						}
						*Result = Weight[0] * Value[7] + Weight[1] * Value[6]
							+ Weight[2] * Value[5] + Weight[3] * Value[4]
							+ Weight[4] * Value[3] + Weight[5] * Value[2]
							+ Weight[6] * Value[1] + Weight[7] * Value[0];
						break;
					case MirrorOffBounds:
					case MirrorOnBounds:
					case Periodic:
						k = Degree;
						for (j = i + (Degree + 1L) / 2L; ((i - (Degree - 1L) / 2L) <= j);
							j--) {
							*Status = GetFoldedIndex(j, &(Index[k--]), LineLength,
								Convention);
							if (*Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
						}
						*Result = Weight[0] * LineCoeff[Index[7]] + Weight[1]
							* LineCoeff[Index[6]] + Weight[2] * LineCoeff[Index[5]]
							+ Weight[3] * LineCoeff[Index[4]] + Weight[4]
							* LineCoeff[Index[3]] + Weight[5] * LineCoeff[Index[2]]
							+ Weight[6] * LineCoeff[Index[1]] + Weight[7]
							* LineCoeff[Index[0]];
						break;
					default:
						*Status = ERROR;
						WRITE_ERROR(SplineInterpolateLine, "Unexpected boundary convention")
						break;
				}
			}
			break;
		case 8L:
			i = ConvertDoubleToLong(Argument);
			w = Argument - (double)i;
			Weight[0] = Bspline08(w - 4.0);
			Weight[1] = Bspline08(w - 3.0);
			Weight[2] = Bspline08(w - 2.0);
			Weight[3] = Bspline08(w - 1.0);
			Weight[4] = Bspline08(w);
			Weight[5] = Bspline08(w + 1.0);
			Weight[6] = Bspline08(w + 2.0);
			Weight[7] = Bspline08(w + 3.0);
			Weight[8] = Bspline08(w + 4.0);
			if ((4L <= i) && (i < (LineLength - 4L))) {
				*Result = Weight[0] * LineCoeff[i + 4L] + Weight[1] * LineCoeff[i + 3L]
					+ Weight[2] * LineCoeff[i + 2L] + Weight[3] * LineCoeff[i + 1L]
					+ Weight[4] * LineCoeff[i] + Weight[5] * LineCoeff[i - 1L]
					+ Weight[6] * LineCoeff[i - 2L] + Weight[7] * LineCoeff[i - 3L]
					+ Weight[8] * LineCoeff[i - 4L];
			}
			else {
				switch (Convention) {
					case AntiMirrorOnBounds:
					case FiniteDataSupport:
						k = Degree;
						for (j = i + Degree / 2L; ((i - Degree / 2L) <= j); j--) {
							*Status = GetFoldedValueDouble(LineCoeff, j, &(Value[k--]),
								LineLength, Convention);
							if (*Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
						}
						*Result = Weight[0] * Value[8] + Weight[1] * Value[7]
							+ Weight[2] * Value[6] + Weight[3] * Value[5]
							+ Weight[4] * Value[4] + Weight[5] * Value[3]
							+ Weight[6] * Value[2] + Weight[7] * Value[1]
							+ Weight[8] * Value[0];
						break;
					case MirrorOffBounds:
					case MirrorOnBounds:
					case Periodic:
						k = Degree;
						for (j = i + Degree / 2L; ((i - Degree / 2L) <= j); j--) {
							*Status = GetFoldedIndex(j, &(Index[k--]), LineLength,
								Convention);
							if (*Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
						}
						*Result = Weight[0] * LineCoeff[Index[8]] + Weight[1]
							* LineCoeff[Index[7]] + Weight[2] * LineCoeff[Index[6]]
							+ Weight[3] * LineCoeff[Index[5]] + Weight[4]
							* LineCoeff[Index[4]] + Weight[5] * LineCoeff[Index[3]]
							+ Weight[6] * LineCoeff[Index[2]] + Weight[7]
							* LineCoeff[Index[1]] + Weight[8] * LineCoeff[Index[0]];
						break;
					default:
						*Status = ERROR;
						WRITE_ERROR(SplineInterpolateLine, "Unexpected boundary convention")
						break;
				}
			}
			break;
		case 9L:
			i = (long)floor(Argument);
			w = Argument - (double)i;
			Weight[0] = Bspline09(w - 5.0);
			Weight[1] = Bspline09(w - 4.0);
			Weight[2] = Bspline09(w - 3.0);
			Weight[3] = Bspline09(w - 2.0);
			Weight[4] = Bspline09(w - 1.0);
			Weight[5] = Bspline09(w);
			Weight[6] = Bspline09(w + 1.0);
			Weight[7] = Bspline09(w + 2.0);
			Weight[8] = Bspline09(w + 3.0);
			Weight[9] = Bspline09(w + 4.0);
			if ((4L <= i) && (i < (LineLength - 5L))) {
				*Result = Weight[0] * LineCoeff[i + 5L] + Weight[1] * LineCoeff[i + 4L]
					+ Weight[2] * LineCoeff[i + 3L] + Weight[3] * LineCoeff[i + 2L]
					+ Weight[4] * LineCoeff[i + 1L] + Weight[5] * LineCoeff[i]
					+ Weight[6] * LineCoeff[i - 1L] + Weight[7] * LineCoeff[i - 2L]
					+ Weight[8] * LineCoeff[i - 3L] + Weight[9] * LineCoeff[i - 4L];
			}
			else {
				switch (Convention) {
					case AntiMirrorOnBounds:
					case FiniteDataSupport:
						k = Degree;
						for (j = i + (Degree + 1L) / 2L; ((i - (Degree - 1L) / 2L) <= j);
							j--) {
							*Status = GetFoldedValueDouble(LineCoeff, j, &(Value[k--]),
								LineLength, Convention);
							if (*Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
						}
						*Result = Weight[0] * Value[9] + Weight[1] * Value[8]
							+ Weight[2] * Value[7] + Weight[3] * Value[6]
							+ Weight[4] * Value[5] + Weight[5] * Value[4]
							+ Weight[6] * Value[3] + Weight[7] * Value[2]
							+ Weight[8] * Value[1] + Weight[9] * Value[0];
						break;
					case MirrorOffBounds:
					case MirrorOnBounds:
					case Periodic:
						k = Degree;
						for (j = i + (Degree + 1L) / 2L; ((i - (Degree - 1L) / 2L) <= j);
							j--) {
							*Status = GetFoldedIndex(j, &(Index[k--]), LineLength,
								Convention);
							if (*Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
						}
						*Result = Weight[0] * LineCoeff[Index[9]] + Weight[1]
							* LineCoeff[Index[8]] + Weight[2] * LineCoeff[Index[7]]
							+ Weight[3] * LineCoeff[Index[6]] + Weight[4]
							* LineCoeff[Index[5]] + Weight[5] * LineCoeff[Index[4]]
							+ Weight[6] * LineCoeff[Index[3]] + Weight[7]
							* LineCoeff[Index[2]] + Weight[8] * LineCoeff[Index[1]]
							+ Weight[9] * LineCoeff[Index[0]];
						break;
					default:
						*Status = ERROR;
						WRITE_ERROR(SplineInterpolateLine, "Unexpected boundary convention")
						break;
				}
			}
			break;
		case 10L:
			i = ConvertDoubleToLong(Argument);
			w = Argument - (double)i;
			Weight[0] = Bspline10(w - 5.0);
			Weight[1] = Bspline10(w - 4.0);
			Weight[2] = Bspline10(w - 3.0);
			Weight[3] = Bspline10(w - 2.0);
			Weight[4] = Bspline10(w - 1.0);
			Weight[5] = Bspline10(w);
			Weight[6] = Bspline10(w + 1.0);
			Weight[7] = Bspline10(w + 2.0);
			Weight[8] = Bspline10(w + 3.0);
			Weight[9] = Bspline10(w + 4.0);
			Weight[10] = Bspline10(w + 5.0);
			if ((5L <= i) && (i < (LineLength - 5L))) {
				*Result = Weight[0] * LineCoeff[i + 5L] + Weight[1] * LineCoeff[i + 4L]
					+ Weight[2] * LineCoeff[i + 3L] + Weight[3] * LineCoeff[i + 2L]
					+ Weight[4] * LineCoeff[i + 1L] + Weight[5] * LineCoeff[i]
					+ Weight[6] * LineCoeff[i - 1L] + Weight[7] * LineCoeff[i - 2L]
					+ Weight[8] * LineCoeff[i - 3L] + Weight[9] * LineCoeff[i - 4L]
					+ Weight[10] * LineCoeff[i - 5L];
			}
			else {
				switch (Convention) {
					case AntiMirrorOnBounds:
					case FiniteDataSupport:
						k = Degree;
						for (j = i + Degree / 2L; ((i - Degree / 2L) <= j); j--) {
							*Status = GetFoldedValueDouble(LineCoeff, j, &(Value[k--]),
								LineLength, Convention);
							if (*Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
						}
						*Result = Weight[0] * Value[10] + Weight[1] * Value[9]
							+ Weight[2] * Value[8] + Weight[3] * Value[7]
							+ Weight[4] * Value[6] + Weight[5] * Value[5]
							+ Weight[6] * Value[4] + Weight[7] * Value[3]
							+ Weight[8] * Value[2] + Weight[9] * Value[1]
							+ Weight[10] * Value[0];
						break;
					case MirrorOffBounds:
					case MirrorOnBounds:
					case Periodic:
						k = Degree;
						for (j = i + Degree / 2L; ((i - Degree / 2L) <= j); j--) {
							*Status = GetFoldedIndex(j, &(Index[k--]), LineLength,
								Convention);
							if (*Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
						}
						*Result = Weight[0] * LineCoeff[Index[10]] + Weight[1]
							* LineCoeff[Index[9]] + Weight[2] * LineCoeff[Index[8]]
							+ Weight[3] * LineCoeff[Index[7]] + Weight[4]
							* LineCoeff[Index[6]] + Weight[5] * LineCoeff[Index[5]]
							+ Weight[6] * LineCoeff[Index[4]] + Weight[7]
							* LineCoeff[Index[3]] + Weight[8] * LineCoeff[Index[2]]
							+ Weight[9] * LineCoeff[Index[1]] + Weight[10]
							* LineCoeff[Index[0]];
						break;
					default:
						*Status = ERROR;
						WRITE_ERROR(SplineInterpolateLine, "Unexpected boundary convention")
						break;
				}
			}
			break;
		case 11L:
			i = (long)floor(Argument);
			w = Argument - (double)i;
			Weight[0] = Bspline11(w - 6.0);
			Weight[1] = Bspline11(w - 5.0);
			Weight[2] = Bspline11(w - 4.0);
			Weight[3] = Bspline11(w - 3.0);
			Weight[4] = Bspline11(w - 2.0);
			Weight[5] = Bspline11(w - 1.0);
			Weight[6] = Bspline11(w);
			Weight[7] = Bspline11(w + 1.0);
			Weight[8] = Bspline11(w + 2.0);
			Weight[9] = Bspline11(w + 3.0);
			Weight[10] = Bspline11(w + 4.0);
			Weight[11] = Bspline11(w + 5.0);
			if ((5L <= i) && (i < (LineLength - 6L))) {
				*Result = Weight[0] * LineCoeff[i + 6L] + Weight[1] * LineCoeff[i + 5L]
					+ Weight[2] * LineCoeff[i + 4L] + Weight[3] * LineCoeff[i + 3L]
					+ Weight[4] * LineCoeff[i + 2L] + Weight[5] * LineCoeff[i + 1L]
					+ Weight[6] * LineCoeff[i] + Weight[7] * LineCoeff[i - 1L]
					+ Weight[8] * LineCoeff[i - 2L] + Weight[9] * LineCoeff[i - 3L]
					+ Weight[10] * LineCoeff[i - 4L] + Weight[11] * LineCoeff[i - 5L];
			}
			else {
				switch (Convention) {
					case AntiMirrorOnBounds:
					case FiniteDataSupport:
						k = Degree;
						for (j = i + (Degree + 1L) / 2L; ((i - (Degree - 1L) / 2L) <= j);
							j--) {
							*Status = GetFoldedValueDouble(LineCoeff, j, &(Value[k--]),
								LineLength, Convention);
							if (*Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
						}
						*Result = Weight[0] * Value[11] + Weight[1] * Value[10]
							+ Weight[2] * Value[9] + Weight[3] * Value[8]
							+ Weight[4] * Value[7] + Weight[5] * Value[6]
							+ Weight[6] * Value[5] + Weight[7] * Value[4]
							+ Weight[8] * Value[3] + Weight[9] * Value[2]
							+ Weight[10] * Value[1] + Weight[11] * Value[0];
						break;
					case MirrorOffBounds:
					case MirrorOnBounds:
					case Periodic:
						k = Degree;
						for (j = i + (Degree + 1L) / 2L; ((i - (Degree - 1L) / 2L) <= j);
							j--) {
							*Status = GetFoldedIndex(j, &(Index[k--]), LineLength,
								Convention);
							if (*Status == ERROR) {
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
						}
						*Result = Weight[0] * LineCoeff[Index[11]] + Weight[1]
							* LineCoeff[Index[10]] + Weight[2] * LineCoeff[Index[9]]
							+ Weight[3] * LineCoeff[Index[8]] + Weight[4]
							* LineCoeff[Index[7]] + Weight[5] * LineCoeff[Index[6]]
							+ Weight[6] * LineCoeff[Index[5]] + Weight[7]
							* LineCoeff[Index[4]] + Weight[8] * LineCoeff[Index[3]]
							+ Weight[9] * LineCoeff[Index[2]] + Weight[10]
							* LineCoeff[Index[1]] + Weight[11] * LineCoeff[Index[0]];
						break;
					default:
						*Status = ERROR;
						WRITE_ERROR(SplineInterpolateLine, "Unexpected boundary convention")
						break;
				}
			}
			break;
		default:
			AllocateLineDouble(&x, Degree + 1L, Status);
			if (*Status == ERROR) {
/**/			DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
				return(*Status);
			}
			if ((Degree & 1L) != 0L) {
				i = (long)floor(Argument);
				w = Argument - (double)i - (double)((Degree - 1L) / 2L);
				p = x;
				for (j = 0L; (j <= Degree); j++) {
					*Status = Bspline(Degree, w++, p++);
					if (*Status == ERROR) {
						FreeLineDouble(&x);
/**/					DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
						return(*Status);
					}
				}
				if ((((Degree - 1L) / 2L) <= i) && (i < (LineLength - (Degree + 1L)
					/ 2L))) {
					*Result = 0.0;
					p = LineCoeff + (ptrdiff_t)(i + (Degree + 1L) / 2L);
					y = x;
					for (j = 0L; (j <= Degree); j++) {
						*Result += *y++ * *p--;
					}
				}
				else {
					switch (Convention) {
						case AntiMirrorOnBounds:
						case FiniteDataSupport:
							AllocateLineDouble(&y, Degree + 1L, Status);
							if (*Status == ERROR) {
								FreeLineDouble(&x);
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
							k = Degree;
							for (j = i + (Degree + 1L) / 2L; ((i - (Degree - 1L) / 2L)
								<= j); j--) {
								*Status = GetFoldedValueDouble(LineCoeff, j, &(y[k--]),
									LineLength, Convention);
								if (*Status == ERROR) {
									FreeLineDouble(&y);
									FreeLineDouble(&x);
/**/								DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
									return(*Status);
								}
							}
							*Result = 0.0;
							k = Degree;
							for (j = 0L; (j <= Degree); j++) {
								*Result += x[j] * y[k--];
							}
							*Status = FreeLineDouble(&y);
							if (*Status == ERROR) {
								FreeLineDouble(&x);
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
							break;
						case MirrorOffBounds:
						case MirrorOnBounds:
						case Periodic:
							n = (long *)malloc((size_t)((Degree + 1L)
								* (long)sizeof(long)));
							if (n == (long *)NULL) {
								*Status = ERROR;
								FreeLineDouble(&x);
								WRITE_ERROR(SplineInterpolateLine,
									"Unable to perform allocation")
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
							k = Degree;
							for (j = i + (Degree + 1L) / 2L; ((i - (Degree - 1L) / 2L)
								<= j); j--) {
								*Status = GetFoldedIndex(j, &(n[k--]), LineLength,
									Convention);
								if (*Status == ERROR) {
									free(n);
									FreeLineDouble(&x);
/**/								DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
									return(*Status);
								}
							}
							*Result = 0.0;
							k = Degree;
							for (j = 0L; (j <= Degree); j++) {
								*Result += x[j] * LineCoeff[n[k--]];
							}
							free(n);
							break;
						default:
							*Status = ERROR;
							FreeLineDouble(&x);
							WRITE_ERROR(SplineInterpolateLine,
								"Unexpected boundary convention")
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
					}
				}
			}
			else {
				i = ConvertDoubleToLong(Argument);
				w = Argument - (double)i - (double)(Degree / 2L);
				p = x;
				for (j = 0L; (j <= Degree); j++) {
					*Status = Bspline(Degree, w++, p++);
					if (*Status == ERROR) {
						FreeLineDouble(&x);
/**/					DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
						return(*Status);
					}
				}
				if (((Degree / 2L) <= i) && (i < (LineLength - Degree / 2L))) {
					*Result = 0.0;
					p = LineCoeff + (ptrdiff_t)(i + Degree / 2L);
					y = x;
					for (j = 0L; (j <= Degree); j++) {
						*Result += *y++ * *p--;
					}
				}
				else {
					switch (Convention) {
						case AntiMirrorOnBounds:
						case FiniteDataSupport:
							AllocateLineDouble(&y, Degree + 1L, Status);
							if (*Status == ERROR) {
								FreeLineDouble(&x);
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
							k = Degree;
							for (j = i + Degree / 2L; ((i - Degree / 2L) <= j); j--) {
								*Status = GetFoldedValueDouble(LineCoeff, j, &(y[k--]),
									LineLength, Convention);
								if (*Status == ERROR) {
									FreeLineDouble(&y);
									FreeLineDouble(&x);
/**/								DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
									return(*Status);
								}
							}
							*Result = 0.0;
							k = Degree;
							for (j = 0L; (j <= Degree); j++) {
								*Result += x[j] * y[k--];
							}
							*Status = FreeLineDouble(&y);
							if (*Status == ERROR) {
								FreeLineDouble(&x);
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
							break;
						case MirrorOffBounds:
						case MirrorOnBounds:
						case Periodic:
							n = (long *)malloc((size_t)((Degree + 1L)
								* (long)sizeof(long)));
							if (n == (long *)NULL) {
								*Status = ERROR;
								FreeLineDouble(&x);
								WRITE_ERROR(SplineInterpolateLine,
									"Unable to perform allocation")
/**/							DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
								return(*Status);
							}
							k = Degree;
							for (j = i + Degree / 2L; ((i - Degree / 2L) <= j); j--) {
								*Status = GetFoldedIndex(j, &(n[k--]), LineLength,
									Convention);
								if (*Status == ERROR) {
									free(n);
									FreeLineDouble(&x);
/**/								DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
									return(*Status);
								}
							}
							*Result = 0.0;
							k = Degree;
							for (j = 0L; (j <= Degree); j++) {
								*Result += x[j] * LineCoeff[n[k--]];
							}
							free(n);
							break;
						default:
							*Status = ERROR;
							FreeLineDouble(&x);
							WRITE_ERROR(SplineInterpolateLine,
								"Unexpected boundary convention")
/**/						DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
							return(*Status);
					}
				}
			}
			*Status = FreeLineDouble(&x);
			break;
	}
/**/DEBUG_WRITE_LEAVING(SplineInterpolateLine, "Done")
	return(*Status);
} /* end SplineInterpolateLine */

/*--------------------------------------------------------------------------*/
extern int		SplineInterpolateVolume
				(
					float	*VolumeCoeff,		/* B-spline coefficients to interpolate */
					long	Nx,					/* width of the volume */
					long	Ny,					/* height of the volume */
					long	Nz,					/* depth of the volume */
					double	Xargument,			/* input X abscissa */
					double	Yargument,			/* input Y abscissa */
					double	Zargument,			/* input Z abscissa */
					double	*Result,			/* output ordinate */
					long	Degree,				/* degree of the spline */
					enum TBoundaryConvention
							Convention,			/* boundary convention */
					int		*Status				/* error management */
				)

{ /* begin SplineInterpolateVolume */

#if defined(CODEWARRIOR)
#pragma unused(VolumeCoeff, Nx, Ny, Nz, Xargument, Yargument, Zargument, Result, Degree)
#endif

	*Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(SplineInterpolateVolume, VolumeCoeff, *Status,
/**/	"Missing VolumeCoeff")
/**/DEBUG_CHECK_NULL_POINTER(SplineInterpolateVolume, Result, *Status,
/**/	"Missing Result")
/**/DEBUG_CHECK_RANGE_LONG(SplineInterpolateVolume, Nx, 1L, LONG_MAX, *Status,
/**/	"Invalid Nx (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SplineInterpolateVolume, Ny, 1L, LONG_MAX, *Status,
/**/	"Invalid Ny (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SplineInterpolateVolume, Nz, 1L, LONG_MAX, *Status,
/**/	"Invalid Nz (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_LONG(SplineInterpolateVolume, Degree, 0L, LONG_MAX, *Status,
/**/	"Invalid Degree (should be positive)")
/**/DEBUG_RETURN_ON_ERROR(SplineInterpolateVolume, *Status)
/**/DEBUG_WRITE_ENTERING(SplineInterpolateVolume,
/**/	"About to execute SplineInterpolateVolume")

	switch (Convention) {
		case AntiMirrorOnBounds:
		case MirrorOffBounds:
		case MirrorOnBounds:
		case Periodic:
			break;
		case FiniteDataSupport:
			if (2L <= Degree) {
				*Status = ERROR;
				WRITE_ERROR(SplineInterpolateVolume, "Invalid boundary convention")
/**/			DEBUG_WRITE_LEAVING(SplineInterpolateVolume, "Done")
				return(*Status);
			}
			break;
		case FiniteCoefficientSupport:
			Convention = FiniteDataSupport;
			break;
		default:
			*Status = ERROR;
			WRITE_ERROR(SplineInterpolateVolume, "Invalid boundary convention")
/**/		DEBUG_WRITE_LEAVING(SplineInterpolateVolume, "Done")
			return(*Status);
	}

	*Status = ERROR;
	WRITE_ERROR(SplineInterpolateVolume, "Not yet implemented")

/**/DEBUG_WRITE_LEAVING(SplineInterpolateVolume, "Done")
	return(*Status);
} /* SplineInterpolateVolume */

