/*****************************************************************************
 *	System includes
 ****************************************************************************/
#include	<stddef.h>
#include	<string.h>

/*****************************************************************************
 *	Toolbox defines
 ****************************************************************************/
#include	"configs.h"
#include	"debug.h"

/*****************************************************************************
 *	Toolbox includes
 ****************************************************************************/
#include	"dht.h"

/*****************************************************************************
 *	Conditional includes
 ****************************************************************************/
#ifdef		DEBUG
#include	<limits.h>
#include	<stdio.h>
#include	"messagedisplay.h"
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
extern int		DiscreteHartleyTransform
				(
					double	Data[],				/* signal */
					double	*ScratchBuffer,		/* scratch buffer */
					double	CaS[],				/* coefficients */
					long	SignalLength,		/* signal length */
					int		Forward				/* direction */
				)

/* computes the Discrete Hartley Transform of a Real signal */
/* the input signal is given by Data (double) */
/* the computation is in-place (the output replaces the input) */
/* ScratchBuffer is a pre-allocated workspace of size SignalLength */
/* the values returned in ScratchBuffer are meaningless */
/* CaS is an input array of coefficients of size SignalLength (see GetCaS function) */
/* CaS is modified internally, but is restored when DHT returns */
/* SignalLength is the length of the signal */
/* no restriction on SignalLength, but best efficiency for radix 2, 3, 4, 5 */
/* performs a DHT transform when (Forward == TRUE) */
/* performs an inverse DHT transform when (Forward == FALSE) */
/* success: return(!ERROR); failure: return(ERROR); */

{ /* begin DiscreteHartleyTransform */

	static double
			*p, *q, *r, *s, *t, *u, *v, *w, *x, *y, *z;
	static double
			*x1, *y1, *x2, *y2, *x3, *y3, *x4, *y4;
	static double
			a, b, c, d;
	long	nHalf, nThird, nFifth;
	static long
			i;
	static int
			Status = !ERROR;

/**/DEBUG_CHECK_NULL_POINTER(DiscreteHartleyTransform, Data, Status,
/**/	"No data")
/**/DEBUG_CHECK_NULL_POINTER(DiscreteHartleyTransform, ScratchBuffer, Status,
/**/	"No buffer")
/**/DEBUG_CHECK_NULL_POINTER(DiscreteHartleyTransform, CaS, Status,
/**/	"No Hartley coefficients")
/**/DEBUG_CHECK_RANGE_LONG(DiscreteHartleyTransform, SignalLength, 1L, LONG_MAX, Status,
/**/	"Invalid length (should be strictly positive)")
/**/DEBUG_CHECK_RANGE_INT(DiscreteHartleyTransform, Forward, FALSE, TRUE, Status,
/**/	"Use TRUE for direct Hartley transform and FALSE for inverse Hartley transform")
/**/DEBUG_RETURN_ON_ERROR(DiscreteHartleyTransform, Status)
/**/DEBUG_WRITE_ENTERING(DiscreteHartleyTransform,
/**/	"About to perform a Hartley transform")

	if (SignalLength == 1L) {
	}
	else if (SignalLength == 2L) {
		a = *Data++;
		b = *Data;
		if (Forward) {
			*Data-- = (a - b) * 0.5;
			*Data = (a + b) * 0.5;
		}
		else {
			*Data-- = a - b;
			*Data = a + b;
		}
	}
	else if (SignalLength == 3L) {
		a = *Data++;
		b = *Data++;
		if (Forward) {
			c = (b - *Data) * 0.2886751345948128822545743902509787278238008756350634380093012;
			b += *Data;
			d = (2.0 * a - b) * (1.0 / 6.0);
			*Data-- = d - c;
			*Data-- = d + c;
			*Data = (a + b) * (1.0 / 3.0);
		}
		else {
			c = (b - *Data) * 0.866025403784438646763723170752936183471402626905190314027903;
			b += *Data;
			d = a - b * 0.5;
			*Data-- = d - c;
			*Data-- = d + c;
			*Data = a + b;
		}
	}
	else if (SignalLength == 4L) {
		p = Data + (ptrdiff_t)2;
		a = *p;
		*p = *Data - *p;
		*Data += a;
		Data++;
		a = *++p;
		*p = *Data - *p;
		*Data += a;
		Data++;
		a = *p;
		if (Forward) {
			*p = (*Data-- - *p) * 0.25;
			b = *Data;
			*Data-- = (*--p + a) * 0.25;
			*p = (*Data - b) * 0.25;
			*Data = (*Data + b) * 0.25;
		}
		else {
			*p = *Data-- - *p;
			b = *Data;
			*Data-- = *--p + a;
			*p = *Data - b;
			*Data += b;
		}
	}
	else if (SignalLength == 5L) {
		ScratchBuffer = (double *)memcpy(ScratchBuffer, Data,
			(size_t)(5L * (long)sizeof(double)));
		p = ScratchBuffer + (ptrdiff_t)5;
		q = CaS + (ptrdiff_t)5;
		r = ScratchBuffer;
		a = 0.0;
		while (r < p) {
			a += *r++;
		}
		*Data++ = (Forward) ? (a * 0.2) : (a);
		for (i = 1L; (i < 5L); i++) {
			r = ScratchBuffer;
			s = CaS;
			a = 0.0;
			while (r < p) {
				a += *r++ * *s;
				s += (ptrdiff_t)i;
				if (q <= s) {
					s -= (ptrdiff_t)5;
				}
			}
			*Data++ = (Forward) ? (a * 0.2) : (a);
		}
	}
	else if ((SignalLength & 1L) == 0L) {
		nHalf = SignalLength / 2L;
		ScratchBuffer = (double *)memcpy(ScratchBuffer, Data,
			(size_t)(SignalLength * (long)sizeof(double)));
		p = ScratchBuffer;
		q = ScratchBuffer + (ptrdiff_t)nHalf;
		r = Data;
		s = Data + (ptrdiff_t)1;
		for (i = nHalf; (0L < i--); r += (ptrdiff_t)2, s += (ptrdiff_t)2) {
			*p++ = *r;
			*q++ = *s;
		}
		Data = (double *)memcpy(Data, CaS, (size_t)(SignalLength * (long)sizeof(double)));
		p = CaS;
		q = CaS + (ptrdiff_t)nHalf;
		r = Data;
		s = Data + (ptrdiff_t)1;
		for (i = nHalf; (0L < i--); r += (ptrdiff_t)2, s += (ptrdiff_t)2) {
			*p++ = *r;
			*q++ = *s;
		}
		Status = DiscreteHartleyTransform(ScratchBuffer, Data, CaS, nHalf, Forward);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(DiscreteHartleyTransform, "Done")
			return(Status);
		}
		Status = DiscreteHartleyTransform(ScratchBuffer + (ptrdiff_t)nHalf,
			Data, CaS, nHalf, Forward);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(DiscreteHartleyTransform, "Done")
			return(Status);
		}
		Data = (double *)memcpy(Data, CaS, (size_t)(SignalLength * (long)sizeof(double)));
		p = CaS;
		r = Data;
		s = Data + (ptrdiff_t)nHalf;
		i = nHalf;
		while (0L < i--) {
			*p++ = *r++;
			*p++ = *s++;
		}
		p = ScratchBuffer + (ptrdiff_t)nHalf;
		q = ScratchBuffer + (ptrdiff_t)SignalLength;
		for (p++, q--; (p < q); p++, q--) {
			a = *p;
			*p = (*p + *q) * 0.5;
			*q = (a - *q) * 0.5;
		}
		p = CaS;
		q = CaS + (ptrdiff_t)SignalLength;
		r = ScratchBuffer + (ptrdiff_t)nHalf;
		s = ScratchBuffer + (ptrdiff_t)SignalLength;
		t = ScratchBuffer;
		if (Forward) {
			*Data++ = (*t++ + *r++) * 0.5;
			p++;
			q--;
			s--;
			while (r < s) {
				*Data++ = (*t++ + *r++ * *p++ + *s-- * *q--) * 0.5;
			}
			if (r == s++) {
				*Data++ = (*t++ + *r * *p++) * 0.5;
				q--;
			}
			r--;
			while (t < r) {
				*Data++ = (*t++ + *r-- * *p++ - *s++ * *q--) * 0.5;
			}
			t = ScratchBuffer;
			*Data++ = (*t++ - *r++) * 0.5;
			p++;
			q--;
			s--;
			while (r < s) {
				*Data++ = (*t++ + *r++ * *p++ + *s-- * *q--) * 0.5;
			}
			if (r == s++) {
				*Data++ = (*t++ + *r * *p++) * 0.5;
				q--;
			}
			r--;
			while (t < r) {
				*Data++ = (*t++ + *r-- * *p++ - *s++ * *q--) * 0.5;
			}
		}
		else {
			*Data++ = *t++ + *r++;
			p++;
			q--;
			s--;
			while (r < s) {
				*Data++ = *t++ + *r++ * *p++ + *s-- * *q--;
			}
			if (r == s++) {
				*Data++ = *t++ + *r * *p++;
				q--;
			}
			r--;
			while (t < r) {
				*Data++ = *t++ + *r-- * *p++ - *s++ * *q--;
			}
			t = ScratchBuffer;
			*Data++ = *t++ - *r++;
			p++;
			q--;
			s--;
			while (r < s) {
				*Data++ = *t++ + *r++ * *p++ + *s-- * *q--;
			}
			if (r == s++) {
				*Data++ = *t++ + *r * *p++;
				q--;
			}
			r--;
			while (t < r) {
				*Data++ = *t++ + *r-- * *p++ - *s++ * *q--;
			}
		}
	}
	else if ((SignalLength % 3L) == 0L) {
		nThird = SignalLength / 3L;
		ScratchBuffer = (double *)memcpy(ScratchBuffer, Data,
			(size_t)(SignalLength * (long)sizeof(double)));
		p = ScratchBuffer;
		q = p + (ptrdiff_t)nThird;
		r = q + (ptrdiff_t)nThird;
		s = Data;
		t = s + (ptrdiff_t)1;
		u = t + (ptrdiff_t)1;
		for (i = nThird; (0L < i--);
			s += (ptrdiff_t)3, t += (ptrdiff_t)3, u += (ptrdiff_t)3) {
			*p++ = *s;
			*q++ = *t;
			*r++ = *u;
		}
		Data = (double *)memcpy(Data, CaS, (size_t)(SignalLength * (long)sizeof(double)));
		p = CaS;
		q = p + (ptrdiff_t)nThird;
		r = q + (ptrdiff_t)nThird;
		s = Data;
		t = s + (ptrdiff_t)1;
		u = t + (ptrdiff_t)1;
		for (i = nThird; (0L < i--);
			s += (ptrdiff_t)3, t += (ptrdiff_t)3, u += (ptrdiff_t)3) {
			*p++ = *s;
			*q++ = *t;
			*r++ = *u;
		}
		Status = DiscreteHartleyTransform(ScratchBuffer, Data, CaS, nThird, Forward);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(DiscreteHartleyTransform, "Done")
			return(Status);
		}
		Status = DiscreteHartleyTransform(ScratchBuffer + (ptrdiff_t)nThird,
			Data, CaS, nThird, Forward);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(DiscreteHartleyTransform, "Done")
			return(Status);
		}
		Status = DiscreteHartleyTransform(ScratchBuffer + (ptrdiff_t)(2L * nThird),
			Data, CaS, nThird, Forward);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(DiscreteHartleyTransform, "Done")
			return(Status);
		}
		Data = (double *)memcpy(Data, CaS, (size_t)(SignalLength * (long)sizeof(double)));
		p = CaS;
		r = Data;
		s = r + (ptrdiff_t)nThird;
		t = s + (ptrdiff_t)nThird;
		i = nThird;
		while (0L < i--) {
			*p++ = *r++;
			*p++ = *s++;
			*p++ = *t++;
		}
		p = ScratchBuffer + (ptrdiff_t)nThird;
		q = p + (ptrdiff_t)nThird;
		r = q;
		s = r + (ptrdiff_t)nThird;
		for (p++, q--, r++, s--; (p < q); p++, q--, r++, s--) {
			a = *p;
			*p = (*p + *q) * 0.5;
			*q = (a - *q) * 0.5;
			a = *r;
			*r = (*r + *s) * 0.5;
			*s = (a - *s) * 0.5;
		}
		p = CaS;
		q = p + (ptrdiff_t)SignalLength;
		x = CaS;
		y = x + (ptrdiff_t)SignalLength;
		t = ScratchBuffer;
		r = ScratchBuffer + (ptrdiff_t)nThird;
		s = r + (ptrdiff_t)nThird;
		u = s;
		v = u + (ptrdiff_t)nThird;
		if (Forward) {
			*Data++ = (*t++ + *r++ + *u++) * (1.0 / 3.0);
			for (p++, q--, x += (ptrdiff_t)2, y -= (ptrdiff_t)2, s--, v--;
				(r < s); x += (ptrdiff_t)2, y -= (ptrdiff_t)2) {
				*Data++ = (*t++ + *r++ * *p++ + *s-- * *q-- + *u++ * *x + *v-- * *y)
					* (1.0 / 3.0);
			}
			for (r--, s++, u--, v++; (t < r); x += (ptrdiff_t)2, y -= (ptrdiff_t)2) {
				*Data++ = (*t++ + *r-- * *p++ - *s++ * *q-- + *u-- * *x - *v++ * *y)
					* (1.0 / 3.0);
			}
			t = ScratchBuffer;
			*Data++ = (*t++ + *r++ * *p++ + *u++ * *x) * (1.0 / 3.0);
			x += (ptrdiff_t)2;
			y -= (ptrdiff_t)2;
			for (q--, s--, v--; (r < s); x += (ptrdiff_t)2, y -= (ptrdiff_t)2) {
				*Data++ = (*t++ + *r++ * *p++ + *s-- * *q-- + *u++ * *x + *v-- * *y)
					* (1.0 / 3.0);
			}
			x -= (ptrdiff_t)2;
			y += (ptrdiff_t)2;
			for (r--, s++, u--, v++; (t < r); x -= (ptrdiff_t)2, y += (ptrdiff_t)2) {
				*Data++ = (*t++ + *r-- * *p++ - *s++ * *q-- + *u-- * *y - *v++ * *x)
					* (1.0 / 3.0);
			}
			t = ScratchBuffer;
			*Data++ = (*t++ + *r++ * *p++ + *u++ * *y) * (1.0 / 3.0);
			x -= (ptrdiff_t)2;
			y += (ptrdiff_t)2;
			for (q--, s--, v--; (r < s); x -= (ptrdiff_t)2, y += (ptrdiff_t)2) {
				*Data++ = (*t++ + *r++ * *p++ + *s-- * *q-- + *u++ * *y + *v-- * *x)
					* (1.0 / 3.0);
			}
			for (r--, s++, u--, v++; (t < r); x -= (ptrdiff_t)2, y += (ptrdiff_t)2) {
				*Data++ = (*t++ + *r-- * *p++ - *s++ * *q-- + *u-- * *y - *v++ * *x)
					* (1.0 / 3.0);
			}
		}
		else {
			*Data++ = *t++ + *r++ + *u++;
			for (p++, q--, x += (ptrdiff_t)2, y -= (ptrdiff_t)2, s--, v--;
				(r < s); x += (ptrdiff_t)2, y -= (ptrdiff_t)2) {
				*Data++ = *t++ + *r++ * *p++ + *s-- * *q-- + *u++ * *x + *v-- * *y;
			}
			for (r--, s++, u--, v++; (t < r); x += (ptrdiff_t)2, y -= (ptrdiff_t)2) {
				*Data++ = *t++ + *r-- * *p++ - *s++ * *q-- + *u-- * *x - *v++ * *y;
			}
			t = ScratchBuffer;
			*Data++ = *t++ + *r++ * *p++ + *u++ * *x;
			x += (ptrdiff_t)2;
			y -= (ptrdiff_t)2;
			for (q--, s--, v--; (r < s); x += (ptrdiff_t)2, y -= (ptrdiff_t)2) {
				*Data++ = *t++ + *r++ * *p++ + *s-- * *q-- + *u++ * *x + *v-- * *y;
			}
			x -= (ptrdiff_t)2;
			y += (ptrdiff_t)2;
			for (r--, s++, u--, v++; (t < r); x -= (ptrdiff_t)2, y += (ptrdiff_t)2) {
				*Data++ = *t++ + *r-- * *p++ - *s++ * *q-- + *u-- * *y - *v++ * *x;
			}
			t = ScratchBuffer;
			*Data++ = *t++ + *r++ * *p++ + *u++ * *y;
			x -= (ptrdiff_t)2;
			y += (ptrdiff_t)2;
			for (q--, s--, v--; (r < s); x -= (ptrdiff_t)2, y += (ptrdiff_t)2) {
				*Data++ = *t++ + *r++ * *p++ + *s-- * *q-- + *u++ * *y + *v-- * *x;
			}
			for (r--, s++, u--, v++; (t < r); x -= (ptrdiff_t)2, y += (ptrdiff_t)2) {
				*Data++ = *t++ + *r-- * *p++ - *s++ * *q-- + *u-- * *y - *v++ * *x;
			}
		}
	}
	else if ((SignalLength % 5L) == 0L) {
		nFifth = SignalLength / 5L;
		ScratchBuffer = (double *)memcpy(ScratchBuffer, Data,
			(size_t)(SignalLength * (long)sizeof(double)));
		p = ScratchBuffer;
		q = p + (ptrdiff_t)nFifth;
		r = q + (ptrdiff_t)nFifth;
		s = r + (ptrdiff_t)nFifth;
		t = s + (ptrdiff_t)nFifth;
		u = Data;
		v = u + (ptrdiff_t)1;
		w = v + (ptrdiff_t)1;
		x = w + (ptrdiff_t)1;
		y = x + (ptrdiff_t)1;
		for (i = nFifth; (0L < i--); u += (ptrdiff_t)5, v += (ptrdiff_t)5,
			w += (ptrdiff_t)5, x += (ptrdiff_t)5, y += (ptrdiff_t)5) {
			*p++ = *u;
			*q++ = *v;
			*r++ = *w;
			*s++ = *x;
			*t++ = *y;
		}
		Data = (double *)memcpy(Data, CaS, (size_t)(SignalLength * (long)sizeof(double)));
		p = CaS;
		q = p + (ptrdiff_t)nFifth;
		r = q + (ptrdiff_t)nFifth;
		s = r + (ptrdiff_t)nFifth;
		t = s + (ptrdiff_t)nFifth;
		u = Data;
		v = u + (ptrdiff_t)1;
		w = v + (ptrdiff_t)1;
		x = w + (ptrdiff_t)1;
		y = x + (ptrdiff_t)1;
		for (i = nFifth; (0L < i--); u += (ptrdiff_t)5, v += (ptrdiff_t)5,
			w += (ptrdiff_t)5, x += (ptrdiff_t)5, y += (ptrdiff_t)5) {
			*p++ = *u;
			*q++ = *v;
			*r++ = *w;
			*s++ = *x;
			*t++ = *y;
		}
		Status = DiscreteHartleyTransform(ScratchBuffer, Data, CaS, nFifth, Forward);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(DiscreteHartleyTransform, "Done")
			return(Status);
		}
		Status = DiscreteHartleyTransform(ScratchBuffer + (ptrdiff_t)nFifth,
			Data, CaS, nFifth, Forward);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(DiscreteHartleyTransform, "Done")
			return(Status);
		}
		Status = DiscreteHartleyTransform(ScratchBuffer + (ptrdiff_t)(2L * nFifth),
			Data, CaS, nFifth, Forward);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(DiscreteHartleyTransform, "Done")
			return(Status);
		}
		Status = DiscreteHartleyTransform(ScratchBuffer + (ptrdiff_t)(3L * nFifth),
			Data, CaS, nFifth, Forward);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(DiscreteHartleyTransform, "Done")
			return(Status);
		}
		Status = DiscreteHartleyTransform(ScratchBuffer + (ptrdiff_t)(4L * nFifth),
			Data, CaS, nFifth, Forward);
		if (Status == ERROR) {
/**/		DEBUG_WRITE_LEAVING(DiscreteHartleyTransform, "Done")
			return(Status);
		}
		Data = (double *)memcpy(Data, CaS, (size_t)(SignalLength * (long)sizeof(double)));
		p = CaS;
		u = Data;
		v = u + (ptrdiff_t)nFifth;
		w = v + (ptrdiff_t)nFifth;
		x = w + (ptrdiff_t)nFifth;
		y = x + (ptrdiff_t)nFifth;
		i = nFifth;
		while (0L < i--) {
			*p++ = *u++;
			*p++ = *v++;
			*p++ = *w++;
			*p++ = *x++;
			*p++ = *y++;
		}
		p = ScratchBuffer + (ptrdiff_t)nFifth;
		q = p + (ptrdiff_t)nFifth;
		r = q;
		s = r + (ptrdiff_t)nFifth;
		t = s;
		u = t + (ptrdiff_t)nFifth;
		v = u;
		w = v + (ptrdiff_t)nFifth;
		for (p++, q--, r++, s--, t++, u--, v++, w--; (p < q);
			p++, q--, r++, s--, t++, u--, v++, w--) {
			a = *p;
			*p = (*p + *q) * 0.5;
			*q = (a - *q) * 0.5;
			a = *r;
			*r = (*r + *s) * 0.5;
			*s = (a - *s) * 0.5;
			a = *t;
			*t = (*t + *u) * 0.5;
			*u = (a - *u) * 0.5;
			a = *v;
			*v = (*v + *w) * 0.5;
			*w = (a - *w) * 0.5;
		}
		x1 = CaS;
		y1 = x1 + (ptrdiff_t)SignalLength;
		x2 = CaS;
		y2 = x2 + (ptrdiff_t)SignalLength;
		x3 = CaS;
		y3 = x3 + (ptrdiff_t)SignalLength;
		x4 = CaS;
		y4 = x4 + (ptrdiff_t)SignalLength;
		z = ScratchBuffer;
		p = z + (ptrdiff_t)nFifth;
		q = p + (ptrdiff_t)nFifth;
		r = q;
		s = r + (ptrdiff_t)nFifth;
		t = s;
		u = t + (ptrdiff_t)nFifth;
		v = u;
		w = v + (ptrdiff_t)nFifth;
		if (Forward) {
			*Data++ = (*z++ + *p++ + *r++ + *t++ + *v++) * 0.2;
			x2 += (ptrdiff_t)2;
			y2 -= (ptrdiff_t)2;
			x3 += (ptrdiff_t)3;
			y3 -= (ptrdiff_t)3;
			x4 += (ptrdiff_t)4;
			y4 -= (ptrdiff_t)4;
			for (x1++, y1--, q--, s--, u--, w--; (p < q);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				*Data++ = (*z++ + *p++ * *x1++ + *q-- * *y1-- + *r++ * *x2 + *s-- * *y2
					+ *t++ * *x3 + *u-- * *y3 + *v++ * *x4 + *w-- * *y4) * 0.2;
			}
			for (p--, q++, r--, s++, t--, u++, v--, w++; (z < p);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				*Data++ = (*z++ + *p-- * *x1++ - *q++ * *y1-- + *r-- * *x2 - *s++ * *y2
					+ *t-- * *x3 - *u++ * *y3 + *v-- * *x4 - *w++ * *y4) * 0.2;
			}
			z = ScratchBuffer;
			*Data++ = (*z++ + *p++ * *x1 + *r++ * *x2 + *t++ * *x3 + *v++ * *x4) * 0.2;
			x2 += (ptrdiff_t)2;
			y2 -= (ptrdiff_t)2;
			x3 += (ptrdiff_t)3;
			y3 -= (ptrdiff_t)3;
			x4 += (ptrdiff_t)4;
			y4 -= (ptrdiff_t)4;
			for (x1++, y1--, q--, s--, u--, w--; (p < q);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				if (y4 < CaS) {
					x4 -= (ptrdiff_t)SignalLength;
					y4 += (ptrdiff_t)SignalLength;
				}
				*Data++ = (*z++ + *p++ * *x1++ + *q-- * *y1-- + *r++ * *x2 + *s-- * *y2
					+ *t++ * *x3 + *u-- * *y3 + *v++ * *x4 + *w-- * *y4) * 0.2;
			}
			for (p--, q++, r--, s++, t--, u++, v--, w++; (z < p);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				if (y3 < CaS) {
					x3 -= (ptrdiff_t)SignalLength;
					y3 += (ptrdiff_t)SignalLength;
				}
				*Data++ = (*z++ + *p-- * *x1++ - *q++ * *y1-- + *r-- * *x2 - *s++ * *y2
					+ *t-- * *x3 - *u++ * *y3 + *v-- * *x4 - *w++ * *y4) * 0.2;
			}
			z = ScratchBuffer;
			*Data++ = (*z++ + *p++ * *x1 + *r++ * *x2 + *t++ * *x3 + *v++ * *x4) * 0.2;
			x2 += (ptrdiff_t)2;
			y2 -= (ptrdiff_t)2;
			x3 += (ptrdiff_t)3;
			y3 -= (ptrdiff_t)3;
			x4 += (ptrdiff_t)4;
			y4 -= (ptrdiff_t)4;
			for (x1++, y1--, q--, s--, u--, w--; (p < q);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				*Data++ = (*z++ + *p++ * *x1++ + *q-- * *y1-- + *r++ * *x2 + *s-- * *y2
					+ *t++ * *x3 + *u-- * *y3 + *v++ * *x4 + *w-- * *y4) * 0.2;
			}
			x2 -= (ptrdiff_t)SignalLength;
			y2 += (ptrdiff_t)SignalLength;
			x4 -= (ptrdiff_t)SignalLength;
			y4 += (ptrdiff_t)SignalLength;
			for (p--, q++, r--, s++, t--, u++, v--, w++; (z < p);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				*Data++ = (*z++ + *p-- * *x1++ - *q++ * *y1-- + *r-- * *x2 - *s++ * *y2
					+ *t-- * *x3 - *u++ * *y3 + *v-- * *x4 - *w++ * *y4) * 0.2;
			}
			z = ScratchBuffer;
			*Data++ = (*z++ + *p++ * *x1 + *r++ * *x2 + *t++ * *x3 + *v++ * *x4) * 0.2;
			x2 += (ptrdiff_t)2;
			y2 -= (ptrdiff_t)2;
			x3 += (ptrdiff_t)3;
			y3 -= (ptrdiff_t)3;
			x4 += (ptrdiff_t)4;
			y4 -= (ptrdiff_t)4;
			for (x1++, y1--, q--, s--, u--, w--; (p < q);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				if (y3 < CaS) {
					x3 -= (ptrdiff_t)SignalLength;
					y3 += (ptrdiff_t)SignalLength;
				}
				*Data++ = (*z++ + *p++ * *x1++ + *q-- * *y1-- + *r++ * *x2 + *s-- * *y2
					+ *t++ * *x3 + *u-- * *y3 + *v++ * *x4 + *w-- * *y4) * 0.2;
			}
			for (p--, q++, r--, s++, t--, u++, v--, w++; (z < p);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				if (y4 < CaS) {
					x4 -= (ptrdiff_t)SignalLength;
					y4 += (ptrdiff_t)SignalLength;
				}
				*Data++ = (*z++ + *p-- * *x1++ - *q++ * *y1-- + *r-- * *x2 - *s++ * *y2
					+ *t-- * *x3 - *u++ * *y3 + *v-- * *x4 - *w++ * *y4) * 0.2;
			}
			z = ScratchBuffer;
			*Data++ = (*z++ + *p++ * *x1 + *r++ * *x2 + *t++ * *x3 + *v++ * *x4) * 0.2;
			x2 += (ptrdiff_t)2;
			y2 -= (ptrdiff_t)2;
			x3 += (ptrdiff_t)3;
			y3 -= (ptrdiff_t)3;
			x4 += (ptrdiff_t)4;
			y4 -= (ptrdiff_t)4;
			for (x1++, y1--, q--, s--, u--, w--; (p < q);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				*Data++ = (*z++ + *p++ * *x1++ + *q-- * *y1-- + *r++ * *x2 + *s-- * *y2
					+ *t++ * *x3 + *u-- * *y3 + *v++ * *x4 + *w-- * *y4) * 0.2;
			}
			for (p--, q++, r--, s++, t--, u++, v--, w++; (z < p);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				*Data++ = (*z++ + *p-- * *x1++ - *q++ * *y1-- + *r-- * *x2 - *s++ * *y2
					+ *t-- * *x3 - *u++ * *y3 + *v-- * *x4 - *w++ * *y4) * 0.2;
			}
		}
		else {
			*Data++ = *z++ + *p++ + *r++ + *t++ + *v++;
			x2 += (ptrdiff_t)2;
			y2 -= (ptrdiff_t)2;
			x3 += (ptrdiff_t)3;
			y3 -= (ptrdiff_t)3;
			x4 += (ptrdiff_t)4;
			y4 -= (ptrdiff_t)4;
			for (x1++, y1--, q--, s--, u--, w--; (p < q);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				*Data++ = *z++ + *p++ * *x1++ + *q-- * *y1-- + *r++ * *x2 + *s-- * *y2
					+ *t++ * *x3 + *u-- * *y3 + *v++ * *x4 + *w-- * *y4;
			}
			for (p--, q++, r--, s++, t--, u++, v--, w++; (z < p);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				*Data++ = *z++ + *p-- * *x1++ - *q++ * *y1-- + *r-- * *x2 - *s++ * *y2
					+ *t-- * *x3 - *u++ * *y3 + *v-- * *x4 - *w++ * *y4;
			}
			z = ScratchBuffer;
			*Data++ = *z++ + *p++ * *x1 + *r++ * *x2 + *t++ * *x3 + *v++ * *x4;
			x2 += (ptrdiff_t)2;
			y2 -= (ptrdiff_t)2;
			x3 += (ptrdiff_t)3;
			y3 -= (ptrdiff_t)3;
			x4 += (ptrdiff_t)4;
			y4 -= (ptrdiff_t)4;
			for (x1++, y1--, q--, s--, u--, w--; (p < q);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				if (y4 < CaS) {
					x4 -= (ptrdiff_t)SignalLength;
					y4 += (ptrdiff_t)SignalLength;
				}
				*Data++ = *z++ + *p++ * *x1++ + *q-- * *y1-- + *r++ * *x2 + *s-- * *y2
					+ *t++ * *x3 + *u-- * *y3 + *v++ * *x4 + *w-- * *y4;
			}
			for (p--, q++, r--, s++, t--, u++, v--, w++; (z < p);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				if (y3 < CaS) {
					x3 -= (ptrdiff_t)SignalLength;
					y3 += (ptrdiff_t)SignalLength;
				}
				*Data++ = *z++ + *p-- * *x1++ - *q++ * *y1-- + *r-- * *x2 - *s++ * *y2
					+ *t-- * *x3 - *u++ * *y3 + *v-- * *x4 - *w++ * *y4;
			}
			z = ScratchBuffer;
			*Data++ = *z++ + *p++ * *x1 + *r++ * *x2 + *t++ * *x3 + *v++ * *x4;
			x2 += (ptrdiff_t)2;
			y2 -= (ptrdiff_t)2;
			x3 += (ptrdiff_t)3;
			y3 -= (ptrdiff_t)3;
			x4 += (ptrdiff_t)4;
			y4 -= (ptrdiff_t)4;
			for (x1++, y1--, q--, s--, u--, w--; (p < q);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				*Data++ = *z++ + *p++ * *x1++ + *q-- * *y1-- + *r++ * *x2 + *s-- * *y2
					+ *t++ * *x3 + *u-- * *y3 + *v++ * *x4 + *w-- * *y4;
			}
			x2 -= (ptrdiff_t)SignalLength;
			y2 += (ptrdiff_t)SignalLength;
			x4 -= (ptrdiff_t)SignalLength;
			y4 += (ptrdiff_t)SignalLength;
			for (p--, q++, r--, s++, t--, u++, v--, w++; (z < p);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				*Data++ = *z++ + *p-- * *x1++ - *q++ * *y1-- + *r-- * *x2 - *s++ * *y2
					+ *t-- * *x3 - *u++ * *y3 + *v-- * *x4 - *w++ * *y4;
			}
			z = ScratchBuffer;
			*Data++ = *z++ + *p++ * *x1 + *r++ * *x2 + *t++ * *x3 + *v++ * *x4;
			x2 += (ptrdiff_t)2;
			y2 -= (ptrdiff_t)2;
			x3 += (ptrdiff_t)3;
			y3 -= (ptrdiff_t)3;
			x4 += (ptrdiff_t)4;
			y4 -= (ptrdiff_t)4;
			for (x1++, y1--, q--, s--, u--, w--; (p < q);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				if (y3 < CaS) {
					x3 -= (ptrdiff_t)SignalLength;
					y3 += (ptrdiff_t)SignalLength;
				}
				*Data++ = *z++ + *p++ * *x1++ + *q-- * *y1-- + *r++ * *x2 + *s-- * *y2
					+ *t++ * *x3 + *u-- * *y3 + *v++ * *x4 + *w-- * *y4;
			}
			for (p--, q++, r--, s++, t--, u++, v--, w++; (z < p);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				if (y4 < CaS) {
					x4 -= (ptrdiff_t)SignalLength;
					y4 += (ptrdiff_t)SignalLength;
				}
				*Data++ = *z++ + *p-- * *x1++ - *q++ * *y1-- + *r-- * *x2 - *s++ * *y2
					+ *t-- * *x3 - *u++ * *y3 + *v-- * *x4 - *w++ * *y4;
			}
			z = ScratchBuffer;
			*Data++ = *z++ + *p++ * *x1 + *r++ * *x2 + *t++ * *x3 + *v++ * *x4;
			x2 += (ptrdiff_t)2;
			y2 -= (ptrdiff_t)2;
			x3 += (ptrdiff_t)3;
			y3 -= (ptrdiff_t)3;
			x4 += (ptrdiff_t)4;
			y4 -= (ptrdiff_t)4;
			for (x1++, y1--, q--, s--, u--, w--; (p < q);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				*Data++ = *z++ + *p++ * *x1++ + *q-- * *y1-- + *r++ * *x2 + *s-- * *y2
					+ *t++ * *x3 + *u-- * *y3 + *v++ * *x4 + *w-- * *y4;
			}
			for (p--, q++, r--, s++, t--, u++, v--, w++; (z < p);
				x2 += (ptrdiff_t)2, y2 -= (ptrdiff_t)2,
				x3 += (ptrdiff_t)3, y3 -= (ptrdiff_t)3,
				x4 += (ptrdiff_t)4, y4 -= (ptrdiff_t)4) {
				*Data++ = *z++ + *p-- * *x1++ - *q++ * *y1-- + *r-- * *x2 - *s++ * *y2
					+ *t-- * *x3 - *u++ * *y3 + *v-- * *x4 - *w++ * *y4;
			}
		}
	}
	else {
		b = 1.0 / (double)SignalLength;
		ScratchBuffer = (double *)memcpy(ScratchBuffer, Data,
			(size_t)(SignalLength * (long)sizeof(double)));
		p = ScratchBuffer + (ptrdiff_t)SignalLength;
		q = CaS + (ptrdiff_t)SignalLength;
		r = ScratchBuffer;
		a = 0.0;
		while (r < p) {
			a += *r++;
		}
		*Data++ = (Forward) ? (a * b) : (a);
		for (i = 1L; (i < SignalLength); i++) {
			r = ScratchBuffer;
			s = CaS;
			a = 0.0;
			while (r < p) {
				a += *r++ * *s;
				s += (ptrdiff_t)i;
				if (q <= s) {
					s -= (ptrdiff_t)SignalLength;
				}
			}
			*Data++ = (Forward) ? (a * b) : (a);
		}
	}
/**/DEBUG_WRITE_LEAVING(DiscreteHartleyTransform, "Done")
	return(Status);
} /* end DiscreteHartleyTransform */

