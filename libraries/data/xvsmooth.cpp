/*
 * xvsmooth.c - smoothing/color dither routines for XV
 *
 *  Author:    John Bradley, University of Pennsylvania
 *                (bradley@cis.upenn.edu)
 *
 *  Contains:
 *            byte *SmoothResize(src8, swide, shigh, dwide, dhigh,
 *                               rmap, gmap, bmap, rdmap, gdmap, bdmap, maplen)
 *            byte *Smooth24(pic824, swide, shigh, dwide, dhigh,
 *                               rmap, gmap, bmap)
 *            byte *DoColorDither(picSrc, pic8, w, h, rmap,gmap,bmap,
 *                                rdisp, gdisp, bdisp, maplen)


 */

/* Copyright Notice
 * ================
 * Copyright 1989, 1990, 1991, 1992, 1993 by John Bradley
 *
 * Permission to use, copy, and distribute XV in its entirety, for
 * non-commercial purposes, is hereby granted without fee, provided that
 * this license information and copyright notice appear in all copies.
 *
 * Note that distributing XV 'bundled' in with ANY product is considered
 * to be a 'commercial purpose'.
 *
 * Also note that any copies of XV that are distributed MUST be built
 * and/or configured to be in their 'unregistered copy' mode, so that it
 * is made obvious to the user that XV is shareware, and that they should
 * consider donating, or at least reading this License Info.
 *
 * The software may be modified for your own purposes, but modified
 * versions may NOT be distributed without prior consent of the author.
 *
 * This software is provided 'as-is', without any express or implied
 * warranty.  In no event will the author be held liable for any damages
 * arising from the use of this software.
 *
 * If you would like to do something with XV that this copyright
 * prohibits (such as distributing it with a commercial product,
 * using portions of the source in some other program, etc.), please
 * contact the author (preferably via email).  Arrangements can
 * probably be worked out.
 *
 * XV is shareware for PERSONAL USE only.  You may use XV for your own
 * amusement, and if you find it nifty, useful, generally cool, or of
 * some value to you, your non-deductable donation would be greatly
 * appreciated.  $25 is the suggested donation, though, of course,
 * larger donations are quite welcome.  Folks who donate $25 or more
 * can receive a Real Nice bound copy of the XV manual for no extra
 * charge.
 *
 * Commercial, government, and institutional users MUST register their
 * copies of XV, for the exceedingly REASONABLE price of just $25 per
 * workstation/X terminal.  Site licenses are available for those who
 * wish to run XV on a large number of machines.  Contact the author
 * for more details.
 *
 * The author may be contacted via:
 *    US Mail:  John Bradley
 *              1053 Floyd Terrace
 *              Bryn Mawr, PA  19010
 *
 *    Phone:    (215) 898-8813
 *    EMail:    bradley@cis.upenn.edu
 */

/* CO: --------------------------------------------------------------------- */
/* RANGE forces a to be in the range b..c (inclusive) */
#define RANGE(a,b,c) { if (a < b) a = b;  if (a > c) a = c; }
typedef unsigned char byte;

/* CO: --------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "xmipp_error.h"

byte *Smooth(byte *picSrc8, size_t swide, size_t shigh, size_t dwide, size_t dhigh);
void DoColorDither(byte *picSmooth, byte *&picDithered, size_t w, size_t h);

#ifdef __STDC__
int SmoothXY(byte *, byte *, int, int, int, int);
#else
int SmoothXY();
#endif

#define xvbzero(s,l) memset(s, 0, l)
/***************************************************/
void SmoothResize(byte *picSrc8, byte *destpic8, size_t swide, size_t shigh, size_t dwide, size_t dhigh) {
	/* generic interface to Smooth and ColorDither code.
	 given an 8-bit-per, swide * shigh image with colormap rmap,gmap,bmap,
	 will generate a new 8-bit-per, dwide * dhigh image, which is dithered
	 using colors found in rdmap, gdmap, bdmap arrays */

	/* returns ptr to a dwide*dhigh array of bytes, or NULL on failure */
	byte * picSmooth = Smooth(picSrc8, swide, shigh, dwide, dhigh);
	if (picSmooth) {
		DoColorDither(picSmooth, destpic8, dwide, dhigh);
		free(picSmooth);
	}
}

/***************************************************/
byte *Smooth(byte *picSrc8, size_t swide, size_t shigh, size_t dwide, size_t dhigh) {
	/* does a SMOOTH resize from pic824 (which is either a swide*shigh, 8-bit
	 pic, with colormap rmap,gmap,bmap OR a swide*shigh, 24-bit image, based
	 on whether 'is24' is set) into a dwide * dhigh 24-bit image

	 returns a dwide*dhigh 24bit image, or NULL on failure (malloc) */
	/* rmap,gmap,bmap should be 'desired' colors */

	byte *picSmooth, *ptrPicSmooth;
	int *cxtab, *pxtab;
	size_t y1Off, cyOff;
	size_t ex, ey, cx, cy, px, py, apx, apy, apx_100, apy_100, x1, y1;
	size_t cA, cB, cC, cD;
	size_t pA, pB, pC, pD;
	int retval;

	cA = cB = cC = cD = 0;
	size_t picSmoothSize = ((size_t) dwide) * dhigh;
	ptrPicSmooth = picSmooth = (byte *) malloc(picSmoothSize);
	if (!picSmooth)
		REPORT_ERROR(ERR_MEM_NOTENOUGH,
				formatString("Unable to alloc: %lu",picSmoothSize));

	/* decide which smoothing routine to use based on type of expansion */
	if (dwide < swide && dhigh < shigh)
		retval = SmoothXY(picSmooth, picSrc8, swide, shigh, dwide, dhigh);
	else {
		/* dwide >= swide && dhigh >= shigh */

		/* cx,cy = original pixel in pic824.  px,py = relative position
		 of pixel ex,ey inside of cx,cy as percentages +-50%, +-50%.
		 0,0 = middle of pixel */

		/* we can save a lot of time by precomputing cxtab[] and pxtab[], both
		 dwide arrays of ints that contain values for the equations:
		 cx = (ex * swide) / dwide;
		 px = ((ex * swide * 100) / dwide) - (cx * 100) - 50; */

		cxtab = (int *) malloc(dwide * sizeof(int));
		if (!cxtab)
			REPORT_ERROR(ERR_MEM_NOTENOUGH, "Unable to alloc for smoothing");

		pxtab = (int *) malloc(dwide * sizeof(int));
		if (!pxtab)
			REPORT_ERROR(ERR_MEM_NOTENOUGH, "Unable to alloc for smoothing");

		for (ex = 0; ex < dwide; ex++) {
			cxtab[ex] = (ex * swide) / dwide;
			pxtab[ex] = (((ex * swide) * 100) / dwide) - (cxtab[ex] * 100) - 50;
		}

		for (ey = 0; ey < dhigh; ey++) {
			cy = (ey * shigh) / dhigh;
			py = (((ey * shigh) * 100) / dhigh) - (cy * 100) - 50;
			if (py < 0) {
				y1 = cy - 1;
				if (y1 < 0)
					y1 = 0;
			} else {
				y1 = cy + 1;
				if (y1 > shigh - 1)
					y1 = shigh - 1;
			}
			apy = abs(py);
			apy_100 = 100 - apy;

			cyOff = (size_t) cy * swide; /* current line */
			y1Off = (size_t) y1 * swide; /* up or down one line, depending */

			/*      if ((ey&15) == 0) WaitCursor(); */

			for (ex = 0; ex < dwide; ex++) {
				cx = cxtab[ex];
				px = pxtab[ex];

				if (px < 0) {
					x1 = cx - 1;
					if (x1 < 0)
						x1 = 0;
				} else {
					x1 = cx + 1;
					if (x1 > swide - 1)
						x1 = swide - 1;
				}

				cA = picSrc8[y1Off + x1]; /* corner pixel */
				cB = picSrc8[y1Off + cx]; /* up/down center pixel */
				cC = picSrc8[cyOff + x1]; /* left/right center pixel */
				cD = picSrc8[cyOff + cx]; /* center pixel */

				/* quick check */
				if (cA == cB && cB == cC && cC == cD)
					/* set this pixel to the same color as in pic8 */
					*ptrPicSmooth++ = cD;
				else {
					/* compute weighting factors */
					apx = abs(px);
					apx_100 = 100 - apx;
					pA = apx * apy;
					pB = apy * apx_100;
					pC = apx * apy_100;
					pD = apx_100*apy_100;

					byte val = (byte) (((pA * cA) + (pB * cB) + (pC * cC) + (pD * cD)) / 10000);
					*ptrPicSmooth++ = val;
				}
			}
		}

		free(cxtab);
		free(pxtab);
		retval = 0; /* okay */
	}

	if (retval) { /* one of the Smooth**() methods failed */
		free(picSmooth);
		picSmooth = (byte *) NULL;
	}

	return picSmooth;
}

/***************************************************/
int SmoothXY(byte *picSmooth, byte *picSrc8, int swide, int shigh, int dwide,
		int dhigh) {
	byte *cptr;
	int i, j;
	int *lbuf;
	int pix;
	int lastline, thisline, lastpix, linecnt, pixcnt;
	int *pixarr, *paptr;

	/* returns '0' if okay, '1' if failed (malloc) */

	/* shrinks pic8 into a dwide * dhigh 24-bit picture.  Only works correctly
	 when swide>=dwide and shigh>=dhigh (ie, the picture is shrunk on both
	 axes) */

	/* malloc some arrays */
	lbuf = (int *) calloc(swide, sizeof(int));
	pixarr = (int *) calloc(swide + 1, sizeof(int));
	if (!lbuf || !pixarr)
		REPORT_ERROR(ERR_MEM_NOTENOUGH, "Cannot allocate memory for smoothing");

	for (j = 0; j <= swide; j++)
		pixarr[j] = (j * dwide + (15 * swide) / 16) / swide;

	lastline = linecnt = pix = 0;
	cptr = picSrc8;

	for (i = 0; i <= shigh; i++) {
		thisline = (i * dhigh + (15 * shigh) / 16) / shigh;

		if ((thisline != lastline)) {
			pix = pixcnt = lastpix = 0;

			for (j = 0, paptr = pixarr; j <= swide; j++, paptr++) {
				if (*paptr != lastpix) { /* write a pixel to pic24 */
					*picSmooth++ = (pix / linecnt) / pixcnt;
					lastpix = *paptr;
					pix = pixcnt = 0;
				}

				if (j < swide) {
					pix += lbuf[j];
					pixcnt++;
				}
			}

			lastline = thisline;
			xvbzero((char *) lbuf, swide * sizeof(int));
			/* clear out line bufs */
			linecnt = 0;
		}

		if (i < shigh) {
			for (j = 0; j < swide; j++, cptr++)
				lbuf[j] += *cptr;

			linecnt++;
		}
	}

	free(lbuf);
	free(pixarr);
	return 0;
}

/********************************************/
void DoColorDither(byte *picSmooth, byte *&picDithered, size_t w, size_t h) {
	/* takes a 24 bit picture, of size w*h, dithers with the colors in
	 rdisp, gdisp, bdisp (which have already been allocated),
	 and generates an 8-bit w*h image, which it returns.
	 ignores input value 'pic8'
	 returns NULL on error

	 note: the rdisp,gdisp,bdisp arrays should be the 'displayed' colors,
	 not the 'desired' colors

	 if picSrc is NULL, uses the passed-in pic8 (an 8-bit image) as
	 the source, and the rmap,gmap,bmap arrays as the desired colors */

	byte *np, *ep;
	short *cache;
	int r2;
	int *thisline, *nextline, *thisptr, *nextptr, *tmpptr;
	size_t i, j;
	int rerr;
	size_t imax, jmax;
	int key;
	long cnt1, cnt2;

	cnt1 = cnt2 = 0;
	imax = h - 1;
	jmax = w - 1;
	ep = picSmooth;

	/* attempt to malloc things */
	cache = (short *) calloc(2 << 14, sizeof(short));
	thisline = (int *) malloc(w * sizeof(int));
	nextline = (int *) malloc(w * sizeof(int));
	if (!cache || !picDithered || !thisline || !nextline)
		REPORT_ERROR(ERR_MEM_NOTENOUGH, "Cannot allocate memory for smoothing");

	np = picDithered;

	/* get first line of picture */

	for (j = w, tmpptr = nextline; j; j--, ep++)
		*tmpptr++ = (int) *ep;

	for (i = 0; i < h; i++) {
		np = picDithered + i * w;
		tmpptr = thisline;
		thisline = nextline;
		nextline = tmpptr; /* swap */

		if (i != imax)
			for (j = w, tmpptr = nextline; j; j--, ep++)
				*tmpptr++ = (int) *ep;

		/* dither a line, doing odd-lines right-to-left (serpentine) */
		thisptr = (i & 1) ? thisline + w - 1 : thisline;
		nextptr = (i & 1) ? nextline + w - 1 : nextline;
		if (i & 1)
			np += w - 1;

		for (j = 0; j < w; j++) {
			int k, d, mind, closest;

			r2 = *thisptr;
			if (i & 1)
				thisptr -= 1; /* move left */
			else
				thisptr += 1; /* move right */

			if (r2 < 0)
				r2 = 0;

			if (r2 > 255)
				r2 = 255;
			key = ((r2 & 0xf8) << 6) | ((r2 & 0xf8) << 1) | (r2 >> 4);
			/*
			const int maxKey=2 << 14;
			if (key >= maxKey)
				REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS,
						"overflow in DoColorDither");
			*/

			if (cache[key]) {
				*np = (byte) (cache[key] - 1);
				cnt1++;
			} else {
				/* not in cache, have to search the colortable */
				cnt2++;

				mind = 10000;
				for (k = closest = 0; k < 256 && mind > 7; k++) {
					d = 3*abs(r2 - k);
					if (d < mind) {
						mind = d;
						closest = k;
					}
				}
				cache[key] = closest + 1;
				*np = closest;
			}

			/* propagate the error */
			rerr = r2 - *np;

			if (j != jmax) { /* adjust LEFT/RIGHT pixel */
				int rerr_2 = rerr / 2;
				thisptr[0] += rerr_2;
				rerr -= rerr_2;
			}

			if (i != imax) { /* adjust BOTTOM pixel */
				nextptr[0] += rerr; /* possibly all err if we're at l/r edge */
			}

			if (i & 1) {
				nextptr -= 1;
				np--;
			} else {
				nextptr += 1;
				np++;
			}
		}
	}

	free(thisline);
	free(nextline);
	free(cache);
}
