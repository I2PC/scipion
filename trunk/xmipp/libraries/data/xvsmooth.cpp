/*
 * xvsmooth.c - smoothing/color dither routines for XV
 *
 *  Author:    John Bradley, University of Pennsylvania
 *                (bradley@cis.upenn.edu)
 *
 *  Contains:
 *            byte *SmoothResize(src8, swide, shigh, dwide, dhigh,
 *                               rmap, gmap, bmap, rdmap, gdmap, bdmap, maplen)
 *            byte *Smooth24(pic824, is24, swide, shigh, dwide, dhigh,
 *                               rmap, gmap, bmap)
 *            byte *DoColorDither(pic24, pic8, w, h, rmap,gmap,bmap,
 *                                rdisp, gdisp, bdisp, maplen)
 *            byte *Do332ColorDither(pic24, pic8, w, h, rmap,gmap,bmap,
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

byte *Smooth24(byte *pic824, int is24, int swide, int shigh,
               int dwide, int dhigh, byte *rmap, byte *gmap, byte *bmap);
byte *DoColorDither(byte *pic24, byte *pic8, int w, int h,
                    byte *rmap, byte *gmap, byte *bmap,
                    byte *rdisp, byte *gdisp, byte *bdisp, int maplen);


#ifdef __STDC__
int SmoothX(byte *, byte *, int, int, int, int, int,
            byte *, byte *, byte *);
int SmoothY(byte *, byte *, int, int, int, int, int,
            byte *, byte *, byte *);
int SmoothXY(byte *, byte *, int, int, int, int, int,
             byte *, byte *, byte *);
#else
int  SmoothX(), SmoothY(), SmoothXY();
#endif


#define xvbzero(s,l) memset(s, 0, l)
/***************************************************/
byte *SmoothResize(byte *srcpic8, int swide, int shigh,
                   int dwide, int dhigh,
                   byte *rmap, byte *gmap, byte *bmap,
                   byte *rdmap, byte *gdmap, byte *bdmap, int maplen)
{
    /* generic interface to Smooth and ColorDither code.
       given an 8-bit-per, swide * shigh image with colormap rmap,gmap,bmap,
       will generate a new 8-bit-per, dwide * dhigh image, which is dithered
       using colors found in rdmap, gdmap, bdmap arrays */

    /* returns ptr to a dwide*dhigh array of bytes, or NULL on failure */

    byte *pic24, *pic8;

    pic24 = Smooth24(srcpic8, 0, swide, shigh, dwide, dhigh, rmap, gmap, bmap);

    if (pic24)
    {
        pic8 = DoColorDither(pic24, NULL, dwide, dhigh, rmap, gmap, bmap,
                             rdmap, gdmap, bdmap, maplen);
        free(pic24);
        return pic8;
    }

    return (byte *) NULL;
}



/***************************************************/
byte *Smooth24(byte *pic824, int is24, int swide, int shigh, int dwide,
               int dhigh, byte *rmap, byte *gmap, byte *bmap)
{
    /* does a SMOOTH resize from pic824 (which is either a swide*shigh, 8-bit
       pic, with colormap rmap,gmap,bmap OR a swide*shigh, 24-bit image, based
       on whether 'is24' is set) into a dwide * dhigh 24-bit image

       returns a dwide*dhigh 24bit image, or NULL on failure (malloc) */
    /* rmap,gmap,bmap should be 'desired' colors */

    byte *pic24, *pp;
    int  *cxtab, *pxtab;
    int   y1Off, cyOff;
    int   ex, ey, cx, cy, px, py, apx, apy, x1, y1;
    int   cA, cB, cC, cD;
    int   pA, pB, pC, pD;
    int   retval, bperpix;

    cA = cB = cC = cD = 0;
    pp = pic24 = (byte *) malloc(dwide * dhigh * 3);
    if (!pic24)
    {
        fprintf(stderr, "unable to malloc pic24 in 'Smooth24()'\n");
        return pic24;
    }

    bperpix = (is24) ? 3 : 1;

    /* decide which smoothing routine to use based on type of expansion */
    if (dwide <  swide && dhigh <  shigh)
        retval = SmoothXY(pic24, pic824, is24, swide, shigh, dwide, dhigh,
                          rmap, gmap, bmap);

    else if (dwide <  swide && dhigh >= shigh)
        retval = SmoothX(pic24, pic824, is24, swide, shigh, dwide, dhigh,
                         rmap, gmap, bmap);

    else if (dwide >= swide && dhigh <  shigh)
        retval = SmoothY(pic24, pic824, is24, swide, shigh, dwide, dhigh,
                         rmap, gmap, bmap);

    else
    {
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
        {
            free(pic24);
            return NULL;
        }

        pxtab = (int *) malloc(dwide * sizeof(int));
        if (!pxtab)
        {
            free(pic24);
            free(cxtab);
            return NULL;
        }

        for (ex = 0; ex < dwide; ex++)
        {
            cxtab[ex] = (ex * swide) / dwide;
            pxtab[ex] = (((ex * swide) * 100) / dwide)
                        - (cxtab[ex] * 100) - 50;
        }

        for (ey = 0; ey < dhigh; ey++)
        {
            cy = (ey * shigh) / dhigh;
            py = (((ey * shigh) * 100) / dhigh) - (cy * 100) - 50;
            if (py < 0)
            {
                y1 = cy - 1;
                if (y1 < 0) y1 = 0;
            }
            else
            {
                y1 = cy + 1;
                if (y1 > shigh - 1) y1 = shigh - 1;
            }

            cyOff = cy * swide * bperpix;    /* current line */
            y1Off = y1 * swide * bperpix;    /* up or down one line, depending */

            /*      if ((ey&15) == 0) WaitCursor(); */

            for (ex = 0; ex < dwide; ex++)
            {
                byte *pptr, rA, gA, bA, rB, gB, bB, rC, gC, bC, rD, gD, bD;

                cx = cxtab[ex];
                px = pxtab[ex];

                if (px < 0)
                {
                    x1 = cx - 1;
                    if (x1 < 0) x1 = 0;
                }
                else
                {
                    x1 = cx + 1;
                    if (x1 > swide - 1) x1 = swide - 1;
                }

                if (is24)
                {
                    pptr = pic824 + y1Off + x1 * bperpix;   /* corner pixel */
                    rA = *pptr++;
                    gA = *pptr++;
                    bA = *pptr++;

                    pptr = pic824 + y1Off + cx * bperpix;   /* up/down center pixel */
                    rB = *pptr++;
                    gB = *pptr++;
                    bB = *pptr++;

                    pptr = pic824 + cyOff + x1 * bperpix;   /* left/right center pixel */
                    rC = *pptr++;
                    gC = *pptr++;
                    bC = *pptr++;

                    pptr = pic824 + cyOff + cx * bperpix;   /* center pixel */
                    rD = *pptr++;
                    gD = *pptr++;
                    bD = *pptr++;
                }
                else
                {
                    /* 8-bit picture */
                    cA = pic824[y1Off + x1];   /* corner pixel */
                    cB = pic824[y1Off + cx];   /* up/down center pixel */
                    cC = pic824[cyOff + x1];   /* left/right center pixel */
                    cD = pic824[cyOff + cx];   /* center pixel */
                }

                /* quick check */
                if (!is24 && cA == cB && cB == cC && cC == cD)
                {
                    /* set this pixel to the same color as in pic8 */
                    *pp++ = rmap[cD];
                    *pp++ = gmap[cD];
                    *pp++ = bmap[cD];
                }

                else
                {
                    /* compute weighting factors */
                    apx = abs(px);
                    apy = abs(py);
                    pA = (apx * apy) / 100;
                    pB = (apy * (100 - apx)) / 100;
                    pC = (apx * (100 - apy)) / 100;
                    pD = 100 - (pA + pB + pC);

                    if (is24)
                    {
                        *pp++ = (pA * rA) / 100 + (pB * rB) / 100 +
                                (pC * rC) / 100 + (pD * rD) / 100;

                        *pp++ = (pA * gA) / 100 + (pB * gB) / 100 +
                                (pC * gC) / 100 + (pD * gD) / 100;

                        *pp++ = (pA * bA) / 100 + (pB * bB) / 100 +
                                (pC * bC) / 100 + (pD * bD) / 100;
                    }
                    else
                    {
                        /* 8-bit pic */
                        *pp++ = (pA * rmap[cA]) / 100 + (pB * rmap[cB]) / 100 +
                                (pC * rmap[cC]) / 100 + (pD * rmap[cD]) / 100;

                        *pp++ = (pA * gmap[cA]) / 100 + (pB * gmap[cB]) / 100 +
                                (pC * gmap[cC]) / 100 + (pD * gmap[cD]) / 100;

                        *pp++ = (pA * bmap[cA]) / 100 + (pB * bmap[cB]) / 100 +
                                (pC * bmap[cC]) / 100 + (pD * bmap[cD]) / 100;
                    }
                }
            }
        }

        free(cxtab);
        free(pxtab);
        retval = 0;    /* okay */
    }

    if (retval)
    {
        /* one of the Smooth**() methods failed */
        free(pic24);
        pic24 = (byte *) NULL;
    }

    return pic24;
}

/***************************************************/
int SmoothX(byte *pic24, byte *pic824, int is24,
            int swide, int shigh, int dwide, int dhigh,
            byte *rmap, byte *gmap, byte *bmap)
{
    byte *cptr, *cptr1;
    int  i, j;
    int  *lbufR, *lbufG, *lbufB;
    int  pixR, pixG, pixB, bperpix;
    int  pcnt0, pcnt1, lastpix, pixcnt, thisline, ypcnt;
    int  *pixarr, *paptr;

    /* returns '0' if okay, '1' if failed (malloc) */

    /* for case where pic8 is shrunk horizontally and stretched vertically
       maps pic8 into an dwide * dhigh 24-bit picture.  Only works correctly
       when swide>=dwide and shigh<=dhigh */


    /* malloc some arrays */
    lbufR = (int *) calloc(swide, sizeof(int));
    lbufG = (int *) calloc(swide, sizeof(int));
    lbufB = (int *) calloc(swide, sizeof(int));
    pixarr = (int *) calloc(swide + 1, sizeof(int));
    if (!lbufR || !lbufG || !lbufB || !pixarr)
    {
        if (lbufR)  free(lbufR);
        if (lbufG)  free(lbufG);
        if (lbufB)  free(lbufB);
        if (pixarr) free(pixarr);
        return 1;
    }

    bperpix = (is24) ? 3 : 1;

    for (j = 0; j <= swide; j++)
        pixarr[j] = (j * dwide + (15 * swide) / 16) / swide;

    cptr = pic824;
    cptr1 = cptr + swide * bperpix;

    for (i = 0; i < dhigh; i++)
    {
        /*    if ((i&15) == 0) WaitCursor();*/

        ypcnt = (((i * shigh) << 6) / dhigh) - 32;
        if (ypcnt < 0) ypcnt = 0;

        pcnt1 = ypcnt & 0x3f;                     /* 64ths of NEXT line to use */
        pcnt0 = 64 - pcnt1;                       /* 64ths of THIS line to use */

        thisline = ypcnt >> 6;

        cptr  = pic824 + thisline * swide * bperpix;
        if (thisline + 1 < shigh) cptr1 = cptr + swide * bperpix;
        else cptr1 = cptr;

        if (is24)
        {
            for (j = 0; j < swide; j++)
            {
                lbufR[j] = ((*cptr++ * pcnt0) + (*cptr1++ * pcnt1)) >> 6;
                lbufG[j] = ((*cptr++ * pcnt0) + (*cptr1++ * pcnt1)) >> 6;
                lbufB[j] = ((*cptr++ * pcnt0) + (*cptr1++ * pcnt1)) >> 6;
            }
        }
        else
        {
            /* 8-bit input pic */
            for (j = 0; j < swide; j++, cptr++, cptr1++)
            {
                lbufR[j] = ((rmap[*cptr] * pcnt0) + (rmap[*cptr1] * pcnt1)) >> 6;
                lbufG[j] = ((gmap[*cptr] * pcnt0) + (gmap[*cptr1] * pcnt1)) >> 6;
                lbufB[j] = ((bmap[*cptr] * pcnt0) + (bmap[*cptr1] * pcnt1)) >> 6;
            }
        }

        pixR = pixG = pixB = pixcnt = lastpix = 0;

        for (j = 0, paptr = pixarr; j <= swide; j++, paptr++)
        {
            if (*paptr != lastpix)
            {
                /* write a pixel to pic24 */
                *pic24++ = pixR / pixcnt;
                *pic24++ = pixG / pixcnt;
                *pic24++ = pixB / pixcnt;
                lastpix = *paptr;
                pixR = pixG = pixB = pixcnt = 0;
            }

            if (j < swide)
            {
                pixR += lbufR[j];
                pixG += lbufG[j];
                pixB += lbufB[j];
                pixcnt++;
            }
        }
    }

    free(lbufR);
    free(lbufG);
    free(lbufB);
    free(pixarr);
    return 0;
}

/***************************************************/
int SmoothY(byte *pic24, byte *pic824, int is24, int swide,
            int shigh, int dwide, int dhigh,
            byte *rmap, byte *gmap, byte *bmap)
{
    byte *clptr, *cptr, *cptr1;
    int  i, j, bperpix;
    int  *lbufR, *lbufG, *lbufB, *pct0, *pct1, *cxarr, *cxptr;
    int  lastline, thisline, linecnt;
    int  retval;


    /* returns '0' if okay, '1' if failed (malloc) */

    /* for case where pic8 is shrunk vertically and stretched horizontally
       maps pic8 into a dwide * dhigh 24-bit picture.  Only works correctly
       when swide<=dwide and shigh>=dhigh */

    retval = 0;   /* no probs, yet... */

    bperpix = (is24) ? 3 : 1;

    lbufR = lbufG = lbufB = pct0 = pct1 = cxarr = NULL;
    lbufR = (int *) calloc(dwide, sizeof(int));
    lbufG = (int *) calloc(dwide, sizeof(int));
    lbufB = (int *) calloc(dwide, sizeof(int));
    pct0  = (int *) calloc(dwide, sizeof(int));
    pct1  = (int *) calloc(dwide, sizeof(int));
    cxarr = (int *) calloc(dwide, sizeof(int));

    if (!lbufR || !lbufG || !lbufB || !pct0 || ! pct1 || !cxarr)
    {
        retval = 1;
        goto smyexit;
    }



    for (i = 0; i < dwide; i++)
    {
        /* precompute some handy tables */
        int cx64;
        cx64 = (((i * swide) << 6) / dwide) - 32;
        if (cx64 < 0) cx64 = 0;
        pct1[i] = cx64 & 0x3f;
        pct0[i] = 64 - pct1[i];
        cxarr[i] = cx64 >> 6;
    }


    lastline = linecnt = 0;

    for (i = 0, clptr = pic824; i <= shigh; i++, clptr += swide * bperpix)
    {
        /*    if ((i&15) == 0) WaitCursor();*/

        thisline = (i * dhigh + (15 * shigh) / 16) / shigh;

        if (thisline != lastline)
        {
            /* copy a line to pic24 */
            for (j = 0; j < dwide; j++)
            {
                *pic24++ = lbufR[j] / linecnt;
                *pic24++ = lbufG[j] / linecnt;
                *pic24++ = lbufB[j] / linecnt;
            }

            xvbzero((char *) lbufR, dwide * sizeof(int));  /* clear out line bufs */
            xvbzero((char *) lbufG, dwide * sizeof(int));
            xvbzero((char *) lbufB, dwide * sizeof(int));
            linecnt = 0;
            lastline = thisline;
        }


        for (j = 0, cxptr = cxarr; j < dwide; j++, cxptr++)
        {
            cptr  = clptr + *cxptr * bperpix;
            if (*cxptr < swide - 1) cptr1 = cptr + 1 * bperpix;
            else cptr1 = cptr;

            if (is24)
            {
                lbufR[j] += (((*cptr++ * pct0[j]) + (*cptr1++ * pct1[j])) >> 6);
                lbufG[j] += (((*cptr++ * pct0[j]) + (*cptr1++ * pct1[j])) >> 6);
                lbufB[j] += (((*cptr++ * pct0[j]) + (*cptr1++ * pct1[j])) >> 6);
            }
            else
            {
                /* 8-bit input pic */
                lbufR[j] += (((rmap[*cptr] * pct0[j]) + (rmap[*cptr1] * pct1[j])) >> 6);
                lbufG[j] += (((gmap[*cptr] * pct0[j]) + (gmap[*cptr1] * pct1[j])) >> 6);
                lbufB[j] += (((bmap[*cptr] * pct0[j]) + (bmap[*cptr1] * pct1[j])) >> 6);
            }
        }

        linecnt++;
    }


smyexit:
    if (lbufR) free(lbufR);
    if (lbufG) free(lbufG);
    if (lbufB) free(lbufB);
    if (pct0)  free(pct0);
    if (pct1)  free(pct1);
    if (cxarr) free(cxarr);

    return retval;
}

/***************************************************/
int SmoothXY(byte *pic24, byte *pic824, int is24,
             int swide, int shigh, int dwide, int dhigh,
             byte *rmap, byte *gmap, byte *bmap)
{
    byte *cptr;
    int  i, j;
    int  *lbufR, *lbufG, *lbufB;
    int  pixR, pixG, pixB, bperpix;
    int  lastline, thisline, lastpix, linecnt, pixcnt;
    int  *pixarr, *paptr;


    /* returns '0' if okay, '1' if failed (malloc) */

    /* shrinks pic8 into a dwide * dhigh 24-bit picture.  Only works correctly
       when swide>=dwide and shigh>=dhigh (ie, the picture is shrunk on both
       axes) */


    /* malloc some arrays */
    lbufR = (int *) calloc(swide, sizeof(int));
    lbufG = (int *) calloc(swide, sizeof(int));
    lbufB = (int *) calloc(swide, sizeof(int));
    pixarr = (int *) calloc(swide + 1, sizeof(int));
    if (!lbufR || !lbufG || !lbufB || !pixarr)
    {
        if (lbufR)  free(lbufR);
        if (lbufG)  free(lbufG);
        if (lbufB)  free(lbufB);
        if (pixarr) free(pixarr);
        return 1;
    }

    bperpix = (is24) ? 3 : 1;

    for (j = 0; j <= swide; j++)
        pixarr[j] = (j * dwide + (15 * swide) / 16) / swide;

    lastline = linecnt = pixR = pixG = pixB = 0;
    cptr = pic824;

    for (i = 0; i <= shigh; i++)
    {
        /*    if ((i&15) == 0) WaitCursor();*/

        thisline = (i * dhigh + (15 * shigh) / 16) / shigh;

        if ((thisline != lastline))
        {
            /* copy a line to pic24 */
            pixR = pixG = pixB = pixcnt = lastpix = 0;

            for (j = 0, paptr = pixarr; j <= swide; j++, paptr++)
            {
                if (*paptr != lastpix)
                {
                    /* write a pixel to pic24 */
                    *pic24++ = (pixR / linecnt) / pixcnt;
                    *pic24++ = (pixG / linecnt) / pixcnt;
                    *pic24++ = (pixB / linecnt) / pixcnt;
                    lastpix = *paptr;
                    pixR = pixG = pixB = pixcnt = 0;
                }

                if (j < swide)
                {
                    pixR += lbufR[j];
                    pixG += lbufG[j];
                    pixB += lbufB[j];
                    pixcnt++;
                }
            }

            lastline = thisline;
            xvbzero((char *) lbufR, swide * sizeof(int));  /* clear out line bufs */
            xvbzero((char *) lbufG, swide * sizeof(int));
            xvbzero((char *) lbufB, swide * sizeof(int));
            linecnt = 0;
        }

        if (i < shigh)
        {
            if (is24)
            {
                for (j = 0; j < swide; j++)
                {
                    lbufR[j] += *cptr++;
                    lbufG[j] += *cptr++;
                    lbufB[j] += *cptr++;
                }
            }
            else
            {
                for (j = 0; j < swide; j++, cptr++)
                {
                    lbufR[j] += rmap[*cptr];
                    lbufG[j] += gmap[*cptr];
                    lbufB[j] += bmap[*cptr];
                }
            }

            linecnt++;
        }
    }

    free(lbufR);
    free(lbufG);
    free(lbufB);
    free(pixarr);
    return 0;
}

/********************************************/
byte *DoColorDither(byte *pic24, byte *pic8, int w, int h, byte *rmap,
                    byte *gmap, byte *bmap, byte *rdisp, byte *gdisp, byte *bdisp, int maplen)
{
    /* takes a 24 bit picture, of size w*h, dithers with the colors in
       rdisp, gdisp, bdisp (which have already been allocated),
       and generates an 8-bit w*h image, which it returns.
       ignores input value 'pic8'
       returns NULL on error

       note: the rdisp,gdisp,bdisp arrays should be the 'displayed' colors,
       not the 'desired' colors

       if pic24 is NULL, uses the passed-in pic8 (an 8-bit image) as
       the source, and the rmap,gmap,bmap arrays as the desired colors */

    byte *np, *ep, *newpic;
    short *cache;
    int r2, g2, b2;
    int *thisline, *nextline, *thisptr, *nextptr, *tmpptr;
    int  i, j, rerr, gerr, berr, pwide3;
    int  imax, jmax;
    int key;
    long cnt1, cnt2;

    cnt1 = cnt2 = 0;
    pwide3 = w * 3;
    imax = h - 1;
    jmax = w - 1;
    ep = (pic24) ? pic24 : pic8;

    /* attempt to malloc things */
    newpic = (byte *) malloc(w * h);
    cache = (short *) calloc(2 << 14, sizeof(short));
    thisline = (int *) malloc(pwide3 * sizeof(int));
    nextline = (int *) malloc(pwide3 * sizeof(int));
    if (!cache || !newpic || !thisline || !nextline)
    {
        if (newpic)   free(newpic);
        if (cache)    free(cache);
        if (thisline) free(thisline);
        if (nextline) free(nextline);

        return (byte *) NULL;
    }

    np = newpic;

    /* get first line of picture */

    if (pic24)
    {
        for (j = pwide3, tmpptr = nextline; j; j--, ep++)
            *tmpptr++ = (int) * ep;
    }
    else
    {
        for (j = w, tmpptr = nextline; j; j--, ep++)
        {
            *tmpptr++ = (int) rmap[*ep];
            *tmpptr++ = (int) gmap[*ep];
            *tmpptr++ = (int) bmap[*ep];
        }
    }


    for (i = 0; i < h; i++)
    {
        np = newpic + i * w;
        /*    if ((i&15) == 0) WaitCursor();*/

        tmpptr = thisline;
        thisline = nextline;
        nextline = tmpptr;   /* swap */

        if (i != imax)
        {
            /* get next line */
            if (!pic24)
                for (j = w, tmpptr = nextline; j; j--, ep++)
                {
                    *tmpptr++ = (int) rmap[*ep];
                    *tmpptr++ = (int) gmap[*ep];
                    *tmpptr++ = (int) bmap[*ep];
                }
            else
                for (j = pwide3, tmpptr = nextline; j; j--, ep++) *tmpptr++ = (int) * ep;
        }

        /* dither a line, doing odd-lines right-to-left (serpentine) */
        thisptr = (i & 1) ? thisline + w * 3 - 3 : thisline;
        nextptr = (i & 1) ? nextline + w * 3 - 3 : nextline;
        if (i&1) np += w - 1;


        for (j = 0; j < w; j++)
        {
            int k, d, mind, closest;

            r2 = *thisptr++;
            g2 = *thisptr++;
            b2 = *thisptr++;
            if (i&1) thisptr -= 6;  /* move left */

            /* map r2,g2,b2 components (could be outside 0..255 range)
            into 0..255 range */

            if (r2 < 0 || g2 < 0 || b2 < 0)
            {
                /* are there any negatives in RGB? */
                if (r2 < g2)
                {
                    if (r2 < b2) k = 0;
                    else k = 2;
                }
                else
                {
                    if (g2 < b2) k = 1;
                    else k = 2;
                }

                switch (k)
                {
                case 0:
                    g2 -= r2;
                    b2 -= r2;
                    d = (abs(r2) * 3) / 2;    /* RED */
                    r2 = 0;
                    g2 = (g2 > d) ? g2 - d : 0;
                    b2 = (b2 > d) ? b2 - d : 0;
                    break;

                case 1:
                    r2 -= g2;
                    b2 -= g2;
                    d = (abs(g2) * 3) / 2;    /* GREEN */
                    r2 = (r2 > d) ? r2 - d : 0;
                    g2 = 0;
                    b2 = (b2 > d) ? b2 - d : 0;
                    break;

                case 2:
                    r2 -= b2;
                    g2 -= b2;
                    d = (abs(b2) * 3) / 2;    /* BLUE */
                    r2 = (r2 > d) ? r2 - d : 0;
                    g2 = (g2 > d) ? g2 - d : 0;
                    b2 = 0;
                    break;
                }
            }

            if (r2 > 255 || g2 > 255 || b2 > 255)
            {
                /* any overflows in RGB? */
                if (r2 > g2)
                {
                    if (r2 > b2) k = 0;
                    else k = 2;
                }
                else
                {
                    if (g2 > b2) k = 1;
                    else k = 2;
                }

                switch (k)
                {
                case 0:
                    g2 = (g2 * 255) / r2;
                    b2 = (b2 * 255) / r2;
                    r2 = 255;
                    break;
                case 1:
                    r2 = (r2 * 255) / g2;
                    b2 = (b2 * 255) / g2;
                    g2 = 255;
                    break;
                case 2:
                    r2 = (r2 * 255) / b2;
                    g2 = (g2 * 255) / b2;
                    b2 = 255;
                    break;
                }
            }

            key = ((r2 & 0xf8) << 6) | ((g2 & 0xf8) << 1) | (b2 >> 4);
            if (key >= (2 << 14))
            {
                fprintf(stderr, "'key' overflow in DoColorDither()");
                exit(-1);
            }

            if (cache[key])
            {
                *np = (byte)(cache[key] - 1);
                cnt1++;
            }
            else
            {
                /* not in cache, have to search the colortable */
                cnt2++;

                mind = 10000;
                for (k = closest = 0; k < maplen && mind > 7; k++)
                {
                    d = abs(r2 - rdisp[k])
                        + abs(g2 - gdisp[k])
                        + abs(b2 - bdisp[k]);
                    if (d < mind)
                    {
                        mind = d;
                        closest = k;
                    }
                }
                cache[key] = closest + 1;
                *np = closest;
            }


            /* propogate the error */
            rerr = r2 - rdisp[*np];
            gerr = g2 - gdisp[*np];
            berr = b2 - bdisp[*np];

            if (j != jmax)
            {
                /* adjust LEFT/RIGHT pixel */
                thisptr[0] += (rerr / 2);
                rerr -= (rerr / 2);
                thisptr[1] += (gerr / 2);
                gerr -= (gerr / 2);
                thisptr[2] += (berr / 2);
                berr -= (berr / 2);
            }

            if (i != imax)
            {
                /* adjust BOTTOM pixel */
                nextptr[0] += rerr;    /* possibly all err if we're at l/r edge */
                nextptr[1] += gerr;
                nextptr[2] += berr;
            }

            if (i&1)
            {
                nextptr -= 3;
                np--;
            }
            else
            {
                nextptr += 3;
                np++;
            }
        }
    }


    free(thisline);
    free(nextline);
    free(cache);

    return newpic;
}

/********************************************/
byte *Do332ColorDither(byte *pic24, byte *pic8, int w, int h, byte *rmap,
                       byte *gmap, byte *bmap, byte *rdisp, byte *gdisp, byte *bdisp, int maplen)
{
    /* some sort of color dither optimized for the 332 std cmap */

    /* takes a 24 bit picture, of size w*h, dithers with the colors in
       rdisp, gdisp, bdisp (which have already been allocated),
       and generates an 8-bit w*h image, which it returns.
       ignores input value 'pic8'
       returns NULL on error

       note: the rdisp,gdisp,bdisp arrays should be the 'displayed' colors,
       not the 'desired' colors

       if pic24 is NULL, uses the passed-in pic8 (an 8-bit image) as
       the source, and the rmap,gmap,bmap arrays as the desired colors */

    byte *np, *ep, *newpic;
    int r2, g2, b2;
    int *thisline, *nextline, *thisptr, *nextptr, *tmpptr;
    int  i, j, rerr, gerr, berr, pwide3;
    int  imax, jmax;
    long cnt1, cnt2;

    cnt1 = cnt2 = 0;
    pwide3 = w * 3;
    imax = h - 1;
    jmax = w - 1;

    /* attempt to malloc things */
    newpic = (byte *) malloc(w * h);
    thisline = (int *) malloc(pwide3 * sizeof(int));
    nextline = (int *) malloc(pwide3 * sizeof(int));
    if (!newpic || !thisline || !nextline)
    {
        if (newpic)   free(newpic);
        if (thisline) free(thisline);
        if (nextline) free(nextline);

        return (byte *) NULL;
    }

    np = newpic;
    ep = (pic24) ? pic24 : pic8;


    /* get first line of picture */

    if (pic24)
    {
        for (j = pwide3, tmpptr = nextline; j; j--, ep++) *tmpptr++ = (int) * ep;
    }
    else
    {
        for (j = w, tmpptr = nextline; j; j--, ep++)
        {
            *tmpptr++ = (int) rmap[*ep];
            *tmpptr++ = (int) gmap[*ep];
            *tmpptr++ = (int) bmap[*ep];
        }
    }


    for (i = 0; i < h; i++)
    {
        np = newpic + i * w;
        /*    if ((i&15) == 0) WaitCursor();*/

        tmpptr = thisline;
        thisline = nextline;
        nextline = tmpptr;   /* swap */

        if (i != imax)
        {
            /* get next line */
            if (!pic24)
                for (j = w, tmpptr = nextline; j; j--, ep++)
                {
                    *tmpptr++ = (int) rmap[*ep];
                    *tmpptr++ = (int) gmap[*ep];
                    *tmpptr++ = (int) bmap[*ep];
                }
            else
                for (j = pwide3, tmpptr = nextline; j; j--, ep++) *tmpptr++ = (int) * ep;
        }


        /* dither a line, doing odd-lines right-to-left (serpentine) */
        thisptr = (i & 1) ? thisline + w * 3 - 3 : thisline;
        nextptr = (i & 1) ? nextline + w * 3 - 3 : nextline;
        if (i&1) np += w - 1;


        for (j = 0; j < w; j++)
        {
            int k, d, mind, closest, rb, gb, bb;

            r2 = *thisptr++;
            g2 = *thisptr++;
            b2 = *thisptr++;
            if (i&1) thisptr -= 6;  /* move left */

            rb = (r2 + 0x10);    /* round top 3 bits */
            RANGE(rb, 0, 255);
            rb = rb & 0xe0;

            gb = (g2 + 0x10);    /* round 3 bits */
            RANGE(gb, 0, 255);
            gb = gb & 0xe0;

            bb = (b2 + 0x20);    /* round 2 bits */
            RANGE(bb, 0, 255);
            bb = bb & 0xc0;


            *np = rb | (gb >> 3) | (bb >> 6);

            /* propogate the error */
            rerr = r2 - rdisp[*np];
            gerr = g2 - gdisp[*np];
            berr = b2 - bdisp[*np];

            if (j != jmax)
            {
                /* adjust LEFT/RIGHT pixel */
                thisptr[0] += (rerr / 2);
                rerr -= (rerr / 2);
                thisptr[1] += (gerr / 2);
                gerr -= (gerr / 2);
                thisptr[2] += (berr / 2);
                berr -= (berr / 2);
            }

            if (i != imax)
            {
                /* adjust BOTTOM pixel */
                nextptr[0] += rerr;    /* possibly all err if we're at l/r edge */
                nextptr[1] += gerr;
                nextptr[2] += berr;
            }

            if (i&1)
            {
                nextptr -= 3;
                np--;
            }
            else
            {
                nextptr += 3;
                np++;
            }
        }
    }


    free(thisline);
    free(nextline);

    return newpic;
}
