/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Copyright (c) 2000 , CSIC .
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 *
 *****************************************************************************/

#ifndef _XV_SMOOTH_H
#define _XV_SMOOTH_H
typedef unsigned char byte;

byte *SmoothResize(byte *srcpic8, int swide, int shigh,
                   int dwide, int dhigh,
                   byte *rmap, byte *gmap, byte *bmap,
                   byte *rdmap, byte *gdmap, byte *bdmap, int maplen);
#endif
