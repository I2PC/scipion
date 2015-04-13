/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef COLOR_H_
#define COLOR_H_

#include <stdio.h>
#include "xmipp_strings.h"
enum colorAttribute
{
    RESET     =  0,
    BRIGHT    =  1,
    DIM       =  2,
    UNDERLINE =  3,
    BLINK     =  4,
    REVERSE   =  7,
    HIDDEN    =  8
};

enum colorCode
{
    BLACK   =  0,
    RED     =  1,
    GREEN   =  2,
    YELLOW  =  3,
    BLUE    =  4,
    MAGENTA =  5,
    CYAN    =  6,
    WHITE   =  7
};

void textcolor(int attr, int fg, int bg);

String colorString(const char * msg, int color, int attribute = BRIGHT, int bgcolor=BLACK);


#endif /* COLOR_H_ */
