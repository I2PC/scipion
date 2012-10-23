/*
 * xmipp_macros.cpp
 *
 *  Created on: Jan 9, 2012
 *      Author: kino
 */



#include "xmipp_macros.h"

#if defined(__APPLE__) || defined(__MINGW32__)
void sincos(double angle, double * sine, double * cosine)
{
	*sine = sin(angle);
	*cosine = cos(angle);
}
#endif

