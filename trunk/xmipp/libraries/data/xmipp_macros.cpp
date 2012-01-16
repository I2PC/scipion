/*
 * xmipp_macros.cpp
 *
 *  Created on: Jan 9, 2012
 *      Author: kino
 */



#include "xmipp_macros.h"

#ifdef __APPLE__
#include <math.h>
void sincos(double angle, double * sine, double * cosine)
{
	*sine = sin(angle);
	*cosine = cos(angle);
}
#endif
