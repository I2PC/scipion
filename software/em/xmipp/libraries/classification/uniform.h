/***************************************************************************
 *
 * Authors:     Alberto Pascual
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
//-----------------------------------------------------------------------------
// RandomUniformGenerator.h
//-----------------------------------------------------------------------------

#ifndef XMIPP_UNIFORM_H
#define XMIPP_UNIFORM_H

#include <stdlib.h>  // RAND_MAX
#include <time.h>    // RAND_MAX


//-----------------------------------------------------------------------------
// uniform return values uniformly distributed over the interval [min, max)
//-----------------------------------------------------------------------------

/**@defgroup UniformRandomGenerator Uniform random generator
   @ingroup ClassificationLibrary */
//@{
/**
* Template class for uniform random number generation
*/
template<class T> class RandomUniformGenerator
{
public:
    /**
    * Constructor
    * Parameter: _min   minimum value for the generated ramdom number
    * Parameter: _max   maximum value for the generated ramdom number
    */
    RandomUniformGenerator(T _min = 0, T _max = 1): thisMin(_min), diff(_max - _min)
    {
        srand((unsigned)time(NULL));
    }

    /**
    * () operator
    * Returns a random number between _min and _max
    *
    */
    T operator()()
    {
        return thisMin + T(diff * rand() / RAND_MAX);
    }

private:
    T thisMin;
    T diff;
};

//@}
#endif//XMIPP_UNIFORM_H
