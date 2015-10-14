/***************************************************************************
 * Authors:     AUTHOR_NAME (jlvilas@cnb.csic.es)
 *
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

#ifndef ASSIGNATION_TILT_PAIR_H_
#define ASSIGNATION_TILT_PAIR_H_
#define PI 3.14159265

#include <data/xmipp_program.h>
#include <math.h>
#include <alglib/src/ap.h>

class ProgassignationTiltPair: public XmippProgram
{


public:
    /** Filenames */
    FileName fntilt, fnuntilt,  fndir, fnmic;

    /** Maximum shift*/
    int mshift;

    /** Particle size*/
    double particle_size;

    /** threshold*/
    double thr;

    MetaData mdPartial;

public:

    void defineParams();

    void readParams();

    void run();

};
#endif
