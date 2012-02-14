/***************************************************************************
 * Authors:     Joaquin Oton (joton@cnb.csic.es)
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

#ifndef MPI_RECONSTRUCT_ART_H_
#define MPI_RECONSTRUCT_ART_H_

#include "xmipp_mpi.h"
#include <reconstruction/reconstruct_art.h>

class ProgMPIReconsArt: public ProgReconsART, public XmippMpiProgram
{
public:
    /** Empty constructor */
    ProgMPIReconsArt();

    /* Constructor */
    ProgMPIReconsArt(int argc, char *argv[]);

    /* constructor providing an MpiNode
     * this is useful for using this programs from others
     */
    ProgMPIReconsArt(MpiNode * node);

    /* Run --------------------------------------------------------------------- */
    void run();
};

/* ------------------------------------------------------------------------- */
/* Time managing stuff                                                       */
/* ------------------------------------------------------------------------- */

typedef struct
{
    double user;  /* User time. */
    double sys;   /* System time. */
    double cpu;   /* CPU time = User + System. */
    double wall;  /* Wall time. */
}
USWtime_t;

// Gets User and System times for use with MPI
void uswtime(USWtime_t *tm);

#endif /* MPI_RECONSTRUCT_ART_H_ */
