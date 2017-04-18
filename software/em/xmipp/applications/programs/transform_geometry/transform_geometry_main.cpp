/***************************************************************************
 *
 * Authors:     Carlos Oscar (           coss@cnb.csic.es (2000)
 *              Roberto Marabini (added fourier option)
 *              Jose Miguel de la Rosa Trevin (fusion of shift, rotate and scale)
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

#include "data/transform_geometry.h"

#include <stdio.h>
#include <time.h>
#include <sys/time.h>

//RUN_XMIPP_PROGRAM(ProgTransformGeometry)

/* retorna "a - b" en segundos */
double timeval_diff(struct timeval *a, struct timeval *b)
{
  return
    (double)(a->tv_sec + (double)a->tv_usec/1000000) -
    (double)(b->tv_sec + (double)b->tv_usec/1000000);
}

int main(int argc, char** argv) {

	int aux;

	struct timeval t_inicpu, t_fincpu;
	double secscpu;
	gettimeofday(&t_inicpu, NULL);

	ProgTransformGeometry program;
	program.read(argc, argv);
    int errorCode = program.tryRun();

    gettimeofday(&t_fincpu, NULL);
    secscpu = timeval_diff(&t_fincpu, &t_inicpu);
    printf("Total Transform Geometry CPU: %.16g milisegundos\n", secscpu * 1000.0);

    return errorCode;

}
