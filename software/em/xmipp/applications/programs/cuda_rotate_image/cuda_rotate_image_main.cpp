/***************************************************************************
 *
 * Authors:     Amaya Jimenez (ajimenez@cnb.csic.es)
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

#include "../../../libraries/reconstruction_adapt_cuda/gpu_rotate_image.h"


//RUN_XMIPP_PROGRAM(ProgGpuRotateImage)

double timeval_diff(struct timeval *a, struct timeval *b)
{
  return
    (double)(a->tv_sec + (double)a->tv_usec/1000000) -
    (double)(b->tv_sec + (double)b->tv_usec/1000000);
}

int main(int argc, char** argv) {

	int aux;

	struct timeval t_inigpu, t_fingpu;
	double secsgpu;
	gettimeofday(&t_inigpu, NULL);

	ProgGpuRotateImage program;
	program.read(argc, argv);
    int errorCode = program.tryRun();

    gettimeofday(&t_fingpu, NULL);
    secsgpu = timeval_diff(&t_fingpu, &t_inigpu);
    printf("Total Rotate Image GPU: %.16g milisegundos\n", secsgpu * 1000.0);

    return errorCode;

}


