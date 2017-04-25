/***************************************************************************
 *
 * Authors:    Amaya Jimenez      ajimenez@cnb.csic.es (2002)
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
#ifndef _PROG_GPU_CORRELATION
#define _PROG_GPU_CORRELATION

#include <data/xmipp_program.h>
#include <data/xmipp_image.h>
#include <data/xmipp_program.h>
#include <data/transformations.h>
#include <data/args.h>
#include <data/matrix2d.h>

#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#include "../reconstruction_cuda/cuda_gpu_correlation.h"


//CUDA functions
void cuda_correlation(float *image, float *rotated_image, float *module, float *angle, size_t Xdim, size_t Ydim, size_t Zdim);


class ProgGpuCorrelation: public XmippProgram
{
private:
    //Input file
	FileName fn_sel;
	//Output file
    FileName fn_out;

public:

    //Input metadata file
    MetaData SF;

    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Define parameters
    void defineParams();

    //void preProcess();

    /** processImage */
    void run();

};
//@}
#endif
