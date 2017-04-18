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
#ifndef _PROG_GPU_ROTATE_IMAGE_TH
#define _PROG_GPU_ROTATE_IMAGE_TH

#include <data/xmipp_program.h>
#include <data/xmipp_image.h>
#include <data/xmipp_program.h>
#include <data/transformations.h>
#include <data/args.h>
#include <data/matrix2d.h>
#include <data/xmipp_threads.h>

#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#include "../reconstruction_cuda/cuda_gpu_rotate_image_bspline.h"
#include "../reconstruction_cuda/cuda_gpu_rotate_image_linear_nearest_threads.h"
#include "../reconstruction_cuda/cuda_gpu_rotate_image_linear_nearest.h"

#define CUDA_NEAREST 0 //xmipp NEAREST
#define CUDA_LINEAR  1 //xmipp linear
#define CUDA_BSPLINE 2 //xmipp_BSPLINE3

#define CUDA_BORDER 3 //cudaAddressModeBorder =  3
#define CUDA_CLAMP  1 //cudaAddressModeClamp =  1
#define CUDA_WRAP 0 //cudaAddressModeWrap =  0

#define USE_MD 0
#define USE_COMMAND 1
#define USE_BOTH 2

#define NUM_TH 2

/**@defgroup AngularDistance angular_distance (Distance between two angular assignments)
   @ingroup ReconsLibrary */
//@{
/** Angular Distance parameters. */


void transformationMatrixCudaThreads(const MDRow &imageGeo, Matrix2D<double> &A, bool isVol);

//CUDA functions
//void cuda_rotate_image_threads(float **image, float **rotated_image, size_t *Xdim, size_t *Ydim, size_t *Zdim, size_t mdInSize, double** transform_matrix, int interp, int wrap, int first_call);


class ProgGpuRotateImageThreads: public XmippProgram
{
private:
    //Input file
	FileName fn_sel;
	//Output file
    FileName fn_out;

    int splineDegree, dim, wrapMode;
    Matrix2D<double> R;
    bool isVol, isMd, flip;
    int first_call;

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
