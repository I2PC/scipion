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

#include "gpu_rotate_image.h"

#include <reconstruction_cuda/cuda_gpu_rotate_image.h>
#include <reconstruction_cuda/cuda_gpu_rotate_image_v2.h>
#include <data/xmipp_image.h>
#include "data/xmipp_program.h"
#include <data/args.h>

//#define PI ((float)3.141592653589793)


// Read arguments ==========================================================
void ProgGpuRotateImage::readParams()
{
	fnRef = getParam("-i");
	fnOut = getParam("-o");
	ang = getParam("--ang");
	interp = getParam("--interp");

}

// Show ====================================================================
void ProgGpuRotateImage::show()
{
    std::cout
	<< "Input:          " << fnRef    << std::endl
	<< "Output:          " << fnOut    << std::endl
    << "Rotate an image by an angle: " << ang       << std::endl
	<< "The user can choose the interpolation method: 0 - Point, 1 - Linear, 2 - Cubic " << interp    << std::endl
    ;
}

// usage ===================================================================
void ProgGpuRotateImage::defineParams()
{
    addUsageLine("Computes the rotation of an image with CUDA in GPU");
    addParamsLine("   -i <Metadata1>        : Input image");
    addParamsLine("   -o <Metadata1>        : Output image");
    addParamsLine("   --ang <Metadata1>        : Rotation angle in degrees");
    addParamsLine("   --interp <Metadata1>        : Interpolation method: 0 - Point, 1 - Linear, 2 - Cubic");

}

//#define DEBUG
// Compute distance --------------------------------------------------------
void ProgGpuRotateImage::run()
{

	int ang2, interpol;
	ang2 = ang.getNumber();
	interpol = interp.getNumber();

	float rad = PI*(float)ang2/180.0;
    std::cout << "Inside run with deg " << ang2 << " and rad " << rad << std::endl;

    Image<float> Iref, Iout;
    Iref.read(fnRef);
    size_t Xdim, Ydim, Zdim, Ndim;
    Iref.getDimensions(Xdim, Ydim, Zdim, Ndim);
    if (Zdim>1 || Ndim>1){
    	REPORT_ERROR(ERR_MATRIX_DIM,"Problem with image dimensions");
    }

    MultidimArray<float> &original_image=Iref();

    MultidimArray<float> rotated_image(Xdim,Ydim);

    float *original_image_gpu = MULTIDIM_ARRAY(original_image);
    float *rotated_image_gpu = MULTIDIM_ARRAY(rotated_image);

    if (interpol<2){
    	cuda_rotate_image(original_image_gpu, rotated_image_gpu, Xdim, Ydim, rad, interpol);
    }else{
    	cuda_rotate_image_v2(original_image_gpu, rotated_image_gpu, Xdim, Ydim, rad);
    }

    Iout() = rotated_image;
    Iout.write(fnOut);


}

