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
#include <data/args.h>


// Read arguments ==========================================================
void ProgGpuRotateImage::readParams()
{
	//ang = getParam("--ang");
	//interp = getParam("--interp");

	(checkParam("--ang"))? getParam("--ang"):  REPORT_ERROR(ERR_ARG_MISSING, (std::string)"Rotation angle unspecified");
	(checkParam("--interp"))? getParam("--interp"): 1;

}

// Show ====================================================================
void ProgGpuRotateImage::show()
{
    std::cout
    << "Rotate an image by an angle: " << ang       << std::endl
	<< "Choose the interpolation method: 0 - Point, 1 - Linear " << interp    << std::endl
    ;
}

// usage ===================================================================
void ProgGpuRotateImage::defineParams()
{
    addUsageLine("Computes the rotation of an image with CUDA in GPU");
    addParamsLine("   --ang <Metadata1>        : Rotation angle in degrees");
    addParamsLine("   --interp <Metadata1>        : Interpolation method");

}

//#define DEBUG
// Compute distance --------------------------------------------------------
void ProgGpuRotateImage::run()
{

	int num=5;
    int ang2 = ang.getNumber();
    int interpol = interp.getNumber();
    float rad = 3.14159*(float)ang2/180.0;
    std::cout << "Inside run with deg " << ang2 << " and rad " << rad << std::endl;
    MultidimArray<float> original_image(num,num);
    MultidimArray<float> rotated_image(num,num);

    /*A2D_ELEM(original_image,0,0) = 38.;
    A2D_ELEM(original_image,1,0) = 39.;
    A2D_ELEM(original_image,0,1) = 118.;
    A2D_ELEM(original_image,1,1) = 13.;*/
    original_image.initRandom(0, 100, RND_UNIFORM);

    float *original_image_gpu = MULTIDIM_ARRAY(original_image);
    float *rotated_image_gpu = MULTIDIM_ARRAY(rotated_image);

    cuda_rotate_image(original_image_gpu, rotated_image_gpu, rad, interpol);

	std::cout << "original_image" << original_image << std::endl;
	std::cout << "rotated_image" << rotated_image << std::endl;


}

