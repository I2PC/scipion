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

#include "gpu_example3.h"

#include <reconstruction_cuda/cuda_gpu_example3.h>
#include <data/args.h>


// Read arguments ==========================================================
void ProgGpuExample3::readParams()
{
	matSize = getParam("--matSize");

}

// Show ====================================================================
void ProgGpuExample3::show()
{
    std::cout
    << "Size of the matrixes: " << matSize       << std::endl
    ;
}

// usage ===================================================================
void ProgGpuExample3::defineParams()
{
    addUsageLine("Computes the sum between two matrixes with CUDA in GPU");
    addParamsLine("   --matSize <Metadata1>        : Size of the matrixes");

}

//#define DEBUG
// Compute distance --------------------------------------------------------
void ProgGpuExample3::run()
{

    int num = matSize.getNumber();
    std::cout << "Inside run with " << num << std::endl;
    MultidimArray<float> m1(num,num);
    MultidimArray<float> m2(num,num);
    MultidimArray<float> mResult(num,num);

    m1.initRandom(0, 10, RND_UNIFORM);
    m2.initRandom(100, 200, RND_UNIFORM);

    float *m1_gpu = MULTIDIM_ARRAY(m1);
    float *m2_gpu = MULTIDIM_ARRAY(m2);
    float *result = MULTIDIM_ARRAY(mResult);


    cuda_funcion(m1_gpu, m2_gpu, result, num);

	std::cout << "m1" << m1 << std::endl;
	std::cout << "m2" << m2 << std::endl;
	std::cout << "mResult" << mResult << std::endl;


}

