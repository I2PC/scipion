/***************************************************************************
 *
 * Authors:    Amaya Jimenez      ajimenez@cnb.csic.es (2017)
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

#include "gpu_correlation.h"


//CUDA functions
void cuda_correlation(float *image, float *fft_image, float *module, float *angle, size_t Xdim, size_t Ydim, size_t Zdim){

	cuda_fft(image, fft_image, module, angle, Xdim, Ydim, Zdim);

}


// Read arguments ==========================================================
void ProgGpuCorrelation::readParams()
{

    fn_sel = getParam("-i");
    fn_out = getParam("-o");

}

// Show ====================================================================

void ProgGpuCorrelation::show()
{
    std::cout
	<< "Input:          " << fn_sel    << std::endl
	<< "Output:          " << fn_out    << std::endl
    ;
}

// usage ===================================================================
void ProgGpuCorrelation::defineParams()
{

	addParamsLine(" -i <input_file>      : Input file.");
	addParamsLine(" -o  <output_file>   : Output file.");
    addUsageLine("Computes the correlation between a set of experimental images with respect "
    		     "to a set of reference images with CUDA in GPU");

}


//#define DEBUG
// Compute distance --------------------------------------------------------
void ProgGpuCorrelation::run()
{

	Image<float> Iref, Iout, mod, ang;
	size_t *Xdim, *Ydim, *Zdim, *Ndim;
	float **original_image_gpu;

	MultidimArray<float> rec_image;
	float **rec_image_gpu;

	MultidimArray<float> module;
	float **module_gpu;

	MultidimArray<float> angle;
	float** angle_gpu;

	//Read input metadataFile
	FileName fnImg, fnImgOut, fnModOut, fnAngOut;
	MDRow rowIn;
	size_t objIndex = 0;

	SF.read(fn_sel,NULL);

	size_t mdInSize = SF.size();
	original_image_gpu = new float *[(int)mdInSize];
	rec_image_gpu = new float *[(int)mdInSize];
	module_gpu = new float *[(int)mdInSize];
	angle_gpu = new float *[(int)mdInSize];

	Xdim = new size_t [(int)mdInSize];
	Ydim = new size_t [(int)mdInSize];
	Zdim = new size_t [(int)mdInSize];
	Ndim = new size_t [(int)mdInSize];

	FOR_ALL_OBJECTS_IN_METADATA(SF){

		++objIndex;
		SF.getRow(rowIn, objIndex);
		rowIn.getValue(MDL_IMAGE, fnImg);
		std::cerr << objIndex << ". Input image: " << fnImg << std::endl;
		Iref.read(fnImg);
		Iref.getDimensions(Xdim[objIndex-1], Ydim[objIndex-1], Zdim[objIndex-1], Ndim[objIndex-1]);
		original_image_gpu[objIndex-1] = new float[Xdim[objIndex-1]*Ydim[objIndex-1]];
		rec_image_gpu[objIndex-1] = new float[Xdim[objIndex-1]*Ydim[objIndex-1]];
		module_gpu[objIndex-1] = new float[Xdim[objIndex-1]*((size_t)(floor(Ydim[objIndex-1]/2)+1))];
		angle_gpu[objIndex-1] = new float[Xdim[objIndex-1]*((size_t)(floor(Ydim[objIndex-1]/2)+1))];
		memcpy(original_image_gpu[objIndex-1], MULTIDIM_ARRAY(Iref()), Xdim[objIndex-1]*Ydim[objIndex-1]*sizeof(float));

		cuda_correlation(original_image_gpu[objIndex-1], rec_image_gpu[objIndex-1], module_gpu[objIndex-1], angle_gpu[objIndex-1], Xdim[objIndex-1], Ydim[objIndex-1], Zdim[objIndex-1]);

		rec_image.coreAllocate(Ndim[objIndex-1], Zdim[objIndex-1], Ydim[objIndex-1], Xdim[objIndex-1]);
		module.coreAllocate(Ndim[objIndex-1], Zdim[objIndex-1], (size_t)(floor(Ydim[objIndex-1]/2)+1), Xdim[objIndex-1]);
		angle.coreAllocate(Ndim[objIndex-1], Zdim[objIndex-1], (size_t)(floor(Ydim[objIndex-1]/2)+1), Xdim[objIndex-1]);
		memcpy(MULTIDIM_ARRAY(rec_image), rec_image_gpu[objIndex-1], Xdim[objIndex-1]*Ydim[objIndex-1]*sizeof(float));
		memcpy(MULTIDIM_ARRAY(module), module_gpu[objIndex-1], Xdim[objIndex-1]*(floor(Ydim[objIndex-1]/2)+1)*sizeof(float));
		memcpy(MULTIDIM_ARRAY(angle), angle_gpu[objIndex-1], Xdim[objIndex-1]*((size_t)(floor(Ydim[objIndex-1]/2)+1))*sizeof(float));
		Iout() = rec_image;
		mod() = module;
		ang() = angle;
		fnImgOut.compose("image", (int)objIndex, "jpg");
		fnModOut.compose("mod", (int)objIndex, "jpg");
		fnAngOut.compose("ang", (int)objIndex, "jpg");
		Iout.write(fnImgOut);
		mod.write(fnModOut);
		ang.write(fnAngOut);
		rec_image.coreDeallocate();
		module.coreDeallocate();
		angle.coreDeallocate();

	}


	for( int i=0 ; i < objIndex; i++){
		delete[] original_image_gpu[i];
		delete[] rec_image_gpu[i];
		delete[] module_gpu[i];
		delete[] angle_gpu[i];
	}
	delete[] original_image_gpu;
	delete[] rec_image_gpu;
	delete[] module_gpu;
	delete[] angle_gpu;
	delete[] Xdim;
	delete[] Ydim;
	delete[] Zdim;
	delete[] Ndim;


}

