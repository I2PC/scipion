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


//CPU functions
void padding(const float* image, float* padded_image, size_t Xdim, size_t Ydim, size_t new_Xdim, size_t new_Ydim)
{
    int im_size = Xdim*Ydim;
    //int new_Xdim = ((Xdim*2)-1);
    //int new_Ydim = ((Ydim*2)-1);

    // Pad image
    memset(padded_image, 0, new_Xdim*new_Ydim * sizeof(float));
    int pointer=(((new_Ydim-Ydim)/2)*new_Xdim)+((new_Xdim-Xdim)/2);
    for (int i=0;i<Ydim;i++){
    	memcpy(&padded_image[pointer], &image[i*Ydim], Ydim*sizeof(float));
	pointer+=new_Xdim;
    }
}



//CUDA functions
void cuda_fft_main(float *image, float *fft_image, float *module, float *angle, size_t Xdim, size_t Ydim, size_t Zdim, size_t batch){

	cuda_fft(image, fft_image, module, angle, Xdim, Ydim, Zdim, batch);


}

void cuda_polar_main(float *image, float *polar, size_t Xdim, size_t Ydim, size_t Zdim, size_t batch){


	cuda_cart2polar(image, polar, Xdim, Ydim, Zdim, batch);

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

	Image<float> Iref, Iout, mod, ang, Ipad, Ipolar;
	size_t Xdim, Ydim, Zdim, Ndim;
	float *original_image_gpu;

	MultidimArray<float> rec_image;
	float *rec_image_gpu;

	MultidimArray<float> module;
	float *module_gpu;

	MultidimArray<float> angle;
	float *angle_gpu;
	
	MultidimArray<float> polar;
	float *polar_gpu;

	MultidimArray<float> padded;
	float *padded_image_gpu;

	//Read input metadataFile
	FileName fnImg, fnImgOut, fnModOut, fnAngOut, fnImgPad, fnPolar;
	MDRow rowIn;
	size_t objIndex = 0;

	SF.read(fn_sel,NULL);

	size_t mdInSize = SF.size();

	size_t pad_Xdim, pad_Ydim;

	getImageSize(SF, Xdim, Ydim, Zdim, Ndim);
	int radio = std::max(Xdim/2, Ydim/2);
	int totalRealDim = Xdim*Ydim*Zdim*mdInSize;
	int totalPolar = radio*360*mdInSize;

	pad_Xdim = ((Xdim*2)-1);
	pad_Ydim = ((Ydim*2)-1);
	int totalPad = pad_Xdim*pad_Ydim*mdInSize;
	int totalComplexDim = pad_Xdim*((pad_Ydim/2)+1)*Zdim*mdInSize;

	original_image_gpu = new float[totalRealDim];
	rec_image_gpu = new float[totalPad];
	module_gpu = new float[totalComplexDim];
	angle_gpu = new float[totalComplexDim];
	polar_gpu = new float[totalPolar];
	padded_image_gpu = new float[totalPad];

	int pointer=0;
	int pointer_pad=0;
	FOR_ALL_OBJECTS_IN_METADATA(SF){

		++objIndex;
		SF.getRow(rowIn, objIndex);
		rowIn.getValue(MDL_IMAGE, fnImg);
		std::cerr << objIndex << ". Input image: " << fnImg << std::endl;
		Iref.read(fnImg);

		memcpy(&original_image_gpu[pointer], MULTIDIM_ARRAY(Iref()), Xdim*Ydim*sizeof(float));
		pointer += Xdim*Ydim;

		//AJ padding
		padded.coreAllocate(Ndim, Zdim, pad_Ydim, pad_Xdim);
		padding(MULTIDIM_ARRAY(Iref()), &padded_image_gpu[pointer_pad], Xdim, Ydim, pad_Xdim, pad_Ydim);
		memcpy(MULTIDIM_ARRAY(padded), &padded_image_gpu[pointer_pad], pad_Xdim*pad_Ydim*sizeof(float));
		fnImgPad.compose("pad", (int)objIndex, "jpg");
		Ipad()=padded;
		Ipad.write(fnImgPad);
		padded.coreDeallocate();
		pointer_pad += pad_Xdim*pad_Ydim;
		//end padding
	}

	cuda_polar_main(original_image_gpu, polar_gpu, Xdim, Ydim, Zdim, mdInSize);
	cuda_fft_main(padded_image_gpu, rec_image_gpu, module_gpu, angle_gpu, pad_Xdim, pad_Ydim, Zdim, mdInSize);

	objIndex=0;

	std::cerr << "Writing output data" << std::endl;
	int pointerPad=0;
	int pointerComplex=0;
	int pointerPolar=0;
	FOR_ALL_OBJECTS_IN_METADATA(SF){

		++objIndex;
		rec_image.coreAllocate(Ndim, Zdim, pad_Ydim, pad_Xdim);
		module.coreAllocate(Ndim, Zdim, ((pad_Ydim/2)+1), pad_Xdim);
		angle.coreAllocate(Ndim, Zdim, ((pad_Ydim/2)+1), pad_Xdim);
		polar.coreAllocate(1, 1, radio, 360);
		memcpy(MULTIDIM_ARRAY(rec_image), &rec_image_gpu[pointerPad], pad_Xdim*pad_Ydim*sizeof(float));
		memcpy(MULTIDIM_ARRAY(module), &module_gpu[pointerComplex], pad_Xdim*((pad_Ydim/2)+1)*sizeof(float));
		memcpy(MULTIDIM_ARRAY(angle), &angle_gpu[pointerComplex], pad_Xdim*((pad_Ydim/2)+1)*sizeof(float));
		memcpy(MULTIDIM_ARRAY(polar), &polar_gpu[pointerPolar], radio*360*sizeof(float));
		pointerPad += pad_Xdim*pad_Ydim;
		pointerComplex += pad_Xdim*((pad_Ydim/2)+1);
		pointerPolar += radio*360;
		Iout() = rec_image;
		mod() = module;
		ang() = angle;
		Ipolar() = polar;
		fnImgOut.compose("rec", (int)objIndex, "jpg");
		fnModOut.compose("mod", (int)objIndex, "jpg");
		fnAngOut.compose("ang", (int)objIndex, "jpg");
		fnPolar.compose("polar", (int)objIndex, "jpg");
		Iout.write(fnImgOut);
		mod.write(fnModOut);
		ang.write(fnAngOut);
		Ipolar.write(fnPolar);
		rec_image.coreDeallocate();
		module.coreDeallocate();
		angle.coreDeallocate();
		polar.coreDeallocate();

	}


	delete[] original_image_gpu;
	delete[] rec_image_gpu;
	delete[] module_gpu;
	delete[] angle_gpu;
	delete[] padded_image_gpu;
	delete[] polar_gpu;
}

