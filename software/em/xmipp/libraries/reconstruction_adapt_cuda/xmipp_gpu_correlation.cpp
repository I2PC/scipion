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

#include <data/xmipp_image.h>
#include <data/mask.h>
#include <data/xmipp_fftw.h>
#include <data/transformations.h>
#include <data/metadata_extension.h>

#include <algorithm>
#include "xmipp_gpu_utils.h"
#include "xmipp_gpu_correlation.h"
#include <reconstruction_cuda/cuda_gpu_correlation.h>

void preprocess_images_reference(MetaData &SF, int numImages, Mask &mask,
		GpuCorrelationAux &d_correlationAux, bool rotation)
{
	size_t Xdim, Ydim, Zdim, Ndim;
	getImageSize(SF,Xdim,Ydim,Zdim,Ndim);
	size_t pad_Xdim=2*Xdim-1;
	size_t pad_Ydim=2*Ydim-1;

	d_correlationAux.Xdim=pad_Xdim;
	d_correlationAux.Ydim=pad_Ydim;

	size_t objIndex = 0;
	MDRow rowIn;
	FileName fnImg;
	Image<double> Iref;
	MultidimArray<double> Iref2, padIref, padIref2, padMask;
	padIref.resizeNoCopy(pad_Ydim,pad_Xdim);
	padIref.setXmippOrigin();
	size_t radius=(size_t)mask.R1;

	GpuMultidimArrayAtCpu<double> original_image_stack(Xdim,Ydim,1,numImages);

	/*GpuMultidimArrayAtCpu<double> padded_image_stack;
	GpuMultidimArrayAtCpu<double> padded_image2_stack;
	GpuMultidimArrayAtCpu<double> padded_mask;
	GpuMultidimArrayAtCpu<double> autocorrelation_mask;
	GpuMultidimArrayAtCpu<double> original_image_stack_aux(Xdim,Ydim,1,numImages);

	if(!rotation){
		padded_image_stack.resize(pad_Xdim,pad_Ydim,1,numImages);
		padded_image2_stack.resize(pad_Xdim,pad_Ydim,1,numImages);
		padded_mask.resize(pad_Xdim,pad_Ydim);
		autocorrelation_mask.resize(pad_Xdim,pad_Ydim);
	}*/

	MDIterator *iter = new MDIterator(SF);

	int pointer=0;
	int pointer_pad=0;
	int available_images=numImages;
	size_t n=0;
	while(available_images && iter->objId!=0){

		objIndex = iter->objId;
		available_images--;

		SF.getRow(rowIn, objIndex);
		rowIn.getValue(MDL_IMAGE, fnImg);
		std::cerr << objIndex << ". Image: " << fnImg << std::endl;
		Iref.read(fnImg);
		original_image_stack.fillImage(n,Iref());

		/*Iref().setXmippOrigin();
		mask.apply_mask(Iref(), Iref());
		original_image_stack.fillImage(n,Iref());*/

		/*if(!rotation){
			Iref2=Iref();
			Iref2*=Iref2;
			Iref().window(padIref, STARTINGY(padIref),STARTINGX(padIref),FINISHINGY(padIref),FINISHINGX(padIref));
			Iref2 .window(padIref2,STARTINGY(padIref),STARTINGX(padIref),FINISHINGY(padIref),FINISHINGX(padIref));
			padded_image_stack.fillImage(n,padIref);
			padded_image2_stack.fillImage(n,padIref2);
		}*/

		if(iter->hasNext())
			iter->moveNext();

		n++;
	}
	delete iter;

	/*if(!rotation){
		mask.get_binary_mask().window(padMask, STARTINGY(padIref),STARTINGX(padIref),FINISHINGY(padIref),FINISHINGX(padIref));
		padded_mask.fillImage(0, padMask);
	}*/
	//MultidimArray<double> dMask, maskAutocorrelation;
	//typeCast(padMask, dMask);
	//auto_correlation_matrix(dMask, maskAutocorrelation);
	//autocorrelation_mask.fillImage(0, maskAutocorrelation);



	//AJ new masking and padding
	original_image_stack.copyToGpu(d_correlationAux.d_original_image);
	GpuMultidimArrayAtGpu<double> image_stack_gpu(Xdim,Ydim,1,numImages);
	original_image_stack.copyToGpu(image_stack_gpu);
	MultidimArray<int> maskArray = mask.get_binary_mask();
	MultidimArray<double> dMask;
	typeCast(maskArray, dMask);
	GpuMultidimArrayAtGpu<double> mask_device(Xdim, Ydim, Zdim, 1);
	mask_device.copyToGpu(MULTIDIM_ARRAY(dMask));

	GpuMultidimArrayAtGpu<double> padded_image_gpu, padded_image2_gpu, padded_mask_gpu;
	//if(!rotation){
		padded_image_gpu.resize(pad_Xdim, pad_Ydim, 1, numImages);
		padded_image2_gpu.resize(pad_Xdim, pad_Ydim, 1, numImages);
		padded_mask_gpu.resize(pad_Xdim, pad_Ydim, 1, 1);
	//}

	padding_masking(image_stack_gpu, mask_device, padded_image_gpu, padded_image2_gpu,
			padded_mask_gpu, rotation, false);
	//END AJ//



	//if(!rotation){
		/*GpuMultidimArrayAtGpu<double> padded_image_gpu, padded_image2_gpu, autocorrelation_mask_gpu, padded_mask_gpu;
		padded_image_stack.copyToGpu(padded_image_gpu); //padded_image_gpu
		padded_image2_stack.copyToGpu(padded_image2_gpu);
		padded_mask.copyToGpu(padded_mask_gpu);*/
		//autocorrelation_mask.copyToGpu(d_correlationAux.maskAutocorrelation);
		//FFT
		padded_image_gpu.fft(d_correlationAux.d_projFFT);
		padded_image2_gpu.fft(d_correlationAux.d_projSquaredFFT);
		padded_mask_gpu.fft(d_correlationAux.d_maskFFT);
	//}


	//Polar transform of the projected images
	//GpuMultidimArrayAtGpu<double> original_image_gpu;
	//if(rotation){
		//original_image_stack.copyToGpu(original_image_gpu);
		d_correlationAux.XdimPolar=360;
		d_correlationAux.YdimPolar=radius;
		GpuMultidimArrayAtGpu<double> polar_gpu(360,radius,1,numImages);
		GpuMultidimArrayAtGpu<double> polar2_gpu(360,radius,1,numImages);
		cuda_cart2polar(d_correlationAux.d_original_image, polar_gpu, polar2_gpu, false);
		//FFT
		polar_gpu.fft(d_correlationAux.d_projPolarFFT);
		polar2_gpu.fft(d_correlationAux.d_projPolarSquaredFFT);
	//}

	/*/AJ for debugging
	size_t xAux= pad_Xdim;
	size_t yAux= pad_Ydim;
	size_t nAux= 5;
	//GpuMultidimArrayAtGpu<double> aux(xAux,yAux,1,nAux);
	//d_correlationAux.d_projFFT.ifft(aux);
	GpuMultidimArrayAtCpu<double> auxCpu(xAux,yAux,1,nAux);
	auxCpu.copyFromGpu(padded_image_gpu);
	pointer=0;
	for(int i=0; i<nAux; i++){
	MultidimArray<double> padded;
	FileName fnImgPad;
	Image<double> Ipad;
	padded.coreAllocate(1, 1, yAux, xAux);
	memcpy(MULTIDIM_ARRAY(padded), &auxCpu.data[pointer], xAux*yAux*sizeof(double));
	fnImgPad.compose("paddedRef", i+1, "jpg");
	Ipad()=padded;
	Ipad.write(fnImgPad);
	padded.coreDeallocate();
	pointer += xAux*yAux;
	}
	//END AJ/*/
	/*/AJ for debugging
	xAux= 360;
	yAux= radius;
	//GpuMultidimArrayAtGpu<double> aux(xAux,yAux,1,nAux);
	//d_correlationAux.d_projFFT.ifft(aux);
	GpuMultidimArrayAtCpu<double> auxCpu2(xAux,yAux,1,nAux);
	auxCpu2.copyFromGpu(polar_gpu);
	pointer=0;
	for(int i=0; i<nAux; i++){
	MultidimArray<double> padded2;
	FileName fnImgPad2;
	Image<double> Ipad2;
	padded2.coreAllocate(1, 1, yAux, xAux);
	memcpy(MULTIDIM_ARRAY(padded2), &auxCpu2.data[pointer], xAux*yAux*sizeof(double));
	fnImgPad2.compose("polarRef", i+1, "jpg");
	Ipad2()=padded2;
	Ipad2.write(fnImgPad2);
	padded2.coreDeallocate();
	pointer += xAux*yAux;
	}
	//END AJ/*/

}



void preprocess_images_experimental(MetaData &SF, int numImages, int numImagesRef, Mask &mask,
		GpuCorrelationAux &d_correlationAux, bool rotation)
{
	size_t Xdim, Ydim, Zdim, Ndim;
	getImageSize(SF,Xdim,Ydim,Zdim,Ndim);
	size_t pad_Xdim=2*Xdim-1;
	size_t pad_Ydim=2*Ydim-1;

	d_correlationAux.Xdim=pad_Xdim;
	d_correlationAux.Ydim=pad_Ydim;

	size_t objIndex = 0;
	MDRow rowIn;
	FileName fnImg;
	Image<double> Iref;
	MultidimArray<double> Iref2, padIref, padIref2, padMask;
	padIref.resizeNoCopy(pad_Ydim,pad_Xdim);
	padIref.setXmippOrigin();
	size_t radius=(size_t)mask.R1;

	GpuMultidimArrayAtCpu<double> original_image_stack(Xdim,Ydim,1,numImagesRef);

	/*GpuMultidimArrayAtCpu<double> padded_image_stack;
	GpuMultidimArrayAtCpu<double> padded_image2_stack;
	GpuMultidimArrayAtCpu<double> padded_mask;
	//GpuMultidimArrayAtCpu<double> autocorrelation_mask;
	GpuMultidimArrayAtCpu<double> original_image_stack_aux(Xdim,Ydim,1,numImagesRef);

	if(!rotation){
		padded_image_stack.resize(pad_Xdim,pad_Ydim);
		padded_image2_stack.resize(pad_Xdim,pad_Ydim);
		padded_mask.resize(pad_Xdim,pad_Ydim);
		//autocorrelation_mask.resize(pad_Xdim,pad_Ydim);
	}*/

	MDIterator *iter = new MDIterator(SF);

	int pointer=0;
	int pointer_pad=0;
	int available_images=numImages;

	objIndex = iter->objId;
	//available_images--;

	SF.getRow(rowIn, objIndex);
	rowIn.getValue(MDL_IMAGE, fnImg);
	std::cerr << objIndex << ". Image: " << fnImg << std::endl;
	Iref.read(fnImg);

	for(size_t i=0; i<numImagesRef; i++)
		original_image_stack.fillImage(i,Iref());

	//Iref().setXmippOrigin();

	//AJ new masking and padding
	original_image_stack.copyToGpu(d_correlationAux.d_original_image);
	GpuMultidimArrayAtGpu<double> image_stack_gpu(Xdim,Ydim,1,numImagesRef);
	original_image_stack.copyToGpu(image_stack_gpu);
	MultidimArray<int> maskArray = mask.get_binary_mask();
	MultidimArray<double> dMask;
	typeCast(maskArray, dMask);
	GpuMultidimArrayAtGpu<double> mask_device(Xdim, Ydim, Zdim, 1);
	mask_device.copyToGpu(MULTIDIM_ARRAY(dMask));

	GpuMultidimArrayAtGpu<double> padded_image_gpu, padded_image2_gpu, padded_mask_gpu;
	if(!rotation){
		padded_image_gpu.resize(pad_Xdim, pad_Ydim, 1, numImagesRef);
		padded_image2_gpu.resize(pad_Xdim, pad_Ydim, 1, numImagesRef);
		padded_mask_gpu.resize(pad_Xdim, pad_Ydim, 1, 1);
	}

	padding_masking(image_stack_gpu, mask_device, padded_image_gpu, padded_image2_gpu,
			padded_mask_gpu, rotation, true);
	//END AJ//

	/*
	mask.apply_mask(Iref(), Iref());
	original_image_stack.fillImage(0,Iref());

	Iref().selfReverseX();
	Iref().selfReverseY();

	if(!rotation){
		Iref2=Iref();
		Iref2*=Iref2;
		Iref().window(padIref, STARTINGY(padIref),STARTINGX(padIref),FINISHINGY(padIref),FINISHINGX(padIref));
		Iref2 .window(padIref2,STARTINGY(padIref),STARTINGX(padIref),FINISHINGY(padIref),FINISHINGX(padIref));
		padded_image_stack.fillImage(0,padIref);
		padded_image2_stack.fillImage(0,padIref2);
	}

	if(iter->hasNext())
		iter->moveNext();

	if(!rotation){
		mask.get_binary_mask().window(padMask, STARTINGY(padIref),STARTINGX(padIref),FINISHINGY(padIref),FINISHINGX(padIref));
		padded_mask.fillImage(0, padMask);
	}
	//MultidimArray<double> dMask, maskAutocorrelation;
	//typeCast(padMask, dMask);
	//auto_correlation_matrix(dMask, maskAutocorrelation);
	//autocorrelation_mask.fillImage(0, maskAutocorrelation);
	*/

	//GpuMultidimArrayAtGpu<double> padded_image_gpu, padded_image2_gpu, autocorrelation_mask_gpu, padded_mask_gpu;
	if(!rotation){
		//padded_image_stack.copyToGpuMultiple(padded_image_gpu, numImagesRef);
		//padded_image2_stack.copyToGpuMultiple(padded_image2_gpu, numImagesRef);
		//padded_mask.copyToGpu(padded_mask_gpu);
		//autocorrelation_mask.copyToGpu(d_correlationAux.maskAutocorrelation);
		//FFT
		padded_image_gpu.fft(d_correlationAux.d_projFFT);
		padded_image2_gpu.fft(d_correlationAux.d_projSquaredFFT);
		padded_mask_gpu.fft(d_correlationAux.d_maskFFT);
	}

	//Polar transform of the projected images
	//GpuMultidimArrayAtGpu<double> original_image_gpu;
	if(rotation){
		//original_image_stack.copyToGpuMultiple(original_image_gpu, numImagesRef);
		d_correlationAux.XdimPolar=360;
		d_correlationAux.YdimPolar=radius;
		GpuMultidimArrayAtGpu<double> polar_gpu(360,radius,1,numImagesRef);
		GpuMultidimArrayAtGpu<double> polar2_gpu(360,radius,1,numImagesRef);
		cuda_cart2polar(d_correlationAux.d_original_image, polar_gpu, polar2_gpu, true); // original_image_gpu
		//FFT
		polar_gpu.fft(d_correlationAux.d_projPolarFFT);
		polar2_gpu.fft(d_correlationAux.d_projPolarSquaredFFT);
	}

	/*/AJ for debugging
	size_t xAux= pad_Xdim;
	size_t yAux= pad_Ydim;
	size_t nAux= 1;
	//GpuMultidimArrayAtGpu<double> aux(xAux,yAux,1,nAux);
	//d_correlationAux.d_maskFFT.ifft(aux);
	GpuMultidimArrayAtCpu<double> auxCpu(xAux,yAux,1,nAux);
	auxCpu.copyFromGpu(padded_mask_gpu);
	pointer=0;
	for(int i=0; i<nAux; i++){
	MultidimArray<double> padded;
	FileName fnImgPad;
	Image<double> Ipad;
	padded.coreAllocate(1, 1, yAux, xAux);
	memcpy(MULTIDIM_ARRAY(padded), &auxCpu.data[pointer], xAux*yAux*sizeof(double));
	fnImgPad.compose("paddedMaskBien", i+1, "jpg");
	Ipad()=padded;
	Ipad.write(fnImgPad);
	padded.coreDeallocate();
	pointer += xAux*yAux;
	}
	//END AJ/*/

}


void preprocess_images_experimental_transform(GpuCorrelationAux &d_correlationAux, Mask &mask, bool rotation)
{

	size_t Xdim = d_correlationAux.d_original_image.Xdim;
	size_t Ydim = d_correlationAux.d_original_image.Ydim;
	size_t Zdim = d_correlationAux.d_original_image.Zdim;
	size_t Ndim = d_correlationAux.d_original_image.Ndim;
	size_t pad_Xdim=2*Xdim-1;
	size_t pad_Ydim=2*Ydim-1;
	size_t radius=(size_t)mask.R1;

	MultidimArray<int> maskArray = mask.get_binary_mask();
	MultidimArray<double> dMask;
	typeCast(maskArray, dMask);
	GpuMultidimArrayAtGpu<double> mask_device(Xdim, Ydim, Zdim, 1);
	mask_device.copyToGpu(MULTIDIM_ARRAY(dMask));

	GpuMultidimArrayAtGpu<double> padded_image_gpu, padded_image2_gpu, padded_mask_gpu;
	if(!rotation){
		padded_image_gpu.resize(pad_Xdim, pad_Ydim, 1, Ndim);
		padded_image2_gpu.resize(pad_Xdim, pad_Ydim, 1, Ndim);
		padded_mask_gpu.resize(pad_Xdim, pad_Ydim, 1, 1);
	}

	padding_masking(d_correlationAux.d_transform_image, mask_device, padded_image_gpu, padded_image2_gpu,
			padded_mask_gpu, rotation, true);

	if(!rotation){
		padded_image_gpu.fft(d_correlationAux.d_projFFT);
		padded_image2_gpu.fft(d_correlationAux.d_projSquaredFFT);
		padded_mask_gpu.fft(d_correlationAux.d_maskFFT);
	}

	//Polar transform of the projected images
	GpuMultidimArrayAtGpu<double> polar_gpu, polar2_gpu;
	if(rotation){
		d_correlationAux.XdimPolar=360;
		d_correlationAux.YdimPolar=radius;
		polar_gpu.resize(360,radius,1,Ndim);
		polar2_gpu.resize(360,radius,1,Ndim);
		cuda_cart2polar(d_correlationAux.d_transform_image, polar_gpu, polar2_gpu, true);
		//FFT
		polar_gpu.fft(d_correlationAux.d_projPolarFFT);
		polar2_gpu.fft(d_correlationAux.d_projPolarSquaredFFT);
	}

	/*/AJ for debugging
	size_t xAux= pad_Xdim;
	size_t yAux= pad_Ydim;
	size_t nAux= 5;
	GpuMultidimArrayAtCpu<double> auxCpu(xAux,yAux,1,nAux);
	auxCpu.copyFromGpu(padded_image_gpu);
	int pointer2=0;
	for(int i=0; i<nAux; i++){
	MultidimArray<double> padded;
	FileName fnImgPad;
	Image<double> Ipad;
	padded.coreAllocate(1, 1, yAux, xAux);
	memcpy(MULTIDIM_ARRAY(padded), &auxCpu.data[pointer2], xAux*yAux*sizeof(double));
	fnImgPad.compose("new", i+1, "jpg");
	Ipad()=padded;
	Ipad.write(fnImgPad);
	padded.coreDeallocate();
	pointer2 += xAux*yAux;
	}
	//END AJ/*/

}



// Read arguments ==========================================================
void ProgGpuCorrelation::readParams()
{

    fn_proj = getParam("-i_proj");
    fn_exp = getParam("-i_exp");

}

// Show ====================================================================

void ProgGpuCorrelation::show()
{
    std::cout
	<< "Input projected:          " << fn_proj    << std::endl
	<< "Input experimental:          " << fn_exp    << std::endl
    ;
}

// usage ===================================================================
void ProgGpuCorrelation::defineParams()
{

	addParamsLine(" -i_proj <input_projected_file>      : Input projected images.");
	addParamsLine(" -i_exp  <input_experimental_file>   : Input experimental images.");
    addUsageLine("Computes the correlation between a set of experimental images with respect "
    		     "to a set of reference images with CUDA in GPU");

}

int check_gpu_memory(size_t Xdim, size_t Ydim, int percent){
	float data[3]={0, 0, 0};
	cuda_check_gpu_memory(data);
	int bytes = 8*(2*((2*Xdim)-1)*((2*Ydim)-1) + 2*(360*(Xdim/2)));
	return (int)((data[1]*percent/100)/bytes);
}

// Compute correlation --------------------------------------------------------
void ProgGpuCorrelation::run()
{

	bool rotation = false;

	//PROJECTION IMAGES PART
	printf("Reference images: \n");
	//Read input metadataFile for projection images
	size_t Xdim, Ydim, Zdim, Ndim;
	SF.read(fn_proj,NULL);
	size_t mdInSize = SF.size();
	getImageSize(SF, Xdim, Ydim, Zdim, Ndim);

	// Generate mask
	Mask mask, maskPolar;
    mask.type = BINARY_CIRCULAR_MASK;
	mask.mode = INNER_MASK;
	mask.R1 = std::min(Xdim*0.45, Ydim*0.45);
	mask.resize(Ydim,Xdim);
	mask.get_binary_mask().setXmippOrigin();
	mask.generate_mask();
	int maskCount = mask.get_binary_mask().sum();

	//MultidimArray<double> dMask, maskAutocorrelation;
	//typeCast(mask.get_binary_mask(), dMask);
	//auto_correlation_matrix(dMask, maskAutocorrelation);

	//AJ check_gpu_memory to know how many images we can copy in the gpu memory
	int percent = 70;
	int numImagesProj = check_gpu_memory(Xdim, Ydim, percent);
	printf("%i images can be copied in the GPU memory \n", numImagesProj);
	int available_images_proj;
	if(numImagesProj>mdInSize)
		available_images_proj = mdInSize;
	else
		available_images_proj = numImagesProj-1;

	GpuCorrelationAux d_referenceAux;
	preprocess_images_reference(SF, available_images_proj, mask, d_referenceAux, rotation);
	//d_referenceAux.maskAutocorrelation.resize(XSIZE(maskAutocorrelation),YSIZE(maskAutocorrelation),ZSIZE(maskAutocorrelation),1);
	//fillImage(d_referenceAux.maskAutocorrelation, maskAutocorrelation);
	//if(!rotation){
		d_referenceAux.maskCount=maskCount;
		d_referenceAux.produceSideInfo();
	//}

	//EXPERIMENTAL IMAGES PART
	printf("Experimental images: \n");
	SFexp.read(fn_exp,NULL);
	size_t mdExpSize = SFexp.size();

	int numImagesExp = check_gpu_memory(Xdim, Ydim, percent);
	printf("%i experimental images can be copied in the GPU memory \n", numImagesExp);
	int available_images_exp;
	if(numImagesExp>mdExpSize)
		available_images_exp = mdExpSize;
	else
		available_images_exp = numImagesExp-1;

	GpuCorrelationAux d_experimentalAux;
	preprocess_images_experimental(SFexp, 1, available_images_proj, mask, d_experimentalAux, rotation);
	//d_experimentalAux.maskAutocorrelation.resize(XSIZE(maskAutocorrelation),YSIZE(maskAutocorrelation),ZSIZE(maskAutocorrelation),1);
	//fillImage(d_experimentalAux.maskAutocorrelation, maskAutocorrelation);
	if(!rotation){
		d_experimentalAux.maskCount=maskCount;
		d_experimentalAux.produceSideInfo();
	}


	//CORRELATION PART
	//TRANSFORMATION MATRIX CALCULATION
	printf("Calculating correlation...\n");
	TransformMatrix<float> transMat(available_images_proj);
	if(!rotation)
		cuda_calculate_correlation(d_referenceAux, d_experimentalAux, transMat);
	else
		cuda_calculate_correlation_rotation(d_referenceAux, d_experimentalAux, transMat);


	//APPLY TRANSFORMATION
	printf("Applying transformation...\n");
	d_experimentalAux.d_transform_image.resize(d_experimentalAux.d_original_image);
	apply_transform(d_experimentalAux.d_original_image, d_experimentalAux.d_transform_image, transMat);

	/*/AJ for debugging
	size_t xAux= Xdim;
	size_t yAux= Ydim;
	size_t nAux= available_images_proj;
	GpuMultidimArrayAtCpu<double> auxCpu(xAux,yAux,1,nAux);
	auxCpu.copyFromGpu(d_experimentalAux.d_transform_image);
	int pointer2=0;
	for(int i=0; i<nAux; i++){
	MultidimArray<double> padded;
	FileName fnImgPad;
	Image<double> Ipad;
	padded.coreAllocate(1, 1, yAux, xAux);
	memcpy(MULTIDIM_ARRAY(padded), &auxCpu.data[pointer2], xAux*yAux*sizeof(double));
	fnImgPad.compose("check", i+1, "jpg");
	Ipad()=padded;
	Ipad.write(fnImgPad);
	padded.coreDeallocate();
	pointer2 += xAux*yAux;
	}
	//END AJ/*/

	//PREPROCESS TO PREPARE DATA TO THE NEXT STEP
	printf("Processing experimental images again...\n");
	rotation = !rotation;
	preprocess_images_experimental_transform(d_experimentalAux, mask, rotation);
	if(!rotation){
		d_experimentalAux.maskCount=maskCount;
		d_experimentalAux.produceSideInfo();
	}

	//CORRELATION PART
	//TRANSFORMATION MATRIX CALCULATION
	printf("Calculating correlation...\n");
	if(!rotation)
		cuda_calculate_correlation(d_referenceAux, d_experimentalAux, transMat);
	else
		cuda_calculate_correlation_rotation(d_referenceAux, d_experimentalAux, transMat);

	/*/AJ for debugging
	size_t xAux= d_experimentalAux.MFrealSpace.Xdim;
	size_t yAux= d_experimentalAux.MFrealSpace.Ydim;
	size_t nAux= 5;
	GpuMultidimArrayAtCpu<double> auxCpu(xAux,yAux,1,nAux);
	auxCpu.copyFromGpu(d_experimentalAux.MFrealSpace);
	int pointer2=0;
	for(int i=0; i<nAux; i++){
	MultidimArray<double> padded;
	FileName fnImgPad;
	Image<double> Ipad;
	padded.coreAllocate(1, 1, yAux, xAux);
	memcpy(MULTIDIM_ARRAY(padded), &auxCpu.data[pointer2], xAux*yAux*sizeof(double));
	fnImgPad.compose("otro", i+1, "jpg");
	Ipad()=padded;
	Ipad.write(fnImgPad);
	padded.coreDeallocate();
	pointer2 += xAux*yAux;
	}
	//END AJ/*/



}

