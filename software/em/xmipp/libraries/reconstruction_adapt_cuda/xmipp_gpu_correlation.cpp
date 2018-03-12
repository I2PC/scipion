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

#include "xmipp_gpu_correlation.h"

#include <data/xmipp_image.h>
#include <data/mask.h>
#include <data/xmipp_fftw.h>
#include <data/transformations.h>
#include <data/metadata_extension.h>
#include <data/filters.h>
#include <data/xmipp_funcs.h>
//#include <data/xmipp_threads.h>

#include "xmipp_gpu_utils.h"
#include <reconstruction_cuda/cuda_gpu_correlation.h>

#include <algorithm>
#include <math.h>
#include <time.h>
#include <sys/time.h>




// A function to print all prime factors of a given number n
void primeFactors(int n, int *out)
{
	int n_orig = n;
	// Print the number of 2s that divide n
	while (n%2 == 0)
	{
		//printf("%d ", 2);
		out[0]++;
		n = n/2;
	}

	// n must be odd at this point. So we can skip
	// one element (Note i = i +2)
	for (int i = 3; i <= sqrt(n_orig); i = i+2)
	{
		// While i divides n, print i and divide n
		while (n%i == 0)
		{
			//printf("%d ", i);
			if (i==3)
				out[1]++;
			else if (i==5)
				out[2]++;
			else if (i==7)
				out[3]++;
			else if(i>7)
				out[4]++;

			n = n/i;
		}
	}

	// This condition is to handle the case when n
	// is a prime number greater than 2
	if (n > 2){
		//printf ("%d ", n);
		out[4]++;
	}
}


void preprocess_images_reference(MetaData &SF, int firstIdx, int numImages, Mask &mask, GpuCorrelationAux &d_correlationAux,
		mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask, mycufftHandle &myhandlePolar, mycufftHandle &myhandleAux,
		StructuresAux &myStructureAux, MDIterator *iter, myStreamHandle myStream)
{
	size_t Xdim, Ydim, Zdim, Ndim;
	getImageSize(SF,Xdim,Ydim,Zdim,Ndim);
	size_t pad_Xdim=d_correlationAux.Xdim;
	size_t pad_Ydim=d_correlationAux.Ydim;

	MDRow rowIn;
	FileName fnImg;
	Image<float> Iref;
	size_t radius = d_correlationAux.YdimPolar;
	size_t angles = d_correlationAux.XdimPolar;

	GpuMultidimArrayAtCpu<float> original_image_stack_ref(Xdim,Ydim,1,numImages);

	//MDIterator *iter = new MDIterator(SF);

	size_t n=0;
	for(int i=firstIdx; i<firstIdx+numImages; i++){

		SF.getValue(MDL_IMAGE,fnImg,iter->objId);
		//std::cerr << iter->objId << ". Image: " << fnImg << std::endl;
		Iref.read(fnImg);
		original_image_stack_ref.fillImage(n,Iref()/8);

		if(iter->hasNext())
			iter->moveNext();

		n++;
	}

	//delete iter;

	GpuMultidimArrayAtGpu<float> image_stack_gpu(Xdim,Ydim,1,numImages);
	original_image_stack_ref.copyToGpu(image_stack_gpu, myStream);

	MultidimArray<int> maskArray = mask.get_binary_mask();
	MultidimArray<float> dMask;
	typeCast(maskArray, dMask);
	d_correlationAux.d_mask.resize(Xdim, Ydim, Zdim, 1);
	float *mask_aux;
	cpuMalloc((void**)&mask_aux, sizeof(float)*Xdim*Ydim*Zdim);
	memcpy(mask_aux, MULTIDIM_ARRAY(dMask), sizeof(float)*Xdim*Ydim*Zdim);
	d_correlationAux.d_mask.copyToGpuStream(mask_aux, myStream);

	padding_masking(image_stack_gpu, d_correlationAux.d_mask, myStructureAux.padded_image_gpu, myStructureAux.padded_image2_gpu,
			myStructureAux.padded_mask_gpu, false, myStream);

	GpuMultidimArrayAtGpu< std::complex<float> > dull;

    myStructureAux.padded_image_gpu.fftStream(d_correlationAux.d_projFFT, myhandlePadded, myStream, false, dull);

    myStructureAux.padded_image2_gpu.fftStream(d_correlationAux.d_projSquaredFFT, myhandlePadded, myStream, false, dull);
    myStructureAux.padded_mask_gpu.fftStream(d_correlationAux.d_maskFFT, myhandleMask, myStream, false, dull);

	//Polar transform of the projected images
	cuda_cart2polar(image_stack_gpu, myStructureAux.polar_gpu, myStructureAux.polar2_gpu, false, myStream);

    myStructureAux.polar_gpu.fftStream(d_correlationAux.d_projPolarFFT, myhandlePolar, myStream, false, dull);

    myStructureAux.polar2_gpu.fftStream(d_correlationAux.d_projPolarSquaredFFT, myhandlePolar, myStream, false, dull);

}



void preprocess_images_experimental(MetaData &SF, FileName &fnImg, int numImagesRef, GpuMultidimArrayAtGpu<float> &mask,
		GpuMultidimArrayAtGpu< std::complex<float> > &d_maskFFT, GpuCorrelationAux &d_correlationAux, bool rotation,
		int firstStep, bool mirror, mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask, mycufftHandle &myhandlePolar,
		StructuresAux &myStructureAux, myStreamHandle myStream)
{
	size_t Xdim, Ydim, Zdim, Ndim;
	getImageSize(SF,Xdim,Ydim,Zdim,Ndim);
	size_t pad_Xdim=d_correlationAux.Xdim;
	size_t pad_Ydim=d_correlationAux.Ydim;
	size_t radius=d_correlationAux.YdimPolar;
	size_t angles = d_correlationAux.XdimPolar;

	GpuMultidimArrayAtCpu<float> original_image_stack(Xdim,Ydim,1,numImagesRef);

	GpuMultidimArrayAtGpu< std::complex<float> > dull;

	if(firstStep==0){

		Image<float> Iref;

		Iref.read(fnImg);

		//AJ mirror of the image
		if(mirror)
			Iref().selfReverseX();
		//END AJ mirror

		for(size_t i=0; i<numImagesRef; i++)
			original_image_stack.fillImage(i,Iref()/8);

	}

		original_image_stack.copyToGpu(d_correlationAux.d_original_image, myStream);

	if(!rotation){
		padding_masking(d_correlationAux.d_original_image, mask, myStructureAux.padded_image_gpu, myStructureAux.padded_image2_gpu,
				myStructureAux.padded_mask_gpu, true, myStream);

		myStructureAux.padded_image_gpu.fftStream(d_correlationAux.d_projFFT, myhandlePadded, myStream, false, dull);

	    myStructureAux.padded_image2_gpu.fftStream(d_correlationAux.d_projSquaredFFT, myhandlePadded, myStream, false, dull);
		d_maskFFT.copyGpuToGpuStream(d_correlationAux.d_maskFFT, myStream);

	}

	if(rotation){
		cuda_cart2polar(d_correlationAux.d_original_image, myStructureAux.polar_gpu, myStructureAux.polar2_gpu, true, myStream);
	    myStructureAux.polar_gpu.fftStream(d_correlationAux.d_projPolarFFT, myhandlePolar, myStream, false, dull);
	    myStructureAux.polar2_gpu.fftStream(d_correlationAux.d_projPolarSquaredFFT, myhandlePolar, myStream, false, dull);
	}

}



void preprocess_images_experimental_two(MetaData &SF, FileName &fnImg, int numImagesRef, GpuMultidimArrayAtGpu<float> &mask,
		GpuMultidimArrayAtGpu< std::complex<float> > &d_maskFFT,
		GpuCorrelationAux &d_correlationAuxTR, GpuCorrelationAux &d_correlationAuxRT,
		bool rotation, int firstStep, bool mirror,
		mycufftHandle &myhandlePaddedTR, mycufftHandle &myhandleMaskTR, mycufftHandle &myhandlePolarTR,
		mycufftHandle &myhandlePaddedRT, mycufftHandle &myhandleMaskRT, mycufftHandle &myhandlePolarRT,
		StructuresAux &myStructureAuxTR, StructuresAux &myStructureAuxRT,
		myStreamHandle &myStreamTR, myStreamHandle &myStreamRT,
		GpuMultidimArrayAtCpu<float> &original_image_stack,
		mycufftHandle &ifftcb)
{


	size_t Xdim, Ydim, Zdim, Ndim;
	getImageSize(SF,Xdim,Ydim,Zdim,Ndim);
	size_t pad_Xdim=d_correlationAuxTR.Xdim;
	size_t pad_Ydim=d_correlationAuxTR.Ydim;
	size_t radius=d_correlationAuxTR.YdimPolar;
	size_t angles = d_correlationAuxTR.XdimPolar;

	original_image_stack.resize(Xdim,Ydim,1,numImagesRef);

	if(firstStep==0){

		Image<float> Iref;

		Iref.read(fnImg);

		//AJ mirror of the image
		if(mirror)
			Iref().selfReverseX();
		//END AJ mirror

		for(size_t i=0; i<numImagesRef; i++)
			original_image_stack.fillImage(i,Iref()/8);

	}

		d_correlationAuxTR.d_original_image.resize(Xdim,Ydim,1,numImagesRef);
		d_correlationAuxRT.d_original_image.resize(Xdim,Ydim,1,numImagesRef);
		d_correlationAuxTR.d_projFFT.resize((pad_Xdim/2)+1, pad_Ydim, 1, numImagesRef);
		d_correlationAuxTR.d_projSquaredFFT.resize((pad_Xdim/2)+1, pad_Ydim, 1, numImagesRef);
		d_correlationAuxRT.d_projPolarFFT.resize((angles/2)+1, radius, 1, numImagesRef);
		d_correlationAuxRT.d_projPolarSquaredFFT.resize((angles/2)+1, radius, 1, numImagesRef);
		d_correlationAuxTR.d_maskFFT.resize(d_maskFFT);

		original_image_stack.copyToGpu(d_correlationAuxTR.d_original_image, myStreamTR);

		padding_masking(d_correlationAuxTR.d_original_image, mask, myStructureAuxTR.padded_image_gpu, myStructureAuxTR.padded_image2_gpu,
				myStructureAuxTR.padded_mask_gpu, true, myStreamTR);

		original_image_stack.copyToGpu(d_correlationAuxRT.d_original_image, myStreamRT);

		cuda_cart2polar(d_correlationAuxRT.d_original_image, myStructureAuxRT.polar_gpu, myStructureAuxRT.polar2_gpu, true, myStreamRT);


		GpuMultidimArrayAtGpu< std::complex<float> > dull;

		myStructureAuxTR.padded_image_gpu.fftStream(d_correlationAuxTR.d_projFFT, myhandlePaddedTR, myStreamTR, false, dull);

		myStructureAuxTR.padded_image2_gpu.fftStream(d_correlationAuxTR.d_projSquaredFFT, myhandlePaddedTR, myStreamTR, false, dull);
		d_maskFFT.copyGpuToGpuStream(d_correlationAuxTR.d_maskFFT, myStreamTR);

	    myStructureAuxRT.polar_gpu.fftStream(d_correlationAuxRT.d_projPolarFFT, myhandlePolarRT, myStreamRT, false, dull);
	    myStructureAuxRT.polar2_gpu.fftStream(d_correlationAuxRT.d_projPolarSquaredFFT, myhandlePolarRT, myStreamRT, false, dull);

}



void preprocess_images_experimental_transform_two(MetaData &SF, FileName &fnImg, int numImagesRef, GpuMultidimArrayAtGpu<float> &mask,
		GpuMultidimArrayAtGpu< std::complex<float> > &d_maskFFT,
		GpuCorrelationAux &d_correlationAuxOne, GpuCorrelationAux &d_correlationAuxTwo,
		bool rotation, int firstStep, bool mirror,
		mycufftHandle &myhandlePaddedOne, mycufftHandle &myhandleMaskOne, mycufftHandle &myhandlePolarOne,
		mycufftHandle &myhandlePaddedTwo, mycufftHandle &myhandleMaskTwo, mycufftHandle &myhandlePolarTwo,
		StructuresAux &myStructureAuxOne, StructuresAux &myStructureAuxTwo,
		myStreamHandle &myStreamOne, myStreamHandle &myStreamTwo, int step)
{

	size_t Xdim = d_correlationAuxOne.d_transform_image.Xdim;
	size_t Ydim = d_correlationAuxOne.d_transform_image.Ydim;
	size_t Zdim = d_correlationAuxOne.d_transform_image.Zdim;
	size_t Ndim = d_correlationAuxOne.d_transform_image.Ndim;
	size_t pad_Xdim=d_correlationAuxOne.Xdim;
	size_t pad_Ydim=d_correlationAuxOne.Ydim;
	size_t radius=d_correlationAuxOne.YdimPolar;
	size_t angles = d_correlationAuxOne.XdimPolar;

	d_correlationAuxOne.d_projFFT.resize((pad_Xdim/2)+1, pad_Ydim, 1, numImagesRef);
	d_correlationAuxOne.d_projSquaredFFT.resize((pad_Xdim/2)+1, pad_Ydim, 1, numImagesRef);
	d_correlationAuxTwo.d_projPolarFFT.resize((angles/2)+1, radius, 1, numImagesRef);
	d_correlationAuxTwo.d_projPolarSquaredFFT.resize((angles/2)+1, radius, 1, numImagesRef);
	d_correlationAuxOne.d_maskFFT.resize(d_maskFFT);


	padding_masking(d_correlationAuxOne.d_transform_image, mask, myStructureAuxOne.padded_image_gpu, myStructureAuxOne.padded_image2_gpu,
			myStructureAuxOne.padded_mask_gpu, true, myStreamOne);

	cuda_cart2polar(d_correlationAuxTwo.d_transform_image, myStructureAuxTwo.polar_gpu, myStructureAuxTwo.polar2_gpu, true, myStreamTwo);

	GpuMultidimArrayAtGpu< std::complex<float> > dull;
	myStructureAuxOne.padded_image_gpu.fftStream(d_correlationAuxOne.d_projFFT, myhandlePaddedOne, myStreamOne, false, dull);
	myStructureAuxOne.padded_image2_gpu.fftStream(d_correlationAuxOne.d_projSquaredFFT, myhandlePaddedOne, myStreamOne, false, dull);
	d_maskFFT.copyGpuToGpuStream(d_correlationAuxOne.d_maskFFT, myStreamOne);

	myStructureAuxTwo.polar_gpu.fftStream(d_correlationAuxTwo.d_projPolarFFT, myhandlePolarTwo, myStreamTwo, false, dull);
	myStructureAuxTwo.polar2_gpu.fftStream(d_correlationAuxTwo.d_projPolarSquaredFFT, myhandlePolarTwo, myStreamTwo, false, dull);

}


void preprocess_images_experimental_transform(GpuCorrelationAux &d_correlationAux, GpuMultidimArrayAtGpu<float> &mask,
		GpuMultidimArrayAtGpu< std::complex<float> > &d_maskFFT, bool rotation, int step, mycufftHandle &myhandlePadded,
		mycufftHandle &myhandleMask, mycufftHandle &myhandlePolar, StructuresAux &myStructureAux, myStreamHandle myStream)
{

	size_t Xdim = d_correlationAux.d_transform_image.Xdim;
	size_t Ydim = d_correlationAux.d_transform_image.Ydim;
	size_t Zdim = d_correlationAux.d_transform_image.Zdim;
	size_t Ndim = d_correlationAux.d_transform_image.Ndim;
	size_t pad_Xdim=d_correlationAux.Xdim;
	size_t pad_Ydim=d_correlationAux.Ydim;
	size_t radius=d_correlationAux.YdimPolar;
	size_t angles = d_correlationAux.XdimPolar;

	GpuMultidimArrayAtGpu< std::complex<float> > dull;

	if(!rotation){
		padding_masking(d_correlationAux.d_transform_image, mask, myStructureAux.padded_image_gpu, myStructureAux.padded_image2_gpu,
				myStructureAux.padded_mask_gpu, true, myStream);

		myStructureAux.padded_image_gpu.fftStream(d_correlationAux.d_projFFT, myhandlePadded, myStream, false, dull);

	    myStructureAux.padded_image2_gpu.fftStream(d_correlationAux.d_projSquaredFFT, myhandlePadded, myStream, false, dull);
		d_maskFFT.copyGpuToGpuStream(d_correlationAux.d_maskFFT, myStream);

	}

	//Polar transform of the projected images
	if(rotation){
		cuda_cart2polar(d_correlationAux.d_transform_image, myStructureAux.polar_gpu, myStructureAux.polar2_gpu, true, myStream);
		myStructureAux.polar_gpu.fftStream(d_correlationAux.d_projPolarFFT, myhandlePolar, myStream, false, dull);
		myStructureAux.polar2_gpu.fftStream(d_correlationAux.d_projPolarSquaredFFT, myhandlePolar, myStream, false, dull);

	}


}

void align_experimental_image(FileName &fnImgExp, GpuCorrelationAux &d_referenceAux,
		GpuCorrelationAux &d_experimentalAuxTR, GpuCorrelationAux &d_experimentalAuxRT,
		TransformMatrix<float> &transMat_tr, TransformMatrix<float> &transMat_rt, float *max_vector_tr, float *max_vector_rt,
		MetaData &SFexp, int available_images_proj, bool mirror, int maxShift,
		mycufftHandle &myhandlePadded_tr, mycufftHandle &myhandleMask_tr, mycufftHandle &myhandlePolar_tr,
		mycufftHandle &myhandlePaddedB_tr, mycufftHandle &myhandleMaskB_tr, mycufftHandle &myhandlePolarB_tr,
		mycufftHandle &myhandlePadded_rt, mycufftHandle &myhandleMask_rt, mycufftHandle &myhandlePolar_rt,
		mycufftHandle &myhandlePaddedB_rt, mycufftHandle &myhandleMaskB_rt, mycufftHandle &myhandlePolarB_rt,
		StructuresAux &myStructureAux_tr, StructuresAux &myStructureAux_rt,
		myStreamHandle &myStreamTR, myStreamHandle &myStreamRT,
		TransformMatrix<float> &resultTR, TransformMatrix<float> &resultRT,
		GpuMultidimArrayAtCpu<float> &original_image_stack, mycufftHandle &ifftcb)
{

	bool rotation;

	////////////////////////
	//FIRST PART FOR TRTRTRT

	int max_step;
	rotation = false;
	//max_vector = max_vector_tr;
	max_step=7;


	preprocess_images_experimental_two(SFexp, fnImgExp, available_images_proj, d_referenceAux.d_mask,
			d_referenceAux.d_maskFFT, d_experimentalAuxTR, d_experimentalAuxRT, true, 0, mirror,
					myhandlePadded_tr, myhandleMask_tr, myhandlePolar_tr,
					myhandlePadded_rt, myhandleMask_rt, myhandlePolar_rt,
					myStructureAux_tr, myStructureAux_rt, myStreamTR, myStreamRT, original_image_stack,
					ifftcb);

	d_experimentalAuxTR.maskCount=d_referenceAux.maskCount;
	d_experimentalAuxTR.produceSideInfo(myhandlePaddedB_tr, myhandleMaskB_tr, myStructureAux_tr,
			d_referenceAux.maskAutocorrelation, myStreamTR);

	d_experimentalAuxTR.d_transform_image.resize(d_experimentalAuxTR.d_original_image);
	d_experimentalAuxRT.d_transform_image.resize(d_experimentalAuxRT.d_original_image);

	//transMat = &transMat_tr;

	for(int step=0; step<6; step++){

		bool saveMaxVector = false;
		if(step==5)
			saveMaxVector = true;

		if(step%2==0){

			//FIRST TRANSLATION AND SECOND ROTATION
			//CORRELATION PART
			//TRANSFORMATION MATRIX CALCULATION
			cuda_calculate_correlation_two(d_referenceAux, d_experimentalAuxTR,
					transMat_tr, max_vector_tr, maxShift,
					myhandlePaddedB_tr, mirror, myStructureAux_tr,
					myStreamTR,
					d_experimentalAuxRT, transMat_rt,
					max_vector_rt, myhandlePolarB_rt,
					myStructureAux_rt, myStreamRT,
					resultTR, resultRT, ifftcb, saveMaxVector);

			//APPLY TRANSFORMATION
			apply_transform(d_experimentalAuxTR.d_original_image, d_experimentalAuxTR.d_transform_image, transMat_tr, myStreamTR);

			apply_transform(d_experimentalAuxRT.d_original_image, d_experimentalAuxRT.d_transform_image, transMat_rt, myStreamRT);

			//PREPROCESS TO PREPARE DATA TO THE NEXT STEP
			preprocess_images_experimental_transform_two(SFexp, fnImgExp, available_images_proj, d_referenceAux.d_mask,
					d_referenceAux.d_maskFFT, d_experimentalAuxRT, d_experimentalAuxTR,	true, 0, mirror,
					myhandlePadded_rt, myhandleMask_rt, myhandlePolar_rt,
					myhandlePadded_tr, myhandleMask_tr, myhandlePolar_tr,
					myStructureAux_rt, myStructureAux_tr, myStreamRT, myStreamTR, 1);

			d_experimentalAuxRT.maskCount=d_referenceAux.maskCount;
			d_experimentalAuxRT.produceSideInfo(myhandlePaddedB_rt, myhandleMaskB_rt, myStructureAux_rt,
					d_referenceAux.maskAutocorrelation, myStreamRT);

		}
		else{

			//FIRST ROTATION AND SECOND TRANSLATION
			//CORRELATION PART
			//TRANSFORMATION MATRIX CALCULATION
			cuda_calculate_correlation_two(d_referenceAux, d_experimentalAuxRT,
					transMat_rt, max_vector_rt, maxShift,
					myhandlePaddedB_rt, mirror, myStructureAux_rt,
					myStreamRT,
					d_experimentalAuxTR, transMat_tr,
					max_vector_tr, myhandlePolarB_tr,
					myStructureAux_tr, myStreamTR,
					resultRT, resultTR, ifftcb, saveMaxVector);


			if(step < 5){

				//APPLY TRANSFORMATION
				apply_transform(d_experimentalAuxRT.d_original_image, d_experimentalAuxRT.d_transform_image, transMat_rt, myStreamRT);

				apply_transform(d_experimentalAuxTR.d_original_image, d_experimentalAuxTR.d_transform_image, transMat_tr, myStreamTR);

				//PREPROCESS TO PREPARE DATA TO THE NEXT STEP
				preprocess_images_experimental_transform_two(SFexp, fnImgExp, available_images_proj, d_referenceAux.d_mask,
						d_referenceAux.d_maskFFT, d_experimentalAuxTR, d_experimentalAuxRT,	true, 0, mirror,
						myhandlePadded_tr, myhandleMask_tr, myhandlePolar_tr,
						myhandlePadded_rt, myhandleMask_rt, myhandlePolar_rt,
						myStructureAux_tr, myStructureAux_rt, myStreamTR, myStreamRT, 2);

				d_experimentalAuxTR.maskCount=d_referenceAux.maskCount;
				d_experimentalAuxTR.produceSideInfo(myhandlePaddedB_tr, myhandleMaskB_tr, myStructureAux_tr,
						d_referenceAux.maskAutocorrelation, myStreamTR);

			}else if(step==5){

				//APPLY TRANSFORMATION
				d_experimentalAuxTR.d_transform_image.resize(d_experimentalAuxTR.d_original_image);
				apply_transform(d_experimentalAuxTR.d_original_image, d_experimentalAuxTR.d_transform_image, transMat_tr, myStreamTR);

				//PREPROCESS TO PREPARE DATA TO THE NEXT STEP
				preprocess_images_experimental_transform(d_experimentalAuxTR, d_referenceAux.d_mask, d_referenceAux.d_maskFFT, false, step,
							myhandlePadded_tr, myhandleMask_tr, myhandlePolar_tr, myStructureAux_tr, myStreamTR);
				d_experimentalAuxTR.maskCount=d_referenceAux.maskCount;
				d_experimentalAuxTR.produceSideInfo(myhandlePaddedB_tr, myhandleMaskB_tr, myStructureAux_tr,
						d_referenceAux.maskAutocorrelation, myStreamTR);

				//CORRELATION PART
				//TRANSFORMATION MATRIX CALCULATION
				cuda_calculate_correlation(d_referenceAux, d_experimentalAuxTR, transMat_tr, max_vector_tr, maxShift, myhandlePaddedB_tr,
							mirror, myStructureAux_tr, myStreamTR, resultTR, saveMaxVector);

			}

		}

	}

}



// Read arguments ==========================================================
void ProgGpuCorrelation::readParams()
{

    fn_ref = getParam("-i_ref");
    fn_exp = getParam("-i_exp");
    fn_out = getParam("-o");
   	generate_out = checkParam("--classify");
   	fn_classes_out = getParam("--classify");
   	significance = checkParam("--significance");
   	simplifiedMd = checkParam("--simplifiedMd");
   	if(significance){
   		alpha=getDoubleParam("--significance");
   		keepN=false;
   	}
   	if(checkParam("--keep_best") && !significance){
   		keepN=true;
   		n_keep=getIntParam("--keep_best");
   	}
   	if(!keepN && !significance){
   		keepN=true;
   		n_keep=getIntParam("--keep_best");
   	}
   	fnDir = getParam("--odir");
   	maxShift = getIntParam("--maxShift");

   	sizePad = getIntParam("--sizePad");

}

// Show ====================================================================

void ProgGpuCorrelation::show()
{
    std::cout
	<< "Input projected:                " << fn_ref    << std::endl
	<< "Input experimental:             " << fn_exp    << std::endl
	<< "Generate output images (y/n):   " << generate_out    << std::endl
    ;
}

// usage ===================================================================
void ProgGpuCorrelation::defineParams()
{

	addParamsLine("   -i_ref  <md_ref_file>                : Metadata file with input reference images");
	addParamsLine("   -i_exp  <md_exp_file>                : Metadata file with input experimental images");
    addParamsLine("   -o      <md_out>                     : Output metadata file");
	addParamsLine("   [--classify <md_classes_out=\"output_classes.xmd\">]	       : To generate the aligned output images and write the associated metadata");
	addParamsLine("   [--keep_best <N=2>]  			       : To keep N aligned images with the highest correlation");
	addParamsLine("   [--significance <alpha=0.2>]  	   : To use significance with the indicated value");
	addParamsLine("   [--odir <outputDir=\".\">]           : Output directory to save the aligned images");
    addParamsLine("   [--maxShift <s=10>]                  : Maximum shift allowed (+-this amount)");
    addParamsLine("   [--simplifiedMd <b=false>]     : To generate a simplified metadata with only the maximum weight image stores");
    addParamsLine("   [--sizePad <pad=100>]    ");
    addUsageLine("Computes the correlation between a set of experimental images with respect "
    		     "to a set of reference images with CUDA in GPU");

}

int check_gpu_memory(size_t Xdim, size_t Ydim, int percent){
	float data[3]={0, 0, 0};
	cuda_check_gpu_memory(data);
	int bytes = 8*(2*((2*Xdim)-1)*((2*Ydim)-1) + 2*(360*(Xdim/2)));
	return (int)((data[1]*percent/100)/bytes);
}


void calculate_weights(MultidimArray<float> &matrixCorrCpu, MultidimArray<float> &matrixCorrCpu_mirror, MultidimArray<float> &corrTotalRow,
		MultidimArray<float> &weights, int Nref, size_t mdExpSize, size_t mdInSize, MultidimArray<float> &weightsMax, bool simplifiedMd){

	MultidimArray<float> colAux;
	for(int i=0; i<2*mdInSize; i++){
		if(i<mdInSize){
			matrixCorrCpu.getRow(i,colAux); //col
			corrTotalRow.setCol(i, colAux);
		}else{
			matrixCorrCpu_mirror.getRow(i-mdInSize,colAux); //col
			corrTotalRow.setCol(i, colAux);
		}
	}
	MultidimArray<float> corrTotalCol(1,1,2*mdExpSize, mdInSize);
	MultidimArray<float> rowAux;
	for(int i=0; i<2*mdExpSize; i++){
		if(i<mdExpSize){
			matrixCorrCpu.getCol(i,rowAux); //row
			corrTotalCol.setRow(i, rowAux);
		}else{
			matrixCorrCpu_mirror.getCol(i-mdExpSize,rowAux); //row
			corrTotalCol.setRow(i, rowAux);
		}
	}

	//Order the correlation matrix by rows and columns
	MultidimArray<float> rowCorr;
	MultidimArray<int> rowIndexOrder;
	MultidimArray<int> corrOrderByRowIndex(1,1,mdExpSize, 2*mdInSize);

	MultidimArray<float> colCorr;
	MultidimArray<int> colIndexOrder;
	MultidimArray<int> corrOrderByColIndex(1,1,2*mdExpSize, mdInSize);

	for (size_t i=0; i<mdExpSize; i++){
		corrTotalRow.getRow(i, rowCorr);
		rowCorr.indexSort(rowIndexOrder);
		corrOrderByRowIndex.setRow(i, rowIndexOrder);
	}
	for (size_t i=0; i<mdInSize; i++){
		corrTotalCol.getCol(i, colCorr);
		colCorr.indexSort(colIndexOrder);
		corrOrderByColIndex.setCol(i, colIndexOrder);
	}
	corrOrderByRowIndex.selfReverseX();
	corrOrderByColIndex.selfReverseY();


	//AJ To calculate the weights of every image
	MultidimArray<float> weights1(1,1,mdExpSize,2*mdInSize);
	MultidimArray<float> weights2(1,1,mdExpSize,2*mdInSize);

	for(int i=0; i<mdExpSize; i++){
		int idxMax = DIRECT_A2D_ELEM(corrOrderByRowIndex,i,0)-1;
		for(int j=0; j<2*mdInSize; j++){
			int idx = DIRECT_A2D_ELEM(corrOrderByRowIndex,i,j)-1;
			float weight;
			if(DIRECT_A2D_ELEM(corrTotalRow,i,idx)<0)
				weight=0.0;
			else
				weight = 1.0 - (j/(float)corrOrderByRowIndex.xdim);
			weight *= DIRECT_A2D_ELEM(corrTotalRow,i,idx) / DIRECT_A2D_ELEM(corrTotalRow,i,idxMax);
			DIRECT_A2D_ELEM(weights1, i, idx) = weight;
		}
	}
	for(int i=0; i<mdInSize; i++){
		int idxMax = DIRECT_A2D_ELEM(corrOrderByColIndex,0,i)-1;
		for(int j=0; j<2*mdExpSize; j++){
			int idx = DIRECT_A2D_ELEM(corrOrderByColIndex,j,i)-1;
			float weight;
			if(DIRECT_A2D_ELEM(corrTotalCol,idx,i)<0)
				weight=0.0;
			else
				weight = 1.0 - (j/(float)corrOrderByColIndex.ydim);
			weight *= DIRECT_A2D_ELEM(corrTotalCol,idx,i) / DIRECT_A2D_ELEM(corrTotalCol,idxMax,i);
			if(idx<mdExpSize){
				DIRECT_A2D_ELEM(weights2, idx, i) = weight;
			}else{
				DIRECT_A2D_ELEM(weights2, idx-mdExpSize, i+mdInSize) = weight;
			}
		}
	}
	weights=weights1*weights2;


	//AJ
	MultidimArray<float> rowWeights;
	MultidimArray<int> rowIndexOrderWeights;
	MultidimArray<int> weightsOrderByRowIndex(1,1,mdExpSize, 2*mdInSize);
	for (size_t i=0; i<mdExpSize; i++){
		weights.getRow(i, rowWeights);
		rowWeights.indexSort(rowIndexOrderWeights);
		weightsOrderByRowIndex.setRow(i, rowIndexOrderWeights);
	}
	weightsOrderByRowIndex.selfReverseX();
	for(int i=0; i<mdExpSize; i++){
		int idxMax = DIRECT_A2D_ELEM(weightsOrderByRowIndex,i,0)-1;
		for(int j=Nref; j<2*mdInSize; j++){
			int idx = DIRECT_A2D_ELEM(weightsOrderByRowIndex,i,j)-1;
			DIRECT_A2D_ELEM(weights, i, idx) = 0;
		}
	}
	//END AJ


	//AJ new to store the maximum weight for every exp image
	if(simplifiedMd && Nref>1){
		weightsMax.resize(mdExpSize);
		for(int i=0; i<mdInSize; i++){
			for(int j=0; j<mdExpSize; j++){
				if(DIRECT_A2D_ELEM(weights,j,i)!=0){
					if(DIRECT_A2D_ELEM(weights,j,i)>DIRECT_A1D_ELEM(weightsMax,j))
						DIRECT_A1D_ELEM(weightsMax,j) = DIRECT_A2D_ELEM(weights,j,i);
				}
				if(DIRECT_A2D_ELEM(weights,j,i+mdInSize)!=0){
					if(DIRECT_A2D_ELEM(weights,j,i+mdInSize)>DIRECT_A1D_ELEM(weightsMax,j))
						DIRECT_A1D_ELEM(weightsMax,j) = DIRECT_A2D_ELEM(weights,j,i+mdInSize);
				}
			}
		}
	}
	//END AJ

}


void generate_metadata(MetaData SF, MetaData SFexp, FileName fnDir, FileName fn_out, size_t mdExpSize, size_t mdInSize, MultidimArray<float> &weights,
		MultidimArray<float> &corrTotalRow, MultidimArray<float> *matrixTransCpu, MultidimArray<float> *matrixTransCpu_mirror, int maxShift,
		MultidimArray<float> &weightsMax, bool simplifiedMd, int Nref){

	double maxShift2 = maxShift*maxShift;
	Matrix2D<double> bestM(3,3);
	MultidimArray<float> out2(3,3);
	Matrix2D<double>out2Matrix(3,3);
	MDRow rowOut;
	MetaData mdOut;
	String nameImg, nameRef;
	bool flip;
	double rot, tilt, psi;
	int idxJ;
	size_t refNum;

	MDIterator *iterExp = new MDIterator(SFexp);
	MDRow rowExp;
	MDIterator *iter = new MDIterator();
	MDRow row;

	for(int i=0; i<mdExpSize; i++){

		iter->init(SF);
		for(int j=0; j<2*mdInSize; j++){

			if(j%mdInSize==0)
				iter->init(SF);
			SF.getRow(row, iter->objId);

			if(DIRECT_A2D_ELEM(weights,i,j)!=0){

				//AJ new to store the maximum weight for every exp image
				if(simplifiedMd && Nref>1){
					if(DIRECT_A2D_ELEM(weights,i,j)!=DIRECT_A1D_ELEM(weightsMax,i)){
						if(iter->hasNext())
							iter->moveNext();
						continue;
					}
				}
				//END AJ

				size_t itemId;
				SFexp.getRow(rowExp, iterExp->objId);
				//rowExp.getValue(MDL_IMAGE, nameImg);
				//rowExp.getValue(MDL_ITEM_ID, itemId);
				//rowOut
				//rowExp.setValue(MDL_ITEM_ID, itemId);
				//rowExp.setValue(MDL_IMAGE,nameImg);
				rowExp.setValue(MDL_WEIGHT, (double)DIRECT_A2D_ELEM(weights, i, j));
				rowExp.setValue(MDL_MAXCC, (double)DIRECT_A2D_ELEM(corrTotalRow, i, j));
				if(j<mdInSize){
					flip = false;
					matrixTransCpu[j].getSlice(i, out2); //matrixTransCpu[i].getSlice(j, out2);
					idxJ = j;
				}else{
					flip = true;
					matrixTransCpu_mirror[j-mdInSize].getSlice(i, out2); //matrixTransCpu_mirror[i].getSlice(j-mdInSize, out2);
					idxJ = j-mdInSize;
				}

				//AJ NEW
				MAT_ELEM(bestM,0,0) = DIRECT_A2D_ELEM(out2,0,0);
				MAT_ELEM(bestM,0,1)=DIRECT_A2D_ELEM(out2,0,1);
				MAT_ELEM(bestM,0,2)=DIRECT_A2D_ELEM(out2,0,2);

				MAT_ELEM(bestM,1,0)=DIRECT_A2D_ELEM(out2,1,0);
				MAT_ELEM(bestM,1,1)=DIRECT_A2D_ELEM(out2,1,1);
				MAT_ELEM(bestM,1,2)=DIRECT_A2D_ELEM(out2,1,2);

				MAT_ELEM(bestM,2,0)=0.0;
				MAT_ELEM(bestM,2,1)=0.0;
				MAT_ELEM(bestM,2,2)=1.0;
				bestM = bestM.inv();
				//FIN AJ NEW

				double shiftX = MAT_ELEM(bestM,0,2);//(double)DIRECT_A2D_ELEM(out2,0,2);
				double shiftY = MAT_ELEM(bestM,1,2);//(double)DIRECT_A2D_ELEM(out2,1,2);
				if (shiftX*shiftX + shiftY*shiftY > maxShift2){
					if(iter->hasNext())
						iter->moveNext();
					continue;
				}

				//rowOut
				rowExp.setValue(MDL_FLIP, flip);

				double scale;
				/*MAT_ELEM(bestM,0,0)=MAT_ELEM(out2Matrix,0,0);//DIRECT_A2D_ELEM(out2,0,0);
				MAT_ELEM(bestM,0,1)=MAT_ELEM(out2Matrix,0,1);//DIRECT_A2D_ELEM(out2,0,1);
				MAT_ELEM(bestM,0,2)=MAT_ELEM(out2Matrix,0,2);//DIRECT_A2D_ELEM(out2,0,2);
				MAT_ELEM(bestM,1,0)=MAT_ELEM(out2Matrix,1,0);//DIRECT_A2D_ELEM(out2,1,0);
				MAT_ELEM(bestM,1,1)=MAT_ELEM(out2Matrix,1,1);//DIRECT_A2D_ELEM(out2,1,1);
				MAT_ELEM(bestM,1,2)=MAT_ELEM(out2Matrix,1,2);//DIRECT_A2D_ELEM(out2,1,2);
				*/

				MAT_ELEM(bestM,2,0)=0.0;
				MAT_ELEM(bestM,2,1)=0.0;
				MAT_ELEM(bestM,2,2)=1.0;
				if(flip){
					MAT_ELEM(bestM,0,0)*=-1; //bestM
					MAT_ELEM(bestM,1,0)*=-1; //bestM
				}
				bestM=bestM.inv(); //bestM

				transformationMatrix2Parameters2D(bestM,flip,scale,shiftX,shiftY,psi); //bestM
				if (flip)
					shiftX*=-1;

				//AJ NEW
				if(flip){
					shiftX*=-1;
					//shiftY*=-1;
					psi*=-1;
				}
				//FIN AJ NEW

				//rowOut
				rowExp.setValue(MDL_SHIFT_X, -shiftX);
				rowExp.setValue(MDL_SHIFT_Y, -shiftY);
				//rowExp.setValue(MDL_SHIFT_Z, 0.0);
				row.getValue(MDL_ANGLE_ROT, rot);
				rowExp.setValue(MDL_ANGLE_ROT, rot);
				row.getValue(MDL_ANGLE_TILT, tilt);
				rowExp.setValue(MDL_ANGLE_TILT, tilt);
				rowExp.setValue(MDL_ANGLE_PSI, psi);
				//rowOut
				if(row.containsLabel(MDL_ITEM_ID))
					row.getValue(MDL_ITEM_ID, refNum);
				else
					refNum = idxJ+1;
				rowExp.setValue(MDL_REF, (int)refNum);
				mdOut.addRow(rowExp);
			}
			if(iter->hasNext())
				iter->moveNext();
		}
		if(iterExp->hasNext())
			iterExp->moveNext();
	}
	String fnFinal=formatString("%s/%s",fnDir.c_str(),fn_out.c_str());
	mdOut.write(fnFinal);

	delete iterExp;

}


void generate_output_classes(MetaData SF, MetaData SFexp, FileName fnDir, size_t mdExpSize, size_t mdInSize,
		MultidimArray<float> &weights, MultidimArray<float> *matrixTransCpu, MultidimArray<float> *matrixTransCpu_mirror,
		int maxShift, FileName fn_classes_out, MultidimArray<float> &weightsMax, bool simplifiedMd, int Nref){

	double maxShift2 = maxShift*maxShift;
	MultidimArray<float> out2(3,3);
	Matrix2D<double> out2Matrix(3,3);
	double rot, tilt, psi;
	int *NexpVector;

	size_t xAux, yAux, zAux, nAux;
	getImageSize(SF,xAux,yAux,zAux,nAux);
	FileName fnImgNew, fnExpNew, fnRoot, fnStackOut, fnOut, fnStackMD, fnClass;
	Image<double> Inew, Iexp_aux, Inew2, Iexp_out;
	Matrix2D<double> E(3,3);
	MultidimArray<float> auxtr(3,3);
	Matrix2D<double> auxtrMatrix(3,3);
	MultidimArray<double> refSum(1, 1, yAux, xAux);
	bool firstTime=true;
	size_t refNum;
	MultidimArray<double> zeros(1, 1, yAux, xAux);

	// Generate mask
	Mask mask;
    mask.type = BINARY_CIRCULAR_MASK;
	mask.mode = INNER_MASK;
	size_t rad = (size_t)std::min(xAux*0.5, yAux*0.5);
	mask.R1 = rad;
	mask.resize(yAux,xAux);
	mask.get_binary_mask().setXmippOrigin();
	mask.generate_mask();

	CorrelationAux auxCenter;
	RotationalCorrelationAux auxCenter2;

	MDIterator *iterSF = new MDIterator(SF);
	MDRow rowSF;
	MDIterator *iterSFexp = new MDIterator();
	MDRow rowSFexp;

	bool read = false;
	int countingClasses=1;
	bool skip_image;
	NexpVector = new int[mdInSize];
	for(int i=0; i<mdInSize; i++){
		NexpVector[i]=0;
		bool change=false;
		double normWeight=0;

		SF.getRow(rowSF, iterSF->objId);
		if(rowSF.containsLabel(MDL_ITEM_ID))
			rowSF.getValue(MDL_ITEM_ID, refNum);
		else
			refNum=countingClasses;

		iterSFexp->init(SFexp);

		refSum.initZeros();

		fnRoot=fn_classes_out.withoutExtension();
		fnStackOut=formatString("%s/%s.stk",fnDir.c_str(),fnRoot.c_str());
		if(fnStackOut.exists() && firstTime)
			fnStackOut.deleteFile();

		firstTime=false;
		for(int j=0; j<mdExpSize; j++){

			read = false;
			skip_image=false;

			long int pointer1=i*xAux*yAux;
			long int pointer2=i*xAux*yAux;

			if(DIRECT_A2D_ELEM(weights,j,i)!=0){

				//AJ new to store the maximum weight for every exp image
				if(simplifiedMd && Nref>1){
					if(DIRECT_A2D_ELEM(weights,j,i)!=DIRECT_A1D_ELEM(weightsMax,j))
						skip_image=true;
				}
				//END AJ

				if(!skip_image){
					matrixTransCpu[i].getSlice(j, auxtr); //matrixTransCpu[j].getSlice(i, auxtr);
					//AJ NEW
					MAT_ELEM(E,0,0)=DIRECT_A2D_ELEM(auxtr,0,0);
					MAT_ELEM(E,0,1)=DIRECT_A2D_ELEM(auxtr,0,1);
					MAT_ELEM(E,0,2)=DIRECT_A2D_ELEM(auxtr,0,2);

					MAT_ELEM(E,1,0)=DIRECT_A2D_ELEM(auxtr,1,0);
					MAT_ELEM(E,1,1)=DIRECT_A2D_ELEM(auxtr,1,1);
					MAT_ELEM(E,1,2)=DIRECT_A2D_ELEM(auxtr,1,2);

					MAT_ELEM(E,2,0)=0.0;
					MAT_ELEM(E,2,1)=0.0;
					MAT_ELEM(E,2,2)=1.0;
					E = E.inv();
					//FIN AJ NEW

					double shiftX = MAT_ELEM(E,0,2);//(double)DIRECT_A2D_ELEM(auxtr,0,2);
					double shiftY = MAT_ELEM(E,1,2);//(double)DIRECT_A2D_ELEM(auxtr,1,2);
					if (shiftX*shiftX + shiftY*shiftY > maxShift2)
						skip_image=true;
				}

				if(!skip_image){

					if(!read){
						SFexp.getRow(rowSFexp, iterSFexp->objId);
						rowSFexp.getValue(MDL_IMAGE, fnExpNew);
						Iexp_aux.read(fnExpNew);
						read = true;
					}

					NexpVector[i]++;

					/*MAT_ELEM(E,0,0)=MAT_ELEM(auxtrMatrix,0,0);//DIRECT_A2D_ELEM(auxtr,0,0);
					MAT_ELEM(E,0,1)=MAT_ELEM(auxtrMatrix,0,1);//DIRECT_A2D_ELEM(auxtr,0,1);
					MAT_ELEM(E,0,2)=MAT_ELEM(auxtrMatrix,0,2);//DIRECT_A2D_ELEM(auxtr,0,2);
					MAT_ELEM(E,1,0)=MAT_ELEM(auxtrMatrix,1,0);//DIRECT_A2D_ELEM(auxtr,1,0);
					MAT_ELEM(E,1,1)=MAT_ELEM(auxtrMatrix,1,1);//DIRECT_A2D_ELEM(auxtr,1,1);
					MAT_ELEM(E,1,2)=MAT_ELEM(auxtrMatrix,1,2);//DIRECT_A2D_ELEM(auxtr,1,2);
					*/

					MAT_ELEM(E,2,0)=0.0;
					MAT_ELEM(E,2,1)=0.0;
					MAT_ELEM(E,2,2)=1.0;

					selfApplyGeometry(LINEAR,Iexp_aux(),E,IS_NOT_INV,DONT_WRAP,0.0); //E
					//applyGeometry(LINEAR,Iexp_out(),Iexp_aux(),auxtrMatrix,IS_NOT_INV,DONT_WRAP,0.0);

					Iexp_aux().resetOrigin();

					refSum += Iexp_aux()*DIRECT_A2D_ELEM(weights,j,i);
					change=true;
					normWeight+=DIRECT_A2D_ELEM(weights,j,i);
				}
			}
			skip_image=false;
			if(DIRECT_A2D_ELEM(weights,j,i+mdInSize)!=0){

				//AJ new to store the maximum weight for every exp image
				if(simplifiedMd && Nref>1){
					if(DIRECT_A2D_ELEM(weights,j,i+mdInSize)!=DIRECT_A1D_ELEM(weightsMax,j))
						skip_image=true;
				}
				//END AJ

				if(!skip_image){
					matrixTransCpu_mirror[i].getSlice(j, auxtr); //matrixTransCpu_mirror[j].getSlice(i, auxtr);
					//AJ NEW
					MAT_ELEM(E,0,0)=DIRECT_A2D_ELEM(auxtr,0,0);
					MAT_ELEM(E,0,1)=DIRECT_A2D_ELEM(auxtr,0,1);
					MAT_ELEM(E,0,2)=DIRECT_A2D_ELEM(auxtr,0,2);

					MAT_ELEM(E,1,0)=DIRECT_A2D_ELEM(auxtr,1,0);
					MAT_ELEM(E,1,1)=DIRECT_A2D_ELEM(auxtr,1,1);
					MAT_ELEM(E,1,2)=DIRECT_A2D_ELEM(auxtr,1,2);

					MAT_ELEM(E,2,0)=0.0;
					MAT_ELEM(E,2,1)=0.0;
					MAT_ELEM(E,2,2)=1.0;
					E = E.inv();
					//FIN AJ NEW

					double shiftX = MAT_ELEM(E,0,2);//(double)DIRECT_A2D_ELEM(auxtr,0,2);
					double shiftY = MAT_ELEM(E,1,2);//(double)DIRECT_A2D_ELEM(auxtr,1,2);
					if (shiftX*shiftX + shiftY*shiftY > maxShift2)
						skip_image=true;
				}

				if(!skip_image){

					if(!read){
						SFexp.getRow(rowSFexp, iterSFexp->objId);
						rowSFexp.getValue(MDL_IMAGE, fnExpNew);
						Iexp_aux.read(fnExpNew);
						read = true;
					}

					NexpVector[i]++;
					Iexp_aux().selfReverseX();

					/*MAT_ELEM(E,0,0)=MAT_ELEM(auxtrMatrix,0,0);//DIRECT_A2D_ELEM(auxtr,0,0);
					MAT_ELEM(E,0,1)=MAT_ELEM(auxtrMatrix,0,1);//DIRECT_A2D_ELEM(auxtr,0,1);
					MAT_ELEM(E,0,2)=MAT_ELEM(auxtrMatrix,0,2);//DIRECT_A2D_ELEM(auxtr,0,2);
					MAT_ELEM(E,1,0)=MAT_ELEM(auxtrMatrix,1,0);//DIRECT_A2D_ELEM(auxtr,1,0);
					MAT_ELEM(E,1,1)=MAT_ELEM(auxtrMatrix,1,1);//DIRECT_A2D_ELEM(auxtr,1,1);
					MAT_ELEM(E,1,2)=MAT_ELEM(auxtrMatrix,1,2);//DIRECT_A2D_ELEM(auxtr,1,2);
					*/

					MAT_ELEM(E,2,0)=0.0;
					MAT_ELEM(E,2,1)=0.0;
					MAT_ELEM(E,2,2)=1.0;

					//AJ NEW
					MAT_ELEM(E,0,2)*=-1; //E
					MAT_ELEM(E,0,1)*=-1; //E
					MAT_ELEM(E,1,0)*=-1; //E
					//FIN AJ NEW//

					selfApplyGeometry(LINEAR,Iexp_aux(),E,IS_NOT_INV,DONT_WRAP,0.0); //E

					Iexp_aux().resetOrigin();

					refSum += Iexp_aux()*DIRECT_A2D_ELEM(weights,j,i+mdInSize);
					change=true;
					normWeight+=DIRECT_A2D_ELEM(weights,j,i+mdInSize);
				}
			}
			if(iterSFexp->hasNext())
				iterSFexp->moveNext();
		}

		FileName fnStackNo;
		fnStackNo.compose(countingClasses, fnStackOut);
		if(change){
			refSum/=normWeight;
			Inew()=refSum;
			centerImage(Inew(), auxCenter, auxCenter2);
			//masking to avoid wrapping in the edges of the image
			mask.apply_mask(Inew(), Inew2());
			Inew2().resetOrigin();
			Inew2.write(fnStackNo,i,true,WRITE_APPEND);
		}else{
			Inew2() = zeros;
			Inew2.write(fnStackNo,i,true,WRITE_APPEND);
		}

		if(iterSF->hasNext())
			iterSF->moveNext();

		countingClasses++;
	}


	iterSFexp->init(SFexp);
	iterSF->init(SF);

	countingClasses=1;
	Matrix2D<double> bestM(3,3);
	MetaData SFout;
	firstTime=true;
	skip_image=false;
	for(int i=0; i<mdInSize; i++){

		//SF.getRow(rowSF, iterSF->objId);
		//rowSF.getValue(MDL_IMAGE, fnImgNew);
		//fnRoot=fnImgNew.withoutExtension().afterLastOf("/").afterLastOf("@");
		SF.getRow(rowSF, iterSF->objId);
		if(rowSF.containsLabel(MDL_ITEM_ID))
			rowSF.getValue(MDL_ITEM_ID, refNum);
		else
			refNum = countingClasses;

		fnRoot=fn_classes_out.withoutExtension();
		fnStackMD=formatString("%s/%s.xmd", fnDir.c_str(), fnRoot.c_str());
		fnClass.compose(countingClasses, fnStackOut);

		if(fnStackMD.exists() && firstTime)
			fnStackMD.deleteFile();

		firstTime=false;
		size_t id = SFout.addObject();
		SFout.setValue(MDL_REF, (int)refNum, id);
		SFout.setValue(MDL_IMAGE, fnClass, id);
		SFout.setValue(MDL_CLASS_COUNT,(size_t)NexpVector[i], id);

		if(iterSF->hasNext())
			iterSF->moveNext();

		countingClasses++;
	}
	SFout.write("classes@"+fnStackMD, MD_APPEND);

	iterSF->init(SF);
	FileName fnExpIm;
	MDRow row;
	for(int i=0; i<mdInSize; i++){
		skip_image=false;
		SF.getRow(rowSF, iterSF->objId);
		if (rowSF.containsLabel(MDL_ITEM_ID))
			rowSF.getValue(MDL_ITEM_ID, refNum);
		else
			refNum=i+1;

		iterSFexp->init(SFexp);

		MetaData SFq;
		for(int j=0; j<mdExpSize; j++){
			read = false;
			skip_image=false;
			//SFexp.getRow(rowSFexp, iterSFexp->objId);
			//rowSFexp.getValue(MDL_IMAGE, fnExpIm);

			if(DIRECT_A2D_ELEM(weights,j,i)!=0){

				//AJ new to store the maximum weight for every exp image
				if(simplifiedMd && Nref>1){
					if(DIRECT_A2D_ELEM(weights,j,i)!=DIRECT_A1D_ELEM(weightsMax,j))
						skip_image=true;
				}
				//END AJ

				if(!skip_image){
					matrixTransCpu[i].getSlice(j, out2); //matrixTransCpu[j].getSlice(i, out2);
					//AJ NEW
					MAT_ELEM(bestM,0,0)=DIRECT_A2D_ELEM(out2,0,0);
					MAT_ELEM(bestM,0,1)=DIRECT_A2D_ELEM(out2,0,1);
					MAT_ELEM(bestM,0,2)=DIRECT_A2D_ELEM(out2,0,2);

					MAT_ELEM(bestM,1,0)=DIRECT_A2D_ELEM(out2,1,0);
					MAT_ELEM(bestM,1,1)=DIRECT_A2D_ELEM(out2,1,1);
					MAT_ELEM(bestM,1,2)=DIRECT_A2D_ELEM(out2,1,2);

					MAT_ELEM(bestM,2,0)=0.0;
					MAT_ELEM(bestM,2,1)=0.0;
					MAT_ELEM(bestM,2,2)=1.0;
					bestM = bestM.inv();
					//FIN AJ NEW

					double sx = MAT_ELEM(bestM,0,2); //(double)DIRECT_A2D_ELEM(out2,0,2);
					double sy = MAT_ELEM(bestM,1,2); //(double)DIRECT_A2D_ELEM(out2,1,2);
					if (sx*sx + sy*sy > maxShift2)
						skip_image=true;
				}

				if(!skip_image){

					size_t itemId;
					if(!read){
						SFexp.getRow(rowSFexp, iterSFexp->objId);
						//rowSFexp.getValue(MDL_IMAGE, fnExpIm);
						//rowSFexp.getValue(MDL_ITEM_ID, itemId);
						read = true;
					}
					//row
					//row.setValue(MDL_ITEM_ID, itemId);
					//row.setValue(MDL_IMAGE, fnExpIm);
					rowSFexp.setValue(MDL_WEIGHT, (double)DIRECT_A2D_ELEM(weights, j, i));
					rowSFexp.setValue(MDL_FLIP, false);

					double scale, shiftX, shiftY, psi;
					bool flip;
					/*MAT_ELEM(bestM,0,0)=MAT_ELEM(out2Matrix,0,0);//DIRECT_A2D_ELEM(out2,0,0);
					MAT_ELEM(bestM,0,1)=MAT_ELEM(out2Matrix,0,1);//DIRECT_A2D_ELEM(out2,0,1);
					MAT_ELEM(bestM,0,2)=MAT_ELEM(out2Matrix,0,2);//DIRECT_A2D_ELEM(out2,0,2);
					MAT_ELEM(bestM,1,0)=MAT_ELEM(out2Matrix,1,0);//DIRECT_A2D_ELEM(out2,1,0);
					MAT_ELEM(bestM,1,1)=MAT_ELEM(out2Matrix,1,1);//DIRECT_A2D_ELEM(out2,1,1);
					MAT_ELEM(bestM,1,2)=MAT_ELEM(out2Matrix,1,2);//DIRECT_A2D_ELEM(out2,1,2);
					*/

					MAT_ELEM(bestM,2,0)=0.0;
					MAT_ELEM(bestM,2,1)=0.0;
					MAT_ELEM(bestM,2,2)=1.0;
					bestM=bestM.inv(); //bestM

					transformationMatrix2Parameters2D(bestM,flip,scale,shiftX,shiftY,psi); //bestM

					//row
					rowSFexp.setValue(MDL_SHIFT_X, -shiftX);
					rowSFexp.setValue(MDL_SHIFT_Y, -shiftY);
					//rowSFexp.setValue(MDL_SHIFT_Z, 0.0);
					rowSF.getValue(MDL_ANGLE_ROT, rot);
					rowSFexp.setValue(MDL_ANGLE_ROT, rot);
					rowSF.getValue(MDL_ANGLE_TILT, tilt);
					rowSFexp.setValue(MDL_ANGLE_TILT, tilt);
					rowSFexp.setValue(MDL_ANGLE_PSI, psi);
					rowSFexp.setValue(MDL_REF,(int)refNum);
					SFq.addRow(rowSFexp);
				}
			}

			skip_image=false;
			if(DIRECT_A2D_ELEM(weights,j,i+mdInSize)!=0){

				//AJ new to store the maximum weight for every exp image
				if(simplifiedMd && Nref>1){
					if(DIRECT_A2D_ELEM(weights,j,i+mdInSize)!=DIRECT_A1D_ELEM(weightsMax,j))
						skip_image=true;
				}
				//END AJ

				if(!skip_image){
					matrixTransCpu_mirror[i].getSlice(j, out2); //matrixTransCpu_mirror[j].getSlice(i, out2);
					//AJ NEW
					MAT_ELEM(bestM,0,0)=DIRECT_A2D_ELEM(out2,0,0);
					MAT_ELEM(bestM,0,1)=DIRECT_A2D_ELEM(out2,0,1);
					MAT_ELEM(bestM,0,2)=DIRECT_A2D_ELEM(out2,0,2);

					MAT_ELEM(bestM,1,0)=DIRECT_A2D_ELEM(out2,1,0);
					MAT_ELEM(bestM,1,1)=DIRECT_A2D_ELEM(out2,1,1);
					MAT_ELEM(bestM,1,2)=DIRECT_A2D_ELEM(out2,1,2);

					MAT_ELEM(bestM,2,0)=0.0;
					MAT_ELEM(bestM,2,1)=0.0;
					MAT_ELEM(bestM,2,2)=1.0;
					bestM = bestM.inv();
					//FIN AJ NEW

					double sx = MAT_ELEM(bestM,0,2); //(double)DIRECT_A2D_ELEM(out2,0,2);
					double sy = MAT_ELEM(bestM,1,2); //(double)DIRECT_A2D_ELEM(out2,1,2);
					if (sx*sx + sy*sy > maxShift2)
						skip_image=true;
				}

				if(!skip_image){

					size_t itemId;
					if(!read){
						SFexp.getRow(rowSFexp, iterSFexp->objId);
						//rowSFexp.getValue(MDL_IMAGE, fnExpIm);
						//rowSFexp.getValue(MDL_ITEM_ID, itemId);
						read = true;
					}
					//row
					//row.setValue(MDL_ITEM_ID, itemId);
					//row.setValue(MDL_IMAGE, fnExpIm);
					rowSFexp.setValue(MDL_WEIGHT, (double)DIRECT_A2D_ELEM(weights, j, i+mdInSize));
					rowSFexp.setValue(MDL_FLIP, true);

					double scale, shiftX, shiftY, psi;
					bool flip;
					/*MAT_ELEM(bestM,0,0)=MAT_ELEM(out2Matrix,0,0);//DIRECT_A2D_ELEM(out2,0,0);
					MAT_ELEM(bestM,0,1)=MAT_ELEM(out2Matrix,0,1);//DIRECT_A2D_ELEM(out2,0,1);
					MAT_ELEM(bestM,0,2)=MAT_ELEM(out2Matrix,0,2);//DIRECT_A2D_ELEM(out2,0,2);
					MAT_ELEM(bestM,1,0)=MAT_ELEM(out2Matrix,1,0);//DIRECT_A2D_ELEM(out2,1,0);
					MAT_ELEM(bestM,1,1)=MAT_ELEM(out2Matrix,1,1);//DIRECT_A2D_ELEM(out2,1,1);
					MAT_ELEM(bestM,1,2)=MAT_ELEM(out2Matrix,1,2);//DIRECT_A2D_ELEM(out2,1,2);
					*/

					MAT_ELEM(bestM,2,0)=0.0;
					MAT_ELEM(bestM,2,1)=0.0;
					MAT_ELEM(bestM,2,2)=1.0;

					MAT_ELEM(bestM,0,0)*=-1; //bestM
					MAT_ELEM(bestM,1,0)*=-1; //bestM
					bestM=bestM.inv(); //bestM

					transformationMatrix2Parameters2D(bestM,flip,scale,shiftX,shiftY,psi); //bestM

					//AJ NEW
					shiftX*=-1;
					psi*=-1;
					//FIN AJ NEW

					shiftX*=-1;
					//row
					rowSFexp.setValue(MDL_SHIFT_X, -shiftX);
					rowSFexp.setValue(MDL_SHIFT_Y, -shiftY);
					//rowSFexp.setValue(MDL_SHIFT_Z, 0.0);
					rowSF.getValue(MDL_ANGLE_ROT, rot);
					rowSFexp.setValue(MDL_ANGLE_ROT, rot);
					rowSF.getValue(MDL_ANGLE_TILT, tilt);
					rowSFexp.setValue(MDL_ANGLE_TILT, tilt);
					rowSFexp.setValue(MDL_ANGLE_PSI, psi);
					rowSFexp.setValue(MDL_REF,(int)refNum);
					SFq.addRow(rowSFexp);
				}
			}
			if(iterSFexp->hasNext())
				iterSFexp->moveNext();
		}
		MetaData SFq_sorted;
		SFq_sorted.sort(SFq, MDL_IMAGE);
		SFq_sorted.write(formatString("class%06d_images@%s",refNum,fnStackMD.c_str()),MD_APPEND);

		if(iterSF->hasNext())
			iterSF->moveNext();
	}


	delete []NexpVector;
	delete iterSF;
	delete iterSFexp;

}

// Compute correlation --------------------------------------------------------
void ProgGpuCorrelation::run()
{

	//PROJECTION IMAGES
	size_t Xdim, Ydim, Zdim, Ndim;
	SF.read(fn_ref,NULL);
	size_t mdInSize = SF.size();
	getImageSize(SF, Xdim, Ydim, Zdim, Ndim);


	//EXPERIMENTAL IMAGES
	SFexp.read(fn_exp,NULL);
	size_t mdExpSize = SFexp.size();

	// Generate mask
	Mask mask;
    mask.type = BINARY_CIRCULAR_MASK;
	mask.mode = INNER_MASK;
	size_t rad = (size_t)std::min(Xdim*0.48, Ydim*0.48);

	int number = rad;
	int *out = new int[5];

	while(true){
		if (number%2!=0){
			number--;
			continue;
		}
		for (int z=0; z<5; z++)
			out[z]=0;
		primeFactors(number, out);
		if ((out[0]!=0 || out[1]!=0 || out[2]!=0 || out[3]!=0) && out[4]==0){
			rad = number;
			break;
		}
		else
			number--;
	}

	mask.R1 = rad;
	mask.resize(Ydim,Xdim);
	mask.get_binary_mask().setXmippOrigin();
	mask.generate_mask();
	int maskCount = mask.get_binary_mask().sum();

	//AJ check the size of the data to avoid exceed the GPU memory
	float memory[3]={0, 0, 0}; //total, free, used
	cuda_check_gpu_memory(memory);

    int maxGridSize[3];
    cuda_check_gpu_properties(maxGridSize);


	//AJ check_gpu_memory to know how many images we can copy in the gpu memory
    float limit=0.4; //0.877; 1.3;
	int available_images_proj = mdExpSize; //mdInSize
	int available1 = mdExpSize;
	int available2 = mdExpSize;
	if(Xdim*Ydim*mdExpSize*4*100/memory[1]>limit){ //mdInSize
		available1 = floor(memory[1]*(limit/100)/(Xdim*Ydim*4));
	}
	if(Xdim*2*Ydim*2*mdExpSize>maxGridSize[0]){ //mdInSize
		available2 = floor((round(maxGridSize[0]*0.9))/(Xdim*Ydim*2*2));
	}
	if (available1<available2)
		available_images_proj = available1;
	else
		available_images_proj = available2;


	//matrix with all the best transformations in CPU
	MultidimArray<float> *matrixTransCpu = new MultidimArray<float> [mdInSize]; //mdExpSize
	for(int i=0; i<mdInSize; i++) //mdExpSize
		matrixTransCpu[i].coreAllocate(1, mdExpSize, 3, 3); //mdInSize
	MultidimArray<float> *matrixTransCpu_mirror = new MultidimArray<float> [mdInSize]; //mdExpSize
	for(int i=0; i<mdInSize; i++) //mdExpSize
		matrixTransCpu_mirror[i].coreAllocate(1, mdExpSize, 3, 3); //mdInSize

	//correlation matrix
	MultidimArray<float> matrixCorrCpu(1, 1, mdInSize, mdExpSize); //mdExpSize, mdInSize
	MultidimArray<float> matrixCorrCpu_mirror(1, 1, mdInSize, mdExpSize); //mdExpSize, mdInSize

	//Aux vectors with maximum values of correlation in RT and TR steps
	float *max_vector_rt;
	float *max_vector_tr;
	float *max_vector_rt_mirror;
	float *max_vector_tr_mirror;

	//Transformation matrix in GPU and CPU
	TransformMatrix<float> transMat_tr;
	TransformMatrix<float> transMat_rt;
	TransformMatrix<float> transMat_tr_mirror;
	TransformMatrix<float> transMat_rt_mirror;

	TransformMatrix<float> resultTR;
    TransformMatrix<float> resultRT;

	int firstIdx=0;
	bool finish=false;

	mycufftHandle myhandlePadded_tr, myhandleMask_tr, myhandlePolar_tr, myhandleAux_tr, myhandlePaddedB_tr, myhandleMaskB_tr, myhandlePolarB_tr, myhandleAuxB_tr;
	mycufftHandle myhandlePadded_rt, myhandleMask_rt, myhandlePolar_rt, myhandleAux_rt, myhandlePaddedB_rt, myhandleMaskB_rt, myhandlePolarB_rt, myhandleAuxB_rt;
	mycufftHandle ifftcb;

	myStreamHandle myStreamTR, myStreamRT;
	myStreamCreate(myStreamTR);
	myStreamCreate(myStreamRT);


	GpuCorrelationAux d_referenceAux;

	size_t pad_Xdim=2*Xdim-1;
	size_t pad_Ydim=2*Ydim-1;

	number = pad_Xdim;
	while(true){
		if (number%2!=0){
			number++;
			continue;
		}
		for (int z=0; z<5; z++)
			out[z]=0;
		primeFactors(number, out);
		if ((out[0]!=0 || out[1]!=0 || out[2]!=0 || out[3]!=0) && out[4]==0){
			pad_Xdim = number;
			break;
		}
		else
			number++;
	}

	pad_Ydim = pad_Xdim;
	d_referenceAux.XdimOrig=Xdim;
	d_referenceAux.YdimOrig=Ydim;
	d_referenceAux.Xdim=pad_Xdim;
	d_referenceAux.Ydim=pad_Ydim;
	d_referenceAux.XdimPolar=360;
	d_referenceAux.YdimPolar=(size_t)mask.R1;


	StructuresAux myStructureAux_tr, myStructureAux_rt;

	MDIterator *iter = new MDIterator(SFexp); //SF

	GpuMultidimArrayAtCpu<float> original_image_stack;

	//Loop over the reference images
	size_t totalWork=mdInSize*mdExpSize;
	size_t workDone=0;
	init_progress_bar(totalWork);
	size_t lastProgressShown=0;
	while(!finish){

		original_image_stack.resize(Xdim,Ydim,1,available_images_proj);

		//Aux vectors with maximum values of correlation in RT and TR steps
		cpuMalloc((void**)&max_vector_tr, sizeof(float)*available_images_proj);
		cpuMalloc((void**)&max_vector_rt, sizeof(float)*available_images_proj);
		cpuMalloc((void**)&max_vector_tr_mirror, sizeof(float)*available_images_proj);
		cpuMalloc((void**)&max_vector_rt_mirror, sizeof(float)*available_images_proj);


		//Transformation matrix in GPU and CPU
		transMat_tr.resize(myStreamTR, available_images_proj);
		transMat_rt.resize(myStreamRT, available_images_proj);
		transMat_tr_mirror.resize(myStreamTR, available_images_proj);
		transMat_rt_mirror.resize(myStreamRT, available_images_proj);

		resultTR.resize(myStreamTR, available_images_proj);
	    resultRT.resize(myStreamRT, available_images_proj);

		//TODO allocate memory with care
		myStructureAux_tr.padded_image_gpu.resize(pad_Xdim, pad_Ydim, 1, available_images_proj);
		myStructureAux_tr.padded_image2_gpu.resize(pad_Xdim, pad_Ydim, 1, available_images_proj);
		myStructureAux_tr.padded_mask_gpu.resize(pad_Xdim, pad_Ydim, 1, 1);
		myStructureAux_tr.polar_gpu.resize(d_referenceAux.XdimPolar,d_referenceAux.YdimPolar,1,available_images_proj);
		myStructureAux_tr.polar2_gpu.resize(d_referenceAux.XdimPolar,d_referenceAux.YdimPolar,1,available_images_proj);

		myStructureAux_rt.padded_image_gpu.resize(pad_Xdim, pad_Ydim, 1, available_images_proj);
		myStructureAux_rt.padded_image2_gpu.resize(pad_Xdim, pad_Ydim, 1, available_images_proj);
		myStructureAux_rt.padded_mask_gpu.resize(pad_Xdim, pad_Ydim, 1, 1);
		myStructureAux_rt.polar_gpu.resize(d_referenceAux.XdimPolar,d_referenceAux.YdimPolar,1,available_images_proj);
		myStructureAux_rt.polar2_gpu.resize(d_referenceAux.XdimPolar,d_referenceAux.YdimPolar,1,available_images_proj);

		//SF
		preprocess_images_reference(SFexp, firstIdx, available_images_proj, mask, d_referenceAux,
				myhandlePadded_tr, myhandleMask_tr, myhandlePolar_tr, myhandleAux_tr, myStructureAux_tr, iter, myStreamTR);

	    d_referenceAux.maskCount=maskCount;
		d_referenceAux.produceSideInfo(myhandlePaddedB_tr, myhandleMaskB_tr, myStructureAux_tr, myStreamTR);

		//AJ calling a cudaDeviceSyncrhonize to be sure that these images are loaded in gpu memory
		// and available for all the streams
		waitGpu(myStreamTR, true);

		//EXPERIMENTAL IMAGES PART
		size_t expIndex = 0;
		MDRow rowExp;
		FileName fnImgExp;
		MDIterator *iterExp = new MDIterator(SF); //SFexp

		GpuCorrelationAux d_experimentalAuxTR, d_experimentalAuxRT;
		d_experimentalAuxTR.XdimOrig=d_referenceAux.XdimOrig;
		d_experimentalAuxTR.YdimOrig=d_referenceAux.YdimOrig;
		d_experimentalAuxTR.Xdim=d_referenceAux.Xdim;
		d_experimentalAuxTR.Ydim=d_referenceAux.Ydim;
		d_experimentalAuxTR.XdimPolar=d_referenceAux.XdimPolar;
		d_experimentalAuxTR.YdimPolar=d_referenceAux.YdimPolar;

		d_experimentalAuxRT.XdimOrig=d_referenceAux.XdimOrig;
		d_experimentalAuxRT.YdimOrig=d_referenceAux.YdimOrig;
		d_experimentalAuxRT.Xdim=d_referenceAux.Xdim;
		d_experimentalAuxRT.Ydim=d_referenceAux.Ydim;
		d_experimentalAuxRT.XdimPolar=d_referenceAux.XdimPolar;
		d_experimentalAuxRT.YdimPolar=d_referenceAux.YdimPolar;

		//TODO: here we can use threads to carry out the alignment of different images in different threads
		size_t n=0;
		int available_images_exp = mdInSize; //mdExpSize
		while(available_images_exp && iterExp->objId!=0){

			transMat_tr.initialize(myStreamTR);
			transMat_rt.initialize(myStreamRT);
			transMat_tr_mirror.initialize(myStreamTR);
			transMat_rt_mirror.initialize(myStreamRT);

			for(int i=0; i<available_images_proj; i++){
				max_vector_tr[i]=-1;
				max_vector_rt[i]=-1;
				max_vector_tr_mirror[i]=-1;
				max_vector_rt_mirror[i]=-1;
			}

			expIndex = iterExp->objId;
			available_images_exp--;

			SF.getRow(rowExp, expIndex); //SFexp
			rowExp.getValue(MDL_IMAGE, fnImgExp);
			//std::cerr << expIndex << ". Image: " << fnImgExp << std::endl;

			//AJ calling the function to align the images
			bool mirror=false;
			//SFexp
			align_experimental_image(fnImgExp, d_referenceAux, d_experimentalAuxTR, d_experimentalAuxRT, transMat_tr, transMat_rt,
					max_vector_tr, max_vector_rt, SF, available_images_proj, mirror, maxShift,
					myhandlePadded_tr, myhandleMask_tr, myhandlePolar_tr, myhandlePaddedB_tr, myhandleMaskB_tr, myhandlePolarB_tr,
					myhandlePadded_rt, myhandleMask_rt, myhandlePolar_rt, myhandlePaddedB_rt, myhandleMaskB_rt, myhandlePolarB_rt,
					myStructureAux_tr, myStructureAux_rt, myStreamTR, myStreamRT,
					resultTR, resultRT, original_image_stack, ifftcb);


			mirror=true;
			//SFexp
			align_experimental_image(fnImgExp, d_referenceAux, d_experimentalAuxTR, d_experimentalAuxRT, transMat_tr_mirror, transMat_rt_mirror,
							max_vector_tr_mirror, max_vector_rt_mirror, SF, available_images_proj, mirror, maxShift,
							myhandlePadded_tr, myhandleMask_tr, myhandlePolar_tr, myhandlePaddedB_tr, myhandleMaskB_tr, myhandlePolarB_tr,
							myhandlePadded_rt, myhandleMask_rt, myhandlePolar_rt, myhandlePaddedB_rt, myhandleMaskB_rt, myhandlePolarB_rt,
							myStructureAux_tr, myStructureAux_rt, myStreamTR, myStreamRT,
							resultTR, resultRT, original_image_stack, ifftcb);

			//AJ to check the best transformation among all the evaluated
			transMat_tr.copyMatrixToCpu(myStreamTR);
			transMat_tr_mirror.copyMatrixToCpu(myStreamRT);
			transMat_rt.copyMatrixToCpu(myStreamTR);
			transMat_rt_mirror.copyMatrixToCpu(myStreamRT);

			waitGpu(myStreamTR, false);
			waitGpu(myStreamRT, false);

			MultidimArray<float> out2(3,3);
			for(int i=0; i<available_images_proj; i++){
				if(max_vector_tr[i]>max_vector_rt[i]){
					memcpy(MULTIDIM_ARRAY(out2), &transMat_tr.h_data[i*9], 9*sizeof(float));
					matrixTransCpu[n].setSlice(firstIdx+i, out2);
					A2D_ELEM(matrixCorrCpu, n, firstIdx+i) = max_vector_tr[i];
				}else{
					memcpy(MULTIDIM_ARRAY(out2), &transMat_rt.h_data[i*9], 9*sizeof(float));
					matrixTransCpu[n].setSlice(firstIdx+i, out2);
					A2D_ELEM(matrixCorrCpu, n, firstIdx+i) = max_vector_rt[i];
				}
				//mirror image
				if(max_vector_tr_mirror[i]>max_vector_rt_mirror[i]){
					memcpy(MULTIDIM_ARRAY(out2), &transMat_tr_mirror.h_data[i*9], 9*sizeof(float));
					matrixTransCpu_mirror[n].setSlice(firstIdx+i, out2);
					A2D_ELEM(matrixCorrCpu_mirror, n, firstIdx+i) = max_vector_tr_mirror[i];
				}else{
					memcpy(MULTIDIM_ARRAY(out2), &transMat_rt_mirror.h_data[i*9], 9*sizeof(float));
					matrixTransCpu_mirror[n].setSlice(firstIdx+i, out2);
					A2D_ELEM(matrixCorrCpu_mirror, n, firstIdx+i) = max_vector_rt_mirror[i];
				}
			}

			if(iterExp->hasNext())
				iterExp->moveNext();

			n++;
            workDone+=available_images_proj;
            if (size_t(workDone/100)>lastProgressShown)
            {
            	progress_bar(workDone);
            	lastProgressShown=size_t(workDone/100);
            }
		}//end while experimental images

		delete iterExp;
		firstIdx +=available_images_proj;
		int aux;
		aux=available_images_proj;
		if(firstIdx+available_images_proj > mdExpSize){ //mdInSize
			aux=available_images_proj;
			available_images_proj=mdExpSize-firstIdx; //mdInSize
		}
		if(firstIdx==mdExpSize){ //mdInSize
			finish=true;
		}
		if(aux!=available_images_proj){
			myhandlePadded_tr.clear();
			myhandleMask_tr.clear();
			myhandlePolar_tr.clear();
			myhandlePaddedB_tr.clear();
			myhandleMaskB_tr.clear();
			myhandlePolarB_tr.clear();

			myhandlePadded_rt.clear();
			myhandleMask_rt.clear();
			myhandlePolar_rt.clear();
			myhandlePaddedB_rt.clear();
			myhandleMaskB_rt.clear();
			myhandlePolarB_rt.clear();
		}


	}//End loop over the reference images while(!finish)
	progress_bar(totalWork);

	delete iter;

	myhandlePadded_tr.clear();
	myhandleMask_tr.clear();
	myhandlePolar_tr.clear();
	myhandlePaddedB_tr.clear();
	myhandleMaskB_tr.clear();
	myhandlePolarB_tr.clear();

	myhandlePadded_rt.clear();
	myhandleMask_rt.clear();
	myhandlePolar_rt.clear();
	myhandlePaddedB_rt.clear();
	myhandleMaskB_rt.clear();
	myhandlePolarB_rt.clear();

	MultidimArray<float> weights(1,1,mdExpSize,2*mdInSize);
	MultidimArray<float> weightsMax;
	MultidimArray<float> corrTotalRow(1,1,mdExpSize, 2*mdInSize);
	int Nref;
	if(keepN){
		Nref=n_keep;
	}else if(significance){
		Nref=round(corrTotalRow.xdim*alpha);
		if(Nref==0)
			Nref=1;
	}

	calculate_weights(matrixCorrCpu, matrixCorrCpu_mirror, corrTotalRow, weights, Nref, mdExpSize, mdInSize, weightsMax, simplifiedMd);

	std::cerr << "Creating output metadatas..." << std::endl;

	generate_metadata(SF, SFexp, fnDir, fn_out, mdExpSize, mdInSize, weights, corrTotalRow, matrixTransCpu,
			matrixTransCpu_mirror, maxShift, weightsMax, simplifiedMd, Nref);

	if(generate_out)
		generate_output_classes(SF, SFexp, fnDir, mdExpSize, mdInSize, weights, matrixTransCpu,
				matrixTransCpu_mirror, maxShift, fn_classes_out, weightsMax, simplifiedMd, Nref);

	//Free memory in CPU
	for(int i=0; i<mdInSize; i++) //mdExpSize
		matrixTransCpu[i].coreDeallocate();
	delete []matrixTransCpu;
	for(int i=0; i<mdInSize; i++) //mdExpSize
		matrixTransCpu_mirror[i].coreDeallocate();
	delete []matrixTransCpu_mirror;

	cpuFree(max_vector_tr);
	cpuFree(max_vector_rt);
	cpuFree(max_vector_tr_mirror);
	cpuFree(max_vector_rt_mirror);



}
