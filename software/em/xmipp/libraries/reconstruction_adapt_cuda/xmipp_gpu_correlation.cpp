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

#include <algorithm>
#include "xmipp_gpu_utils.h"
#include <reconstruction_cuda/cuda_gpu_correlation.h>

#include <math.h>
#include <time.h>
#include <sys/time.h>

void preprocess_images_reference(MetaData &SF, int firstIdx, int numImages, Mask &mask, GpuCorrelationAux &d_correlationAux,
		mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask, mycufftHandle &myhandlePolar, mycufftHandle &myhandleAux,
		StructuresAux &myStructureAux)
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

	GpuMultidimArrayAtCpu<float> original_image_stack(Xdim,Ydim,1,numImages);

	size_t n=0;
	for(int i=firstIdx; i<firstIdx+numImages; i++){

		SF.getValue(MDL_IMAGE,fnImg,i+1);
		std::cerr << i << ". Image: " << fnImg << std::endl;
		Iref.read(fnImg);
		original_image_stack.fillImage(n,Iref()/8);

		n++;
	}

	//AJ new masking and padding
	//original_image_stack.copyToGpu(d_correlationAux.d_original_image);
	GpuMultidimArrayAtGpu<float> image_stack_gpu(Xdim,Ydim,1,numImages);
	original_image_stack.copyToGpu(image_stack_gpu);
	MultidimArray<int> maskArray = mask.get_binary_mask();
	MultidimArray<float> dMask;
	typeCast(maskArray, dMask);
	d_correlationAux.d_mask.resize(Xdim, Ydim, Zdim, 1);
	d_correlationAux.d_mask.copyToGpu(MULTIDIM_ARRAY(dMask));

	/*GpuMultidimArrayAtGpu<float> padded_image_gpu, padded_image2_gpu, padded_mask_gpu;
	padded_image_gpu.resize(pad_Xdim, pad_Ydim, 1, numImages);
	padded_image2_gpu.resize(pad_Xdim, pad_Ydim, 1, numImages);
	padded_mask_gpu.resize(pad_Xdim, pad_Ydim, 1, 1);*/

	padding_masking(image_stack_gpu, d_correlationAux.d_mask, myStructureAux.padded_image_gpu, myStructureAux.padded_image2_gpu,
			myStructureAux.padded_mask_gpu, false);

	//AJ TIME
	timeval start, end;
	double secs;
    gettimeofday(&start, NULL);

    myStructureAux.padded_image_gpu.fft(d_correlationAux.d_projFFT, myhandlePadded);

    //AJ TIME
    gettimeofday(&end, NULL);
    secs = timeval_diff(&end, &start);
    //printf("REF padded_image_gpu.fft: %.16g miliseconds\n", secs * 1000.0);

    myStructureAux.padded_image2_gpu.fft(d_correlationAux.d_projSquaredFFT, myhandlePadded);
    myStructureAux.padded_mask_gpu.fft(d_correlationAux.d_maskFFT, myhandleMask);

	//Polar transform of the projected images
	/*GpuMultidimArrayAtGpu<float> polar_gpu(angles,radius,1,numImages);
	GpuMultidimArrayAtGpu<float> polar2_gpu(angles,radius,1,numImages);*/
	cuda_cart2polar(image_stack_gpu, myStructureAux.polar_gpu, myStructureAux.polar2_gpu, false); //d_correlationAux.d_original_image
	//FFT

	//AJ TIME
	timeval start2, end2;
	double secs2;
    gettimeofday(&start2, NULL);


    myStructureAux.polar_gpu.fft(d_correlationAux.d_projPolarFFT, myhandlePolar);

    //AJ TIME
    gettimeofday(&end2, NULL);
    secs2 = timeval_diff(&end2, &start2);
    //printf("REF polar_gpu.fft: %.16g miliseconds\n", secs2 * 1000.0);

    myStructureAux.polar2_gpu.fft(d_correlationAux.d_projPolarSquaredFFT, myhandlePolar);

	/*/AJ for debugging
	mycufftHandle myhandle;
	size_t xAux1= padded_image_gpu.Xdim;
	size_t yAux1= padded_image_gpu.Ydim;
	size_t nAux1= padded_image_gpu.Ndim;
	GpuMultidimArrayAtGpu<float> aux(xAux1,yAux1,1,nAux1);
	d_correlationAux.d_projFFT.ifft(aux, myhandle);
	GpuMultidimArrayAtCpu<float> auxCpu1(xAux1,yAux1,1,nAux1);
	auxCpu1.copyFromGpu(aux);
	int pointer1=0;
	for(int i=0; i<nAux1; i++){
	MultidimArray<float> padded1;
	FileName fnImgPad1;
	Image<float> Ipad1;
	padded1.coreAllocate(1, 1, yAux1, xAux1);
	memcpy(MULTIDIM_ARRAY(padded1), &auxCpu1.data[pointer1], xAux1*yAux1*sizeof(float));
	fnImgPad1.compose("original", i+1, "mrc");
	Ipad1()=padded1;
	Ipad1.write(fnImgPad1);
	padded1.coreDeallocate();
	pointer1 += xAux1*yAux1;
	}
	//END AJ/*/

}



void preprocess_images_experimental(MetaData &SF, FileName &fnImg, int numImagesRef, GpuMultidimArrayAtGpu<float> &mask,
		GpuMultidimArrayAtGpu< std::complex<float> > &d_maskFFT, GpuCorrelationAux &d_correlationAux, bool rotation,
		int firstStep, bool mirror, mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask, mycufftHandle &myhandlePolar,
		StructuresAux &myStructureAux)
{
	size_t Xdim, Ydim, Zdim, Ndim;
	getImageSize(SF,Xdim,Ydim,Zdim,Ndim);
	size_t pad_Xdim=d_correlationAux.Xdim;
	size_t pad_Ydim=d_correlationAux.Ydim;
	size_t radius=d_correlationAux.YdimPolar;
	size_t angles = d_correlationAux.XdimPolar;

	//GpuMultidimArrayAtGpu<float> image_stack_gpu(Xdim,Ydim,1,numImagesRef);

	if(firstStep==0){

		Image<float> Iref;

		GpuMultidimArrayAtCpu<float> original_image_stack(Xdim,Ydim,1,numImagesRef);

		Iref.read(fnImg);

		//AJ adding noise to the image
		//Iref().addNoise(0, 2000, "gaussian");
		//END AJ adding noise

		//AJ mirror of the image
		if(mirror)
			Iref().selfReverseX();
		//END AJ mirror

		for(size_t i=0; i<numImagesRef; i++)
			original_image_stack.fillImage(i,Iref()/8);
		original_image_stack.copyToGpu(d_correlationAux.d_original_image);
		//original_image_stack.copyToGpu(image_stack_gpu);

	}//else{
		//d_correlationAux.d_original_image.copyGpuToGpu(image_stack_gpu);
	//}

	//AJ new masking and padding
	/*MultidimArray<int> maskArray = mask.get_binary_mask();
	MultidimArray<float> dMask;
	typeCast(maskArray, dMask);
	GpuMultidimArrayAtGpu<float> mask_device(Xdim, Ydim, Zdim, 1);
	mask_device.copyToGpu(MULTIDIM_ARRAY(dMask));*/

	//GpuMultidimArrayAtGpu<float> padded_image_gpu, padded_image2_gpu, padded_mask_gpu;
	if(!rotation){
		//padded_image_gpu.resize(pad_Xdim, pad_Ydim, 1, numImagesRef);
		//padded_image2_gpu.resize(pad_Xdim, pad_Ydim, 1, numImagesRef);
		//padded_mask_gpu.resize(pad_Xdim, pad_Ydim, 1, 1);

		padding_masking(d_correlationAux.d_original_image, mask, myStructureAux.padded_image_gpu, myStructureAux.padded_image2_gpu,
				myStructureAux.padded_mask_gpu, true);

		//FFT
		//AJ TIME
		timeval start, end;
		double secs;
	    gettimeofday(&start, NULL);

	    myStructureAux.padded_image_gpu.fft(d_correlationAux.d_projFFT, myhandlePadded);

	    //AJ TIME
	    gettimeofday(&end, NULL);
	    secs = timeval_diff(&end, &start);
	    //printf("EXP padded_image_gpu.fft: %.16g miliseconds\n", secs * 1000.0);


	    myStructureAux.padded_image2_gpu.fft(d_correlationAux.d_projSquaredFFT, myhandlePadded);
		//padded_mask_gpu.fft(d_correlationAux.d_maskFFT, myhandleMask);
		d_maskFFT.copyGpuToGpu(d_correlationAux.d_maskFFT);

	}

	//GpuMultidimArrayAtGpu<float> polar_gpu, polar2_gpu;
	if(rotation){
		//polar_gpu.resize(angles,radius,1,numImagesRef);
		//polar2_gpu.resize(angles,radius,1,numImagesRef);
		cuda_cart2polar(d_correlationAux.d_original_image, myStructureAux.polar_gpu, myStructureAux.polar2_gpu, true);
		//FFT
		//AJ TIME
		timeval start2, end2;
		double secs2;
	    gettimeofday(&start2, NULL);


	    myStructureAux.polar_gpu.fft(d_correlationAux.d_projPolarFFT, myhandlePolar);

	    //AJ TIME
	    gettimeofday(&end2, NULL);
	    secs2 = timeval_diff(&end2, &start2);
	    //printf("EXP polar_gpu.fft: %.16g miliseconds\n", secs2 * 1000.0);

	    myStructureAux.polar2_gpu.fft(d_correlationAux.d_projPolarSquaredFFT, myhandlePolar);
	}

	/*/AJ for debugging
	if(!rotation){
	mycufftHandle myhandleAux;
	size_t xAux1= padded_image_gpu.Xdim;
	size_t yAux1= padded_image_gpu.Ydim;
	size_t nAux1= padded_image_gpu.Ndim;
	GpuMultidimArrayAtGpu<float> aux(xAux1,yAux1,1,nAux1);
	d_correlationAux.d_projFFT.ifft(aux, myhandleAux);
	GpuMultidimArrayAtCpu<float> auxCpu1(xAux1,yAux1,1,nAux1);
	auxCpu1.copyFromGpu(aux);
	int pointer1=0;
	for(int i=0; i<nAux1; i++){
	MultidimArray<float> padded1;
	FileName fnImgPad1;
	Image<float> Ipad1;
	padded1.coreAllocate(1, 1, yAux1, xAux1);
	memcpy(MULTIDIM_ARRAY(padded1), &auxCpu1.data[pointer1], xAux1*yAux1*sizeof(float));
	fnImgPad1.compose("condemor", i+1, "mrc");
	Ipad1()=padded1;
	Ipad1.write(fnImgPad1);
	padded1.coreDeallocate();
	pointer1 += xAux1*yAux1;
	}
	}
	//END AJ/*/

}


void preprocess_images_experimental_transform(GpuCorrelationAux &d_correlationAux, GpuMultidimArrayAtGpu<float> &mask,
		GpuMultidimArrayAtGpu< std::complex<float> > &d_maskFFT, bool rotation, int step, mycufftHandle &myhandlePadded,
		mycufftHandle &myhandleMask, mycufftHandle &myhandlePolar, StructuresAux &myStructureAux)
{

	size_t Xdim = d_correlationAux.d_transform_image.Xdim;
	size_t Ydim = d_correlationAux.d_transform_image.Ydim;
	size_t Zdim = d_correlationAux.d_transform_image.Zdim;
	size_t Ndim = d_correlationAux.d_transform_image.Ndim;
	size_t pad_Xdim=d_correlationAux.Xdim;
	size_t pad_Ydim=d_correlationAux.Ydim;
	size_t radius=d_correlationAux.YdimPolar;
	size_t angles = d_correlationAux.XdimPolar;

	/*MultidimArray<int> maskArray = mask.get_binary_mask();
	MultidimArray<float> dMask;
	typeCast(maskArray, dMask);
	GpuMultidimArrayAtGpu<float> mask_device(Xdim, Ydim, Zdim, 1);
	mask_device.copyToGpu(MULTIDIM_ARRAY(dMask));*/

	//GpuMultidimArrayAtGpu<float> padded_image_gpu, padded_image2_gpu, padded_mask_gpu;
	if(!rotation){
		//AJ TIME
		timeval start0, end0;
		double secs0;
	    gettimeofday(&start0, NULL);

		//padded_image_gpu.resize(pad_Xdim, pad_Ydim, 1, Ndim);
		//padded_image2_gpu.resize(pad_Xdim, pad_Ydim, 1, Ndim);
		//padded_mask_gpu.resize(pad_Xdim, pad_Ydim, 1, 1);

	    //AJ TIME
	    gettimeofday(&end0, NULL);
	    secs0 = timeval_diff(&end0, &start0);
	    //printf("TRANS resize: %.16g miliseconds\n", secs0 * 1000.0);

		//AJ TIME
		timeval start, end;
		double secs;
	    gettimeofday(&start, NULL);

		padding_masking(d_correlationAux.d_transform_image, mask, myStructureAux.padded_image_gpu, myStructureAux.padded_image2_gpu,
				myStructureAux.padded_mask_gpu, true);

	    //AJ TIME
	    gettimeofday(&end, NULL);
	    secs = timeval_diff(&end, &start);
	    //printf("TRANS padding_masking: %.16g miliseconds\n", secs * 1000.0);

		//AJ TIME
		timeval start1, end1;
		double secs1;
	    gettimeofday(&start1, NULL);

		//FFT
	    myStructureAux.padded_image_gpu.fft(d_correlationAux.d_projFFT, myhandlePadded);
	    myStructureAux.padded_image2_gpu.fft(d_correlationAux.d_projSquaredFFT, myhandlePadded);
		//padded_mask_gpu.fft(d_correlationAux.d_maskFFT, myhandleMask);
		d_maskFFT.copyGpuToGpu(d_correlationAux.d_maskFFT);

	    //AJ TIME
	    gettimeofday(&end1, NULL);
	    secs1 = timeval_diff(&end1, &start1);
	    //printf("TRANS FFTs: %.16g miliseconds\n", secs1 * 1000.0);

	}

	//Polar transform of the projected images
	//GpuMultidimArrayAtGpu<float> polar_gpu, polar2_gpu;
	if(rotation){
		//polar_gpu.resize(angles,radius,1,Ndim);
		//polar2_gpu.resize(angles,radius,1,Ndim);
		cuda_cart2polar(d_correlationAux.d_transform_image, myStructureAux.polar_gpu, myStructureAux.polar2_gpu, true);
		//FFT
		myStructureAux.polar_gpu.fft(d_correlationAux.d_projPolarFFT, myhandlePolar);
		myStructureAux.polar2_gpu.fft(d_correlationAux.d_projPolarSquaredFFT, myhandlePolar);
	}

	/*/AJ for debugging
	if(rotation){
	mycufftHandle myhandleAux;
	size_t xAux1= d_correlationAux.d_transform_image.Xdim;
	size_t yAux1= d_correlationAux.d_transform_image.Ydim;
	size_t nAux1= d_correlationAux.d_transform_image.Ndim;
	//GpuMultidimArrayAtGpu<float> aux(xAux1,yAux1,1,nAux1);
	//d_correlationAux.d_projFFT.ifft(aux, myhandleAux);
	GpuMultidimArrayAtCpu<float> auxCpu1(xAux1,yAux1,1,nAux1);
	auxCpu1.copyFromGpu(d_correlationAux.d_transform_image);
	int pointer1=0;
	for(int i=0; i<nAux1; i++){
	MultidimArray<float> padded1;
	FileName fnImgPad1;
	Image<float> Ipad1;
	padded1.coreAllocate(1, 1, yAux1, xAux1);
	memcpy(MULTIDIM_ARRAY(padded1), &auxCpu1.data[pointer1], xAux1*yAux1*sizeof(float));
	fnImgPad1.compose("condemor", i+1, "mrc");
	Ipad1()=padded1;
	Ipad1.write(fnImgPad1);
	padded1.coreDeallocate();
	pointer1 += xAux1*yAux1;
	}
	}
	//END AJ/*/

}

void align_experimental_image(FileName &fnImgExp, GpuCorrelationAux &d_referenceAux, GpuCorrelationAux &d_experimentalAux,
		TransformMatrix<float> &transMat_tr, TransformMatrix<float> &transMat_rt, float *max_vector_tr, float *max_vector_rt,
		MetaData &SFexp, int available_images_proj, bool mirror, int maxShift,
		mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask, mycufftHandle &myhandlePolar,
		mycufftHandle &myhandlePaddedB, mycufftHandle &myhandleMaskB, mycufftHandle &myhandlePolarB, StructuresAux &myStructureAux)
{

	bool rotation;
	TransformMatrix<float> *transMat;
	float *max_vector;

	for(int firstStep=0; firstStep<2; firstStep++){ //2

		if (firstStep==0){
			rotation = false;
			max_vector = max_vector_tr;
		}else{
			rotation = true;
			max_vector = max_vector_rt;
		}

		//AJ TIME
		timeval start2, end2;
		double secs2;
	    gettimeofday(&start2, NULL);

		preprocess_images_experimental(SFexp, fnImgExp, available_images_proj, d_referenceAux.d_mask, d_referenceAux.d_maskFFT,
				d_experimentalAux, rotation, firstStep, mirror, myhandlePadded, myhandleMask, myhandlePolar, myStructureAux);

	    //AJ TIME
	    gettimeofday(&end2, NULL);
	    secs2 = timeval_diff(&end2, &start2);
	    //printf("preprocess_images_experimental: %.16g miliseconds\n", secs2 * 1000.0);

		//AJ TIME
		timeval start3, end3;
		double secs3;
	    gettimeofday(&start3, NULL);

		if(!rotation){
			d_experimentalAux.maskCount=d_referenceAux.maskCount;
			d_experimentalAux.produceSideInfo(myhandlePaddedB, myhandleMaskB, myStructureAux, d_referenceAux.maskAutocorrelation);
		}

	    //AJ TIME
	    gettimeofday(&end3, NULL);
	    secs3 = timeval_diff(&end3, &start3);
	    //printf("d_experimentalAux.produceSideInfo: %.16g miliseconds\n", secs3 * 1000.0);

		if(firstStep==0)
			transMat = &transMat_tr;
		else
			transMat = &transMat_rt;



		int max_step=6;
		char stepchar[20]="";
		for(int step=0; step<max_step; step++){ //loop over consecutive translations and rotations (TRTRTR or RTRTRT) 6

			/*if(!rotation){
				stepchar[step]='T';
				printf("step %i of %i %s\n",step+1, max_step, stepchar);
			}else{
				stepchar[step]='R';
				printf("step %i of %i %s\n",step+1, max_step, stepchar);
			}*/


			//AJ TIME
			timeval start4, end4;
			double secs4;
		    gettimeofday(&start4, NULL);

			//CORRELATION PART
			//TRANSFORMATION MATRIX CALCULATION
			//printf("Calculating correlation...\n");
			if(!rotation)
				cuda_calculate_correlation(d_referenceAux, d_experimentalAux, *transMat, max_vector, maxShift, myhandlePaddedB, mirror, myStructureAux);
			else
				cuda_calculate_correlation_rotation(d_referenceAux, d_experimentalAux, *transMat, max_vector, maxShift, myhandlePolarB, mirror, myStructureAux);

		    //AJ TIME
		    gettimeofday(&end4, NULL);
		    secs4 = timeval_diff(&end4, &start4);
		    //printf("cuda_calculate_correlation: %.16g miliseconds\n", secs4 * 1000.0);

			//AJ TIME
			timeval start5, end5;
			double secs5;
		    gettimeofday(&start5, NULL);

			//APPLY TRANSFORMATION
			if(step<max_step-1){
				//printf("Applying transformation...\n");
				d_experimentalAux.d_transform_image.resize(d_experimentalAux.d_original_image);
				apply_transform(d_experimentalAux.d_original_image, d_experimentalAux.d_transform_image, *transMat);
			}

		    //AJ TIME
		    gettimeofday(&end5, NULL);
		    secs5 = timeval_diff(&end5, &start5);
		    //printf("apply_transform: %.16g miliseconds\n", secs5 * 1000.0);

			/*/AJ for debugging
			if(!rotation){
			mycufftHandle myhandleAux;
			size_t xAux1= d_experimentalAux.d_transform_image.Xdim;
			size_t yAux1= d_experimentalAux.d_transform_image.Ydim;
			size_t nAux1= d_experimentalAux.d_transform_image.Ndim;
			//GpuMultidimArrayAtGpu<float> aux(xAux1,yAux1,1,nAux1);
			//d_correlationAux.d_projFFT.ifft(aux, myhandleAux);
			GpuMultidimArrayAtCpu<float> auxCpu1(xAux1,yAux1,1,nAux1);
			auxCpu1.copyFromGpu(d_experimentalAux.d_transform_image);
			int pointer1=0;
			for(int i=0; i<nAux1; i++){
			MultidimArray<float> padded1;
			FileName fnImgPad1;
			Image<float> Ipad1;
			padded1.coreAllocate(1, 1, yAux1, xAux1);
			memcpy(MULTIDIM_ARRAY(padded1), &auxCpu1.data[pointer1], xAux1*yAux1*sizeof(float));
			fnImgPad1.compose("transform", i+1, "mrc");
			Ipad1()=padded1;
			Ipad1.write(fnImgPad1);
			padded1.coreDeallocate();
			pointer1 += xAux1*yAux1;
			}
			}
			//END AJ/*/

			/*/AJ TIME
			timeval start, end;
			double secs;
		    gettimeofday(&start, NULL);*/

			//AJ TIME
			timeval start6, end6;
			double secs6;
		    gettimeofday(&start6, NULL);

			//PREPROCESS TO PREPARE DATA TO THE NEXT STEP
			if(step<max_step-1){
				//printf("Re-processing experimental images...\n");
				rotation = !rotation;
				preprocess_images_experimental_transform(d_experimentalAux, d_referenceAux.d_mask, d_referenceAux.d_maskFFT, rotation, step,
						myhandlePadded, myhandleMask, myhandlePolar, myStructureAux);

				//AJ TIME
				timeval start6, end6;
				double secs6;
				  gettimeofday(&start6, NULL);

				if(!rotation){
					d_experimentalAux.maskCount=d_referenceAux.maskCount;
					d_experimentalAux.produceSideInfo(myhandlePaddedB, myhandleMaskB, myStructureAux, d_referenceAux.maskAutocorrelation);
				}

			    //AJ TIME
			    gettimeofday(&end6, NULL);
			    secs6 = timeval_diff(&end6, &start6);
			    //printf("TRANS produceSideInfo: %.16g miliseconds\n", secs6 * 1000.0);

			}

		    //AJ TIME
		    gettimeofday(&end6, NULL);
		    secs6 = timeval_diff(&end6, &start6);
		    //printf("preprocess_images_experimental_transform: %.16g miliseconds\n", secs6 * 1000.0);

		    /*/AJ TIME
		    gettimeofday(&end, NULL);
		    secs = timeval_diff(&end, &start);
		    //printf("preprocess_images_experimental_transform: %.16g miliseconds\n", secs * 1000.0);*/

		}//end for(int step=0; step<6; step++)

	}//end for(int firstStep=0; firstStep<2; firstStep++)

}



// Read arguments ==========================================================
void ProgGpuCorrelation::readParams()
{

    fn_ref = getParam("-i_ref");
    fn_exp = getParam("-i_exp");
    fn_out = getParam("-o");
   	generate_out = checkParam("--out");
   	significance = checkParam("--significance");
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
	addParamsLine("   [--out]  						       : To generate the aligned output images");
	addParamsLine("   [--keep_best <N=2>]  			       : To keep N aligned images with the highest correlation");
	addParamsLine("   [--significance <alpha=0.2>]  	   : To use significance with the indicated value");
	addParamsLine("   [--odir <outputDir=\".\">]           : Output directory to save the aligned images");
    addParamsLine("   [--maxShift <s=20>]                  : Maximum shift allowed (+-this amount)");
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
		MultidimArray<float> &weights, int Nref, size_t mdExpSize, size_t mdInSize){

	//Concatenate corr matrix of the original and mirror images in both directions
	//v2 = v1.sort();
	//v2 = v1.indexSort(); //indexes start at 1
	MultidimArray<float> colAux;
	for(int i=0; i<2*mdInSize; i++){
		if(i<mdInSize){
			matrixCorrCpu.getCol(i,colAux);
			corrTotalRow.setCol(i, colAux);
		}else{
			matrixCorrCpu_mirror.getCol(i-mdInSize,colAux);
			corrTotalRow.setCol(i, colAux);
		}
	}
	MultidimArray<float> corrTotalCol(1,1,2*mdExpSize, mdInSize);
	MultidimArray<float> rowAux;
	for(int i=0; i<2*mdExpSize; i++){
		if(i<mdExpSize){
			matrixCorrCpu.getRow(i,rowAux);
			corrTotalCol.setRow(i, rowAux);
		}else{
			matrixCorrCpu_mirror.getRow(i-mdExpSize,rowAux);
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
		for(int j=0; j<Nref; j++){
			int idx = DIRECT_A2D_ELEM(corrOrderByRowIndex,i,j)-1;
			if(DIRECT_A2D_ELEM(corrTotalRow,i,idx)<0)
				break;
			float weight = 1.0 - (j/(float)corrOrderByRowIndex.xdim);
			DIRECT_A2D_ELEM(weights1, i, idx) = weight;
		}
	}
	for(int i=0; i<mdInSize; i++){
		for(int j=0; j<2*mdExpSize; j++){
			int idx = DIRECT_A2D_ELEM(corrOrderByColIndex,j,i)-1;
			float weight = 1.0 - (j/(float)corrOrderByColIndex.ydim);
			if(idx<mdExpSize){
				DIRECT_A2D_ELEM(weights2, idx, i) = weight;
			}else{
				DIRECT_A2D_ELEM(weights2, idx-mdExpSize, i+mdInSize) = weight;
			}
		}
	}
	weights=weights1*weights2;


}


void generate_metadata(MetaData SF, MetaData SFexp, FileName fnDir, FileName fn_out, size_t mdExpSize, size_t mdInSize, MultidimArray<float> weights,
		MultidimArray<float> corrTotalRow, MultidimArray<float> *matrixTransCpu, MultidimArray<float> *matrixTransCpu_mirror){

	//AJ To generate the output metadata
	MultidimArray<float> out2(3,3);
	MDRow rowOut;
	MetaData mdOut;
	String nameImg, nameRef;
	bool flip;
	double rot, tilt, psi;
	size_t count=0;
	int idxJ;
	for(int i=0; i<mdExpSize; i++){
		for(int j=0; j<2*mdInSize; j++){
			if(DIRECT_A2D_ELEM(weights,i,j)!=0){
				rowOut.setValue(MDL_ITEM_ID, count);
				SFexp.getValue(MDL_IMAGE,nameImg,i+1);
				//printf("nameExp %s \n", nameImg.c_str());
				rowOut.setValue(MDL_IMAGE,nameImg);
				rowOut.setValue(MDL_WEIGHT, (double)DIRECT_A2D_ELEM(weights, i, j));
				rowOut.setValue(MDL_MAXCC, (double)DIRECT_A2D_ELEM(corrTotalRow, i, j));
				if(j<mdInSize){
					flip = false;
					matrixTransCpu[i].getSlice(j, out2);
					idxJ = j;
				}else{
					flip = true;
					matrixTransCpu_mirror[i].getSlice(j-mdInSize, out2);
					idxJ = j-mdInSize;
				}
				rowOut.setValue(MDL_SHIFT_X, (double)DIRECT_A2D_ELEM(out2,0,2));
				rowOut.setValue(MDL_SHIFT_Y, (double)DIRECT_A2D_ELEM(out2,1,2));
				psi = RAD2DEG(atan2(DIRECT_A2D_ELEM(out2,1,0), DIRECT_A2D_ELEM(out2,0,0)));
				SF.getValue(MDL_IMAGE,nameRef,idxJ+1);
				//printf("nameRef %s \n", nameRef.c_str());
				SF.getValue(MDL_ANGLE_ROT,rot,idxJ+1);
				rowOut.setValue(MDL_ANGLE_ROT, rot);
				SF.getValue(MDL_ANGLE_TILT,tilt,idxJ+1);
				rowOut.setValue(MDL_ANGLE_TILT, tilt);
				rowOut.setValue(MDL_ANGLE_PSI, psi);
				rowOut.setValue(MDL_FLIP, flip);
				//printf("rot %lf tilt %lf psi %lf shiftX %lf shiftY %lf flip %i corr %lf weight %lf \n", rot, tilt, psi, (double)DIRECT_A2D_ELEM(out2,0,2), (double)DIRECT_A2D_ELEM(out2,1,2), (int)flip, (double)(DIRECT_A2D_ELEM(corrTotalRow, i, j)), (double)(DIRECT_A2D_ELEM(weights, i, j)));
				mdOut.addRow(rowOut);
				count++;
			}
		}
	}
	String fnFinal=formatString("%s/%s",fnDir.c_str(),fn_out.c_str());
	mdOut.write(fnFinal);


}




// Compute correlation --------------------------------------------------------
void ProgGpuCorrelation::run()
{


	//AJ TIME
	timeval start, end;
	double secs;
    gettimeofday(&start, NULL);


	//PROJECTION IMAGES PART
	//printf("Reference images: \n");
	//Read input metadataFile for projection images
	size_t Xdim, Ydim, Zdim, Ndim;
	SF.read(fn_ref,NULL);
	size_t mdInSize = SF.size();
	getImageSize(SF, Xdim, Ydim, Zdim, Ndim);

	// Generate mask
	Mask mask, maskPolar;
    mask.type = BINARY_CIRCULAR_MASK;
	mask.mode = INNER_MASK;
	size_t rad = (size_t)std::min(Xdim*0.45, Ydim*0.45);
	rad = (rad%2==0) ? rad : (rad+1);
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
    float limit=1.3; //0.877;
	int available_images_proj = mdInSize;
	if(Xdim*Ydim*mdInSize*4*100/memory[1]>limit){ //TODO revisar en otras GPUs
		available_images_proj = floor(memory[1]*(limit/100)/(Xdim*Ydim*4));
	}
	if(Xdim*2*Ydim*2*mdInSize>maxGridSize[0]){
		available_images_proj = floor((round(maxGridSize[0]*0.9))/(Xdim*Ydim));
	}

	SFexp.read(fn_exp,NULL);
	size_t mdExpSize = SFexp.size();

	//matrix with all the best transformations in CPU
	MultidimArray<float> *matrixTransCpu = new MultidimArray<float> [mdExpSize];
	for(int i=0; i<mdExpSize; i++)
		matrixTransCpu[i].coreAllocate(1, mdInSize, 3, 3);
	MultidimArray<float> *matrixTransCpu_mirror = new MultidimArray<float> [mdExpSize];
	for(int i=0; i<mdExpSize; i++)
		matrixTransCpu_mirror[i].coreAllocate(1, mdInSize, 3, 3);

	//correlation matrix
	MultidimArray<float> matrixCorrCpu(1, 1, mdExpSize, mdInSize);
	MultidimArray<float> matrixCorrCpu_mirror(1, 1, mdExpSize, mdInSize);

	//Aux vectors with maximum values of correlation in RT and TR steps
	float *max_vector_rt;
	float *max_vector_tr;
	float *max_vector_rt_mirror;
	float *max_vector_tr_mirror;

	//Transformation matrix in GPU
	TransformMatrix<float> transMat_tr;
	TransformMatrix<float> transMat_rt;
	TransformMatrix<float> transMat_tr_mirror;
	TransformMatrix<float> transMat_rt_mirror;


	int firstIdx=0;
	bool finish=false;

	mycufftHandle myhandlePadded, myhandleMask, myhandlePolar, myhandleAux;
	mycufftHandle myhandlePaddedB, myhandleMaskB, myhandlePolarB, myhandleAuxB;

	GpuCorrelationAux d_referenceAux;

	size_t pad_Xdim=2*Xdim-1;
	pad_Xdim = (pad_Xdim%2==0) ? pad_Xdim:(pad_Xdim+1);
	size_t pad_Ydim=2*Ydim-1;
	pad_Ydim = (pad_Ydim%2==0) ? pad_Ydim:(pad_Ydim+1);
	//printf("pad_Xdim %i, pad_Ydim %i, radius %i \n", pad_Xdim, pad_Ydim, (int)mask.R1);
	d_referenceAux.XdimOrig=Xdim;
	d_referenceAux.YdimOrig=Ydim;
	d_referenceAux.Xdim=pad_Xdim;
	d_referenceAux.Ydim=pad_Ydim;
	d_referenceAux.XdimPolar=360;
	d_referenceAux.YdimPolar=(size_t)mask.R1;

	StructuresAux myStructureAux;

	//Loop over the reference images
	while(!finish){

		//Aux vectors with maximum values of correlation in RT and TR steps
		max_vector_rt = new float [available_images_proj];
		max_vector_tr = new float [available_images_proj];
		max_vector_rt_mirror = new float [available_images_proj];
		max_vector_tr_mirror = new float [available_images_proj];

		//Transformation matrix in GPU
		transMat_tr.resize(available_images_proj);
		transMat_rt.resize(available_images_proj);
		transMat_tr_mirror.resize(available_images_proj);
		transMat_rt_mirror.resize(available_images_proj);

		myStructureAux.padded_image_gpu.resize(pad_Xdim, pad_Ydim, 1, available_images_proj);
		myStructureAux.padded_image2_gpu.resize(pad_Xdim, pad_Ydim, 1, available_images_proj);
		myStructureAux.padded_mask_gpu.resize(pad_Xdim, pad_Ydim, 1, 1);
		myStructureAux.polar_gpu.resize(d_referenceAux.XdimPolar,d_referenceAux.YdimPolar,1,available_images_proj);
		myStructureAux.polar2_gpu.resize(d_referenceAux.XdimPolar,d_referenceAux.YdimPolar,1,available_images_proj);

		//AJ TIME
		timeval start1, end1;
		double secs1;
	    gettimeofday(&start1, NULL);

		preprocess_images_reference(SF, firstIdx, available_images_proj, mask, d_referenceAux,
				myhandlePadded, myhandleMask, myhandlePolar, myhandleAux, myStructureAux);

	    //AJ TIME
	    gettimeofday(&end1, NULL);
	    secs1 = timeval_diff(&end1, &start1);
	    //printf("preprocess_images_reference: %.16g miliseconds\n", secs1 * 1000.0);

		//AJ TIME
		timeval start2, end2;
		double secs2;
	    gettimeofday(&start2, NULL);

	    d_referenceAux.maskCount=maskCount;
		d_referenceAux.produceSideInfo(myhandlePaddedB, myhandleMaskB, myStructureAux);

		//AJ TIME
	    gettimeofday(&end2, NULL);
	    secs2 = timeval_diff(&end2, &start2);
	    //printf("d_referenceAux.produceSideInfo: %.16g miliseconds\n", secs2 * 1000.0);

		//EXPERIMENTAL IMAGES PART
		//printf("Experimental image: \n");

		size_t expIndex = 0;
		MDRow rowExp;
		FileName fnImgExp;
		MDIterator *iterExp = new MDIterator(SFexp);

		GpuCorrelationAux d_experimentalAux;
		d_experimentalAux.XdimOrig=d_referenceAux.XdimOrig;
		d_experimentalAux.YdimOrig=d_referenceAux.YdimOrig;
		d_experimentalAux.Xdim=d_referenceAux.Xdim;
		d_experimentalAux.Ydim=d_referenceAux.Ydim;
		d_experimentalAux.XdimPolar=d_referenceAux.XdimPolar;
		d_experimentalAux.YdimPolar=d_referenceAux.YdimPolar;

		size_t n=0;
		int available_images_exp = mdExpSize;
		while(available_images_exp && iterExp->objId!=0){

			transMat_tr.initialize();
			transMat_rt.initialize();
			transMat_tr_mirror.initialize();
			transMat_rt_mirror.initialize();

			for(int i=0; i<available_images_proj; i++){
				max_vector_tr[i]=-1;
				max_vector_rt[i]=-1;
				max_vector_tr_mirror[i]=-1;
				max_vector_rt_mirror[i]=-1;
			}

			expIndex = iterExp->objId;
			available_images_exp--;

			SFexp.getRow(rowExp, expIndex);
			rowExp.getValue(MDL_IMAGE, fnImgExp);
			std::cerr << expIndex << ". Image: " << fnImgExp << std::endl;

			//AJ calling the function to align the images
			bool mirror=false;
			align_experimental_image(fnImgExp, d_referenceAux, d_experimentalAux, transMat_tr, transMat_rt,
					max_vector_tr, max_vector_rt, SFexp, available_images_proj, mirror, maxShift,
					myhandlePadded, myhandleMask, myhandlePolar, myhandlePaddedB, myhandleMaskB, myhandlePolarB, myStructureAux);


			////printf("Repeating process with mirror image...\n");
			mirror=true;
			align_experimental_image(fnImgExp, d_referenceAux, d_experimentalAux, transMat_tr_mirror, transMat_rt_mirror,
							max_vector_tr_mirror, max_vector_rt_mirror, SFexp, available_images_proj, mirror, maxShift,
							myhandlePadded, myhandleMask, myhandlePolar, myhandlePaddedB, myhandleMaskB, myhandlePolarB, myStructureAux);

			//AJ to check the best transformation among all the evaluated
			for(int i=0; i<available_images_proj; i++){
				if(max_vector_tr[i]>max_vector_rt[i]){
					transMat_tr.copyOneMatrixToCpu(MULTIDIM_ARRAY(matrixTransCpu[n]), firstIdx+i, i);
					A2D_ELEM(matrixCorrCpu, n, firstIdx+i) = max_vector_tr[i];
				}else{
					transMat_rt.copyOneMatrixToCpu(MULTIDIM_ARRAY(matrixTransCpu[n]), firstIdx+i, i);
					A2D_ELEM(matrixCorrCpu, n, firstIdx+i) = max_vector_rt[i];
				}
				//mirror image
				if(max_vector_tr_mirror[i]>max_vector_rt_mirror[i]){
					transMat_tr_mirror.copyOneMatrixToCpu(MULTIDIM_ARRAY(matrixTransCpu_mirror[n]), firstIdx+i, i);
					A2D_ELEM(matrixCorrCpu_mirror, n, firstIdx+i) = max_vector_tr_mirror[i];
				}else{
					transMat_rt_mirror.copyOneMatrixToCpu(MULTIDIM_ARRAY(matrixTransCpu_mirror[n]), firstIdx+i, i);
					A2D_ELEM(matrixCorrCpu_mirror, n, firstIdx+i) = max_vector_rt_mirror[i];
				}
			}
			/*std::cerr << "Trans Matrix = " << matrixTransCpu[n] << std::endl;
			std::cerr << "Corr Matrix = " << matrixCorrCpu << std::endl;
			std::cerr << "Mirror image - Trans Matrix = " << matrixTransCpu_mirror[n] << std::endl;
			std::cerr << "Mirror image - Corr Matrix = " << matrixCorrCpu_mirror << std::endl;*/


			if(iterExp->hasNext())
				iterExp->moveNext();

			n++;

		}//end while experimental images

		delete iterExp;
		firstIdx +=available_images_proj;
		int aux;
		if(firstIdx+available_images_proj > mdInSize){
			aux=available_images_proj;
			available_images_proj=mdInSize-firstIdx;
		}
		if(firstIdx==mdInSize){
			finish=true;
		}
		if(aux!=available_images_proj){
			myhandlePadded.clear();
			myhandleMask.clear();
			myhandlePolar.clear();
			myhandlePaddedB.clear();
			myhandleMaskB.clear();
			myhandlePolarB.clear();
		}

	}//End loop over the reference images while(!finish)

	myhandlePadded.clear();
	myhandleMask.clear();
	myhandlePolar.clear();
	myhandlePaddedB.clear();
	myhandleMaskB.clear();
	myhandlePolarB.clear();


	MultidimArray<float> weights(1,1,mdExpSize,2*mdInSize);
	MultidimArray<float> corrTotalRow(1,1,mdExpSize, 2*mdInSize);
	int Nref;
	if(keepN){
		Nref=n_keep;
	}else if(significance){
		Nref=round(corrTotalRow.xdim*alpha);
		if(Nref==0)
			Nref=1;
	}
	calculate_weights(matrixCorrCpu, matrixCorrCpu_mirror, corrTotalRow, weights, Nref, mdExpSize, mdInSize);
	std::cerr << "Correlation matrix = " << corrTotalRow << std::endl;
	std::cerr << "Weights = " << weights << std::endl;

    //AJ TIME
    gettimeofday(&end, NULL);
    secs = timeval_diff(&end, &start);
    printf("Parte CUDA: %.16g miliseconds\n", secs * 1000.0);

	generate_metadata(SF, SFexp, fnDir, fn_out, mdExpSize, mdInSize, weights, corrTotalRow, matrixTransCpu, matrixTransCpu_mirror);


	MultidimArray<float> out2(3,3);
	double rot, tilt, psi;

	int *NexpVector;
	if(generate_out){
		size_t xAux, yAux, zAux, nAux;
		getImageSize(SF,xAux,yAux,zAux,nAux);
		FileName fnImgNew, fnExpNew, fnRoot, fnStackOut, fnOut;
		Image<double> Inew, Iexp_aux;
		Matrix2D<double> E(3,3);
		MultidimArray<float> auxtr(3,3);
		MultidimArray<double> refSum(1, 1, yAux, xAux);
		bool firstTime=true;

		CorrelationAux auxCenter;
		RotationalCorrelationAux auxCenter2;

		NexpVector = new int[mdInSize];
		for(int i=0; i<mdInSize; i++){
			NexpVector[i]=0;
			bool change=false;
			double normWeight=0;
			SF.getValue(MDL_IMAGE,fnImgNew,i+1);
			printf("REF: %s\n",fnImgNew.c_str());
			refSum.initZeros();

			fnRoot=fnImgNew.withoutExtension().afterLastOf("/").afterLastOf("@");
			fnStackOut=formatString("%s/%s_aligned.stk",fnDir.c_str(),fnRoot.c_str());
			if(fnStackOut.exists() && firstTime)
				fnStackOut.deleteFile();

			//fnOut = formatString("%s/%s_classes.stk",fnDir.c_str(),fnRoot.c_str());
			firstTime=false;
			for(int j=0; j<mdExpSize; j++){
				long int pointer1=i*xAux*yAux;
				long int pointer2=i*xAux*yAux;

				if(DIRECT_A2D_ELEM(weights,j,i)!=0){
					NexpVector[i]++;
					SFexp.getValue(MDL_IMAGE,fnExpNew,j+1);
					printf("EXP: %s\n",fnExpNew.c_str());
					Iexp_aux.read(fnExpNew);

					matrixTransCpu[j].getSlice(i, auxtr);
					for(int n=0; n<9; n++)
						E.mdata[n]=(double)auxtr.data[n];

					std::cerr << "auxtr = " << auxtr << std::endl;
					std::cerr << "E = " << E << std::endl;

					selfApplyGeometry(LINEAR,Iexp_aux(),E,IS_NOT_INV,WRAP,0.0);



					centerImage(Iexp_aux(), auxCenter, auxCenter2);
					Iexp_aux().resetOrigin();
					refSum += Iexp_aux()*DIRECT_A2D_ELEM(weights,j,i);
					change=true;
					normWeight+=DIRECT_A2D_ELEM(weights,j,i);
				}
				if(DIRECT_A2D_ELEM(weights,j,i+mdInSize)!=0){
					NexpVector[i]++;
					SFexp.getValue(MDL_IMAGE,fnExpNew,j+1);
					Iexp_aux.read(fnExpNew);
					Iexp_aux().selfReverseX();

					matrixTransCpu_mirror[j].getSlice(i, auxtr);
					for(int n=0; n<9; n++)
						E.mdata[n]=(double)auxtr.data[n];

					selfApplyGeometry(LINEAR,Iexp_aux(),E,IS_NOT_INV,WRAP,0.0);

					centerImage(Iexp_aux(), auxCenter, auxCenter2);
					Iexp_aux().resetOrigin();
					refSum += Iexp_aux()*DIRECT_A2D_ELEM(weights,j,i+mdInSize);
					change=true;
					normWeight+=DIRECT_A2D_ELEM(weights,j,i+mdInSize);
				}
			}
			if(change){
				refSum/=normWeight;
				Inew()=refSum;
				Inew.write(fnStackOut,i+1,true,WRITE_APPEND);
			}
		}
	}

	if(generate_out){
		MetaData SFout;
		FileName fnImgNew, fnRoot, fnStackOut, fnStackMD, fnClass;
		bool firstTime=true;
		for(int i=0; i<mdInSize; i++){
			SF.getValue(MDL_IMAGE,fnImgNew,i+1);
			fnRoot=fnImgNew.withoutExtension().afterLastOf("/").afterLastOf("@");
			fnStackOut=formatString("%s/%s_aligned.stk",fnDir.c_str(),fnRoot.c_str());
			fnStackMD=formatString("%s/%s_aligned.xmd",fnDir.c_str(),fnRoot.c_str());
			fnClass.compose(i + 1, fnStackOut);

			if(fnStackMD.exists() && firstTime)
				fnStackMD.deleteFile();

			firstTime=false;
			if(NexpVector[i]==0)
				continue;

			size_t id = SFout.addObject();
			SFout.setValue(MDL_REF, i + 1, id);
			SFout.setValue(MDL_IMAGE, fnClass, id);
			SFout.setValue(MDL_CLASS_COUNT,(size_t)NexpVector[i], id);
		}
		SFout.write("classes@"+fnStackMD, MD_APPEND);

		FileName fnExpIm;
		MDRow row;
		for(int i=0; i<mdInSize; i++){//
			if(NexpVector[i]==0)
				continue;
			MetaData SFq;
			for(int j=0; j<mdExpSize; j++){
				if(DIRECT_A2D_ELEM(weights,j,i)!=0){
					SFexp.getValue(MDL_IMAGE,fnExpIm,j+1);
					row.setValue(MDL_IMAGE, fnExpIm);
					row.setValue(MDL_WEIGHT, (double)DIRECT_A2D_ELEM(weights, j, i));
					matrixTransCpu[j].getSlice(i, out2);
					row.setValue(MDL_FLIP, false);
					row.setValue(MDL_SHIFT_X, (double)DIRECT_A2D_ELEM(out2,0,2));
					row.setValue(MDL_SHIFT_Y, (double)DIRECT_A2D_ELEM(out2,1,2));
					psi = RAD2DEG(atan2(DIRECT_A2D_ELEM(out2,1,0), DIRECT_A2D_ELEM(out2,0,0)));
					SF.getValue(MDL_ANGLE_ROT,rot,i+1);
					row.setValue(MDL_ANGLE_ROT, rot);
					SF.getValue(MDL_ANGLE_TILT,tilt,i+1);
					row.setValue(MDL_ANGLE_TILT, tilt);
					row.setValue(MDL_ANGLE_PSI, psi);
					SFq.addRow(row);
				}
				if(DIRECT_A2D_ELEM(weights,j,i+mdInSize)!=0){
					SFexp.getValue(MDL_IMAGE,fnExpIm,j+1);
					row.setValue(MDL_IMAGE, fnExpIm);
					row.setValue(MDL_WEIGHT, (double)DIRECT_A2D_ELEM(weights, j, i+mdInSize));
					matrixTransCpu_mirror[j].getSlice(i, out2);
					row.setValue(MDL_FLIP, true);
					row.setValue(MDL_SHIFT_X, (double)DIRECT_A2D_ELEM(out2,0,2));
					row.setValue(MDL_SHIFT_Y, (double)DIRECT_A2D_ELEM(out2,1,2));
					psi = RAD2DEG(atan2(DIRECT_A2D_ELEM(out2,1,0), DIRECT_A2D_ELEM(out2,0,0)));
					SF.getValue(MDL_ANGLE_ROT,rot,i+1);
					row.setValue(MDL_ANGLE_ROT, rot);
					SF.getValue(MDL_ANGLE_TILT,tilt,i+1);
					row.setValue(MDL_ANGLE_TILT, tilt);
					row.setValue(MDL_ANGLE_PSI, psi);
					SFq.addRow(row);
				}
			}
			MetaData SFq_sorted;
			SFq_sorted.sort(SFq, MDL_IMAGE);
			SFq_sorted.write(formatString("class%06d_images@%s",i+1,fnStackMD.c_str()),MD_APPEND);

		}

	}



	//Free memory in CPU
	for(int i=0; i<mdExpSize; i++)
		matrixTransCpu[i].coreDeallocate();
	delete []matrixTransCpu;
	delete []max_vector_tr;
	delete []max_vector_rt;
	for(int i=0; i<mdExpSize; i++)
		matrixTransCpu_mirror[i].coreDeallocate();
	delete []matrixTransCpu_mirror;
	delete []max_vector_tr_mirror;
	delete []max_vector_rt_mirror;
	delete []NexpVector;



}
