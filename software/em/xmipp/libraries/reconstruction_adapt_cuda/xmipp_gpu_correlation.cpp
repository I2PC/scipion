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

#include <math.h>


void preprocess_images_reference(MetaData &SF, int firstIdx, int numImages, Mask &mask,
		GpuCorrelationAux &d_correlationAux)
{
	size_t Xdim, Ydim, Zdim, Ndim;
	getImageSize(SF,Xdim,Ydim,Zdim,Ndim);
	size_t pad_Xdim=2*Xdim-1;
	size_t pad_Ydim=2*Ydim-1;

	d_correlationAux.Xdim=pad_Xdim;
	d_correlationAux.Ydim=pad_Ydim;

	MDRow rowIn;
	FileName fnImg;
	Image<double> Iref;
	size_t radius=(size_t)mask.R1;

	GpuMultidimArrayAtCpu<double> original_image_stack(Xdim,Ydim,1,numImages);

	size_t n=0;
	for(int i=firstIdx; i<firstIdx+numImages; i++){

		SF.getValue(MDL_IMAGE,fnImg,i+1);
		std::cerr << i << ". Image: " << fnImg << std::endl;
		Iref.read(fnImg);
		original_image_stack.fillImage(n,Iref());

		n++;
	}


	//AJ new masking and padding
	//original_image_stack.copyToGpu(d_correlationAux.d_original_image);
	GpuMultidimArrayAtGpu<double> image_stack_gpu(Xdim,Ydim,1,numImages);
	original_image_stack.copyToGpu(image_stack_gpu);
	MultidimArray<int> maskArray = mask.get_binary_mask();
	MultidimArray<double> dMask;
	typeCast(maskArray, dMask);
	GpuMultidimArrayAtGpu<double> mask_device(Xdim, Ydim, Zdim, 1);
	mask_device.copyToGpu(MULTIDIM_ARRAY(dMask));

	GpuMultidimArrayAtGpu<double> padded_image_gpu, padded_image2_gpu, padded_mask_gpu;
	padded_image_gpu.resize(pad_Xdim, pad_Ydim, 1, numImages);
	padded_image2_gpu.resize(pad_Xdim, pad_Ydim, 1, numImages);
	padded_mask_gpu.resize(pad_Xdim, pad_Ydim, 1, 1);

	padding_masking(image_stack_gpu, mask_device, padded_image_gpu, padded_image2_gpu,
			padded_mask_gpu, NULL, false);

	padded_image_gpu.fft(d_correlationAux.d_projFFT);
	padded_image2_gpu.fft(d_correlationAux.d_projSquaredFFT);
	padded_mask_gpu.fft(d_correlationAux.d_maskFFT);

	//Polar transform of the projected images
	d_correlationAux.XdimPolar=360;
	d_correlationAux.YdimPolar=radius;
	GpuMultidimArrayAtGpu<double> polar_gpu(360,radius,1,numImages);
	GpuMultidimArrayAtGpu<double> polar2_gpu(360,radius,1,numImages);
	cuda_cart2polar(image_stack_gpu, polar_gpu, polar2_gpu, false); //d_correlationAux.d_original_image
	//FFT
	polar_gpu.fft(d_correlationAux.d_projPolarFFT);
	polar2_gpu.fft(d_correlationAux.d_projPolarSquaredFFT);

	/*/AJ for debugging
	size_t xAux1= 360;
	size_t yAux1= radius;
	size_t nAux1= d_correlationAux.d_projPolarFFT.Ndim;
	GpuMultidimArrayAtGpu<double> aux(xAux1,yAux1,1,nAux1);
	d_correlationAux.d_projPolarFFT.ifft(aux);
	GpuMultidimArrayAtCpu<double> auxCpu1(xAux1,yAux1,1,nAux1);
	auxCpu1.copyFromGpu(aux);
	int pointer1=0;
	for(int i=0; i<nAux1; i++){
	MultidimArray<double> padded1;
	FileName fnImgPad1;
	Image<double> Ipad1;
	padded1.coreAllocate(1, 1, yAux1, xAux1);
	memcpy(MULTIDIM_ARRAY(padded1), &auxCpu1.data[pointer1], xAux1*yAux1*sizeof(double));
	fnImgPad1.compose("./kk/70set", i+1, "mrc");
	Ipad1()=padded1;
	Ipad1.write(fnImgPad1);
	padded1.coreDeallocate();
	pointer1 += xAux1*yAux1;
	}
	//END AJ/*/

}



void preprocess_images_experimental(MetaData &SF, FileName &fnImg, int numImagesRef, Mask &mask,
		GpuCorrelationAux &d_correlationAux, bool rotation, int firstStep, bool mirror)
{
	size_t Xdim, Ydim, Zdim, Ndim;
	getImageSize(SF,Xdim,Ydim,Zdim,Ndim);
	size_t pad_Xdim=2*Xdim-1;
	size_t pad_Ydim=2*Ydim-1;
	size_t radius=(size_t)mask.R1;

	GpuMultidimArrayAtGpu<double> image_stack_gpu(Xdim,Ydim,1,numImagesRef);

	if(firstStep==0){

		d_correlationAux.Xdim=pad_Xdim;
		d_correlationAux.Ydim=pad_Ydim;

		Image<double> Iref;

		GpuMultidimArrayAtCpu<double> original_image_stack(Xdim,Ydim,1,numImagesRef);

		Iref.read(fnImg);

		//AJ adding noise to the image
		//Iref().addNoise(0, 2000, "gaussian");
		//END AJ adding noise

		//AJ mirror of the image
		if(mirror)
			Iref().selfReverseX();
		//END AJ mirror

		for(size_t i=0; i<numImagesRef; i++)
			original_image_stack.fillImage(i,Iref());
		original_image_stack.copyToGpu(d_correlationAux.d_original_image);
		original_image_stack.copyToGpu(image_stack_gpu);

	}else{
		d_correlationAux.d_original_image.copyGpuToGpu(image_stack_gpu);
	}

	//AJ new masking and padding
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

		padding_masking(image_stack_gpu, mask_device, padded_image_gpu, padded_image2_gpu,
			padded_mask_gpu, rotation, true);

		//FFT
		padded_image_gpu.fft(d_correlationAux.d_projFFT);
		padded_image2_gpu.fft(d_correlationAux.d_projSquaredFFT);
		padded_mask_gpu.fft(d_correlationAux.d_maskFFT);
	}

	GpuMultidimArrayAtGpu<double> polar_gpu, polar2_gpu;
	if(rotation){
		d_correlationAux.XdimPolar=360;
		d_correlationAux.YdimPolar=radius;
		polar_gpu.resize(360,radius,1,numImagesRef);
		polar2_gpu.resize(360,radius,1,numImagesRef);
		cuda_cart2polar(d_correlationAux.d_original_image, polar_gpu, polar2_gpu, true);
		//FFT
		polar_gpu.fft(d_correlationAux.d_projPolarFFT);
		polar2_gpu.fft(d_correlationAux.d_projPolarSquaredFFT);
	}

}


void preprocess_images_experimental_transform(GpuCorrelationAux &d_correlationAux, Mask &mask, bool rotation, int step)
{

	size_t Xdim = d_correlationAux.d_transform_image.Xdim;
	size_t Ydim = d_correlationAux.d_transform_image.Ydim;
	size_t Zdim = d_correlationAux.d_transform_image.Zdim;
	size_t Ndim = d_correlationAux.d_transform_image.Ndim;
	size_t pad_Xdim=2*Xdim-1;
	size_t pad_Ydim=2*Ydim-1;
	size_t radius=(size_t)mask.R1;

	//GpuMultidimArrayAtGpu<double> image_stack_gpu(Xdim,Ydim,1,Ndim);
	//d_correlationAux.d_transform_image.copyGpuToGpu(image_stack_gpu);

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

		padding_masking(d_correlationAux.d_transform_image, mask_device, padded_image_gpu, padded_image2_gpu,
				padded_mask_gpu, rotation, true);

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
	if(rotation){
	size_t xAux= 360;
	size_t yAux= radius;
	size_t nAux= 1;
	//GpuMultidimArrayAtGpu<double> aux(xAux,yAux,1,nAux);
	//d_correlationAux.d_projPolarFFT.ifft(aux);
	GpuMultidimArrayAtCpu<double> auxCpu(xAux,yAux,1,nAux);
	auxCpu.copyFromGpu(polar_gpu);
	int pointer2=0;
	for(int i=0; i<nAux; i++){
	MultidimArray<double> padded;
	FileName fnImgPad;
	Image<double> Ipad;
	padded.coreAllocate(1, 1, yAux, xAux);
	memcpy(MULTIDIM_ARRAY(padded), &auxCpu.data[pointer2], xAux*yAux*sizeof(double));
	fnImgPad.compose("polarBad", step+1, "mrc");
	Ipad()=padded;
	Ipad.write(fnImgPad);
	padded.coreDeallocate();
	pointer2 += xAux*yAux;
	}
	}
	//END AJ/*/


}

void align_experimental_image(FileName &fnImgExp, GpuCorrelationAux &d_referenceAux, GpuCorrelationAux &d_experimentalAux,
		TransformMatrix<float> &transMat_tr, TransformMatrix<float> &transMat_rt, double *max_vector_tr, double *max_vector_rt,
		MetaData &SFexp, int available_images_proj, Mask &mask, bool mirror, double maxShift)
{

	bool rotation;
	TransformMatrix<float> *transMat;
	double *max_vector;

	for(int firstStep=0; firstStep<2; firstStep++){ //2

		if (firstStep==0){
			rotation = false;
			printf("First step\n");
			max_vector = max_vector_tr;
		}else{
			rotation = true;
			printf("Second step\n");
			max_vector = max_vector_rt;
		}

		preprocess_images_experimental(SFexp, fnImgExp, available_images_proj, mask, d_experimentalAux, rotation, firstStep, mirror);

		if(!rotation){
			d_experimentalAux.maskCount=d_referenceAux.maskCount;
			d_experimentalAux.produceSideInfo();
		}

		if(firstStep==0)
			transMat = &transMat_tr;
		else
			transMat = &transMat_rt;

		int max_step=6;
		char stepchar[20]="";
		for(int step=0; step<max_step; step++){ //loop over consecutive translations and rotations (TRTRTR or RTRTRT) 6

			if(!rotation){
				stepchar[step]='T';
				printf("step %i of %i %s\n",step+1, max_step, stepchar);
			}else{
				stepchar[step]='R';
				printf("step %i of %i %s\n",step+1, max_step, stepchar);
			}

			//CORRELATION PART
			//TRANSFORMATION MATRIX CALCULATION
			//printf("Calculating correlation...\n");
			if(!rotation)
				cuda_calculate_correlation(d_referenceAux, d_experimentalAux, *transMat, max_vector, maxShift);
			else
				cuda_calculate_correlation_rotation(d_referenceAux, d_experimentalAux, *transMat, max_vector);


			//APPLY TRANSFORMATION
			if(step<max_step-1){
				//printf("Applying transformation...\n");
				d_experimentalAux.d_transform_image.resize(d_experimentalAux.d_original_image);
				apply_transform(d_experimentalAux.d_original_image, d_experimentalAux.d_transform_image, *transMat);
			}


			/*/AJ for debugging
			if(firstStep==0 && !mirror && step==0){
			size_t xAux= d_experimentalAux.debug.Xdim;
			size_t yAux= d_experimentalAux.debug.Ydim;
			size_t nAux= d_experimentalAux.debug.Ndim;
			GpuMultidimArrayAtCpu<double> auxCpu(xAux,yAux,1,nAux);
			auxCpu.copyFromGpu(d_experimentalAux.debug);
			int pointer2=0;
			for(int i=0; i<nAux; i++){
			MultidimArray<double> padded;
			FileName fnImgPad;
			Image<double> Ipad;
			padded.coreAllocate(1, 1, yAux, xAux);
			memcpy(MULTIDIM_ARRAY(padded), &auxCpu.data[pointer2], xAux*yAux*sizeof(double));
			fnImgPad.compose("./kk/70set", i+1, "mrc");
			Ipad()=padded;
			Ipad.write(fnImgPad);
			padded.coreDeallocate();
			pointer2 += xAux*yAux;
			}
			}
			//END AJ/*/

			//PREPROCESS TO PREPARE DATA TO THE NEXT STEP
			if(step<max_step-1){
				//printf("Re-processing experimental images...\n");
				rotation = !rotation;
				preprocess_images_experimental_transform(d_experimentalAux, mask, rotation, step);
				if(!rotation){
					d_experimentalAux.maskCount=d_referenceAux.maskCount;
					d_experimentalAux.produceSideInfo();
				}
			}

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
   	maxShift = getDoubleParam("--maxShift");

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

// Compute correlation --------------------------------------------------------
void ProgGpuCorrelation::run()
{


	//PROJECTION IMAGES PART
	printf("Reference images: \n");
	//Read input metadataFile for projection images
	size_t Xdim, Ydim, Zdim, Ndim;
	SF.read(fn_ref,NULL);
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

	//AJ check the size of the data to avoid exceed the GPU memory
	float memory[3]={0, 0, 0}; //total, free, used
	cuda_check_gpu_memory(memory);
	//printf("START Total %f   Free %f   Used %f\n", memory[0], memory[1], memory[2]);

    int maxGridSize[3];
    cuda_check_gpu_properties(maxGridSize);
    //printf("maxGridSize[0] %i, maxGridSize[1] %i, maxGridSize[2] %i \n",maxGridSize[0],maxGridSize[1],maxGridSize[2]);


	//AJ check_gpu_memory to know how many images we can copy in the gpu memory
	int available_images_proj = mdInSize;
	if(Xdim*Ydim*mdInSize*8*100/memory[1]>2.4){ //TODO debe depender de la memoria disponible en la GPU
		available_images_proj = floor(memory[1]*0.024/(Xdim*Ydim*8));
	}
	if(Xdim*2*Ydim*2*mdInSize>maxGridSize[0]){
		available_images_proj = floor((round(maxGridSize[0]*0.9))/(Xdim*Ydim));
	}
	//available_images_proj=2;
	printf("available_images_proj %i\n",available_images_proj);

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
	MultidimArray<double> matrixCorrCpu(1, 1, mdExpSize, mdInSize);
	MultidimArray<double> matrixCorrCpu_mirror(1, 1, mdExpSize, mdInSize);

	//Auxiliar vectors with maximum values of correlation in RT and TR steps
	double *max_vector_rt = new double [available_images_proj];
	double *max_vector_tr = new double [available_images_proj];
	double *max_vector_rt_mirror = new double [available_images_proj];
	double *max_vector_tr_mirror = new double [available_images_proj];

	//Transformation matrix in GPU
	TransformMatrix<float> transMat_tr(available_images_proj);
	TransformMatrix<float> transMat_rt(available_images_proj);
	TransformMatrix<float> transMat_tr_mirror(available_images_proj);
	TransformMatrix<float> transMat_rt_mirror(available_images_proj);


	int firstIdx=0;
	bool finish=false;

	//Loop over the reference images
	while(!finish){


		GpuCorrelationAux d_referenceAux;
		preprocess_images_reference(SF, firstIdx, available_images_proj, mask, d_referenceAux);
		d_referenceAux.maskCount=maskCount;
		d_referenceAux.produceSideInfo();

		//AJ check the size of the data to avoid exceed the CUDA FFT size
		//float memory[3]={0, 0, 0}; //total, free, used
		//cuda_check_gpu_memory(memory);
		//printf("Total %f   Free %f   Used %f\n", memory[0], memory[1], memory[2]);


		//EXPERIMENTAL IMAGES PART
		printf("Experimental image: \n");

		size_t expIndex = 0;
		MDRow rowExp;
		FileName fnImgExp;
		MDIterator *iterExp = new MDIterator(SFexp);

		GpuCorrelationAux d_experimentalAux;

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
					max_vector_tr, max_vector_rt, SFexp, available_images_proj, mask, mirror, maxShift);


			printf("Repeating process with mirror image...\n");
			mirror=true;
			align_experimental_image(fnImgExp, d_referenceAux, d_experimentalAux, transMat_tr_mirror, transMat_rt_mirror,
							max_vector_tr_mirror, max_vector_rt_mirror, SFexp, available_images_proj, mask, mirror, maxShift);

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
			std::cerr << "Trans Matrix = " << matrixTransCpu[n] << std::endl;
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
		if(firstIdx+available_images_proj > mdInSize){
			available_images_proj=mdInSize-firstIdx;
		}
		if(firstIdx==mdInSize){
			finish=true;
		}

	}//End loop over the reference images while(!finish)


	//Concatenate corr matrix of the original and mirror images in both directions
	//v2 = v1.sort();
	//v2 = v1.indexSort(); //indexes start at 1
	MultidimArray<double> corrTotalRow(1,1,mdExpSize, 2*mdInSize);
	MultidimArray<double> colAux;
	for(int i=0; i<2*mdInSize; i++){
		if(i<mdInSize){
			matrixCorrCpu.getCol(i,colAux);
			corrTotalRow.setCol(i, colAux);
		}else{
			matrixCorrCpu_mirror.getCol(i-mdInSize,colAux);
			corrTotalRow.setCol(i, colAux);
		}
	}
	std::cerr << "Correlation matrix = " << corrTotalRow << std::endl;

	MultidimArray<double> corrTotalCol(1,1,2*mdExpSize, mdInSize);
	MultidimArray<double> rowAux;
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
	MultidimArray<double> rowCorr;
	MultidimArray<int> rowIndexOrder;
	MultidimArray<int> corrOrderByRowIndex(1,1,mdExpSize, 2*mdInSize);

	MultidimArray<double> colCorr;
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

	int Nref;
	if(keepN){
		Nref=n_keep;
	}else if(significance){
		Nref=round(corrOrderByRowIndex.xdim*alpha);
		if(Nref==0)
			Nref=1;
	}

	//AJ To calculate the weights of every image
	//MultidimArray<double> out(3,3);
	MultidimArray<double> weights1(1,1,mdExpSize,2*mdInSize);
	MultidimArray<double> weights2(1,1,mdExpSize,2*mdInSize);
	MultidimArray<double> weights(1,1,mdExpSize,2*mdInSize);

	for(int i=0; i<mdExpSize; i++){
		for(int j=0; j<Nref; j++){
			int idx = DIRECT_A2D_ELEM(corrOrderByRowIndex,i,j)-1;
			if(DIRECT_A2D_ELEM(corrTotalRow,i,idx)<0)
				break;
			double weight = 1.0 - (j/(double)corrOrderByRowIndex.xdim);
			DIRECT_A2D_ELEM(weights1, i, idx) = weight;
		}
	}

	for(int i=0; i<mdInSize; i++){
		for(int j=0; j<2*mdExpSize; j++){
			int idx = DIRECT_A2D_ELEM(corrOrderByColIndex,j,i)-1;
			double weight = 1.0 - (j/(double)corrOrderByColIndex.ydim);
			if(idx<mdExpSize){
				DIRECT_A2D_ELEM(weights2, idx, i) = weight;
			}else{
				DIRECT_A2D_ELEM(weights2, idx-mdExpSize, i+mdInSize) = weight;
			}
		}
	}
	weights=weights1*weights2;
	std::cerr << "Weights = " << weights << std::endl;

	//AJ To generate the output metadata
	MultidimArray<double> out2(3,3);
	MDRow rowOut;
	MetaData mdOut;
	String nameImg;
	bool flip;
	double rot, tilt, psi;
	size_t count=0;
	for(int i=0; i<mdExpSize; i++){
		for(int j=0; j<2*mdInSize; j++){
			if(DIRECT_A2D_ELEM(weights,i,j)!=0){
				rowOut.setValue(MDL_ITEM_ID, count);
				SFexp.getValue(MDL_IMAGE,nameImg,i+1);
				rowOut.setValue(MDL_IMAGE,nameImg);
				rowOut.setValue(MDL_WEIGHT, DIRECT_A2D_ELEM(weights, i, j));
				if(j<mdInSize){
					flip = false;
					matrixTransCpu[i].getSlice(j, out2);
				}else{
					flip = true;
					matrixTransCpu_mirror[i].getSlice(j-mdInSize, out2);
				}
				rowOut.setValue(MDL_SHIFT_X, DIRECT_A2D_ELEM(out2,0,2));
				rowOut.setValue(MDL_SHIFT_Y, DIRECT_A2D_ELEM(out2,1,2));
				psi = RAD2DEG(atan2(DIRECT_A2D_ELEM(out2,1,0), DIRECT_A2D_ELEM(out2,0,0)));
				rowOut.setValue(MDL_ANGLE_PSI, psi);
				if(significance){
					SF.getValue(MDL_ANGLE_ROT,rot,j+1);
					rowOut.setValue(MDL_ANGLE_ROT, rot);
					SF.getValue(MDL_ANGLE_TILT,tilt,j+1);
					rowOut.setValue(MDL_ANGLE_TILT, tilt);
				}
				rowOut.setValue(MDL_FLIP, flip);
				mdOut.addRow(rowOut);
				count++;
			}
		}
	}
	String fnFinal=formatString("%s/%s",fnDir.c_str(),fn_out.c_str());
	mdOut.write(fnFinal);


	if(generate_out){
		size_t xAux, yAux, zAux, nAux;
		getImageSize(SF,xAux,yAux,zAux,nAux);
		String fnImgNew, fnExpNew;
		Image<double> Inew, Iexp_aux;
		Matrix2D<double> E(3,3);
		MultidimArray<float> auxtr(3,3);
		MultidimArray<double> refSum(1, 1, yAux, xAux);

		for(int i=0; i<mdInSize; i++){
			bool change=false;
			double normWeight=0;
			SF.getValue(MDL_IMAGE,fnImgNew,i+1);
			refSum.initZeros();

			for(int j=0; j<mdExpSize; j++){
				int pointer1=i*xAux*yAux;
				int pointer2=i*xAux*yAux;

				if(DIRECT_A2D_ELEM(weights,j,i)!=0){
					SFexp.getValue(MDL_IMAGE,fnExpNew,j+1);
					Iexp_aux.read(fnExpNew);

					matrixTransCpu[j].getSlice(i, auxtr);
					for(int n=0; n<9; n++)
						E.mdata[n]=(double)auxtr.data[n];

					selfApplyGeometry(LINEAR,Iexp_aux(),E,IS_NOT_INV,WRAP,0.0);

					refSum += Iexp_aux()*DIRECT_A2D_ELEM(weights,j,i);
					change=true;
					normWeight+=DIRECT_A2D_ELEM(weights,j,i);
				}
				if(DIRECT_A2D_ELEM(weights,j,i+mdInSize)!=0){
					SFexp.getValue(MDL_IMAGE,fnExpNew,j+1);
					Iexp_aux.read(fnExpNew);
					Iexp_aux().selfReverseX();
					Iexp_aux().resetOrigin(); //TODO AJ preguntar esto

					matrixTransCpu_mirror[j].getSlice(i, auxtr);
					for(int n=0; n<9; n++)
						E.mdata[n]=(double)auxtr.data[n];

					selfApplyGeometry(LINEAR,Iexp_aux(),E,IS_NOT_INV,WRAP,0.0);

					refSum += Iexp_aux()*DIRECT_A2D_ELEM(weights,j,i+mdInSize);
					change=true;
					normWeight+=DIRECT_A2D_ELEM(weights,j,i+mdInSize);
				}
			}
			if(change){
				refSum/=normWeight;
				Inew()=refSum;
				fnImgNew.insert(fnImgNew.rfind("."),"_aligned");
				String fnImFinal=formatString("%s/%s",fnDir.c_str(),fnImgNew.substr(fnImgNew.rfind("/")).c_str());
				Inew.write(fnImFinal);
			}
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

}

