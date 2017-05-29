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

void preprocess_images(MetaData &SF, int numImages, Mask &mask,
		GpuCorrelationAux &d_correlationAux, bool rotate)
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
	GpuMultidimArrayAtCpu<double> padded_image_stack(pad_Xdim,pad_Ydim,1,numImages);
	GpuMultidimArrayAtCpu<double> padded_image2_stack(pad_Xdim,pad_Ydim,1,numImages);
	GpuMultidimArrayAtCpu<double> padded_mask(pad_Xdim,pad_Ydim);
	GpuMultidimArrayAtCpu<double> autocorrelation_mask(pad_Xdim,pad_Ydim);

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
		Iref().setXmippOrigin();
		mask.apply_mask(Iref(), Iref());
		original_image_stack.fillImage(n,Iref());

		if (rotate)
		{
			Iref().selfReverseX();
			Iref().selfReverseY();
		}

		Iref2=Iref();
		Iref2*=Iref2;
		Iref().window(padIref, STARTINGY(padIref),STARTINGX(padIref),FINISHINGY(padIref),FINISHINGX(padIref));
		Iref2 .window(padIref2,STARTINGY(padIref),STARTINGX(padIref),FINISHINGY(padIref),FINISHINGX(padIref));
		padded_image_stack.fillImage(n,padIref);
		padded_image2_stack.fillImage(n,padIref2);

		if(iter->hasNext())
			iter->moveNext();

		n++;
	}
	delete iter;

	mask.get_binary_mask().window(padMask, STARTINGY(padIref),STARTINGX(padIref),FINISHINGY(padIref),FINISHINGX(padIref));
	padded_mask.fillImage(0, padMask);
	//MultidimArray<double> dMask, maskAutocorrelation;
	//typeCast(padMask, dMask);
	//auto_correlation_matrix(dMask, maskAutocorrelation);
	//autocorrelation_mask.fillImage(0, maskAutocorrelation);

	GpuMultidimArrayAtGpu<double> original_image_gpu, padded_image_gpu, padded_image2_gpu, autocorrelation_mask_gpu, padded_mask_gpu;
	original_image_stack.copyToGpu(original_image_gpu);
	padded_image_stack.copyToGpu(padded_image_gpu);
	padded_image2_stack.copyToGpu(padded_image2_gpu);
	padded_mask.copyToGpu(padded_mask_gpu);
	//autocorrelation_mask.copyToGpu(d_correlationAux.maskAutocorrelation);

	//Polar transform of the projected images
	d_correlationAux.XdimPolar=360;
	d_correlationAux.YdimPolar=radius;
	GpuMultidimArrayAtGpu<double> polar_gpu(360,radius,1,numImages);
	GpuMultidimArrayAtGpu<double> polar2_gpu(360,radius,1,numImages);

	cuda_cart2polar(original_image_gpu, polar_gpu, polar2_gpu, rotate);

	//FFT
	padded_image_gpu.fft(d_correlationAux.d_projFFT);
	padded_image2_gpu.fft(d_correlationAux.d_projSquaredFFT);
	polar_gpu.fft(d_correlationAux.d_projPolarFFT);
	polar2_gpu.fft(d_correlationAux.d_projPolarSquaredFFT);
	padded_mask_gpu.fft(d_correlationAux.d_maskFFT);

	/*/AJ for debugging
	size_t xAux= 360;
	size_t yAux= radius;
	size_t nAux= 1;
	GpuMultidimArrayAtGpu<double> aux(xAux,yAux,1,numImages);
	d_correlationAux.d_maskFFTPolar.ifft(aux);
	GpuMultidimArrayAtCpu<double> auxCpu(xAux,yAux,1,numImages);
	auxCpu.copyFromGpu(aux);
	pointer=0;
	for(int i=0; i<nAux; i++){
	MultidimArray<double> padded;
	FileName fnImgPad;
	Image<double> Ipad;
	padded.coreAllocate(1, 1, yAux, xAux);
	memcpy(MULTIDIM_ARRAY(padded), &auxCpu.data[pointer], xAux*yAux*sizeof(double));
	fnImgPad.compose("mask", i+1, "jpg");
	Ipad()=padded;
	Ipad.write(fnImgPad);
	padded.coreDeallocate();
	pointer += xAux*yAux;
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

//#define DEBUG
// Compute correlation --------------------------------------------------------
void ProgGpuCorrelation::run()
{

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
	//auto_correlation_matrix(dMask, maskAutocorrelation); //AJ con la padded??????????????

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
	preprocess_images(SF, available_images_proj, mask, d_referenceAux, false);
	//d_referenceAux.maskAutocorrelation.resize(XSIZE(maskAutocorrelation),YSIZE(maskAutocorrelation),ZSIZE(maskAutocorrelation),1);
	//fillImage(d_referenceAux.maskAutocorrelation, maskAutocorrelation);
	d_referenceAux.maskCount=maskCount;
	d_referenceAux.produceSideInfo();



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
	preprocess_images(SFexp, available_images_exp, mask, d_experimentalAux, true);
	//d_experimentalAux.maskAutocorrelation.resize(XSIZE(maskAutocorrelation),YSIZE(maskAutocorrelation),ZSIZE(maskAutocorrelation),1);
	//fillImage(d_experimentalAux.maskAutocorrelation, maskAutocorrelation);
	d_experimentalAux.maskCount=maskCount;
	d_experimentalAux.produceSideInfo();

	/*/AJ for debugging
	size_t xAux= d_experimentalAux.d_denom.Xdim;
	size_t yAux= d_experimentalAux.d_denom.Ydim;
	size_t nAux= d_experimentalAux.d_denom.Ndim;
	GpuMultidimArrayAtCpu<double> auxCpu(xAux,yAux,1,nAux);
	auxCpu.copyFromGpu(d_experimentalAux.d_denom);
	int pointer2=0;
	for(int i=0; i<nAux; i++){
	MultidimArray<double> padded;
	FileName fnImgPad;
	Image<double> Ipad;
	padded.coreAllocate(1, 1, yAux, xAux);
	memcpy(MULTIDIM_ARRAY(padded), &auxCpu.data[pointer2], xAux*yAux*sizeof(double));
	fnImgPad.compose("nuevo", i+1, "jpg");
	Ipad()=padded;
	Ipad.write(fnImgPad);
	padded.coreDeallocate();
	pointer2 += xAux*yAux;
	}
	//END AJ/*/


	//CORRELATION PART

	//Translational part
	printf("Calculating correlation...\n");

	bool rotation = false;

	double **NCC;

	if (!rotation)
		NCC = cuda_calculate_correlation(d_referenceAux, d_experimentalAux);
	else
		NCC = cuda_calculate_correlation_rotation(d_referenceAux, d_experimentalAux);

	MultidimArray<double> MDAncc;
	MultidimArray<double> NCC_matrix(available_images_exp, available_images_proj);

	MDIterator *iterExp = new MDIterator(SFexp);
	size_t objIndex = 0;
	size_t objIndexExp = 0;
	int	pointer=0;
	int count=1;

	size_t sizeX, sizeY;
	if (!rotation){
		sizeX = d_experimentalAux.Xdim;
		sizeY = d_experimentalAux.Ydim;
	}else{
		sizeX = d_experimentalAux.XdimPolar;
		sizeY = d_experimentalAux.YdimPolar;
	}

	FileName fnImgPad;
	Image<double> Ipad;
	while(available_images_exp && iterExp->objId!=0)
	{
		objIndexExp = iterExp->objId;
		available_images_exp--;

		MDIterator *iterProj = new MDIterator(SF);
		objIndex = 0;
		pointer = 0;

		int aux_available_images_proj = available_images_proj;

		while(available_images_proj && iterProj->objId!=0){

			objIndex = iterProj->objId;
			available_images_proj--;

			MDAncc.coreAllocate(Ndim, Zdim, sizeY, sizeX);
			memcpy(MULTIDIM_ARRAY(MDAncc), &NCC[iterExp->objId-1][pointer], sizeX*sizeY*sizeof(double));
			fnImgPad.compose("NCC", count, "jpg");
			Ipad()=MDAncc;
			Ipad.write(fnImgPad);

			pointer += (sizeY*sizeX);
			double max = MDAncc.computeMax();
			NCC_matrix(objIndexExp-1, objIndex-1) = max;
			printf("Max=%lf\n",max);

			int posX, posY;
			MDAncc.maxIndex(posY, posX);
			if(posX>sizeX/2 && posY>sizeY/2){
				posX = sizeX-1-posX;
				posY = sizeY-1-posY;
			}else if(posX<sizeX/2 && posY>sizeY/2){
				posX = -(posX+1);
				posY = sizeY-1-posY;
			}else if(posX<sizeX/2 && posY<sizeY/2){
				posX = -(posX+1);
				posY = -(posY+1);
			}else if(posX>sizeX/2 && posY<sizeY/2){
				posX = sizeX-1-posX;
				posY = -(posY+1);
			}
			//TODO: be careful with the sign!!!
			//TODO: los signos en rotacion salen al contrario... hay que establecer una forma de ponerlos...
			printf("Max x-index=%i\n", posX);
			printf("Max y-index=%i\n", posY);

			if(iterProj->hasNext())
				iterProj->moveNext();

			MDAncc.coreDeallocate();
			count++;

		}
		delete iterProj;
		available_images_proj = aux_available_images_proj;

		if(iterExp->hasNext())
			iterExp->moveNext();

	}

	std::cout << "NCC matrix = " << NCC_matrix << std::endl;

	for(int n = 0; n < available_images_exp; n++)
		delete[] NCC[n];
	delete[] NCC;


}

