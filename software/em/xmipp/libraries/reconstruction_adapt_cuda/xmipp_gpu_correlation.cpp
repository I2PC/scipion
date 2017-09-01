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

	MDIterator *iter = new MDIterator(SF);

	size_t n=0;
	for(int i=firstIdx; i<firstIdx+numImages; i++){

		SF.getValue(MDL_IMAGE,fnImg,iter->objId);
		std::cerr << i << ". Image: " << fnImg << std::endl;
		Iref.read(fnImg);
		original_image_stack.fillImage(n,Iref()/8);

		if(iter->hasNext())
			iter->moveNext();

		n++;
	}

	delete iter;

	GpuMultidimArrayAtGpu<float> image_stack_gpu(Xdim,Ydim,1,numImages);
	original_image_stack.copyToGpu(image_stack_gpu);
	MultidimArray<int> maskArray = mask.get_binary_mask();
	MultidimArray<float> dMask;
	typeCast(maskArray, dMask);
	d_correlationAux.d_mask.resize(Xdim, Ydim, Zdim, 1);
	d_correlationAux.d_mask.copyToGpu(MULTIDIM_ARRAY(dMask));

	padding_masking(image_stack_gpu, d_correlationAux.d_mask, myStructureAux.padded_image_gpu, myStructureAux.padded_image2_gpu,
			myStructureAux.padded_mask_gpu, false);

    myStructureAux.padded_image_gpu.fft(d_correlationAux.d_projFFT, myhandlePadded);

    myStructureAux.padded_image2_gpu.fft(d_correlationAux.d_projSquaredFFT, myhandlePadded);
    myStructureAux.padded_mask_gpu.fft(d_correlationAux.d_maskFFT, myhandleMask);

	//Polar transform of the projected images
	cuda_cart2polar(image_stack_gpu, myStructureAux.polar_gpu, myStructureAux.polar2_gpu, false);
    myStructureAux.polar_gpu.fft(d_correlationAux.d_projPolarFFT, myhandlePolar);
    myStructureAux.polar2_gpu.fft(d_correlationAux.d_projPolarSquaredFFT, myhandlePolar);

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

	if(firstStep==0){

		Image<float> Iref;

		GpuMultidimArrayAtCpu<float> original_image_stack(Xdim,Ydim,1,numImagesRef);

		Iref.read(fnImg);

		//AJ mirror of the image
		if(mirror)
			Iref().selfReverseX();
		//END AJ mirror

		for(size_t i=0; i<numImagesRef; i++)
			original_image_stack.fillImage(i,Iref()/8);
		original_image_stack.copyToGpu(d_correlationAux.d_original_image);

	}

	if(!rotation){
		padding_masking(d_correlationAux.d_original_image, mask, myStructureAux.padded_image_gpu, myStructureAux.padded_image2_gpu,
				myStructureAux.padded_mask_gpu, true);
	    myStructureAux.padded_image_gpu.fft(d_correlationAux.d_projFFT, myhandlePadded);
	    myStructureAux.padded_image2_gpu.fft(d_correlationAux.d_projSquaredFFT, myhandlePadded);
		d_maskFFT.copyGpuToGpu(d_correlationAux.d_maskFFT);
	}

	if(rotation){
		cuda_cart2polar(d_correlationAux.d_original_image, myStructureAux.polar_gpu, myStructureAux.polar2_gpu, true);
	    myStructureAux.polar_gpu.fft(d_correlationAux.d_projPolarFFT, myhandlePolar);
	    myStructureAux.polar2_gpu.fft(d_correlationAux.d_projPolarSquaredFFT, myhandlePolar);
	}

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

	if(!rotation){
		padding_masking(d_correlationAux.d_transform_image, mask, myStructureAux.padded_image_gpu, myStructureAux.padded_image2_gpu,
				myStructureAux.padded_mask_gpu, true);
	    myStructureAux.padded_image_gpu.fft(d_correlationAux.d_projFFT, myhandlePadded);
	    myStructureAux.padded_image2_gpu.fft(d_correlationAux.d_projSquaredFFT, myhandlePadded);
		d_maskFFT.copyGpuToGpu(d_correlationAux.d_maskFFT);
	}

	//Polar transform of the projected images
	if(rotation){
		cuda_cart2polar(d_correlationAux.d_transform_image, myStructureAux.polar_gpu, myStructureAux.polar2_gpu, true);
		myStructureAux.polar_gpu.fft(d_correlationAux.d_projPolarFFT, myhandlePolar);
		myStructureAux.polar2_gpu.fft(d_correlationAux.d_projPolarSquaredFFT, myhandlePolar);
	}

}

void align_experimental_image(FileName &fnImgExp, GpuCorrelationAux &d_referenceAux, GpuCorrelationAux &d_experimentalAux,
		TransformMatrix<float> &transMat_tr, TransformMatrix<float> &transMat_rt, float *max_vector_tr, float *max_vector_rt,
		MetaData &SFexp, int available_images_proj, bool mirror, int maxShift,
		mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask, mycufftHandle &myhandlePolar,
		mycufftHandle &myhandlePaddedB, mycufftHandle &myhandleMaskB, mycufftHandle &myhandlePolarB,
		StructuresAux &myStructureAux)
{

	bool rotation;
	TransformMatrix<float> *transMat, *transMat_md;
	float *max_vector;

	for(int firstStep=0; firstStep<2; firstStep++){

		int max_step;
		if (firstStep==0){
			rotation = false;
			max_vector = max_vector_tr;
			max_step=7;
		}else{
			rotation = true;
			max_vector = max_vector_rt;
			max_step=6;
		}

		preprocess_images_experimental(SFexp, fnImgExp, available_images_proj, d_referenceAux.d_mask, d_referenceAux.d_maskFFT,
				d_experimentalAux, rotation, firstStep, mirror, myhandlePadded, myhandleMask, myhandlePolar, myStructureAux);

		if(!rotation){
			d_experimentalAux.maskCount=d_referenceAux.maskCount;
			d_experimentalAux.produceSideInfo(myhandlePaddedB, myhandleMaskB, myStructureAux, d_referenceAux.maskAutocorrelation);
		}

		if(firstStep==0)
			transMat = &transMat_tr;
		else
			transMat = &transMat_rt;

		for(int step=0; step<max_step; step++){ //loop over consecutive translations and rotations (TRTRTRT or RTRTRT)

			//CORRELATION PART
			//TRANSFORMATION MATRIX CALCULATION
			if(!rotation){
				cuda_calculate_correlation(d_referenceAux, d_experimentalAux, *transMat, max_vector, maxShift, myhandlePaddedB,
						mirror, myStructureAux);
			}else{
				cuda_calculate_correlation_rotation(d_referenceAux, d_experimentalAux, *transMat, max_vector, maxShift, myhandlePolarB,
						mirror, myStructureAux);
			}

			//APPLY TRANSFORMATION
			if(step<max_step-1){
				d_experimentalAux.d_transform_image.resize(d_experimentalAux.d_original_image);
				apply_transform(d_experimentalAux.d_original_image, d_experimentalAux.d_transform_image, *transMat);
			}

			//PREPROCESS TO PREPARE DATA TO THE NEXT STEP
			if(step<max_step-1){
				rotation = !rotation;
				preprocess_images_experimental_transform(d_experimentalAux, d_referenceAux.d_mask, d_referenceAux.d_maskFFT, rotation, step,
						myhandlePadded, myhandleMask, myhandlePolar, myStructureAux);

				if(!rotation){
					d_experimentalAux.maskCount=d_referenceAux.maskCount;
					d_experimentalAux.produceSideInfo(myhandlePaddedB, myhandleMaskB, myStructureAux, d_referenceAux.maskAutocorrelation);
				}

			}

		}//end for(int step=0; step<max_step; step++)

	}//end for(int firstStep=0; firstStep<2; firstStep++)

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
	addParamsLine("   [--classify <md_classes_out>]	       : To generate the aligned output images and write the associated metadata");
	addParamsLine("   [--keep_best <N=2>]  			       : To keep N aligned images with the highest correlation");
	addParamsLine("   [--significance <alpha=0.2>]  	   : To use significance with the indicated value");
	addParamsLine("   [--odir <outputDir=\".\">]           : Output directory to save the aligned images");
    addParamsLine("   [--maxShift <s=10>]                  : Maximum shift allowed (+-this amount)");
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
		int idxMax = DIRECT_A2D_ELEM(corrOrderByRowIndex,i,0)-1;
		for(int j=0; j<Nref; j++){
			int idx = DIRECT_A2D_ELEM(corrOrderByRowIndex,i,j)-1;
			if(DIRECT_A2D_ELEM(corrTotalRow,i,idx)<0)
				break;
			float weight = 1.0 - (j/(float)corrOrderByRowIndex.xdim);
			weight *= DIRECT_A2D_ELEM(corrTotalRow,i,idx) / DIRECT_A2D_ELEM(corrTotalRow,i,idxMax);
			DIRECT_A2D_ELEM(weights1, i, idx) = weight;
		}
	}
	for(int i=0; i<mdInSize; i++){
		int idxMax = DIRECT_A2D_ELEM(corrOrderByColIndex,0,i)-1;
		for(int j=0; j<2*mdExpSize; j++){
			int idx = DIRECT_A2D_ELEM(corrOrderByColIndex,j,i)-1;
			float weight = 1.0 - (j/(float)corrOrderByColIndex.ydim);
			weight *= DIRECT_A2D_ELEM(corrTotalCol,idx,i) / DIRECT_A2D_ELEM(corrTotalCol,idxMax,i);
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
		MultidimArray<float> corrTotalRow, MultidimArray<float> *matrixTransCpu, MultidimArray<float> *matrixTransCpu_mirror, int maxShift){

	double maxShift2 = maxShift*maxShift;
	Matrix2D<double> bestM(3,3);
	MultidimArray<float> out2(3,3);
	MDRow rowOut;
	MetaData mdOut;
	String nameImg, nameRef;
	bool flip;
	double rot, tilt, psi;
	int idxJ;

	MDIterator *iterExp = new MDIterator(SFexp);
	MDRow rowExp;
	MDIterator *iter = new MDIterator();
	MDRow row;

	for(int i=0; i<mdExpSize; i++){
		SFexp.getRow(rowExp, iterExp->objId);
		rowExp.getValue(MDL_IMAGE, nameImg);
		iter->init(SF);
		for(int j=0; j<2*mdInSize; j++){

			if(j%mdInSize==0)
				iter->init(SF);
			SF.getRow(row, iter->objId);

			if(DIRECT_A2D_ELEM(weights,i,j)!=0){
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
				double shiftX = (double)DIRECT_A2D_ELEM(out2,0,2);
				double shiftY = (double)DIRECT_A2D_ELEM(out2,1,2);
				if (shiftX*shiftX + shiftY*shiftY > maxShift2){
					if(iter->hasNext())
						iter->moveNext();
					continue;
				}

				rowOut.setValue(MDL_FLIP, flip);

				double scale;
				MAT_ELEM(bestM,0,0)=DIRECT_A2D_ELEM(out2,0,0);
				MAT_ELEM(bestM,0,1)=DIRECT_A2D_ELEM(out2,0,1);
				MAT_ELEM(bestM,0,2)=DIRECT_A2D_ELEM(out2,0,2);

				MAT_ELEM(bestM,1,0)=DIRECT_A2D_ELEM(out2,1,0);
				MAT_ELEM(bestM,1,1)=DIRECT_A2D_ELEM(out2,1,1);
				MAT_ELEM(bestM,1,2)=DIRECT_A2D_ELEM(out2,1,2);

				MAT_ELEM(bestM,2,0)=0.0;
				MAT_ELEM(bestM,2,1)=0.0;
				MAT_ELEM(bestM,2,2)=1.0;
				if(flip){
					MAT_ELEM(bestM,0,0)*=-1;
					MAT_ELEM(bestM,1,0)*=-1;
				}
				bestM=bestM.inv();

				transformationMatrix2Parameters2D(bestM,flip,scale,shiftX,shiftY,psi);
				if (flip)
					shiftX*=-1;
				rowOut.setValue(MDL_SHIFT_X, -shiftX);
				rowOut.setValue(MDL_SHIFT_Y, -shiftY);
				row.getValue(MDL_ANGLE_ROT, rot);
				rowOut.setValue(MDL_ANGLE_ROT, rot);
				row.getValue(MDL_ANGLE_TILT, tilt);
				rowOut.setValue(MDL_ANGLE_TILT, tilt);
				rowOut.setValue(MDL_ANGLE_PSI, psi);

				rowOut.setValue(MDL_REF, idxJ);
				mdOut.addRow(rowOut);
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
		MultidimArray<float> weights, MultidimArray<float> *matrixTransCpu, MultidimArray<float> *matrixTransCpu_mirror,
		int maxShift, FileName fn_classes_out){

	double maxShift2 = maxShift*maxShift;
	MultidimArray<float> out2(3,3);
	double rot, tilt, psi;
	int *NexpVector;

	size_t xAux, yAux, zAux, nAux;
	getImageSize(SF,xAux,yAux,zAux,nAux);
	FileName fnImgNew, fnExpNew, fnRoot, fnStackOut, fnOut, fnStackMD, fnClass;
	Image<double> Inew, Iexp_aux, Inew2;
	Matrix2D<double> E(3,3);
	MultidimArray<float> auxtr(3,3);
	MultidimArray<double> refSum(1, 1, yAux, xAux);
	bool firstTime=true;

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

	NexpVector = new int[mdInSize];
	for(int i=0; i<mdInSize; i++){
		NexpVector[i]=0;
		bool change=false;
		double normWeight=0;

		SF.getRow(rowSF, iterSF->objId);
		rowSF.getValue(MDL_IMAGE, fnImgNew);
		iterSFexp->init(SFexp);

		refSum.initZeros();

		//fnRoot=fnImgNew.withoutExtension().afterLastOf("/").afterLastOf("@");
		fnRoot=fn_classes_out.withoutExtension();
		fnStackOut=formatString("%s/%s.stk",fnDir.c_str(),fnRoot.c_str());
		if(fnStackOut.exists() && firstTime)
			fnStackOut.deleteFile();

		firstTime=false;
		for(int j=0; j<mdExpSize; j++){

			SFexp.getRow(rowSFexp, iterSFexp->objId);
			rowSFexp.getValue(MDL_IMAGE, fnExpNew);
			Iexp_aux.read(fnExpNew);

			long int pointer1=i*xAux*yAux;
			long int pointer2=i*xAux*yAux;

			if(DIRECT_A2D_ELEM(weights,j,i)!=0){

				matrixTransCpu[j].getSlice(i, auxtr);

				double shiftX = (double)DIRECT_A2D_ELEM(auxtr,0,2);
				double shiftY = (double)DIRECT_A2D_ELEM(auxtr,1,2);
				if (shiftX*shiftX + shiftY*shiftY > maxShift2){
					if(iterSFexp->hasNext())
						iterSFexp->moveNext();
					continue;
				}

				NexpVector[i]++;

				MAT_ELEM(E,0,0)=DIRECT_A2D_ELEM(auxtr,0,0);
				MAT_ELEM(E,0,1)=DIRECT_A2D_ELEM(auxtr,0,1);
				MAT_ELEM(E,0,2)=DIRECT_A2D_ELEM(auxtr,0,2);

				MAT_ELEM(E,1,0)=DIRECT_A2D_ELEM(auxtr,1,0);
				MAT_ELEM(E,1,1)=DIRECT_A2D_ELEM(auxtr,1,1);
				MAT_ELEM(E,1,2)=DIRECT_A2D_ELEM(auxtr,1,2);

				MAT_ELEM(E,2,0)=0.0;
				MAT_ELEM(E,2,1)=0.0;
				MAT_ELEM(E,2,2)=1.0;

				selfApplyGeometry(LINEAR,Iexp_aux(),E,IS_NOT_INV,DONT_WRAP,0.0);

				Iexp_aux().resetOrigin();
				refSum += Iexp_aux()*DIRECT_A2D_ELEM(weights,j,i);
				change=true;
				normWeight+=DIRECT_A2D_ELEM(weights,j,i);
			}
			if(DIRECT_A2D_ELEM(weights,j,i+mdInSize)!=0){

				matrixTransCpu_mirror[j].getSlice(i, auxtr);

				double shiftX = (double)DIRECT_A2D_ELEM(auxtr,0,2);
				double shiftY = (double)DIRECT_A2D_ELEM(auxtr,1,2);
				if (shiftX*shiftX + shiftY*shiftY > maxShift2){
					if(iterSFexp->hasNext())
						iterSFexp->moveNext();
					continue;
				}

				NexpVector[i]++;
				Iexp_aux().selfReverseX();

				MAT_ELEM(E,0,0)=DIRECT_A2D_ELEM(auxtr,0,0);
				MAT_ELEM(E,0,1)=DIRECT_A2D_ELEM(auxtr,0,1);
				MAT_ELEM(E,0,2)=DIRECT_A2D_ELEM(auxtr,0,2);

				MAT_ELEM(E,1,0)=DIRECT_A2D_ELEM(auxtr,1,0);
				MAT_ELEM(E,1,1)=DIRECT_A2D_ELEM(auxtr,1,1);
				MAT_ELEM(E,1,2)=DIRECT_A2D_ELEM(auxtr,1,2);

				MAT_ELEM(E,2,0)=0.0;
				MAT_ELEM(E,2,1)=0.0;
				MAT_ELEM(E,2,2)=1.0;

				selfApplyGeometry(LINEAR,Iexp_aux(),E,IS_NOT_INV,DONT_WRAP,0.0);

				Iexp_aux().resetOrigin();
				refSum += Iexp_aux()*DIRECT_A2D_ELEM(weights,j,i+mdInSize);
				change=true;
				normWeight+=DIRECT_A2D_ELEM(weights,j,i+mdInSize);
			}
			if(iterSFexp->hasNext())
				iterSFexp->moveNext();
		}

		if(change){
			refSum/=normWeight;
			Inew()=refSum;
			centerImage(Inew(), auxCenter, auxCenter2);
			//masking to avoid wrapping in the edges of the image
			mask.apply_mask(Inew(), Inew2());
			Inew2().resetOrigin();
			Inew2.write(fnStackOut,i+1,true,WRITE_APPEND);
		}

		if(iterSF->hasNext())
			iterSF->moveNext();
	}

	iterSFexp->init(SFexp);
	iterSF->init(SF);

	Matrix2D<double> bestM(3,3);
	MetaData SFout;
	firstTime=true;
	for(int i=0; i<mdInSize; i++){

		SF.getRow(rowSF, iterSF->objId);
		rowSF.getValue(MDL_IMAGE, fnImgNew);
		//fnRoot=fnImgNew.withoutExtension().afterLastOf("/").afterLastOf("@");
		fnRoot=fn_classes_out.withoutExtension();
		fnStackMD=formatString("%s/%s.xmd",fnDir.c_str(),fnRoot.c_str());
		fnClass.compose(i + 1, fnStackOut);

		if(fnStackMD.exists() && firstTime)
			fnStackMD.deleteFile();

		firstTime=false;
		if(NexpVector[i]==0){
			if(iterSF->hasNext())
				iterSF->moveNext();
			continue;
		}

		size_t id = SFout.addObject();
		SFout.setValue(MDL_REF, i + 1, id);
		SFout.setValue(MDL_IMAGE, fnClass, id);
		SFout.setValue(MDL_CLASS_COUNT,(size_t)NexpVector[i], id);

		if(iterSF->hasNext())
			iterSF->moveNext();
	}
	SFout.write("classes@"+fnStackMD, MD_APPEND);

	iterSF->init(SF);
	FileName fnExpIm;
	MDRow row;
	for(int i=0; i<mdInSize; i++){
		SF.getRow(rowSF, iterSF->objId);
		iterSFexp->init(SFexp);
		if(NexpVector[i]==0){
			if(iterSF->hasNext())
				iterSF->moveNext();
			continue;
		}
		MetaData SFq;
		for(int j=0; j<mdExpSize; j++){
			SFexp.getRow(rowSFexp, iterSFexp->objId);
			rowSFexp.getValue(MDL_IMAGE, fnExpIm);
			if(DIRECT_A2D_ELEM(weights,j,i)!=0){

				matrixTransCpu[j].getSlice(i, out2);

				double sx = (double)DIRECT_A2D_ELEM(out2,0,2);
				double sy = (double)DIRECT_A2D_ELEM(out2,1,2);
				if (sx*sx + sy*sy > maxShift2){
					if(iterSFexp->hasNext())
						iterSFexp->moveNext();
					continue;
				}

				row.setValue(MDL_IMAGE, fnExpIm);
				row.setValue(MDL_WEIGHT, (double)DIRECT_A2D_ELEM(weights, j, i));
				row.setValue(MDL_FLIP, false);

				double scale, shiftX, shiftY, psi;
				bool flip;
				MAT_ELEM(bestM,0,0)=DIRECT_A2D_ELEM(out2,0,0);
				MAT_ELEM(bestM,0,1)=DIRECT_A2D_ELEM(out2,0,1);
				MAT_ELEM(bestM,0,2)=DIRECT_A2D_ELEM(out2,0,2);

				MAT_ELEM(bestM,1,0)=DIRECT_A2D_ELEM(out2,1,0);
				MAT_ELEM(bestM,1,1)=DIRECT_A2D_ELEM(out2,1,1);
				MAT_ELEM(bestM,1,2)=DIRECT_A2D_ELEM(out2,1,2);

				MAT_ELEM(bestM,2,0)=0.0;
				MAT_ELEM(bestM,2,1)=0.0;
				MAT_ELEM(bestM,2,2)=1.0;
				bestM=bestM.inv();

				transformationMatrix2Parameters2D(bestM,flip,scale,shiftX,shiftY,psi);

				row.setValue(MDL_SHIFT_X, -shiftX);
				row.setValue(MDL_SHIFT_Y, -shiftY);
				rowSF.getValue(MDL_ANGLE_ROT, rot);
				row.setValue(MDL_ANGLE_ROT, rot);
				rowSF.getValue(MDL_ANGLE_TILT, tilt);
				row.setValue(MDL_ANGLE_TILT, tilt);
				row.setValue(MDL_ANGLE_PSI, psi);
				SFq.addRow(row);
			}
			if(DIRECT_A2D_ELEM(weights,j,i+mdInSize)!=0){

				matrixTransCpu_mirror[j].getSlice(i, out2);

				double sx = (double)DIRECT_A2D_ELEM(out2,0,2);
				double sy = (double)DIRECT_A2D_ELEM(out2,1,2);
				if (sx*sx + sy*sy > maxShift2){
					if(iterSFexp->hasNext())
						iterSFexp->moveNext();
					continue;
				}

				row.setValue(MDL_IMAGE, fnExpIm);
				row.setValue(MDL_WEIGHT, (double)DIRECT_A2D_ELEM(weights, j, i+mdInSize));
				row.setValue(MDL_FLIP, true);

				double scale, shiftX, shiftY, psi;
				bool flip;
				MAT_ELEM(bestM,0,0)=DIRECT_A2D_ELEM(out2,0,0);
				MAT_ELEM(bestM,0,1)=DIRECT_A2D_ELEM(out2,0,1);
				MAT_ELEM(bestM,0,2)=DIRECT_A2D_ELEM(out2,0,2);

				MAT_ELEM(bestM,1,0)=DIRECT_A2D_ELEM(out2,1,0);
				MAT_ELEM(bestM,1,1)=DIRECT_A2D_ELEM(out2,1,1);
				MAT_ELEM(bestM,1,2)=DIRECT_A2D_ELEM(out2,1,2);

				MAT_ELEM(bestM,2,0)=0.0;
				MAT_ELEM(bestM,2,1)=0.0;
				MAT_ELEM(bestM,2,2)=1.0;

				MAT_ELEM(bestM,0,0)*=-1;
				MAT_ELEM(bestM,1,0)*=-1;
				bestM=bestM.inv();

				transformationMatrix2Parameters2D(bestM,flip,scale,shiftX,shiftY,psi);

				shiftX*=-1;
				row.setValue(MDL_SHIFT_X, -shiftX);
				row.setValue(MDL_SHIFT_Y, -shiftX);
				rowSF.getValue(MDL_ANGLE_ROT, rot);
				row.setValue(MDL_ANGLE_ROT, rot);
				rowSF.getValue(MDL_ANGLE_TILT, tilt);
				row.setValue(MDL_ANGLE_TILT, tilt);
				row.setValue(MDL_ANGLE_PSI, psi);
				SFq.addRow(row);
			}
			if(iterSFexp->hasNext())
				iterSFexp->moveNext();
		}
		MetaData SFq_sorted;
		SFq_sorted.sort(SFq, MDL_IMAGE);
		SFq_sorted.write(formatString("class%06d_images@%s",i+1,fnStackMD.c_str()),MD_APPEND);

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

	//PROJECTION IMAGES PART
	size_t Xdim, Ydim, Zdim, Ndim;
	SF.read(fn_ref,NULL);
	size_t mdInSize = SF.size();
	getImageSize(SF, Xdim, Ydim, Zdim, Ndim);

	// Generate mask
	Mask mask;
    mask.type = BINARY_CIRCULAR_MASK;
	mask.mode = INNER_MASK;
	size_t rad = (size_t)std::min(Xdim*0.48, Ydim*0.48);
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
	if(Xdim*Ydim*mdInSize*4*100/memory[1]>limit){
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


		preprocess_images_reference(SF, firstIdx, available_images_proj, mask, d_referenceAux,
				myhandlePadded, myhandleMask, myhandlePolar, myhandleAux, myStructureAux);

	    d_referenceAux.maskCount=maskCount;
		d_referenceAux.produceSideInfo(myhandlePaddedB, myhandleMaskB, myStructureAux);


		//EXPERIMENTAL IMAGES PART
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
					myhandlePadded, myhandleMask, myhandlePolar, myhandlePaddedB, myhandleMaskB, myhandlePolarB,
					myStructureAux);


			mirror=true;
			align_experimental_image(fnImgExp, d_referenceAux, d_experimentalAux, transMat_tr_mirror, transMat_rt_mirror,
							max_vector_tr_mirror, max_vector_rt_mirror, SFexp, available_images_proj, mirror, maxShift,
							myhandlePadded, myhandleMask, myhandlePolar, myhandlePaddedB, myhandleMaskB, myhandlePolarB,
							myStructureAux);


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

	generate_metadata(SF, SFexp, fnDir, fn_out, mdExpSize, mdInSize, weights, corrTotalRow, matrixTransCpu,
			matrixTransCpu_mirror, maxShift);

	if(generate_out)
		generate_output_classes(SF, SFexp, fnDir, mdExpSize, mdInSize, weights, matrixTransCpu,
				matrixTransCpu_mirror, maxShift, fn_classes_out);

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
