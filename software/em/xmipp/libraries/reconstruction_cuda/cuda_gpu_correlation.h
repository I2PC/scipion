
#ifndef CUDA_GPU_CORRELATION_H
#define CUDA_GPU_CORRELATION_H
#include <complex>
#include "cuda_xmipp_utils.h"

void cuda_cart2polar(GpuMultidimArrayAtGpu<float> &image, GpuMultidimArrayAtGpu<float> &polar_image, GpuMultidimArrayAtGpu<float> &polar2_image, bool rotate, myStreamHandle &myStream);
void waitGpu(myStreamHandle &myStream, bool allStreams);
void calculateAbs (std::complex<float> *data, float *out, int size, myStreamHandle &myStream);


class StructuresAux {
public:
	GpuMultidimArrayAtGpu<float> padded_image_gpu;
	GpuMultidimArrayAtGpu<float> padded_image2_gpu;
	GpuMultidimArrayAtGpu<float> padded_mask_gpu;
	GpuMultidimArrayAtGpu<float> polar_gpu;
	GpuMultidimArrayAtGpu<float> polar2_gpu;
	GpuMultidimArrayAtGpu< std::complex<float> > MF, MF2;
	GpuMultidimArrayAtGpu< std::complex<float> > RefExpFourier;
	GpuMultidimArrayAtGpu< float > RefExpRealSpace;
	GpuMultidimArrayAtGpu< std::complex<float> > RefExpFourierPolar;
	GpuMultidimArrayAtGpu< float > RefExpRealSpacePolar;
	GpuMultidimArrayAtGpu<float> d_NCC, d_NCCPolar, d_NCCPolar1D, auxMax, auxZero;
	GpuMultidimArrayAtGpu<float> d_out_max, d_pos_max, d_out_polar_max, d_pos_polar_max;
	GpuMultidimArrayAtGpu<float> maxGpu;

};

class GpuCorrelationAux {
public:
	GpuMultidimArrayAtGpu<float> d_original_image;
	GpuMultidimArrayAtGpu<float> d_transform_image;
	GpuMultidimArrayAtGpu< std::complex<float> > d_projFFT;
	GpuMultidimArrayAtGpu< std::complex<float> > d_projSquaredFFT;
	GpuMultidimArrayAtGpu< std::complex<float> > d_projPolarFFT;
	GpuMultidimArrayAtGpu< std::complex<float> > d_projPolarSquaredFFT;
	GpuMultidimArrayAtGpu< std::complex<float> > d_maskFFT;
	GpuMultidimArrayAtGpu<float> maskAutocorrelation;
	int maskCount;
	size_t Xdim, Ydim, XdimOrig, YdimOrig, XdimPolar, YdimPolar;
	//AJ Xdim and Ydim are the padded sizes

	GpuMultidimArrayAtGpu< float > d_denom;
	GpuMultidimArrayAtGpu< float > MFrealSpace, MF2realSpace;

	GpuMultidimArrayAtGpu<float> d_mask;

	void produceSideInfo(mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask, StructuresAux &myStructureAux, myStreamHandle &myStream);
	void produceSideInfo(mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask, StructuresAux &myStructureAux, GpuMultidimArrayAtGpu<float> &maskAutocorr, myStreamHandle &myStream);
	void produceSideInfo(mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask, StructuresAux &myStructureAux, GpuMultidimArrayAtGpu<float> &maskAutocorr, myStreamHandle &myStream, bool skip, mycufftHandle &ifftcb);
};


void cuda_calculate_correlation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux, TransformMatrix<float> &transMat, float *max_vector, int maxShift, mycufftHandle &myhandlePadded, bool mirror, StructuresAux &myStructureAux, myStreamHandle &myStream, TransformMatrix<float> &resultTR, bool saveMaxVector);
void cuda_calculate_correlation_rotation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux, TransformMatrix<float> &transMat, float *max_vector, int maxShift, mycufftHandle &myhandlePadded, bool mirror, StructuresAux &myStructureAux, myStreamHandle &myStream, TransformMatrix<float> &resultRT);
void apply_transform(GpuMultidimArrayAtGpu<float> &d_original_image, GpuMultidimArrayAtGpu<float> &d_transform_image, TransformMatrix<float> &transMat, myStreamHandle &myStream);
void padding_masking(GpuMultidimArrayAtGpu<float> &d_orig_image, GpuMultidimArrayAtGpu<float> &mask, GpuMultidimArrayAtGpu<float> &padded_image_gpu, GpuMultidimArrayAtGpu<float> &padded_image2_gpu, GpuMultidimArrayAtGpu<float> &padded_mask_gpu, bool experimental, myStreamHandle &myStream);
/*void calculateFFTPlanSize(mycufftHandle &myhandle);
void produceSideInfoCuda(GpuCorrelationAux &aux, mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask, mycufftHandle &myhandlePolar, size_t memoryUsed);
void cuda_align_experimental_image(float *expImages, GpuCorrelationAux &d_referenceAux, GpuCorrelationAux &d_experimentalAux,
		TransformMatrix<float> &transMat_tr, TransformMatrix<float> &transMat_rt, float *max_vector_tr, float *max_vector_rt,
		int numImagesRef, GpuMultidimArrayAtGpu<float> &mask, double maxShift,
		mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask, mycufftHandle &myhandlePolar, size_t mdExpSize);
*/
void cuda_calculate_correlation_two(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAuxTR,
		TransformMatrix<float> &transMatTR, float *max_vectorTR, int maxShift,
		mycufftHandle &myhandlePaddedTR, bool mirror, StructuresAux &myStructureAuxTR,
		myStreamHandle &myStreamTR,
		GpuCorrelationAux &experimentalAuxRT, TransformMatrix<float> &transMatRT,
		float *max_vectorRT, mycufftHandle &myhandlePaddedRT,
		StructuresAux &myStructureAuxRT, myStreamHandle &myStreamRT,
		TransformMatrix<float> &resultTR, TransformMatrix<float> &resultRT,
		mycufftHandle &ifftcb, bool saveMaxVector);
#endif
