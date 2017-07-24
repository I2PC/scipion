
#ifndef CUDA_GPU_CORRELATION_H
#define CUDA_GPU_CORRELATION_H
#include <complex>
#include "cuda_xmipp_utils.h"

void cuda_cart2polar(GpuMultidimArrayAtGpu<double> &image, GpuMultidimArrayAtGpu<double> &polar_image, GpuMultidimArrayAtGpu<double> &polar2_image, bool rotate);

class GpuCorrelationAux {
public:
	GpuMultidimArrayAtGpu<double> d_original_image;
	GpuMultidimArrayAtGpu<double> d_transform_image;
	GpuMultidimArrayAtGpu< std::complex<double> > d_projFFT;
	GpuMultidimArrayAtGpu< std::complex<double> > d_projSquaredFFT;
	GpuMultidimArrayAtGpu< std::complex<double> > d_projPolarFFT;
	GpuMultidimArrayAtGpu< std::complex<double> > d_projPolarSquaredFFT;
	GpuMultidimArrayAtGpu< std::complex<double> > d_maskFFT;
	GpuMultidimArrayAtGpu<double> maskAutocorrelation;
	int maskCount;
	int Xdim, Ydim, XdimPolar, YdimPolar;

	GpuMultidimArrayAtGpu< double > d_denom;
	GpuMultidimArrayAtGpu< double > MFrealSpace;
	GpuMultidimArrayAtGpu<double> debug;
	//GpuMultidimArrayAtGpu< std::complex<double> > debug_complex;

	//AJ new
	GpuMultidimArrayAtGpu<double> d_padded_image;
	GpuMultidimArrayAtGpu<double> d_padded2_image;
	GpuMultidimArrayAtGpu<double> d_mask;
	GpuMultidimArrayAtGpu<double> d_polar_image;
	GpuMultidimArrayAtGpu<double> d_polar2_image;

	void produceSideInfo(mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask);
};


void cuda_calculate_correlation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux, TransformMatrix<float> &transMat, double *max_vector, double maxShift, mycufftHandle &myhandlePadded);
void cuda_calculate_correlation_rotation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux, TransformMatrix<float> &transMat, double *max_vector, double maxShift, mycufftHandle &myhandlePadded);
void apply_transform(GpuMultidimArrayAtGpu<double> &d_original_image, GpuMultidimArrayAtGpu<double> &d_transform_image, TransformMatrix<float> &transMat);
void padding_masking(GpuMultidimArrayAtGpu<double> &d_orig_image, GpuMultidimArrayAtGpu<double> &mask, GpuMultidimArrayAtGpu<double> &padded_image_gpu, GpuMultidimArrayAtGpu<double> &padded_image2_gpu, GpuMultidimArrayAtGpu<double> &padded_mask_gpu, bool rotation, bool experimental);
void fft_v2(size_t Xdim, size_t Ydim, size_t Ndim, size_t Zdim, double *in, std::complex<double> *out, mycufftHandle &myhandle);
void ifft_v2(size_t Xdim, size_t Ydim, size_t Ndim, size_t Zdim, std::complex<double> *in, double *out, mycufftHandle &myhandle);
//void createPlanFFT(size_t Xdim, size_t Ydim, size_t Ndim, size_t Zdim, mycufftHandle &myhandle, size_t &workSize);
void calculateFFTPlanSize(mycufftHandle &myhandle);
void produceSideInfoCuda(GpuCorrelationAux &aux, mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask, mycufftHandle &myhandlePolar, size_t memoryUsed);
void cuda_align_experimental_image(double *expImages, GpuCorrelationAux &d_referenceAux, GpuCorrelationAux &d_experimentalAux,
		TransformMatrix<float> &transMat_tr, TransformMatrix<float> &transMat_rt, double *max_vector_tr, double *max_vector_rt,
		int numImagesRef, GpuMultidimArrayAtGpu<double> &mask, double maxShift,
		mycufftHandle &myhandlePadded, mycufftHandle &myhandleMask, mycufftHandle &myhandlePolar, size_t mdExpSize);

#endif
