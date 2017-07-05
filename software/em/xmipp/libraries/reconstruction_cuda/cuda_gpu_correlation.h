
#ifndef CUDA_GPU_CORRELATION_H
#define CUDA_GPU_CORRELATION_H
#include <complex>
#include "cuda_xmipp_utils.h"

void cuda_cart2polar(GpuMultidimArrayAtGpu<double> &image, GpuMultidimArrayAtGpu<double> &polar_image, GpuMultidimArrayAtGpu<double> &polar2_image, bool rotate);

//void cuda_masking(float *image, size_t Xdim, size_t Ydim, size_t Zdim, size_t batch, int radius);
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

	void produceSideInfo();
};

void cuda_calculate_correlation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux, TransformMatrix<float> &transMat, double *max_vector);
void cuda_calculate_correlation_rotation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux, TransformMatrix<float> &transMat, double *max_vector);
void apply_transform(GpuMultidimArrayAtGpu<double> &d_original_image, GpuMultidimArrayAtGpu<double> &d_transform_image, TransformMatrix<float> &transMat);
void padding_masking(GpuMultidimArrayAtGpu<double> &d_orig_image, GpuMultidimArrayAtGpu<double> &mask, GpuMultidimArrayAtGpu<double> &padded_image_gpu, GpuMultidimArrayAtGpu<double> &padded_image2_gpu, GpuMultidimArrayAtGpu<double> &padded_mask_gpu, bool rotation, bool experimental);

#endif
