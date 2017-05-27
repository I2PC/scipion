
#ifndef CUDA_GPU_CORRELATION_H
#define CUDA_GPU_CORRELATION_H
#include <complex>
#include "cuda_xmipp_utils.h"

void cuda_cart2polar(GpuMultidimArrayAtGpu<double> &image, GpuMultidimArrayAtGpu<double> &polar_image, GpuMultidimArrayAtGpu<double> &polar2_image, bool rotate);

//void cuda_masking(float *image, size_t Xdim, size_t Ydim, size_t Zdim, size_t batch, int radius);
class GpuCorrelationAux {
public:
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

	GpuMultidimArrayAtGpu< std::complex<double> > d_maskFFTPolar;

	void produceSideInfo();
};

double** cuda_calculate_correlation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux);
double** cuda_calculate_correlation_rotation(GpuCorrelationAux &referenceAux, GpuCorrelationAux &experimentalAux);

#endif
