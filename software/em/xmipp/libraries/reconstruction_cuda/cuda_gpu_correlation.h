
#include <time.h>
#include <sys/time.h>
#include <complex>
#include "cuda_utils.h"


void cuda_check_gpu_memory(float* data);

std::complex<double>* cuda_fft(double* &image, double* recovered_image, double *module, double *angle, size_t Xdim, size_t Ydim, size_t Zdim, size_t batch);

void cuda_cart2polar(GpuMultidimArrayAtGpu<double> &image, GpuMultidimArrayAtGpu<double> &polar_image, GpuMultidimArrayAtGpu<double> &polar2_image, bool rotate);

//void cuda_masking(float *image, size_t Xdim, size_t Ydim, size_t Zdim, size_t batch, int radius);

double** cuda_calculate_correlation(std::complex<double>* d_projFFTPointer, std::complex<double>* d_projSquaredFFTPointer, std::complex<double>* d_expFFTPointer, std::complex<double>* d_expSquaredFFTPointer, std::complex<double>* d_maskFFTPointer, size_t Xdim, size_t Ydim, size_t Zdim, size_t numProj, size_t numExp, int counting);


