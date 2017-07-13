#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

#include <complex>
#include <stdio.h>
#include <cuda_runtime.h>
#include <cufft.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true);

#define gpuErrchkFFT(ans) { gpuAssertFFT((ans), __FILE__, __LINE__); }
void gpuAssertFFT(cufftResult_t code, const char *file, int line, bool abort=true);

cudaPitchedPtr CopyVolumeHostToDevice(const float* host, uint width, uint height, uint depth);
void CopyVolumeDeviceToHost(float* host, const cudaPitchedPtr device, uint width, uint height, uint depth);

#endif
