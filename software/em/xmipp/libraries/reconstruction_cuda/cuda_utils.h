#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

#include <complex>
#include <stdio.h>
#include <cuda_runtime.h>
#include <cufft.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true);

// cuFFT API errors
static const char *_cudaGetErrorEnum(cufftResult error)
{
    switch (error)
    {
        case CUFFT_SUCCESS:
            return "CUFFT_SUCCESS";

        case CUFFT_INVALID_PLAN:
            return "CUFFT_INVALID_PLAN";

        case CUFFT_ALLOC_FAILED:
            return "CUFFT_ALLOC_FAILED";

        case CUFFT_INVALID_TYPE:
            return "CUFFT_INVALID_TYPE";

        case CUFFT_INVALID_VALUE:
            return "CUFFT_INVALID_VALUE";

        case CUFFT_INTERNAL_ERROR:
            return "CUFFT_INTERNAL_ERROR";

        case CUFFT_EXEC_FAILED:
            return "CUFFT_EXEC_FAILED";

        case CUFFT_SETUP_FAILED:
            return "CUFFT_SETUP_FAILED";

        case CUFFT_INVALID_SIZE:
            return "CUFFT_INVALID_SIZE";

        case CUFFT_UNALIGNED_DATA:
            return "CUFFT_UNALIGNED_DATA";
    }

    return "<unknown>";
}

#define gpuErrchkFFT(ans) { gpuAssertFFT((ans), __FILE__, __LINE__); }
void gpuAssertFFT(cufftResult_t code, const char *file, int line, bool abort=true);

template<typename T>
void hostRegister(T* srcArray, size_t bytes, unsigned int flags=0) {
	if (NULL != srcArray && (0 != bytes)) {
		cudaHostRegister(srcArray, bytes, flags);
		gpuErrchk( cudaPeekAtLastError() );
	}
}

template<typename T>
void hostUnregister(T* srcArray) {
	if (NULL != srcArray) {
		cudaHostUnregister(srcArray);
		gpuErrchk( cudaPeekAtLastError() );
	}
}

cudaPitchedPtr CopyVolumeHostToDevice(const float* host, uint width, uint height, uint depth);
void CopyVolumeDeviceToHost(float* host, const cudaPitchedPtr device, uint width, uint height, uint depth);

#endif
