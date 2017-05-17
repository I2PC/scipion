
#include "cuda_utils.h"

void gpuAssert(cudaError_t code, const char *file, int line, bool abort)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

void gpuAssertFFT(cufftResult_t code, const char *file, int line, bool abort)
{
   if (code != CUFFT_SUCCESS)
   {
	   fprintf(stderr, "CUFFT error");
      if (abort) exit(code);
   }
}

cudaPitchedPtr CopyVolumeHostToDevice(const float* host, uint width, uint height, uint depth)
{
	cudaPitchedPtr device = {0};
	const cudaExtent extent = make_cudaExtent(width * sizeof(float), height, depth);
	gpuErrchk(cudaMalloc3D(&device, extent));
	cudaMemcpy3DParms p = {0};
	p.srcPtr = make_cudaPitchedPtr((void*)host, width * sizeof(float), width, height);
	p.dstPtr = device;
	p.extent = extent;
	p.kind = cudaMemcpyHostToDevice;
	gpuErrchk(cudaMemcpy3D(&p));
	return device;
}

void CopyVolumeDeviceToHost(float* host, const cudaPitchedPtr device, uint width, uint height, uint depth)
{
	const cudaExtent extent = make_cudaExtent(width * sizeof(float), height, depth);
	cudaMemcpy3DParms p = {0};
	p.srcPtr = device;
	p.dstPtr = make_cudaPitchedPtr((void*)host, width * sizeof(float), width, height);
	p.extent = extent;
	p.kind = cudaMemcpyDeviceToHost;
	gpuErrchk(cudaMemcpy3D(&p));
	gpuErrchk(cudaFree(device.ptr));  //free the GPU volume
}
