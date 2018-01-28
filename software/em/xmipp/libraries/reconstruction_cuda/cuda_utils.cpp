
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
	   fprintf(stderr, "CUFFT error: %s %s %d\n", _cudaGetErrorEnum(code), file, line);
      if (abort) exit(code);
   }
}
