

#ifndef CUDA_MEM_CHK_H
#define CUDA_MEM_CHK_H

//CUDA includes
#include <cuda_runtime.h>
#include "cuda_check_errors.h"

inline __host__ __device__ float3 getAvailableMem()
{
size_t free_byte, total_byte;

gpuErrchk(cudaMemGetInfo( &free_byte, &total_byte ));

double free_db = (double)free_byte;
double total_db = (double)total_byte;
double used_db = total_db - free_db;

float3 memData;
memData.x=total_db;
memData.y=free_db;
memData.z=used_db;

return memData;
}

#endif
