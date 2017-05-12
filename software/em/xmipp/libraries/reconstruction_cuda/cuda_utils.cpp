
#include "cuda_utils.h"
#include "cuda_check_errors.h"
#include <cuda_runtime.h>


struct ioTime *mytimes;

void gpuMalloc(void** d_data, size_t Nbytes)
{
	gpuErrchk(cudaMalloc(d_data, Nbytes));
}

void gpuFree(void* d_data)
{
	gpuErrchk(cudaFree(d_data));
}

void gpuCopyFromCPUToGPU(void* data, void* d_data, size_t Nbytes)
{
	gpuErrchk(cudaMemcpy(d_data, data, Nbytes, cudaMemcpyHostToDevice));
}

void gpuCopyFromGPUToCPU(void* d_data, void* data, size_t Nbytes)
{
	gpuErrchk(cudaMemcpy(data, d_data, Nbytes, cudaMemcpyDeviceToHost));
}

int gridFromBlock(int tasks, int Nthreads)
{
    int numBlk = tasks/Nthreads;
    if(tasks%Nthreads>0)
    	numBlk++;
    return numBlk;
}

