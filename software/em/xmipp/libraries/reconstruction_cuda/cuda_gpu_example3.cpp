
//Host includes
#include "cuda_gpu_example3.h"
#include <iostream>
#include <stdio.h>
//CUDA includes
#include <cuda_runtime.h>


//CUDA functions
__global__ void
vecAdd(const float *A, const float *B, float *C, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        C[i] = A[i] + B[i];
        printf("elem %d, m1 = %f, m2= %f, suma= %f \n", i, A[i], B[i], C[i]);
    }
}

void cuda_funcion(float *m1, float* m2, float *result, int num){

	std::cerr  << "Inside CUDA function " << num << std::endl;


	//CUDA code
	float *d_m1, *d_m3, *d_m2;

	size_t matSize=num*num*sizeof(float);
	cudaMalloc((void **)&d_m1, matSize);
	cudaMalloc((void **)&d_m2, matSize);
	cudaMalloc((void **)&d_m3, matSize);


	cudaMemcpy(d_m1, m1, matSize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_m2, m2, matSize, cudaMemcpyHostToDevice);

	int numTh = 1024;
	int numBlk = num*num/numTh;
	if ((num*num)%numTh >0){
			numBlk++;
	}

	printf("CUDA kernel launch with %d blocks of %d threads\n", numBlk, numTh);
	vecAdd<<<numBlk, numTh>>>(d_m1, d_m2, d_m3, num*num);

	cudaMemcpy(result, d_m3, matSize, cudaMemcpyDeviceToHost);

	cudaFree(d_m1);
	cudaFree(d_m2);
	cudaFree(d_m3);

}
