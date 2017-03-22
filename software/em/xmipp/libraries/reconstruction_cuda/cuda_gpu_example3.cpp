
//Host includes
#include "cuda_gpu_example3.h"
#include <iostream>
#include <data/multidim_array.h>
//CUDA includes
#include <cuda_runtime.h>


//CUDA functions
__global__ void
vecAdd(const double *A, const double *B, double *C, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
        C[i] = A[i] + B[i];
	//printf("elem %d, m1 = %f, m2= %f, suma= %f \n", i, A[i], B[i], C[i]);
    }
}

void cuda_funcion(){

	std::cerr  << "Inside CUDA function" << std::endl;
	double m3_out_device[4], m3_out_host[4];

	MultidimArray<double> m1(2,2);
	MultidimArray<double> m2(2,2);
	MultidimArray<double> mResult(2,2);
	double * result;
	result = MULTIDIM_ARRAY(mResult);

	A2D_ELEM(m1,0,0) = 1.;
	A2D_ELEM(m1,1,0) = 2.;
	A2D_ELEM(m1,0,1) = 3.;
	A2D_ELEM(m1,1,1) = 4.;

	A2D_ELEM(m2,0,0) = 11.;
	A2D_ELEM(m2,1,0) = 22.;
	A2D_ELEM(m2,0,1) = 1033.;
	A2D_ELEM(m2,1,1) = 44.;

	//CUDA code
	double *d_m1, *d_m3, *d_m2;
	cudaMalloc((void **)&d_m1, 4*sizeof(double));
	cudaMalloc((void **)&d_m2, 4*sizeof(double));
	cudaMalloc((void **)&d_m3, 4*sizeof(double));


	cudaMemcpy(d_m1, MULTIDIM_ARRAY(m1), 4*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_m2, MULTIDIM_ARRAY(m2), 4*sizeof(double), cudaMemcpyHostToDevice);

	int numTh = 1024;
	int numBlk = 1;
	printf("CUDA kernel launch with %d blocks of %d threads\n", numBlk, numTh);
	vecAdd<<<numBlk, numTh>>>(d_m1, d_m2, d_m3, 4);

	cudaMemcpy(result, d_m3, 4*sizeof(double), cudaMemcpyDeviceToHost);

	cudaFree(d_m1);
	cudaFree(d_m2);
	cudaFree(d_m3);

	std::cout << "m1" << m1 << std::endl;
	std::cout << "m2" << m2 << std::endl;
	std::cout << "mResult" << mResult << std::endl;


}
