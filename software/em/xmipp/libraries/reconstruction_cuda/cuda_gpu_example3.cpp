
//Host includes
#include "cuda_gpu_example3.h"
#include <iostream>
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

void cuda_funcion(int num){

	std::cerr  << "Inside CUDA function " << num << std::endl;

	/*MultidimArray<double> m1(num,num);
	MultidimArray<double> m2(num,num);
	MultidimArray<double> mResult(num,num);
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
	m1.initRandom(0, 10, RND_UNIFORM);
	m2.initRandom(100, 200, RND_UNIFORM);

	//CUDA code
	double *d_m1, *d_m3, *d_m2;
	size_t matSize=num*num*sizeof(double);
	cudaMalloc((void **)&d_m1, matSize);
	cudaMalloc((void **)&d_m2, matSize);
	cudaMalloc((void **)&d_m3, matSize);


	cudaMemcpy(d_m1, MULTIDIM_ARRAY(m1), matSize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_m2, MULTIDIM_ARRAY(m2), matSize, cudaMemcpyHostToDevice);

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

	std::cout << "m1" << m1 << std::endl;
	std::cout << "m2" << m2 << std::endl;
	std::cout << "mResult" << mResult << std::endl;*/


}
