

//Host includes
#include <data/multidim_array.h>
#include <data/matrix2d.h>
//#include <data/xmipp_fftw.h>

//CUDA includes
#include <cuda_runtime.h>

/*__global__ void
vecAdd(const float *A, MultidimArray<float> *B, float *C, int numElements)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numElements)
    {
	printf("Soy %i\n",i);
        C[i] = A[i] + (*B).data[i];
	printf("elem %d, m1 = %f, m2= %f, suma= %f \n", i, A[i], (*B).data[i], C[i]);
    }
}*/

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


int main(void)
{

    //Host code
/*
    float m1[4];
	m1[0]=1.;
	m1[1]=2.;
	m1[2]=3.;
	m1[3]=4.;

    float m2[4];
	m2[0]=11.;
	m2[1]=22.;
	m2[2]=33.;
	m2[3]=44.;
*/
    double m3_out_device[4], m3_out_host[4];

    MultidimArray<double> m1(2,2);
    MultidimArray<double> m2(2,2);
    MultidimArray<double> mResult(2,2);
    double * result;
    result = MULTIDIM_ARRAY(mResult);
    //MultidimArray<float> m3_out_host(4);
    //MultidimArray<float> m3_out_device(4);

    A2D_ELEM(m1,0,0) = 1.;
    A2D_ELEM(m1,1,0) = 2.;
    A2D_ELEM(m1,0,1) = 3.;
    A2D_ELEM(m1,1,1) = 4.;

    A2D_ELEM(m2,0,0) = 11.;
    A2D_ELEM(m2,1,0) = 22.;
    A2D_ELEM(m2,0,1) = 333.;
    A2D_ELEM(m2,1,1) = 44.;


    //m3_out_host = m1 + m2;
/*
    for(int i=0; i< 4;i++){
           m3_out_host[i] = m1[i]+m2[i];	
     }
*/



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

    
    //Print values in host and device to check
    /*
    for(int i=0; i< 4;i++){
           printf("m1[%i] = %f\n", i, m1[i]);	
    }
    for(int i=0; i< 4;i++){
           printf("m2[%i] = %f\n", i, m2[i]);	
    }

   
    for(int i=0; i< 4;i++){
           printf("m3_out_host[%i] = %f\n", i, m3_out_host[i]);	
    }
*/
    std::cout << "m1" << m1 << std::endl;
    std::cout << "m2" << m2 << std::endl;
    std::cout << "mResult" << mResult << std::endl;

/*
    for(int i=0; i< 4;i++){
           printf("m3_out_device[%i] = %f\n", i, m3_out_device[i]);	
    }

*/
    /*MultidimArray< std::complex< double > > FFT1, coeffs;
    FourierTransformer transformer1;

    //void FourierTransform(T& v, T1& V, bool getCopy=true)
    transformer1.FourierTransform(m2, FFT1, true);
    //void getFourierCopy(T& V)

    std::cerr << "FFT1" << FFT1 << std::endl;

    //void inverseFourierTransform(T& V, T1& v)
    MultidimArray<double> recover(2,2);
    transformer1.inverseFourierTransform(FFT1, recover);

    std::cerr << "FFT1" << FFT1 << std::endl;
    //std::cerr << "What????" << std::endl;
    std::cerr << "recover" << recover << std::endl;

    transformer1.cleanup();*/



    return 0;
}

