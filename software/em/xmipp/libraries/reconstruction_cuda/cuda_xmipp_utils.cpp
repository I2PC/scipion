
#include "cuda_xmipp_utils.h"
#include "cuda_utils.h"

#include <cufft.h>
#include <cuComplex.h>

#include <thrust/extrema.h>
#include <thrust/device_ptr.h>

#define PI 3.14159265

struct ioTime *mytimes;

void gpuMalloc(void** d_data, size_t Nbytes)
{
	gpuErrchk(cudaMalloc(d_data, Nbytes));
}

void gpuFree(void* d_data)
{
	gpuErrchk(cudaFree(d_data));
}


void initializeIdentity(float* d_data, size_t Ndim)
{
	float identity[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
	for(int i=0; i<Ndim; i++)
		gpuErrchk(cudaMemcpy((void*)&d_data[i*9], &identity[0], 9*sizeof(float), cudaMemcpyHostToDevice));
}

void setTranslationMatrix(float* d_data, float posX, float posY, int n)
{
	float matrix[9] = {1, 0, posX, 0, 1, posY, 0, 0, 1};
	gpuErrchk(cudaMemcpy((void*)&d_data[n*9], &matrix[0], 9*sizeof(float), cudaMemcpyHostToDevice));
}

void setRotationMatrix(float* d_data, float ang, int n)
{
	float rad = (float)(ang*PI/180);
	float matrix[9] = {cosf(rad), -sinf(rad), 0, sinf(rad), cosf(rad), 0, 0, 0, 1};
	gpuErrchk(cudaMemcpy((void*)&d_data[n*9], &matrix[0], 9*sizeof(float), cudaMemcpyHostToDevice));
}

void gpuCopyFromCPUToGPU(void* data, void* d_data, size_t Nbytes)
{
	gpuErrchk(cudaMemcpy(d_data, data, Nbytes, cudaMemcpyHostToDevice));
}

void gpuCopyFromGPUToCPU(void* d_data, void* data, size_t Nbytes)
{
	gpuErrchk(cudaMemcpy(data, d_data, Nbytes, cudaMemcpyDeviceToHost));
}

void gpuCopyFromGPUToGPU(void* d_dataFrom, void* d_dataTo, size_t Nbytes)
{
	gpuErrchk(cudaMemcpy(d_dataTo, d_dataFrom, Nbytes, cudaMemcpyDeviceToDevice));
}

int gridFromBlock(int tasks, int Nthreads)
{
    int numBlk = tasks/Nthreads;
    if(tasks%Nthreads>0)
    	numBlk++;
    return numBlk;
}
void cuda_check_gpu_memory(float* data)
{
	size_t free_byte, total_byte;
	gpuErrchk(cudaMemGetInfo( &free_byte, &total_byte ));

	double free_db = (double)free_byte;
	double total_db = (double)total_byte;
	double used_db = total_db - free_db;

	data[0]=total_db;
	data[1]=free_db;
	data[2]=used_db;
}


template<>
template<>
void GpuMultidimArrayAtGpu<double>::fft(GpuMultidimArrayAtGpu< std::complex<double> > &fourierTransform)
{

	//printf("EN FFT\n");

	int Xfdim=(Xdim/2)+1;
	fourierTransform.resize(Xfdim,Ydim,Zdim,Ndim);

	//TODO AJ check the size of the data to avoid exceed the CUDA FFT size
	//TODO depende de como de llena este la GPU
	int maxPoints=51000000;
	int positionReal=0;
	int positionFFT=0;
	size_t NdimNew;
	if(Xdim*Ydim*Ndim>maxPoints){
		NdimNew=floor(maxPoints/(Xdim*Ydim));
	}else
		NdimNew = Ndim;
	size_t aux=Ndim;

	int nr1[] = {Xdim};   // --- Size of the image in real space
	int nr2[] = {Ydim, Xdim};   // --- Size of the image in real space
	int nr3[] = {Zdim, Ydim, Xdim};   // --- Size of the image in real space

	int nf1[] = {Xfdim};   // --- Size of the Fourier transform
	int nf2[] = {Ydim, Xfdim};   // --- Size of the Fourier transform
	int nf3[] = {Zdim, Ydim, Xfdim};   // --- Size of the Fourier transform
	int *nr=NULL, *nf=NULL;
	int NRANK; // 1D, 2D or 3D FFTs
	if (Ydim==1 && Zdim==1)
	{
		NRANK=1;
		nr=nr1;
		nf=nf1;
	}
	else if (Zdim==1)
	{
		NRANK=2;
		nr=nr2;
		nf=nf2;
	}
	else
	{
		NRANK=3;
		nr=nr3;
		nf=nf3;
	}

	int rstride = 1;				// --- Distance between two successive input/output elements
	int fstride = 1;
	int rdist = Xdim*Ydim*Zdim;	    // --- Distance between batches
	int fdist = Xfdim*Ydim*Zdim;

	while(aux){

		//printf("NdimNew EN FFT %i\n",(int)NdimNew);
		//printf("Number of points %i \n",(int)NdimNew*Xdim*Ydim);
		//printf("%i %i %i %i %i %i %i \n",NRANK, *nr, rstride, rdist, *nf, fstride, fdist);
		//printf("positionReal EN FFT %i\n",(int)positionReal);
		//printf("positionFFT EN FFT %i\n",(int)positionFFT);

		cufftHandle planF;
		gpuErrchkFFT(cufftPlanMany(&planF, NRANK, nr, nr, rstride, rdist, nf, fstride, fdist, CUFFT_D2Z, NdimNew));
		gpuErrchkFFT(cufftExecD2Z(planF, (cufftDoubleReal*)&d_data[positionReal], (cufftDoubleComplex*)&fourierTransform.d_data[positionFFT]));
		gpuErrchk(cudaDeviceSynchronize());
		cufftDestroy(planF);

		positionReal+=(NdimNew*Xdim*Ydim*Zdim);
		positionFFT+=(NdimNew*Xfdim*Ydim*Zdim);
		aux-=NdimNew;
		if(aux*Xdim*Ydim <= maxPoints){
			NdimNew=aux;
		}
	}//AJ end while
}

template<>
template<>
void GpuMultidimArrayAtGpu< std::complex<double> >::ifft(GpuMultidimArrayAtGpu<double> &realSpace)
{

	//printf("EN IFFT\n");

	//TODO AJ check the size of the data to avoid exceed the CUDA FFT size
	int maxPoints=51000000;
	int positionReal=0;
	int positionFFT=0;
	size_t NdimNew;
	if(realSpace.Xdim*realSpace.Ydim*realSpace.Ndim>maxPoints){
		NdimNew=floor(maxPoints/(realSpace.Xdim*realSpace.Ydim));
	}else
		NdimNew = realSpace.Ndim;
	size_t aux=realSpace.Ndim;


	int nr1[] = {realSpace.Xdim};   // --- Size of the image in real space
	int nr2[] = {realSpace.Ydim, realSpace.Xdim};   // --- Size of the image in real space
	int nr3[] = {realSpace.Zdim, realSpace.Ydim, realSpace.Xdim};   // --- Size of the image in real space

	int Xfdim=(realSpace.Xdim/2)+1;
	int nf1[] = {Xfdim};   // --- Size of the Fourier transform
	int nf2[] = {realSpace.Ydim, Xfdim};   // --- Size of the Fourier transform
	int nf3[] = {realSpace.Zdim, realSpace.Ydim, Xfdim};   // --- Size of the Fourier transform
	int *nr=NULL, *nf=NULL;
	int NRANK; // 1D, 2D or 3D FFTs
	if (realSpace.Ydim==1 && realSpace.Zdim==1)
	{
		NRANK=1;
		nr=nr1;
		nf=nf1;
	}
	else if (realSpace.Zdim==1)
	{
		NRANK=2;
		nr=nr2;
		nf=nf2;
	}
	else
	{
		NRANK=3;
		nr=nr3;
		nf=nf3;
	}

    int rstride = 1;				// --- Distance between two successive input/output elements
	int fstride = 1;
    int rdist = realSpace.Xdim*realSpace.Ydim*realSpace.Zdim;	    // --- Distance between batches
	int fdist = Xfdim*realSpace.Ydim*realSpace.Zdim;

	while(aux){

		//printf("NdimNew EN IFFT %i\n",(int)NdimNew);
		//printf("Number of points %i \n",(int)NdimNew*realSpace.Xdim*realSpace.Ydim);
		cufftHandle planB;
		gpuErrchkFFT(cufftPlanMany(&planB, NRANK, nr, nf, fstride, fdist, nr, rstride, rdist, CUFFT_Z2D, NdimNew));
		gpuErrchkFFT(cufftExecZ2D(planB, (cufftDoubleComplex *)&d_data[positionFFT], (cufftDoubleReal*)&realSpace.d_data[positionReal]));
		gpuErrchk(cudaDeviceSynchronize());
		cufftDestroy(planB);

		positionReal+=(NdimNew*realSpace.Xdim*realSpace.Ydim*realSpace.Zdim);
		positionFFT+=(NdimNew*Xfdim*realSpace.Ydim*realSpace.Zdim);
		aux-=NdimNew;
		if(aux*realSpace.Xdim*realSpace.Ydim <= maxPoints){
			NdimNew=aux;
		}

	}//AJ end while
}


template<>
void GpuMultidimArrayAtGpu<double>::calculateMax(double *max_values, float *posX, float *posY){

	int index = 0;

	for(int i=0; i<Ndim; i++){
		thrust::device_ptr<double> dev_ptr = thrust::device_pointer_cast(&d_data[index]);
		thrust::device_ptr<double> max_ptr = thrust::max_element(dev_ptr, dev_ptr + (int)yxdim);
		unsigned int position = &max_ptr[0] - &dev_ptr[0];
		max_values[i] = max_ptr[0];
		//printf("max_thrust %lf   ", max_values[i]);

		float posX_aux = (float)(position%Xdim);
		float posY_aux = (float)(position/Xdim);
		float Xdim2 = (float)(Xdim/2);
		float Ydim2 = (float)(Ydim/2);

		if(posX_aux>=Xdim2 && posY_aux>=Ydim2){
			posX[i] = Xdim-1-posX_aux;
			posY[i] = Ydim-1-posY_aux;
		}else if(posX_aux<Xdim2 && posY_aux>=Ydim2){
			posX[i] = -(posX_aux+1);
			posY[i] = Ydim-1-posY_aux;
		}else if(posX_aux<Xdim2 && posY_aux<Ydim2){
			posX[i] = -(posX_aux+1);
			posY[i] = -(posY_aux+1);
		}else if(posX_aux>=Xdim2 && posY_aux<Ydim2){
			posX[i] = Xdim-1-posX_aux;
			posY[i] = -(posY_aux+1);
		}
		//printf("Max x-thr=%f   ", posX[i]);
		//printf("Max y-thr=%f\n", posY[i]);

		index = index+(int)yxdim;
	}

}


