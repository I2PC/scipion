
#include "cuda_xmipp_utils.h"
#include "cuda_utils.h"

#include <cufft.h>
#include <cuComplex.h>

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
	int nr1[] = {Xdim};   // --- Size of the image in real space
	int nr2[] = {Ydim, Xdim};   // --- Size of the image in real space
	int nr3[] = {Zdim, Ydim, Xdim};   // --- Size of the image in real space

	int Xfdim=(Xdim/2)+1;
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

    fourierTransform.resize(Xfdim,Ydim,Zdim,Ndim);
    cufftHandle planF;
    gpuErrchkFFT(cufftPlanMany(&planF, NRANK, nr, nr, rstride, rdist, nf, fstride, fdist, CUFFT_D2Z, Ndim));
    gpuErrchkFFT(cufftExecD2Z(planF, (cufftDoubleReal*)d_data, (cufftDoubleComplex*) fourierTransform.d_data));
    gpuErrchk(cudaDeviceSynchronize());
    cufftDestroy(planF);
}

template<>
template<>
void GpuMultidimArrayAtGpu< std::complex<double> >::ifft(GpuMultidimArrayAtGpu<double> &realSpace)
{
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

    cufftHandle planB;
	gpuErrchkFFT(cufftPlanMany(&planB, NRANK, nr, nf, fstride, fdist, nr, rstride, rdist, CUFFT_Z2D, realSpace.Ndim));
	gpuErrchkFFT(cufftExecZ2D(planB, (cufftDoubleComplex *) d_data, (cufftDoubleReal*) realSpace.d_data));
	gpuErrchk(cudaDeviceSynchronize());
	cufftDestroy(planB);
}

