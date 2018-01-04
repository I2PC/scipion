

#include "cuda_xmipp_utils.h"
#include "cuda_utils.h"

#include <cufft.h>
#include <cuComplex.h>

#define PI 3.14159265

void mycufftDestroy(void* ptr)
{
	cufftHandle *planPtr = (cufftHandle *)ptr;
	cufftDestroy(*planPtr);
	delete planPtr;
}

void calculateFFTPlanSize(mycufftHandle &myhandle){
	printf("calculateFFTPlanSize  myhandle.ptr: %p\n",myhandle.ptr);
	size_t ws2;
	cufftHandle *planFptr=(cufftHandle *)myhandle.ptr;
	cufftGetSize(*planFptr, &ws2);
	printf("calculateFFTPlanSize  size %i\n", (int)ws2);
}


void createPlanFFT(size_t Xdim, size_t Ydim, size_t Ndim, size_t Zdim, bool forward, cufftHandle *plan){

	int Xfdim=(Xdim/2)+1;

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

	if(forward){
		gpuErrchkFFT(cufftPlanMany(plan, NRANK, nr, nr, rstride, rdist, nf, fstride, fdist, CUFFT_R2C, Ndim));
	}else{
		gpuErrchkFFT(cufftPlanMany(plan, NRANK, nr, nf, fstride, fdist, nr, rstride, rdist, CUFFT_C2R, Ndim));
	}

}


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

	float free_db = (float)free_byte;
	float total_db = (float)total_byte;
	float used_db = total_db - free_db;

	data[0]=total_db;
	data[1]=free_db;
	data[2]=used_db;
}

void cuda_check_gpu_properties(int* grid)
{
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	grid[0] = prop.maxGridSize[0];
	grid[1] = prop.maxGridSize[1];
	grid[2] = prop.maxGridSize[2];
}

template<>
template<>
void GpuMultidimArrayAtGpu<float>::fft(GpuMultidimArrayAtGpu< std::complex<float> > &fourierTransform, mycufftHandle &myhandle)
{

	int Xfdim=(Xdim/2)+1;
	fourierTransform.resize(Xfdim,Ydim,Zdim,Ndim);

	int positionReal=0;
	int positionFFT=0;
	size_t NdimNew, auxNdim;
	NdimNew = Ndim;
	size_t aux=Ndim;

	auxNdim=Ndim;

	if(myhandle.ptr!=NULL) {
		NdimNew = Ndim;
	}

	while(aux>0){

		GpuMultidimArrayAtGpu<cufftReal> auxInFFT;
		GpuMultidimArrayAtGpu<cufftComplex> auxOutFFT;
		if(NdimNew!=Ndim){
			auxInFFT.resize(Xdim,Ydim,Zdim,NdimNew);
			gpuCopyFromGPUToGPU((cufftReal*)&d_data[positionReal], auxInFFT.d_data, Xdim*Ydim*Zdim*NdimNew*sizeof(cufftReal));
			auxOutFFT.resize(Xfdim,Ydim,Zdim,NdimNew);
		}

		cufftHandle *planFptr = NULL;
		cufftHandle *planAuxFptr = NULL;
		if(auxNdim!=NdimNew){
			planAuxFptr = new cufftHandle;
			createPlanFFT(Xdim, Ydim, NdimNew, Zdim, true, planAuxFptr);
		}else{
			if(myhandle.ptr == NULL){
				myhandle.ptr = planFptr = new cufftHandle;
				createPlanFFT(Xdim, Ydim, NdimNew, Zdim, true, planFptr);
			}
			planFptr=(cufftHandle *)myhandle.ptr;
		}

		if(auxNdim==NdimNew){
			if(NdimNew!=Ndim){
				gpuErrchkFFT(cufftExecR2C(*planFptr, auxInFFT.d_data, auxOutFFT.d_data));
			}else{
				gpuErrchkFFT(cufftExecR2C(*planFptr, (cufftReal*)&d_data[positionReal], (cufftComplex*)&fourierTransform.d_data[positionFFT]));
			}
		}else{
			if(NdimNew!=Ndim){
				gpuErrchkFFT(cufftExecR2C(*planAuxFptr, auxInFFT.d_data, auxOutFFT.d_data));
			}else{
				gpuErrchkFFT(cufftExecR2C(*planAuxFptr, (cufftReal*)&d_data[positionReal], (cufftComplex*)&fourierTransform.d_data[positionFFT]));
			}
		}

		gpuErrchk(cudaDeviceSynchronize());

		if(NdimNew!=Ndim){
			gpuCopyFromGPUToGPU(auxOutFFT.d_data, (cufftComplex*)&fourierTransform.d_data[positionFFT], Xfdim*Ydim*Zdim*NdimNew*sizeof(cufftComplex));
			auxOutFFT.clear();
			auxInFFT.clear();
		}

		auxNdim=NdimNew;

		positionReal+=(NdimNew*Xdim*Ydim*Zdim);
		positionFFT+=(NdimNew*Xfdim*Ydim*Zdim);
		aux-=NdimNew;
		if(aux<NdimNew)
			NdimNew=aux;

		if (NULL != planAuxFptr) {
			cufftDestroy(*planAuxFptr); // destroy if created
		}

	}//AJ end while

}

template<>
template<>
void GpuMultidimArrayAtGpu< std::complex<float> >::ifft(GpuMultidimArrayAtGpu<float> &realSpace, mycufftHandle &myhandle)
{

	int Xfdim=(realSpace.Xdim/2)+1;

	int positionReal=0;
	int positionFFT=0;
	size_t NdimNew, auxNdim;
	NdimNew = realSpace.Ndim;
	size_t aux=realSpace.Ndim;

	auxNdim=realSpace.Ndim;

	if(myhandle.ptr!=NULL)
		NdimNew = Ndim;

	while(aux>0){

		GpuMultidimArrayAtGpu<cufftComplex> auxInFFT;
		GpuMultidimArrayAtGpu<cufftReal> auxOutFFT;
		if(NdimNew!=Ndim){
			auxInFFT.resize(Xfdim,realSpace.Ydim,realSpace.Zdim,NdimNew);
			gpuCopyFromGPUToGPU((cufftComplex*)&d_data[positionFFT], auxInFFT.d_data, Xfdim*realSpace.Ydim*realSpace.Zdim*NdimNew*sizeof(cufftComplex));
			auxOutFFT.resize(realSpace.Xdim,realSpace.Ydim,realSpace.Zdim, NdimNew);
		}

		cufftHandle *planBptr = new cufftHandle;
		cufftHandle *planAuxBptr = new cufftHandle;
		if(auxNdim!=NdimNew){
			createPlanFFT(Xdim, Ydim, NdimNew, Zdim, false, planAuxBptr);
		}else{
			if(myhandle.ptr == NULL){
				createPlanFFT(realSpace.Xdim, realSpace.Ydim, NdimNew, Zdim, false, planBptr);
				myhandle.ptr = (void *)planBptr;
				planBptr=(cufftHandle *)myhandle.ptr;
			}else{
				planBptr=(cufftHandle *)myhandle.ptr;
			}
		}

		if(auxNdim==NdimNew){
			if(NdimNew!=Ndim){
				gpuErrchkFFT(cufftExecC2R(*planBptr, auxInFFT.d_data, auxOutFFT.d_data));
			}else{
				gpuErrchkFFT(cufftExecC2R(*planBptr, (cufftComplex *)&d_data[positionFFT], (cufftReal*)&realSpace.d_data[positionReal]));
			}
		}else{
			if(NdimNew!=Ndim){
				gpuErrchkFFT(cufftExecC2R(*planAuxBptr, auxInFFT.d_data, auxOutFFT.d_data));
			}else{
				gpuErrchkFFT(cufftExecC2R(*planAuxBptr, (cufftComplex *)&d_data[positionFFT], (cufftReal*)&realSpace.d_data[positionReal]));
			}
		}


		gpuErrchk(cudaDeviceSynchronize());

		if(NdimNew!=Ndim){
			gpuCopyFromGPUToGPU(auxOutFFT.d_data, (cufftReal*)&realSpace.d_data[positionReal], realSpace.Xdim*realSpace.Ydim*realSpace.Zdim*NdimNew*sizeof(cufftReal));
			auxOutFFT.clear();
			auxInFFT.clear();
		}

		auxNdim=NdimNew;

		positionReal+=(NdimNew*realSpace.Xdim*realSpace.Ydim*realSpace.Zdim);
		positionFFT+=(NdimNew*Xfdim*realSpace.Ydim*realSpace.Zdim);
		aux-=NdimNew;
		if(aux<NdimNew)
			NdimNew=aux;

		if(auxNdim!=NdimNew && NdimNew!=0)
			cufftDestroy(*planAuxBptr);

	}//AJ end while

}
