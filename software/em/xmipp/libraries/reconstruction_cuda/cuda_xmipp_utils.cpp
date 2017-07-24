

#include "cuda_xmipp_utils.h"
#include "cuda_utils.h"

#include <cufft.h>
#include <cuComplex.h>

#include <thrust/extrema.h>
#include <thrust/device_ptr.h>

#include <time.h>
#include <sys/time.h>

#define PI 3.14159265

struct ioTime *mytimes;

void mycufftDestroy(void *ptr)
{
	cufftHandle *planPtr = (cufftHandle *)ptr;
	cufftDestroy(*planPtr);
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

	//cufftHandle planF;
	//cufftCreate(&planF);
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

/*
template<>
template<>
void GpuMultidimArrayAtGpu<double>::fftNew(GpuMultidimArrayAtGpu< std::complex<double> > &fourierTransform, GpuPlanFFT plan)
{

	int Xfdim=(Xdim/2)+1;
	fourierTransform.resize(Xfdim,Ydim,Zdim,Ndim);

	//AJ check the size of the data to avoid exceed the CUDA FFT size
	float memory[3]={0, 0, 0}; //total, free, used
	cuda_check_gpu_memory(memory);

	int positionReal=0;
	int positionFFT=0;
	size_t NdimNew;
	if(Xdim*Ydim*Ndim*64>memory[1]*0.8){
		float sizeAuxIm = Xdim*Ydim*Zdim*sizeof(cufftDoubleReal) + Xfdim*Ydim*Zdim*sizeof(cufftDoubleComplex);
		NdimNew=floor((memory[1]*0.6)/(64*Xdim*Ydim + sizeAuxIm));
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

	while(aux>0){

		GpuMultidimArrayAtGpu<cufftDoubleReal> auxInFFT;
		GpuMultidimArrayAtGpu<cufftDoubleComplex> auxOutFFT;
		if(NdimNew!=Ndim){
			auxInFFT.resize(Xdim,Ydim,Zdim,NdimNew);
			gpuCopyFromGPUToGPU((cufftDoubleReal*)&d_data[positionReal], auxInFFT.d_data, Xdim*Ydim*Zdim*NdimNew*sizeof(cufftDoubleReal));
			auxOutFFT.resize(Xfdim,Ydim,Zdim,NdimNew);
		}


		//AJ TIME
		timeval start, end;
		double secs;
	    gettimeofday(&start, NULL);

    	cufftHandle planF = plan.planF;
	    //gpuErrchkFFT(cufftPlanMany(&planF, NRANK, nr, nr, rstride, rdist, nf, fstride, fdist, CUFFT_D2Z, NdimNew));

	    //AJ TIME
	    gettimeofday(&end, NULL);
	    secs = timeval_diff(&end, &start);
	    printf("cufftPlanMany: %.16g miliseconds\n", secs * 1000.0);

		//AJ TIME
		timeval start2, end2;
		double secs2;
	    gettimeofday(&start2, NULL);


		if(NdimNew!=Ndim){
			gpuErrchkFFT(cufftExecD2Z(planF, auxInFFT.d_data, auxOutFFT.d_data));
		}else{
			gpuErrchkFFT(cufftExecD2Z(planF, (cufftDoubleReal*)&d_data[positionReal], (cufftDoubleComplex*)&fourierTransform.d_data[positionFFT]));
		}

	    //AJ TIME
	    gettimeofday(&end2, NULL);
	    secs2 = timeval_diff(&end2, &start2);
	    printf("cufftExecD2Z: %.16g miliseconds\n", secs2 * 1000.0);

		gpuErrchk(cudaDeviceSynchronize());
		//cufftDestroy(planF);

		if(NdimNew!=Ndim){
			gpuCopyFromGPUToGPU(auxOutFFT.d_data, (cufftDoubleComplex*)&fourierTransform.d_data[positionFFT], Xfdim*Ydim*Zdim*NdimNew*sizeof(cufftDoubleComplex));
			auxOutFFT.clear();
			auxInFFT.clear();
		}

		positionReal+=(NdimNew*Xdim*Ydim*Zdim);
		positionFFT+=(NdimNew*Xfdim*Ydim*Zdim);
		aux-=NdimNew;
		if(aux<NdimNew)
			NdimNew=aux;


	}//AJ end while

}
*/

template<>
template<>
void GpuMultidimArrayAtGpu<float>::fft(GpuMultidimArrayAtGpu< std::complex<float> > &fourierTransform, mycufftHandle &myhandle)
{

	int Xfdim=(Xdim/2)+1;
	fourierTransform.resize(Xfdim,Ydim,Zdim,Ndim);

	//AJ check the size of the data to avoid exceed the CUDA FFT size
	float memory[3]={0, 0, 0}; //total, free, used
	cuda_check_gpu_memory(memory);

	int positionReal=0;
	int positionFFT=0;
	size_t NdimNew, auxNdim;
	if(Xdim*Ydim*Ndim*64>memory[1]*0.8){
		float sizeAuxIm = Xdim*Ydim*Zdim*sizeof(cufftReal) + Xfdim*Ydim*Zdim*sizeof(cufftComplex);
		NdimNew=floor((memory[1]*0.6)/(64*Xdim*Ydim + sizeAuxIm));
	}else
		NdimNew = Ndim;
	size_t aux=Ndim;

	auxNdim=NdimNew;

	while(aux>0){

		GpuMultidimArrayAtGpu<cufftReal> auxInFFT;
		GpuMultidimArrayAtGpu<cufftComplex> auxOutFFT;
		if(NdimNew!=Ndim){
			auxInFFT.resize(Xdim,Ydim,Zdim,NdimNew);
			gpuCopyFromGPUToGPU((cufftReal*)&d_data[positionReal], auxInFFT.d_data, Xdim*Ydim*Zdim*NdimNew*sizeof(cufftReal));
			auxOutFFT.resize(Xfdim,Ydim,Zdim,NdimNew);
		}


		/*/AJ TIME
		timeval start, end;
		double secs;
	    gettimeofday(&start, NULL);*/

		cufftHandle *planFptr = new cufftHandle;
		cufftHandle *planAuxFptr = new cufftHandle;
		if(myhandle.ptr == NULL){
			printf("creo el plan \n");
			createPlanFFT(Xdim, Ydim, NdimNew, Zdim, true, planFptr);
			myhandle.ptr = (void *)planFptr;
			planFptr=(cufftHandle *)myhandle.ptr;
		}else{
			planFptr=(cufftHandle *)myhandle.ptr;
		}
		if(auxNdim!=NdimNew){
			printf("CUIDADO \n");
			createPlanFFT(Xdim, Ydim, NdimNew, Zdim, true, planAuxFptr);
		}


	    /*/AJ TIME
	    gettimeofday(&end, NULL);
	    secs = timeval_diff(&end, &start);
	    printf("cufftPlanMany: %.16g miliseconds\n", secs * 1000.0);*/

		/*/AJ TIME
		timeval start2, end2;
		double secs2;
	    gettimeofday(&start2, NULL);*/

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

	    /*/AJ TIME
	    gettimeofday(&end2, NULL);
	    secs2 = timeval_diff(&end2, &start2);
	    printf("cufftExecD2Z: %.16g miliseconds\n", secs2 * 1000.0);*/

		gpuErrchk(cudaDeviceSynchronize());
		if(auxNdim!=NdimNew)
			cufftDestroy(*planAuxFptr);


		if(NdimNew!=Ndim){
			gpuCopyFromGPUToGPU(auxOutFFT.d_data, (cufftComplex*)&fourierTransform.d_data[positionFFT], Xfdim*Ydim*Zdim*NdimNew*sizeof(cufftComplex));
			auxOutFFT.clear();
			auxInFFT.clear();
		}

		positionReal+=(NdimNew*Xdim*Ydim*Zdim);
		positionFFT+=(NdimNew*Xfdim*Ydim*Zdim);
		aux-=NdimNew;
		if(aux<NdimNew)
			NdimNew=aux;

		auxNdim=NdimNew;

	}//AJ end while

}

template<>
template<>
void GpuMultidimArrayAtGpu< std::complex<float> >::ifft(GpuMultidimArrayAtGpu<float> &realSpace, mycufftHandle &myhandle)
{

	int Xfdim=(realSpace.Xdim/2)+1;

	//AJ check the size of the data to avoid exceed the CUDA FFT size
	float memory[3]={0, 0, 0}; //total, free, used
	cuda_check_gpu_memory(memory);

	int positionReal=0;
	int positionFFT=0;
	size_t NdimNew, auxNdim;
	if(realSpace.Xdim*realSpace.Ydim*realSpace.Ndim*64>memory[1]*0.8){
		float sizeAuxIm = realSpace.Xdim*realSpace.Ydim*realSpace.Zdim*sizeof(cufftReal) + Xfdim*realSpace.Ydim*realSpace.Zdim*sizeof(cufftComplex);
		NdimNew=floor((memory[1]*0.6)/(64*realSpace.Xdim*realSpace.Ydim + sizeAuxIm));
	}else
		NdimNew = realSpace.Ndim;
	size_t aux=realSpace.Ndim;

	auxNdim=NdimNew;

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
		if(myhandle.ptr == NULL){
			printf("creo el plan \n");
			createPlanFFT(realSpace.Xdim, realSpace.Ydim, NdimNew, Zdim, false, planBptr);
			myhandle.ptr = (void *)planBptr;
			planBptr=(cufftHandle *)myhandle.ptr;
		}else{
			planBptr=(cufftHandle *)myhandle.ptr;
		}
		if(auxNdim!=NdimNew){
			printf("CUIDADO \n");
			createPlanFFT(Xdim, Ydim, NdimNew, Zdim, false, planAuxBptr);
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
		if(auxNdim!=NdimNew)
			cufftDestroy(*planAuxBptr);

		if(NdimNew!=Ndim){
			gpuCopyFromGPUToGPU(auxOutFFT.d_data, (cufftReal*)&realSpace.d_data[positionReal], realSpace.Xdim*realSpace.Ydim*realSpace.Zdim*NdimNew*sizeof(cufftReal));
			auxOutFFT.clear();
			auxInFFT.clear();
		}

		positionReal+=(NdimNew*realSpace.Xdim*realSpace.Ydim*realSpace.Zdim);
		positionFFT+=(NdimNew*Xfdim*realSpace.Ydim*realSpace.Zdim);
		aux-=NdimNew;
		if(aux<NdimNew)
			NdimNew=aux;

		auxNdim=NdimNew;

	}//AJ end while

}
/*
void GpuPlanFFT::createGpuPlanFFT()
{
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
    gpuErrchkFFT(cufftPlanMany(&planF, NRANK, nr, nr, rstride, rdist, nf, fstride, fdist, CUFFT_D2Z, Ndim));
}

void GpuPlanFFT::destroyGpuPlanFFT()
{
	cufftDestroy(planF);
}*/

template<>
void GpuMultidimArrayAtGpu<float>::calculateMax(float *max_values, float *posX, float *posY){

	int index = 0;

	for(int i=0; i<Ndim; i++){
		thrust::device_ptr<float> dev_ptr = thrust::device_pointer_cast(&d_data[index]);
		thrust::device_ptr<float> max_ptr = thrust::max_element(dev_ptr, dev_ptr + (int)yxdim);
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


