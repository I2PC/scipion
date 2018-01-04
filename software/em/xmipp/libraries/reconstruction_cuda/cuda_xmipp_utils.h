
#ifndef CUDA_XMIPP_UTILS_H
#define CUDA_XMIPP_UTILS_H

#include <stdio.h>

void mycufftDestroy(void *ptr);
void gpuMalloc(void** d_data, size_t Nbytes);
void gpuFree(void* d_data);
void gpuCopyFromGPUToGPU(void* d_dataFrom, void* d_dataTo, size_t Nbytes);
void gpuCopyFromCPUToGPU(void* data, void* d_data, size_t Nbytes);
void gpuCopyFromGPUToCPU(void* d_data, void* data, size_t Nbytes);
int gridFromBlock(int tasks, int Nthreads);

class mycufftHandle {
public:
	void *ptr;

	mycufftHandle(){
			ptr=NULL;
	}

	void clear()
	{
		if (ptr!=NULL) {
			mycufftDestroy(ptr);
		}
		ptr=NULL;
	}


};

class XmippDim3 {
public:
	int x;
	int y;
	int z;

	XmippDim3(int _x, int _y, int _z)
	{
		x=_x;
		y=_y;
		z=_z;
	}

	XmippDim3()
	{
		x=y=z=0;
	}
};

template<typename T>
class GpuMultidimArrayAtGpu
{
public:
	size_t Xdim, Ydim, Zdim, Ndim, yxdim, zyxdim, nzyxdim;
    T* d_data;

	GpuMultidimArrayAtGpu()
    {
		Xdim=Ydim=Zdim=Ndim=yxdim=zyxdim=nzyxdim=0;
		d_data=NULL;
    }

	GpuMultidimArrayAtGpu(size_t _Xdim, size_t _Ydim=1, size_t _Zdim=1, size_t _Ndim=1)
    {
		Xdim=Ydim=Zdim=Ndim=yxdim=zyxdim=nzyxdim=0;
		d_data=NULL;
		resize(_Xdim, _Ydim, _Zdim, _Ndim);
    }

	GpuMultidimArrayAtGpu(size_t _Xdim, size_t _Ydim, size_t _Zdim, size_t _Ndim, T* deviceData)
    {
		setDims(_Xdim, _Ydim, _Zdim, _Ndim);
		d_data = deviceData;
    }

	template<typename T1>
	void resize(const GpuMultidimArrayAtGpu<T1>& array)
	{

		resize(array.Xdim, array.Ydim, array.Zdim, array.Ndim);
	}

	void resize(size_t _Xdim, size_t _Ydim=1, size_t _Zdim=1, size_t _Ndim=1)
    {

		if (_Xdim*_Ydim*_Zdim*_Ndim==nzyxdim){

			return;
		}

		clear();
		setDims(_Xdim, _Ydim, _Zdim, _Ndim);

        gpuMalloc((void**) &d_data,nzyxdim*sizeof(T));

    }

	bool isEmpty()
	{
		return d_data==NULL;
	}

	void clear()
	{
		if (d_data!=NULL){
			gpuFree((void*) d_data);

		}
		Xdim=Ydim=Zdim=Ndim=yxdim=zyxdim=nzyxdim=0;
		d_data=NULL;
	}

	~GpuMultidimArrayAtGpu()
	{
		clear();
	}

	void copyToGpu(T* data)
	{
		gpuCopyFromCPUToGPU((void *)data, (void *)d_data, nzyxdim*sizeof(T));
	}

	void fillImageToGpu(T* data, size_t n=0)
	{
		gpuCopyFromCPUToGPU((void *)data, (void *)&d_data[n*zyxdim], zyxdim*sizeof(T));
	}

	void copyGpuToGpu(GpuMultidimArrayAtGpu<T> &gpuArray)
	{
		if (gpuArray.isEmpty())
			gpuArray.resize(Xdim,Ydim,Zdim,Ndim);

		gpuCopyFromGPUToGPU(d_data, gpuArray.d_data, nzyxdim*sizeof(T));
	}

	void calculateGridSize(const XmippDim3 &blockSize, XmippDim3 &gridSize) const
	{
		gridSize.x=gridFromBlock(Xdim,blockSize.x);
		gridSize.y=gridFromBlock(Ydim,blockSize.y);
		gridSize.z=gridFromBlock(Zdim,blockSize.z);
	}

	void calculateGridSizeVectorized(const XmippDim3 &blockSize, XmippDim3 &gridSize) const
	{
		gridSize.x=gridFromBlock(nzyxdim,blockSize.x);
		gridSize.y=1;
		gridSize.z=1;
	}

	template <typename T1>
	void fft(GpuMultidimArrayAtGpu<T1> &fourierTransform, mycufftHandle &myhandle);

	// RealSpace must already be resized
	template <typename T1>
	void ifft(GpuMultidimArrayAtGpu<T1> &realSpace, mycufftHandle &myhandle);

	void calculateMax(float *max_values, float *posX, float *posY, int fixPadding);


private:
	void setDims(size_t _Xdim, size_t _Ydim=1, size_t _Zdim=1, size_t _Ndim=1) {
		Xdim=_Xdim;
		Ydim=_Ydim;
		Zdim=_Zdim;
		Ndim=_Ndim;
		yxdim=(size_t)_Ydim*_Xdim;
		zyxdim=yxdim*_Zdim;
		nzyxdim=zyxdim*_Ndim;
	}
};

#endif
