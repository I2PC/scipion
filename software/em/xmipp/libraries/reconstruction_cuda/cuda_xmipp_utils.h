
#ifndef CUDA_XMIPP_UTILS_H
#define CUDA_XMIPP_UTILS_H

#include <stdio.h>

void gpuMalloc(void** d_data, size_t Nbytes);
void gpuFree(void* d_data);
void gpuCopyFromGPUToGPU(void* d_dataFrom, void* d_dataTo, size_t Nbytes);
void initializeIdentity(float* d_data, size_t Ndim);
void setTranslationMatrix(float* d_data, float posX, float posY, int n);
void setRotationMatrix(float* d_data, float ang, int n);
void gpuCopyFromCPUToGPU(void* data, void* d_data, size_t Nbytes);
void gpuCopyFromGPUToCPU(void* d_data, void* data, size_t Nbytes);
int gridFromBlock(int tasks, int Nthreads);

struct ioTime
{
	//AJ times for CPU processing
	struct timeval t_ini_cpu_mem_in, t_fin_cpu_mem_in, t_ini_cpu_mem_out, t_fin_cpu_mem_out;
	double secs_cpu_mem_in, secs_cpu_mem_out;

	//AJ times for GPU with linear and nearest interpolation
	struct timeval t_ini_gpu_mem_in, t_fin_gpu_mem_in, t_ini_gpu_kernel, t_fin_gpu_kernel, t_ini_gpu_mem_out, t_fin_gpu_mem_out;
	double secs_gpu_mem_in, secs_gpu_kernel, secs_gpu_mem_out;

	//AJ times for GPU with bspline interpolation
	struct timeval t_ini_gpu_mem_inS, t_fin_gpu_mem_inS, t_ini_gpu_filtS, t_fin_gpu_filtS, t_ini_gpu_kernelS, t_fin_gpu_kernelS, t_ini_gpu_mem_outS, t_fin_gpu_mem_outS;
	double secs_gpu_mem_inS, secs_gpu_filtS, secs_gpu_kernelS, secs_gpu_mem_outS;

	//AJ to know how many images are processed
	int size;

	//AJ some flags
	bool calcTime;
	bool spline;
};

inline double timeval_diff(struct timeval *a, struct timeval *b)
{
  return
    (double)(a->tv_sec + (double)a->tv_usec/1000000) -
    (double)(b->tv_sec + (double)b->tv_usec/1000000);
}

extern struct ioTime *mytimes;

void cuda_check_gpu_memory(float* data);

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

#define CONVERT2DIM3(d) (dim3((d).x,(d).y,(d).z))

template<typename T>
class TransformMatrix
{
public:
	size_t Xdim, Ydim, Zdim, Ndim, yxdim, zyxdim, nzyxdim;
    T* d_data;

    TransformMatrix()
    {
		Xdim=Ydim=Zdim=Ndim=yxdim=zyxdim=nzyxdim=0;
		d_data=NULL;
    }

    TransformMatrix(size_t _Ndim, size_t _Xdim=3, size_t _Ydim=3, size_t _Zdim=1)
    {
		Xdim=Ydim=Zdim=Ndim=yxdim=zyxdim=nzyxdim=0;
		d_data=NULL;
		resize(_Ndim, _Xdim, _Ydim, _Zdim);
    }

	template<typename T1>
	void resize(const TransformMatrix<T1>& array)
	{

		resize(array.Ndim, array.Xdim, array.Ydim, array.Zdim);
	}

	void resize(size_t _Ndim, size_t _Xdim=3, size_t _Ydim=3, size_t _Zdim=1)
    {
		if (_Xdim*_Ydim*_Zdim*_Ndim==nzyxdim)
			return;

		clear();

		Xdim=_Xdim;
		Ydim=_Ydim;
		Zdim=_Zdim;
		Ndim=_Ndim;
        yxdim=(size_t)_Ydim*_Xdim;
        zyxdim=yxdim*_Zdim;
        nzyxdim=zyxdim*_Ndim;
        gpuMalloc((void**) &d_data,nzyxdim*sizeof(T));
        initializeIdentity(d_data, Ndim);
    }

	bool isEmpty()
	{
		return d_data==NULL;
	}

	void clear()
	{
		if (d_data!=NULL)
			gpuFree((void*) d_data);
		Xdim=Ydim=Zdim=Ndim=yxdim=zyxdim=nzyxdim=0;
		d_data=NULL;
	}

	~TransformMatrix()
	{
		clear();
	}

	void setTranslation(const float* posX, const float* posY)
	{
		for(int i=0; i<Ndim; i++)
			setTranslationMatrix(d_data, posX[i], posY[i], i);
	}

	void setRotation(const float* ang)
	{
		for(int i=0; i<Ndim; i++)
			setRotationMatrix(d_data, ang[i], i);
	}

	void copyMatrix(TransformMatrix<float> &lastMatrix)
	{
		if (lastMatrix.isEmpty())
			lastMatrix.resize(Ndim);

		gpuCopyFromGPUToGPU(d_data, lastMatrix.d_data, nzyxdim*sizeof(float));
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

	template<typename T1>
	void resize(const GpuMultidimArrayAtGpu<T1>& array)
	{

		resize(array.Xdim, array.Ydim, array.Zdim, array.Ndim);
	}

	void resize(size_t _Xdim, size_t _Ydim=1, size_t _Zdim=1, size_t _Ndim=1)
    {
		if (_Xdim*_Ydim*_Zdim*_Ndim==nzyxdim)
			return;

		clear();

		Xdim=_Xdim;
		Ydim=_Ydim;
		Zdim=_Zdim;
		Ndim=_Ndim;
        yxdim=(size_t)_Ydim*_Xdim;
        zyxdim=yxdim*_Zdim;
        nzyxdim=zyxdim*_Ndim;
        gpuMalloc((void**) &d_data,nzyxdim*sizeof(T));

    }

	bool isEmpty()
	{
		return d_data==NULL;
	}

	void clear()
	{
		if (d_data!=NULL)
			gpuFree((void*) d_data);
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
	void fft(GpuMultidimArrayAtGpu<T1> &fourierTransform);

	// RealSpace must already be resized
	template <typename T1>
	void ifft(GpuMultidimArrayAtGpu<T1> &realSpace);

	void calculateMax(double *max_values, float *posX, float *posY);

};



#endif
