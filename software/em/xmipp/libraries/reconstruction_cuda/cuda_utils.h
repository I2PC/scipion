
#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

#include <complex>
#include <stdio.h>


std::complex<double>* cuda_fft(double* &image, double* recovered_image, double *module, double *angle, size_t Xdim, size_t Ydim, size_t Zdim, size_t batch);
//double* cuda_ifft(std::complex<double>* &data, size_t Xdim, size_t Ydim, size_t Zdim, size_t batch);

struct ioTime {

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

extern struct ioTime *mytimes;
struct dim3;

inline double timeval_diff(struct timeval *a, struct timeval *b)
{
  return
    (double)(a->tv_sec + (double)a->tv_usec/1000000) -
    (double)(b->tv_sec + (double)b->tv_usec/1000000);
}

void gpuMalloc(void** d_data, size_t Nbytes);
void gpuFree(void* d_data);
void gpuCopyFromCPUToGPU(void* data, void* d_data, size_t Nbytes);
void gpuCopyFromGPUToCPU(void* d_data, void* data, size_t Nbytes);
int gridFromBlock(int tasks, int Nthreads);

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
        resize(_Xdim, _Ydim, _Zdim, _Ndim);
    }

	void resize(size_t _Xdim, size_t _Ydim=1, size_t _Zdim=1, size_t _Ndim=1)
    {
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

	~GpuMultidimArrayAtGpu()
	{
		if (d_data!=NULL)
			gpuFree((void*) d_data);
	}

	void calculateGridSize(const XmippDim3 &blockSize, XmippDim3 &gridSize)
	{
		gridSize.x=gridFromBlock(Xdim,blockSize.x);
		gridSize.y=gridFromBlock(Ydim,blockSize.y);
		gridSize.z=gridFromBlock(Zdim,blockSize.z);
	}

	void calculateGridSizeVectorized(const XmippDim3 &blockSize, XmippDim3 &gridSize)
	{
		gridSize.x=gridFromBlock(nzyxdim,blockSize.x);
		gridSize.y=1;
		gridSize.z=1;
	}

	std::complex<double>* fft()
	{
		return cuda_fft(d_data, NULL, NULL, NULL, Xdim, Ydim, Zdim, Ndim);
	}

	//AJ only for debugging purposes
	std::complex<double>* fft(double* output)
	{
		return cuda_fft(d_data, output, NULL, NULL, Xdim, Ydim, Zdim, Ndim);
	}
	//END AJ

	/*double* ifft()
	{
		return cuda_ifft(d_data, Xdim, Ydim, Zdim, Ndim);
	}*/

};


#endif
