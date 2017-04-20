
#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

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

inline double timeval_diff(struct timeval *a, struct timeval *b)
{
  return
    (double)(a->tv_sec + (double)a->tv_usec/1000000) -
    (double)(b->tv_sec + (double)b->tv_usec/1000000);
}

#endif
