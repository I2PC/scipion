/***************************************************************************
 *
 * Authors:    Amaya Jimenez      ajimenez@cnb.csic.es (2017)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
#ifndef _PROG_GPU_UTILS
#define _PROG_GPU_UTILS

#include <data/multidim_array.h>
#include <reconstruction_cuda/cuda_gpu_correlation.h>
#include <reconstruction_cuda/cuda_xmipp_utils.h>

template<typename T>
class GpuMultidimArrayAtCpu
{
public:
	int Xdim, Ydim, Zdim, Ndim, yxdim, zyxdim, nzyxdim;
    T* data;

	GpuMultidimArrayAtCpu()
    {
		Xdim=Ydim=Zdim=Ndim=yxdim=zyxdim=nzyxdim=0;
		data=NULL;
    }

	GpuMultidimArrayAtCpu(int _Xdim, int _Ydim=1, int _Zdim=1, int _Ndim=1)
    {
		Xdim=Ydim=Zdim=Ndim=yxdim=zyxdim=nzyxdim=0;
		data=NULL;
		resize(_Xdim, _Ydim, _Zdim, _Ndim);
    }

	void resize(int _Xdim, int _Ydim=1, int _Zdim=1, int _Ndim=1)
    {
		//FIXME what happens if x and y swaps?
		if (_Xdim*_Ydim*_Zdim*_Ndim==nzyxdim){

			return;
		}

		clear();

		Xdim=_Xdim;
		Ydim=_Ydim;
		Zdim=_Zdim;
		Ndim=_Ndim;
        yxdim=(size_t)_Ydim*_Xdim;
        zyxdim=yxdim*_Zdim;
        nzyxdim=zyxdim*_Ndim;
        if (nzyxdim>0){
        	//data=new T[nzyxdim];
        	cpuMalloc((void**)&data, nzyxdim*sizeof(T));
        }else
        	data=NULL;
    }

	void fillImage(int n, const MultidimArray<T> &from)
	{
		//TODO: preguntar esto, en el host puedo usar un memcpy??
		memcpy(data+n*zyxdim, MULTIDIM_ARRAY(from), MULTIDIM_SIZE(from)*sizeof(T));
	}

	bool isEmpty()
	{
		return data==NULL;
	}

	void copyToGpu(GpuMultidimArrayAtGpu<T> &gpuArray, myStreamHandle &myStream)
	{
		if (gpuArray.isEmpty()){
			gpuArray.resize(Xdim,Ydim,Zdim,Ndim);
		}

		gpuCopyFromCPUToGPUStream(data, gpuArray.d_data, nzyxdim*sizeof(T), myStream);
	}

	void copyFromGpu(GpuMultidimArrayAtGpu<T> &gpuArray, myStreamHandle myStream)
	{
		gpuCopyFromGPUToCPUStream(gpuArray.d_data, data, nzyxdim*sizeof(T), myStream);

	}

	void copyToGpuMultiple(GpuMultidimArrayAtGpu<T> &gpuArray, int numCopy, myStreamHandle myStream)
	{
		if (gpuArray.isEmpty())
			gpuArray.resize(Xdim,Ydim,Zdim,numCopy);

		int index = 0;
		for(int i=0; i<numCopy; i++){
			gpuCopyFromCPUToGPUStream(data, &gpuArray.d_data[index], nzyxdim*sizeof(T), myStream);
			index=index+nzyxdim;
		}
	}

	void giveOneImage(GpuMultidimArrayAtCpu<T> out, int idx){
		if (out.isEmpty())
			out.resize(Xdim,Ydim,Zdim,1);
		memcpy(out.data, data[idx*zyxdim], zyxdim*sizeof(T));
	}

	void setOneImage(GpuMultidimArrayAtCpu<T> in, int idx){
		memcpy(data[idx*zyxdim], in.data, zyxdim*sizeof(T));
	}

	void clear()
	{
		if (data!=NULL){;
			cpuFree(data);
		}
		data=NULL;

	}

	~GpuMultidimArrayAtCpu()
	{
		clear();
	}
};

template<class T>
void fillImage(GpuMultidimArrayAtGpu<T> &to, const MultidimArray<T> &from, myStreamHandle myStream, int n=0)
{
	gpuCopyFromCPUToGPUStream(MULTIDIM_ARRAY(from), to.d_data+n*MULTIDIM_SIZE(from), MULTIDIM_SIZE(from)*sizeof(T), myStream);
}

//@}
#endif






