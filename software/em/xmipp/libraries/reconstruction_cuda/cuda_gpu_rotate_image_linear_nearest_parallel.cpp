
//Host includes
#include "cuda_gpu_rotate_image_linear_nearest_parallel.h"

#include <iostream>
#include <stdio.h>
//CUDA includes
#include "cuda_utils.h"
#include "cuda_interpolation2D_basic_rotation_parallel.h"

#define NUM_TEX 5
cudaTextureObject_t texImage[NUM_TEX];
cudaStream_t stream[NUM_TEX] ;

void freeRotatedImage (float **rotated_image){

	for(int i=0; i< NUM_TEX; i++){
		cudaFreeHost(rotated_image[i]);
	}
	cudaFreeHost(rotated_image);
}


void CreateTexture(float **image, size_t *Xdim, size_t *Ydim, size_t *Zdim, int interp, int wrap, int i)
{

	//for (int i = 0; i < NUM_TEX; i++){

		std::cerr << "Create texture: " << i << std::endl;

		float *d_inputIm;//Device Array
		size_t pitch;
		gpuErrchk(cudaMallocPitch(&d_inputIm, &pitch, sizeof(float)*Xdim[i], Ydim[i]));
		gpuErrchk(cudaMemcpy2D(d_inputIm, pitch, image[i], sizeof(float)*Xdim[i], sizeof(float)*Xdim[i], Ydim[i], cudaMemcpyHostToDevice));

		cudaResourceDesc resDesc;
		memset(&resDesc, 0, sizeof(resDesc));
		resDesc.resType = cudaResourceTypePitch2D;
		resDesc.res.pitch2D.devPtr = d_inputIm;
		resDesc.res.pitch2D.pitchInBytes =  pitch;
		resDesc.res.pitch2D.width = Xdim[i];
		resDesc.res.pitch2D.height = Ydim[i];
		resDesc.res.pitch2D.desc = cudaCreateChannelDesc<float>(); // bits per channel

		cudaTextureDesc texDesc;
		memset(&texDesc, 0, sizeof(texDesc));
		texDesc.readMode = cudaReadModeElementType;
		texDesc.normalizedCoords = true;
		if (interp==0){
			texDesc.filterMode = cudaFilterModePoint;
		}else{
		    texDesc.filterMode = cudaFilterModeLinear;
		}
		texDesc.addressMode[0] = (cudaTextureAddressMode)wrap;
		texDesc.addressMode[1] = (cudaTextureAddressMode)wrap;
		gpuErrchk(cudaCreateTextureObject(&texImage[i], &resDesc, &texDesc, NULL));


	    /*float *d_inputIm;//Device Array
        gpuErrchk(cudaMalloc((void **)&d_inputIm, Xdim[i]*Ydim[i]*sizeof(float)));
        gpuErrchk(cudaMemcpy(d_inputIm, &image[i], Xdim[i]*Ydim[i]*sizeof(float), cudaMemcpyHostToDevice));

        //cudaArray Descriptor
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
        //cuda Array
        cudaArray *d_cuArr;
        gpuErrchk( cudaMalloc3DArray(&d_cuArr, &channelDesc, make_cudaExtent(Xdim[i]*sizeof(float), Ydim[i], Zdim[i])) );
        cudaMemcpy3DParms copyParams = {0};

        //Array creation
        copyParams.srcPtr   = make_cudaPitchedPtr(d_inputIm, Xdim[i]*sizeof(float), Xdim[i], Ydim[i]);
        copyParams.dstArray = d_cuArr;
        copyParams.extent   = make_cudaExtent(Xdim[i], Ydim[i], Zdim[i]);
        copyParams.kind     = cudaMemcpyDeviceToDevice;
        gpuErrchk(cudaMemcpy3D(&copyParams));
        //Array creation End

        cudaResourceDesc    texRes;
        memset(&texRes, 0, sizeof(cudaResourceDesc));
        texRes.resType = cudaResourceTypeArray;
        texRes.res.array.array  = d_cuArr;
        cudaTextureDesc     texDescr;
        memset(&texDescr, 0, sizeof(cudaTextureDesc));
        texDescr.normalizedCoords = true;
        if (interp==0){
        	texDescr.filterMode = cudaFilterModePoint;
        }else{
         	texDescr.filterMode = cudaFilterModeLinear;
        }
        texDescr.addressMode[0] = (cudaTextureAddressMode)wrap;
        texDescr.addressMode[1] = (cudaTextureAddressMode)wrap;
        texDescr.readMode = cudaReadModeElementType;

        gpuErrchk(cudaCreateTextureObject(&texImage[i], &texRes, &texDescr, NULL));*/



	//}

}


float** cuda_rotate_image_linear_nearest_parallel(float **image, size_t *Xdim,
		size_t *Ydim, size_t *Zdim, double** ang, int interp, int wrap, int first_call){

	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	if(prop.concurrentKernels==1)
		printf("Device support concurrent kernels\n");
	else if (prop.concurrentKernels==0)
		printf("Device does not support concurrent kernels\n");

    size_t free_byte, total_byte;

    gpuErrchk(cudaMemGetInfo( &free_byte, &total_byte ));

    double free_db = (double)free_byte ;
    double total_db = (double)total_byte ;
    double used_db = total_db - free_db ;

    printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n", used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);

	int mdSizeInfo = NUM_TEX;
	int angle_size = 9;
    float *d_output0, *d_output1, *d_output2, *d_output3, *d_output4;
	double *d_angle0, *d_angle1, *d_angle2, *d_angle3, *d_angle4;

	float **rotated_image;
	gpuErrchk(cudaMallocHost(&rotated_image, mdSizeInfo*sizeof(float*)));

	int maxXdim=0;
	int maxYdim=0;


	for (int i=0; i<mdSizeInfo; i++){

		gpuErrchk(cudaMallocHost(&rotated_image[i], Xdim[i]*Ydim[i]*sizeof(float)));

		CreateTexture(image, Xdim, Ydim, Zdim, interp, wrap, i);
		gpuErrchk(cudaStreamCreate(&stream[i]));

		if(i==0){
			gpuErrchk(cudaMalloc((void **)&d_output0, Xdim[i]*Ydim[i]*sizeof(float)));
			gpuErrchk(cudaMalloc((void**)&d_angle0, angle_size * sizeof(double)));
			gpuErrchk(cudaMemcpy(d_angle0, ang[i], angle_size * sizeof(double), cudaMemcpyHostToDevice));
		}else if(i==1){
			gpuErrchk(cudaMalloc((void **)&d_output1, Xdim[i]*Ydim[i]*sizeof(float)));
			gpuErrchk(cudaMalloc((void**)&d_angle1, angle_size * sizeof(double)));
			gpuErrchk(cudaMemcpy(d_angle1, ang[i], angle_size * sizeof(double), cudaMemcpyHostToDevice));
		}else if(i==2){
			gpuErrchk(cudaMalloc((void **)&d_output2, Xdim[i]*Ydim[i]*sizeof(float)));
			gpuErrchk(cudaMalloc((void**)&d_angle2, angle_size * sizeof(double)));
			gpuErrchk(cudaMemcpy(d_angle2, ang[i], angle_size * sizeof(double), cudaMemcpyHostToDevice));
		}else if(i==3){
			gpuErrchk(cudaMalloc((void **)&d_output3, Xdim[i]*Ydim[i]*sizeof(float)));
			gpuErrchk(cudaMalloc((void**)&d_angle3, angle_size * sizeof(double)));
			gpuErrchk(cudaMemcpy(d_angle3, ang[i], angle_size * sizeof(double), cudaMemcpyHostToDevice));
		}else if(i==4){
			gpuErrchk(cudaMalloc((void **)&d_output4, Xdim[i]*Ydim[i]*sizeof(float)));
			gpuErrchk(cudaMalloc((void**)&d_angle4, angle_size * sizeof(double)));
			gpuErrchk(cudaMemcpy(d_angle4, ang[i], angle_size * sizeof(double), cudaMemcpyHostToDevice));
		}

		if(Xdim[i]>maxXdim)
			maxXdim=Xdim[i];
		if(Ydim[i]>maxYdim)
			maxYdim=Ydim[i];

	}

	int numTh = 32;
	const dim3 blockSize(numTh, numTh, 1);
	int numBlkx = (int)(maxXdim)/numTh;
	if((maxXdim)%numTh>0){
		numBlkx++;
	}
	int numBlky = (int)(maxYdim)/numTh;
	if((maxYdim)%numTh>0){
		numBlky++;
	}
	const dim3 gridSize(numBlkx, numBlky, 1);


	rotate_kernel_normalized_2D_parallel<<<gridSize, blockSize, 0, stream[0]>>>(d_output0, Xdim[0], Ydim[0], d_angle0, texImage[0]);
	gpuErrchk(cudaMemcpyAsync(rotated_image[0], d_output0, Xdim[0]*Ydim[0]*sizeof(float), cudaMemcpyDeviceToHost, stream[0]));
	rotate_kernel_normalized_2D_parallel<<<gridSize, blockSize, 0, stream[1]>>>(d_output1, Xdim[1], Ydim[1], d_angle1, texImage[1]);
	gpuErrchk(cudaMemcpyAsync(rotated_image[1], d_output1, Xdim[1]*Ydim[1]*sizeof(float), cudaMemcpyDeviceToHost, stream[1]));
	rotate_kernel_normalized_2D_parallel<<<gridSize, blockSize, 0, stream[2]>>>(d_output2, Xdim[2], Ydim[2], d_angle2, texImage[2]);
	gpuErrchk(cudaMemcpyAsync(rotated_image[2], d_output2, Xdim[2]*Ydim[2]*sizeof(float), cudaMemcpyDeviceToHost, stream[2]));
	rotate_kernel_normalized_2D_parallel<<<gridSize, blockSize, 0, stream[3]>>>(d_output3, Xdim[3], Ydim[3], d_angle3, texImage[3]);
	gpuErrchk(cudaMemcpyAsync(rotated_image[3], d_output3, Xdim[3]*Ydim[3]*sizeof(float), cudaMemcpyDeviceToHost, stream[3]));
	rotate_kernel_normalized_2D_parallel<<<gridSize, blockSize, 0, stream[4]>>>(d_output4, Xdim[4], Ydim[4], d_angle4, texImage[4]);
	gpuErrchk(cudaMemcpyAsync(rotated_image[4], d_output4, Xdim[4]*Ydim[4]*sizeof(float), cudaMemcpyDeviceToHost, stream[4]));

	//AJ NO FUNCIONA
	//rotate_kernel_normalized_2D_parallel_inside<<<gridSize, blockSize>>>(d_output0, d_output1, d_output2, d_output3, d_output4, Xdim, Ydim, d_angle0, d_angle1, d_angle2, d_angle3, d_angle4, texImage);


	gpuErrchk(cudaDeviceSynchronize());
	for (int i=0; i<mdSizeInfo; i++){
		if(i==0){
			//gpuErrchk(cudaMemcpy(rotated_image[i], d_output0, Xdim[i]*Ydim[i]*sizeof(float), cudaMemcpyDeviceToHost));
			gpuErrchk(cudaFree(d_output0));
			gpuErrchk(cudaFree(d_angle0));
		}else if(i==1){
			//gpuErrchk(cudaMemcpy(rotated_image[i], d_output1, Xdim[i]*Ydim[i]*sizeof(float), cudaMemcpyDeviceToHost));
			gpuErrchk(cudaFree(d_output1));
			gpuErrchk(cudaFree(d_angle1));
		}else if(i==2){
			//gpuErrchk(cudaMemcpy(rotated_image[i], d_output2, Xdim[i]*Ydim[i]*sizeof(float), cudaMemcpyDeviceToHost));
			gpuErrchk(cudaFree(d_output2));
			gpuErrchk(cudaFree(d_angle2));
		}else if(i==3){
			//gpuErrchk(cudaMemcpy(rotated_image[i], d_output3, Xdim[i]*Ydim[i]*sizeof(float), cudaMemcpyDeviceToHost));
			gpuErrchk(cudaFree(d_output3));
			gpuErrchk(cudaFree(d_angle3));
		}else if(i==4){
			//gpuErrchk(cudaMemcpy(rotated_image[i], d_output4, Xdim[i]*Ydim[i]*sizeof(float), cudaMemcpyDeviceToHost));
			gpuErrchk(cudaFree(d_output4));
			gpuErrchk(cudaFree(d_angle4));
		}

	}

    return rotated_image;

}
