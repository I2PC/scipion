#define TILE_DIM 32
#define TILE_DIMH (TILE_DIM/2)
#define TILE_DIM2 (TILE_DIM*2)
#define BLOCK_ROWS 8
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }


#define FILL_GPU(d_array, value, Xdim, Ydim) { \
	block.x=floor(((Xdim)+(TILE_DIM-1))/TILE_DIM); block.y=floor(((Ydim)+TILE_DIMH-1)/(TILE_DIMH)); \
	thread.x=TILE_DIM; thread.y=TILE_DIMH; \
	Kfill<<< block, thread >>>((d_array), (value), (Xdim), (Ydim));\
	cudaThreadSynchronize(); \
	gpuErrchk(cudaGetLastError()); \
	}

// transpose output (note that the desination has Ydim columns!!!)
#define TRANSPOSE_GPU(d_dest, d_source, Xdim, Ydim) { \
	block.x=floor(((Xdim)+(TILE_DIM-1))/TILE_DIM); block.y=floor(((Ydim)+(TILE_DIM-1))/TILE_DIM); \
	thread.x=TILE_DIM; thread.y=BLOCK_ROWS; \
	Ktranspose<<< block, thread >>>((d_dest), (d_source), (Xdim), (Ydim)); \
	cudaThreadSynchronize(); \
	gpuErrchk(cudaGetLastError()); \
	}




inline void gpuAssert(cudaError_t code, const char *file, int line,
			 bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

/****
 sumall - sum no_imgs images in matrix sumObjs, starting from firstImg
****/
__global__
void Ksumall( float* sumObjs, int** array_Img, int no_imgs, int Xdim, int Ydim){

   int x = blockIdx.x*blockDim.x + threadIdx.x;
   int y = blockIdx.y*blockDim.y + threadIdx.y;
   int offset = x+y*Xdim;
   int px_img = Xdim*Ydim;
  
   if ((x<Xdim)&&(y<Ydim)){
//   	printf("no_imgs=%d, (%d, %d) th=(%d,%d) blk=(%d,%d)\n", no_imgs, x, y, threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);
	double sum=0;
	int* ptr;
	for (int img=0; img<no_imgs; img++){
	    	ptr=array_Img[img];
		sum+= (double)ptr[offset];
	}
        sum*=2;
        sumObjs[x+y*Xdim]=(float)sum;
//	printf("data %f\n", sumObjs[x+y+y*Xdim]);

   } //end-if Xdim, Ydim
}//sumall

/****
 mult - c=a+b
****/
__global__
void Kmult( int* a, float* b, int* c, int Xdim, int Ydim){

   int x = blockIdx.x*blockDim.x + threadIdx.x;
   int y = blockIdx.y*blockDim.y + threadIdx.y;
   int offset = x+y*Xdim;

   if ((x<Xdim)&&(y<Ydim))
	c[offset]=int((double)a[offset]*b[offset]);
}//mult

/****
 fill - element=value
****/
template<typename T>
__global__
void Kfill( T* a, T value, int Xdim, int Ydim){

   int x = blockIdx.x*blockDim.x + threadIdx.x;
   int y = blockIdx.y*blockDim.y + threadIdx.y;
   int offset = x+y*Xdim;

   if ((x<Xdim)&&(y<Ydim))
	a[offset]=value;
}//fill



/****
 smooth1: first kernel of constructSmoothHistogramColumn/Row
          * if used for constructSmoothHistogramByRow, rowH (colrow) 
            must be transposed before calling the kernel, and smooth 
            after calling it
	  * Each thread computes a semicolumn

****/

/*********** IT CAN BE OPTIMIZED ************/
/*********** by adding more threads in the Y dim ************/
__global__ void Ksmooth1(float *smooth, int *colrow, const float *listOfWeights, int width, int Xdim, int Ydim)
{
   int x = blockIdx.x*TILE_DIM2 + threadIdx.x;
   int y = blockIdx.y*TILE_DIM2 + threadIdx.y;

   if (x<Xdim){
	   double sumWeightsC = 0;
	   int tmp=0;
	   for(int k = -width; k<=width; ++k){
		if (x+k<0 || x+k>=Xdim)
			continue;
		tmp++;
	// DIRECT_A2D_ELEM(smoothColumnH,i,j) += actualWeightC * DIRECT_A2D_ELEM(columnH,i,j+k);
			float actualWeightC = listOfWeights[abs(k)];
		sumWeightsC += actualWeightC;
		for (size_t i=0; i<TILE_DIM2; ++i){
			if ((y+i)<Ydim){
				smooth[x+(y+i)*Xdim] += actualWeightC * (colrow[(x+k)+(y+i)*Xdim]);
			}	
		}	
	  }	
	  if ((width==0)&&(sumWeightsC==0)){
		printf("tmp=%d sumWeightsC=%f width=%d\n", tmp, sumWeightsC, width);	
  		printf("BLK(%d,%d), TH(%d,%d)\n", blockIdx.x, blockIdx.y, threadIdx.x, threadIdx.y);
	  }
	 double iSumWeightsC=1./sumWeightsC;
	  for (size_t i=0; i<TILE_DIM2; ++i)
		if ((y+i)<Ydim)
			smooth[x+(y+i)*Xdim] *= (float)iSumWeightsC;
   }//end-if(x<=Xdim) 
}

/****
 Ksmooth2: average of rows into the first row
          * There is a thread per column. 
            A thread averages a whole column into the first element

            A thread copies all the values of the first element into 
            the rest of the elements of the column
 	  (Not the best implementation but easy to code)

          * For constructSmoothHistogramsByCol, 
            the smooth matrix must be first transposed

****/
/*********** IT CAN BE OPTIMIZED ************/
/*********** adding more threads in the Y dim ************/
__global__ void Ksmooth2(float *smooth, int Xdim, int Ydim)
{
   int x = blockIdx.x*TILE_DIM2 + threadIdx.x;

   if (x<Xdim){
	   double sum = 0.;
 	   for (size_t y=0; y<Ydim; ++y){
		sum += smooth[x+y*Xdim];
	  }	

	  sum /= Ydim;
 	   for (size_t y=0; y<Ydim; ++y){ 
		smooth[x+y*Xdim]=(float)sum;
	  }	
   }//end-if(x<=Xdim) 
}

/****
K transformGray: each thread deals with a whole column
      1. read value from Iframe - pixval
      2. look for position of first value bigger than pixval - pos
      3. write in IframeTransformedColRow smoothColRos(pos-1) 

****/
__global__ void KtransformGray(const int *Iframe, const int *colrowH, int *IframeTransformColRow, const float* smoothColRow, int Xdim, int Ydim)
{
   int x = blockIdx.x*TILE_DIM2 + threadIdx.x;
 
   if (x<Xdim){
	for (size_t y0=0; y0<Ydim; ++y0){
	   int pixval = Iframe[x+y0*Xdim];
		// upperbounds
	   size_t y1, y2;
	   y2=Ydim-1; // It fails if Ydim==1
   	   for (y1=0; y1<Ydim; ++y1){
		if (colrowH[x+y1*Xdim]>pixval){
			y2=y1-1;
			break;
		}
	   } 
           if (y2<0)
		y2=0;

//	__syncthreads();   
	   IframeTransformColRow[x+y0*Xdim]=(int)smoothColRow[x+y2*Xdim];
   	}//end-for-y0	
   }//end-if-Xdim
}



// Kernel to transpose a matrix

// The kernel assumes that each block deals with a 32x32 tile
// each block has 32x8 threads 
// IMPORTANT: The input has Xdim columns X Ydim rows
//            The output has Ydim columns X Xdim rows
//
//	Optimization: move conditional outside the loop
//
 
template <typename T>
__global__ void Ktranspose(T *odata, const T *idata, int Xdim, int Ydim)
{
__shared__ float tile[TILE_DIM+1][TILE_DIM];
//__shared__ float tile[TILE_DIM][TILE_DIM];
  int x = blockIdx.x * TILE_DIM + threadIdx.x;
  int y = blockIdx.y * TILE_DIM + threadIdx.y;

  if (x<Xdim)
	for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
		if (y+j<Ydim)
		     tile[threadIdx.y+j][threadIdx.x] = idata[(y+j)*Xdim + x];

  __syncthreads();

  x = blockIdx.y * TILE_DIM + threadIdx.x;  // transpose block offset
  y = blockIdx.x * TILE_DIM + threadIdx.y;

  if (x<Ydim)
	  for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
		if (y+j<Xdim)
		     odata[(y+j)*Ydim + x] = tile[threadIdx.x][threadIdx.y + j];
}


/***************************
	KcomputeTVCols

	Each block has 64x1 threads, but is linked to a 64x64 tile 
							(TILE_DIM2xTILE_DIM2)
	* Each thread deals with  two columns in a tile
        * The column-wise partial accumulations are stored in shared memory
        * Eventually, the block-wise partial accumulation is performed 
	  by reduction
	  The reduction assumes that the number of threads is multiple of 2	
          Thus, the shared memory array is filled with zeros for
          threads with idx outside the image
	  The block-wise partial acc. is stored in global memory (d_avgTV)
          The final reduction is perfomed in C
	
	* Note that shared memory is dynamic, so it is necessary to pass 
          the size of partialACC using a third parameter in the kernal 
          configuration: <<grid, block, size_of_partialACC_in_bytes>>Kcompute...

	* Optimization: move conditional outside the loop
*****************************/
template <typename T>
__global__ void KcomputeTVcolumns(T* array,  double* avgTV, int Xdim, int Ydim)
{
	extern __shared__ double partACC[];
	
	int x = blockIdx.x*TILE_DIM2 + threadIdx.x;
	int y = blockIdx.y*TILE_DIM2 + threadIdx.y;
        int tileDim = blockDim.x;
        
	double retvalC=0;
	if (x<Xdim-1){ //two consecutive conlumns are accessed
		for (int j = 0; j < tileDim; ++j)
			if (y+j<Ydim)
		     		retvalC+=abs(array[(y+j)*Xdim + x]-array[(y+j)*Xdim +(x+1)]);
		     	//	retvalC+=abs(((y+j)*Xdim + x)-((y+j)*Xdim +(x+1)));
	}
	partACC[threadIdx.x]=retvalC;
/*	if (blockIdx.x==3 && blockIdx.y==3)
		printf("retvalC=%f",retvalC);
*/
	__syncthreads();
	
	for (int st=tileDim/2; st>0; st>>=1){
		if (threadIdx.x<st)
		    partACC[threadIdx.x] += partACC[threadIdx.x+st];
		__syncthreads();
	}
//	__syncthreads();

	if (threadIdx.x==0){
		avgTV[blockIdx.x+blockIdx.y*gridDim.x]=partACC[0];
//		if (blockIdx.x==3 && blockIdx.y==3)
//			printf("partACC=%f",partACC[0]);
	}
		
	return;
}

