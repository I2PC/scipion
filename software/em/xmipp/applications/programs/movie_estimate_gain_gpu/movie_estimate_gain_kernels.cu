#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
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
void sumall( float* sumObjs, float** array_Img, int no_imgs, int Xdim, int Ydim){

   int x = blockIdx.x*blockDim.x + threadIdx.x;
   int y = blockIdx.y*blockDim.y + threadIdx.y;
   int offset = x+y*Xdim;
   int px_img = Xdim*Ydim;

   if ((x<Xdim)&&(y<Ydim)){
//   	printf("no_imgs=%d, (%d, %d) th=(%d,%d) blk=(%d,%d)\n", no_imgs, x, y, threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);
	float sum=0;
	float* ptr;
	for (int img=0; img<no_imgs; img++){
	    	ptr=array_Img[img];
		sum+= ptr[offset];
	}
        sum*=2;
        sumObjs[x+y*Xdim]=sum;
//	printf("data %f\n", sumObjs[x+y+y*Xdim]);

   } //end-if Xdim, Ydim
} 
