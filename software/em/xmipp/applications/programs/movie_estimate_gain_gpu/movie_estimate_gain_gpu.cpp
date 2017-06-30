#include <cuda_runtime.h>

#include <data/xmipp_program.h>

#include <time.h> // GCF
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/host_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/generate.h>
#include <thrust/equal.h>
#include <thrust/sequence.h>
#include <thrust/for_each.h>
#include "movie_estimate_gain_kernels.cu"


// GPU version part of the class

// Used for nested sort
struct sort_functor
{
  thrust::device_ptr<int> data;
  int dsize;
  __host__ __device__
  void operator()(int start_idx)
  {
    thrust::sort(thrust::device, data+(dsize*start_idx), data+(dsize*(start_idx+1)));
  }
};

int my_mod_start;
int no_el;

//used for vectorized sort
int my_modC(){
  return (my_mod_start++)%no_el;
}

int my_modR(){
  return (my_mod_start++)/no_el;
}



class ProgMovieEstimateGainGPU: public XmippProgram
{
public:
	FileName fnIn; // Set of input images
	FileName fnRoot; // Correction image
	int Niter; // Number of iterations
	float maxSigma, sigmaStep;  //it was double
	bool singleReference;
public:
    void defineParams();
    void readParams();
    void show();
    void run();

    void produceSideInfo();
    void computeHistogramsGPU(const int* d_IframeIdeal);
    void constructSmoothHistogramsByColumnGPU(const float *listOfWeights, int width);
    void constructSmoothHistogramsByRowGPU(const float *listOfWeights, int width);
    void transformGrayValuesColumnGPU(const int* d_Iframe, int* d_IframeTransformedColumn);
    void transformGrayValuesRowGPU(const int* d_Iframe, int* d_IframeTransformed, bool out_trans=false);
    size_t selectBestSigmaByColumnGPU(const int* d_Iframe);
    size_t selectBestSigmaByRowGPU(const int* d_Iframe);
    double computeTVColumnsGPU(int* d_I); 
    double computeTVRowsGPU(int* d_I); 
  
    void averageCorrectionGPU(void);

public:
	MetaData mdIn;
	MultidimArray<float> sumObs; // it was double
	Image<float> ICorrection;         // it was double
	std::vector<float> listOfSigmas;  // it was double
	std::vector<int> listOfWidths; // it was double
	std::vector<float *> listOfWeights; // it was double

	int Xdim, Ydim;

	// Note that most GPU data is float instead of double
	// d_ means GPU data, h_ means host (uP) data 
	cudaDeviceProp GPUprop;
        int  GPUdevice; 
	int* d_columnH, *d_rowH; //matrix
	int* d_Ht; // Temporary col/rowH (used for transposition) 
	float* d_smoothColumnH, * d_smoothRowH, * d_sumObs; //matrix
	float* d_smoothHt; // temporary smooth matrix(used for transposition)
	float* d_sumIdeal;

	float* d_ICorrection; 		//Image
	int* d_IframeIdeal; 		//Image
	int* d_Iframet1;		//Image (used for transposition)
	int* d_IframeTransformed; 	//Image
	int* d_Iframet2;		//Image (used for transposition)
	float** d_listOfWeights;	//vector of pointers
	float* d_listOfSigmas;		//vector
	float* d_weights;                 	// vector
	int** Iframe_vecGPU;			// vector of images
	std::vector<float *> listOfWeightsGPU; 	// it was double
	double* d_avgTVcol, *d_avgTVrow;	// vector (d_avgTV[0] is eq. to avgTV)
	double* avgTVcolTmp, *avgTVrowTmp;	// temporal vector used for reduction
	double* d_avgCorr;			// vector (d_avgTV[0] is eq. to avgTV)
	double* avgCorrTmp;			// temporal vector used for reduction


	clock_t start, end;
	dim3 block, thread;
	size_t sz_imgFL, sz_imgINT;
	size_t sz_avgTVcol, nel_avgTVcol, sz_avgTVrow, nel_avgTVrow;
	size_t sz_avgCorr, nel_avgCorr;
};

void ProgMovieEstimateGainGPU::defineParams()
{
    addUsageLine("Estimate the gain image for a movie");
    addParamsLine(" -i <movie>: Input movie");
    addParamsLine(" [--oroot <fn=\"estimated\">]: Estimated corrections and gains");
    addParamsLine("                             :+(Ideal=Observed*Corr, Observed=Ideal*Gain)");
    addParamsLine(" [--iter <N=3>]: Number of iterations");
    addParamsLine(" [--maxSigma <s=3>]: Maximum number of neighbour rows/columns to analyze");
    addParamsLine(" [--sigmaStep <s=0.5>]: Step size for sigma");
    addParamsLine(" [--singleRef] : Use a single histogram reference");
    addParamsLine("               :+This assumes that there is noi image contamination or carbon holes");
}

void ProgMovieEstimateGainGPU::readParams()
{
	fnIn=getParam("-i");
	fnRoot=getParam("--oroot");
	Niter=getIntParam("--iter");
	maxSigma=getDoubleParam("--maxSigma");
	sigmaStep=getDoubleParam("--sigmaStep");
	singleReference=checkParam("--singleRef");
}

void ProgMovieEstimateGainGPU::produceSideInfo()
{
	mdIn.read(fnIn);
	mdIn.removeDisabled();
	if (mdIn.size()==0)
		exit(0);
	Image<float> Iframe; // it was double
	FileName fnFrame;
	mdIn.getValue(MDL_IMAGE,fnFrame,mdIn.firstObject());
	Iframe.read(fnFrame);
	Xdim=XSIZE(Iframe());
	Ydim=YSIZE(Iframe());
	sz_imgFL=Xdim*Ydim*sizeof(float);
	sz_imgINT=Xdim*Ydim*sizeof(int);

	ICorrection().resizeNoCopy(Ydim,Xdim);
	ICorrection().initConstant(1);
	sumObs.initZeros(Ydim,Xdim);

	// Send ICorrection to GPU
	gpuErrchk(cudaMalloc(&d_ICorrection, sz_imgFL));
	FILL_GPU(d_ICorrection, (float)1., Xdim, Ydim)
		
	// Allocate and intialize GPU data	
	gpuErrchk(cudaMalloc( &d_sumIdeal, sz_imgFL)); // no need to initialized
	gpuErrchk(cudaMalloc( &d_sumObs, sz_imgFL)); // no need to initialized
	gpuErrchk(cudaMalloc( &d_rowH, sz_imgINT)); // no need to initialized
	gpuErrchk(cudaMalloc( &d_columnH, sz_imgINT)); // no need to initialized
	gpuErrchk(cudaMalloc( &d_Ht, sz_imgINT)); // temporal vector, no need to initialized
	gpuErrchk(cudaMalloc( &d_smoothColumnH, sz_imgFL));
	gpuErrchk(cudaMalloc( &d_smoothRowH, sz_imgFL));
	gpuErrchk(cudaMalloc( &d_smoothHt, sz_imgFL));
	gpuErrchk(cudaMalloc( &d_Iframet1, sz_imgINT));
	gpuErrchk(cudaMalloc( &d_Iframet2, sz_imgINT));
 	block.x=block.y=block.z=1;
 	thread.x=thread.y=thread.z=1;
 
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		mdIn.getValue(MDL_IMAGE,fnFrame,__iter.objId);
		Iframe.read(fnFrame);
		sumObs+=Iframe();
	}
	sumObs*=2;
	
	gpuErrchk(cudaMemcpy(d_sumObs, &sumObs(0,0), sz_imgFL, cudaMemcpyHostToDevice));
	
	nel_avgTVcol = floor((Xdim+(TILE_DIM2-1))/TILE_DIM2)*floor((Ydim+(TILE_DIM2-1))/TILE_DIM2);
	sz_avgTVcol = sizeof(double)*nel_avgTVcol;
	gpuErrchk(cudaMalloc( &d_avgTVcol, sz_avgTVcol)); // no need to be initialized
        avgTVcolTmp = (double*)malloc(sz_avgTVcol);	
	nel_avgTVrow = floor((Ydim+(TILE_DIM2-1))/TILE_DIM2)*floor((Xdim+(TILE_DIM2-1))/TILE_DIM2);
	sz_avgTVrow = sizeof(double)*nel_avgTVrow;
	gpuErrchk(cudaMalloc( &d_avgTVrow, sz_avgTVrow)); // no need to initialized
        avgTVrowTmp = (double*)malloc(sz_avgTVrow);	

	nel_avgCorr = (floor((Xdim+(TILE_DIM2-1))/TILE_DIM2))*(floor((Ydim+(TILE_DIM2-1))/TILE_DIM2));
	sz_avgCorr = sizeof(double)*nel_avgCorr;
	gpuErrchk(cudaMalloc( &d_avgCorr, sz_avgCorr)); // no need to initialized
        avgCorrTmp = (double*)malloc(sz_avgCorr);	

	//****

	// Initialize sigma values
	for (float sigma=0; sigma<=maxSigma; sigma+=sigmaStep)
		listOfSigmas.push_back(sigma);
	
	for (size_t i=0; i<listOfSigmas.size(); ++i)
	{
		int jmax=ceil(3*listOfSigmas[i]); 
		listOfWidths.push_back(jmax+1);
		float *weights=new float[jmax+1];    //GCF:it was double and (jmax)
		float K=1; //GCF: it was double
                if (listOfSigmas[i]>0)
                   K=-0.5/(listOfSigmas[i]*listOfSigmas[i]);
		for (int j=0; j<=jmax; ++j)
			weights[j]=exp(K*j*j);
		listOfWeights.push_back(weights);
	}
	// Copy to GPU
	size_t sz_listS = sizeof(float)*listOfSigmas.size();
	gpuErrchk(cudaMalloc(&d_listOfSigmas, sz_listS));
	gpuErrchk(cudaMemcpy(d_listOfSigmas, &listOfSigmas[0], sz_listS, cudaMemcpyHostToDevice));
	size_t sz_listW = sizeof(float*)*listOfSigmas.size();
	gpuErrchk(cudaMalloc(&d_listOfWeights, sz_listW));
	
	for (size_t i=0; i<listOfSigmas.size(); ++i)
	{
		int jmax=ceil(3*listOfSigmas[i]);
		size_t sz_weights=sizeof(float)*(jmax+1); //GCF: it was (jmax)
		gpuErrchk(cudaMalloc(&d_weights, sz_weights));
	        gpuErrchk(cudaMemcpy(d_weights, listOfWeights[i], sz_weights, cudaMemcpyHostToDevice));
		listOfWeightsGPU.push_back(d_weights);
	}
		
	
}

void ProgMovieEstimateGainGPU::show()
{
	if (verbose==0)
		return;
	std::cout
	<< "Input movie:     " << fnIn            << std::endl
	<< "Output rootname: " << fnRoot          << std::endl
	<< "N. Iterations:   " << Niter           << std::endl
	<< "Max. Sigma:      " << maxSigma        << std::endl
	<< "Sigma step:      " << sigmaStep       << std::endl
	<< "Single ref:      " << singleReference << std::endl
	;
}


void ProgMovieEstimateGainGPU::run()
{
	// Select GPU and show
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	GPUdevice = -1;
	for (int device = 1; device < deviceCount; ++device) {
	    cudaDeviceProp deviceProp;
	    if (!cudaSetDevice(device)){
		    cudaGetDeviceProperties(&GPUprop, device);
		    printf("Device %d has compute capability %d.%d.\n",
        		   device, GPUprop.major, GPUprop.minor);
		    GPUdevice=device;
		    break;
	    }
	}

	if (GPUdevice==-1){
		std::cout << "No GPU found" << std::endl;
		while(1);
	}
	else if (GPUdevice!=0)
		std::cout << "GPU device 0 not selected, risk of not choosing the fastest one" << std::endl;
	std::cout << "* Total global mem = " << GPUprop. totalGlobalMem << " bytes" << std::endl;

	
	produceSideInfo();

	FileName fnFrame;
	Image<int> Iframe;
	MultidimArray<int> IframeTransformed, IframeIdeal;
	MultidimArray<float> &mICorrection=ICorrection();   // it was double. used in the final step. CPU since is only once


	// Copy all int images onto GPU
	int im=0;
	size_t sz_Iframe_vec=sizeof(int*)*mdIn.size();
	size_t px_img=Xdim*Ydim;

	Iframe_vecGPU = (int**)malloc(sz_Iframe_vec);
 
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		mdIn.getValue(MDL_IMAGE,fnFrame,__iter.objId);
		Iframe.read(fnFrame);

		gpuErrchk(cudaMalloc(&Iframe_vecGPU[im], sz_imgINT));
		gpuErrchk(cudaMemcpy(Iframe_vecGPU[im++], &Iframe(0,0), sz_imgINT, cudaMemcpyHostToDevice)); 
	}

 
	// Malloc several data
	gpuErrchk(cudaMalloc(&d_IframeIdeal, sz_imgINT));
	gpuErrchk(cudaMalloc(&d_IframeTransformed, sz_imgINT));
	
	
//*********** GPU Computation
	
	for (int n=0; n<Niter; n++)
	{
		std::cout << "Iteration " << n << std::endl;
        	FILL_GPU(d_sumIdeal, (float)0., Xdim, Ydim);
		
		im=0;
		FOR_ALL_OBJECTS_IN_METADATA(mdIn)
		{
			mdIn.getValue(MDL_IMAGE,fnFrame,__iter.objId);
			std::cout << "   Frame " << fnFrame << std::endl;
			
			dim3 block1(floor((Xdim+(TILE_DIM-1))/TILE_DIM),floor((Ydim+(TILE_DIMH-1))/(TILE_DIMH)),1);
			dim3 thread1( TILE_DIM, TILE_DIMH);
			Kmult<<< block1, thread1 >>>(Iframe_vecGPU[im++], d_ICorrection, d_IframeIdeal, Xdim, Ydim);
			cudaThreadSynchronize();
			gpuErrchk(cudaGetLastError());
			
		        computeHistogramsGPU(d_IframeIdeal); // GPU test	

			size_t bestSigmaColGPU = selectBestSigmaByColumnGPU(d_IframeIdeal);
			std::cout << "      sigmaColGPU: " << listOfSigmas[bestSigmaColGPU] << std::endl;
		 	
			size_t bestSigmaRowGPU = selectBestSigmaByRowGPU(d_IframeIdeal);
			std::cout << "      sigmaRowGPU: " << listOfSigmas[bestSigmaRowGPU] << std::endl;

		constructSmoothHistogramsByRowGPU(listOfWeightsGPU[bestSigmaRowGPU],listOfWidths[bestSigmaRowGPU]);   //smooth
			transformGrayValuesRowGPU(d_IframeIdeal, d_IframeTransformed, false); // IframeTransformed is not transposed
			dim3 block2(floor((Xdim+(TILE_DIM-1))/TILE_DIM),floor((Ydim+(TILE_DIMH-1))/(TILE_DIMH)),1);
			dim3 thread2( TILE_DIM, TILE_DIMH);
			Ksum<<< block2, thread2 >>>(d_sumIdeal,  d_IframeTransformed, d_sumIdeal,Xdim, Ydim);
			cudaThreadSynchronize();
			gpuErrchk(cudaGetLastError());

			constructSmoothHistogramsByColumnGPU(listOfWeightsGPU[bestSigmaColGPU],listOfWidths[bestSigmaColGPU]);
			transformGrayValuesColumnGPU(d_IframeIdeal, d_IframeTransformed);
			Ksum<<< block2, thread2 >>>(d_sumIdeal,  d_IframeTransformed, d_sumIdeal,Xdim, Ydim);
			cudaThreadSynchronize();
			gpuErrchk(cudaGetLastError());
			
	
		}//for-ALL_OBJECTS

	
		// kernel with partial average
		averageCorrectionGPU();
		
	}//for-iter


	gpuErrchk(cudaMemcpy(&DIRECT_A2D_ELEM(mICorrection,0,0), d_ICorrection, sz_imgFL, cudaMemcpyDeviceToHost));	
	
	ICorrection.write(fnRoot+"_correction.xmp");
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(mICorrection)
		if (DIRECT_A2D_ELEM(mICorrection,i,j)>1e-5)
			DIRECT_A2D_ELEM(mICorrection,i,j)=1.0/DIRECT_A2D_ELEM(mICorrection,i,j);
		else
			DIRECT_A2D_ELEM(mICorrection,i,j)=1;
	ICorrection.write(fnRoot+"_gain.xmp");

	// Free memory
	for (int i=0; i< mdIn.size(); i++){
		cudaFree(Iframe_vecGPU[i]);
	}

	for (int i=0; i< listOfWeightsGPU.size(); i++){
		cudaFree(listOfWeightsGPU[i]);
	}
		
	cudaFree(d_rowH);
	cudaFree(d_columnH);
	cudaFree(d_Ht);
	cudaFree(d_smoothColumnH);
	cudaFree(d_smoothRowH);
	cudaFree(d_smoothHt);
	cudaFree(d_sumObs);
	cudaFree(d_sumIdeal);
	cudaFree(d_ICorrection);
	cudaFree(d_IframeIdeal);
	cudaFree(d_IframeTransformed);
	cudaFree(d_Iframet1);
	cudaFree(d_Iframet2);
	cudaFree(d_avgTVcol);
	cudaFree(d_avgTVrow);
	free(Iframe_vecGPU);
    	cudaFree(d_avgCorr);
    	free(avgCorrTmp);
    	free(avgTVcolTmp);
    	free(avgTVrowTmp);
	
}

// compute colH and rowH of image in GPU (Iframe)
// it's always applied to d_IframeIdeal
// rowH is transposed (it's faster his way when used by other kernels)
void ProgMovieEstimateGainGPU::computeHistogramsGPU(const int* d_IframeIdeal)
{

////		clock_t start = clock();
		gpuErrchk(cudaMemcpy(d_Ht, d_IframeIdeal, sz_imgINT, cudaMemcpyDeviceToDevice));	
		thrust::device_ptr<int> th_d_colHt = thrust::device_pointer_cast(d_Ht);
		thrust::host_vector<int> h_segmentsC(Xdim*Ydim);
		my_mod_start = 0;
		no_el = Xdim;
		thrust::generate(h_segmentsC.begin(), h_segmentsC.end(), my_modC);
		thrust::device_vector<int> d_segmentsC = h_segmentsC;
		thrust::stable_sort_by_key(th_d_colHt, th_d_colHt+Xdim*Ydim, d_segmentsC.begin());
		thrust::stable_sort_by_key(d_segmentsC.begin(), d_segmentsC.end(), th_d_colHt);
		cudaDeviceSynchronize();


		// transpose output (note that the input has Ydim columns!!!)
		TRANSPOSE_GPU(d_columnH, d_Ht, Ydim, Xdim);

// VECTORIZED SORT: much better than CPU
		gpuErrchk(cudaMemcpy(d_Ht, d_IframeIdeal, sz_imgINT, cudaMemcpyDeviceToDevice));	
		thrust::device_ptr<int> th_d_rowH = thrust::device_pointer_cast(d_Ht);
		thrust::host_vector<int> h_segmentsR(Xdim*Ydim);
		my_mod_start = 0;
		no_el = Xdim;
		thrust::generate(h_segmentsR.begin(), h_segmentsR.end(), my_modR);
		thrust::device_vector<int> d_segmentsR= h_segmentsR;
		thrust::stable_sort_by_key(th_d_rowH, th_d_rowH+Xdim*Ydim, d_segmentsR.begin());
		thrust::stable_sort_by_key(d_segmentsR.begin(), d_segmentsR.end(), th_d_rowH);
		cudaDeviceSynchronize();
		
        	// transpose output (note that the input has Ydim columns!!!)
		TRANSPOSE_GPU(d_rowH, d_Ht, Xdim, Ydim);
}

// argument weights is a pointer to a device array with the list of weights
void ProgMovieEstimateGainGPU::constructSmoothHistogramsByColumnGPU(const float *weights, int width)
{

	//smoothColumnH.initZeros(columnH);
        FILL_GPU(d_smoothColumnH, (float)0., Xdim, Ydim);
	
        dim3 block1(floor((Xdim+(TILE_DIM2-1))/TILE_DIM2),floor((Ydim+(TILE_DIM2-1))/(TILE_DIM2)),1);
	dim3 thread1( TILE_DIM2, 1, 1);
	Ksmooth1<<< block1, thread1>>>(d_smoothColumnH, d_columnH, weights, width, Xdim, Ydim);
	cudaThreadSynchronize();
	gpuErrchk(cudaGetLastError());


	
	if (singleReference)
	{
		// transpose smoothColumnH
		TRANSPOSE_GPU(d_smoothHt, d_smoothColumnH, Xdim, Ydim);

		// Average all rows and replicate in all rows
		dim3 block4(floor((Ydim+(TILE_DIM2-1))/TILE_DIM2),1,1);
		dim3 thread4( TILE_DIM2, 1);
		Ksmooth2<<< block4, thread4 >>>(d_smoothHt, Ydim, Xdim);
		cudaThreadSynchronize();
		gpuErrchk(cudaGetLastError());

		// transpose back smoothColumnH
		TRANSPOSE_GPU(d_smoothColumnH, d_smoothHt, Ydim, Xdim);

	}

}

// argument weights is a pointer to a device array with the list of weights
// d_smoothRowH is transposed (this is the best way to use it later)
void ProgMovieEstimateGainGPU::constructSmoothHistogramsByRowGPU(const float *weights, int width)
{

	//smoothColumnH.initZeros(columnH);
        FILL_GPU(d_smoothRowH, (float)0., Xdim, Ydim);
	
        //TRANSPOSE_GPU(d_Ht, d_rowH,  Xdim, Ydim);
	//TRANSPOSE_GPU(d_smoothHt, d_smoothRowH, Xdim, Ydim);
        
	dim3 block1(floor((Ydim+(TILE_DIM2-1))/TILE_DIM2),floor((Xdim+(TILE_DIM2-1))/(TILE_DIM2)),1);
	dim3 thread1( TILE_DIM2, 1, 1);
	Ksmooth1<<< block1, thread1>>>(d_smoothRowH, d_rowH, weights, width, Ydim, Xdim);
	cudaThreadSynchronize();
	gpuErrchk(cudaGetLastError());


	
	if (singleReference)
	{
		// transpose smoothColumnH
		TRANSPOSE_GPU(d_smoothHt, d_smoothRowH, Ydim, Xdim);

		// Average all rows and replicate in all rows
		dim3 block4(floor((Ydim+(TILE_DIM2-1))/TILE_DIM2),1,1);
		dim3 thread4( TILE_DIM2, 1);
		Ksmooth2<<< block4, thread4 >>>(d_smoothHt, Ydim, Xdim);
		cudaThreadSynchronize();
		gpuErrchk(cudaGetLastError());
	

		// transpose back smoothColumnH
		TRANSPOSE_GPU(d_smoothRowH, d_smoothHt, Xdim, Ydim);

	}

	
}


void ProgMovieEstimateGainGPU::transformGrayValuesColumnGPU(const int* d_Iframe, int* d_IframeTransformedColumn)
{

	FILL_GPU(d_IframeTransformedColumn, 0, Xdim, Ydim); // for debuging purposes

	dim3 block1(floor((Xdim+(TILE_DIM-1))/TILE_DIM),floor((Ydim+(TILE_DIMH-1))/TILE_DIMH),1);
	dim3 thread1( TILE_DIM, TILE_DIMH);
		
      	KtransformGray3<<<block1,thread1>>>(d_Iframe, d_columnH, d_IframeTransformedColumn, d_smoothColumnH, Xdim, Ydim);

	cudaThreadSynchronize();
	gpuErrchk(cudaGetLastError());


}

// Note that it is posible to select if the output is transposed or not
// When computeTV is called after, the ouput should remain transposed (true)
// When it is called from the main program it must not be transposed (false)
// Note that d_smoothRowH comes from ComputeHistogramByRow already transposed
void ProgMovieEstimateGainGPU::transformGrayValuesRowGPU(const int* d_Iframe, int* d_IframeTransformed, bool out_trans)
{
	// transpose smoothRowH, d_Iframe, d_colH
	
	TRANSPOSE_GPU(d_Iframet1, d_Iframe, Xdim, Ydim);

	

	int* p_IframeTransformed=(!out_trans?d_Iframet2:d_IframeTransformed);

	FILL_GPU(p_IframeTransformed, 0, Ydim, Xdim);


	dim3 block1(floor((Ydim+(TILE_DIM-1))/TILE_DIM),floor((Xdim+(TILE_DIMH-1))/TILE_DIMH),1);
	dim3 thread1( TILE_DIM, TILE_DIMH);
		
	KtransformGray3<<<block1,thread1>>>(d_Iframet1, d_rowH,   p_IframeTransformed, d_smoothRowH, Ydim, Xdim);
	cudaThreadSynchronize();
	gpuErrchk(cudaGetLastError());

		
// When out_trans==true we do nothing, since the output of KtransformGray
//   is transposed
	if (!out_trans){
	   TRANSPOSE_GPU(d_IframeTransformed, d_Iframet2, Ydim, Xdim);

	}
}

double ProgMovieEstimateGainGPU::computeTVColumnsGPU(int* d_I){

	
		dim3 block1(floor((Xdim+(TILE_DIM2-1))/TILE_DIM2),floor((Ydim+(TILE_DIM2-1))/TILE_DIM2),1);
		dim3 thread1( TILE_DIM2, 1, 1);
		KcomputeTVcolumns<int><<<block1, thread1, sizeof(double)*TILE_DIM2>>>(d_I, d_avgTVcol, Xdim, Ydim);

		cudaThreadSynchronize();
		gpuErrchk(cudaGetLastError());
		double avgTVgpu=0.;
		gpuErrchk(cudaMemcpy(avgTVcolTmp, d_avgTVcol, sz_avgTVcol, cudaMemcpyDeviceToHost));
		for (size_t h=0; h<nel_avgTVcol; ++h){
			avgTVgpu+=avgTVcolTmp[h];
		}
		avgTVgpu/=(Xdim-1)*Ydim;

		return avgTVgpu;
}

double ProgMovieEstimateGainGPU::computeTVRowsGPU(int* d_I){
		dim3 block1(floor((Ydim+(TILE_DIM2-1))/TILE_DIM2),floor((Xdim+(TILE_DIM2-1))/TILE_DIM2),1);
		dim3 thread1( TILE_DIM2, 1, 1);
		KcomputeTVcolumns<int><<<block1, thread1, sizeof(double)*TILE_DIM2>>>(d_I, d_avgTVrow, Ydim, Xdim);

		cudaThreadSynchronize();
		gpuErrchk(cudaGetLastError());
		double avgTVgpu=0.;
		gpuErrchk(cudaMemcpy(avgTVrowTmp, d_avgTVrow, sz_avgTVrow, cudaMemcpyDeviceToHost));
		for (size_t h=0; h<nel_avgTVrow; ++h){
			avgTVgpu+=avgTVrowTmp[h];
		}
		avgTVgpu/=(Ydim-1)*Xdim;
		return avgTVgpu;
}


size_t ProgMovieEstimateGainGPU::selectBestSigmaByRowGPU(const int* d_Iframe)
{
	double bestAvgTV=1e38;
	size_t best_s=0;

	for(size_t s = 0; s< listOfWeights.size(); ++s)
	{
		constructSmoothHistogramsByRowGPU(listOfWeightsGPU[s],listOfWidths[s]);
		transformGrayValuesRowGPU(d_Iframe,d_IframeTransformed,true); // note that d_IframeTransformed is tranposed (so rows are columns)
	        // no need to create/using a KcomputeTVrow kernel since the transformed image is transposed
		double avgTVgpu=computeTVRowsGPU(d_IframeTransformed);
		if (avgTVgpu<bestAvgTV){
		
			bestAvgTV=avgTVgpu;
			best_s=s;
		}
	}//for-s

	return best_s;
}

size_t ProgMovieEstimateGainGPU::selectBestSigmaByColumnGPU(const int* d_Iframe)
{
	double bestAvgTV=1e38;
	size_t best_s=0;

	for(size_t s = 0; s< listOfWeights.size(); ++s)
	{
		constructSmoothHistogramsByColumnGPU(listOfWeightsGPU[s],listOfWidths[s]);
		transformGrayValuesColumnGPU(d_Iframe,d_IframeTransformed);
		
		double avgTVgpu=computeTVColumnsGPU(d_IframeTransformed);
		if (avgTVgpu<bestAvgTV)
		{
			bestAvgTV=avgTVgpu;
			best_s=s;
		}
	}//for-s
	return best_s;
}

void ProgMovieEstimateGainGPU::averageCorrectionGPU(void){
		// Average mIcorrection with sumObs
		// leave partial summation in avgCorriTmp
		dim3 block4(floor((Xdim+(TILE_DIM2-1))/TILE_DIM2),floor((Ydim+(TILE_DIM2-1))/TILE_DIM2),1);
		dim3 thread4( TILE_DIM2, 1);
		Kavg1<<< block4, thread4, sizeof(double)*TILE_DIM2>>>(d_ICorrection, d_sumObs, d_sumIdeal, d_avgCorr, Xdim, Ydim);
		cudaThreadSynchronize();
		gpuErrchk(cudaGetLastError());
		
		// Finish average computation
		double avgCorr=0.;
		gpuErrchk(cudaMemcpy(avgCorrTmp, d_avgCorr, sz_avgCorr, cudaMemcpyDeviceToHost));
		for (size_t h=0; h<nel_avgCorr; ++h){
			avgCorr+=avgCorrTmp[h];
//			std::cout << "avgTVrowTmp[" << h << "]=" << avgTVrowTmp[h] << std::endl;
		}
		avgCorr/=Ydim*Xdim;

		// Normalize mIcorrection
		dim3 block1(floor((Xdim+(TILE_DIM-1))/TILE_DIM),floor((Ydim+(TILE_DIMH-1))/(TILE_DIMH)),1);
		dim3 thread1( TILE_DIM, TILE_DIMH);
		KdivC<<< block1, thread1 >>>(d_ICorrection, (float)avgCorr, Xdim, Ydim);
		cudaThreadSynchronize();
		gpuErrchk(cudaGetLastError());

	

}



RUN_XMIPP_PROGRAM(ProgMovieEstimateGainGPU)

/*	std::cout << "Test on small matrix" << std::endl;

	int xd=3;
	int yd=4;
  	size_t sz_t= sizeof(int)*xd*yd;	
	int h_t[]={11,21,51, 1002,32,22, 13,53,33, 44,84,34};

	int *h_t2=(int*)malloc(sz_t);
	int* d_t;
	gpuErrchk(cudaMalloc(&d_t,sz_t));
	gpuErrchk(cudaMemcpy(d_t, h_t, sz_t, cudaMemcpyHostToDevice));
	thrust::device_ptr<int> th_d_t = thrust::device_pointer_cast(d_t);
	thrust::host_vector<int> h_segR(xd*yd);
	my_mod_start = 0;
	no_el = xd;

	thrust::generate(h_segR.begin(), h_segR.end(), my_modC);
	for (int i=0; i<xd*yd; i++)
1	std::cout<<h_segR[i] << " ";
	std::cout << std::endl;
	thrust::device_vector<int> d_segR= h_segR;
	thrust::stable_sort_by_key(th_d_t, th_d_t+xd*yd, d_segR.begin());
	thrust::stable_sort_by_key(d_segR.begin(), d_segR.end(), th_d_t);
	cudaDeviceSynchronize();
	std::cout<<"after sorting"<<std::endl;
	gpuErrchk(cudaMemcpy(h_t2, d_t, sz_t, cudaMemcpyDeviceToHost));
	std::cout<<"after copying"<<std::endl;
	isequal(h_t,h_t2, xd*yd, xd*yd);



 std::cout << "avgTV=" << avgTV << "; GPU avgTV=" << avgTVgpu << std::endl;
*/
