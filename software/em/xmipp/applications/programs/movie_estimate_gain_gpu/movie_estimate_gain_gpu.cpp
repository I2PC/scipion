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

double computeTVColumns(MultidimArray<int> &I);
double computeTVRows(MultidimArray<int> &I);


template <typename T>
 bool isequal(const T* A, const T* B, int size, int show=10, float th_perc=.1)  { 
	bool ok=true;
	float max_perc=0.;
	float tmp_perc;
	for (int i=0; i<size; i++){
		if (A[i]!=B[i]){
			if (A[i]==0.){
				ok=false;
				std::cout << "A[i]=0!!!" << std::endl;
				break;
			}
			tmp_perc = 100.*fabs(1.-B[i]/A[i]);
			if (tmp_perc>max_perc)
				max_perc=tmp_perc;
			if (max_perc > th_perc)
				ok=false;
		}	
	}	

	if (!ok){
		std::cout << "MAX ERR = " << max_perc << "\%" << std::endl;
		for (int i=0; i<show; i++){
			std::cout << "A[" <<i<<"]="<< A[i];
			std::cout << " B[" <<i<<"]="<< B[i] << std::endl;
		}
        }
	return ok;
}

// Check if matrix A with Xdim columns is the transposed of B (with Ydim columns)
template <typename T>
 bool isequalT(const T* A, const T* B, int xdim, int ydim, int show=10, float th_perc=.1)  { 
	bool ok=true;
	float max_perc=0.;
	float tmp_perc;
		for (int i=0; i<ydim; i++){
		for (int j=0; j<xdim; j++){
			if (A[j+i*xdim]!=B[i+j*ydim]){
				if (A[j+i*xdim]==0.){
					ok=false;
					std::cout << "A[j+i*ydim]=0!!!" << std::endl;
					break;
				}
			tmp_perc = 100.*fabs(1-B[i+j*ydim]/A[j+i*xdim]);
			if (tmp_perc>max_perc)
				max_perc=tmp_perc;
			if (max_perc > th_perc)
				ok=false;
			}	
		}	
	}

	if (!ok){
		std::cout << "MAX ERR = " << max_perc << "\%" << std::endl;
		for (int i=0; i<show; i++){
			std::cout << "A[" <<i<<"]="<< A[i];
			std::cout << " B'[" <<i<<"]="<< B[i*ydim] << std::endl;
		}
        }
	return ok;
} 

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
    void computeHistograms(const MultidimArray<int> &Iframe);
    void computeHistogramsGPU(const int* d_IframeIdeal);
    void normalizeHistograms();
    void invertHistograms();

    void constructSmoothHistogramsByColumn(const float *listOfWeights, int width);
    void constructSmoothHistogramsByRow(const float *listOfWeights, int width);
    void constructSmoothHistogramsByColumnGPU(const float *listOfWeights, int width);
    //*** REMOVE
    void constructSmoothHistogramsByColumnTEST(const float *listOfWeights, int width);
    //***********
    void constructSmoothHistogramsByRowGPU(const float *listOfWeights, int width);
    void transformGrayValuesColumn(const MultidimArray<int> &Iframe, MultidimArray<int> &IframeTransformedColumn);
    void transformGrayValuesColumnGPU(const int* d_Iframe, int* d_IframeTransformedColumn);
    void transformGrayValuesRow(const MultidimArray<int> &Iframe, MultidimArray<int> &IframeTransformedRow);
    void computeTransformedHistograms(const MultidimArray<int> &Iframe);

    size_t selectBestSigmaByColumn(const MultidimArray<int> &Iframe);
    size_t selectBestSigmaByColumnGPU_TEST(const MultidimArray<int> &Iframe, const int* d_Iframe);
    size_t selectBestSigmaByRow(const MultidimArray<int> &Iframe);

  //  int my_modXdim(void);	// used for vectorized sort



public:
	MetaData mdIn;
	MultidimArray<int> columnH,rowH, aSingleColumnH, aSingleRowH;
	MultidimArray<float> smoothColumnH, smoothRowH;
	MultidimArray<double> sumObs; // it was double
	Image<float> ICorrection;         // it was double
	std::vector<float> listOfSigmas;  // it was double
	std::vector<int> listOfWidths; // it was double
	std::vector<float *> listOfWeights; // it was double

	int Xdim, Ydim;
//	int my_mod_start;  //used for vectorized sort

	// Note that most GPU data is float instead of double
	// d_ means GPU data, h_ means host (uP) data 
	cudaDeviceProp GPUprop;
        int  GPUdevice; 
	int* d_columnH, *d_rowH, *d_aSingleColumnH, *d_aSingleRowH; //matrix
	int* d_columnHt; //, *d_rowH, *d_aSingleColumnH, *d_aSingleRowH; //matrix
	float* d_smoothColumnH, * d_smoothRowH, * d_sumObs; //matrix
	float* d_smoothT; // temporary smooth matrix
	// d_Iframes and h_Iframes contains the same info
	// d_Iframes is used internally by the GPU to access the images
	// h_Iframes is used if a kernel is using a single image
	float* d_ICorrection; 		//Image
	int* d_IframeIdeal; 		//Image
	int* d_IframeTransformed; 	//Image
//	int* d_listOfWidths;		//vector
	float** d_listOfWeights;	//vector of pointers
	float* d_listOfSigmas;		//vector
	float* d_weights;                 // vector
	int** Iframe_vecGPU;		// vector of images
	std::vector<float *> listOfWeightsGPU; // it was double
	clock_t start, end;
	dim3 block, thread;
	size_t sz_imgFL, sz_imgINT;
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
	Image<double> Iframe;
	FileName fnFrame;
	mdIn.getValue(MDL_IMAGE,fnFrame,mdIn.firstObject());
	Iframe.read(fnFrame);
	Xdim=XSIZE(Iframe());
	Ydim=YSIZE(Iframe());
	sz_imgFL=Xdim*Ydim*sizeof(float);
	sz_imgINT=Xdim*Ydim*sizeof(int);

	columnH.initZeros(Ydim,Xdim);
	rowH.initZeros(Ydim,Xdim);
	ICorrection().resizeNoCopy(Ydim,Xdim);
	ICorrection().initConstant(1);
	sumObs.initZeros(Ydim,Xdim);

	// Send ICorrection to GPU
	gpuErrchk(cudaMalloc(&d_ICorrection, sz_imgFL));
	gpuErrchk(cudaMemcpy(d_ICorrection, &ICorrection(0,0), sz_imgFL, cudaMemcpyHostToDevice));        
		
	// Allocate and intialize GPU data	
	gpuErrchk(cudaMalloc( &d_sumObs, sz_imgFL)); // no need to initialized
	gpuErrchk(cudaMalloc( &d_rowH, sz_imgINT)); // no need to initialized
	gpuErrchk(cudaMalloc( &d_columnH, sz_imgINT)); // no need to initialized
	gpuErrchk(cudaMalloc( &d_columnHt, sz_imgINT)); // temporal vector, no need to initialized
	gpuErrchk(cudaMalloc( &d_smoothColumnH, sz_imgFL));
	gpuErrchk(cudaMalloc( &d_smoothT, sz_imgFL));
 	block.x=block.y=block.z=1;
 	thread.x=thread.y=thread.z=1;
 
	MultidimArray<float> sumObsFL; // it was double
	sumObsFL.initZeros(Ydim,Xdim);
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		mdIn.getValue(MDL_IMAGE,fnFrame,__iter.objId);
		Iframe.read(fnFrame);
		sumObs+=Iframe();
	}
	sumObs*=2;
	
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(sumObsFL)
		DIRECT_A2D_ELEM(sumObsFL,i,j)=(float)(DIRECT_A2D_ELEM(sumObs,i,j));
	gpuErrchk(cudaMemcpy(d_sumObs, &sumObsFL(0,0), sz_imgFL, cudaMemcpyHostToDevice));
	
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
		std::cout <<" i=" << i  <<"; p_weights=" << d_weights << std::endl;
		listOfWeightsGPU.push_back(d_weights);
	}
//	gpuErrchk(cudaMemcpy(d_listOfWeightGPU, listOfWeights, sz_listW, cudaMemcpyHostToDevice));
	
	
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
	std::cout << "Total global mem = " << GPUprop. totalGlobalMem << " bytes" << std::endl;

	
	produceSideInfo();

	FileName fnFrame;
	Image<int> Iframe;
	MultidimArray<int> IframeTransformed, IframeIdeal;
	//MultidimArray<double> sumIdeal;
	//MultidimArray<double> &mICorrection=ICorrection();
	MultidimArray<float> sumIdeal;
	MultidimArray<float> &mICorrection=ICorrection();


	// Copy all int images onto GPU
	int im=0;
	size_t sz_Iframe_vec=sizeof(int*)*mdIn.size();
	size_t px_img=Xdim*Ydim;
	//size_t sz_imgINT=px_img*sizeof(int);

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
		sumIdeal.initZeros(Ydim,Xdim);
		im=0;
		FOR_ALL_OBJECTS_IN_METADATA(mdIn)
		{
			mdIn.getValue(MDL_IMAGE,fnFrame,__iter.objId);
			std::cout << "   Frame " << fnFrame << std::endl;
			Iframe.read(fnFrame);
			IframeIdeal = Iframe();
		
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(IframeIdeal)
				DIRECT_A2D_ELEM(IframeIdeal,i,j)=(int)(DIRECT_A2D_ELEM(IframeIdeal,i,j)*DIRECT_A2D_ELEM(mICorrection,i,j));
			//GPU kernel
		        gpuErrchk(cudaMemcpy(d_IframeIdeal, Iframe_vecGPU[im], sz_imgINT, cudaMemcpyDeviceToDevice));
			dim3 block1(floor((Xdim+(TILE_DIM-1))/TILE_DIM),floor((Ydim+(TILE_DIMH-1))/(TILE_DIMH)),1);
			dim3 thread1( TILE_DIM, TILE_DIMH);
			Kmult<<< block1, thread1 >>>(Iframe_vecGPU[im], d_ICorrection, d_IframeIdeal, Xdim, Ydim);
			cudaThreadSynchronize();
			gpuErrchk(cudaGetLastError());

	
			int* IframeIdeal_tmp=(int*)malloc(sz_imgINT); // remove
			gpuErrchk(cudaMemcpy(IframeIdeal_tmp, d_IframeIdeal, sz_imgINT, cudaMemcpyDeviceToHost));
			if (!isequal(&IframeIdeal(0,0), IframeIdeal_tmp, Xdim*Ydim))	
		 		std::cout << "IframeIdeal is not equal!!!!!!!!!!!!" << std::endl;
			im++;

			computeHistograms(IframeIdeal); //rowH columnH
		        computeHistogramsGPU(d_IframeIdeal); // GPU test	
		
std::cout << "1st CHECK" << std::endl;

			int* colH_tmp=(int*)malloc(sz_imgINT); // remove
std::cout << "2nd CHECK" << std::endl;
			gpuErrchk(cudaMemcpy(colH_tmp, d_columnH, sz_imgINT, cudaMemcpyDeviceToHost));
std::cout << "3rd CHECK" << std::endl;
			if (!isequal(&DIRECT_A2D_ELEM(columnH,0,0), colH_tmp, Xdim*Ydim))	
		 		std::cout << "columnH is not equal!!!!!!!!!!!!" << std::endl;
			int* rowH_tmp=(int*)malloc(sz_imgINT); // remove
			gpuErrchk(cudaMemcpy(rowH_tmp, d_rowH, sz_imgINT, cudaMemcpyDeviceToHost));
			if (!isequal(&DIRECT_A2D_ELEM(rowH,0,0), rowH_tmp, Xdim*Ydim))	
		 		std::cout << "rowH is not equal!!!!!!!!!!!!" << std::endl;
			
std::cout << "FINISHED" << std::endl;
	
			size_t bestSigmaCol = selectBestSigmaByColumn(IframeIdeal);
			std::cout << "      sigmaCol: " << listOfSigmas[bestSigmaCol] << std::endl;
	
			size_t bestSigmaColGPU = selectBestSigmaByColumnGPU_TEST(IframeIdeal, d_IframeIdeal);
			std::cout << "      sigmaCol: " << listOfSigmas[bestSigmaCol] << std::endl;
			size_t bestSigmaRow = selectBestSigmaByRow(IframeIdeal);
			std::cout << "      sigmaRow: " << listOfSigmas[bestSigmaRow] << std::endl;

			constructSmoothHistogramsByRow(listOfWeights[bestSigmaRow],listOfWidths[bestSigmaRow]);   //smooth
			transformGrayValuesRow(IframeIdeal,IframeTransformed);
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(IframeTransformed)
				DIRECT_A2D_ELEM(sumIdeal,i,j)+=DIRECT_A2D_ELEM(IframeTransformed,i,j);
			constructSmoothHistogramsByColumn(listOfWeights[bestSigmaCol],listOfWidths[bestSigmaCol]);
			transformGrayValuesColumn(IframeIdeal,IframeTransformed);
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(IframeTransformed)
			DIRECT_A2D_ELEM(sumIdeal,i,j)+=DIRECT_A2D_ELEM(IframeTransformed,i,j);
		}
		//***
		// kernel with partial average
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(mICorrection)
		{
			double den=DIRECT_A2D_ELEM(sumObs,i,j);
			if (fabs(den)<1e-6)
				DIRECT_A2D_ELEM(mICorrection,i,j)=1.0;
			else
				DIRECT_A2D_ELEM(mICorrection,i,j)=DIRECT_A2D_ELEM(sumIdeal,i,j)/den;
		}
		mICorrection/=mICorrection.computeAvg();
		//****

	}
	
//***********
	for (int n=0; n<Niter; n++)
	{
		std::cout << "Iteration " << n << std::endl;
		sumIdeal.initZeros(Ydim,Xdim);
		FOR_ALL_OBJECTS_IN_METADATA(mdIn)
		{
			mdIn.getValue(MDL_IMAGE,fnFrame,__iter.objId);
			std::cout << "   Frame " << fnFrame << std::endl;
			Iframe.read(fnFrame);
			IframeIdeal = Iframe();
                     	std::cout << " XSIZE = " << XSIZE(Iframe());  // gacaffe: just' checking, REMOVE
                     	std::cout << " YSIZE = " << YSIZE(Iframe()) << std::endl;  // gacaffe: just' checking, REMOVE

			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(IframeIdeal)
				DIRECT_A2D_ELEM(IframeIdeal,i,j)=(int)(DIRECT_A2D_ELEM(IframeIdeal,i,j)*DIRECT_A2D_ELEM(mICorrection,i,j));
			computeHistograms(IframeIdeal);

			size_t bestSigmaCol = selectBestSigmaByColumn(IframeIdeal);
			std::cout << "      sigmaCol: " << listOfSigmas[bestSigmaCol] << std::endl;
			size_t bestSigmaRow = selectBestSigmaByRow(IframeIdeal);
			std::cout << "      sigmaRow: " << listOfSigmas[bestSigmaRow] << std::endl;

			constructSmoothHistogramsByRow(listOfWeights[bestSigmaRow],listOfWidths[bestSigmaRow]);
			transformGrayValuesRow(IframeIdeal,IframeTransformed);
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(IframeTransformed)
				DIRECT_A2D_ELEM(sumIdeal,i,j)+=DIRECT_A2D_ELEM(IframeTransformed,i,j);
			constructSmoothHistogramsByColumn(listOfWeights[bestSigmaCol],listOfWidths[bestSigmaCol]);
			transformGrayValuesColumn(IframeIdeal,IframeTransformed);
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(IframeTransformed)
			DIRECT_A2D_ELEM(sumIdeal,i,j)+=DIRECT_A2D_ELEM(IframeTransformed,i,j);
		}

		//***
		// kernel with partial average
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(mICorrection)
		{
			double den=DIRECT_A2D_ELEM(sumObs,i,j);
			if (fabs(den)<1e-6)
				DIRECT_A2D_ELEM(mICorrection,i,j)=1.0;
			else
				DIRECT_A2D_ELEM(mICorrection,i,j)=DIRECT_A2D_ELEM(sumIdeal,i,j)/den;
		}
		mICorrection/=mICorrection.computeAvg();
		//****

	}
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
	cudaFree(d_rowH);
	cudaFree(d_columnH);
	cudaFree(d_columnHt);
	cudaFree(d_smoothColumnH);
	cudaFree(d_IframeIdeal);
	cudaFree(d_IframeTransformed);
	free(Iframe_vecGPU);

}

void ProgMovieEstimateGainGPU::computeHistograms(const MultidimArray<int> &Iframe)
{
		int* auxElemC=new int[Ydim];
		int* auxElemR=new int[Xdim];


		
		start = clock();

		for(size_t j=0; j<XSIZE(columnH); j++)
		{
			for(size_t i=0; i<Ydim; i++)
				auxElemC[i]=A2D_ELEM(Iframe,i,j);
			std::sort(auxElemC, auxElemC+Ydim);
			for(size_t i=0; i<Ydim; i++)
				A2D_ELEM(columnH,i,j)=auxElemC[i];
		}
		delete auxElemC;

		end = clock();
		std::cout << "colH " <<  (float)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl;

		start=clock();
		
		for(size_t i=0; i<YSIZE(rowH); i++)
		{
			memcpy(auxElemR,&A2D_ELEM(Iframe,i,0),Xdim*sizeof(int));
			std::sort(auxElemR, auxElemR+Xdim);
			memcpy(&A2D_ELEM(rowH,i,0),auxElemR,Xdim*sizeof(int));
		}
		delete auxElemR;


		end = clock();
		std::cout << "rowH " <<  (float)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl;
}

// compute colH and rowH of image in GPU (Iframe)
// it's always applied to d_IframeIdeal
void ProgMovieEstimateGainGPU::computeHistogramsGPU(const int* d_IframeIdeal)
{
		//So here we order the rows of the transposed columnH matrix
/*		clock_t	start = clock();
	
	
		thrust::device_ptr<int> th_d_columnH = thrust::device_pointer_cast(d_columnH);
		thrust::host_vector<int> h_segmentsC(Xdim*Ydim);
		my_mod_start = 0;
		no_el = Ydim;
		thrust::generate(h_segmentsC.begin(), h_segmentsC.end(), my_modR);
		thrust::device_vector<int> d_segmentsC= h_segmentsC;
		thrust::stable_sort_by_key(th_d_columnH, th_d_columnH+Xdim*Ydim, d_segmentsC.begin());
		thrust::stable_sort_by_key(d_segmentsC.begin(), d_segmentsC.end(), th_d_columnH);
		cudaDeviceSynchronize();
		clock_t	end = clock();
		std::cout << "(VECTORIZED) GPU colH " <<  (float)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl;
*/

		clock_t start = clock();
		gpuErrchk(cudaMemcpy(d_columnHt, d_IframeIdeal, sz_imgINT, cudaMemcpyDeviceToDevice));	
		thrust::device_ptr<int> th_d_colHt = thrust::device_pointer_cast(d_columnHt);
		thrust::host_vector<int> h_segmentsC(Xdim*Ydim);
		my_mod_start = 0;
		no_el = Xdim;
		thrust::generate(h_segmentsC.begin(), h_segmentsC.end(), my_modC);
		thrust::device_vector<int> d_segmentsC = h_segmentsC;
		thrust::stable_sort_by_key(th_d_colHt, th_d_colHt+Xdim*Ydim, d_segmentsC.begin());
		thrust::stable_sort_by_key(d_segmentsC.begin(), d_segmentsC.end(), th_d_colHt);
		cudaDeviceSynchronize();


		// transpose output (note that the input has Ydim columns!!!)
		TRANSPOSE_GPU(d_columnH, d_columnHt, Ydim, Xdim);
			 	
		end = clock();
		std::cout << "(VECTORIZED and transposed) GPU colH " <<  (float)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl;


// VECTORIZED SORT: much better than CPU
		start = clock();
		gpuErrchk(cudaMemcpy(d_rowH, d_IframeIdeal, sz_imgINT, cudaMemcpyDeviceToDevice));	
		thrust::device_ptr<int> th_d_rowH = thrust::device_pointer_cast(d_rowH);
		thrust::host_vector<int> h_segmentsR(Xdim*Ydim);
		my_mod_start = 0;
		no_el = Xdim;
		thrust::generate(h_segmentsR.begin(), h_segmentsR.end(), my_modR);
		thrust::device_vector<int> d_segmentsR= h_segmentsR;
		thrust::stable_sort_by_key(th_d_rowH, th_d_rowH+Xdim*Ydim, d_segmentsR.begin());
		thrust::stable_sort_by_key(d_segmentsR.begin(), d_segmentsR.end(), th_d_rowH);
		cudaDeviceSynchronize();
		end = clock();
		std::cout << "(VECTORIZED) GPU rowH " <<  (float)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl;


}


void ProgMovieEstimateGainGPU::constructSmoothHistogramsByColumn(const float *listOfWeights, int width)
{

	smoothColumnH.initZeros(columnH);

	for (size_t j=0; j<XSIZE(columnH); ++j)
	{
		double sumWeightsC = 0;
		for(int k = -width; k<=width; ++k)
		{
			if (j+k<0 || j+k>=XSIZE(columnH))
				continue;
			float actualWeightC = listOfWeights[abs(k)];
			sumWeightsC += actualWeightC;
			for (size_t i=0; i<Ydim; ++i)
				DIRECT_A2D_ELEM(smoothColumnH,i,j) += actualWeightC * DIRECT_A2D_ELEM(columnH,i,j+k);
		}

		double iSumWeightsC=1/sumWeightsC;
		for (size_t i=0; i<Ydim; ++i)
			DIRECT_A2D_ELEM(smoothColumnH,i,j) *= iSumWeightsC;
	}

	if (singleReference)
	{
		// Compute the average of all column histograms
		for (size_t i=0; i<Ydim; ++i)
			for (size_t j=1; j<Xdim; ++j)
				DIRECT_A2D_ELEM(smoothColumnH,i,0)+=DIRECT_A2D_ELEM(smoothColumnH,i,j);

		double iXdim=1.0/Xdim;
		for (size_t i=0; i<Ydim; ++i)
		{
			DIRECT_A2D_ELEM(smoothColumnH,i,0)*=iXdim;
			double aux=DIRECT_A2D_ELEM(smoothColumnH,i,0);
			for (size_t j=1; j<Xdim; ++j)
				DIRECT_A2D_ELEM(smoothColumnH,i,j)=aux;
		}
	}
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
		TRANSPOSE_GPU(d_smoothT, d_smoothColumnH, Xdim, Ydim);

		// Average all rows and replicate in all rows
		dim3 block4(floor((Xdim+(TILE_DIM2-1))/TILE_DIM2),1,1);
		dim3 thread4( TILE_DIM2, 1);
		Ksmooth2<<< block4, thread4 >>>(d_smoothT, Ydim, Xdim);
		cudaThreadSynchronize();
		gpuErrchk(cudaGetLastError());

		// transpose back smoothColumnH
		TRANSPOSE_GPU(d_smoothColumnH, d_smoothT, Ydim, Xdim);

	}

}

void ProgMovieEstimateGainGPU::constructSmoothHistogramsByRow(const float *listOfWeights, int width)
{
	smoothRowH.initZeros(rowH);
	for(size_t i = 0; i<YSIZE(rowH); ++i)
	{
		double sumWeightsR = 0;
		for(int k = -width; k<=width; ++k)
		{
			if (i+k<0 || i+k>=YSIZE(rowH))
				continue;
			double actualWeightR = listOfWeights[abs(k)];
			sumWeightsR += actualWeightR;
			for (size_t j=0; j< Xdim; ++j)
				DIRECT_A2D_ELEM(smoothRowH,i,j) += actualWeightR * DIRECT_A2D_ELEM(rowH,i+k,j);
		}
		double iSumWeightsR=1/sumWeightsR;
		for (size_t j=0; j<Xdim; ++j)
			DIRECT_A2D_ELEM(smoothRowH,i,j) *= iSumWeightsR;
	}

	if (singleReference)
	{
		// Compute the average of all row histograms
		for (size_t j=0; j<Xdim; ++j)
			for (size_t i=1; i<Ydim; ++i)
				DIRECT_A2D_ELEM(smoothRowH,0,j)+=DIRECT_A2D_ELEM(smoothRowH,i,j);

		double iYdim=1.0/Ydim;
		for (size_t j=0; j<Xdim; ++j)
		{
			DIRECT_A2D_ELEM(smoothRowH,0,j)*=iYdim;
			double aux=DIRECT_A2D_ELEM(smoothRowH,0,j);
			for (size_t i=1; i<Ydim; ++i)
				DIRECT_A2D_ELEM(smoothRowH,i,j)=aux;
		}
	}


}


void ProgMovieEstimateGainGPU::transformGrayValuesColumn(const MultidimArray<int> &Iframe, MultidimArray<int> &IframeTransformedColumn)
{
	IframeTransformedColumn.initZeros(Ydim,Xdim);
	aSingleColumnH.resizeNoCopy(Ydim);
	int *aSingleColumnH0=&DIRECT_A1D_ELEM(aSingleColumnH,0);
	int *aSingleColumnHF=(&DIRECT_A1D_ELEM(aSingleColumnH,Ydim-1))+1;


	for (size_t j=0; j<Xdim; ++j)
	{
		for (size_t i=0; i<Ydim; ++i)
			DIRECT_A1D_ELEM(aSingleColumnH,i)=DIRECT_A2D_ELEM(columnH,i,j);

		for (size_t i=0; i<Ydim; ++i)
		{
			int pixval=DIRECT_A2D_ELEM(Iframe,i,j);
			int *pixvalPtr=std::upper_bound(aSingleColumnH0,aSingleColumnHF,pixval);
			pixvalPtr-=1;
			int idx=pixvalPtr-aSingleColumnH0;
			DIRECT_A2D_ELEM(IframeTransformedColumn,i,j)+=(int)DIRECT_A2D_ELEM(smoothColumnH,idx,j);
		}
	}


}

void ProgMovieEstimateGainGPU::transformGrayValuesColumnGPU(const int* d_Iframe, int* d_IframeTransformedColumn)
{
	//IframeTransformedColumn.initZeros(Ydim,Xdim);

	FILL_GPU(d_IframeTransformedColumn, 0, Xdim, Ydim);

	//aSingleColumnH.resizeNoCopy(Ydim);
	//int *aSingleColumnH0=&DIRECT_A1D_ELEM(aSingleColumnH,0);
	//int *aSingleColumnHF=(&DIRECT_A1D_ELEM(aSingleColumnH,Ydim-1))+1;

	dim3 block1(floor((Xdim+(TILE_DIM2-1))/TILE_DIM2),1,1);
	dim3 thread1( TILE_DIM2, 1);
		
      	KtransformGray<<<block1,thread1>>>(d_Iframe, d_columnH, d_IframeTransformedColumn, d_smoothColumnH, Xdim, Ydim);


}
void ProgMovieEstimateGainGPU::transformGrayValuesRow(const MultidimArray<int> &Iframe, MultidimArray<int> &IframeTransformedRow)
{
	IframeTransformedRow.initZeros(Ydim,Xdim);
		aSingleRowH.resizeNoCopy(Xdim);
		int *aSingleRowH0=&DIRECT_A1D_ELEM(aSingleRowH,0);
		int *aSingleRowHF=(&DIRECT_A1D_ELEM(aSingleRowH,Xdim-1))+1;

	for (size_t i=0; i<Ydim; ++i)
	{
		memcpy(aSingleRowH0,&DIRECT_A2D_ELEM(rowH,i,0),Xdim*sizeof(int));

			for (size_t j=0; j<Xdim; ++j)
			{
				int pixvalR=DIRECT_A2D_ELEM(Iframe,i,j);
				int *pixvalPtrR=std::upper_bound(aSingleRowH0,aSingleRowHF,pixvalR);
				pixvalPtrR-=1;
				int idxR=pixvalPtrR-aSingleRowH0;
				DIRECT_A2D_ELEM(IframeTransformedRow,i,j)+=(int)DIRECT_A2D_ELEM(smoothRowH,i,idxR);
			}
	}

}


double computeTVColumns(MultidimArray<int> &I)
{
	double retvalC=0;
	for (size_t i=0; i<YSIZE(I); ++i)
	    for (size_t j=0; j<XSIZE(I)-1; ++j)
	    	retvalC+=std::abs(DIRECT_A2D_ELEM(I,i,j)-DIRECT_A2D_ELEM(I,i,j+1));

	return retvalC/((XSIZE(I)-1)*YSIZE(I));
}

double computeTVRows(MultidimArray<int> &I)
{
	double retvalR=0;
	for (size_t i=0; i<YSIZE(I)-1; ++i)
	    for (size_t j=0; j<XSIZE(I); ++j)
	    	retvalR+=std::abs(DIRECT_A2D_ELEM(I,i,j)-DIRECT_A2D_ELEM(I,i+1,j));

	return retvalR/((YSIZE(I)-1)*XSIZE(I));
}


size_t ProgMovieEstimateGainGPU::selectBestSigmaByColumn(const MultidimArray<int> &Iframe)
{
	double bestAvgTV=1e38;
	size_t best_s=0;
	MultidimArray<int> IframeTransformed;

	for(size_t s = 0; s< listOfWeights.size(); ++s)
	{
		constructSmoothHistogramsByColumn(listOfWeights[s],listOfWidths[s]);

		transformGrayValuesColumn(Iframe,IframeTransformed);
		double avgTV=computeTVColumns(IframeTransformed);
		if (avgTV<bestAvgTV)
		{
			bestAvgTV=avgTV;
			best_s=s;
		}
	}
	return best_s;
}


size_t ProgMovieEstimateGainGPU::selectBestSigmaByColumnGPU_TEST(const MultidimArray<int> &Iframe, const int* d_Iframe)
{
	double bestAvgTV=1e38;
	size_t best_s=0;
	MultidimArray<int> IframeTransformed;

	for(size_t s = 0; s< listOfWeights.size(); ++s)
	{
//		constructSmoothHistogramsByColumn(listOfWeights[s],listOfWidths[s]);
		start = clock();
		constructSmoothHistogramsByColumn(listOfWeights[s],listOfWidths[s]);
		end = clock();
		std::cout << "SMOOTH_HBCTEST " <<  (float)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl;

		start = clock();
		constructSmoothHistogramsByColumnGPU(listOfWeightsGPU[s],listOfWidths[s]);
		end = clock();
		std::cout << "SMOOTH_HBCGPU " <<  (float)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl;

		
		//size_t sz_imgFL = sizeof(float)*Xdim*Ydim;
		float* smoothColumnH_tmp=(float*)malloc(sz_imgFL); // remove
		gpuErrchk(cudaMemcpy(smoothColumnH_tmp, d_smoothColumnH, sz_imgFL, cudaMemcpyDeviceToHost));
		
		if (!isequal(&DIRECT_A2D_ELEM(smoothColumnH,0,0), smoothColumnH_tmp, Xdim*Ydim))	
		 	std::cout << "SmoothColumn  is not equal!!!!!!!!!!!!" << std::endl;

		start = clock();
		transformGrayValuesColumn(Iframe,IframeTransformed);
		end = clock();
		std::cout << "TransformColCPU " <<  (float)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl;

		start = clock();
		transformGrayValuesColumnGPU(d_Iframe,d_IframeTransformed);
       		end = clock();
		std::cout << "TransformCoGPU " <<  (float)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl;

        	int* IframeTransformed_tmp=(int*)malloc(sz_imgINT); // remove
		gpuErrchk(cudaMemcpy(IframeTransformed_tmp, d_IframeTransformed, sz_imgINT, cudaMemcpyDeviceToHost));
		
		if (!isequal(&DIRECT_A2D_ELEM(IframeTransformed,0,0), IframeTransformed_tmp, Xdim*Ydim))	
		 	std::cout << "IframeTransformed  is not equal!!!!!!!!!!!!" << std::endl;

 
		double avgTV=computeTVColumns(IframeTransformed);
		if (avgTV<bestAvgTV)
		{
			bestAvgTV=avgTV;
			best_s=s;
		}
	}
	return best_s;
}

size_t ProgMovieEstimateGainGPU::selectBestSigmaByRow(const MultidimArray<int> &Iframe)
{
	double bestAvgTV=1e38;
	size_t best_s=0;
	MultidimArray<int> IframeTransformed;

	for(size_t s = 0; s< listOfWeights.size(); ++s)
	{
		constructSmoothHistogramsByRow(listOfWeights[s],listOfWidths[s]);
		transformGrayValuesRow(Iframe,IframeTransformed);
		double avgTV=computeTVRows(IframeTransformed);
		if (avgTV<bestAvgTV)
		{
			bestAvgTV=avgTV;
			best_s=s;
		}
	}

	return best_s;
}

/*********** REMOVE **************/
void ProgMovieEstimateGainGPU::constructSmoothHistogramsByColumnTEST(const float *listOfWeights, int width)
{

	smoothColumnH.initZeros(columnH);

	for (size_t j=0; j<XSIZE(columnH); ++j)
	{
		double sumWeightsC = 0;
		for(int k = -width; k<=width; ++k)
		{
			if (j+k<0 || j+k>=XSIZE(columnH))
				continue;
			double actualWeightC = listOfWeights[abs(k)];
			sumWeightsC += actualWeightC;
			for (size_t i=0; i<Ydim; ++i)
				DIRECT_A2D_ELEM(smoothColumnH,i,j) += actualWeightC * DIRECT_A2D_ELEM(columnH,i,j+k);
		}

		double iSumWeightsC=1./sumWeightsC; //GCF: it wass 1/sum...
		for (size_t i=0; i<Ydim; ++i)
			DIRECT_A2D_ELEM(smoothColumnH,i,j) *= iSumWeightsC;
	}

//	if (singleReference)
/*	if (0)

	{
		// Compute the average of all column histograms
		for (size_t i=0; i<Ydim; ++i)
			for (size_t j=1; j<Xdim; ++j)
				DIRECT_A2D_ELEM(smoothColumnH,i,0)+=DIRECT_A2D_ELEM(smoothColumnH,i,j);

		double iXdim=1.0/Xdim;
		for (size_t i=0; i<Ydim; ++i)
		{
			DIRECT_A2D_ELEM(smoothColumnH,i,0)*=iXdim;
			double aux=DIRECT_A2D_ELEM(smoothColumnH,i,0);
			for (size_t j=1; j<Xdim; ++j)
				DIRECT_A2D_ELEM(smoothColumnH,i,j)=aux;
		}
	}
*/
}


//*********** REMOVE **************/


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
		std::cout<<h_segR[i] << " ";
	std::cout << std::endl;
	thrust::device_vector<int> d_segR= h_segR;
	thrust::stable_sort_by_key(th_d_t, th_d_t+xd*yd, d_segR.begin());
	thrust::stable_sort_by_key(d_segR.begin(), d_segR.end(), th_d_t);
	cudaDeviceSynchronize();
	std::cout<<"after sorting"<<std::endl;
	gpuErrchk(cudaMemcpy(h_t2, d_t, sz_t, cudaMemcpyDeviceToHost));
	std::cout<<"after copying"<<std::endl;
	isequal(h_t,h_t2, xd*yd, xd*yd);



*/
