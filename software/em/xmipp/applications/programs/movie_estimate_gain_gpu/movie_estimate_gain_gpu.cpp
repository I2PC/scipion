#include <cuda_runtime.h>

#include <data/xmipp_program.h>

#include <time.h> // GCF
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "movie_estimate_gain_kernels.cu"

double computeTVColumns(MultidimArray<int> &I);
double computeTVRows(MultidimArray<int> &I);


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
    void normalizeHistograms();
    void invertHistograms();

    void constructSmoothHistogramsByColumn(const float *listOfWeights, int width);
    void constructSmoothHistogramsByRow(const float *listOfWeights, int width);
    void transformGrayValuesColumn(const MultidimArray<int> &Iframe, MultidimArray<int> &IframeTransformedColumn);
    void transformGrayValuesRow(const MultidimArray<int> &Iframe, MultidimArray<int> &IframeTransformedRow);
    void computeTransformedHistograms(const MultidimArray<int> &Iframe);

    size_t selectBestSigmaByColumn(const MultidimArray<int> &Iframe);
    size_t selectBestSigmaByRow(const MultidimArray<int> &Iframe);



public:
	MetaData mdIn;
	MultidimArray<int> columnH,rowH, aSingleColumnH, aSingleRowH;
	MultidimArray<float> smoothColumnH, smoothRowH, sumObs; // it was double
	Image<float> ICorrection;         // it was double
	std::vector<float> listOfSigmas;  // it was double
	std::vector<float> listOfWidths; // it was double
	std::vector<float *> listOfWeights; // it was double

	int Xdim, Ydim;

	// Note that most GPU data is float instead of double
	// d_ means GPU data, h_ means host (uP) data 
	cudaDeviceProp GPUprop;
        int  GPUdevice; 
	int* d_columnH, d_rowH, d_aSingleColumnH, d_aSingleRowH; //matrix
	float* d_smoothColumnH, * d_smoothRowH, * d_sumObs; //matrix
	// d_Iframes and h_Iframes contains the same info
	// d_Iframes is used internally by the GPU to access the images
	// h_Iframes is used if a kernel is using a single image
	float** d_IframesFL;		// vector of images
	float** h_IframesFL;		// vector of images
	int** d_IframesINT;		// vector of images
	int** h_IframesINT;		// vector of images
	float* d_ICorrection; 		//Image
	int* d_IframeIdeal; 		//Image
	float* d_listOfWidths;		//vector
	float** d_listOfWeights;	//vector
	float* d_listOfSigmas;		//vector

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
    addParamsLine("               :+This assumes that there is no image contamination or carbon holes");
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
//	Image<double> Iframe;
	Image<float> Iframe;	// GPU => float
	FileName fnFrame;
	mdIn.getValue(MDL_IMAGE,fnFrame,mdIn.firstObject());
	Iframe.read(fnFrame);
	Xdim=XSIZE(Iframe());
	Ydim=YSIZE(Iframe());

	columnH.initZeros(Ydim,Xdim);
	rowH.initZeros(Ydim,Xdim);
	ICorrection().resizeNoCopy(Ydim,Xdim);
	ICorrection().initConstant(1);
	sumObs.initZeros(Ydim,Xdim);

	clock_t start=clock();
	// Copy all images to GPU
	int im=0;
	size_t sz_IframesFL=sizeof(float*)*mdIn.size();
	h_IframesFL = (float**)malloc(sz_IframesFL);
	size_t sz_img=Xdim*Ydim*sizeof(float);
 
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		std::cout << "im="<< im <<"; Image size in bytes = " << sz_img << std::endl;
		mdIn.getValue(MDL_IMAGE,fnFrame,__iter.objId);
		Iframe.read(fnFrame);

		gpuErrchk(cudaMalloc(&h_IframesFL[im], sz_img));
		gpuErrchk(cudaMemcpy(h_IframesFL[im++], &Iframe(0,0), sz_img, cudaMemcpyHostToDevice)); 
	}
	// send vector of images to GPU
	gpuErrchk(cudaMalloc(&d_IframesFL, sz_IframesFL));
	gpuErrchk(cudaMemcpy(d_IframesFL, h_IframesFL, sz_IframesFL, cudaMemcpyHostToDevice));
	// Allocate and intialize GPU data	gpuErrchk(cudaMalloc( &d_sumObs, sz_img)); // no need to initialize (it's done in kernel)
	//****
	// Kernel
	std::cout << "CALLING KERNEL" << std::endl;

        dim3 block(floor((Xdim+31)/32),floor((Ydim+17)/16),1);
	dim3 thread( 32, 16);
	sumall<<< block, thread >>>(d_sumObs, d_IframesFL, mdIn.size(), Xdim, Ydim);
 	cudaThreadSynchronize();
	gpuErrchk(cudaGetLastError());

	float* sumObs_tmp=(float*)malloc(sz_img); // remove
	
	std::cout << "READING RESULTS " << sz_img << std::endl;
	cudaMemcpy(sumObs_tmp, d_sumObs, sz_img, cudaMemcpyDeviceToHost);
	std::cout << "RESULTS READ" << std::endl;
	clock_t end=clock();
	std::cout << "SEND+KERNEL+SEND= " <<  (float)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl;
	
	start=clock();
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		mdIn.getValue(MDL_IMAGE,fnFrame,__iter.objId);
		Iframe.read(fnFrame);
		sumObs+=Iframe();
	}
	sumObs*=2;
	end=clock();	
	std::cout << "SEQUENTIAL= " <<  (float)(end - start) / CLOCKS_PER_SEC << "secs" << std::endl;
	
	// check if sumObs and sumObs_tmp are the same
	for (int i=0; i<Xdim*Ydim; i++){
		if (DIRECT_A1D_ELEM(sumObs,i)!=sumObs_tmp[i]){
			std::cout << "ERROR: sumObs"<< std::endl;
			break;
		}	
	}	
	std::cout << "sumObs[0]=" << DIRECT_A2D_ELEM(sumObs,0,0) << " -> d_sumObs[0]=" << sumObs_tmp[0] << std::endl;

	// free stuff // REMOVE
	free(sumObs_tmp);
	cudaFree(d_sumObs);
/*	for (int i=0; i< mdIn.size(); i++){
		cudaFree(d_Iframes[i]);
		free(h_Iframes[i]);
	}
	cudaFree(d_Iframes);
	free(h_Iframes);
*/ 
	//****

	// Initialize sigma values
	for (float sigma=0; sigma<=maxSigma; sigma+=sigmaStep)
		listOfSigmas.push_back(sigma);
	
	for (size_t i=0; i<listOfSigmas.size(); ++i)
	{
		int jmax=ceil(3*listOfSigmas[i]);
		listOfWidths.push_back(jmax);
		float *weights=new float[jmax];    //it was double
		float K=-0.5/(listOfSigmas[i]*listOfSigmas[i]); //it was double
//********* MODIFED BY GCF
		//for (int j=1; j<=jmax; ++j)
		//	weights[j-1]=exp(K*j*j);
		for (int j=2; j<=jmax; ++j)
			weights[j-1]=exp(K*j*j);
		weights[0]=1.;			// this is always 1. when used 
//***********************/Weights
		listOfWeights.push_back(weights);
	}
	// Copy to GPU
	size_t sz_list = sizeof(float)*listOfSigmas.size();
	gpuErrchk(cudaMalloc(&d_listOfSigmas, sz_list));
	gpuErrchk(cudaMemcpy(d_listOfSigmas, &listOfSigmas[0], sz_list, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMalloc(&d_listOfWeights, sz_list));
	
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
	size_t sz_IframesINT=sizeof(int*)*mdIn.size();
	h_IframesINT = (int**)malloc(sz_IframesINT);
	size_t px_img=Xdim*Ydim;
	size_t sz_img=px_img*sizeof(int);

 
	FOR_ALL_OBJECTS_IN_METADATA(mdIn)
	{
		std::cout << "im="<< im <<"; Image size in bytes = " << sz_img << std::endl;
		mdIn.getValue(MDL_IMAGE,fnFrame,__iter.objId);
		Iframe.read(fnFrame);

		gpuErrchk(cudaMalloc(&h_IframesINT[im], sz_img));
		gpuErrchk(cudaMemcpy(h_IframesINT[im++], &Iframe(0,0), sz_img, cudaMemcpyHostToDevice)); 
	}
	// send vector of images to GPU
	gpuErrchk(cudaMalloc(&d_IframesINT, sz_IframesINT));
	gpuErrchk(cudaMemcpy(d_IframesINT, h_IframesINT, sz_IframesINT, cudaMemcpyHostToDevice));
	// Malloc several data
	gpuErrchk(cudaMalloc(&d_IframeIdeal, sz_img));
	
	
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
			thrust::device_ptr<int> thptr_d_IframeIdeal = thrust::device_pointer_cast(d_IframeIdeal);
			thrust::device_ptr<int> thptr_d_Iframe = thrust::device_pointer_cast(h_IframesINT[im]);
			thrust::copy(thptr_d_Iframe, thptr_d_Iframe+px_img, thptr_d_IframeIdeal);
			thrust::copy(thptr_d_Iframe, thptr_d_Iframe+px_img, thptr_d_IframeIdeal);
			thrust::transform(thptr_d_IframeIdeal, thptr_d_IframeIdeal+px_img, thptr_d_ICorrection, thrust::multiplies());
			im++;

			computeHistograms(IframeIdeal); //rowH columnH

			size_t bestSigmaCol = selectBestSigmaByColumn(IframeIdeal);
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

#ifdef NEVER_DEFINED
	ICorrection.write(fnCorr);
	Image<double> save;
	typeCast(sumIdeal,save());
	save.write("PPPSumIdeal.xmp");
	typeCast(sumObs,save());
	save.write("PPPSumObs.xmp");
	//std::cout << "Press any key\n";
	//char c; std::cin >> c;
#endif
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

#ifdef NEVER_DEFINED
	ICorrection.write(fnCorr);
	Image<double> save;
	typeCast(sumIdeal,save());
	save.write("PPPSumIdeal.xmp");
	typeCast(sumObs,save());
	save.write("PPPSumObs.xmp");
	//std::cout << "Press any key\n";
	//char c; std::cin >> c;
#endif
	}
	ICorrection.write(fnRoot+"_correction.xmp");
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(mICorrection)
		if (DIRECT_A2D_ELEM(mICorrection,i,j)>1e-5)
			DIRECT_A2D_ELEM(mICorrection,i,j)=1.0/DIRECT_A2D_ELEM(mICorrection,i,j);
		else
			DIRECT_A2D_ELEM(mICorrection,i,j)=1;
	ICorrection.write(fnRoot+"_gain.xmp");
}

void ProgMovieEstimateGainGPU::computeHistograms(const MultidimArray<int> &Iframe)
{
		int* auxElemC=new int[Ydim];
		int* auxElemR=new int[Xdim];


		
		clock_t start = clock();

		for(size_t j=0; j<XSIZE(columnH); j++)
		{
			for(size_t i=0; i<Ydim; i++)
				auxElemC[i]=A2D_ELEM(Iframe,i,j);
			std::sort(auxElemC, auxElemC+Ydim);
			for(size_t i=0; i<Ydim; i++)
				A2D_ELEM(columnH,i,j)=auxElemC[i];
		}
		delete auxElemC;

		clock_t end = clock();
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
#ifdef NEVER_DEFINED
	Image<double> save;
	typeCast(columnH,save());
	save.write("PPPcolumnH.xmp");
	typeCast(rowH,save());
	save.write("PPProwH.xmp");
#endif
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
			//**** MODIFED BY GCF
			//double actualWeightC = k==0? 1:listOfWeights[abs(k)];
			double actualWeightC = listOfWeights[abs(k)];
			//********	
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
#ifdef NEVER_DEFINED
	Image<double> save;
	typeCast(smoothColumnH,save());
	save.write("PPPsmoothColumnH.xmp");
#endif
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
			//**** MODIFED BY GCF
			//double actualWeightC = k==0? 1:listOfWeights[abs(k)];
			double actualWeightR = listOfWeights[abs(k)];
			//********	
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

#ifdef NEVER_DEFINED
	Image<double> save;
	typeCast(smoothRowH,save());
	save.write("PPPsmoothRowH.xmp");
#endif
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


#ifdef NEVER_DEFINED
	Image<double> save;
	typeCast(IframeTransformedColumn,save());
	save.write("PPPIframeTransformedColumn.xmp");
	typeCast(Iframe,save());
	save.write("PPPIframe.xmp");
	std::cout << "Press any key to continue\n";
	char c; std::cin >> c;
#endif
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

#ifdef NEVER_DEFINED
	Image<double> save;
	typeCast(IframeTransformedRow,save());
	save.write("PPPIframeTransformedRow.xmp");
	typeCast(Iframe,save());
	save.write("PPPIframe.xmp");
	std::cout << "Press any key to continue\n";
	char c; std::cin >> c;
#endif

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

RUN_XMIPP_PROGRAM(ProgMovieEstimateGainGPU)
