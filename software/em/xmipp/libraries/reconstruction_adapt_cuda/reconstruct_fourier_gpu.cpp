/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
 *              Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Jose Roman Bilbao-Castro (jrbcast@ace.ual.es)
 *              Vahid Abrishami (vabrishami@cnb.csic.es)
 *              David Strelak (davidstrelak@gmail.com)
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

#include "reconstruct_fourier_gpu.h"

// Define params
void ProgRecFourierGPU::defineParams()
{
    //usage
    addUsageLine("Generate 3D reconstructions from projections using direct Fourier interpolation with arbitrary geometry.");
    addUsageLine("Kaisser-windows are used for interpolation in Fourier space.");
    //params
    addParamsLine("   -i <md_file>                : Metadata file with input projections");
    addParamsLine("  [-o <volume_file=\"rec_fourier.vol\">]  : Filename for output volume");
    addParamsLine("  [--sym <symfile=c1>]              : Enforce symmetry in projections");
    addParamsLine("  [--padding <proj=2.0> <vol=2.0>]  : Padding used for projections and volume");
    addParamsLine("  [--max_resolution <p=0.5>]     : Max resolution (Nyquist=0.5)");
    addParamsLine("  [--weight]                     : Use weights stored in the image metadata");
    addParamsLine("  [--blob <radius=1.9> <order=0> <alpha=15>] : Blob parameters");
    addParamsLine("                                 : radius in pixels, order of Bessel function in blob and parameter alpha");
    addParamsLine("  [--fast]                       : Do the blobing at the end of the computation.");
    addParamsLine("                                 : Gives slightly different results, but is faster.");
    addParamsLine("  [--fftOnGPU]                   : Perform the FFT conversion of the input images on GPU.");
    addParamsLine("                                 : Requires more memory on GPU (see also --bufferSize).");
    addParamsLine("  [--useCTF]                     : Use CTF information if present");
    addParamsLine("  [--sampling <Ts=1>]            : sampling rate of the input images in Angstroms/pixel");
    addParamsLine("                                 : It is only used when correcting for the CTF");
    addParamsLine("  [--phaseFlipped]               : Give this flag if images have been already phase flipped");
    addParamsLine("  [--minCTF <ctf=0.01>]          : Minimum value of the CTF that will be inverted");
    addParamsLine("                                 : CTF values (in absolute value) below this one will not be corrected");
    addParamsLine("  [--bufferSize <size=25>]        : Number of projection loaded in memory (will be actually 2x as much.");
    addParamsLine("                                 : This will require up to 4*size*projSize*projSize*16B, e.g.");
    addParamsLine("                                 : 100MB for projection of 256x256 or 400MB for projection of 512x512");
    addExampleLine("For reconstruct enforcing i3 symmetry and using stored weights:", false);
    addExampleLine("   xmipp_reconstruct_fourier  -i reconstruction.sel --sym i3 --weight");
}

// Read arguments ==========================================================
void ProgRecFourierGPU::readParams()
{
    fn_in = getParam("-i");
    fn_out = getParam("-o");
    fn_sym = getParam("--sym");
    do_weights = checkParam("--weight");
    padding_factor_proj = getDoubleParam("--padding", 0);
    padding_factor_vol = getDoubleParam("--padding", 1);
    blob.radius   = getDoubleParam("--blob", 0);
    blob.order    = getIntParam("--blob", 1);
    blob.alpha    = getDoubleParam("--blob", 2);
    useFast		  = checkParam("--fast");
    fftOnGPU	  = checkParam("--fftOnGPU");
    maxResolution = getDoubleParam("--max_resolution");
    useCTF = checkParam("--useCTF");
    isPhaseFlipped = checkParam("--phaseFlipped");
    minCTF = getDoubleParam("--minCTF");
    if (useCTF)
        iTs = 1 / getDoubleParam("--sampling");
    bufferSize = getIntParam("--bufferSize");
}

// Show ====================================================================
void ProgRecFourierGPU::show()
{
    if (verbose > 0)
    {
        std::cout << " =====================================================================" << std::endl;
        std::cout << " Direct 3D reconstruction method using Kaiser windows as interpolators" << std::endl;
        std::cout << " =====================================================================" << std::endl;
        std::cout << " Input selfile             : "  << fn_in << std::endl;
        std::cout << " padding_factor_proj       : "  << padding_factor_proj << std::endl;
        std::cout << " padding_factor_vol        : "  << padding_factor_vol << std::endl;
        std::cout << " Output volume             : "  << fn_out << std::endl;
        if (fn_sym != "")
            std::cout << " Symmetry file for projections : "  << fn_sym << std::endl;
        if (do_weights)
            std::cout << " Use weights stored in the image headers or doc file" << std::endl;
        else
            std::cout << " Do NOT use weights" << std::endl;
        if (useCTF)
            std::cout << "Using CTF information" << std::endl
            << "Sampling rate: " << 1 / iTs << std::endl
            << "Phase flipped: " << isPhaseFlipped << std::endl
            << "Minimum CTF: " << minCTF << std::endl;
        std::cout << "\n Interpolation Function"
        << "\n   blrad                 : "  << blob.radius
        << "\n   blord                 : "  << blob.order
        << "\n   blalpha               : "  << blob.alpha
        << "\n max_resolution          : "  << maxResolution
        << "\n -----------------------------------------------------------------" << std::endl;
    }
}

// Main routine ------------------------------------------------------------
void ProgRecFourierGPU::run()
{
    show();
    produceSideinfo();
    initProgress();
    //Computing interpolated volume
    processImages(0, SF.size() - 1);
    // Get data from GPU
    getGPUData();
    // remove complex conjugate of the intermediate result
    mirrorAndCropTempSpaces();
    //Saving the volume
    finishComputations(fn_out);
}

void ProgRecFourierGPU::createWorkThread(int gpuStream, int startIndex, int endIndex, RecFourierWorkThread& thread) {
	thread.parent = this;
	thread.startImageIndex = startIndex;
	thread.endImageIndex = endIndex;
	thread.gpuStream = gpuStream;
	pthread_create( &thread.id , NULL, threadRoutine, (void *)(&thread) );
}


void ProgRecFourierGPU::produceSideinfo()
{
    // Translate the maximum resolution to digital frequency
    // maxResolution=sampling_rate/maxResolution;
    maxResolutionSqr=maxResolution*maxResolution;

    // Read the input images
    SF.read(fn_in);
    SF.removeDisabled();
    SF.getDatabase()->activateThreadMuting();

    // Ask for memory for the output volume and its Fourier transform
    size_t objId = SF.firstObject();
    FileName fnImg;
    SF.getValue(MDL_IMAGE,fnImg,objId);
    Image<double> I;
    I.read(fnImg, HEADER);
    int Ydim=YSIZE(I());
    int Xdim=XSIZE(I());
    if (Ydim!=Xdim)
        REPORT_ERROR(ERR_MULTIDIM_SIZE,"This algorithm only works for squared images");
    imgSize=Xdim;
    paddedImgSize = Xdim*padding_factor_vol;
	size_t conserveRows = (size_t) ceil((double) paddedImgSize * maxResolution * 2.0);
	conserveRows = (size_t) ceil((double) conserveRows / 2.0);
	maxVolumeIndexX = maxVolumeIndexYZ = 2 * conserveRows;

    // Build a table of blob values
    Fourier_blob_table.resize(BLOB_TABLE_SIZE_SQRT);

    struct blobtype blobFourier,blobnormalized;
    blobFourier=blob;
    blobFourier.radius/=(padding_factor_vol*Xdim);
    blobnormalized=blob;
    blobnormalized.radius/=((double)padding_factor_proj/padding_factor_vol);
    double deltaSqrt     = (blob.radius*blob.radius) /(BLOB_TABLE_SIZE_SQRT-1);
    double deltaFourier  = (sqrt(3.)*Xdim/2.)/(BLOB_TABLE_SIZE_SQRT-1);

    // The interpolation kernel must integrate to 1
    iw0 = 1.0 / blob_Fourier_val(0.0, blobnormalized);
    double padXdim3 = padding_factor_vol * Xdim;
    padXdim3 = padXdim3 * padXdim3 * padXdim3;
    double blobTableSize = blob.radius*sqrt(1./ (BLOB_TABLE_SIZE_SQRT-1));
    for (int i = 0; i < BLOB_TABLE_SIZE_SQRT; i++)
    {
		//use a r*r sample instead of r
		//DIRECT_VEC_ELEM(blob_table,i)         = blob_val(delta*i, blob)  *iw0;
		blobTableSqrt[i] = blob_val(blobTableSize*sqrt((double)i), blob)  *iw0;
        //***
        //DIRECT_VEC_ELEM(fourierBlobTableSqrt,i) =
        //     blob_Fourier_val(fourierBlobTableSize*sqrt(i), blobFourier)*padXdim3  *iw0;
        VEC_ELEM(Fourier_blob_table,i) =
            blob_Fourier_val(deltaFourier*i, blobFourier)*padXdim3  *iw0;
        //#define DEBUG
#ifdef DEBUG

        std::cout << VEC_ELEM(Fourier_blob_table,i)
        << " " << VEC_ELEM(fourierBlobTableSqrt,i)
        << std::endl;
#endif
  #undef DEBUG

    }
    //iDelta        = 1/delta;
    iDeltaSqrt    = 1/deltaSqrt;
    iDeltaFourier = 1/deltaFourier;

    // Get symmetries
    Matrix2D<double>  Identity(3,3);
    Identity.initIdentity();
    R_repository.push_back(Identity);
    if (fn_sym != "")
    {
        SymList SL;
        SL.readSymmetryFile(fn_sym);
        for (int isym = 0; isym < SL.symsNo(); isym++)
        {
            Matrix2D<double>  L(4, 4), R(4, 4);
            SL.getMatrices(isym, L, R);
            R.resize(3, 3);
            R_repository.push_back(R);
        }
    }

    // query the system on the number of cores
	noOfCores = 1;//std::max(1l, sysconf(_SC_NPROCESSORS_ONLN));
}

void ProgRecFourierGPU::cropAndShift(
		MultidimArray<std::complex<double> >& paddedFourier,
		ProgRecFourierGPU* parent,
		RecFourierBufferData* buffer,
		float* dest) {
	int sizeX = buffer->fftSizeX;
	int sizeY = buffer->fftSizeY;

	// convert image (shift to center and remove high frequencies)
	std::complex<double> paddedFourierTmp;
	int halfY = paddedFourier.ydim / 2;
	double tempMyPadd[2];
	for (size_t i = 0; i < paddedFourier.ydim; i++) {
		for (size_t j = 0; j < sizeX; j++) {
			if (i < sizeX || i >= (paddedFourier.ydim - sizeX)) {
				// check the frequency
				paddedFourierTmp = DIRECT_A2D_ELEM(paddedFourier, i, j);
				FFT_IDX2DIGFREQ(j, parent->paddedImgSize, tempMyPadd[0]);
				FFT_IDX2DIGFREQ(i, parent->paddedImgSize, tempMyPadd[1]);
				if (tempMyPadd[0] * tempMyPadd[0] + tempMyPadd[1] * tempMyPadd[1]> parent->maxResolutionSqr) {
					continue;
				}
				// do the shift
				int myPadI = (i < halfY) ?	i + sizeX : i - paddedFourier.ydim + sizeX;
				int index = myPadI * sizeX + j;
				dest[index * 2] = paddedFourierTmp.real();
				dest[index * 2 + 1] = paddedFourierTmp.imag();
			}
		}
	}
}

void ProgRecFourierGPU::prepareBuffer(RecFourierWorkThread* threadParams,
		ProgRecFourierGPU* parent,
		bool hasCTF, std::vector<size_t>& objId)
{
    ApplyGeoParams params;
    params.only_apply_shifts = true;
    MultidimArray<float> localPaddedImgFloat;
    MultidimArray<double> localPaddedImgDouble;
	FourierTransformer localTransformerImg;
	// FIXME cannot be complex float due to build errors
	MultidimArray< std::complex<double> > localPaddedFourier;

	RecFourierBufferData* buffer = threadParams->buffer;
	buffer->noOfImages = 0; // 'clean' buffer
	for (int projIndex = 0; projIndex < parent->bufferSize; projIndex++) {
		double rot, tilt, psi, weight;
		Projection proj;
		Matrix2D<double> Ainv(3, 3);
		int imgIndex = threadParams->startImageIndex + projIndex;

		if (imgIndex >= threadParams->endImageIndex) {
			continue;
		}
		//Read projection from selfile, read also angles and shifts if present
		//but only apply shifts (set above)
		// FIXME following line is a current bottleneck, as it calls BSpline interpolation
		proj.readApplyGeo(*(threadParams->selFile), objId[imgIndex], params);
		if (parent->do_weights && proj.weight() == 0.f) {
			continue;
		}

		// Compute the coordinate axes associated to this projection
		rot = proj.rot();
		tilt = proj.tilt();
		psi = proj.psi();
		Euler_angles2matrix(rot, tilt, psi, Ainv);

		Ainv = Ainv.transpose();
		float transf[3][3];
		float transfInv[3][3];

		// prepare transforms for all  symmetries
		int travSpaceOffset = buffer->noOfImages * buffer->noOfSymmetries;
		for (int j = 0; j < parent->R_repository.size(); j++) {
			RecFourierProjectionTraverseSpace* space = &buffer->spaces[travSpaceOffset + j];

			Matrix2D<double> A_SL=parent->R_repository[j]*(Ainv);
			Matrix2D<double> A_SLInv=A_SL.inv();
			A_SL.convertTo(transf);
			A_SLInv.convertTo(transfInv);

			parent->computeTraverseSpace(buffer->fftSizeX, buffer->fftSizeY, projIndex,
				transf, transfInv, space);
			space->weight = (parent->do_weights) ? proj.weight() : 1.0;
		}

		// Copy the projection to the center of the padded image
		// and compute its Fourier transform, if requested
		proj().setXmippOrigin();
		const MultidimArray<double> &mProj = proj();
		if (buffer->hasFFTs) {
			localPaddedImgDouble.initZeros(buffer->paddedImgSize, buffer->paddedImgSize);
			localPaddedImgDouble.setXmippOrigin();
			FOR_ALL_ELEMENTS_IN_ARRAY2D(mProj)
				A2D_ELEM(localPaddedImgDouble,i,j) = A2D_ELEM(mProj, i, j);
			CenterFFT(localPaddedImgDouble, true);

			// Fourier transformer for the images
			localTransformerImg.setReal(localPaddedImgDouble);
			localTransformerImg.FourierTransform();
			localTransformerImg.getFourierAlias(localPaddedFourier);
			cropAndShift(localPaddedFourier, parent,
					buffer, buffer->getNthItem(buffer->FFTs, projIndex));
		} else {
			localPaddedImgFloat.initZeros(parent->paddedImgSize, parent->paddedImgSize);
			localPaddedImgFloat.setXmippOrigin();
			FOR_ALL_ELEMENTS_IN_ARRAY2D(mProj)
				A2D_ELEM(localPaddedImgFloat,i,j) = A2D_ELEM(mProj, i, j);
			CenterFFT(localPaddedImgFloat, true);

			// add image at the end of the stack (that is already long enough)
			memcpy(buffer->getNthItem(buffer->paddedImages, projIndex),
					localPaddedImgFloat.data,
					buffer->getPaddedImgByteSize());
		}

		if (hasCTF) {
			computeCTFCorrection(threadParams, objId[imgIndex],parent, buffer, projIndex);
		}

		buffer->noOfImages++; // new image added to buffer
	}
}


void* ProgRecFourierGPU::threadRoutine(void* threadArgs) {
	XMIPP_TRY // in case some method throws a xmipp exception

    RecFourierWorkThread* threadParams = (RecFourierWorkThread *) threadArgs;
    ProgRecFourierGPU* parent = threadParams->parent;
    std::vector<size_t> objId;

    threadParams->selFile = &parent->SF;
    threadParams->selFile->findObjects(objId);
    bool hasCTF = parent->useCTF
    		&& (threadParams->selFile->containsLabel(MDL_CTF_MODEL)
    				|| threadParams->selFile->containsLabel(MDL_CTF_DEFOCUSU));

    // allocate buffer
    threadParams->buffer = new RecFourierBufferData( ! parent->fftOnGPU, hasCTF,
    		parent->maxVolumeIndexX / 2, parent->maxVolumeIndexYZ, parent->paddedImgSize,
			parent->bufferSize, (int)parent->R_repository.size());
    // allocate GPU
    allocateWrapper(threadParams->buffer, threadParams->gpuStream);


	threadParams->endImageIndex = std::min(threadParams->endImageIndex+1, threadParams->startImageIndex+parent->bufferSize);
	prepareBuffer(threadParams, parent, false, objId);
	int noOfSpaces = threadParams->buffer->getNoOfElements(threadParams->buffer->spaces);

    size_t deviceIndex = 0;
	std::string kernelFile = "/home/david/GIT/Scipion/software/em/xmipp/libraries/reconstruction_cuda/reconstruct_fourier.cu";
	std::string referenceKernelFile = "/home/david/GIT/Scipion/software/em/xmipp/libraries/reconstruction_cuda/reconstruct_fourier_ref.cu";

	// Create tuner object for specified device, platform index is ignored in case of CUDA API usage
	ktt::Tuner tuner(0, deviceIndex, ktt::ComputeApi::Cuda);

	// Add new kernel to tuner, specify kernel name, grid dimensions and block dimensions
	ktt::KernelId kernelId = tuner.addKernelFromFile(kernelFile, "processBufferKernel", ktt::DimensionVector(1), ktt::DimensionVector(1));
	int localSize = 16;
	int size2D = parent->maxVolumeIndexX + 1;
	int globalSize = ceil(size2D/(float)localSize);
	ktt::KernelId referenceKernelId = tuner.addKernelFromFile(referenceKernelFile, "processBufferKernel", ktt::DimensionVector(globalSize, globalSize), ktt::DimensionVector(localSize, localSize));

	int volumeSize = std::pow(parent->maxVolumeIndexYZ + 1, 3);
	parent->tempVolumeGPU = new float[volumeSize * 2];
	parent->tempWeightsGPU = new float[volumeSize];

	// Add new arguments to tuner, argument data is copied from std::vector containers
	ktt::ArgumentId volId = tuner.addArgumentVector(std::vector<float>(parent->tempVolumeGPU, parent->tempVolumeGPU+volumeSize*2), ktt::ArgumentAccessType::ReadWrite);
	ktt::ArgumentId weightId = tuner.addArgumentVector(std::vector<float>(parent->tempWeightsGPU, parent->tempWeightsGPU +volumeSize), ktt::ArgumentAccessType::ReadWrite);
	ktt::ArgumentId spaceId = tuner.addArgumentVector(std::vector<RecFourierProjectionTraverseSpace>(threadParams->buffer->spaces, threadParams->buffer->spaces + noOfSpaces), ktt::ArgumentAccessType::ReadOnly);
	ktt::ArgumentId spaceNoId = tuner.addArgumentScalar(noOfSpaces);
	ktt::ArgumentId FFTsId = tuner.addArgumentVector(std::vector<float>(threadParams->buffer->FFTs, threadParams->buffer->FFTs+threadParams->buffer->getNoOfElements(threadParams->buffer->FFTs)), ktt::ArgumentAccessType::ReadOnly);
	ktt::ArgumentId fftSizeXId = tuner.addArgumentScalar(threadParams->buffer->fftSizeX);
	ktt::ArgumentId fftSizeYId = tuner.addArgumentScalar(threadParams->buffer->fftSizeY);
	ktt::ArgumentId blobTableId = tuner.addArgumentVector(std::vector<float>((float*)parent->blobTableSqrt, (float*)parent->blobTableSqrt+BLOB_TABLE_SIZE_SQRT), ktt::ArgumentAccessType::ReadOnly);
	ktt::ArgumentId imgCacheId = tuner.addArgumentScalar(0);
	ktt::ArgumentId sharedMemId = tuner.addArgumentLocal<std::complex<float> >(1); // will be set eventually
	ktt::ArgumentId maxVolIndexXId = tuner.addArgumentScalar(parent->maxVolumeIndexX);
	ktt::ArgumentId maxVolIndexYZId = tuner.addArgumentScalar(parent->maxVolumeIndexYZ);


	tuner.addParameter(kernelId, "BLOCK_DIM_X", {8,12,16,20,24,28,32}, ktt::ThreadModifierType::Local, ktt::ThreadModifierAction::Multiply, ktt::Dimension::X);
	tuner.addParameter(kernelId, "BLOCK_DIM_Y", {8,12,16,20,24,28,32}, ktt::ThreadModifierType::Local, ktt::ThreadModifierAction::Multiply, ktt::Dimension::Y);
	tuner.addParameter(kernelId, "TILE", {1,2,4,8});

	tuner.addParameter(kernelId, "GRID_DIM_Z", {1,4,8,16}, ktt::ThreadModifierType::Global, ktt::ThreadModifierAction::Multiply, ktt::Dimension::Z);

	tuner.addParameter(kernelId, "SHARED_BLOB_TABLE", {0,1});
	tuner.addParameter(kernelId, "SHARED_IMG", {0, 1});
	tuner.addParameter(kernelId, "USE_ATOMICS", {0, 1});
	tuner.addParameter(kernelId, "BLOB_TABLE_SIZE_SQRT", {BLOB_TABLE_SIZE_SQRT});
	tuner.addParameter(kernelId, "PRECOMPUTE_BLOB_VAL", {0,1});
	tuner.addParameter(kernelId, "cMaxVolumeIndexX", {parent->maxVolumeIndexX});
	tuner.addParameter(kernelId, "cMaxVolumeIndexYZ", {parent->maxVolumeIndexYZ});
	tuner.addParameter(kernelId, "blobOrder", {parent->blob.order});

//	tuner.addParameter(kernelId, "cBlobRadius", {parent->blob.radius});
//	tuner.addParameter(kernelId, "cBlobAlpha", {parent->blob.alpha});
//	tuner.addParameter(kernelId, "cIw0", {parent->iw0});
//	tuner.addParameter(kernelId, "cIDeltaSqrt", {parent->iDeltaSqrt});

	auto blocksDimEqConstr = [](std::vector<size_t> vector) {return vector.at(0)== vector.at(1);};
	tuner.addConstraint(kernelId, blocksDimEqConstr, std::vector<std::string>{"BLOCK_DIM_X", "BLOCK_DIM_Y"});

	auto tileMultXConstr = [](std::vector<size_t> vector) {return vector.at(1) == 1 || (vector.at(0) % vector.at(1) == 0);};
	tuner.addConstraint(kernelId, tileMultXConstr, std::vector<std::string>{"BLOCK_DIM_X", "TILE"});

	auto tileMultYConstr = [](std::vector<size_t> vector) {return vector.at(1) == 1 || (vector.at(0) % vector.at(1) == 0);};
	tuner.addConstraint(kernelId, tileMultYConstr, std::vector<std::string>{"BLOCK_DIM_Y", "TILE"});

	auto tileSharedImgConstr = [](std::vector<size_t> vector) {return vector.at(0) == 0 || vector.at(1) == 1;};
	tuner.addConstraint(kernelId, tileMultYConstr, std::vector<std::string>{"SHARED_IMG", "TILE"});

	auto useAtomicsZDimConstr = [](std::vector<size_t> vector) {return !(vector.at(0) == 0 && vector.at(1) != 1);};
	tuner.addConstraint(kernelId, tileMultYConstr, std::vector<std::string>{"USE_ATOMICS", "GRID_DIM_Z"});

	auto tileSmallerThanBlockConstr = [](std::vector<size_t> vector) {return vector.at(0) > vector.at(2) && vector.at(1) > vector.at(2);};
	tuner.addConstraint(kernelId, tileSmallerThanBlockConstr, std::vector<std::string>{"BLOCK_DIM_X", "BLOCK_DIM_Y", "TILE"});

	auto tooMuchSharedMemConstr = [](std::vector<size_t> vector) {return !(vector.at(0)==1 && vector.at(1)==1);};
	tuner.addConstraint(kernelId, tooMuchSharedMemConstr, std::vector<std::string>{"SHARED_BLOB_TABLE", "SHARED_IMG"});

	auto blobTableConstr = [](std::vector<size_t> vector) {return vector.at(0)==0 || (vector.at(0)==1 && vector.at(1)==1);};
	tuner.addConstraint(kernelId, blobTableConstr, std::vector<std::string>{"SHARED_BLOB_TABLE", "PRECOMPUTE_BLOB_VAL"});

	tuner.setTuningManipulator(kernelId, std::make_unique<Manipulator>(parent,objId,threadParams->buffer,
			imgCacheId,spaceId,spaceNoId,FFTsId, sharedMemId, threadParams->startImageIndex, threadParams->endImageIndex));

	// Set kernel arguments by providing corresponding argument ids returned by addArgument() method, order of arguments is important
	tuner.setKernelArguments(kernelId, std::vector<ktt::ArgumentId>{volId, weightId,
		spaceId, spaceNoId,
		FFTsId,
		fftSizeXId, fftSizeYId, blobTableId, imgCacheId, sharedMemId});
	tuner.setKernelArguments(referenceKernelId, std::vector<ktt::ArgumentId>{volId, weightId,
		spaceId, spaceNoId,
		FFTsId,
		fftSizeXId, fftSizeYId, blobTableId, maxVolIndexXId, maxVolIndexYZId});


	// Specify custom tolerance threshold for validation of floating point arguments. Default threshold is 1e-4.
	tuner.setValidationMethod(ktt::ValidationMethod::SideBySideRelativeComparison, 0.01f);
    tuner.setReferenceKernel(kernelId, referenceKernelId, {}, std::vector<ktt::ArgumentId>{volId,
        weightId});


        // Set reference class, which implements C++ version of kernel computation in order to validate results provided by kernel,
        // provide list of arguments which will be validated
    //    tuner.setReferenceClass(kernelId, std::make_unique<SimpleReferenceClass>(data, dataId), std::vector<ktt::ArgumentId>{dataId});



    //    tuner.setArgumentComparator(dataId, compareData);


	// Print tuning results to standard output and to output.csv file
	// Launch kernel tuning
//	for (auto i = 0; i < 2; i++) {
//		std::cout << "run no " << i << std::endl;
//	tuner.tuneKernelByStep(kernelId, {
//			ktt::ArgumentOutputDescriptor(volId, parent->tempVolumeGPU),
//			ktt::ArgumentOutputDescriptor(weightId, parent->tempWeightsGPU)});
//	}
	tuner.tuneKernel(kernelId);
	tuner.printResult(kernelId, std::cout, ktt::PrintFormat::Verbose);
	size_t lastindex = parent->fn_out.getString().find_last_of(".");
	std::string rawname = parent->fn_out.getString().substr(0, lastindex);
	tuner.printResult(kernelId, rawname + "_results.csv", ktt::PrintFormat::CSV);

	// TODO prekopirovat pamet z CPU na CPU


    // clean after itself
    releaseWrapper(threadParams->gpuStream);
    delete threadParams->buffer;
    threadParams->buffer = NULL;
    threadParams->selFile = NULL;
    barrier_wait( &parent->barrier );// notify that thread finished

	XMIPP_CATCH // catch possible exception
	return NULL;
}

inline void ProgRecFourierGPU::multiply(const float transform[3][3], Point3D<float>& inOut) {
	float tmp0 = transform[0][0] * inOut.x + transform[0][1] * inOut.y + transform[0][2] * inOut.z;
	float tmp1 = transform[1][0] * inOut.x + transform[1][1] * inOut.y + transform[1][2] * inOut.z;
	float tmp2 = transform[2][0] * inOut.x + transform[2][1] * inOut.y + transform[2][2] * inOut.z;
	inOut.x = tmp0;
	inOut.y = tmp1;
	inOut.z = tmp2;
}

void ProgRecFourierGPU::createProjectionCuboid(Point3D<float>* cuboid, float sizeX, float sizeY, float blobSize)
{
	float halfY = sizeY / 2.0f;
	cuboid[3].x = cuboid[2].x = cuboid[7].x = cuboid[6].x = 0.f - blobSize;
	cuboid[0].x = cuboid[1].x = cuboid[4].x = cuboid[5].x = sizeX + blobSize;

	cuboid[3].y = cuboid[0].y = cuboid[7].y = cuboid[4].y = -(halfY + blobSize);
	cuboid[1].y = cuboid[2].y = cuboid[5].y = cuboid[6].y = halfY + blobSize;

	cuboid[3].z = cuboid[0].z = cuboid[1].z = cuboid[2].z = 0.f + blobSize;
	cuboid[7].z = cuboid[4].z = cuboid[5].z = cuboid[6].z = 0.f - blobSize;
}

inline void ProgRecFourierGPU::translateCuboid(Point3D<float>* cuboid, Point3D<float> vector) {
	for (int i = 0; i < 8; i++) {
		cuboid[i].x += vector.x;
		cuboid[i].y += vector.y;
		cuboid[i].z += vector.z;
	}
}

void ProgRecFourierGPU::computeAABB(Point3D<float>* AABB, Point3D<float>* cuboid,
	float minX, float minY, float minZ,
	float maxX, float maxY, float maxZ) {
	AABB[0].x = AABB[0].y = AABB[0].z = std::numeric_limits<float>::max();
	AABB[1].x = AABB[1].y = AABB[1].z = std::numeric_limits<float>::min();
	Point3D<float> tmp;
	for (int i = 0; i < 8; i++) {
		tmp = cuboid[i];
		if (AABB[0].x > tmp.x) AABB[0].x = tmp.x;
		if (AABB[0].y > tmp.y) AABB[0].y = tmp.y;
		if (AABB[0].z > tmp.z) AABB[0].z = tmp.z;
		if (AABB[1].x < tmp.x) AABB[1].x = tmp.x;
		if (AABB[1].y < tmp.y) AABB[1].y = tmp.y;
		if (AABB[1].z < tmp.z) AABB[1].z = tmp.z;
	}
	// limit to max size
	if (AABB[0].x < minX) AABB[0].x = minX;
	if (AABB[0].y < minY) AABB[0].y = minY;
	if (AABB[0].z < minZ) AABB[0].z = minZ;
	if (AABB[1].x > maxX) AABB[1].x = maxX;
	if (AABB[1].y > maxY) AABB[1].y = maxY;
	if (AABB[1].z > maxZ) AABB[1].z = maxZ;
}

void static printAABB(Point3D<float> AABB[]) {
	std::cout
	// one base
		<< AABB[0].x << " " << AABB[0].y << " " << AABB[0].z << "\n"
		<< AABB[1].x << " " << AABB[0].y << " " << AABB[0].z << "\n"
		<< AABB[1].x << " " << AABB[1].y << " " << AABB[0].z << "\n"
		<< AABB[0].x << " " << AABB[1].y << " " << AABB[0].z << "\n"
		<< AABB[0].x << " " << AABB[0].y << " " << AABB[0].z << "\n"
	// other base with one connection
		<< AABB[0].x << " " << AABB[0].y << " " << AABB[1].z << "\n"
		<< AABB[1].x << " " << AABB[0].y << " " << AABB[1].z << "\n"
		<< AABB[1].x << " " << AABB[1].y << " " << AABB[1].z << "\n"
		<< AABB[0].x << " " << AABB[1].y << " " << AABB[1].z << "\n"
		<< AABB[0].x << " " << AABB[0].y << " " << AABB[1].z << "\n"
	// lines between bases
		<< AABB[1].x << " " << AABB[0].y << " " << AABB[1].z << "\n"
		<< AABB[1].x << " " << AABB[0].y << " " << AABB[0].z << "\n"
		<< AABB[1].x << " " << AABB[1].y << " " << AABB[0].z << "\n"
		<< AABB[1].x << " " << AABB[1].y << " " << AABB[1].z << "\n"
		<< AABB[0].x << " " << AABB[1].y << " " << AABB[1].z << "\n"
		<< AABB[0].x << " " << AABB[1].y << " " << AABB[0].z
		<< std::endl;
}

inline void ProgRecFourierGPU::computeCTFCorrection(RecFourierWorkThread* threadParams,
		size_t imgIndex,
		ProgRecFourierGPU* parent,
		RecFourierBufferData* buffer,
		int storeIndex) {
	CTFDescription ctf;
	ctf.readFromMetadataRow(*(threadParams->selFile), imgIndex);
	ctf.produceSideInfo();
	float freqX, freqY;
	float CTFVal, modulatorVal;
	for (int y = 0; y < buffer->fftSizeY; y++) {
		// since Y axis is shifted to center, we have to use different calculation
		freqY = (y - (parent->paddedImgSize / 2.f)) / (float) parent->paddedImgSize;
		for (int x = 0; x < buffer->fftSizeY; x++) {
			CTFVal = modulatorVal = 1.f;
			// get respective frequency
			FFT_IDX2DIGFREQ(x, parent->paddedImgSize, freqX);
			ctf.precomputeValues(freqX * parent->iTs, freqY * parent->iTs);
			CTFVal = ctf.getValuePureNoKAt();
			if (std::isnan(CTFVal)) {
				if ((x == 0) && (y == 0)) {
					modulatorVal = CTFVal = 1.0;
				}
				else {
					modulatorVal = CTFVal = 0.0;
				}
			}
			if (fabs(CTFVal) < parent->minCTF) {
				modulatorVal = fabs(CTFVal);
				CTFVal = SGN(CTFVal);
			} else {
				CTFVal = 1.0 / CTFVal;
			}
			if (parent->isPhaseFlipped)
				CTFVal = fabs(CTFVal);

			int index = y * buffer->fftSizeX + x;
			buffer->getNthItem(buffer->CTFs, storeIndex)[index] = CTFVal;
			buffer->getNthItem(buffer->modulators, storeIndex)[index] = modulatorVal;
		}
	}
}

template<typename T>
T*** ProgRecFourierGPU::allocate(T***& where, int xSize, int ySize, int zSize) {
	where = new T**[zSize];
	for (int z = 0; z < zSize; z++) {
		where[z] = new T*[ySize];
		for (int y = 0; y < ySize; y++) {
			where[z][y] = new T[xSize];
			for (int x = 0; x < xSize; x++) {
				where[z][y][x] = (T) 0;
			}
		}
	}
	return where;
}

template<typename T>
void ProgRecFourierGPU::release(T***& array, int ySize, int zSize) {
	for(int z = 0; z < zSize; z++) {
		for(int y = 0; y < ySize; y++) {
			delete[] array[z][y];
		}
		delete[] array[z];
	}
	delete[] array;
	array = NULL;
}

template<typename T>
T*** ProgRecFourierGPU::applyBlob(T***& input, float blobSize,
		float* blobTableSqrt, float iDeltaSqrt) {
	float blobSizeSqr = blobSize * blobSize;
	int blob = floor(blobSize); // we are using integer coordinates, so we cannot hit anything further
	T tmp;
	T*** output;
	// create new storage
	allocate(output, maxVolumeIndexX+1, maxVolumeIndexYZ+1, maxVolumeIndexYZ+1);

	// traverse new storage
	for (int i = 0; i <= maxVolumeIndexYZ; i++) {
		for (int j = 0; j <= maxVolumeIndexYZ; j++) {
			for (int k = 0; k <= maxVolumeIndexX; k++) {
				// traverse input storage
				tmp = (T) 0;
				for (int z = std::max(0, i-blob); z <= std::min(maxVolumeIndexYZ, i+blob); z++) {
					float dZSqr = (i - z) * (i - z);
					for (int y = std::max(0, j-blob); y <= std::min(maxVolumeIndexYZ, j+blob); y++) {
						float dYSqr = (j - y) * (j - y);
						for (int x = std::max(0, k-blob); x <= std::min(maxVolumeIndexX, k+blob); x++) {
							float dXSqr = (k - x) * (k - x);
							float distanceSqr = dZSqr + dYSqr + dXSqr;
							if (distanceSqr > blobSizeSqr) {
								continue;
							}
							int aux = (int) ((distanceSqr * iDeltaSqrt + 0.5)); //Same as ROUND but avoid comparison
							float tmpWeight = blobTableSqrt[aux];
							tmp += tmpWeight * input[z][y][x];
						}
					}
				}
				output[i][j][k] = tmp;
			}
		}
	}
	// free original data
	release(input, maxVolumeIndexYZ+1, maxVolumeIndexYZ+1);
	return output;
}

template<typename T, typename U>
void ProgRecFourierGPU::convertToExpectedSpace(T*** input, int size,
	MultidimArray<U>& VoutFourier) {
	int halfSize = size / 2;
	for (int z = 0; z <= size; z++) {
    	for (int y = 0; y <= size; y++) {
			for (int x = 0; x <= halfSize; x++) {
    			int newPos[3];
    			// shift FFT from center to corners
				newPos[0] = x; // no need to move
    			newPos[1] = (y < halfSize) ? VoutFourier.ydim - halfSize + y : y - halfSize ;
    			newPos[2] = (z < halfSize) ? VoutFourier.zdim - halfSize + z : z - halfSize ;
    			// store to output array
    			// += in necessary as VoutFourier might be used multiple times when used with MPI
				DIRECT_A3D_ELEM(VoutFourier, newPos[2], newPos[1], newPos[0]) += input[z][y][x];
    		}
    	}
    }
}

void ProgRecFourierGPU::getGPUData() {
	if (NULL == tempVolume) {
		allocate(tempVolume, maxVolumeIndexYZ + 1, maxVolumeIndexYZ + 1,
				maxVolumeIndexYZ + 1);
	}
	if (NULL == tempWeights) {
		allocate(tempWeights, maxVolumeIndexYZ + 1, maxVolumeIndexYZ + 1,
				maxVolumeIndexYZ + 1);
	}
	int size =  maxVolumeIndexYZ + 1;
	for (int z = 0; z < size; z++) {
		for (int y = 0; y < size; y++) {
			int index = (z * size * size) + (y * size);
			memcpy(tempVolume[z][y], &tempVolumeGPU[2 * index], 2 * size * sizeof(float));
			memcpy(tempWeights[z][y] , &tempWeightsGPU[index], size * sizeof(float));
		}
	}



	delete[] tempVolumeGPU;
	delete[] tempWeightsGPU;
}

void ProgRecFourierGPU::mirrorAndCropTempSpaces() {
	maxVolumeIndexX = maxVolumeIndexYZ/2; // just half of the space is necessary, the rest is complex conjugate
	mirrorAndCrop(tempWeights, &identity<float>);
	mirrorAndCrop(tempVolume, &conjugate);
}

template<typename T>
void ProgRecFourierGPU::mirrorAndCrop(T***& input, T (*f)(T)) {
	T*** output;
	// create new storage, notice that just 'right hand side - X axis' of the input will be preserved, left will be converted to its complex conjugate
	allocate(output, maxVolumeIndexX+1, maxVolumeIndexYZ+1, maxVolumeIndexYZ+1);
	// traverse old storage
	for (int z = 0; z <= maxVolumeIndexYZ; z++) {
		for (int y = 0; y <= maxVolumeIndexYZ; y++) {
			for (int x = 0; x <= maxVolumeIndexYZ; x++) {
				if (x < maxVolumeIndexX) {
					int newPos[3];
					// mirror against center of the volume, e.g. [0,0,0]->[size,size,size]. It will fit as the input space is one voxel bigger
					newPos[0] = maxVolumeIndexYZ - x;
					newPos[1] = maxVolumeIndexYZ - y;
					newPos[2] = maxVolumeIndexYZ - z;
					output[newPos[2]][newPos[1]][newPos[0]-maxVolumeIndexX] += f(input[z][y][x]);
				} else {
					// copy with X shifted by (-halfSize)
					output[z][y][x-maxVolumeIndexX] += input[z][y][x];
				}
			}
		}
	}
	// free original data
	release(input, maxVolumeIndexYZ+1, maxVolumeIndexYZ+1);
	// set new data
	input = output;
}

void ProgRecFourierGPU::forceHermitianSymmetry() {
	int x = 0;
	for (int z = 0; z <= maxVolumeIndexYZ; z++) {
		for (int y = 0; y <= maxVolumeIndexYZ/2; y++) {
			int newPos[3];
			// mirror against center of the volume, e.g. [0,0,0]->[size,size,size]. It will fit as the input space is one voxel biger
			newPos[0] = x;
			newPos[1] = maxVolumeIndexYZ - y;
			newPos[2] = maxVolumeIndexYZ - z;
			std::complex<float> tmp1 = 0.5f * (tempVolume[newPos[2]][newPos[1]][newPos[0]] + conj(tempVolume[z][y][x]));
			float tmp2 = 0.5f * (tempWeights[newPos[2]][newPos[1]][newPos[0]] + tempWeights[z][y][x]);

			tempVolume[newPos[2]][newPos[1]][newPos[0]] = tmp1;
			tempVolume[z][y][x] = conj(tmp1);
			tempWeights[newPos[2]][newPos[1]][newPos[0]] = tempWeights[z][y][x] = tmp2;
		}
	}
}

void ProgRecFourierGPU::processWeights() {

//	Image<double> Vout;
//	Vout().initZeros(maxVolumeIndexYZ+1,maxVolumeIndexYZ+1,maxVolumeIndexX);
//	for (int z = 0; z <= maxVolumeIndexYZ; z++ ) {
//		for (int y = 0; y <= maxVolumeIndexYZ; y++ ) {
//			for (int x = 0; x <= maxVolumeIndexX; x++ ) {
//				DIRECT_A3D_ELEM(Vout.data, z, y, x) = tempVolume[z][y][x].real();
//			}
//		}
//	}
//	Vout.write("test_objemu.vol");

    // Get a first approximation of the reconstruction
    float corr2D_3D=pow(padding_factor_proj,2.)/
                     (imgSize* pow(padding_factor_vol,3.));
	for (int z = 0; z <= maxVolumeIndexYZ; z++) {
		for (int y = 0; y <= maxVolumeIndexYZ; y++) {
			for (int x = 0; x <= maxVolumeIndexX; x++) {
				float weight = tempWeights[z][y][x];
				std::complex<float> val = tempVolume[z][y][x];
				if (fabs(weight) > 1e-3) {
					weight = 1.f/weight;
				}

				if (1.0/weight > ACCURACY)
					tempVolume[z][y][x] *= corr2D_3D*weight;
				else
					tempVolume[z][y][x] = 0;
			}
		}
	}
}

void ProgRecFourierGPU::computeTraverseSpace(int imgSizeX, int imgSizeY, int projectionIndex,
		MATRIX& transform, MATRIX& transformInv, RecFourierProjectionTraverseSpace* space) {
	Point3D<float> cuboid[8];
	Point3D<float> AABB[2];
	Point3D<float> origin = {maxVolumeIndexX/2.f, maxVolumeIndexYZ/2.f, maxVolumeIndexYZ/2.f};
	createProjectionCuboid(cuboid, imgSizeX, imgSizeY, useFast ? 0.f : blob.radius);
	rotateCuboid(cuboid, transform);
	translateCuboid(cuboid, origin);
	computeAABB(AABB, cuboid, 0, 0, 0, maxVolumeIndexX, maxVolumeIndexYZ, maxVolumeIndexYZ);

	// store data
	space->projectionIndex = projectionIndex;
	space->minZ = floor(AABB[0].z);
	space->minY = floor(AABB[0].y);
	space->minX = floor(AABB[0].x);
	space->maxZ = ceil(AABB[1].z);
	space->maxY = ceil(AABB[1].y);
	space->maxX = ceil(AABB[1].x);
	space->topOriginX = cuboid[4].x;
	space->topOriginY = cuboid[4].y;
	space->topOriginZ = cuboid[4].z;
	space->bottomOriginX = cuboid[0].x;
	space->bottomOriginY = cuboid[0].y;
	space->bottomOriginZ = cuboid[0].z;
	space->maxDistanceSqr = (imgSizeX+(useFast ? 0.f : blob.radius))
			* (imgSizeX+(useFast ? 0.f : blob.radius));

	space->transformInv00 = transformInv[0][0];
	space->transformInv01 = transformInv[0][1];
	space->transformInv02 = transformInv[0][2];

	space->transformInv10 = transformInv[1][0];
	space->transformInv11 = transformInv[1][1];
	space->transformInv12 = transformInv[1][2];

	space->transformInv20 = transformInv[2][0];
	space->transformInv21 = transformInv[2][1];
	space->transformInv22 = transformInv[2][2];

	// calculate best traverse direction
	Point3D<float> normal;
	normal.x = normal.y = 0.f;
	normal.z = 1.f;
	multiply(transform, normal);
	space->unitNormalX = normal.x;
	space->unitNormalY = normal.y;
	space->unitNormalZ = normal.z;
	float nX = std::abs(normal.x);
	float nY = std::abs(normal.y);
	float nZ = std::abs(normal.z);

	// biggest vector indicates ideal direction
	if (nX >= nY  && nX >= nZ) { // iterate YZ plane
		space->dir = space->YZ;
	} else if (nY >= nX && nY >= nZ) { // iterate XZ plane
		space->dir = space->XZ;
	} else if (nZ >= nX && nZ >= nY) { // iterate XY plane
		space->dir = space->XY;
	}
}

void ProgRecFourierGPU::logProgress(int increment, bool reset) {
	static int repaintAfter = (int)ceil((double)SF.size()/60);
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static int noOfDone = 0;
	static int noOfLogs = 0;
	if (reset) {
		noOfDone = noOfLogs = 0;
	}
	if (verbose) {
		pthread_mutex_lock(&mutex);
		noOfDone += increment;
		while (noOfLogs < (noOfDone / repaintAfter)) {
			progress_bar(noOfDone);
			noOfLogs++;
		}
		pthread_mutex_unlock(&mutex);
	}
}

void ProgRecFourierGPU::processImages( int firstImageIndex, int lastImageIndex)
{


	// initialize GPU
//    if (NULL == tempVolumeGPU) {
//    	allocateTempVolumeGPU(tempVolumeGPU, maxVolumeIndexYZ+1, sizeof(std::complex<float>));
//    }
//    if (NULL == tempWeightsGPU) {
//    	allocateTempVolumeGPU(tempWeightsGPU, maxVolumeIndexYZ+1, sizeof(float));
//    }
    createStreams(noOfCores);
//    copyConstants(maxVolumeIndexX, maxVolumeIndexYZ,
//    		blob.radius, blob.alpha, iDeltaSqrt, iw0);
    if ( ! useFast) {
    	copyBlobTable(blobTableSqrt, BLOB_TABLE_SIZE_SQRT);
    }

	// create threads
	int imgPerThread = ceil((lastImageIndex-firstImageIndex+1) / (float)noOfCores);
	workThreads = new RecFourierWorkThread[noOfCores];
	barrier_init( &barrier, noOfCores + 1 ); // + main thread
	for (int i = 0; i < noOfCores; i++) {
		int sIndex = firstImageIndex + i*imgPerThread;
		int eIndex = std::min(lastImageIndex, sIndex + imgPerThread-1);
		createWorkThread(i, sIndex, eIndex, workThreads[i]);
	}

	// Waiting for threads to finish
	barrier_wait( &barrier );

	// clean threads and GPU
	for (int i = 0; i < noOfCores; i++) {
		pthread_join(workThreads[i].id, NULL);
	}
	barrier_destroy( &barrier );
	delete[] workThreads;
	releaseBlobTable();
	deleteStreams(noOfCores);
}

void ProgRecFourierGPU::releaseTempSpaces() {
	release(tempWeights, maxVolumeIndexYZ+1, maxVolumeIndexYZ+1);
	release(tempVolume, maxVolumeIndexYZ+1, maxVolumeIndexYZ+1);
}

void ProgRecFourierGPU::finishComputations( const FileName &out_name )
{
	if (useFast) {
		tempVolume = applyBlob(tempVolume, blob.radius, blobTableSqrt, iDeltaSqrt);
		tempWeights = applyBlob(tempWeights, blob.radius, blobTableSqrt, iDeltaSqrt);
	}

	forceHermitianSymmetry();
	processWeights();
	release(tempWeights, maxVolumeIndexYZ+1, maxVolumeIndexYZ+1);
	MultidimArray< std::complex<double> > VoutFourier;
	allocateVoutFourier(VoutFourier);
	convertToExpectedSpace(tempVolume, maxVolumeIndexYZ, VoutFourier);
	release(tempVolume, maxVolumeIndexYZ+1, maxVolumeIndexYZ+1);

    // Output volume
    Image<double> Vout;
    Vout().initZeros(paddedImgSize,paddedImgSize,paddedImgSize);
    FourierTransformer transformerVol;
    transformerVol.setThreadsNumber(noOfCores);
    transformerVol.fReal = &(Vout.data);
    transformerVol.setFourierAlias(VoutFourier);
    transformerVol.recomputePlanR2C();

    transformerVol.inverseFourierTransform();
    transformerVol.clear();
    CenterFFT(Vout(),false);

    // Correct by the Fourier transform of the blob
    Vout().setXmippOrigin();
    Vout().selfWindow(FIRST_XMIPP_INDEX(imgSize),FIRST_XMIPP_INDEX(imgSize),
                      FIRST_XMIPP_INDEX(imgSize),LAST_XMIPP_INDEX(imgSize),
                      LAST_XMIPP_INDEX(imgSize),LAST_XMIPP_INDEX(imgSize));
    double pad_relation= ((double)padding_factor_proj/padding_factor_vol);
    pad_relation = (pad_relation * pad_relation * pad_relation);

    MultidimArray<double> &mVout=Vout();
    double ipad_relation=1.0/pad_relation;
    double meanFactor2=0;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mVout)
    {
        double radius=sqrt((double)(k*k+i*i+j*j));
        double aux=radius*iDeltaFourier;
        double factor = Fourier_blob_table(ROUND(aux));
        double factor2=(pow(Sinc(radius/(2*(imgSize))),2));
		A3D_ELEM(mVout,k,i,j) /= (ipad_relation*factor2*factor);
		meanFactor2+=factor2;
    }
	meanFactor2/=MULTIDIM_SIZE(mVout);
	FOR_ALL_ELEMENTS_IN_ARRAY3D(mVout)
	A3D_ELEM(mVout,k,i,j) *= meanFactor2;
    Vout.write(out_name);
    Vout.clear();
}

void ProgRecFourierGPU::setIO(const FileName &fn_in, const FileName &fn_out)
{
    this->fn_in = fn_in;
    this->fn_out = fn_out;
}
