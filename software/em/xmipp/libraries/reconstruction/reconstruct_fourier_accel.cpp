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

#include "reconstruct_fourier_accel.h"

// Define params
void ProgRecFourierAccel::defineParams()
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
    addExampleLine("   xmipp_reconstruct_fourier_accel  -i reconstruction.sel --sym i3 --weight");
}

// Read arguments ==========================================================
void ProgRecFourierAccel::readParams()
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
    maxResolution = getDoubleParam("--max_resolution");
    useCTF = checkParam("--useCTF");
    isPhaseFlipped = checkParam("--phaseFlipped");
    minCTF = getDoubleParam("--minCTF");
    if (useCTF)
        iTs = 1 / getDoubleParam("--sampling");
    bufferSize = getIntParam("--bufferSize");
}

// Show ====================================================================
void ProgRecFourierAccel::show()
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
void ProgRecFourierAccel::run()
{
    show();
    produceSideinfo();
    // Process all images in the selfile
    if (verbose) {
		init_progress_bar(SF.size());
    }
    // Create loading thread stuff
    createLoadingThread();
    //Computing interpolated volume
    processImages(0, SF.size() - 1);
    // remove complex conjugate of the intermediate result
    mirrorAndCropTempSpaces();
    //Saving the volume
    finishComputations(fn_out);
    cleanLoadingThread();
}

void ProgRecFourierAccel::createLoadingThread() {
	barrier_init( &barrier, 2 ); // two barries - for main and loading thread
	loadThread.parent = this;
	loadThread.selFile = &SF;
	pthread_create( &loadThread.id , NULL, loadImageThread, (void *)(&loadThread) );
	threadOpCode = PRELOAD_IMAGE;
}

void ProgRecFourierAccel::cleanLoadingThread() {
	threadOpCode = EXIT_THREAD;
	// Waiting for thread to finish
	barrier_wait( &barrier );
	pthread_join(*(&loadThread.id), NULL);
	barrier_destroy( &barrier );
}

void ProgRecFourierAccel::produceSideinfo()
{
    // Translate the maximum resolution to digital frequency
    // maxResolution=sampling_rate/maxResolution;
    maxResolutionSqr=maxResolution*maxResolution;

    // Read the input images
    SF.read(fn_in);
    SF.removeDisabled();

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
    double iw0 = 1.0 / blob_Fourier_val(0.0, blobnormalized);
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
}

inline void ProgRecFourierAccel::getVectors(const Point3D* plane, Point3D& u, Point3D& v) {
	float x0 = plane[0].x;
	float y0 = plane[0].y;
	float z0 = plane[0].z;
	u.x = plane[1].x - x0;
	u.y = plane[1].y - y0;
	u.z = plane[1].z - z0;
	v.x = plane[3].x - x0;
	v.y = plane[3].y - y0;
	v.z = plane[3].z - z0;
}

Array2D<std::complex<float> >* ProgRecFourierAccel::cropAndShift(MultidimArray<std::complex<double> >& paddedFourier,
		ProgRecFourierAccel * parent) {
	int sizeX = parent->maxVolumeIndexX / 2; // input Fourier contains just one half of the space, second is complex conjugate
	int sizeY = parent->maxVolumeIndexYZ;

	Array2D<std::complex<float> >* result = new Array2D<std::complex<float> >(sizeX, sizeY);
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
				(*result)(j, myPadI) = std::complex<float>(paddedFourierTmp.real(), paddedFourierTmp.imag());
			}
		}
	}
	return result;
}

void ProgRecFourierAccel::preloadBuffer(LoadThreadParams* threadParams,
		ProgRecFourierAccel* parent,
		bool hasCTF, std::vector<size_t>& objId)
{
    ApplyGeoParams params;
    MultidimArray<double> localPaddedImg;
	FourierTransformer localTransformerImg;
	MultidimArray< std::complex<double> > localPaddedFourier;
    params.only_apply_shifts = true;
	if (0 == threadParams->buffer1) {
		threadParams->buffer1 = new ProjectionData[parent->bufferSize];
	}
	for (int bIndex = 0; bIndex < parent->bufferSize; bIndex++) {
		double rot, tilt, psi, weight;
		Projection proj;
		Matrix2D<double> localA(3, 3);
		int imgIndex = threadParams->startImageIndex + bIndex;
		ProjectionData* data = &threadParams->buffer1[bIndex];
		if (imgIndex >= threadParams->endImageIndex) {
			continue;
		}
		//Read projection from selfile, read also angles and shifts if present
		//but only apply shifts
		proj.readApplyGeo(*(threadParams->selFile), objId[imgIndex], params);
		rot = proj.rot();
		tilt = proj.tilt();
		psi = proj.psi();
		if (parent->do_weights && proj.weight() == 0.f) {
			continue;
		}
		data->weight = (parent->do_weights) ? proj.weight() : 1.0;

		// Copy the projection to the center of the padded image
		// and compute its Fourier transform
		proj().setXmippOrigin();
		localPaddedImg.initZeros(parent->paddedImgSize, parent->paddedImgSize);
		localPaddedImg.setXmippOrigin();
		const MultidimArray<double> &mProj = proj();
		FOR_ALL_ELEMENTS_IN_ARRAY2D(mProj)
			A2D_ELEM(localPaddedImg,i,j) = A2D_ELEM(mProj, i, j);
		CenterFFT(localPaddedImg, true);

		// Fourier transformer for the images
		localTransformerImg.setReal(localPaddedImg);
		localTransformerImg.FourierTransform();
		localTransformerImg.getFourierAlias(localPaddedFourier);

		// Compute the coordinate axes associated to this image
		Euler_angles2matrix(rot, tilt, psi, localA);

		data->localAInv = localA.transpose();
		data->img = cropAndShift(localPaddedFourier, parent);
		data->imgIndex = imgIndex;
		if (hasCTF) {
			Array2D<float>* CTF = new Array2D<float>(data->img->getXSize(), data->img->getYSize());
			Array2D<float>* modulator = new Array2D<float>(data->img->getXSize(), data->img->getYSize());
			preloadCTF(threadParams, objId[imgIndex],parent, CTF, modulator);
		}
		// set data as usable
		data->skip = false;
		//#define DEBUG22
#ifdef DEBUG22
                //CORRECTO

			if(threadParams->myThreadID%1==0)
			{
				proj.write((std::string) integerToString(threadParams->myThreadID) + "_" +\
 integerToString(threadParams->imageIndex) + "proj.spi");

				ImageXmipp save44;
				save44()=localPaddedImg;
				save44.write((std::string) integerToString(threadParams->myThreadID) + "_" +\
 integerToString(threadParams->imageIndex) + "local_padded_img.spi");

				FourierImage save33;
				save33()=localPaddedFourier;
				save33.write((std::string) integerToString(threadParams->myThreadID) + "_" +\
 integerToString(threadParams->imageIndex) + "local_padded_fourier.spi");
				FourierImage save22;
				//save22()=*paddedFourier;
				save22().alias(*(threadParams->localPaddedFourier));
				save22.write((std::string) integerToString(threadParams->myThreadID) + "_" +\
 integerToString(threadParams->imageIndex) + "_padded_fourier.spi");
			}

#endif
#undef DEBUG22
	}
}

void * ProgRecFourierAccel::loadImageThread( void * threadArgs )
{
    LoadThreadParams * threadParams = (LoadThreadParams *) threadArgs;
    ProgRecFourierAccel * parent = threadParams->parent;
    barrier_t * barrier = &(parent->barrier);

    std::vector<size_t> objId;
    threadParams->selFile->findObjects(objId);
    bool hasCTF = parent->useCTF
    		&& (threadParams->selFile->containsLabel(MDL_CTF_MODEL)
    				|| threadParams->selFile->containsLabel(MDL_CTF_DEFOCUSU));
    do
    {
        barrier_wait( barrier );

        switch ( parent->threadOpCode )
        {
        case PRELOAD_IMAGE:
            {
            	preloadBuffer(threadParams, parent, hasCTF, objId);
                break;
            }
        case EXIT_THREAD:
            return NULL;
        default:
            break;
        }

        barrier_wait( barrier );
    }
    while ( true );
}

template<typename T, typename U>
inline U ProgRecFourierAccel::clamp(U val, T min, T max) {
	U res = val;
	res = (res > max) ? max : res;
	res = (res < min) ? min : res;
	return res;
}

inline void ProgRecFourierAccel::multiply(const float transform[3][3], Point3D& inOut) {
	float tmp0 = transform[0][0] * inOut.x + transform[0][1] * inOut.y + transform[0][2] * inOut.z;
	float tmp1 = transform[1][0] * inOut.x + transform[1][1] * inOut.y + transform[1][2] * inOut.z;
	float tmp2 = transform[2][0] * inOut.x + transform[2][1] * inOut.y + transform[2][2] * inOut.z;
	inOut.x = tmp0;
	inOut.y = tmp1;
	inOut.z = tmp2;
}

void ProgRecFourierAccel::createProjectionCuboid(Point3D* cuboid, float sizeX, float sizeY, float blobSize)
{
	float halfY = sizeY / 2.0f;
	cuboid[0].x = cuboid[3].x = cuboid[4].x = cuboid[7].x = 0.f - blobSize;
	cuboid[1].x = cuboid[2].x = cuboid[5].x = cuboid[6].x = sizeX + blobSize;

	cuboid[0].y = cuboid[1].y = cuboid[4].y = cuboid[5].y = -(halfY + blobSize);
	cuboid[2].y = cuboid[3].y = cuboid[6].y = cuboid[7].y = halfY + blobSize;

	cuboid[0].z = cuboid[1].z = cuboid[2].z = cuboid[3].z = 0.f + blobSize;
	cuboid[4].z = cuboid[5].z = cuboid[6].z = cuboid[7].z = 0.f - blobSize;
}

inline bool ProgRecFourierAccel::getZ(float x, float y, float& z, const Point3D& a, const Point3D& b, const Point3D& p0) {
	// from parametric eq. of the plane
	float x0 = p0.x;
	float y0 = p0.y;
	float z0 = p0.z;

	float u = ((y-y0)*a.x + (x0-x)*a.y) / (a.x * b.y - b.x * a.y);
	float t = (-x0 + x - u*b.x) / (a.x);

	z = z0 + t*a.z + u*b.z;
	return inRange(t, 0.f, 1.f) && inRange(u, 0.f, 1.f);
}

inline bool ProgRecFourierAccel::getY(float x, float& y, float z, const Point3D& a, const Point3D& b, const Point3D& p0) {
	// from parametric eq. of the plane
	float x0 = p0.x;
	float y0 = p0.y;
	float z0 = p0.z;

	float u = ((z-z0)*a.x + (x0-x)*a.z) / (a.x * b.z - b.x * a.z);
	float t = (-x0 + x - u*b.x) / (a.x);

	y = y0 + t*a.y + u*b.y;
	return inRange(t, 0.f, 1.f) && inRange(u, 0.f, 1.f);
}

inline bool ProgRecFourierAccel::getX(float& x, float y, float z, const Point3D& a, const Point3D& b, const Point3D& p0) {
	// from parametric eq. of the plane
	float x0 = p0.x;
	float y0 = p0.y;
	float z0 = p0.z;

	float u = ((z-z0)*a.y + (y0-y)*a.z) / (a.y * b.z - b.y * a.z);
	float t = (-y0 + y - u*b.y) / (a.y);

	x = x0 + t*a.x + u*b.x;
	return inRange(t, 0.f, 1.f) && inRange(u, 0.f, 1.f);
}

inline void ProgRecFourierAccel::translateCuboid(Point3D* cuboid, Point3D vector) {
	for (int i = 0; i < 8; i++) {
		cuboid[i].x += vector.x;
		cuboid[i].y += vector.y;
		cuboid[i].z += vector.z;
	}
}

void ProgRecFourierAccel::computeAABB(Point3D* AABB, Point3D* cuboid,
	float minX, float minY, float minZ,
	float maxX, float maxY, float maxZ) {
	AABB[0].x = AABB[0].y = AABB[0].z = std::numeric_limits<float>::max();
	AABB[1].x = AABB[1].y = AABB[1].z = std::numeric_limits<float>::min();
	Point3D tmp;
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

void ProgRecFourierAccel::printAABB(Point3D* AABB) {
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

inline void ProgRecFourierAccel::preloadCTF(LoadThreadParams* threadParams,
		size_t imgIndex,
		ProgRecFourierAccel* parent,
		Array2D<float>* CTF,
		Array2D<float>* modulator)
{
	CTFDescription ctf;
	ctf.readFromMetadataRow(*(threadParams->selFile), imgIndex);
	ctf.produceSideInfo();
	float freqX, freqY;
	float CTFVal, modulatorVal;
	for (int y = 0; y < CTF->getYSize(); y++) {
		// since Y axis is shifted to center, we have to use different calculation
		freqY = (y - (parent->paddedImgSize / 2.f)) / (float) parent->paddedImgSize;
		for (int x = 0; x < CTF->getXSize(); x++) {
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

			(*CTF)(x, y) = CTFVal;
			(*modulator)(x, y) = modulatorVal;
		}
	}



}


inline void ProgRecFourierAccel::processVoxel(int x, int y, int z, const float transform[3][3], float maxDistanceSqr,
		ProjectionData* const data) {
	Point3D imgPos;
	float wBlob = 1.f;
	float wCTF = 1.f;
	float wModulator = 1.f;

	// transform current point to center
	imgPos.x = x - maxVolumeIndexX/2;
	imgPos.y = y - maxVolumeIndexYZ/2;
	imgPos.z = z - maxVolumeIndexYZ/2;
	if (imgPos.x*imgPos.x + imgPos.y*imgPos.y + imgPos.z*imgPos.z > maxDistanceSqr) {
		return; // discard iterations that would access pixel with too high frequency
	}
	// rotate around center
	multiply(transform, imgPos);
	// transform back and round
	// just Y coordinate needs adjusting, since X now matches to picture and Z is irrelevant
	int imgX = clamp((int)(imgPos.x + 0.5f), 0, data->img->getXSize() - 1);
	int imgY = clamp((int)(imgPos.y + 0.5f + maxVolumeIndexYZ / 2), 0, data->img->getYSize() - 1);

	if (0 != data->CTF) {
		wCTF = (*data->CTF)(imgX, imgY);
		wModulator = (*data->modulator)(imgX, imgY);
	}

	float weight = wBlob * wModulator * data->weight;

	tempVolume[z][y][x] += (*data->img)(imgX, imgY) * weight * wCTF;
	tempWeights[z][y][x] += weight;
}

inline void ProgRecFourierAccel::processVoxelBlob(int x, int y, int z, const float transform[3][3], float maxDistanceSqr,
		ProjectionData* const data) {
	Point3D imgPos;
	// transform current point to center
	imgPos.x = x - maxVolumeIndexX/2;
	imgPos.y = y - maxVolumeIndexYZ/2;
	imgPos.z = z - maxVolumeIndexYZ/2;
	if ((imgPos.x*imgPos.x + imgPos.y*imgPos.y + imgPos.z*imgPos.z) > maxDistanceSqr) {
		return; // discard iterations that would access pixel with too high frequency
	}
	// rotate around center
	multiply(transform, imgPos);
	// transform back just Y coordinate, since X now matches to picture and Z is irrelevant
	imgPos.y += maxVolumeIndexYZ / 2;

	// check that we don't want to collect data from far far away ...
	float radiusSqr = blob.radius * blob.radius;
	float zSqr = imgPos.z * imgPos.z;
	if (zSqr > radiusSqr) return;

	// create blob bounding box
	int minX = std::ceil(imgPos.x - blob.radius);
	int maxX = std::floor(imgPos.x + blob.radius);
	int minY = std::ceil(imgPos.y - blob.radius);
	int maxY = std::floor(imgPos.y + blob.radius);
	minX = std::max(minX, 0);
	minY = std::max(minY, 0);
	maxX = std::min(maxX, data->img->getXSize()-1);
	maxY = std::min(maxY, data->img->getYSize()-1);
	std::complex<float>* targetVolume = &tempVolume[z][y][x];
	float* targetWeight = &tempWeights[z][y][x];
	// ugly spaghetti code, but improves performance by app. 10%
	if (0 != data->CTF) {
		// check which pixel in the vicinity that should contribute
		for (int i = minY; i <= maxY; i++) {
			float ySqr = (imgPos.y - i) * (imgPos.y - i);
			float yzSqr = ySqr + zSqr;
			if (yzSqr > radiusSqr) continue;
			for (int j = minX; j <= maxX; j++) {
				float xD = imgPos.x - j;
				float distanceSqr = xD*xD + yzSqr;
				if (distanceSqr > radiusSqr) continue;

				float wCTF = (*data->CTF)(j, i);
				float wModulator = (*data->modulator)(j, i);
				int aux = (int) ((distanceSqr * iDeltaSqrt + 0.5f)); //Same as ROUND but avoid comparison
				float wBlob = blobTableSqrt[aux];
				float weight = wBlob * wModulator * data->weight;
				*targetWeight += weight;
				*targetVolume += (*data->img)(j, i) * weight * wCTF;
			}
		}
	} else {
		// check which pixel in the vicinity that should contribute
		for (int i = minY; i <= maxY; i++) {
			float ySqr = (imgPos.y - i) * (imgPos.y - i);
			float yzSqr = ySqr + zSqr;
			if (yzSqr > radiusSqr) continue;
			for (int j = minX; j <= maxX; j++) {
				float xD = imgPos.x - j;
				float distanceSqr = xD*xD + yzSqr;
				if (distanceSqr > radiusSqr) continue;

				int aux = (int) ((distanceSqr * iDeltaSqrt + 0.5f)); //Same as ROUND but avoid comparison
				float wBlob = blobTableSqrt[aux];

				float weight = wBlob * data->weight;
				*targetWeight += weight;
				*targetVolume += (*data->img)(j, i) * weight;
			}
		}
	}

}

inline void ProgRecFourierAccel::convert(Matrix2D<double>& in, float out[3][3]) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			out[i][j] = in(i, j);
		}
	}
}

void ProgRecFourierAccel::processProjection(
	ProjectionData* projectionData,
	const float transform[3][3],
	const float transformInv[3][3])
{
	int imgSizeX = projectionData->img->getXSize();
	int imgSizeY = projectionData->img->getYSize();
	const float maxDistanceSqr = (imgSizeX+(useFast ? 0.f : blob.radius)) * (imgSizeX+(useFast ? 0.f : blob.radius));
	Point3D origin = {maxVolumeIndexX/2.f, maxVolumeIndexYZ/2.f, maxVolumeIndexYZ/2.f};
	Point3D u, v, normal;
	Point3D AABB[2];
	Point3D cuboid[8];

	// calculate affected space
	createProjectionCuboid(cuboid, imgSizeX, imgSizeY, useFast ? 0.f : blob.radius);
	rotateCuboid(cuboid, transform);
	translateCuboid(cuboid, origin);
	computeAABB(AABB, cuboid, 0, 0, 0, maxVolumeIndexX, maxVolumeIndexYZ, maxVolumeIndexYZ);
	getVectors(cuboid, u, v);
	normal = getNormal(u, v);

	// prepare traversing
	int minY, minX, minZ;
	int maxY, maxX, maxZ;
	minZ = floor(AABB[0].z);
	minY = floor(AABB[0].y);
	minX = floor(AABB[0].x);
	maxZ = ceil(AABB[1].z);
	maxY = ceil(AABB[1].y);
	maxX = ceil(AABB[1].x);
	// iterate along the longest axes, because it has the shortest projection to traverse plane
	int nX, nY, nZ;
	nX = std::abs(nX);
	nY = std::abs(nY);
	nZ = std::abs(nZ);

	if (nZ >= nX && nZ >= nY) { // iterate XY plane
		for(int y = minY; y <= maxY; y++) {
			for(int x = minX; x <= maxX; x++) {
				if (useFast) {
					float hitZ;
					if (getZ(x, y, hitZ, u, v, *cuboid)) {
						int z = (int)(hitZ + 0.5f); // rounding
						processVoxel(x, y, z, transformInv, maxDistanceSqr, projectionData);
					}
				} else {
					float z1, z2;
					bool hit1 = getZ(x, y, z1, u, v, *cuboid); // lower plane
					bool hit2 = getZ(x, y, z2, u, v, *(cuboid + 4)); // upper plane
					if (hit1 || hit2) {
						z1 = clamp(z1, 0, maxVolumeIndexYZ);
						z2 = clamp(z2, 0, maxVolumeIndexYZ);
						float lower = std::min(z1, z2);
						float upper = std::max(z1, z2);
						for (int z = std::floor(lower); z <= std::ceil(upper); z++) {
							processVoxelBlob(x, y, z, transformInv, maxDistanceSqr, projectionData);
						}
					}
				}
			}
		}
	} else if (nY >= nX && nY >= nZ) { // iterate XZ plane
		for(int z = minZ; z <= maxZ; z++) {
			for(int x = minX; x <= maxX; x++) {
				if (useFast) {
					float hitY;
					if (getY(x, hitY, z, u, v, *cuboid)) {
						int y = (int)(hitY + 0.5f); // rounding
						processVoxel(x, y, z, transformInv, maxDistanceSqr, projectionData);
					}
				} else {
					float y1, y2;
					bool hit1 = getY(x, y1, z, u, v, *cuboid); // lower plane
					bool hit2 = getY(x, y2, z, u, v, *(cuboid + 4)); // upper plane
					if (hit1 || hit2) {
						y1 = clamp(y1, 0, maxVolumeIndexYZ);
						y2 = clamp(y2, 0, maxVolumeIndexYZ);
						float lower = std::min(y1, y2);
						float upper = std::max(y1, y2);
						for (int y = std::floor(lower); y <= std::ceil(upper); y++) {
							processVoxelBlob(x, y, z, transformInv, maxDistanceSqr, projectionData);
						}
					}
				}
			}
		}
	} else if(nX >= nY  && nX >= nZ) { // iterate YZ plane
		for(int z = minZ; z <= maxZ; z++) {
			for(int y = minY; y <= maxY; y++) {
				if (useFast) {
					float hitX;
					if (getX(hitX, y, z, u, v, *cuboid)) {
						int x = (int)(hitX + 0.5f); // rounding
						processVoxel(x, y, z, transformInv, maxDistanceSqr, projectionData);
					}
				} else {
					float x1, x2;
					bool hit1 = getX(x1, y, z, u, v, *cuboid); // lower plane
					bool hit2 = getX(x2, y, z, u, v, *(cuboid + 4)); // upper plane
					if (hit1 || hit2) {
						x1 = clamp(x1, 0, maxVolumeIndexX);
						x2 = clamp(x2, 0, maxVolumeIndexX);
						float lower = std::min(x1, x2);
						float upper = std::max(x1, x2);
						for (int x = std::floor(lower); x <= std::ceil(upper); x++) {
							processVoxelBlob(x, y, z, transformInv, maxDistanceSqr, projectionData);
						}
					}
				}
			}
		}
	}
}

template<typename T>
T*** ProgRecFourierAccel::allocate(T***& where, int xSize, int ySize, int zSize) {
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
void ProgRecFourierAccel::release(T***& array, int ySize, int zSize) {
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
T*** ProgRecFourierAccel::applyBlob(T***& input, float blobSize,
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
							int aux = (int) ((distanceSqr * iDeltaSqrt + 0.5f)); //Same as ROUND but avoid comparison
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
void ProgRecFourierAccel::convertToExpectedSpace(T*** input, int size,
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

void ProgRecFourierAccel::mirrorAndCropTempSpaces() {
	maxVolumeIndexX = maxVolumeIndexYZ/2; // just half of the space is necessary, the rest is complex conjugate
	mirrorAndCrop(tempWeights, &identity<float>);
	mirrorAndCrop(tempVolume, &conjugate);
}

template<typename T>
void ProgRecFourierAccel::mirrorAndCrop(T***& input, T (*f)(T)) {

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

void ProgRecFourierAccel::forceHermitianSymmetry() {
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

void ProgRecFourierAccel::processWeights() {
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

void ProgRecFourierAccel::loadImages(int startIndex, int endIndex) {
	loadThread.startImageIndex = startIndex;
	loadThread.endImageIndex = endIndex;
	// Awaking sleeping threads
	barrier_wait( &barrier );
}

void ProgRecFourierAccel::swapLoadBuffers() {
	ProjectionData* tmp = loadThread.buffer2;
	loadThread.buffer2 = loadThread.buffer1;
	loadThread.buffer1 = tmp;
}

void ProgRecFourierAccel::processBuffer(ProjectionData* buffer)
{
	int repaint = (int)ceil((double)SF.size()/60);
	for ( int i = 0 ; i < bufferSize; i++ ) {
		ProjectionData* projData = &buffer[i];
		Array2D<std::complex<float> >* myPaddedFourier = projData->img;
		if (projData->skip) {
			continue;
		}
		if (verbose && projData->imgIndex%repaint==0) {
			progress_bar(projData->imgIndex);
		}

		Matrix2D<double> *Ainv = &projData->localAInv;
		// Loop over all symmetries
		for (size_t isym = 0; isym < R_repository.size(); isym++)
		{
			// Compute the coordinate axes of the symmetrized projection
			Matrix2D<double> A_SL=R_repository[isym]*(*Ainv);
			Matrix2D<double> A_SLInv=A_SL.inv();
			float transf[3][3];
			float transfInv[3][3];
			convert(A_SL, transf);
			convert(A_SLInv, transfInv);
			processProjection(//tempVolume, tempWeights, size,
					projData, transf, transfInv);
		}
		projData->clean();
	}
}

void ProgRecFourierAccel::processImages( int firstImageIndex, int lastImageIndex)
{
    int loops = ceil((lastImageIndex-firstImageIndex+1)/(float)bufferSize);

	// the +1 is to prevent outOfBound reading when mirroring the result (later)
    if (NULL == tempVolume) {
    	allocate(tempVolume, maxVolumeIndexYZ+1, maxVolumeIndexYZ+1, maxVolumeIndexYZ+1);
    }
    if (NULL == tempWeights) {
    	allocate(tempWeights, maxVolumeIndexYZ+1, maxVolumeIndexYZ+1, maxVolumeIndexYZ+1);
    }

    int startLoadIndex = firstImageIndex;

    loadImages(startLoadIndex, std::min(lastImageIndex+1, startLoadIndex+bufferSize));
	barrier_wait( &barrier );
    for(int i = 0; i < loops; i++) {
    	swapLoadBuffers();
    	startLoadIndex += bufferSize;
    	loadImages(startLoadIndex, std::min(lastImageIndex+1, startLoadIndex+bufferSize));
    	processBuffer(loadThread.buffer2);
    	barrier_wait( &barrier );
    }
	delete[] loadThread.buffer1;
	delete[] loadThread.buffer2;
	loadThread.buffer1 = loadThread.buffer2 = NULL;
}

void ProgRecFourierAccel::releaseTempSpaces() {
	release(tempWeights, maxVolumeIndexYZ+1, maxVolumeIndexYZ+1);
	release(tempVolume, maxVolumeIndexYZ+1, maxVolumeIndexYZ+1);
}

void ProgRecFourierAccel::finishComputations( const FileName &out_name )
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
    transformerVol.setThreadsNumber(2); // use just main and 'loading' thread
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

void ProgRecFourierAccel::setIO(const FileName &fn_in, const FileName &fn_out)
{
    this->fn_in = fn_in;
    this->fn_out = fn_out;
}
