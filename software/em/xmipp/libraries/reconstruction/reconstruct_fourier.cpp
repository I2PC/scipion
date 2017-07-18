/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
 *              Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Jose Roman Bilbao-Castro (jrbcast@ace.ual.es)
 *              Vahid Abrishami (vabrishami@cnb.csic.es)
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

#include "reconstruct_fourier.h"
#include <fstream>
#include <set>
#include <sstream>
#include <math.h>
#include <limits>
#include <data/xmipp_fft.h>

#include <time.h>

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()


#define DEBUG_DUMP 0

// Define params
void ProgRecFourier::defineParams()
{
    //usage
    addUsageLine("Generate 3D reconstructions from projections using direct Fourier interpolation with arbitrary geometry.");
    addUsageLine("Kaisser-windows are used for interpolation in Fourier space.");
    //params
    addParamsLine("   -i <md_file>                : Metadata file with input projections");
    addParamsLine("  [-o <volume_file=\"rec_fourier.vol\">]  : Filename for output volume");
    addParamsLine("  [--iter <iterations=1>]      : Number of iterations for weight correction");
    addParamsLine("  [--sym <symfile=c1>]              : Enforce symmetry in projections");
    addParamsLine("  [--padding <proj=2.0> <vol=2.0>]  : Padding used for projections and volume");
    addParamsLine("  [--prepare_fsc <fscfile>]      : Filename root for FSC files");
    addParamsLine("  [--max_resolution <p=0.5>]     : Max resolution (Nyquist=0.5)");
    addParamsLine("  [--weight]                     : Use weights stored in the image metadata");
    addParamsLine("  [--thr <threads=1> <rows=1>]   : Number of concurrent threads and rows processed at time by a thread");
    addParamsLine("  [--blob <radius=1.9> <order=0> <alpha=15>] : Blob parameters");
    addParamsLine("                                 : radius in pixels, order of Bessel function in blob and parameter alpha");
    addParamsLine("  [--useCTF]                     : Use CTF information if present");
    addParamsLine("  [--sampling <Ts=1>]            : sampling rate of the input images in Angstroms/pixel");
    addParamsLine("                                 : It is only used when correcting for the CTF");
    addParamsLine("  [--phaseFlipped]               : Give this flag if images have been already phase flipped");
    addParamsLine("  [--minCTF <ctf=0.01>]          : Minimum value of the CTF that will be inverted");
    addParamsLine("                                 : CTF values (in absolute value) below this one will not be corrected");
    addExampleLine("For reconstruct enforcing i3 symmetry and using stored weights:", false);
    addExampleLine("   xmipp_reconstruct_fourier  -i reconstruction.sel --sym i3 --weight");
}

// Read arguments ==========================================================
void ProgRecFourier::readParams()
{
    fn_sel = getParam("-i");
    fn_out = getParam("-o");
    fn_sym = getParam("--sym");
    if(checkParam("--prepare_fsc"))
        fn_fsc = getParam("--prepare_fsc");
    do_weights = checkParam("--weight");
    padding_factor_proj = getDoubleParam("--padding", 0);
    padding_factor_vol = getDoubleParam("--padding", 1);
    blob.radius   = getDoubleParam("--blob", 0);
    blob.order    = getIntParam("--blob", 1);
    blob.alpha    = getDoubleParam("--blob", 2);
    maxResolution = getDoubleParam("--max_resolution");
    numThreads = getIntParam("--thr");
    thrWidth = getIntParam("--thr", 1);
    NiterWeight = getIntParam("--iter");
    useCTF = checkParam("--useCTF");
    phaseFlipped = checkParam("--phaseFlipped");
    minCTF = getDoubleParam("--minCTF");
    if (useCTF)
        Ts=getDoubleParam("--sampling");
}

// Show ====================================================================
void ProgRecFourier::show()
{
    if (verbose > 0)
    {
        std::cout << " =====================================================================" << std::endl;
        std::cout << " Direct 3D reconstruction method using Kaiser windows as interpolators" << std::endl;
        std::cout << " =====================================================================" << std::endl;
        std::cout << " Input selfile             : "  << fn_sel << std::endl;
        std::cout << " padding_factor_proj       : "  << padding_factor_proj << std::endl;
        std::cout << " padding_factor_vol        : "  << padding_factor_vol << std::endl;
        std::cout << " Output volume             : "  << fn_out << std::endl;
        if (fn_sym != "")
            std::cout << " Symmetry file for projections : "  << fn_sym << std::endl;
        if (fn_fsc != "")
            std::cout << " File root for FSC files: " << fn_fsc << std::endl;
        if (do_weights)
            std::cout << " Use weights stored in the image headers or doc file" << std::endl;
        else
            std::cout << " Do NOT use weights" << std::endl;
        if (useCTF)
            std::cout << "Using CTF information" << std::endl
            << "Sampling rate: " << Ts << std::endl
            << "Phase flipped: " << phaseFlipped << std::endl
            << "Minimum CTF: " << minCTF << std::endl;
        std::cout << "\n Interpolation Function"
        << "\n   blrad                 : "  << blob.radius
        << "\n   blord                 : "  << blob.order
        << "\n   blalpha               : "  << blob.alpha
        //<< "\n sampling_rate           : "  << sampling_rate
        << "\n max_resolution          : "  << maxResolution
        << "\n -----------------------------------------------------------------" << std::endl;
    }
}

// Main routine ------------------------------------------------------------
void ProgRecFourier::run()
{
    show();
    produceSideinfo();
    // Process all images in the selfile
    if (verbose)
    {
        if (NiterWeight!=0)
            init_progress_bar(NiterWeight*SF.size());
        else
            init_progress_bar(SF.size());
    }
    // Create threads stuff
    barrier_init( &barrier, numThreads+1 );
    pthread_mutex_init( &workLoadMutex, NULL );
    statusArray = NULL;
    th_ids = (pthread_t *)malloc( numThreads * sizeof( pthread_t));
    th_args = (ImageThreadParams *) malloc ( numThreads * sizeof( ImageThreadParams ) );

    // Create threads
    for ( int nt = 0 ; nt < numThreads ; nt ++ )
    {
        // Passing parameters to each thread
        th_args[nt].parent = this;
        th_args[nt].myThreadID = nt;
        th_args[nt].selFile = new MetaData(SF);
        pthread_create( (th_ids+nt) , NULL, processImageThread, (void *)(th_args+nt) );
    }

    //Computing interpolated volume
    processImages(0, SF.size() - 1, !fn_fsc.empty(), false);

    // Correcting the weights
    correctWeight();

    //Saving the volume
    finishComputations(fn_out);

    threadOpCode = EXIT_THREAD;

    // Waiting for threads to finish
    barrier_wait( &barrier );
    for ( int nt = 0 ; nt < numThreads ; nt ++ )
        pthread_join(*(th_ids+nt), NULL);
    barrier_destroy( &barrier );
}


void ProgRecFourier::produceSideinfo()
{
    // Translate the maximum resolution to digital frequency
    // maxResolution=sampling_rate/maxResolution;
    maxResolution2=maxResolution*maxResolution;

    // Read the input images
    SF.read(fn_sel);
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
    volPadSizeX = volPadSizeY = volPadSizeZ=(int)(Xdim*padding_factor_vol);
    Vout().initZeros(volPadSizeZ,volPadSizeY,volPadSizeX);

    //use threads for volume inverse fourier transform, plan is created in setReal()
    transformerVol.setThreadsNumber(numThreads);
    transformerVol.setReal(Vout());

    Vout().clear(); // Free the memory so that it is available for FourierWeights
    transformerVol.getFourierAlias(VoutFourier);
    VoutFourier.initZeros();
    FourierWeights.initZeros(VoutFourier);


    // Ask for memory for the padded images
    size_t paddedImgSize=(size_t)(Xdim*padding_factor_proj);
    paddedImg.resize(paddedImgSize,paddedImgSize);
    paddedImg.setXmippOrigin();
    transformerImg.setReal(paddedImg);

    // Build a table of blob values
    blobTableSqrt.resize(BLOB_TABLE_SIZE_SQRT);
    fourierBlobTableSqrt.resize(BLOB_TABLE_SIZE_SQRT);
    Fourier_blob_table.resize(BLOB_TABLE_SIZE_SQRT);

    struct blobtype blobFourier,blobnormalized;
    blobFourier=blob;
    //Sjors 18aug10 blobFourier.radius/=(padding_factor_proj*Xdim);
    blobFourier.radius/=(padding_factor_vol*Xdim);
    blobnormalized=blob;
    blobnormalized.radius/=((double)padding_factor_proj/padding_factor_vol);
    double deltaSqrt     = (blob.radius*blob.radius) /(BLOB_TABLE_SIZE_SQRT-1);
    double deltaFourier  = (sqrt(3.)*Xdim/2.)/(BLOB_TABLE_SIZE_SQRT-1);

    // The interpolation kernel must integrate to 1
    double iw0 = 1.0 / blob_Fourier_val(0.0, blobnormalized);
    //Sjors 18aug10 double padXdim3 = padding_factor_proj * Xdim;
    double padXdim3 = padding_factor_vol * Xdim;
    padXdim3 = padXdim3 * padXdim3 * padXdim3;
    double blobTableSize = blob.radius*sqrt(1./ (BLOB_TABLE_SIZE_SQRT-1));
    //***
    //Following commented line seems to be the right thing but I do not understand it
    //double fourierBlobTableSize = (sqrt(3.)*Xdim*Xdim/2.)*blobFourier.radius *sqrt(1./ (BLOB_TABLE_SIZE_SQRT-1));
    FOR_ALL_ELEMENTS_IN_MATRIX1D(blobTableSqrt)

    {
        //use a r*r sample instead of r
        //DIRECT_VEC_ELEM(blob_table,i)         = blob_val(delta*i, blob)  *iw0;
        VEC_ELEM(blobTableSqrt,i)    = blob_val(blobTableSize*sqrt((double)i), blob)  *iw0;
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

void to3D( int i ) {
	int z = i % 720;
	int y = (i / 720) % 720;
	int x = i / (720 * 720);
    std::cout << x << " " << y << " " << z;
}

template <class T>
Matrix1D<T> revert(T i1, T i2, T i3, Matrix2D<T>& A_SL, int volPadSizeX, int volPadSizeY, int volPadSizeZ, int xdim, int ydim) {
	Matrix2D<T> A_SL_inv = A_SL.inv();
	Matrix1D<T> real_pos(3), img_pos(3);
	FFT_IDX2DIGFREQ_DOUBLE(i1, volPadSizeX, real_pos[0]);
	FFT_IDX2DIGFREQ_DOUBLE(i2, volPadSizeY, real_pos[1]);
	FFT_IDX2DIGFREQ_DOUBLE(i3, volPadSizeZ, real_pos[2]);
	// invert map
	SPEED_UP_temps012;
	M3x3_BY_V3x1(real_pos, A_SL_inv, real_pos);
	DIGFREQ2FFT_IDX_DOUBLE(real_pos[0], xdim, img_pos[0]);
	DIGFREQ2FFT_IDX_DOUBLE(real_pos[1], ydim, img_pos[1]);

	img_pos[2] = real_pos[2];

	return img_pos;
}

template <typename T>
bool ellipseHitTest(T majorAxis, T minorAxis, T x0, T y0, T xp, T yp, T alpha)
{
	T tmp1 = cos(alpha) * (xp - x0) + sin(alpha) * (yp-y0);
	T tmp2 = sin(alpha) * (xp - x0) - cos(alpha) * (yp-y0);

	T distance = (tmp1 * tmp1)/(majorAxis * majorAxis) + (tmp2 * tmp2) / (minorAxis * minorAxis);

	return distance <= 1;
}

void ProgRecFourier::processCube(
		int i,
		int j,
		const Matrix1D<int>& corner1,
		const Matrix1D<int>& corner2,
		const MultidimArray<float>& z2precalculated,
		const MultidimArray<int>& zWrapped,
		const MultidimArray<int>& zNegWrapped,
		const MultidimArray<float>& y2precalculated, float blobRadiusSquared,
		const MultidimArray<int>& yWrapped,
		const MultidimArray<int>& yNegWrapped,
		const MultidimArray<float>& x2precalculated, float iDeltaSqrt,
		float wModulator, const MultidimArray<int>& xWrapped, int xsize_1,
		const MultidimArray<int>& xNegWrapped, bool reprocessFlag, float wCTF,
		MultidimArray<std::complex<double> >& VoutFourier,
		Matrix1D<double>& blobTableSqrt, ImageThreadParams* threadParams,
		MultidimArray<double>& fourierWeights, double* ptrIn,
		float weight,
		ProgRecFourier * parent,
		Matrix1D<double>& real_position) {
	// Actually compute
	for (int intz = corner1[2]; intz <= corner2[2]; ++intz) {
		float z2 = z2precalculated(intz);
		int iz = zWrapped(intz);
		int izneg = zNegWrapped(intz);
		for (int inty = corner1[1]; inty <= corner2[1]; ++inty) {
			float y2z2 = y2precalculated(inty) + z2;
			if (y2z2 > blobRadiusSquared)
				continue;

			int iy = yWrapped(inty);
			int iyneg = yNegWrapped(inty);
			int size1 = VoutFourier.yxdim * (izneg)
					+ ((iyneg) * VoutFourier.xdim);
			int size2 = VoutFourier.yxdim * (iz) + ((iy) * VoutFourier.xdim);
			int fixSize = 0;
			for (int intx = corner1[0]; intx <= corner2[0]; ++intx) {
				// Compute distance to the center of the blob
				// Compute blob value at that distance
				float d2 = x2precalculated(intx) + y2z2;
				if (d2 > blobRadiusSquared)
					continue;

				int aux = (int) ((d2 * iDeltaSqrt + 0.5)); //Same as ROUND but avoid comparison
				float w = blobTableSqrt[aux] * weight
						* wModulator;
				int ix = xWrapped(intx);
				bool conjugate = false;
				int izp, iyp, ixp;
				if (ix > xsize_1) {
					izp = izneg;
					iyp = iyneg;
					ixp = xNegWrapped(intx);
					conjugate = true;
					fixSize = size1;
				} else {
					izp = iz;
					iyp = iy;
					ixp = ix;
					fixSize = size2;
				}


				// Add the weighted coefficient
				if (reprocessFlag) {
					// Use VoutFourier as temporary to save the memory
					double* ptrOut = (double*) (&(DIRECT_A3D_ELEM(VoutFourier,
							izp, iyp, ixp)));
					DIRECT_A3D_ELEM(fourierWeights, izp,iyp,ixp) += (w
							* ptrOut[0]);
				} else {
					float wEffective = w * wCTF;
					size_t memIdx = fixSize + ixp; //YXSIZE(VoutFourier)*(izp)+((iyp)*XSIZE(VoutFourier))+(ixp);


//						double tmp[3];
//						FFT_IDX2DIGFREQ(intx, parent->volPadSizeX,
//								tmp[0]);
//						FFT_IDX2DIGFREQ(inty, parent->volPadSizeX,
//								tmp[1]);
//						FFT_IDX2DIGFREQ(intz, parent->volPadSizeX,
//								tmp[2]);

//					if (iyp > 64) {
						std::cout
							<< ixp
							<< " "
							<< iyp
							<< " "
							<< izp
							<< " -> "
//							<< tmp[0]
//							<< " "
//							<< tmp[1]
//							<< " "
//							<< tmp[2]
//							   << " -> "
							<< j
							<< " "
							<< ((i < 64) ?	i + 64 : i - 128 + 64)
//							<< " "
//							<< wEffective
//							<< " "
//							<< (conjugate ? "conjugate" : "no_conjugate")
							<< std::endl;
//					}
					double* ptrOut = (double*) (&(VoutFourier[memIdx]));
					ptrOut[0] += wEffective * ptrIn[0];
					fourierWeights[memIdx] += w;
					if (conjugate)
						ptrOut[1] -= wEffective * ptrIn[1];
					else
						ptrOut[1] += wEffective * ptrIn[1];
				}
			}
		}
	}
}
struct Point3D {
	float x, y, z;
};


inline bool exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

inline Point3D getNormal(const Point3D& u, const Point3D& v) {
	float x = u.y*v.z - u.z*v.y;
	float y = u.z*v.x - u.x*v.z;
	float z = u.x*v.y - u.y*v.x;
	return Point3D {x, y, z};
}

inline float getDot(const Point3D&u, const Point3D& v) {
	return u.x*v.x + u.y*v.y + u.z*v.z;
}

inline float getLength(const Point3D& v) {
	return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

inline float getCosTheta(const Point3D& u, const Point3D& v) {
	return getDot(u, v) / (getLength(u)*getLength(v));
}

inline void getVectors(const Point3D* plane, Point3D& u, Point3D& v) {
	float x0 = plane[0].x;
	float y0 = plane[0].y;
	float z0 = plane[0].z;
	u = { plane[1].x - x0, plane[1].y - y0, plane[1].z - z0 };
	v = { plane[3].x - x0, plane[3].y - y0, plane[3].z - z0 };
}

template <typename T>
void print(MultidimArray<T>& VoutFourier, bool full = true) {
	std::string fileName = "VoutFourier";
	std::string ext = ".txt";
	int index = 0;
	std::string fullName;
	do {
	   fullName = fileName + SSTR(++index) + ext;
	}
	while (exists(fullName));

	std::cout << "about to write " << fullName << std::endl;
	std::ofstream myfile;
	myfile.open(fullName.c_str());
	myfile << std::fixed;
	myfile.precision(7);
	myfile << std::showpos;
	for (int z = 0; z < (full ? VoutFourier.zdim : 1); z++) {
//		myfile << "Slice " << z << std::endl;
		for (int y = 0; y < (full ? VoutFourier.ydim : 100); y++) {
//			myfile << "Row " << y << ": ";
			for (int x = 0; x < VoutFourier.xdim; x++) {
				T data = VoutFourier.data[(z * VoutFourier.xdim * VoutFourier.ydim) + (y * VoutFourier.xdim) + x];
				if (data == (T)0){
					continue;
				}
				myfile << x << " " << y << " " << z << " -> " << data << std::endl;
			}
//			myfile << std::endl;
		}
//		myfile << std::endl;
	}
	myfile.close();
}

void * ProgRecFourier::processImageThread( void * threadArgs )
{

    ImageThreadParams * threadParams = (ImageThreadParams *) threadArgs;
    ProgRecFourier * parent = threadParams->parent;
    barrier_t * barrier = &(parent->barrier);

    int minSeparation;

    if ( (int)ceil(parent->blob.radius) > parent->thrWidth )
        minSeparation = (int)ceil(parent->blob.radius);
    else
        minSeparation = parent->thrWidth;

    minSeparation+=1;

    Matrix2D<double>  localA(3, 3), localAinv(3, 3);
    MultidimArray< std::complex<double> > localPaddedFourier;
    MultidimArray<double> localPaddedImg;
    FourierTransformer localTransformerImg;

    std::vector<size_t> objId;

    threadParams->selFile->findObjects(objId);
    ApplyGeoParams params;
    params.only_apply_shifts = true;
    MultidimArray<int> zWrapped(3*parent->volPadSizeZ),yWrapped(3*parent->volPadSizeY),xWrapped(3*parent->volPadSizeX),
    zNegWrapped, yNegWrapped, xNegWrapped;
    zWrapped.initConstant(-1);
    yWrapped.initConstant(-1);
    xWrapped.initConstant(-1);
    zWrapped.setXmippOrigin();
    yWrapped.setXmippOrigin();
    xWrapped.setXmippOrigin();
    zNegWrapped=zWrapped;
    yNegWrapped=yWrapped;
    xNegWrapped=xWrapped;

    MultidimArray<float> x2precalculated(xWrapped.xdim), y2precalculated(yWrapped.xdim), z2precalculated(zWrapped.xdim);
    x2precalculated.initConstant(-1);
    y2precalculated.initConstant(-1);
    z2precalculated.initConstant(-1);
    x2precalculated.setXmippOrigin();
    y2precalculated.setXmippOrigin();
    z2precalculated.setXmippOrigin();



    bool hasCTF=(threadParams->selFile->containsLabel(MDL_CTF_MODEL) || threadParams->selFile->containsLabel(MDL_CTF_DEFOCUSU)) &&
                parent->useCTF;
    if (hasCTF)
    {
        threadParams->ctf.enable_CTF=true;
        threadParams->ctf.enable_CTFnoise=false;
    }
    do
    {
        barrier_wait( barrier );

        switch ( parent->threadOpCode )
        {
        case PRELOAD_IMAGE:
            {

                threadParams->read = 0;

                if ( threadParams->imageIndex >= 0 )
                {
                    // Read input image
                    double rot, tilt, psi, weight;
                    Projection proj;

                    //Read projection from selfile, read also angles and shifts if present
                    //but only apply shifts

                    proj.readApplyGeo(*(threadParams->selFile), objId[threadParams->imageIndex], params);
                    rot  = proj.rot();
                    tilt = proj.tilt();
                    psi  = proj.psi();
                    weight = proj.weight();
                    if (hasCTF)
                    {
                        threadParams->ctf.readFromMetadataRow(*(threadParams->selFile),objId[threadParams->imageIndex]);
                        threadParams->ctf.Tm=threadParams->parent->Ts;
                        threadParams->ctf.produceSideInfo();
                    }

                    threadParams->weight = 1.;

                    if(parent->do_weights)
                        threadParams->weight = weight;
                    else if (!parent->do_weights)
                    {
                        weight=1.0;
                    }
                    else if (weight==0.0)
                    {
                        threadParams->read = 2;
                        break;
                    }

                    // Copy the projection to the center of the padded image
                    // and compute its Fourier transform
                    proj().setXmippOrigin();
                    size_t localPaddedImgSize=(size_t)(parent->imgSize*parent->padding_factor_proj);
                    if (threadParams->reprocessFlag)
                        localPaddedFourier.initZeros(localPaddedImgSize,localPaddedImgSize/2+1);
                    else
                    {
                        localPaddedImg.initZeros(localPaddedImgSize,localPaddedImgSize);
                        localPaddedImg.setXmippOrigin();
                        const MultidimArray<double> &mProj=proj();
                        FOR_ALL_ELEMENTS_IN_ARRAY2D(mProj)
                        A2D_ELEM(localPaddedImg,i,j)=A2D_ELEM(mProj,i,j);
                        // COSS A2D_ELEM(localPaddedImg,i,j)=weight*A2D_ELEM(mProj,i,j);
                        CenterFFT(localPaddedImg,true);

                        /*localPaddedImg.initZeros(6,6);
                            DIRECT_A2D_ELEM(localPaddedImg,1,0)=1;
                        std::cout << localPaddedImg << std::endl;*/

                        // Fourier transformer for the images
                        localTransformerImg.setReal(localPaddedImg);
                        localTransformerImg.FourierTransform();
                        localTransformerImg.getFourierAlias(localPaddedFourier);
                        /*std::cout << localPaddedFourier << std::endl;
                        MultidimArray< std::complex<double> > intermediate;
                        FourierTransform(localPaddedImg,intermediate);
                        std::cout << intermediate << std::endl;
                        CenterFFT(intermediate,true);
                        std::cout << intermediate << std::endl;*/
                    }

                    // Compute the coordinate axes associated to this image
                    Euler_angles2matrix(rot, tilt, psi, localA);
                    localAinv=localA.transpose();

                    threadParams->localweight = weight;
                    threadParams->localAInv = &localAinv;
                    threadParams->localPaddedFourier = &localPaddedFourier;
                    //#define DEBUG22
#ifdef DEBUG22

                    {//CORRECTO

                        if(threadParams->myThreadID%1==0)
                        {
                            proj.write((std::string) integerToString(threadParams->myThreadID)  + "_" +\
                                       integerToString(threadParams->imageIndex) + "proj.spi");

                            ImageXmipp save44;
                            save44()=localPaddedImg;
                            save44.write((std::string) integerToString(threadParams->myThreadID)  + "_" +\
                                         integerToString(threadParams->imageIndex) + "local_padded_img.spi");

                            FourierImage save33;
                            save33()=localPaddedFourier;
                            save33.write((std::string) integerToString(threadParams->myThreadID)  + "_" +\
                                         integerToString(threadParams->imageIndex) + "local_padded_fourier.spi");
                            FourierImage save22;
                            //save22()=*paddedFourier;
                            save22().alias(*(threadParams->localPaddedFourier));
                            save22.write((std::string) integerToString(threadParams->myThreadID)  + "_" +\
                                         integerToString(threadParams->imageIndex) + "_padded_fourier.spi");
                        }

                    }
#endif
                    #undef DEBUG22

                    threadParams->read = 1;
                }
                break;
            }
        case EXIT_THREAD:
            return NULL;
        case PROCESS_WEIGHTS:
            {

                // Get a first approximation of the reconstruction
                double corr2D_3D=pow(parent->padding_factor_proj,2.)/
                                 (parent->imgSize* pow(parent->padding_factor_vol,3.));
                // Divide by Zdim because of the
                // the extra dimension added
                // and padding differences
                MultidimArray<double> &mFourierWeights=parent->FourierWeights;
                for (int k=threadParams->myThreadID; k<=FINISHINGZ(mFourierWeights); k+=parent->numThreads)
                    for (int i=STARTINGY(mFourierWeights); i<=FINISHINGY(mFourierWeights); i++)
                        for (int j=STARTINGX(mFourierWeights); j<=FINISHINGX(mFourierWeights); j++)
                        {
                        	if (parent->NiterWeight==0)
                        		A3D_ELEM(parent->VoutFourier,k,i,j)*=corr2D_3D;
                        	else
                        	{
                        		double weight_kij=A3D_ELEM(mFourierWeights,k,i,j);
                        		if (1.0/weight_kij>ACCURACY)
                        			A3D_ELEM(parent->VoutFourier,k,i,j)*=corr2D_3D*A3D_ELEM(mFourierWeights,k,i,j);
                        		else
                        			A3D_ELEM(parent->VoutFourier,k,i,j)=0;
                        	}
                        }
                break;
            }
        case PROCESS_IMAGE:
            {
            	std::cout << "processing image " + SSTR(threadParams->imageIndex) << " from thread " + SSTR(threadParams->myThreadID) << std::endl; //+ " with symmetry " << *threadParams->symmetry << ": ";
            				MultidimArray< std::complex<double> > *paddedFourier = threadParams->paddedFourier;
            				if (threadParams->weight==0.0)
            					break;
            				bool reprocessFlag = threadParams->reprocessFlag;
            				int * statusArray = parent->statusArray;

            				// Get the inverse of the sampling rate
            				double iTs=1.0/parent->Ts; // The padding factor is not considered here, but later when the indexes
            										   // are converted to digital frequencies


            				Matrix2D<double> * A_SL = threadParams->symmetry;

            				// Loop over all Fourier coefficients in the padded image
            				Matrix1D<double> freq(3), gcurrent(3), real_position(3), contFreq(3);
            				Matrix1D<int> corner1(3), corner2(3);

            				// Some alias and calculations moved from heavy loops
            				float wCTF=1, wModulator=1.0;
            				float blobRadiusSquared = parent->blob.radius * parent->blob.radius;
            				float iDeltaSqrt = parent->iDeltaSqrt;
            				Matrix1D<double> & blobTableSqrt = parent->blobTableSqrt;
            				int xsize_1 = parent->VoutFourier.xdim - 1;
            				int zsize_1 = parent->VoutFourier.zdim - 1;
            				MultidimArray< std::complex<double> > &VoutFourier=parent->VoutFourier;
            				MultidimArray<double> &fourierWeights = parent->FourierWeights;

            				for (int i = 0; i < paddedFourier->ydim; i++) { // for each line
            					parent->rowsProcessed++;
            					if (statusArray[i] != 0) {
            						continue; // process only lines marked as '0', i.e. those should be processed
            					}

            					// otherwise process current line
            					for (int j = paddedFourier->xinit;
            							j <= paddedFourier->xinit + paddedFourier->xdim - 1;
            							j++) // for each point
            							{
            						// Compute the frequency of this coefficient in the
            						// universal coordinate system
            						FFT_IDX2DIGFREQ(j, parent->paddedImg.xdim, freq[0]);
            						FFT_IDX2DIGFREQ(i, parent->paddedImg.ydim, freq[1]);
            						freq[2] = 0;

            						if (freq[0] * freq[0] + freq[1] * freq[1] > parent->maxResolution2) {
            							continue;
            						}
            						wModulator = 1.0;
            						if (hasCTF && !reprocessFlag) {
            							contFreq[0] = freq[0] * iTs;
            							contFreq[1] = freq[1] * iTs;
            							threadParams->ctf.precomputeValues(contFreq[0],
            									contFreq[1]);
            							//wCTF=threadParams->ctf.getValueAt();
            							wCTF = threadParams->ctf.getValuePureNoKAt();
            							//wCTF=threadParams->ctf.getValuePureWithoutDampingAt();

            							if (std::isnan(wCTF)) {
            								if (i == 0 && j == 0)
            									wModulator = wCTF = 1.0;
            								else
            									wModulator = wCTF = 0.0;
            							}
            							if (fabs(wCTF) < parent->minCTF) {
            								wModulator = fabs(wCTF);
            								wCTF = SGN(wCTF);
            							} else
            								wCTF = 1.0 / wCTF;
            							if (parent->phaseFlipped)
            								wCTF = fabs(wCTF);
            						}

            	//					std::cout << "orig " << freq[0]
            	//							<< " "
            	//							<< freq[1]
            	//							<< " "
            	//							<< freq[2]
            	//							<< std::endl;


            						SPEED_UP_temps012;
            						M3x3_BY_V3x1(freq, *A_SL, freq);

            	//                    if (!((freq[0] > 0.0f) && (freq[1] > 0.0f) && (freq[2] > 0.0f)))
            	//                        continue;
            	#if DEBUG_DUMP >= 1
            				std::cout << "pred " << freq[0]
            						<< " "
            						<< freq[1]
            						<< " "
            						<< freq[2]
            	                    << " x "
            	                    << j
            	                    << " "
            	                    << i
            						<< std::endl;
            	#endif

            						// Look for the corresponding index in the volume Fourier transform
            						DIGFREQ2FFT_IDX_DOUBLE(freq[0], parent->volPadSizeX,
            								real_position[0]);
            						DIGFREQ2FFT_IDX_DOUBLE(freq[1], parent->volPadSizeY,
            								real_position[1]);
            						DIGFREQ2FFT_IDX_DOUBLE(freq[2], parent->volPadSizeZ,
            								real_position[2]);

            	//					{
            	//						double tmp[3];
            	//						FFT_IDX2DIGFREQ_DOUBLE(real_position[0], 720,
            	//								tmp[0]);
            	//						FFT_IDX2DIGFREQ_DOUBLE(real_position[1], 720,
            	//								tmp[1]);
            	//						FFT_IDX2DIGFREQ_DOUBLE(real_position[2], 720,
            	//								tmp[2]);
            	//					std::cout << "po " << tmp[0]
            	//							<< " "
            	//							<< tmp[1]
            	//							<< " "
            	//							<< tmp[2]
            	//							<< std::endl;
            	//					}

            						int roundedPosition[3];
            						for (int a = 0; a < 3; a++){
            							roundedPosition[a] = (int)(real_position[a] + 0.5f);
            						}

            						// Put a box around that coefficient
            						corner1[0] = roundedPosition[0];
            						corner1[1] = roundedPosition[1];
            						corner1[2] = roundedPosition[2];
            						corner2[0] = roundedPosition[0];
            						corner2[1] = roundedPosition[1];
            						corner2[2] = roundedPosition[2];

            	//						std::cout << freq[0]
            	//							<< " "
            	//							<< freq[1]
            	//							<< " "
            	//							<< freq[2]
            	//							<< std::endl;

            	#ifdef DEBUG

            						std::cout << "Idx Img=(0," << i << "," << j << ") -> Freq Img=("
            						<< freq.transpose() << ") ->\n    Idx Vol=("
            						<< real_position.transpose() << ")\n"
            						<< "   Corner1=" << corner1.transpose() << std::endl
            						<< "   Corner2=" << corner2.transpose() << std::endl;
            	#endif
            						// Loop within the box
            						double *ptrIn = (double *) &(*paddedFourier)(i, j);

            						// Some precalculations
            						for (int intz = corner1[2]; intz <= corner2[2]; ++intz) {
            							float z = 0;
            							z2precalculated(intz) = z * z;
            							if (zWrapped(intz) < 0) {
            								int iz, izneg;
            								fastIntWRAP(iz, intz, 0, zsize_1);
            								zWrapped(intz) = iz;
            								int miz = -iz;
            								fastIntWRAP(izneg, miz, 0, zsize_1);
            								zNegWrapped(intz) = izneg;
            							}
            						}
            						for (int inty = corner1[1]; inty <= corner2[1]; ++inty) {
            							float y = 0;
            							y2precalculated(inty) = y * y;
            							if (yWrapped(inty) < 0) {
            								int iy, iyneg;
            								fastIntWRAP(iy, inty, 0, zsize_1);
            								yWrapped(inty) = iy;
            								int miy = -iy;
            								fastIntWRAP(iyneg, miy, 0, zsize_1);
            								yNegWrapped(inty) = iyneg;
            							}
            						}
            						for (int intx = corner1[0]; intx <= corner2[0]; ++intx) {
            							float x = 0;
            							x2precalculated(intx) = x * x;
            							if (xWrapped(intx) < 0) {
            								int ix, ixneg;
            								fastIntWRAP(ix, intx, 0, zsize_1);
            								xWrapped(intx) = ix;
            								int mix = -ix;
            								fastIntWRAP(ixneg, mix, 0, zsize_1);
            								xNegWrapped(intx) = ixneg;
            							}
            						}

            						// Actually compute
            						processCube(i, j, corner1, corner2, z2precalculated, zWrapped,
            								zNegWrapped, y2precalculated, blobRadiusSquared,
            								yWrapped, yNegWrapped, x2precalculated, iDeltaSqrt,
            								wModulator, xWrapped, xsize_1, xNegWrapped,
            								reprocessFlag, wCTF, VoutFourier, blobTableSqrt,
            								threadParams, fourierWeights, ptrIn, threadParams->weight,
            								parent,
            								real_position);

            						statusArray[i] = -1;
            					}

            				}
            				break;
            	}
        default:
            break;
        }

        barrier_wait( barrier );
    }
    while ( 1 );
//    std::cout << std::endl;
}


void testMapping(MultidimArray<double> paddedImg,
		Matrix2D<double> * A_SL,
		int volPadSizeX,
		int volPadSizeY,
		int volPadSizeZ) {

	int fx = 23;
	int fy = 0;

	Matrix1D<double> ffreq(3), freal_pos(3), rreal_pos(3), rfreq(3);
	FFT_IDX2DIGFREQ(fx, paddedImg.xdim, ffreq[0]);
	FFT_IDX2DIGFREQ(fy, paddedImg.ydim, ffreq[1]);
	ffreq[2] = 0;
	std::cout << ffreq << std::endl;
	SPEED_UP_temps012;
	M3x3_BY_V3x1(ffreq, *A_SL, ffreq);
	std::cout << ffreq << std::endl;
	// Look for the corresponding index in the volume Fourier transform
	DIGFREQ2FFT_IDX_DOUBLE(ffreq[0], volPadSizeX,
			freal_pos[0]);
	DIGFREQ2FFT_IDX_DOUBLE(ffreq[1], volPadSizeY,
			freal_pos[1]);
	DIGFREQ2FFT_IDX_DOUBLE(ffreq[2], volPadSizeZ,
			freal_pos[2]);
	std::cout << freal_pos << std::endl;

	FFT_IDX2DIGFREQ_DOUBLE(freal_pos[0], volPadSizeX, rreal_pos[0]);
	FFT_IDX2DIGFREQ_DOUBLE(freal_pos[1], volPadSizeY, rreal_pos[1]);
	FFT_IDX2DIGFREQ_DOUBLE(freal_pos[2], volPadSizeZ, rreal_pos[2]);
	std::cout << rreal_pos << std::endl;
	// invert map
	Matrix2D<double> A_SL_inv = A_SL->inv();
	std::cout << A_SL_inv << std::endl;
	M3x3_BY_V3x1(rreal_pos, A_SL_inv, rreal_pos);
	std::cout << rreal_pos << std::endl;
//	rreal_pos[2] = 0;
	DIGFREQ2FFT_IDX_DOUBLE(rreal_pos[0], paddedImg.xdim, rfreq[0]);
	DIGFREQ2FFT_IDX_DOUBLE(rreal_pos[1], paddedImg.ydim, rfreq[1]);
	DIGFREQ2FFT_IDX_DOUBLE(rreal_pos[2], paddedImg.zdim, rfreq[2]);
	std::cout << rfreq << std::endl;


}

template <typename T>
inline std::complex<T>
BilinearInterpolation(std::complex<T> bl, std::complex<T> br, std::complex<T> tl, std::complex<T> tr,
		float x, float y)
{
	std::complex<T> top = (T)x * (tl - tr) + tr;
	std::complex<T> bottom = (T)x * (bl - br) + br;
	return (T)y * (top - bottom) + bottom;
}

template<typename T, typename U>
inline U clamp(U val, T min, T max) {
	U res = val;
	res = (res > max) ? max : res;
	res = (res < min) ? min : res;
	return res;
}

inline std::complex<float>
getPixelBilinear(std::complex<float>** img, float x, float y, int imgSizeX, int imgSizeY)
{
	float posX, posY;
	float fractX, fractY;
	fractX = std::modf(x , &posX);
	fractY = std::modf(y , &posY);
	std::complex<float> bl,br,tl,tr;
	bl = br = tl = tr = 0; // values outside of the boundaries will be zero
	if (posY >= 0) {
		if (posX >= 0) {
			bl = img[(int)posY][(int)posX];
		}
		if ((posX + 1)< imgSizeX) {
			br = img[(int)posY][(int)posX + 1];
		}
	}
	if ((posY+1) < imgSizeY) {
		if (posX >= 0) {
			tl = img[(int)posY+1][(int)posX];
		}
		if ((posX + 1) < imgSizeX) {
			tr = img[(int)posY+1][(int)posX + 1];
		}
	}
	return BilinearInterpolation(bl, br, tl, tr, 1. - fractX, 1. - fractY);
}

inline std::complex<float>
getPixelClamped(std::complex<float>** img, float x, float y, int imgSizeX, int imgSizeY)
{
	int imgX = clamp((int)(x + 0.5f), 0, imgSizeX - 1);
	int imgY = clamp((int)(y + 0.5f), 0, imgSizeY - 1);
	return img[imgY][imgX];
}


#if 0
template <typename T>
bool pad(size_t maxX, int offsetX, int offsetY, T& x, T& y){
    /*int ox = x;
    int oy = y;
    if (y < (float)(offsetY)/2.0f)
        y += (float)(offsetY)/2.0f;
    else
        y -= (float)(offsetY)/2.0f;

    if (x > (float)(offsetX-1)/2.0f)
        x -= (float)(offsetX-1);
    bool ret = false;
    if (x < 0.0f) {
        x = -x;
        y = (float)(offsetY-1) - y;
        //if (x > (float)(offsetX-1) || y > (float)(offsetY-1) || x < 0.0f || y < 0.0f) std::cout << "HAF HAF\n";
        ret = true;
    }
    //if (x > (float)(offsetX-1) || y > (float)(offsetY-1) || x < 0.0f || y < 0.0f) std::cout << "HAF HAF\n";
    if (y < (float)(offsetY)/2.0f)
        y += (float)(offsetY)/2.0f;
    else
        y -= (float)(offsetY)/2.0f;

    if (x > (float)(offsetX-1) || y > (float)(offsetY-1) || x < 0.0f || y < 0.0f) 
        std::cout << "HAF HAF << " << ox << " " << oy << " -> " << x << " " << y << "\n";

    return ret;*/

	if (x > (float)(offsetX/2-1)) {
        std::cout << "A\n";
		x = (float)offsetX - x;
        y = (float)offsetY - y;
	}

    if ((x < 0.0f) && (y < 0.0f)) {
        std::cout << "B\n";
        x = -x;
        y = -y;
    }

	if (x < 0.f) {
        std::cout << "C\n";
		//x += (float)offsetX;
        x = -x;
        y = offsetY-1 - y;
	}
	if (y > (float)(offsetY-1)) {
        std::cout << "D\n";
		y -= (float)offsetY;
	}
	if (y < 0.f) {
        std::cout << "A\n";
		y += (float)offsetY;
	}
	if ((int)x >= maxX) { // Careful here, we must not get under -1. -0.9999 is fine, as cast to int will convert it to zero(0)
        std::cout << "(int)x >= maxX\n";
		x = - x + offsetX/2;
		y = - y + offsetY;
        return true;
	}
	if (y > maxX) {// FIXME this is probably wrong
        std::cout << "(y > maxX)\n";
		return true;
    }
	return false;
}
#endif

static clock_t ticks_pad = 0;

template <typename T>
bool pad(size_t maxX, int offsetX, int offsetY, T& x, T& y, size_t imgSizeX, size_t imgSizeY){

    clock_t start = clock();

    bool conjugate = false;
    float halfX = (float)offsetX/2.0f;
    float halfY = (float)offsetY/2.0f;

    if (x > halfX-1) {
#if DEBUG_DUMP >= 2
        std::cout << "x > halfX-1\n";
#endif
        x = (float)offsetX - x;
        y = (float)offsetY - y;
        conjugate = true;
    }

    if ((x < 0.0f) && (y < 0.0f)) {
#if DEBUG_DUMP >= 2
        std::cout << "x < 0 && y < 0\n";
#endif
        x = -x;
        y = -y;
        conjugate = true;
    }
    if (x < 0.f && y > offsetY-1) {
#if DEBUG_DUMP >= 2
        std::cout << "x < 0 && y > offsetY-1\n";
#endif
        x = -x;
        y = 2*(offsetY-1) - y;
        conjugate = true;
    }
    if (x < 0.f && y >= 0.f) {
#if DEBUG_DUMP >= 2
        std::cout << "x < 0 && y > 0\n";
#endif
        //x += (float)offsetX;
        x = -x;
        y = offsetY-1 - y;
        conjugate = true;
    }
    if (x >= 0.f && y < 0.f) {
#if DEBUG_DUMP >= 2
        std::cout << "x > 0 && y < 0\n";
#endif
        y += (float)offsetY;
    }
    if (y > (float)(offsetY-1)) {
#if DEBUG_DUMP >= 2
        std::cout << "y > offsetY-1\n";
#endif
        y -= (float)offsetY;
    }

    ticks_pad += clock() - start;

    return conjugate;
}

inline void multiply(Matrix2D<double>& transform, float inOut[3]) {
	float tmp0 = transform(0, 0) * inOut[0] + transform(0, 1) * inOut[1] + transform(0, 2) * inOut[2];
	float tmp1 = transform(1, 0) * inOut[0] + transform(1, 1) * inOut[1] + transform(1, 2) * inOut[2];
	float tmp2 = transform(2, 0) * inOut[0] + transform(2, 1) * inOut[1] + transform(2, 2) * inOut[2];
	inOut[0] = tmp0;
	inOut[1] = tmp1;
	inOut[2] = tmp2;
}

inline void multiply(Matrix2D<double>& transform, Point3D& inOut) {
	float tmp0 = transform(0, 0) * inOut.x + transform(0, 1) * inOut.y + transform(0, 2) * inOut.z;
	float tmp1 = transform(1, 0) * inOut.x + transform(1, 1) * inOut.y + transform(1, 2) * inOut.z;
	float tmp2 = transform(2, 0) * inOut.x + transform(2, 1) * inOut.y + transform(2, 2) * inOut.z;
	inOut.x = tmp0;
	inOut.y = tmp1;
	inOut.z = tmp2;
}

inline void multiply(float transform[3][3], Point3D& inOut) {
	float tmp0 = transform[0][0] * inOut.x + transform[0][1] * inOut.y + transform[0][2] * inOut.z;
	float tmp1 = transform[1][0] * inOut.x + transform[1][1] * inOut.y + transform[1][2] * inOut.z;
	float tmp2 = transform[2][0] * inOut.x + transform[2][1] * inOut.y + transform[2][2] * inOut.z;
	inOut.x = tmp0;
	inOut.y = tmp1;
	inOut.z = tmp2;
}

inline void multiply(float transform[3][3], float inOut[3]) {
	float tmp0 = transform[0][0] * inOut[0] + transform[0][1] * inOut[1] + transform[0][2] * inOut[2];
	float tmp1 = transform[1][0] * inOut[0] + transform[1][1] * inOut[1] + transform[1][2] * inOut[2];
	float tmp2 = transform[2][0] * inOut[0] + transform[2][1] * inOut[1] + transform[2][2] * inOut[2];
	inOut[0] = tmp0;
	inOut[1] = tmp1;
	inOut[2] = tmp2;
}


template <typename T>
inline std::complex<T>
getPixelValue( MultidimArray< std::complex<T> > *img, int offsetX, int offsetY, float x, float y)
{
	T posX, posY;
	T fractX;
	T fractY;
	fractX = std::modf(x , &posX);
	fractY = std::modf(y , &posY);

	pad(img->xdim, offsetX, offsetY, posX, posY);
	std::complex<T> bottom_left = DIRECT_A2D_ELEM(*img, (int)posY, (int)posX);
	pad(img->xdim, offsetX, offsetY, ++posX, posY);
	std::complex<T> bottom_right = DIRECT_A2D_ELEM(*img, (int)posY, (int)posX);
	pad(img->xdim, offsetX, offsetY, --posX, ++posY);
	std::complex<T> top_left = DIRECT_A2D_ELEM(*img, (int)posY, (int)posX);
	pad(img->xdim, offsetX, offsetY, ++posX, posY);
	std::complex<T> top_right = DIRECT_A2D_ELEM(*img, (int)posY, (int)posX);

	return BilinearInterpolation(bottom_left, bottom_right, top_left, top_right, 1. - fractX, 1. - fractY);
}

template <typename T>
std::complex<T>
getVoxelValue( MultidimArray< std::complex<T> > *img, int offsetX, int offsetY,
		float x, float y, float z,
		blobtype blob, T& weight,
		int xdim, int ydim, double maxResolutionSqr,
		int k, int l, int m,
		Matrix1D<double>& blobTableSqrt, float iDeltaSqrt)
{
	return getVoxelValue(img, offsetX, offsetY, x, y, z, blob.radius, weight, xdim, ydim, maxResolutionSqr,
			k, l, m, blobTableSqrt, iDeltaSqrt);
}

std::complex<float>
getVoxelValue(std::complex<float>** img, int offsetX, int offsetY,
		float x, float y, float z,
		float blobRadius, float& weight,
		int xdim, int ydim, double maxResolutionSqr,
		float k, float l, float m,
		Matrix1D<double>& blobTableSqrt, float iDeltaSqrt,
		size_t imgSizeX, size_t imgSizeY)
{
    bool dumpujeme = false;
    //std::cout << "YYY " << k << " " << l << " " << m << std::endl;
    /*if ((fabsf(k-0.0f) < 0.00005) && (fabsf(l-0.0f) < 0.00005) && (fabsf(m-0.0f) < 0.00005)) {
        std::cout << "XXX dumpujeme\n";
        dumpujeme = true;
    }*/
	weight = 0.;
	float radiusSqr = blobRadius * blobRadius;
	std::complex<float> result = (0, 0);
	int minX = CEIL(x - blobRadius);
	int maxX = FLOOR(x + blobRadius);
	int minY = CEIL(y - blobRadius);
	int maxY = FLOOR(y + blobRadius);
#if DEBUG_DUMP >= 2
    std::cout << "center: " << x << " " << y << std::endl;
#endif
	float zSqr = z * z;
	for (int i = minY; i <= maxY; i++) {
		float ySqr = (y - i) * (y - i);
		for (int j = minX; j <= maxX; j++) {
			float xD = x - j;
			float distanceSqr = ySqr + xD*xD + zSqr;
			if (distanceSqr > radiusSqr) {
				continue;
			}
//			T freq[2];
//			FFT_IDX2DIGFREQ(j, xdim, freq[0]);
//			FFT_IDX2DIGFREQ(i, ydim, freq[1]);
//			if (freq[0] * freq[0] + freq[1] * freq[1] > maxResolutionSqr) {
//				continue;
//			}
			// all tests passed, you can use the pixel value
			int tmpX = j;
			int tmpY = i;
			int aux = (int) ((distanceSqr * iDeltaSqrt + 0.5)); //Same as ROUND but avoid comparison
			float tmpWeight = blobTableSqrt[aux];
    //if (x > (float)(offsetX-1) || y > (float)(offsetY-1) || x < 0.0f || y < 0.0f) std::cout << "HAF HAF\n";
//			T tmpWeight = blob_val(std::sqrt(distanceSqr), blob);
			bool conjugate = pad(imgSizeX, offsetX, offsetY, tmpX, tmpY, imgSizeX, imgSizeY); // FIXME use X dimension of the image
			float conj = conjugate ? -1. : 1.;
			std::complex<float> pixVal = img[tmpY][tmpX];
			result.real() += tmpWeight * pixVal.real();
			result.imag() += tmpWeight * conj * pixVal.imag();
			weight += tmpWeight;
            /*if (dumpujeme) {
                std::cout << "XXX " << j << " " << i << ", " << tmpX << " " << tmpY << ", " << pixVal << " " << tmpWeight << std::endl;
            }*/


//			float dest[3];
//			if (k < 0) {
//				dest[0] = -k;
//				dest[1] = -l;
//				dest[2] = -m;
//			}

			float tmp[3];
//			          				if (k < 0) {
//			          					tmp[0] = -k;
//			          					tmp[1] = -l;
//			          					tmp[2] = -m;
//			          				} else {
//			          					tmp[0] = k;
//			          					tmp[1] = l;
//			          					tmp[2] = m;
//			          				}
//			          				DIGFREQ2FFT_IDX(tmp[0], 720, tmp[0]);
//			          				DIGFREQ2FFT_IDX(tmp[1], 720, tmp[1]);
//			          				DIGFREQ2FFT_IDX(tmp[2], 720, tmp[2]);


#if DEBUG_DUMP >= 1
			std::cout << tmp[0]
					<< " "
					<< tmp[1]
					<< " "
					<< tmp[2]
					<< " -> "
					<< k
					<< " "
					<< l
					<< " "
					<< m
					<< " -> "
					<< tmpX
					<< " "
					<< tmpY
					<< " "
					<< tmpWeight
					<< " "
					<< (conjugate ? "conjugate" : "no_conjugate")
					<< std::endl;
#endif
		}
	}
	return result;

}



template <typename T>
inline bool inRange(T x, T min, T max) {
	return (x > min) && (x < max);
}



// 7____6
//3/__2//
// ||  ||
// ||4 ||5
//0|/__|/1
Point3D* createBoundingCuboid(size_t paddedVolumeSize,
		size_t conserveRows,
		float blobRadius)
{
	static Point3D result[8];
	float sizeLength = (float)conserveRows / (float)paddedVolumeSize;
	float height = ceil(blobRadius) / (float)paddedVolumeSize;
	result[0].x = result[3].x = result[4].x = result[7].x = 0.f - height;
	result[1].x = result[2].x = result[5].x = result[6].x = sizeLength + height;

	result[0].y = result[1].y = result[4].y = result[5].y = -sizeLength - height;
	result[2].y = result[3].y = result[6].y = result[7].y = sizeLength + height;

	result[0].z = result[1].z = result[2].z = result[3].z = height;
	result[4].z = result[5].z = result[6].z = result[7].z = -height;

	return result;
}

// 7____6
//3/__2//
// ||  ||
// ||4 ||5
//0|/__|/1
// origin is in the center [0, 0, 0]
Point3D* createBoundingCuboid(float sizeX, float sizeY,	float sizeZ)
{
	static Point3D result[8];
	float halfX = sizeX/2;
	float halfY = sizeY/2;
	float halfZ = sizeZ/2;
	result[0].x = result[3].x = result[4].x = result[7].x = -halfX;
	result[1].x = result[2].x = result[5].x = result[6].x = halfX;

	result[0].y = result[1].y = result[4].y = result[5].y = -halfY;
	result[2].y = result[3].y = result[6].y = result[7].y = halfY;

	result[0].z = result[1].z = result[2].z = result[3].z = -halfZ;
	result[4].z = result[5].z = result[6].z = result[7].z = halfZ;

	return result;
}



//       3___2
//   +   |   |
// [0,0] |   |y
//   -  0|___|1  sizes are padded with blob-radius
//         x
Point3D* createProjectionPlane(size_t paddedVolumeSize,
		size_t conserveRows,
		float blobRadius)
{
	static Point3D result[4];
	float sizeLength = (float)conserveRows / (float)paddedVolumeSize;
	float padding = ceil(blobRadius) / (float)paddedVolumeSize;
	result[0].x = result[3].x = - (0.f + padding);
	result[1].x = result[2].x = sizeLength + padding;

	result[0].y = result[1].y = -(sizeLength + padding);
	result[2].y = result[3].y = sizeLength + padding;

	result[0].z = result[1].z = result[2].z = result[3].z = 0.f;

	return result;
}

//       3___2
//   +   |   |
// [0,0] |   |y
//   -  0|___|1  sizes are padded with blob-radius
//         x
// [0,0] is in the middle of the left side (point [0] and [3])
Point3D* createProjectionPlane(float sizeX,	float sizeY)
{
	static Point3D result[4];
	float halfY = sizeY / 2.0f;
	result[0].x = result[3].x = 0.f;
	result[1].x = result[2].x = sizeX;

	result[0].y = result[1].y = -halfY;
	result[2].y = result[3].y = halfY;

	result[0].z = result[1].z = result[2].z = result[3].z = 0.f;

	return result;
}

inline bool getZ(float x, float y, float& z, Point3D* plane) {
	// from parametric eq. of the plane
	float x0 = plane[0].x;
	float y0 = plane[0].y;
	float z0 = plane[0].z;
	Point3D a = { plane[1].x - x0, plane[1].y - y0, plane[1].z - z0 };
	Point3D b = { plane[3].x - x0, plane[3].y - y0, plane[3].z - z0 };

	float u = ((y-y0)*a.x + (x0-x)*a.y) / (a.x * b.y - b.x * a.y);
	float t = (-x0 + x - u*b.x) / (a.x);

	z = z0 + t*a.z + u*b.z;
	return inRange(t, 0.f, 1.f) && inRange(u, 0.f, 1.f);
}

inline bool getZ(float x, float y, float& z, const Point3D& a, const Point3D& b, const Point3D& p0) {
	// from parametric eq. of the plane
	float x0 = p0.x;
	float y0 = p0.y;
	float z0 = p0.z;

	float u = ((y-y0)*a.x + (x0-x)*a.y) / (a.x * b.y - b.x * a.y);
	float t = (-x0 + x - u*b.x) / (a.x);

	z = z0 + t*a.z + u*b.z;
	return inRange(t, 0.f, 1.f) && inRange(u, 0.f, 1.f);
}

inline bool getY(float x, float& y, float z, const Point3D& a, const Point3D& b, const Point3D& p0) {
	// from parametric eq. of the plane
	float x0 = p0.x;
	float y0 = p0.y;
	float z0 = p0.z;

	float u = ((z-z0)*a.x + (x0-x)*a.z) / (a.x * b.z - b.x * a.z);
	float t = (-x0 + x - u*b.x) / (a.x);

	y = y0 + t*a.y + u*b.y;
	return inRange(t, 0.f, 1.f) && inRange(u, 0.f, 1.f);
}

inline bool getX(float& x, float y, float z, const Point3D& a, const Point3D& b, const Point3D& p0) {
	// from parametric eq. of the plane
	float x0 = p0.x;
	float y0 = p0.y;
	float z0 = p0.z;

	float u = ((z-z0)*a.y + (y0-y)*a.z) / (a.y * b.z - b.y * a.z);
	float t = (-y0 + y - u*b.y) / (a.y);

	x = x0 + t*a.x + u*b.x;
	return inRange(t, 0.f, 1.f) && inRange(u, 0.f, 1.f);
}

Point3D* pad(Point3D* cuboid) {
	for (int i = 0; i < 8; i++) {
		if (cuboid[i].x < 0) {
			cuboid[i].x *= -1;
			cuboid[i].y *= -1;
		}
	}
	return cuboid;
}

void print(Point3D* cuboid) {
	for (int i = 0; i < 9; i++) {
		std::cout << cuboid[i%8].x << " " << cuboid[i%8].y << " " << cuboid[i%8].z << std::endl;
	}
}

void printPlane(Point3D* plane) {
	for (int i = 0; i < 5; i++) {
		std::cout << plane[i%4].x << " " << plane[i%4].y << " " << plane[i%4].z << std::endl;
	}
}

void print1(Point3D* cuboid) {
	std::cout << cuboid[0].x << " " << cuboid[0].y << " " << cuboid[0].z << std::endl;
	std::cout << cuboid[1].x << " " << cuboid[1].y << " " << cuboid[1].z << std::endl;
	std::cout << cuboid[2].x << " " << cuboid[2].y << " " << cuboid[2].z << std::endl;
	std::cout << cuboid[3].x << " " << cuboid[3].y << " " << cuboid[3].z << std::endl;
	std::cout << cuboid[0].x << " " << cuboid[0].y << " " << cuboid[0].z << std::endl;
	std::cout << cuboid[2].x << " " << cuboid[2].y << " " << cuboid[2].z << std::endl;
}

Point3D* transformCuboid(Point3D* cuboid, Matrix2D<double>& transform) {
	float tmp[3];
	for (int i = 0; i < 8; i++) {
		tmp[0] = cuboid[i].x;
		tmp[1] = cuboid[i].y;
		tmp[2] = cuboid[i].z;
		cuboid[i].x = transform(0,0)*tmp[0] + transform(0,1)*tmp[1] + transform(0,2)*tmp[2];
		cuboid[i].y = transform(1,0)*tmp[0] + transform(1,1)*tmp[1] + transform(1,2)*tmp[2];
		cuboid[i].z = transform(2,0)*tmp[0] + transform(2,1)*tmp[1] + transform(2,2)*tmp[2];
	}
	return cuboid;
}

Point3D* rotateCuboid(Point3D* cuboid, Matrix2D<double>& transform) {
	float tmp[3];
	for (int i = 0; i < 8; i++) {
		tmp[0] = cuboid[i].x;
		tmp[1] = cuboid[i].y;
		tmp[2] = cuboid[i].z;
		cuboid[i].x = transform(0,0)*tmp[0] + transform(0,1)*tmp[1] + transform(0,2)*tmp[2];
		cuboid[i].y = transform(1,0)*tmp[0] + transform(1,1)*tmp[1] + transform(1,2)*tmp[2];
		cuboid[i].z = transform(2,0)*tmp[0] + transform(2,1)*tmp[1] + transform(2,2)*tmp[2];
	}
	return cuboid;
}

//Point3D* rotateCuboid(Point3D* cuboid, Matrix2D<double>& transform, int tmp) {
//	Point3D center = {};
//	for (int i = 0; i < 8; i++) {
//		center.x += cuboid[i].x;
//		center.y += cuboid[i].y;
//		center.z += cuboid[i].z;
//	}
//	center.x /= 8.f;
//	center.y /= 8.f;
//	center.z /= 8.f;
//	std::cout << "center " << center.x << " " << center.y << " " <<center.z << std::endl;
// 	float tmp[3];
//	for (int i = 0; i < 8; i++) {
//		tmp[0] = cuboid[i].x - center.x;
//		tmp[1] = cuboid[i].y - center.y;
//		tmp[2] = cuboid[i].z - center.z;
//		multiply(transform, tmp);
//		cuboid[i].x = tmp[0] + center.x;
//		cuboid[i].y = tmp[1] + center.y;
//		cuboid[i].z = tmp[2] + center.z;
//	}
//	return cuboid;
//}

Point3D* translateCuboid(Point3D* cuboid, Point3D vector) {
	for (int i = 0; i < 8; i++) {
		cuboid[i].x += vector.x;
		cuboid[i].y += vector.y;
		cuboid[i].z += vector.z;
	}
	return cuboid;
}

Point3D* transformPlane(Point3D* plane, Matrix2D<double>& transform) {
	float tmp[3];
	for (int i = 0; i < 4; i++) {
		tmp[0] = plane[i].x;
		tmp[1] = plane[i].y;
		tmp[2] = plane[i].z;
		plane[i].x = transform(0,0)*tmp[0] + transform(0,1)*tmp[1] + transform(0,2)*tmp[2];
		plane[i].y = transform(1,0)*tmp[0] + transform(1,1)*tmp[1] + transform(1,2)*tmp[2];
		plane[i].z = transform(2,0)*tmp[0] + transform(2,1)*tmp[1] + transform(2,2)*tmp[2];
	}
	return plane;
}

Point3D* rotatePlane(Point3D* plane, float transform[3][3]) {
	for (int i = 0; i < 4; i++) {
		multiply(transform, plane[i]);
	}
	return plane;
}

Point3D* translatePlane(Point3D* plane, Point3D vector) {
	for (int i = 0; i < 4; i++) {
		plane[i].x += vector.x;
		plane[i].y += vector.y;
		plane[i].z += vector.z;
	}
	return plane;
}

Point3D* getAABB(Point3D* cuboid) {
	static Point3D AABB[2];
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
	return AABB;
}

Point3D* getPlaneAABB(Point3D* plane) {
	static Point3D AABB[2];
	AABB[0].x = AABB[0].y = AABB[0].z = std::numeric_limits<float>::max();
	AABB[1].x = AABB[1].y = AABB[1].z = std::numeric_limits<float>::min();
	Point3D tmp;
	for (int i = 0; i < 4; i++) {
		tmp = plane[i];
		if (AABB[0].x > tmp.x) AABB[0].x = tmp.x;
		if (AABB[0].y > tmp.y) AABB[0].y = tmp.y;
		if (AABB[0].z > tmp.z) AABB[0].z = tmp.z;
		if (AABB[1].x < tmp.x) AABB[1].x = tmp.x;
		if (AABB[1].y < tmp.y) AABB[1].y = tmp.y;
		if (AABB[1].z < tmp.z) AABB[1].z = tmp.z;
	}
	return AABB;
}

Point3D* getAABB(Point3D* cuboid,
	float minX, float minY, float minZ,
	float maxX, float maxY, float maxZ) {
	Point3D* AABB = getAABB(cuboid);
	// limit to max size
	if (AABB[0].x < minX) AABB[0].x = minX;
	if (AABB[0].y < minY) AABB[0].y = minY;
	if (AABB[0].z < minZ) AABB[0].z = minZ;
	if (AABB[1].x > maxX) AABB[1].x = maxX;
	if (AABB[1].y > maxY) AABB[1].y = maxY;
	if (AABB[1].z > maxZ) AABB[1].z = maxZ;
	
	return AABB;
}

Point3D* getPlaneAABB(Point3D* plane,
	float minX, float minY, float minZ,
	float maxX, float maxY, float maxZ) {
	Point3D* AABB = getPlaneAABB(plane);
	// limit to max size
	if (AABB[0].x < minX) AABB[0].x = minX;
	if (AABB[0].y < minY) AABB[0].y = minY;
	if (AABB[0].z < minZ) AABB[0].z = minZ;
	if (AABB[1].x > maxX) AABB[1].x = maxX;
	if (AABB[1].y > maxY) AABB[1].y = maxY;
	if (AABB[1].z > maxZ) AABB[1].z = maxZ;

	return AABB;
}

void druhyPokus(std::complex<float>*** outputVolume,
		float*** outputWeight,
		size_t conserveRows,
		std::complex<float>** paddedFourier,
		size_t paddedVolumeSize,
		float maxResolutionSqr,
		float blobRadius,
		Matrix2D<double>& transform,
		Matrix1D<double>& blobTableSqrt, float iDeltaSqrt,
		size_t imgSizeX, size_t imgSizeY) {
//std::cout << "plane:" << std::endl;
//	Point3D* plane = transformPlane(createProjectionPlane(paddedVolumeSize, conserveRows, blobRadius), transform);
//	printPlane(plane);
//std::cout << "cuboid:" << std::endl;
	Point3D* cuboid = transformCuboid(createBoundingCuboid(paddedVolumeSize, conserveRows, blobRadius), transform);
//	print(cuboid);
//std::cout << "cuboid:" << std::endl;
//	cuboid = createBoundingCuboid(paddedVolumeSize, conserveRows, blobRadius);
//	print(cuboid);
//std::cout << "test:" << std::endl;

	float upperIndex = (float)conserveRows / (float)paddedVolumeSize;
	float step = upperIndex / (float)(conserveRows);
	int size = 2* conserveRows;

	Matrix2D<double> transfInv = transform.inv();
	float freq[3];
	float img_pos[3];

	// test orientation
	Point3D u;
	Point3D v;
	getVectors(cuboid, u, v);
	Point3D normal = getNormal(u, v);
	float XY = std::abs(getCosTheta(normal, Point3D {0.f, 0.f, 1.f})); // iterate XY plane
	float XZ = std::abs(getCosTheta(normal, Point3D {0.f, 1.f, 0.f})); // iterate XZ plane
	float YZ = std::abs(getCosTheta(normal, Point3D {1.f, 0.f, 0.f})); // iterate YZ plane
	if ((XY >= XZ) && (XY >= YZ)) {
//		std::cout << "XY" << std::endl;
	}
	if ((XZ >= XY) && (XZ >= YZ)) {
//		std::cout << "XZ" << std::endl;
	}
	if ((YZ >= XY) && (YZ >= XZ)) {
//		std::cout << "YZ" << std::endl; // FIXME implement proper space traversing
	}

	Point3D* AABB = getAABB(cuboid);
	int minY, minX;
	int maxY, maxX;
	minY = floor((AABB[0].y + upperIndex) * paddedVolumeSize);
	minX = floor((AABB[0].x + upperIndex) * paddedVolumeSize);
	maxY = ceil((AABB[1].y + upperIndex) * paddedVolumeSize);
	maxX = ceil((AABB[1].x + upperIndex) * paddedVolumeSize);
    minX = std::max(0, minX);
    minY = std::max(0, minY);
    maxX = std::min((int)size-1, maxX);
    maxY = std::min((int)size-1, maxY);
    

//	for (int m = 0; m < size; m++ ) {
		for (int l = minY; l < maxY; l++) {
			for (int k = minX; k < maxX; k++) {
				float x = k*step - upperIndex;
				float y = l*step - upperIndex;
				float z1;
				float z2;
				bool hit1 = getZ(x, y, z1, u, v, *cuboid); // FIXME once iterating in proper plane, can be reduced to one calculation +- blob distance in 'Z' axis
				bool hit2 = getZ(x, y, z2, u, v, *(cuboid + 4));
				if (hit1 || hit2) {
					float mLower = (std::min(z1, z2) + upperIndex) * paddedVolumeSize ;
                    mLower = std::max(0.0f, mLower);
					float mUpper = (std::max(z1, z2) + upperIndex) * paddedVolumeSize ;
                    mUpper = std::min((float)(size-1), mUpper); //FIXME hardcoded value
//					float tmp = std::pow(ceil(blobRadius), 1/3.)( *  / (float)paddedVolumeSize;
					for (int m = floor(mLower); m <= ceil(mUpper); m++) {
						float z = m*step - upperIndex;
						freq[0] = x;
						freq[1] = y;
						freq[2] = z;

                        //if (!((x > 0) && (y > 0) && (z >0)))
                        //    continue;

					// invert map
						multiply(transfInv, freq);

						DIGFREQ2FFT_IDX_DOUBLE(freq[0], paddedVolumeSize, img_pos[0]); // pixel index
						DIGFREQ2FFT_IDX_DOUBLE(freq[1], paddedVolumeSize, img_pos[1]); // line index
						img_pos[2] = std::abs(freq[2] * paddedVolumeSize); // XXX HACK store transformed Z value

						// check restrictions
						if (inRange((size_t) img_pos[1], conserveRows, paddedVolumeSize - conserveRows) // or the line has too high frequencies
								|| (img_pos[2] > blobRadius) // or projection plane is too far
								|| ((freq[0] * freq[0] + freq[1] * freq[1]) > maxResolutionSqr) // or the pixel frequency is too high
								) {
							continue;
						}
//													std::cout << x
//															<< " "
//															<< y
//															<< " "
//															<< z
//															<< std::endl;

/*														std::cout << x
																<< " "
    															<< y
																<< " "
																<< z
																<< " -> "
																<< img_pos[0]
																<< " "
																<< img_pos[1]
																<< " "
																<< img_pos[2]
					//											<< " "
					//											<< (conjugate ? "conjugate" : "no_conjugate")
																<< std::endl;*/

						float weight;



						std::complex<float> val = getVoxelValue(paddedFourier,
								paddedVolumeSize, paddedVolumeSize, img_pos[0], img_pos[1],
								img_pos[2], blobRadius, weight, paddedVolumeSize,
								paddedVolumeSize, maxResolutionSqr, x, y, z,
								blobTableSqrt, iDeltaSqrt,
								imgSizeX, imgSizeY);
                        //if ((m == l == 64) && (k == 63)) std::cout << "YYY " << val << " " << weight << " " << outputVolume[m][l][k] << std::endl;
                        //if (val.imag() < 100.0) std::cout << "XXX " << m << " " << l << " " << k << std::endl;

						outputWeight[m][l][k] += weight;
						outputVolume[m][l][k] += val;


//
//					DIRECT_A3D_ELEM(VoutFourier_muj, m, l, k) += val;
//					//									DIRECT_A3D_ELEM(VoutFourier_muj, m, l, k).real() += weight * val.real();
//					//									DIRECT_A3D_ELEM(VoutFourier_muj, m, l, k).imag() += weight * conj * val.imag();
//					DIRECT_A3D_ELEM(FourierWrights_moje, m, l, k) += weight;



					}
				}
			}
		}
//	}

//	print(transformCuboid(createBoundingCuboid(paddedVolumeSize, conserveRows, blobRadius), transform));
//	print(createBoundingCuboid(paddedVolumeSize, conserveRows, blobRadius));
//	std::cout << "BB" << std::endl;
//
//	Matrix1D<float> freq(3);
//
//	for (size_t y = 0; y < paddedFourier->ydim; y++) { // for each line
//		if (inRdruhyPokusange(y, conserveRows, paddedFourier->ydim - conserveRows)) {
//			continue; // process only lines with small frequencies
//		}
//		// otherwise process current line
//		for (int x = 0;	x < conserveRows; x++) { // for each pixel (with low frequency)
//			// Compute the frequency of this coefficient in the
//			// universal coordinate system
//			FFT_IDX2DIGFREQ(x, paddedVolumeSize, freq[0]);
//			FFT_IDX2DIGFREQ(y, paddedVolumeSize, freq[1]);
//			freq[2] = 0;
//
//			if (freq[0] * freq[0] + freq[1] * freq[1] > maxResolutionSqr) {
//				continue;
//			}
////			SPEED_UP_temps012;
////			M3x3_BY_V3x1(freq, transform, freq);
//			std::cout << freq[0] << " " << freq[1] << " " << freq[2] << " -> " << x << " " << y << std::endl;
//
//		}
//	}

}

void printAABB(Point3D* AABB) {
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

void getVoxelValue(std::complex<float>& outputVal, float& outputWeight,
	float imgPos[3],
	std::complex<float>** img, int imgSizeX, int imgSizeY,
	Matrix1D<double>& blobTableSqrt, float iDeltaSqrt, float blobRadius) {
	float maxDistanceSqr = blobRadius*blobRadius + imgSizeX*imgSizeX;
	outputWeight = 0.f;
	outputVal = 0.f;
	int minX = ceil(imgPos[0] - blobRadius);
	int maxX = floor(imgPos[0] + blobRadius);
	int minY = ceil(imgPos[1] - blobRadius);
	int maxY = floor(imgPos[1] + blobRadius);
	float zSqr = imgPos[2] * imgPos[2];
	for (int i = minY; i <= maxY; i++) {
		float ySqr = (imgPos[1] - i) * (imgPos[1] - i);
		for (int j = minX; j <= maxX; j++) {
			float xD = imgPos[0] - j;
			float distanceSqr = ySqr + xD*xD + zSqr;
			if (distanceSqr > maxDistanceSqr) {
				continue;
			}
			// get weight
			int aux = (int) ((distanceSqr * iDeltaSqrt + 0.5)); //Same as ROUND but avoid comparison
			float tmpWeight = blobTableSqrt[aux];
		}
	}

}	



inline void processVoxel(int x, int y, int z, int halfVolSize, float transform[3][3], float maxDistanceSqr,
		std::complex<float>*** outputVolume, float*** outputWeights, std::complex<float>** inputImage, int imgSizeX, int imgSizeY) {
	float imgPos[3];
	int imgX, imgY;
	// transform current point to center
	imgPos[0] = x - halfVolSize;
	imgPos[1] = y - halfVolSize;
	imgPos[2] = z - halfVolSize;
	if (imgPos[0]*imgPos[0] + imgPos[1]*imgPos[1] + imgPos[2]*imgPos[2] > maxDistanceSqr) {
		return; // discard iterations that would access pixel with too high frequency
	}
	// rotate around center
	multiply(transform, imgPos);
	// transform back
	imgPos[1] += halfVolSize; // just Y coordinate, since X now match to picture and Z is irrelevant

	outputVolume[z][y][x] += getPixelClamped(inputImage, imgPos[0], imgPos[1], imgSizeX, imgSizeY);
	outputWeights[z][y][x] += 1.f;
}

inline void convert(Matrix2D<double>& in, float out[3][3]) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			out[i][j] = in(i, j);
		}
	}
}

void thirdAttempt(std::complex<float>*** outputVolume, float*** outputWeights, int outputVolumeSize,
	std::complex<float>** inputImage, int imgSizeX, int imgSizeY,
	float transform[3][3],
	float transformInv[3][3],
	Matrix1D<double>& blobTableSqrt, float iDeltaSqrt,
	float blobRadius)
{
	float halfOutVolSize = outputVolumeSize / 2.f;
	const int maxOutputIndex = outputVolumeSize-1;
	const float maxDistanceSqr = imgSizeX * imgSizeX;
	Point3D origin = {halfOutVolSize, halfOutVolSize, halfOutVolSize};

	// prepare traversing
	Point3D* plane = createProjectionPlane(imgSizeX, imgSizeY);
	plane = rotatePlane(plane, transform);
	plane = translatePlane(plane, origin);
	Point3D u, v;
	getVectors(plane, u, v);
	Point3D* AABB = getPlaneAABB(plane, 0, 0, 0, maxOutputIndex, maxOutputIndex, maxOutputIndex);
	int minY, minX, minZ;
	int maxY, maxX, maxZ;
	minZ = floor(AABB[0].z);
	minY = floor(AABB[0].y);
	minX = floor(AABB[0].x);
	maxZ = ceil(AABB[1].z);
	maxY = ceil(AABB[1].y);
	maxX = ceil(AABB[1].x);
	int dX, dY, dZ;
	dX = maxX - minX;
	dY = maxY - minY;
	dZ = maxZ - minZ;

	if (dZ <= dX && dZ <= dY) { // iterate XY plane
		for(int y = minY; y <= maxY; y++) {
			for(int x = minX; x <= maxX; x++) {
				float hitZ;
				if (getZ(x, y, hitZ, u, v, *plane)) {
					int z = (int)(hitZ + 0.5f); // rounding
					processVoxel(x, y, z, halfOutVolSize, transformInv, maxDistanceSqr, outputVolume, outputWeights, inputImage, imgSizeX, imgSizeY);
				}
			}
		}
	} else if (dY <= dX && dY <= dZ) { // iterate XZ plane
		for(int z = minZ; z <= maxZ; z++) {
			for(int x = minX; x <= maxX; x++) {
				float hitY;
				if (getY(x, hitY, z, u, v, *plane)) {
					int y = (int)(hitY + 0.5f); // rounding
					processVoxel(x, y, z, halfOutVolSize, transformInv, maxDistanceSqr, outputVolume, outputWeights, inputImage, imgSizeX, imgSizeY);
				}
			}
		}
	} else if(dX <= dY && dX <= dZ) { // iterate YZ plane
		for(int z = minZ; z <= maxZ; z++) {
			for(int y = minY; y <= maxY; y++) {
				float hitX;
				if (getX(hitX, y, z, u, v, *plane)) {
					int x = (int)(hitX + 0.5f); // rounding
					processVoxel(x, y, z, halfOutVolSize, transformInv, maxDistanceSqr, outputVolume, outputWeights, inputImage, imgSizeX, imgSizeY);
				}
			}
		}
	} else {
		std::cout << "ALERT" << std::endl;
	}
}

template<typename T>
T*** allocate(T***& where, int xSize, int ySize, int zSize) {
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
T** allocate(T**& where, int xSize, int ySize) {
	where = new T*[ySize];
	for (int y = 0; y < ySize; y++) {
		where[y] = new T[xSize];
		for (int x = 0; x < xSize; x++) {
			where[y][x] = (T) 0;
		}
	}
	return where;
}

template<typename T>
void release(T***& array, int xSize, int ySize) {
	for(int y = 0; y < ySize; y++) {
		for(int x = 0; x < xSize; x++) {
			delete[] array[y][x];
		}
		delete[] array[y];
	}
	delete[] array;
	array = NULL;
}

template<typename T>
void release(T**& array, int xSize) {
	for(int x = 0; x < xSize; x++) {
		delete[] array[x];
	}
	delete[] array;
	array = NULL;
}


template<typename T>
T*** applyBlob(T***& input, int size, float blobSize,
		Matrix1D<double>& blobTableSqrt, float iDeltaSqrt) {
	float blobSizeSqr = blobSize * blobSize;
	int xSize = size/2 + 1; // just half of the space is necessary, the rest is complex conjugate
	int blob = floor(blobSize); // we are using integer coordinates, so we cannot hit anything further
	T tmp;
	T*** output;
	// create new storage, notice that just 'right hand side' of the input will be used
	allocate(output, xSize, size+1, size+1);

	// traverse new storage
	for (int i = 0; i <= size; i++) {
		for (int j = 0; j <= size; j++) {
			for (int k = 0; k < xSize; k++) {
				// traverse input storage
				tmp = (T) 0;
				for (int z = std::max(0, i-blob); z < std::min(size+1, i+blob); z++) {
					float dZSqr = (i - z) * (i - z);
					for (int y = std::max(0, j-blob); y < std::min(size+1, j+blob); y++) {
						float dYSqr = (j - y) * (j - y);
						for (int x = std::max(0, k-blob); x < std::min(xSize, k+blob); x++) {
							float dXSqr = (k - x) * (k - x);
							float distanceSqr = dZSqr + dYSqr + dXSqr;
							if (distanceSqr > blobSizeSqr) {
								continue;
							}
							int aux = (int) ((distanceSqr * iDeltaSqrt + 0.5)); //Same as ROUND but avoid comparison
							float tmpWeight = blobTableSqrt[aux];
							tmp += tmpWeight * input[z][y][x+size/2];
						}
					}
				}
				output[i][j][k] = tmp;
			}
		}
	}
	// free original data
	release(input, size+1, size+1);
	return output;
}

void convertToExpectedSpace(std::complex<float>*** outputVolume, float*** outputWeight, int size,
	MultidimArray< std::complex<double> >& VoutFourier,
	MultidimArray<double>& FourierWeights) {
	int halfSize = size / 2;
	for (int z = 0; z <= size; z++) {
    	for (int y = 0; y <= size; y++) {
			for (int x = 0; x <= halfSize; x++) {
				if (outputWeight[z][y][x] == 0.f) {
					continue;
				}
    			int newPos[3];
    			// shift FFT from center to corners
				newPos[0] = x; // no need to move
    			newPos[1] = (y < halfSize) ? VoutFourier.ydim - halfSize + y : y - halfSize ;
    			newPos[2] = (z < halfSize) ? VoutFourier.zdim - halfSize + z : z - halfSize ;
    			// store to output array
				DIRECT_A3D_ELEM(VoutFourier, newPos[2], newPos[1], newPos[0]).real() += outputVolume[z][y][x].real();
				DIRECT_A3D_ELEM(VoutFourier, newPos[2], newPos[1], newPos[0]).imag() += outputVolume[z][y][x].imag();
    			DIRECT_A3D_ELEM(FourierWeights, newPos[2], newPos[1], newPos[0]) += outputWeight[z][y][x];
    		}
    	}
    }
}

void mirror(std::complex<float>*** outputVolume, float*** outputWeight, int size) {
	for (int z = 0; z < size; z++) {
		for (int y = 0; y < size; y++) {
			for (int x = 0; x < size/2; x++) {
				int newPos[3];
				// mirror against center of the volume, e.g. [0,0,0]->[size,size,size]. It will fit as the input space is one voxel biger
				newPos[0] = size - x;
				newPos[1] = size - y;
				newPos[2] = size - z;

				outputVolume[newPos[2]][newPos[1]][newPos[0]].real() += outputVolume[z][y][x].real();
				outputVolume[newPos[2]][newPos[1]][newPos[0]].imag() -= outputVolume[z][y][x].imag();
				outputWeight[newPos[2]][newPos[1]][newPos[0]] += outputWeight[z][y][x];
			}
		}
	}
}


//#define DEBUG
void ProgRecFourier::processImages( int firstImageIndex, int lastImageIndex, bool saveFSC, bool reprocessFlag)
{
    MultidimArray< std::complex<double> > *paddedFourier;

    int repaint = (int)ceil((double)SF.size()/60);

    bool processed;
    int imgno = 0;
    int imgIndex = firstImageIndex;

    // This index tells when to save work for later FSC usage
    int FSCIndex = (firstImageIndex + lastImageIndex)/2;

    // Index of the image that has just been processed. Used for
    // FSC purposes
    int current_index;
	std::complex<float>*** outputVolume;
	float*** outputWeight;
	std::complex<float>** myPaddedFourier;
	size_t size;
	size_t conserveRows;
    bool initialized = false;

    do
    {
        threadOpCode = PRELOAD_IMAGE;

        for ( int nt = 0 ; nt < numThreads ; nt ++ )
        {
            if ( imgIndex <= lastImageIndex )
            {
                th_args[nt].imageIndex = imgIndex;
                th_args[nt].reprocessFlag = reprocessFlag;
                imgIndex++;
            }
            else
            {
                th_args[nt].imageIndex = -1;
            }
        }

        // Awaking sleeping threads
        barrier_wait( &barrier );
        // here each thread is reading a different image and compute fft
        // Threads are working now, wait for them to finish
        // processing current projection
        barrier_wait( &barrier );

        // each threads have read a different image and now
        // all the thread will work in a different part of a single image.
        threadOpCode = PROCESS_IMAGE;

        processed = false;


        for ( int nt = 0 ; nt < numThreads ; nt ++ )
        {
            if ( th_args[nt].read == 2 )
                processed = true;
            else if ( th_args[nt].read == 1 )
            {
                processed = true;
                if (verbose && imgno++%repaint==0)
                    progress_bar(imgno);

                double weight = th_args[nt].localweight;
                paddedFourier = th_args[nt].localPaddedFourier;
                current_index = th_args[nt].imageIndex;
                Matrix2D<double> *Ainv = th_args[nt].localAInv;

                //#define DEBUG22
#ifdef DEBUG22

                {
                    static int ii=0;
                    if(ii%1==0)
                    {
                        FourierImage save22;
                        //save22()=*paddedFourier;
                        save22().alias(*paddedFourier);
                        save22.write((std::string) integerToString(ii)  + "_padded_fourier.spi");
                    }
                    ii++;
                }
#endif
                #undef DEBUG22

                // Initialized just once
                if ( statusArray == NULL )
                {
                    statusArray = (int *) malloc ( sizeof(int) * paddedFourier->ydim );
                }

                // Determine how many rows of the fourier
				// transform are of interest for us. This is because
				// the user can avoid to explore at certain resolutions
				conserveRows = (size_t) ceil((double) paddedFourier->ydim * maxResolution * 2.0);
				conserveRows = (size_t) ceil((double) conserveRows / 2.0);
				size = 2 * conserveRows;

				// allocate space for input image
				allocate(myPaddedFourier, conserveRows, size);
				// convert image (shift to center and remove high frequencies)
				std::complex<double> paddedFourierTmp;
				int halfY = paddedFourier->ydim / 2;
				double tempMyPadd[2];
				for (size_t i = 0; i < paddedFourier->ydim; i++) {
					for (size_t j = 0; j < conserveRows; j++) {
						if (i < conserveRows || i >= (paddedFourier->ydim - conserveRows)) {
							paddedFourierTmp = DIRECT_A2D_ELEM(*paddedFourier, i, j);
							FFT_IDX2DIGFREQ(j, paddedImg.xdim, tempMyPadd[0]);
							FFT_IDX2DIGFREQ(i, paddedImg.ydim, tempMyPadd[1]);
							if (tempMyPadd[0] * tempMyPadd[0] + tempMyPadd[1] * tempMyPadd[1]> maxResolution2) {
								continue;
							}
							int myPadI = (i < halfY) ?	i + conserveRows : i - paddedFourier->ydim + conserveRows;
							myPaddedFourier[myPadI][j].real() =	paddedFourierTmp.real();
							myPaddedFourier[myPadI][j].imag() =	paddedFourierTmp.imag();
						}
					}
				}

                if (!initialized) {
                	initialized = true;
					// the +1 is to prevent outOfBound reading when mirroring the result (later)
					allocate(outputVolume, size+1, size+1, size+1);
					allocate(outputWeight, size+1, size+1, size+1);
                }

                // Loop over all symmetries
                for (size_t isym = 0; isym < R_repository.size(); isym++)
                {
                    rowsProcessed = 0;

					// Compute the coordinate axes of the symmetrized projection
					Matrix2D<double> A_SL=R_repository[isym]*(*Ainv);
					Matrix2D<double> A_SLInv=A_SL.inv();
					float transf[3][3];
					float transfInv[3][3];
					convert(A_SL, transf);
					convert(A_SLInv, transfInv);

 ////////////////////////////////////////////////

#if DEBUG_DUMP > 0
					std::cout << "Img " << imgIndex << " symmetry " << isym << std::endl;
                    if (isym > 0 || current_index != 1) {
                       	continue;
                    }
#endif

//////////////////////////////////////////////////

                    // Fill the thread arguments for each thread
                    for ( int th = 0 ; th < numThreads ; th ++ )
                    {
                        // Passing parameters to each thread
                        th_args[th].symmetry = &A_SL;
                        th_args[th].paddedFourier = paddedFourier;
                        th_args[th].weight = weight;
                        th_args[th].reprocessFlag = reprocessFlag;
                    }

                    // Init status array
                    for (size_t i = 0 ; i < paddedFourier->ydim ; i ++ )
                    {
                        if ( i >= conserveRows && i < (paddedFourier->ydim-conserveRows))
                        {
                            // -2 means "discarded"
                            statusArray[i] = -2;
                            rowsProcessed++;
                        }
                        else
                        {
                            statusArray[i] = 0;
                        }
                    }

					thirdAttempt(outputVolume, outputWeight, size,
							myPaddedFourier, conserveRows, size, transf, transfInv,
							blobTableSqrt, iDeltaSqrt, blob.radius);

//                    // Awaking sleeping threads
//                    barrier_wait( &barrier );
//                    // Threads are working now, wait for them to finish
//                    // processing current projection
//                    barrier_wait( &barrier );

                    //#define DEBUG2
#ifdef DEBUG2

                    {
                        static int ii=0;
                        if(ii%1==0)
                        {
                            Image<double> save;
                            save().alias( FourierWeights );
                            save.write((std::string) integerToString(ii)  + "_1_Weights.vol");

                            Image< std::complex<double> > save2;
                            save2().alias( VoutFourier );
                            save2.write((std::string) integerToString(ii)  + "_1_Fourier.vol");
                        }
                        ii++;
                    }
#endif
                    #undef DEBUG2

                }

                if ( current_index == FSCIndex && saveFSC )
                {
                    // Save Current Fourier, Reconstruction and Weights
                    Image<double> save;
                    save().alias( FourierWeights );
                    save.write((std::string)fn_fsc + "_1_Weights.vol");

                    Image< std::complex<double> > save2;
                    save2().alias( VoutFourier );
                    save2.write((std::string) fn_fsc + "_1_Fourier.vol");

                    finishComputations(FileName((std::string) fn_fsc + "_1_recons.vol"));
                    Vout().initZeros(volPadSizeZ, volPadSizeY, volPadSizeX);
                    transformerVol.setReal(Vout());
                    Vout().clear();
                    transformerVol.getFourierAlias(VoutFourier);
                    FourierWeights.initZeros(VoutFourier);
                    VoutFourier.initZeros();
                }

                release(myPaddedFourier, size);
            }
        }



    }
    while ( processed );

#if DEBUG_DUMP > 0
    std::cout << "about to mirror" << std::endl;
#endif
    // get rid of the replicated data
    mirror(outputVolume, outputWeight, size);
#if DEBUG_DUMP > 0
    std::cout << "about to apply blob" << std::endl;
#endif
    outputVolume = applyBlob(outputVolume, size, blob.radius, blobTableSqrt, iDeltaSqrt);
    outputWeight = applyBlob(outputWeight, size, blob.radius, blobTableSqrt, iDeltaSqrt);
#if DEBUG_DUMP > 0
    std::cout << "about to convert to expected format" << std::endl;
#endif
    convertToExpectedSpace(outputVolume, outputWeight, size, VoutFourier, FourierWeights);
    release(outputVolume, size/2 + 1, size);
    release(outputWeight, size/2 + 1, size);
#if DEBUG_DUMP > 0
	std::cout << "done" << std::endl;
    print(VoutFourier, true);
    print(VoutFourier_muj, true);
#endif

    if( saveFSC )
    {
        // Save Current Fourier, Reconstruction and Weights
        Image<double> auxVolume;
        auxVolume().alias( FourierWeights );
        auxVolume.write((std::string)fn_fsc + "_2_Weights.vol");

        Image< std::complex<double> > auxFourierVolume;
        auxFourierVolume().alias( VoutFourier );
        auxFourierVolume.write((std::string) fn_fsc + "_2_Fourier.vol");

        finishComputations(FileName((std::string) fn_fsc + "_2_recons.vol"));

        Vout().initZeros(volPadSizeZ, volPadSizeY, volPadSizeX);
        transformerVol.setReal(Vout());
        Vout().clear();
        transformerVol.getFourierAlias(VoutFourier);
        FourierWeights.initZeros(VoutFourier);
        VoutFourier.initZeros();

        auxVolume.sumWithFile(fn_fsc + "_1_Weights.vol");
        auxVolume.sumWithFile(fn_fsc + "_2_Weights.vol");
        auxFourierVolume.sumWithFile(fn_fsc + "_1_Fourier.vol");
        auxFourierVolume.sumWithFile(fn_fsc + "_2_Fourier.vol");
        remove((fn_fsc + "_1_Weights.vol").c_str());
        remove((fn_fsc + "_2_Weights.vol").c_str());
        remove((fn_fsc + "_1_Fourier.vol").c_str());
        remove((fn_fsc + "_2_Fourier.vol").c_str());

        /*
        //Save SUM
                                    //this is an image but not an xmipp image
                                    auxFourierVolume.write((std::string)fn_fsc + "_all_Fourier.vol",
                                            false,VDOUBLE);
                                    auxVolume.write((std::string)fn_fsc + "_all_Weight.vol",
                                            false,VDOUBLE);
        //
        */
    }
}

void ProgRecFourier::correctWeight()
{
    // If NiterWeight=0 then set the weights to one
	forceWeightSymmetry(FourierWeights);
    if (NiterWeight==0)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(FourierWeights)
        DIRECT_A3D_ELEM(FourierWeights, k,i,j)=1;

    }
    else
    {
        // Temporary save the Fourier of the volume
        MultidimArray< std::complex<double> > VoutFourierTmp;
        VoutFourierTmp=VoutFourier;
        forceWeightSymmetry(FourierWeights);
        // Prepare the VoutFourier
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(VoutFourier)
        {
            double *ptrOut=(double *)&(DIRECT_A3D_ELEM(VoutFourier, k,i,j));
            if (fabs(A3D_ELEM(FourierWeights,k,i,j))>1e-3)
                ptrOut[0] = 1.0/DIRECT_A3D_ELEM(FourierWeights, k,i,j);
        }

        for (int i=1;i<NiterWeight;i++)
        {
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(FourierWeights)
            A3D_ELEM(FourierWeights,k,i,j)=0;
            processImages(0, SF.size() - 1, !fn_fsc.empty(), true);
            forceWeightSymmetry(FourierWeights);
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(VoutFourier)
            {
                double *ptrOut=(double *)&(DIRECT_A3D_ELEM(VoutFourier, k,i,j));
                if (fabs(A3D_ELEM(FourierWeights,k,i,j))>1e-3)
                    ptrOut[0] /= A3D_ELEM(FourierWeights,k,i,j);
            }
        }
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(VoutFourier)
        {
            // Put back the weights to FourierWeights from temporary variable VoutFourier
            double *ptrOut=(double *)&(DIRECT_A3D_ELEM(VoutFourier, k,i,j));
            A3D_ELEM(FourierWeights,k,i,j) = ptrOut[0];
        }
        VoutFourier = VoutFourierTmp;
    }
}

void ProgRecFourier::finishComputations( const FileName &out_name )
{
    //#define DEBUG_VOL
#ifdef DEBUG_VOL
    {
        VolumeXmipp save;
        save().alias( FourierWeights );
        save.write((std::string) fn_out + "Weights.vol");

        FourierVolume save2;
        save2().alias( VoutFourier );
        save2.write((std::string) fn_out + "FourierVol.vol");
    }
#endif

    // Enforce symmetry in the Fourier values as well as the weights
    // Sjors 19aug10 enforceHermitianSymmetry first checks ndim...
    Vout().initZeros(volPadSizeZ,volPadSizeY,volPadSizeX);
    transformerVol.setReal(Vout());
    transformerVol.enforceHermitianSymmetry();
    //forceWeightSymmetry(preFourierWeights);

    // Tell threads what to do
    //#define DEBUG_VOL1
#ifdef DEBUG_VOL1

    {
        Image<double> save;
        save().alias( FourierWeights );
        save.write((std::string) fn_out + "hermiticWeights.vol");

        Image< std::complex<double> > save2;
        save2().alias( VoutFourier );
        save2.write((std::string) fn_out + "hermiticFourierVol.vol");
    }
#endif
    threadOpCode = PROCESS_WEIGHTS;
    // Awake threads
    barrier_wait( &barrier );
    // Threads are working now, wait for them to finish
    barrier_wait( &barrier );

    transformerVol.inverseFourierTransform();
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
        if (NiterWeight!=0)
        {
            A3D_ELEM(mVout,k,i,j) /= (ipad_relation*factor2*factor);
            meanFactor2+=factor2;
        }
        else
            A3D_ELEM(mVout,k,i,j) /= (ipad_relation*factor);
    }
    if (NiterWeight!=0)
    {
        meanFactor2/=MULTIDIM_SIZE(mVout);
        FOR_ALL_ELEMENTS_IN_ARRAY3D(mVout)
        A3D_ELEM(mVout,k,i,j) *= meanFactor2;
    }
    Vout.write(out_name);
}

void ProgRecFourier::setIO(const FileName &fn_in, const FileName &fn_out)
{
    this->fn_sel = fn_in;
    this->fn_out = fn_out;
}

void ProgRecFourier::forceWeightSymmetry(MultidimArray<double> &FourierWeights)
{
    int yHalf=YSIZE(FourierWeights)/2;
    if (YSIZE(FourierWeights)%2==0)
        yHalf--;
    int zHalf=ZSIZE(FourierWeights)/2;
    if (ZSIZE(FourierWeights)%2==0)
        zHalf--;
    int zsize=(int)ZSIZE(FourierWeights);
    int zsize_1=zsize-1;
    int ysize_1=(int)YSIZE(FourierWeights)-1;
    for (int k=0; k<zsize; k++)
    {
        int ksym=intWRAP(-k,0,zsize_1);
        for (int i=1; i<=yHalf; i++)
        {
            int isym=intWRAP(-i,0,ysize_1);
            double mean=0.5*(
                            DIRECT_A3D_ELEM(FourierWeights,k,i,0)+
                            DIRECT_A3D_ELEM(FourierWeights,ksym,isym,0));
            DIRECT_A3D_ELEM(FourierWeights,k,i,0)=
                DIRECT_A3D_ELEM(FourierWeights,ksym,isym,0)=mean;
        }
    }
    for (int k=1; k<=zHalf; k++)
    {
        int ksym=intWRAP(-k,0,zsize_1);
        double mean=0.5*(
                        DIRECT_A3D_ELEM(FourierWeights,k,0,0)+
                        DIRECT_A3D_ELEM(FourierWeights,ksym,0,0));
        DIRECT_A3D_ELEM(FourierWeights,k,0,0)=
            DIRECT_A3D_ELEM(FourierWeights,ksym,0,0)=mean;
    }
}
