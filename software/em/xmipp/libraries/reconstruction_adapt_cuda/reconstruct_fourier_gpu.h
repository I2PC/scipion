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

#ifndef __RECONSTRUCT_FOURIER_GPU_H
#define __RECONSTRUCT_FOURIER_GPU_H

#include <iostream>
#include <limits>
#include <data/xmipp_fftw.h>
#include <data/xmipp_funcs.h>
#include <data/xmipp_image.h>
#include <data/projection.h>
#include <data/xmipp_threads.h>
#include <data/blobs.h>
#include <data/metadata.h>
#include <data/ctf.h>
#include <data/array_2D.h>
#include <data/args.h>
#include <data/xmipp_fft.h>
#include <sys/time.h>
#include <data/metadata.h>
#include <reconstruction/recons.h>
#include <reconstruction/directions.h>
#include <reconstruction/symmetrize.h>
#include <data/point3D.h>
#include <reconstruction/reconstruct_fourier_defines.h>
#include <reconstruction/reconstruct_fourier_projection_traverse_space.h>
#include <reconstruction_cuda/cuda_gpu_reconstruct_fourier.h>

/**@defgroup FourierReconstruction Fourier reconstruction
   @ingroup ReconsLibrary */
//@{
class ProgRecFourierGPU;

/** Struct representing the working thread */
struct RecFourierWorkThread
{
    pthread_t id;
    ProgRecFourierGPU * parent;
    int startImageIndex; // index of the first projection to process
    int endImageIndex; // index of the last projection to process
    MetaData* selFile; // used for loading data
    RecFourierBufferData* buffer; // where data are loaded
    int gpuStream; // index of stream on GPU device
};

class ProgRecFourierGPU : public ProgReconsBase
{
public:
    /**
     * Run the image processing.
     * Method will load data, process them and store result to final destination.
     */
    void run();

    /**
     * Method will take data stored in tempVolume and tempWeights
     * (which should be cropped in X axis), calculate IFFT and store result
     * to final destination.
     */
    void finishComputations( const FileName &out_name );

    /** Functions of common reconstruction interface */
    virtual void setIO(const FileName &fn_in, const FileName &fn_out);

protected:
// 	FIELDS

    /** Thread used for loading input data */
    RecFourierWorkThread* workThreads;

    /** SelFile containing all projections */
    MetaData SF;

    /** Output file name */
    FileName fn_out;

    /** Input file name */
    FileName fn_in;

    /** maximal index in the temporal volumes, Y and Z axis */
	int maxVolumeIndexYZ;

    /** maximal index in the temporal volumes, X axis */
	int maxVolumeIndexX;

	/**
	 * 3D volume holding the cropped (without high frequencies) Fourier space.
	 * Lowest frequencies are in the center, i.e. Fourier space creates a
	 * sphere within a cube.
	 */
	std::complex<float>*** tempVolume = NULL;

	/**
	 * 3D volume holding the weights of the Fourier coefficients stored
	 * in tempVolume.
	 */
	float*** tempWeights = NULL;

	/**
	 * 3D volume holding the cropped (without high frequencies) Fourier space.
	 * Lowest frequencies are in the center, i.e. Fourier space creates a
	 * sphere within a cube.
	 * Valid for GPU, i.e. it is equivalent to tempVolume, which has been flattened
	 * Each two successive values represent one complex number
	 */
	float* tempVolumeGPU = NULL;

	/**
	 * 3D volume holding the weights of the Fourier coefficients stored
	 * in tempVolume.
	 * Valid for GPU, i.e. it is equivalent to tempWeights, which has been flattened
	 */
	float* tempWeightsGPU = NULL;


    /** Holds number of cores available at the host system */
    int noOfThreads;

    /** CUDA device to use */
    int device;

//	METHODS

	/**
	 * Method checks that there is no logical problem in the defines used by
	 * program. If found, error is thrown.
	 */
	void checkDefines();

	/**
	 * Method will take temp spaces (containing complex conjugate values
	 * in the 'right X side'), transfer them to 'left X side' and remove
	 * the 'right X side'. As a result, the X dimension of the temp spaces
	 * will be half of the original.
	 */
    void mirrorAndCropTempSpaces();

	/**
	 * Method will fill temporal spaces with data stored at the GPU.
	 */
    void getGPUData();

    /**
     * Method will enforce Hermitian symmetry, i.e will make sure
     * that the values in temporal space at X=0 are complex conjugate of
     * in respect to center of the space
     */
    void forceHermitianSymmetry();

    /**
     * Method will in effect do the point-wise division of
     * tempVolume and tempWeights
     * (i.e. correct Fourier coefficients by proper weight)
     */
    void processWeights();

    /**
     * Method will create thread used for loading and processing images
     * Thread starts the work immediately
     * gpuStream - stream to use
     * startIndex - index of the first projection to process
     * endIndex - index of the last projection to process
     * thread - used for referencing the thread
     */
    void createWorkThread(int gpuStream, int startIndex, int endIndex, RecFourierWorkThread& thread);

    /** Read arguments from command line */
    void readParams();

    /** Specify supported command line arguments */
    void defineParams();

    /** Show basic info to standard output */
    void show();

    /** Method will fill other help structures and variables */
    void produceSideinfo();

    /**
     * Method will initialize temporal storage (if necessary),
     *  processes images in given interval (including boundary indexes)
     *  and store results in temporal storage (@see tempVolume and tempWeights)
     */
    void processImages(int firstImageIndex, int lastImageIndex);

    /** Method will release temporal spaces for weights and Fourier coefs. */
    void releaseTempSpaces();

private:
//    FIELDS

    /** variable used for blob table values calculation */
    double iw0;

    /** File with symmetries */
    FileName fn_sym;

    /** Flag whether to use the weights in the image metadata */
    bool do_weights;

    /** If true, blobing is done at the end of the computation */
    bool useFast;

    /** Use CTF */
    bool useCTF;

    /** True if the images have been already phase flipped */
    bool isPhaseFlipped;

    /** Minimum CTF value to invert */
    double minCTF;

    /** Inverse of the sampling rate */
    double iTs;

    /** Projection padding Factor */
    double padding_factor_proj;

    /** Volume padding Factor */
    double padding_factor_vol;

    /** Max resolution in Angstroms */
    float maxResolution;

    /** Maximum interesting resolution squared */
    float maxResolutionSqr;

    /** Barrier synchronization for threads */
    barrier_t barrier;

    /** Table with blob values, linear sampling */
    Matrix1D<double>  Fourier_blob_table;

    /** Table with blob values, squared sampling */
    float blobTableSqrt[BLOB_TABLE_SIZE_SQRT];

    /** Inverse of the delta and deltaFourier used in the tables */
    float iDeltaFourier, iDeltaSqrt;

    /** Definition of the blob */
    struct blobtype blob;

    /** Vector with R symmetry matrices */
    std::vector <Matrix2D<double> > R_repository;

    /** Size of the original projection, must be a square */
    int imgSize;

    /** Size of the image including padding. Image must be a square */
    int paddedImgSize;

    /** Size of loading buffer (i.e. max number of projection loaded in one buffer) */
    int bufferSize;

    /** If set to true, FFT of the input projections shall be done on GPU */
    bool fftOnGPU;

// STATIC METHODS

    /** Method to allocate 3D array (not continuous) of given size */
    template<typename T>
    static T*** allocate(T***& where, int xSize, int ySize, int zSize);

    /** Method to release 3D array of given size */
    template<typename T>
    void release(T***& array, int ySize, int zSize);

    /** Method running on separate thread, loading images and processing them on GPU */
	static void* threadRoutine(void* threadArgs);

	/** Function behaving like an identity, i.e returning passed value */
	template<typename T>
	static T identity(T val) { return val;}; // if used with some big type, use reference

	/** Function returning conjugate of a complex number */
	template<typename T>
	static std::complex<T> conjugate(std::complex<T> f) { return conj(f);};

	static float getBessiOrderAlpha(blobtype blob) {
		switch (blob.order) {
		case 0: return bessi0(blob.alpha);
		case 1: return bessi1(blob.alpha);
		case 2: return bessi2(blob.alpha);
		case 3: return bessi3(blob.alpha);
		case 4: return bessi4(blob.alpha);
		REPORT_ERROR(ERR_VALUE_INCORRECT,"Order must be in interval [0..4]");
		}
	}

    /**
     * Method will process the 'paddedFourier' (not shifted, i.e. low frequencies are in corners)
     * in the following way:
     * high frequencies are skipped (replaced by zero (0))
     * space is shifted, so that low frequencies are in the middle of the Y axis
     * resulting space is cropped.
     * Method returns a 2D array with Fourier coefficients, shifted so that low frequencies are
     * in the center of the Y axis (i.e. semicircle)
     */
    static void cropAndShift(
    		MultidimArray<std::complex<double> >& paddedFourier,
    		ProgRecFourierGPU* parent,
			RecFourierBufferData* buffer,
			float* dest);

    /**
    *          6____5
    *         2/___1/
    *    +    | |  ||   y
    * [0,0,0] |*|7 ||4
    *    -    |/___|/  z  sizes are padded with blob-radius
    *         3  x  0
    * [0,0] is in the middle of the left side (point [2] and [3]), provided the blobSize is 0
    * otherwise the values can go to negative values
    * origin(point[0]) in in 'high' frequencies so that possible numerical instabilities are moved to high frequencies
    */
    static void createProjectionCuboid(Point3D<float>* cuboid, float sizeX, float sizeY, float blobSize);

    //* Apply rotation transform to cuboid */
    static void rotateCuboid(Point3D<float>* cuboid, const float transform[3][3]) {
    	for (int i = 0; i < 8; i++) {
    		multiply(transform, cuboid[i]);
    	}
    }

    /** Do 3x3 x 1x3 matrix-vector multiplication */
    static void multiply(const float transform[3][3], Point3D<float>& inOut);

    /** Add 'vector' to each element of 'cuboid' */
    static void translateCuboid(Point3D<float>* cuboid, Point3D<float> vector);

    /**
     * Method will calculate Axis Aligned Bound Box of the cuboid and restrict
     * its maximum size
     */
    static void computeAABB(Point3D<float>* AABB, Point3D<float>* cuboid,
    		float minX, float minY, float minZ,
    		float maxX, float maxY, float maxZ);

    /** DEBUG ONLY method, prints AABB to std::cout. Output can be used in e.g. GNUPLOT */
    static void printAABB(Point3D<float> AABB[]);

    /** Method to convert temporal space to expected (original) format */
    template<typename T, typename U>
    static void convertToExpectedSpace(T*** input, int size,
    	MultidimArray<U>& VoutFourier);

    /** Method to load a buffer of images from input file */
    static void prepareBuffer(RecFourierWorkThread * threadParams,
    		ProgRecFourierGPU* parent,
    		bool hasCTF, std::vector<size_t>& objId);

    /**
     * Method computes CTF and weight modulator for each pixel in the image
     */
    static void computeCTFCorrection(RecFourierWorkThread* threadParams,
    		size_t imgIndex,
			ProgRecFourierGPU* parent,
			RecFourierBufferData* buffer,
    		int storeIndex);

// METHODS

    /** method will parse number of threads from CMD and set accordingly
     *  or default value will be used (all available threads)
     */
    void parseNoOfThreads();

    /**
	 * Method will take input array (of size
	 * maxVolumeIndexYZ*maxVolumeIndexYZ*maxVolumeIndexYZ
	 * and transfer data to newly created array of size
	 * maxVolumeIndexX*maxVolumeIndexYZ*maxVolumeIndexYZ (i.e. X side is half).
	 * so that data on the 'left side' are mirrored against the origin of the
	 * array, and 'f' function is applied on them.
	 * 'right hand side' is only transfered.
	 * As a result, input will contain an array with X side half of the original
	 * size, with all data transfered from 'missing' side to 'preserved' side
	 */
    template<typename T>
    void mirrorAndCrop(T***& input,T (*f)(T));

    /**
     * Method will apply blob to input 3D array.
     * Original array will be released and new (blurred) will be returned
     */
    template<typename T>
    T*** applyBlob(T***& input, float blobSize,
    		float* blobTableSqrt, float iDeltaSqrt);

    /**
     * Method will allocate space for output Fourier transformation.
     * If space is already allocated, method will have no effect
     */
    void allocateVoutFourier(MultidimArray<std::complex<double> >&VoutFourier) {
    	if ((NULL == VoutFourier.data) || (0 == VoutFourier.getSize())) {
    		VoutFourier.initZeros(paddedImgSize, paddedImgSize, paddedImgSize/2 +1);
    	}
    }

    /**
     * Method calculates a traversal space information for specific projection
     * imgSizeX - X size of the projection
     * imgSizeY - Y size of the projection
     * projectionIndex - index to some array where the respective projection is stored
     * transform - forward rotation that should be applied to the projection
     * transformInv - inverse transformation
     * space - which will be filled
     */
    void computeTraverseSpace(int imgSizeX, int imgSizeY, int projectionIndex,
    		MATRIX& transform, MATRIX& transformInv, RecFourierProjectionTraverseSpace* space);

    /**
     * Method logs than 'increment' more pictures were processed
     * Thread safe
     */
    void logProgress(int increment);
};
//@}
#endif
