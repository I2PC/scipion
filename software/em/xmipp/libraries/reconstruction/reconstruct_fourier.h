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

#ifndef __RECONSTRUCT_FOURIER_H
#define __RECONSTRUCT_FOURIER_H

#include <iostream>
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
#include "recons.h"
#include <reconstruction/directions.h>
#include <reconstruction/symmetrize.h>
#define BLOB_TABLE_SIZE 5000
#define BLOB_TABLE_SIZE_SQRT 10000

#define MINIMUMWEIGHT 0.001
#define ACCURACY 0.001

#define EXIT_THREAD 0
#define PRELOAD_IMAGE 1

/**@defgroup FourierReconstruction Fourier reconstruction
   @ingroup ReconsLibrary */
//@{
class ProgRecFourier;

/** Struct representing all data regarding one projection */
struct ProjectionData
{
	Array2D<std::complex<float> >* img;
	CTFDescription ctf;
	int imgIndex;
	float weight;
	Matrix2D<double> localAInv;
	bool skip;
};

/** Struct represents a point in 3D */
struct Point3D {
	float x, y, z;
};

/** Struct holding information for loading thread */
struct LoadThreadParams
{
    pthread_t id;
    ProgRecFourier * parent;
    int startImageIndex;
    int endImageIndex;
    MetaData* selFile;
    ProjectionData* buffer1 = NULL;
    ProjectionData* buffer2 = NULL;
};

class ProgRecFourier : public ProgReconsBase
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
    LoadThreadParams loadThread;
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

//	METHODS
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
	 * Method will take temp spaces (containing complex conjugate values
	 * in the 'right X side'), transfer them to 'left X side' and remove
	 * the 'right X side'. As a result, the X dimension of the temp spaces
	 * will be half of the original.
	 */
    void mirrorAndCropTempSpaces();
    /**
     * Method will enforce Hermitian symmetry, i.e will remove make sure
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
     * Method will create thread used for loading images
     * and set all necessary values/variables
     */
    void createLoadingThread();
    /**
     * Method will release all resources allocated for loading images
     */
    void cleanLoadingThread();
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

private:
//    FIELDS
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
    /** Sampling rate in Angstroms/pixel */
    double sampling_rate;
    /** Max resolution in Angstroms */
    double maxResolution;
    /** Maximum interesting resolution squared */
    double maxResolutionSqr;
    /** Tells the loading thread what to do next */
    int threadOpCode;
    /** Barrier synchronization for threads */
    barrier_t barrier;
    /** Table with blob values, linear sampling */
    Matrix1D<double>  Fourier_blob_table;
    /** Table with blob values, squared sampling */
    Matrix1D<double> blobTableSqrt;
    /** Inverse of the delta and deltaFourier used in the tables */
    double iDeltaFourier, iDeltaSqrt;
    /** Definition of the blob */
    struct blobtype blob;
    /** Vector with R symmetry matrices */
    std::vector <Matrix2D<double> > R_repository;
    /** Size of the original projection, must be a square */
    int imgSize;
    /** Size of the image including padding. Image must be a square */
    int paddedImgSize;
    /** Size of loading buffer (i.e. number of projection loaded in one buffer) */
    int bufferSize;

// STATIC METHODS
    /** Method to allocate 3D array (not continuous) of given size */
    template<typename T>
    static T*** allocate(T***& where, int xSize, int ySize, int zSize);
    /** Method to allocate 2D array (not continuous) of given size */
    template<typename T>
    static T** allocate(T**& where, int xSize, int ySize);
    /** Method running on separate thread */
	static void* loadImageThread(void* threadArgs);
	/** Function behaving like an identity, i.e returning passed value */
	template<typename T>
	static T identity(T val) { return val;}; // if used with some big type, use reference
	/** Function returning conjugate of a complex number */
	template<typename T>
	static std::complex<T> conjugate(std::complex<T> f) { return conj(f);};
    /**
     * Method will process the 'paddedFourier' (not shifted, i.e. low frequencies are in corners)
     * in the following way:
     * high frequencies are skipped (replaced by zero (0))
     * space is shifted, so that low frequencies are in the middle of the Y axis
     * resulting space is cropped.
     * Method returns a 2D array with Fourier coefficients, shifted so that low frequencies are
     * in the center of the Y axis (i.e. semicircle)
     */
    static Array2D<std::complex<float> >* cropAndShift(
    		MultidimArray<std::complex<double> >& paddedFourier,
    		ProgRecFourier* parent);
    /** Returns value within the range (included) */
    template<typename T, typename U>
    static U clamp(U val, T min, T max);
    /**
    *          7____6
    *         3/___2/
    *    +    | |  ||   y
    * [0,0,0] |*|4 ||5
    *    -    |/___|/  z  sizes are padded with blob-radius
    *        0   x  1
    * [0,0] is in the middle of the left side (point [0] and [3]), provided the blobSize is 0
    * otherwise the values can go to negative values
    */
    static void createProjectionCuboid(Point3D* cuboid, float sizeX, float sizeY, float blobSize);
    //* Apply rotation transform to cuboid */
    static void rotateCuboid(Point3D* cuboid, const float transform[3][3]) {
    	for (int i = 0; i < 8; i++) {
    		multiply(transform, cuboid[i]);
    	}
    }
    /** Do 3x3 x 1x3 matrix-vector multiplication */
    static void multiply(const float transform[3][3], Point3D& inOut);
    /** Add 'vector' to each element of 'cuboid' */
    static void translateCuboid(Point3D* cuboid, Point3D vector);
    /**
     * Method will calculate Axis Aligned Bound Box of the cuboid and restrict
     * its maximum size
     */
    static void computeAABB(Point3D* AABB, Point3D* cuboid,
    		float minX, float minY, float minZ,
    		float maxX, float maxY, float maxZ);

// METHODS
	/** Method will set indexes of the images to load and open sync barrier */
    void loadImages(int startIndex, int endIndex);
    /** Method swaps buffers of the loading thread */
    void swapLoadBuffers();
    /**
     * Method will use data stored in the buffer and update temporal
     * storages appropriately.
     */
    void processBuffer(ProjectionData* buffer);
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
    		Matrix1D<double>& blobTableSqrt, float iDeltaSqrt);
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
     * Method will process one projection image and add result to temporal
     * spaces. Method also shows progress of the calculation.
     */
    void processProjection(
    	ProjectionData* projectionData,
    	const float transform[3][3],
    	const float transformInv[3][3]);
    /**
     * Method will map one voxel from the temporal
     * spaces to the given projection and update temporal spaces
     * using the pixel value of the projection.
     */
    void processVoxel(
    		int x, int y, int z,
			const float transform[3][3], float maxDistanceSqr,
    		ProjectionData* const data);


    void processVoxelBlob(int x, int y, int z, const float transform[3][3], float maxDistanceSqr,
    		ProjectionData* const data);


    /**
     * Method sets wCTF and wModulator based on position
     * in the image
     */
    void processCTF(ProjectionData* const data,
    		int imgX, int imgY,
			float& wCTF, float& wModulator);

//    /**
//     * Method will return a pixel on given position (rounded to closest int,
//     * clamped to img size)
//     */
//    std::complex<float>
//    getPixelClamped(Array2D<std::complex<float> >* img, float x, float y);
    //    inline void allocateFourierWeights() {
//    	if (NULL == FourierWeights.data) {
//    		FourierWeights.initZeros(paddedImgSize, paddedImgSize, paddedImgSize/2 +1);
//    	}
//    }
//
//    void saveIntermediateFSC(std::complex<float>*** outputVolume, float*** outputWeight, int volSize,
//    		const std::string &weightsFileName, const std::string &fourierFileName, const FileName &resultFileName,
//			bool storeResult);
//    void saveFinalFSC(std::complex<float>*** outputVolume, float*** outputWeight, int volSize);
//
//    template<typename T>
//    void crop(T***& inOut, int size);
//
//	static void processCube(
//			int i,
//			int j,
//			const Matrix1D<int>& corner1,
//			const Matrix1D<int>& corner2,
//			const MultidimArray<float>& z2precalculated,
//			const MultidimArray<int>& zWrapped,
//			const MultidimArray<int>& zNegWrapped,
//			const MultidimArray<float>& y2precalculated,
//			float blobRadiusSquared,
//			const MultidimArray<int>& yWrapped,
//			const MultidimArray<int>& yNegWrapped,
//			const MultidimArray<float>& x2precalculated,
//			float iDeltaSqrt,
//			float wModulator,
//			const MultidimArray<int>& xWrapped,
//			int xsize_1,
//			const MultidimArray<int>& xNegWrapped,
//			bool reprocessFlag,
//			float wCTF,
//			MultidimArray<std::complex<double> >& VoutFourier,
//			Matrix1D<double>& blobTableSqrt,
//			LoadThreadParams* threadParams,
//			MultidimArray<double>& fourierWeights,
//			double* ptrIn,
//			float weight,
//			ProgRecFourier * parent,
//			Matrix1D<double>& real_position);
};
//@}
#endif
