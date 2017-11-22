#include "/home/david/GIT/Scipion/software/em/xmipp/libraries/data/reconstruct_fourier_projection_traverse_space.h"
#include "/home/david/GIT/Scipion/software/em/xmipp/libraries/data/point3D.h"

__device__ __constant__ float cBlobRadius = 1.900000f;
__device__ __constant__ float cIDeltaSqrt = 2769.806152f;


/** Returns value within the range (included) */
template<typename T, typename U>
__device__
static U clamp(U val, T min, T max) {
	U res = (val > max) ? max : val;
	return (res < min) ? min : res;
}

__device__
const float* getNthItem(const float* array, int itemIndex, int fftSizeX, int fftSizeY, bool isFFT) {
		return array + (fftSizeX * fftSizeY * itemIndex * (isFFT ? 2 : 1)); // *2 since it's complex
}


/**
 * Calculates Z coordinate of the point [x, y] on the plane defined by p0 (origin) and normal
 */
__device__
float getZ(float x, float y, float nX, float nY, float nZ, float p0X, float p0Y, float p0Z) {
	// from a(x-x0)+b(y-y0)+c(z-z0)=0
	return (-nX*(x-p0X)-nY*(y-p0Y))/nZ + p0Z;
}

/**
 * Calculates Y coordinate of the point [x, z] on the plane defined by p0 (origin) and normal
 */
__device__
float getY(float x, float z, float nX, float nY, float nZ, float p0X, float p0Y, float p0Z) {
	// from a(x-x0)+b(y-y0)+c(z-z0)=0
	return (-nX*(x-p0X)-nZ*(z-p0Z))/nY + p0Y;
}

/**
 * Calculates X coordinate of the point [y, z] on the plane defined by p0 (origin) and normal
 */
__device__
float getX(float y, float z, float nX, float nY, float nZ, float p0X, float p0Y, float p0Z) {
	// from a(x-x0)+b(y-y0)+c(z-z0)=0
	return (-nY*(y-p0Y)-nZ*(z-p0Z))/nX + p0X;
}

/** Do 3x3 x 1x3 matrix-vector multiplication */
__device__
void multiply(float t00, float t01, float t02, float t10, float t11, float t12, float t20, float t21, float t22, Point3D<float>& inOut) {
	float tmp0 = t00 * inOut.x + t01 * inOut.y + t02 * inOut.z;
	float tmp1 = t10 * inOut.x + t11 * inOut.y + t12 * inOut.z;
	float tmp2 = t20 * inOut.x + t21 * inOut.y + t22 * inOut.z;
	inOut.x = tmp0;
	inOut.y = tmp1;
	inOut.z = tmp2;
}

/** Do 3x3 x 1x3 matrix-vector multiplication */
__device__
void multiply(const RecFourierProjectionTraverseSpace* tSpace, Point3D<float>& inOut) {
	multiply(tSpace->transformInv00, tSpace->transformInv01, tSpace->transformInv02,
			tSpace->transformInv10, tSpace->transformInv11, tSpace->transformInv12,
			tSpace->transformInv20, tSpace->transformInv21, tSpace->transformInv22, inOut);
}

/**
 * Method will map one voxel from the temporal
 * spaces to the given projection and update temporal spaces
 * using the pixel values of the projection withing the blob distance.
 */
__device__
void processVoxelBlob(
	float* tempVolumeGPU, float *tempWeightsGPU,
	const float* FFTs,
	int fftSizeX, int fftSizeY,
	int x, int y, int z,
	const RecFourierProjectionTraverseSpace* const space,
	const float* blobTableSqrt,
	int cMaxVolumeIndexX, int cMaxVolumeIndexYZ)
{
	Point3D<float> imgPos;
	// transform current point to center
	imgPos.x = x - cMaxVolumeIndexX/2;
	imgPos.y = y - cMaxVolumeIndexYZ/2;
	imgPos.z = z - cMaxVolumeIndexYZ/2;
	if ((imgPos.x*imgPos.x + imgPos.y*imgPos.y + imgPos.z*imgPos.z) > space->maxDistanceSqr) {
		return; // discard iterations that would access pixel with too high frequency
	}
	// rotate around center
	multiply(space, imgPos);
	if (imgPos.x < -cBlobRadius) return; // reading outside of the image boundary. Z is always correct and Y is checked by the condition above
	// transform back just Y coordinate, since X now matches to picture and Z is irrelevant
	imgPos.y += cMaxVolumeIndexYZ / 2;

	// check that we don't want to collect data from far far away ...
	float radiusSqr = cBlobRadius * cBlobRadius;
	float zSqr = imgPos.z * imgPos.z;
	if (zSqr > radiusSqr) return;

	// create blob bounding box
	int minX = ceilf(imgPos.x - cBlobRadius);
	int maxX = floorf(imgPos.x + cBlobRadius);
	int minY = ceilf(imgPos.y - cBlobRadius);
	int maxY = floorf(imgPos.y + cBlobRadius);
	minX = fmaxf(minX, 0);
	minY = fmaxf(minY, 0);
	maxX = fminf(maxX, fftSizeX-1);
	maxY = fminf(maxY, fftSizeY-1);

	int index3D = z * (cMaxVolumeIndexYZ+1) * (cMaxVolumeIndexX+1) + y * (cMaxVolumeIndexX+1) + x;
	float volReal, volImag, w;
	volReal = volImag = w = 0.f;
	const float* __restrict__ img = getNthItem(FFTs, space->projectionIndex, fftSizeX, fftSizeY, true);
	float dataWeight = space->weight;

	// check which pixel in the vicinity should contribute
	for (int i = minY; i <= maxY; i++) {
		float ySqr = (imgPos.y - i) * (imgPos.y - i);
		float yzSqr = ySqr + zSqr;
		if (yzSqr > radiusSqr) continue;
		for (int j = minX; j <= maxX; j++) {
			float xD = imgPos.x - j;
			float distanceSqr = xD*xD + yzSqr;
			if (distanceSqr > radiusSqr) continue;

			int index2D = i * fftSizeX + j;
			int aux = (int) ((distanceSqr * cIDeltaSqrt + 0.5f));
			float wBlob = blobTableSqrt[aux];
			float weight = wBlob * dataWeight;

			w += weight;
			volReal += img[2*index2D] * weight;
			volImag += img[2*index2D + 1] * weight;
		}
	}
	// use atomic as two blocks can write to same voxel
	atomicAdd(&tempVolumeGPU[2*index3D], volReal);
	atomicAdd(&tempVolumeGPU[2*index3D + 1], volImag);
	atomicAdd(&tempWeightsGPU[index3D], w);
}

/**
  * Method will process one projection image and add result to temporal
  * spaces.
  */
__device__
void processProjection(
	float* tempVolumeGPU, float *tempWeightsGPU,
	const float* FFTs,
	int fftSizeX, int fftSizeY,
	const RecFourierProjectionTraverseSpace* const tSpace,
	const float* devBlobTableSqrt,
	int cMaxVolumeIndexX, int cMaxVolumeIndexYZ)
{
	// map thread to each (2D) voxel
	volatile int idx = blockIdx.x*blockDim.x + threadIdx.x;
	volatile int idy = blockIdx.y*blockDim.y + threadIdx.y;

	if (tSpace->XY == tSpace->dir) { // iterate XY plane
		if (idy >= tSpace->minY && idy <= tSpace->maxY) {
			if (idx >= tSpace->minX && idx <= tSpace->maxX) {
				float z1 = getZ(idx, idy, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ); // lower plane
				float z2 = getZ(idx, idy, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ); // upper plane
				z1 = clamp(z1, 0, cMaxVolumeIndexYZ);
				z2 = clamp(z2, 0, cMaxVolumeIndexYZ);
				int lower = floorf(fminf(z1, z2));
				int upper = ceilf(fmaxf(z1, z2));
				for (int z = lower; z <= upper; z++) {
					processVoxelBlob(tempVolumeGPU, tempWeightsGPU, FFTs, fftSizeX, fftSizeY, idx, idy, z, tSpace, devBlobTableSqrt, cMaxVolumeIndexX, cMaxVolumeIndexYZ);
				}
			}
		}
	} else if (tSpace->XZ == tSpace->dir) { // iterate XZ plane
		if (idy >= tSpace->minZ && idy <= tSpace->maxZ) { // map z -> y
			if (idx >= tSpace->minX && idx <= tSpace->maxX) {
				float y1 = getY(idx, idy, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ); // lower plane
				float y2 = getY(idx, idy, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ); // upper plane
				y1 = clamp(y1, 0, cMaxVolumeIndexYZ);
				y2 = clamp(y2, 0, cMaxVolumeIndexYZ);
				int lower = floorf(fminf(y1, y2));
				int upper = ceilf(fmaxf(y1, y2));
				for (int y = lower; y <= upper; y++) {
					processVoxelBlob(tempVolumeGPU, tempWeightsGPU, FFTs, fftSizeX, fftSizeY, idx, y, idy, tSpace, devBlobTableSqrt, cMaxVolumeIndexX, cMaxVolumeIndexYZ);
				}
			}
		}
	} else { // iterate YZ plane
		if (idy >= tSpace->minZ && idy <= tSpace->maxZ) { // map z -> y
			if (idx >= tSpace->minY && idx <= tSpace->maxY) { // map y > x
				float x1 = getX(idx, idy, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ); // lower plane
				float x2 = getX(idx, idy, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ); // upper plane
				x1 = clamp(x1, 0, cMaxVolumeIndexX);
				x2 = clamp(x2, 0, cMaxVolumeIndexX);
				int lower = floorf(fminf(x1, x2));
				int upper = ceilf(fmaxf(x1, x2));
				for (int x = lower; x <= upper; x++) {
					processVoxelBlob(tempVolumeGPU, tempWeightsGPU, FFTs, fftSizeX, fftSizeY, x, idx, idy, tSpace, devBlobTableSqrt, cMaxVolumeIndexX, cMaxVolumeIndexYZ);
				}
			}
		}
	}
}


/**
 * Method will use data stored in the buffer and update temporal
 * storages appropriately.
 */
extern "C" __global__
void processBufferKernel(
		float* tempVolumeGPU, float *tempWeightsGPU,
		RecFourierProjectionTraverseSpace* spaces, int noOfSpaces,
		const float* FFTs,
		int fftSizeX, int fftSizeY,
		float* devBlobTableSqrt,
		int cMaxVolumeIndexX, int cMaxVolumeIndexYZ) {
	for (int i = 0; i < noOfSpaces; i++) {
		RecFourierProjectionTraverseSpace* space = &spaces[i];
		processProjection(
			tempVolumeGPU, tempWeightsGPU,
			FFTs, fftSizeX, fftSizeY,
			space,
			devBlobTableSqrt,
			cMaxVolumeIndexX, cMaxVolumeIndexYZ);
		__syncthreads(); // sync threads to avoid write after read problems
	}
}
