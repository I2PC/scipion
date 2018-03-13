#include "/home/david/GIT/Scipion/software/em/xmipp/libraries/data/reconstruct_fourier_projection_traverse_space.h"
#include "/home/david/GIT/Scipion/software/em/xmipp/libraries/data/point3D.h"

__device__ __constant__ float cBlobRadius = 1.900000f;
__device__ __constant__ float cIDeltaSqrt = 2769.806152f;

__device__
const float* getNthItem(const float* array, int itemIndex, int fftSizeX, int fftSizeY, bool isFFT) {
		return array + (fftSizeX * fftSizeY * itemIndex * (isFFT ? 2 : 1)); // *2 since it's complex
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

__device__
void processPixelBlob(
	float* tempVolumeGPU, float* tempWeightsGPU,
	const float* FFTs,
	int fftSizeX, int fftSizeY,
	const RecFourierProjectionTraverseSpace* const space,
	const float* blobTableSqrt,
	int cMaxVolumeIndexX, int cMaxVolumeIndexYZ,
	int x, int y)
{
	Point3D<float> imgPos;
	const float* __restrict__ img = getNthItem(FFTs, space->projectionIndex, fftSizeX, fftSizeY, true);

	float dataWeight = space->weight;
	int index2D = y * fftSizeX + x;

	if (img[2*index2D] == 0.f && img[2*index2D + 1] == 0.f) {
		return;
	}

	imgPos.x = x;
	imgPos.y = y - cMaxVolumeIndexYZ/2;
	imgPos.z = 0.f;
	// rotate around center
	multiply(space, imgPos);
	// transform to center of the volume
	float volX = imgPos.x + cMaxVolumeIndexX/2;
	float volY = imgPos.y + cMaxVolumeIndexYZ/2;
	float volZ = imgPos.z + cMaxVolumeIndexYZ/2;

	// create blob bounding box
	int minX = ceilf(volX - cBlobRadius);
	int maxX = floorf(volX + cBlobRadius);
	int minY = ceilf(volY - cBlobRadius);
	int maxY = floorf(volY + cBlobRadius);
	int minZ = ceilf(volZ - cBlobRadius);
	int maxZ = floorf(volZ + cBlobRadius);
	minX = fmaxf(minX, 0);
	minY = fmaxf(minY, 0);
	minZ = fmaxf(minZ, 0);
	maxX = fminf(maxX, cMaxVolumeIndexX);
	maxY = fminf(maxY, cMaxVolumeIndexYZ);
	maxZ = fminf(maxZ, cMaxVolumeIndexYZ);
	float radiusSqr = cBlobRadius * cBlobRadius;
	for (int l = minZ; l <= maxZ; l++) {
		float zSqr = (volZ - l) * (volZ - l);
		for (int i = minY; i <= maxY; i++) {
			float ySqr = (volY - i) * (volY - i);
			float yzSqr = ySqr + zSqr;
			if (yzSqr > radiusSqr) continue;
			for (int j = minX; j <= maxX; j++) {
				float xD = volX - j;
				float distanceSqr = xD*xD + yzSqr;
				if (distanceSqr > radiusSqr) continue;

				int index3D = l * (cMaxVolumeIndexYZ+1) * (cMaxVolumeIndexX+1) + i * (cMaxVolumeIndexX+1) + j;

				int aux = (int) ((distanceSqr * cIDeltaSqrt + 0.5f));
				float wBlob = blobTableSqrt[aux];
				float weight = wBlob * dataWeight;

				atomicAdd(&tempVolumeGPU[2*index3D], img[2*index2D] * weight );
				atomicAdd(&tempVolumeGPU[2*index3D + 1], img[2*index2D + 1] * weight );
				atomicAdd(&tempWeightsGPU[index3D], weight);
			}
		}
	}
}

/**
 * Method will use data stored in the buffer and update temporal
 * storages appropriately.
 */

extern "C" __global__
void processBufferKernelInverseReference(float* tempVolumeGPU, float* tempWeightsGPU,
		RecFourierProjectionTraverseSpace* spaces, int noOfSpaces,
		const float* FFTs,
		int fftSizeX, int fftSizeY,
		float* devBlobTableSqrt,
		int cMaxVolumeIndexX, int cMaxVolumeIndexYZ) {
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int idy = blockIdx.y*blockDim.y + threadIdx.y;
	if(idx >= fftSizeX  || idy >= fftSizeY) {
		return;
	}

	for (int i = 0; i < noOfSpaces; i++) {
			RecFourierProjectionTraverseSpace* space = &spaces[i];

			processPixelBlob(tempVolumeGPU, tempWeightsGPU,
				FFTs,
				fftSizeX, fftSizeY,
				space,
				devBlobTableSqrt,
				cMaxVolumeIndexX, cMaxVolumeIndexYZ, idx, idy);
		__syncthreads(); // sync threads to avoid write after read problems
	}
}
