#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime_api.h>
#include "cuda_gpu_reconstruct_fourier.h"

#define BLOCK_DIM 16
#define USE_SHARED_MEM 0

#if USE_SHARED_MEM
#define BLOB_TABLE_SIZE_SQRT 10000 // keep consistent with reconstruct_fourier_gpu.h
__shared__ float BLOB_TABLE[BLOB_TABLE_SIZE_SQRT];
#endif

//#define SHARED_MEM_SIZE (int) (2.449489743*(BLOCK_DIM + 2*1.9) + 1) // FIXME block size must be settable. +1 is 'ceil' .. (block+2*blob)*sqrt(2)*sqrt(3)
//#define SHARED_MEM_SIZE_SQR (SHARED_MEM_SIZE*SHARED_MEM_SIZE)

__shared__ Point3D<float> SHARED_AABB[2];
extern __shared__ float2 INPUT_IMG[];

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort =
		true) {
	if (code != cudaSuccess) {
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
				line);
		if (abort)
			exit(code);
	}
}


__device__
void printAABB(Point3D<float>* AABB) {
	// one base
	printf("\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n",
		AABB[0].x, AABB[0].y, AABB[0].z,
		AABB[1].x, AABB[0].y, AABB[0].z,
		AABB[1].x, AABB[1].y, AABB[0].z,
		AABB[0].x, AABB[1].y, AABB[0].z,
		AABB[0].x, AABB[0].y, AABB[0].z);
	// other base with one connection
	printf("%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n",
		AABB[0].x, AABB[0].y, AABB[1].z,
		AABB[1].x, AABB[0].y, AABB[1].z,
		AABB[1].x, AABB[1].y, AABB[1].z,
		AABB[0].x, AABB[1].y, AABB[1].z,
		AABB[0].x, AABB[0].y, AABB[1].z);
	// lines between bases
	printf("%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n\n",
		AABB[1].x, AABB[0].y, AABB[1].z,
		AABB[1].x, AABB[0].y, AABB[0].z,
		AABB[1].x, AABB[1].y, AABB[0].z,
		AABB[1].x, AABB[1].y, AABB[1].z,
		AABB[0].x, AABB[1].y, AABB[1].z,
		AABB[0].x, AABB[1].y, AABB[0].z);
}

float* allocateGPU(float*& where, int size) {
	std::cout << "allocating " << size << "x" << sizeof(float) << std::endl;
	cudaMalloc((void**)&where, size * size * size * sizeof(float));
	gpuErrchk( cudaPeekAtLastError() );
	std::cout << "allocated " << where << std::endl;
	cudaMemset(where, 0.f, size * size * size * sizeof(float));
	gpuErrchk( cudaPeekAtLastError() );

	return where;
}


__device__ __constant__ int cMaxVolumeIndexX = 0;
__device__ __constant__ int cMaxVolumeIndexYZ = 0;
__device__ __constant__ bool cUseFast = false;
__device__ __constant__ float cBlobRadius = 0.f;
__device__ __constant__ float cIDeltaSqrt = 0.f;


__global__ void test_update(float* test) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	test[2*i] += i;
	test[2*i+1] += i/10.f;
	printf("%d", cMaxVolumeIndexX);
}

//void processBufferGPU(float* tempVolumeGPU,
//		float* tempWeightsGPU,
//		ProjectionData* data, int N, int bufferSize,
//		MATRIX* symmetries, MATRIX* symmetriesInv, int symSize) {
//
//// First update with respect to initial values
//   test_init<<<1,27>>>(tempVolumeGPU);
//   		gpuErrchk(cudaPeekAtLastError());
//		gpuErrchk(cudaDeviceSynchronize());
//
//// New update of device variable with respect
//// to last update
//test_update<<<1,27>>>(tempVolumeGPU);
//   		gpuErrchk(cudaPeekAtLastError());
//		gpuErrchk(cudaDeviceSynchronize());
//	}



void copyTempSpaces(std::complex<float>*** tempVol, float*** tempWeights,
		float* tempVolGPU, float* tempWeightsGPU,
		int size) {
	for (int z = 0; z < size; z++) {
		for (int y = 0; y < size; y++) {
			int index = (z * size * size) + (y * size);
			cudaMemcpy(tempVol[z][y], &tempVolGPU[2 * index], 2 * size * sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(tempWeights[z][y] , &tempWeightsGPU[index], size * sizeof(float), cudaMemcpyDeviceToHost);
		}
	}
	gpuErrchk(cudaPeekAtLastError());
}

void releaseGPU(float*& where) {
	cudaFree(where);
	gpuErrchk(cudaPeekAtLastError());
	where = NULL;
}


//#include <cuda_runtime.h>
//#include <limits>
////
////// FIXME remove
//////#include <iostream>
//////#include <stdio.h>
////
////static cuFloatComplex *d_m1;
//

//
//__device__
//float* tempVolume = NULL;
//__device__
//float* tempWeights = NULL;
////
////
////
////
__device__
inline void multiply(const MATRIX& a, const MATRIX& b, MATRIX& c) {
	c[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0];
	c[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1];
	c[0][2] = a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2];

	c[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0];
	c[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1];
	c[1][2] = a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2];

	c[2][0] = a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0];
	c[2][1] = a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1];
	c[2][2] = a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2];
}
//
//
////
////
////
////
//////CUDA functions
////__global__ void
////vecAdd(cuFloatComplex *A)
////{
////	int i = blockDim.x * blockIdx.x + threadIdx.x;
////	A[i] = cuCaddf(A[i], make_cuFloatComplex (i, 0));
////    printf("%f ", cuCrealf(A[i]));
////}
////
////__global__ void test(ProjectionDataGPU* data) {
////	int i = blockDim.x * blockIdx.x + threadIdx.x;
//////	data->img[i] += i;
////	data->img[i] = cuCaddf(data->img[i], make_cuFloatComplex (i, 0));
////    printf("%f ", cuCrealf(data->img[i]));
////}
////
////#define MUL(a,b) cuCmulf(a, make_cuFloatComplex(b, 0))
//
//
//
struct ProjectionDataGPU
{
	bool skip;
	float* img;
	float* CTF;
	float* modulator;
	int xSize;
	int ySize;
	int imgIndex;
	float weight;
	float localAInv[3][3];
	float localA[3][3];
public:
	ProjectionDataGPU() {
			setDefault();
		}
	ProjectionDataGPU(const ProjectionData& data) {
		skip = data.skip;
		if (skip) {
			setDefault();
			return;
		}
		copy(*data.img, img);
		if (NULL != data.CTF) {
			copy(*data.CTF, CTF);
		} else {
			CTF = NULL;
		}
		if (NULL != data.modulator) {
			copy(*data.modulator, modulator);
		} else {
			modulator = NULL;
		}
		xSize = data.img->getXSize();
		ySize = data.img->getYSize();
		imgIndex = data.imgIndex;
		weight = data.weight;
		data.localAInv.convertTo(localAInv);
		data.localA.convertTo(localA);
	}
	/** Remove stored data and set to skip */
	void clean() {
		cudaFree(img);
		gpuErrchk(cudaPeekAtLastError());
		cudaFree(CTF);
		gpuErrchk(cudaPeekAtLastError());
		cudaFree(modulator);
		gpuErrchk(cudaPeekAtLastError());
		setDefault();
	}
	void setDefault() {
		skip = true;
		img = NULL;
		CTF = NULL;
		modulator = NULL;
		xSize = 0;
		ySize = 0;
		imgIndex = -1;
		weight = 0;
	}
private:
	template<typename T, typename U>
	void copy(const Array2D<T>& from, U& to) {
		int xSize = from.getXSize();
		int ySize = from.getYSize();
		int N = xSize * ySize;
		// flatten the input array
		T* tmp = new T[N];
		for (int y = 0; y < ySize; y++) {
			memcpy(&tmp[y * xSize],
					from.getRow(y),
					sizeof(T)*xSize);
		}
		// Allocate device pointer.
		cudaMalloc((void**) &(to), sizeof(T) * N);
		// Copy content from host to device.
		cudaMemcpy(to, tmp, sizeof(T)*N, cudaMemcpyHostToDevice);
		gpuErrchk(cudaPeekAtLastError());
		delete[] tmp;
	}
};

__global__ void test_init(float* test,
		ProjectionDataGPU* buffer, int bufferSize,
		MATRIX* devSymmetries) {
//	if (NULL == test) {
//		printf("setting\n");
//		test = input;
//	}
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (1 == i) {
		printf("kernel %d %s %f\n", buffer[0].imgIndex, buffer[0].skip ? "true" : "false", buffer[0].weight);
		printf("kernel %d %s %f\n", buffer[1].imgIndex, buffer[1].skip ? "true" : "false", buffer[1].weight);

		for (int a = 0; a < 12; a++) {
			printf("(%f,%f) ", buffer[0].img[2*a], buffer[0].img[(2*a) + 1]);
		}

		printf("vlakno 1\n");
		for(int a = 0; a < bufferSize; a++) {
			if (buffer[a].skip) {
//				printf("preskakuju\n");
				continue;
			}
			printf("%f %f %f\n%f %f %f\n%f %f %f\n", devSymmetries[0][0][0], devSymmetries[0][0][1], devSymmetries[0][0][2],
					devSymmetries[0][1][0], devSymmetries[0][1][1], devSymmetries[0][1][2],
					devSymmetries[0][2][0], devSymmetries[0][2][1], devSymmetries[0][2][2]);
//			int index = ((buffer[a].ySize/2) * buffer[a].xSize) + buffer[a].xSize/2;
//			printf("img kernel %d: (%f,%f)\n",
//					buffer[a].imgIndex,
//					buffer[a].img[index],
//					buffer[a].img[index+1]);
		}
	}
	test[2*i] += 0.f;
	test[2*i+1] += 0.f;
}

static ProjectionDataGPU* copyProjectionData(ProjectionDataGPU* hostBuffer,
		ProjectionData* data, int bufferSize) {
	for (int i = 0; i < bufferSize; i++) {
		hostBuffer[i] = *new ProjectionDataGPU(data[i]);
	}
	ProjectionDataGPU* devBuffer;
	int size = bufferSize * sizeof(ProjectionDataGPU);
	cudaMalloc((void **) &devBuffer, size);
	gpuErrchk( cudaPeekAtLastError() );
	cudaMemcpy(devBuffer, hostBuffer, size, cudaMemcpyHostToDevice);
	gpuErrchk( cudaPeekAtLastError() );
	return devBuffer;
}
/** Returns true if x is in (min, max) interval */

template <typename T>
__device__ static bool inRange(T x, T min, T max) {
	return (x > min) && (x < max);
}

template<typename T, typename U>
__device__ inline U clamp(U val, T min, T max) {
	U res = (val > max) ? max : val;
	return (res < min) ? min : res;
}
__device__
bool getZ(float x, float y, float& z, const Point3D<float>& a, const Point3D<float>& b, const Point3D<float>& p0) {
	// from parametric eq. of the plane
	float x0 = p0.x;
	float y0 = p0.y;
	float z0 = p0.z;

	float u = ((y-y0)*a.x + (x0-x)*a.y) / (a.x * b.y - b.x * a.y);
	float t = (-x0 + x - u*b.x) / (a.x);

	z = z0 + t*a.z + u*b.z;
	return inRange(t, 0.f, 1.f) && inRange(u, 0.f, 1.f);
}

__device__
bool getY(float x, float& y, float z, const Point3D<float>& a, const Point3D<float>& b, const Point3D<float>& p0) {
	// from parametric eq. of the plane
	float x0 = p0.x;
	float y0 = p0.y;
	float z0 = p0.z;

	float u = ((z-z0)*a.x + (x0-x)*a.z) / (a.x * b.z - b.x * a.z);
	float t = (-x0 + x - u*b.x) / (a.x);

	y = y0 + t*a.y + u*b.y;
	return inRange(t, 0.f, 1.f) && inRange(u, 0.f, 1.f);
}

__device__
bool getX(float& x, float y, float z, const Point3D<float>& a, const Point3D<float>& b, const Point3D<float>& p0) {
	// from parametric eq. of the plane
	float x0 = p0.x;
	float y0 = p0.y;
	float z0 = p0.z;

	float u = ((z-z0)*a.y + (y0-y)*a.z) / (a.y * b.z - b.y * a.z);
	float t = (-y0 + y - u*b.y) / (a.y);

	x = x0 + t*a.x + u*b.x;
	return inRange(t, 0.f, 1.f) && inRange(u, 0.f, 1.f);
}

//__device__
//void createProjectionCuboid(Point3D* cuboid, float sizeX, float sizeY, float blobSize)
//{
//	float halfY = sizeY / 2.0f;
//	cuboid[0].x = cuboid[3].x = cuboid[4].x = cuboid[7].x = 0.f - blobSize;
//	cuboid[1].x = cuboid[2].x = cuboid[5].x = cuboid[6].x = sizeX + blobSize;
//
//	cuboid[0].y = cuboid[1].y = cuboid[4].y = cuboid[5].y = -(halfY + blobSize);
//	cuboid[2].y = cuboid[3].y = cuboid[6].y = cuboid[7].y = halfY + blobSize;
//
//	cuboid[0].z = cuboid[1].z = cuboid[2].z = cuboid[3].z = 0.f + blobSize;
//	cuboid[4].z = cuboid[5].z = cuboid[6].z = cuboid[7].z = 0.f - blobSize;
//}
//
__device__
void multiply(const float transform[3][3], Point3D<float>& inOut) {
	float tmp0 = transform[0][0] * inOut.x + transform[0][1] * inOut.y + transform[0][2] * inOut.z;
	float tmp1 = transform[1][0] * inOut.x + transform[1][1] * inOut.y + transform[1][2] * inOut.z;
	float tmp2 = transform[2][0] * inOut.x + transform[2][1] * inOut.y + transform[2][2] * inOut.z;
	inOut.x = tmp0;
	inOut.y = tmp1;
	inOut.z = tmp2;
}
//
//__device__
////* Apply rotation transform to cuboid */
//void rotateCuboid(Point3D* cuboid, const float transform[3][3]) {
//	for (int i = 0; i < 8; i++) {
//		multiply(transform, cuboid[i]);
//	}
//}
//
//__device__
//void translateCuboid(Point3D* cuboid, Point3D vector) {
//	for (int i = 0; i < 8; i++) {
//		cuboid[i].x += vector.x;
//		cuboid[i].y += vector.y;
//		cuboid[i].z += vector.z;
//	}
//}
//__device__
//void getVectors(const Point3D* plane, Point3D& u, Point3D& v) {
//	float x0 = plane[0].x;
//	float y0 = plane[0].y;
//	float z0 = plane[0].z;
//	u.x = plane[1].x - x0;
//	u.y = plane[1].y - y0;
//	u.z = plane[1].z - z0;
//	v.x = plane[3].x - x0;
//	v.y = plane[3].y - y0;
//	v.z = plane[3].z - z0;
//}
__device__
void computeAABB(Point3D<float>* AABB, Point3D<float>* cuboid) {
	AABB[0].x = AABB[0].y = AABB[0].z = INFINITY;
	AABB[1].x = AABB[1].y = AABB[1].z = -INFINITY;
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
	AABB[0].x = ceilf(AABB[0].x);
	AABB[0].y = ceilf(AABB[0].y);
	AABB[0].z = ceilf(AABB[0].z);

	AABB[1].x = floorf(AABB[1].x);
	AABB[1].y = floorf(AABB[1].y);
	AABB[1].z = floorf(AABB[1].z);
}


__device__
void processVoxel(float* tempVolumeGPU, float *tempWeightsGPU,
		int x, int y, int z, const float transform[3][3], float maxDistanceSqr,
		const ProjectionDataGPU* data) {
	Point3D<float> imgPos;
	float wBlob = 1.f;
	float wCTF = 1.f;
	float wModulator = 1.f;
	const float* __restrict__ img = data->img;
	const float* __restrict__ CTF = data->CTF;
	const float* __restrict__ modulator = data->modulator;
	float dataWeight = data->weight;
	int xSize = data->xSize;
	int ySize = data->ySize;

	// transform current point to center
	imgPos.x = x - cMaxVolumeIndexX/2;
	imgPos.y = y - cMaxVolumeIndexYZ/2;
	imgPos.z = z - cMaxVolumeIndexYZ/2;
	if (imgPos.x*imgPos.x + imgPos.y*imgPos.y + imgPos.z*imgPos.z > maxDistanceSqr) {
		return; // discard iterations that would access pixel with too high frequency
	}
	// rotate around center
	multiply(transform, imgPos);
	// transform back and round
	// just Y coordinate needs adjusting, since X now matches to picture and Z is irrelevant
	int imgX = clamp((int)(imgPos.x + 0.5f), 0, xSize - 1);
	int imgY = clamp((int)(imgPos.y + 0.5f + cMaxVolumeIndexYZ / 2), 0, ySize - 1);

	int index3D = z * (cMaxVolumeIndexYZ+1) * (cMaxVolumeIndexX+1) + y * (cMaxVolumeIndexX+1) + x;
	int index2D = imgY * xSize + imgX;

	if (0 != CTF) {
		wCTF = CTF[index2D];
		wModulator = modulator[index2D];
	}

	float weight = wBlob * wModulator * dataWeight;

//	if (data->imgIndex == 12 || data->imgIndex == 56) {
//		printf("%d %d %d -> %d %d (%f, %f)\n",
//				x, y, z,
//				imgX, imgY,
//				data->img[2*index2D], data->img[2*index2D+1]);
//	}


	tempVolumeGPU[2*index3D] += img[2*index2D] * weight * wCTF;
	tempVolumeGPU[2*index3D + 1] += img[2*index2D + 1] * weight * wCTF;
	tempWeightsGPU[index3D] += weight;
}

__device__
void processVoxelBlob(
		float* tempVolumeGPU, float *tempWeightsGPU,
		float* blobTableSqrt,
		int x, int y, int z, const float transform[3][3], float maxDistanceSqr,
		const ProjectionDataGPU* data, bool print, int UUID,int imgCacheDim) {
	Point3D<float> imgPos;
	int xSize = data->xSize;
	int ySize = data->ySize;
	// transform current point to center
	imgPos.x = x - cMaxVolumeIndexX/2;
	imgPos.y = y - cMaxVolumeIndexYZ/2;
	imgPos.z = z - cMaxVolumeIndexYZ/2;
	if ((imgPos.x*imgPos.x + imgPos.y*imgPos.y + imgPos.z*imgPos.z) > maxDistanceSqr) {
		return; // discard iterations that would access pixel with too high frequency
	}
	// rotate around center
	multiply(transform, imgPos);
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
	maxX = fminf(maxX, xSize-1);
	maxY = fminf(maxY, ySize-1);

	int index3D = z * (cMaxVolumeIndexYZ+1) * (cMaxVolumeIndexX+1) + y * (cMaxVolumeIndexX+1) + x;
	float volReal, volImag, w;
	volReal = volImag = w = 0.f;
	const float* __restrict__ img = data->img;
	const float* __restrict__ CTF = data->CTF;
	const float* __restrict__ modulator = data->modulator;
	float dataWeight = data->weight;

	// ugly spaghetti code, but improves performance by app. 10%
	if (0 != CTF) {
		// check which pixel in the vicinity that should contribute
		for (int i = minY; i <= maxY; i++) {
			float ySqr = (imgPos.y - i) * (imgPos.y - i);
			float yzSqr = ySqr + zSqr;
			if (yzSqr > radiusSqr) continue;
			for (int j = minX; j <= maxX; j++) {
				float xD = imgPos.x - j;
				float distanceSqr = xD*xD + yzSqr;
				if (distanceSqr > radiusSqr) continue;

				int index2D = i * xSize + j;

				float wCTF = CTF[index2D];
				float wModulator = modulator[index2D];
				int aux = (int) ((distanceSqr * cIDeltaSqrt + 0.5f)); //Same as ROUND but avoid comparison
#if USE_SHARED_MEM
				float wBlob = BLOB_TABLE[aux];
#else
				float wBlob = blobTableSqrt[aux];
#endif
				float weight = wBlob * wModulator * dataWeight;
				w += weight;
				volReal += img[2*index2D] * weight * wCTF;
				volImag += img[2*index2D + 1] * weight * wCTF;
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

				int index2D = (i - SHARED_AABB[0].y) * imgCacheDim + (j-SHARED_AABB[0].x);

				int aux = (int) ((distanceSqr * cIDeltaSqrt + 0.5f)); //Same as ROUND but avoid comparison
#if USE_SHARED_MEM
				float wBlob = BLOB_TABLE[aux];
#else
				float wBlob = blobTableSqrt[aux];
#endif

				float weight = wBlob * dataWeight;
				w += weight;
				volReal += INPUT_IMG[index2D].x * weight;
				volImag += INPUT_IMG[index2D].y * weight;
			}
		}
	}
	tempVolumeGPU[2*index3D] += volReal;
	tempVolumeGPU[2*index3D + 1] += volImag;
	tempWeightsGPU[index3D] += w;

}


//__device__
//void print(Point3D* cuboid) {
//	printf("\n");
//	for (int i = 0; i < 9; i++) {
//		printf("%f %f %f\n", cuboid[i%8].x, cuboid[i%8].y, cuboid[i%8].z);
//	}
//	printf("\n");
//}

//void AABBToCube(Point3D<float>& AABB[], Point3D<float>& cube[]) {
//	cube[0].x = AABB[0].x;
//	cube[0].y = AABB[0].y;
//	cube[0].z = AABB[0].z;
//
//	cube[1].x = AABB[0].x;
//	cube[1].y = AABB[0].y;
//	cube[1].z = AABB[0].z;
//
//	cube[2].x = AABB[0].x;
//	cube[2].y = AABB[0].y;
//	cube[2].z = AABB[0].z;
//
//	cube[3].x = AABB[0].x;
//	cube[3].y = AABB[0].y;
//	cube[3].z = AABB[0].z;
//
//	cube[4].x = AABB[0].x;
//	cube[4].y = AABB[0].y;
//	cube[4].z = AABB[0].z;
//
//	cube[5].x = AABB[0].x;
//	cube[5].y = AABB[0].y;
//	cube[5].z = AABB[0].z;
//
//	cube[6].x = AABB[0].x;
//	cube[6].y = AABB[0].y;
//	cube[6].z = AABB[0].z;
//
//	cube[7].x = AABB[0].x;
//	cube[7].y = AABB[0].y;
//	cube[7].z = AABB[0].z;
//}

__device__
void mapToImage(Point3D<float>* box, int xSize, int ySize, const float transform[3][3]) {
	for (int i = 0; i < 8; i++) {
		Point3D<float> imgPos;
		// transform current point to center
		imgPos.x = box[i].x - cMaxVolumeIndexX/2;
		imgPos.y = box[i].y - cMaxVolumeIndexYZ/2;
		imgPos.z = box[i].z - cMaxVolumeIndexYZ/2;
//		if ((imgPos.x*imgPos.x + imgPos.y*imgPos.y + imgPos.z*imgPos.z) > maxDistanceSqr) {
//			return; // discard iterations that would access pixel with too high frequency
//		}
		// rotate around center
		multiply(transform, imgPos);
		// transform back just Y coordinate, since X now matches to picture and Z is irrelevant
		imgPos.y += cMaxVolumeIndexYZ / 2;

		box[i] = imgPos;
	}
}

__device__
void calculateAABB(const TraverseSpace& tSpace, ProjectionDataGPU* projectionData) {
//	printf("img size: %d %d\n", projectionData->xSize, projectionData->ySize);
	Point3D<float> box[8];
	if (tSpace.XY == tSpace.dir) { // iterate XY plane
		box[0].x = box[3].x = box[4].x = box[7].x = blockIdx.x*blockDim.x - cBlobRadius;
		box[1].x = box[2].x = box[5].x = box[6].x = (blockIdx.x+1)*blockDim.x + cBlobRadius - 1.f;

		box[2].y = box[3].y = box[6].y = box[7].y = (blockIdx.y+1)*blockDim.y + cBlobRadius - 1.f;
		box[0].y = box[1].y = box[4].y = box[5].y = blockIdx.y*blockDim.y- cBlobRadius;

		getZ(box[0].x, box[0].y, box[0].z, tSpace.u, tSpace.v, tSpace.bottomOrigin); // lower plane
		getZ(box[4].x, box[4].y, box[4].z, tSpace.u, tSpace.v, tSpace.topOrigin); // upper plane

		getZ(box[3].x, box[3].y, box[3].z, tSpace.u, tSpace.v, tSpace.bottomOrigin); // lower plane
		getZ(box[7].x, box[7].y, box[7].z, tSpace.u, tSpace.v, tSpace.topOrigin); // upper plane

		getZ(box[2].x, box[2].y, box[2].z, tSpace.u, tSpace.v, tSpace.bottomOrigin); // lower plane
		getZ(box[6].x, box[6].y, box[6].z, tSpace.u, tSpace.v, tSpace.topOrigin); // upper plane

		getZ(box[1].x, box[1].y, box[1].z, tSpace.u, tSpace.v, tSpace.bottomOrigin); // lower plane
		getZ(box[5].x, box[5].y, box[5].z, tSpace.u, tSpace.v, tSpace.topOrigin); // upper plane

		mapToImage(box, projectionData->xSize, projectionData->ySize, tSpace.transformInv);

		computeAABB(SHARED_AABB, box);

//		if (blockIdx.x == 2 && blockIdx.y == 6) {
//			printf("blocksize: %f %f\n", AABB[1].x - AABB[0].x, AABB[1].y - AABB[0].y);
//		printf("blockid: %d %d\n", blockIdx.x, blockIdx.y);
//			printAABB(AABB);
//		}

//		printAABB(AABB);
//		printf ("projekce %d pouziva %d\n", tSpace.projectionIndex, tSpace.UUID );


	} else if (tSpace.XZ == tSpace.dir) { // iterate XZ plane
		box[0].x = box[3].x = box[4].x = box[7].x = blockIdx.x*blockDim.x - cBlobRadius;
		box[1].x = box[2].x = box[5].x = box[6].x = (blockIdx.x+1)*blockDim.x + cBlobRadius - 1.f;

		box[2].z = box[3].z = box[6].z = box[7].z = (blockIdx.y+1)*blockDim.y + cBlobRadius - 1.f;
		box[0].z = box[1].z = box[4].z = box[5].z = blockIdx.y*blockDim.y- cBlobRadius;

		getY(box[0].x, box[0].y, box[0].z, tSpace.u, tSpace.v, tSpace.bottomOrigin); // lower plane
		getY(box[4].x, box[4].y, box[4].z, tSpace.u, tSpace.v, tSpace.topOrigin); // upper plane

		getY(box[3].x, box[3].y, box[3].z, tSpace.u, tSpace.v, tSpace.bottomOrigin); // lower plane
		getY(box[7].x, box[7].y, box[7].z, tSpace.u, tSpace.v, tSpace.topOrigin); // upper plane

		getY(box[2].x, box[2].y, box[2].z, tSpace.u, tSpace.v, tSpace.bottomOrigin); // lower plane
		getY(box[6].x, box[6].y, box[6].z, tSpace.u, tSpace.v, tSpace.topOrigin); // upper plane

		getY(box[1].x, box[1].y, box[1].z, tSpace.u, tSpace.v, tSpace.bottomOrigin); // lower plane
		getY(box[5].x, box[5].y, box[5].z, tSpace.u, tSpace.v, tSpace.topOrigin); // upper plane

		mapToImage(box, projectionData->xSize, projectionData->ySize, tSpace.transformInv);

		computeAABB(SHARED_AABB, box);
	} else { // iterate YZ plane
		box[0].y = box[3].y = box[4].y = box[7].y = blockIdx.x*blockDim.x - cBlobRadius;
		box[1].y = box[2].y = box[5].y = box[6].y = (blockIdx.x+1)*blockDim.x + cBlobRadius - 1.f;

		box[2].z = box[3].z = box[6].z = box[7].z = (blockIdx.y+1)*blockDim.y + cBlobRadius - 1.f;
		box[0].z = box[1].z = box[4].z = box[5].z = blockIdx.y*blockDim.y- cBlobRadius;

		getX(box[0].x, box[0].y, box[0].z, tSpace.u, tSpace.v, tSpace.bottomOrigin); // lower plane
		getX(box[4].x, box[4].y, box[4].z, tSpace.u, tSpace.v, tSpace.topOrigin); // upper plane

		getX(box[3].x, box[3].y, box[3].z, tSpace.u, tSpace.v, tSpace.bottomOrigin); // lower plane
		getX(box[7].x, box[7].y, box[7].z, tSpace.u, tSpace.v, tSpace.topOrigin); // upper plane

		getX(box[2].x, box[2].y, box[2].z, tSpace.u, tSpace.v, tSpace.bottomOrigin); // lower plane
		getX(box[6].x, box[6].y, box[6].z, tSpace.u, tSpace.v, tSpace.topOrigin); // upper plane

		getX(box[1].x, box[1].y, box[1].z, tSpace.u, tSpace.v, tSpace.bottomOrigin); // lower plane
		getX(box[5].x, box[5].y, box[5].z, tSpace.u, tSpace.v, tSpace.topOrigin); // upper plane

//		if (tSpace.UUID == 1058) {
//			computeAABB(SHARED_AABB, box);
////			printAABB(SHARED_AABB);
//			return;
//		}

		mapToImage(box, projectionData->xSize, projectionData->ySize, tSpace.transformInv);

		computeAABB(SHARED_AABB, box);

		//		if (blockIdx.x == 2 && blockIdx.y == 6) {
		//		printf("blockid: %d %d\n", blockIdx.x, blockIdx.y);
		//			printAABB(AABB);
		//		}

		//		printAABB(AABB);
		//		printf ("projekce %d pouziva %d\n", tSpace.projectionIndex, tSpace.UUID );

	}

//	if ((int)SHARED_AABB[0].x == 101
//			&& (int)SHARED_AABB[0].y == 106
//			&& (int)SHARED_AABB[0].z == -2) {
//		printf("mam te: %d %d blok %d %d space %d\n", threadIdx.x,threadIdx.y, blockIdx.x, blockIdx.y, tSpace.UUID);
//	}
}

__device__
void processProjection(
		float* tempVolumeGPU, float *tempWeightsGPU,
	ProjectionDataGPU* projectionData,
	const TraverseSpace& tSpace,
	const float transformInv[3][3],
	float* devBlobTableSqrt,
	int imgCacheDim)
{
//	int imgSizeX = projectionData->xSize;
//	int imgSizeY = projectionData->ySize;
//	const float maxDistanceSqr = (imgSizeX+(cUseFast ? 0.f : cBlobRadius)) * (imgSizeX+(cUseFast ? 0.f : cBlobRadius));
//	Point3D origin = {cMaxVolumeIndexX/2.f, cMaxVolumeIndexYZ/2.f, cMaxVolumeIndexYZ/2.f};
//	Point3D u, v;
//	Point3D AABB[2];
//	Point3D cuboid[8];

	// calculate affected space
//	createProjectionCuboid(cuboid, imgSizeX, imgSizeY, cUseFast ? 0.f : cBlobRadius);
//	rotateCuboid(cuboid, transform);
//	translateCuboid(cuboid, origin);
//	computeAABB(AABB, cuboid, 0, 0, 0, cMaxVolumeIndexX, cMaxVolumeIndexYZ, cMaxVolumeIndexYZ);

//	getVectors(cuboid, u, v);

	// prepare traversing
//	int minY, minX, minZ;
//	int maxY, maxX, maxZ;
//	minZ = floorf(AABB[0].z);
//	minY = floorf(AABB[0].y);
//	minX = floorf(AABB[0].x);
//	maxZ = ceilf(AABB[1].z);
//	maxY = ceilf(AABB[1].y);
//	maxX = ceilf(AABB[1].x);
//	int dX, dY, dZ;
//	dX = maxX - minX;
//	dY = maxY - minY;
//	dZ = maxZ - minZ;





	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int idy = blockIdx.y*blockDim.y + threadIdx.y;
//	printf("%dx%d zpracovava %d\n", idx, idy, projectionData->imgIndex);

	if (tSpace.XY == tSpace.dir) { // iterate XY plane
		if (idy >= tSpace.minY && idy <= tSpace.maxY) {
			if (idx >= tSpace.minX && idx <= tSpace.maxX) {
				if (cUseFast) {
					float hitZ;
					if (getZ(idx, idy, hitZ, tSpace.u, tSpace.v, tSpace.bottomOrigin)) {
						int z = (int)(hitZ + 0.5f); // rounding
						processVoxel(tempVolumeGPU, tempWeightsGPU, idx, idy, z, transformInv, tSpace.maxDistanceSqr, projectionData);
					}
				} else {

					bool print = (tSpace.UUID == 1295);// && (blockIdx.x == 2) && (blockIdx.y == 6) && (threadIdx.y == 0 || threadIdx.y == 15);

					float z1, z2;
					bool hit1 = getZ(idx, idy, z1, tSpace.u, tSpace.v, tSpace.bottomOrigin); // lower plane
					bool hit2 = getZ(idx, idy, z2, tSpace.u, tSpace.v, tSpace.topOrigin); // upper plane
					if (hit1 || hit2) {
						z1 = clamp(z1, 0, cMaxVolumeIndexYZ);
						z2 = clamp(z2, 0, cMaxVolumeIndexYZ);
						int lower = floorf(fminf(z1, z2));
						int upper = ceilf(fmaxf(z1, z2));
						for (int z = lower; z <= upper; z++) {
							processVoxelBlob(tempVolumeGPU, tempWeightsGPU, devBlobTableSqrt, idx, idy, z, transformInv, tSpace.maxDistanceSqr, projectionData, print, tSpace.UUID, imgCacheDim );
						}
					}
				}
			}
		}
	} else if (tSpace.XZ == tSpace.dir) { // iterate XZ plane
		if (idy >= tSpace.minZ && idy <= tSpace.maxZ) { // map z -> y
			if (idx >= tSpace.minX && idx <= tSpace.maxX) {
				if (cUseFast) {
					float hitY;
					if (getY(idx, hitY, idy, tSpace.u, tSpace.v, tSpace.bottomOrigin)) {
						int y = (int)(hitY + 0.5f); // rounding
						processVoxel(tempVolumeGPU, tempWeightsGPU, idx, y, idy, transformInv, tSpace.maxDistanceSqr, projectionData);
					}
				} else {
					float y1, y2;
					bool hit1 = getY(idx, y1, idy, tSpace.u, tSpace.v, tSpace.bottomOrigin); // lower plane
					bool hit2 = getY(idx, y2, idy, tSpace.u, tSpace.v, tSpace.topOrigin); // upper plane
					if (hit1 || hit2) {
						y1 = clamp(y1, 0, cMaxVolumeIndexYZ);
						y2 = clamp(y2, 0, cMaxVolumeIndexYZ);
						int lower = floorf(fminf(y1, y2));
						int upper = ceilf(fmaxf(y1, y2));
						for (int y = lower; y <= upper; y++) {
							processVoxelBlob(tempVolumeGPU, tempWeightsGPU, devBlobTableSqrt, idx, y, idy, transformInv, tSpace.maxDistanceSqr, projectionData, false, tSpace.UUID, imgCacheDim);
						}
					}
				}
			}
		}
	} else { // iterate YZ plane
		if (idy >= tSpace.minZ && idy <= tSpace.maxZ) { // map z -> y
			if (idx >= tSpace.minY && idx <= tSpace.maxY) { // map y > x
				if (cUseFast) {
					float hitX;
					if (getX(hitX, idx, idy, tSpace.u, tSpace.v, tSpace.bottomOrigin)) {
						int x = (int)(hitX + 0.5f); // rounding
						processVoxel(tempVolumeGPU, tempWeightsGPU, x, idx, idy, transformInv, tSpace.maxDistanceSqr, projectionData);
					}
				} else {
					float x1, x2;
					bool hit1 = getX(x1, idx, idy, tSpace.u, tSpace.v, tSpace.bottomOrigin); // lower plane
					bool hit2 = getX(x2, idx, idy, tSpace.u, tSpace.v, tSpace.topOrigin); // upper plane
					if (hit1 || hit2) {
						x1 = clamp(x1, 0, cMaxVolumeIndexX);
						x2 = clamp(x2, 0, cMaxVolumeIndexX);
						int lower = floorf(fminf(x1, x2));
						int upper = ceilf(fmaxf(x1, x2));
						for (int x = lower; x <= upper; x++) {
							processVoxelBlob(tempVolumeGPU, tempWeightsGPU, devBlobTableSqrt, x, idx, idy, transformInv, tSpace.maxDistanceSqr, projectionData, false,tSpace.UUID, imgCacheDim);
						}
					}
				}
			}
		}
	}

//	if (projectionData->imgIndex == 12 || projectionData->imgIndex == 56) {
//		printf("\n----------------\nimg %d:", projectionData->imgIndex);
//		print(cuboid);
//		printf("------------------\n");
//		printAABB(AABB);
//	}

}

__device__
bool blockHasWork(int imgXSize, int imgYSize) {
	bool hitX = (SHARED_AABB[0].x < imgXSize) && (SHARED_AABB[1].x >= 0);
	bool hitY = (SHARED_AABB[0].y < imgYSize) && (SHARED_AABB[1].y >= 0);
	return hitX && hitY;
}

__device__
void getImgData(int tXindex, int tYindex, ProjectionDataGPU* data,float& vReal, float& vImag) {
	int imgXindex = tXindex + SHARED_AABB[0].x;
	int imgYindex = tYindex + SHARED_AABB[0].y;
	if ((imgXindex >=0)
			&& (imgXindex < data->xSize)
			&& (imgYindex >=0)
			&& (imgYindex < data->ySize))
	{
		int index = imgYindex * data->xSize + imgXindex;
		vReal = data->img[2*index];
		vImag = data->img[2*index + 1];

	} else {
		vReal = vImag = 0.f;
	}
}


__device__
void copyData(ProjectionDataGPU* data, int imgCacheDim) {
	float xAABBSize = SHARED_AABB[1].x - SHARED_AABB[0].x + 1;
	float yAABBSize = SHARED_AABB[1].y - SHARED_AABB[0].y + 1;
	int xLoops = ceilf(xAABBSize / blockDim.x);
	int yLoops = ceilf(yAABBSize / blockDim.y);
//					if (threadIdx.x == 0 && threadIdx.y == 0 && blockIdx.x == 3 && blockIdx.y == 1 && UUID == 16751) {
//						printAABB(SHARED_AABB);
//					}

	for (int y = threadIdx.y; y < imgCacheDim; y += blockDim.y) {
//		int tYindex = threadIdx.y + i*blockDim.y;
//		if (tYindex < SHARED_MEM_SIZE) {
			for (int x = threadIdx.x; x < imgCacheDim; x += blockDim.x) {
//				int tXindex = threadIdx.x + j*blockDim.x;
//				if (tXindex < SHARED_MEM_SIZE) {
					int memIndex = y * imgCacheDim + x;
					float vReal, vImag;
					getImgData(x, y, data, vReal, vImag);


					INPUT_IMG[memIndex].x = vReal;
					INPUT_IMG[memIndex].y = vImag;
//				}
			}
//		}
	}
//
//
//	int id = threadIdx.y*blockDim.x + threadIdx.x;
//	int blockSize = blockDim.x * blockDim.y;
//	for (int i = id; i < BLOB_TABLE_SIZE_SQRT; i+= blockSize) {
//		BLOB_TABLE[i] = devBlobTableSqrt[i];
//	}
}

__global__
void processBufferKernel(
		float* tempVolumeGPU, float *tempWeightsGPU,
		ProjectionDataGPU* buffer, int bufferSize,
		TraverseSpace* traverseSpaces, int noOfTransforms,
		float* devBlobTableSqrt,
		int imgCacheDim) {
#if USE_SHARED_MEM
	if ( ! cUseFast) {
		int id = threadIdx.y*blockDim.x + threadIdx.x;
		int blockSize = blockDim.x * blockDim.y;
		for (int i = id; i < BLOB_TABLE_SIZE_SQRT; i+= blockSize)
			BLOB_TABLE[i] = devBlobTableSqrt[i];
		__syncthreads();
	}
#endif



	for (int i = 0; i < noOfTransforms; i++) {
		TraverseSpace* space = &traverseSpaces[i];
		ProjectionDataGPU* data = &buffer[space->projectionIndex];

		if (threadIdx.x == 0 && threadIdx.y == 0 && data->imgIndex == 0 && blockIdx.x == 0 && blockIdx.y == 0) {
								printf("velikost %d \n", imgCacheDim);
		//							threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);
							}



		if ( ! cUseFast) {
			if ((threadIdx.x == 0) && (threadIdx.y == 0)) {
				calculateAABB(*space, data);
//					printAABB(SHARED_AABB);
//					printf("vetev %d\n", space->dir);
//				}
			}
			__syncthreads();
			if (blockHasWork(data->xSize, data->ySize)) {
//					if (space->UUID == 3916 && threadIdx.x == 0 && threadIdx.y == 0) {
//						printAABB(SHARED_AABB);
//						printf("vetev %d\n", space->dir);
//					}
				copyData(data, imgCacheDim);
//				if (blockIdx.x ==11 && blockIdx.y == 12 && space->UUID == 1058 && threadIdx.x == 0 && threadIdx.y == 0) {
//					printf("nakopiroval jsem data\n");
//					printAABB(SHARED_AABB);
//				}
				__syncthreads();
			} else {
//				__syncthreads();
				continue; // whole block can exit
			}


		}

//		if (blockIdx.x ==11 && blockIdx.y == 12 && space->UUID == 1058 && threadIdx.x == 0 && threadIdx.y == 0) {
//							printf("jdu zpracovat projekci\n");
//							printAABB(SHARED_AABB);
//						}
//		__syncthreads();
		processProjection(
			tempVolumeGPU, tempWeightsGPU,
			data, *space, space->transformInv,
			devBlobTableSqrt,
			imgCacheDim);

		__syncthreads();
	}
}

void processBufferGPU(float* tempVolumeGPU,
		float* tempWeightsGPU,
		ProjectionData* data, int N, int bufferSize,
		TraverseSpace* traverseSpaces, int noOfTransforms,
		int maxVolIndexX, int maxVolIndexYZ,
		bool useFast, float blobRadius,
		float iDeltaSqrt,
		float* blobTableSqrt, int blobTableSize) {

	ProjectionDataGPU* hostBuffer = new ProjectionDataGPU[bufferSize];
	ProjectionDataGPU* devBuffer = copyProjectionData(hostBuffer, data, bufferSize);

//	size_t currentLimit;
//	cudaDeviceGetLimit(&currentLimit, cudaLimitPrintfFifoSize);
//	cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 1048576 * 2);
//	gpuErrchk( cudaPeekAtLastError() );



//	printf("vstup %f %f %f\n%f %f %f\n%f %f %f\n",
//			data[0].localA(0, 0), data[0].localA(0,1), data[0].localA(0,2),
//			data[0].localA(1, 0), data[0].localA(1,1), data[0].localA(1,2),
//			data[0].localA(2, 0), data[0].localA(2,1), data[0].localA(2,2));
//	printf("vstup1 %f %f %f\n%f %f %f\n%f %f %f\n",
//			data[0].localAInv(0, 0), data[0].localAInv(0,1), data[0].localAInv(0,2),
//			data[0].localAInv(1, 0), data[0].localAInv(1,1), data[0].localAInv(1,2),
//			data[0].localAInv(2, 0), data[0].localAInv(2,1), data[0].localAInv(2,2));

	cudaMemcpyToSymbol(cMaxVolumeIndexX, &maxVolIndexX,sizeof(maxVolIndexX));
	gpuErrchk( cudaPeekAtLastError() );
	cudaMemcpyToSymbol(cMaxVolumeIndexYZ, &maxVolIndexYZ,sizeof(maxVolIndexYZ));
	gpuErrchk( cudaPeekAtLastError() );
	cudaMemcpyToSymbol(cUseFast, &useFast,sizeof(useFast));
	gpuErrchk( cudaPeekAtLastError() );
	cudaMemcpyToSymbol(cBlobRadius, &blobRadius,sizeof(blobRadius));
	gpuErrchk( cudaPeekAtLastError() );
	cudaMemcpyToSymbol(cIDeltaSqrt, &iDeltaSqrt,sizeof(iDeltaSqrt));
	gpuErrchk( cudaPeekAtLastError() );


//		// New update of device variable with respect
//		// to last update
//		test_update<<<1,27>>>(tempVolumeGPU);
//		   		gpuErrchk(cudaPeekAtLastError());
//				gpuErrchk(cudaDeviceSynchronize());
//				return;

	TraverseSpace* devTravSpaces;
	cudaMalloc((void **) &devTravSpaces, noOfTransforms*sizeof(TraverseSpace));
	gpuErrchk( cudaPeekAtLastError() );
	float* devBlobTableSqrt;
	cudaMalloc((void **) &devBlobTableSqrt, blobTableSize*sizeof(float));
	gpuErrchk( cudaPeekAtLastError() );

	cudaMemcpy(devTravSpaces, traverseSpaces, noOfTransforms*sizeof(TraverseSpace), cudaMemcpyHostToDevice);
	gpuErrchk( cudaPeekAtLastError() );
	cudaMemcpy(devBlobTableSqrt, blobTableSqrt, blobTableSize*sizeof(float), cudaMemcpyHostToDevice);
	gpuErrchk( cudaPeekAtLastError() );

	int size2D = maxVolIndexYZ + 1;

	dim3 dimBlock(BLOCK_DIM, BLOCK_DIM);
	dim3 dimGrid((int)ceil(size2D/dimBlock.x),(int)ceil(size2D/dimBlock.y));

	int imgCacheDim = ceil(sqrt(2.f) * sqrt(3.f) *(BLOCK_DIM + 2*blobRadius));
	printf ("imgCacheDim: %d\n", imgCacheDim);

	processBufferKernel<<<dimGrid, dimBlock, imgCacheDim*imgCacheDim*sizeof(float2)>>>(
			tempVolumeGPU, tempWeightsGPU,
			devBuffer, bufferSize,
			devTravSpaces, noOfTransforms,
			devBlobTableSqrt,
			imgCacheDim);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
//	ProjectionDataGPU*	hostBuffer = new ProjectionDataGPU[bufferSize];
//	cudaMemcpy(devBuffer, hostBuffer, bufferSize*sizeof(ProjectionDataGPU), cudaMemcpyDeviceToHost);
//	for (int i = 0; i < bufferSize; i++) {
//		hostBuffer[i].clean();
//		data[i].clean();
//	}

//	for(int a = 0; a < bufferSize; a++) {
//		Array2D<std::complex<float> >* imgt = data[a].img;
//		if (data[a].skip) continue;
//		std::cout << "img " << data[a].imgIndex <<
//				": " << (*imgt)(imgt->getXSize()/2, imgt->getYSize()/2) <<
//				std::endl;
//	}



//	// First update with respect to initial values
//		   test_init<<<1,27>>>(tempVolumeGPU, devBuffer, bufferSize, devSymmetries);
//		   		gpuErrchk(cudaPeekAtLastError());
//				gpuErrchk(cudaDeviceSynchronize());
//

	for (int i = 0 ; i < bufferSize; i++) {
		hostBuffer[i].clean();
	}
	cudaFree(devBuffer);
	gpuErrchk( cudaPeekAtLastError() );
	delete[] hostBuffer;
	cudaFree(devTravSpaces);
	gpuErrchk( cudaPeekAtLastError() );
	cudaFree(devBlobTableSqrt);
	gpuErrchk( cudaPeekAtLastError() );
}
//
//void getTempSpaces(int size, std::complex<float>***& volume, float***& tempWeights) {
//	std::cout << size << std::endl;
//	std::complex<float>* tmp1 = new std::complex<float>[20];
//	std::complex<float>* tmp2;
//	cudaGetSymbolAddress((void **)&tmp2, tempVolume);
//	gpuErrchk( cudaPeekAtLastError() );
////	std::cout << cuCrealf(tempVolume[0]) << std::endl;
////	cudaMemcpy(tmp1, tmp2, 20*sizeof(std::complex<float>), cudaMemcpyDeviceToHost);
//	cudaMemcpyFromSymbol(tmp1, tempVolume, 20*sizeof(std::complex<float>));
//	gpuErrchk( cudaPeekAtLastError() );
//	for(int z = 0; z < 20; z++) {
//	std::cout << tmp1[z] << " " << std::endl;
////		for(int y = 0; y < 20; y++) {
////			for(int x = 0; x < 20; x++) {
////				volume[z][y][x] = tmp1[z*20*20+ y * size+ x];
////			}
////		}
//	}
//	delete[] tmp1;
//
//
////	float* tmp2 = new float[20];
////	cudaMemcpy(tmp2, tempWeights, 20*sizeof(float), cudaMemcpyDeviceToHost);
////	gpuErrchk( cudaPeekAtLastError() );
////	for(int z = 0; z < size; z++) {
////		for(int y = 0; y < size; y++) {
////			for(int x = 0; x < size; x++) {
////				tempWeights[z][y][x] = tmp2[z*size*size+ y * size+ x];
////			}
////		}
////	}
////	delete[] tmp2;
//}
//
//void allocateTempSpaces(int size) {
//	float* tmp1;
//	cudaGetSymbolAddress((void **)&tmp1, tempVolume);
//	gpuErrchk( cudaPeekAtLastError() );
//	if (NULL == tmp1) {
//		std::cout << "allocating " << size << std::endl;
//		cudaMalloc((void **)&tmp1, 20*sizeof(std::complex<float>));
//		gpuErrchk( cudaPeekAtLastError() );
//		cudaMemset(tmp1, 0.1,20 * sizeof(std::complex<float>));
//		gpuErrchk( cudaPeekAtLastError() );
////		cudaMalloc((void **)&tempWeights, 20*sizeof(float));
////		gpuErrchk( cudaPeekAtLastError() );
//	}
//}
//
////void copyBuffer(ProjectionData* data, int size) {
//////	if (NULL == buffer) {
//////		cudaMalloc((void **)&buffer, size*sizeof(ProjectionDataGPU));
////	ProjectionDataGPU*	hostBuffer = new ProjectionDataGPU[size];
//////	}
////	for (int i = 0; i < size; i++) {
////		hostBuffer[i] = *copyProjectionData(data[i]);
////	}
////	cudaMemcpy(devBuffer, hostBuffer, size*sizeof(ProjectionDataGPU), cudaMemcpyHostToDevice);
////}
////
////
////void a(std::complex<float> *arr) {
////	std::cerr  << "A"<< std::endl;
////
////	for(int i = 0; i < num; i++) {
////		printf("%f ", arr[i].real());
////	}
////	printf("\n");
////
//////	//CUDA code
//////	float *d_m1;
//////
////	size_t size=num*sizeof(std::complex<float>);
////	cudaMalloc((void **)&d_m1, size);
////	cudaMemcpy(d_m1, arr, size, cudaMemcpyHostToDevice);
////	std::cerr  << "A hotovo"<< std::endl;
////}
////
////void b() {
////	std::cerr  << "B"<< std::endl;
////	vecAdd<<<1, 20>>>(d_m1);
////	std::cerr  << "B hotovo"<< std::endl;
////}
////
////void c(std::complex<float>* arr) {
////	std::cerr  << "C"<< std::endl;
////	size_t size=num*sizeof(std::complex<float>);
////	cudaMemcpy(arr, d_m1, size, cudaMemcpyDeviceToHost);
////	cudaFree(d_m1);
////	for(int i = 0; i < num; i++) {
////		printf("%f ", arr[i].real());
////	}
////	printf("\n");
////	std::cerr  << "C hotovo"<< std::endl;
////}
////
////void moje(){
//////
//////	int num = 20;
//////	std::cerr  << "Inside CUDA function "<< std::endl;
//////
//////	float* array = new float[num];
//////	memset(array, 1, num*sizeof(float));
//////
//////	for(int i = 0; i < num; i++) {
//////		printf("%f ", array[i]);
//////	}
//////	printf("\n");
//////
////////	//CUDA code
//////
////////
//////	size_t size=num*sizeof(float);
//////	cudaMalloc((void **)&d_m1, size);
////////	cudaMalloc((void **)&d_m2, matSize);
////////	cudaMalloc((void **)&d_m3, matSize);
////////
////////
//////	cudaMemcpy(d_m1, array, size, cudaMemcpyHostToDevice);
////////	cudaMemcpy(d_m2, m2, matSize, cudaMemcpyHostToDevice);
////////
////////	int numTh = 1024;
////////	int numBlk = num*num/numTh;
////////	if ((num*num)%numTh >0){
////////			numBlk++;
////////	}
//////
//////	vecAdd<<<1, 20>>>(d_m1);
//////	cudaMemcpy(array, d_m1, size, cudaMemcpyDeviceToHost);
////////
//////
//////
//////	cudaFree(d_m1);
////////	cudaFree(d_m2);
////////	cudaFree(d_m3);
//////	std::cerr  << "CUDA function done"<< std::endl;
////}
