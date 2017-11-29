#include "/home/david/GIT/Scipion/software/em/xmipp/libraries/data/reconstruct_fourier_projection_traverse_space.h"
#include "/home/david/GIT/Scipion/software/em/xmipp/libraries/data/point3D.h"



#if SHARED_BLOB_TABLE
__shared__ float BLOB_TABLE[BLOB_TABLE_SIZE_SQRT];
#endif

#if SHARED_IMG
__shared__ Point3D<float> SHARED_AABB[2];
extern __shared__ float2 IMG[];
#endif


// FIELDS



// Holding blob coefficient table. Present on GPU
//float* devBlobTableSqrt = NULL;

//__device__ __constant__ int cMaxVolumeIndexX = 0;
//__device__ __constant__ int cMaxVolumeIndexYZ = 0;
__device__ __constant__ float cBlobRadius = 1.900000f;
__device__ __constant__ float cBlobAlpha = 15.000000f;
__device__ __constant__ float cIw0 = 0.581191f;
__device__ __constant__ float cIDeltaSqrt = 2769.806152f;


__device__ void operator+=(float2 &a, float2 b)
{
    a.x += b.x;
    a.y += b.y;
}

__device__ void operator*=(float2 &a, float2 b)
{
    a.x *= b.x;
    a.y *= b.y;
}

__device__ float2 operator*(float2 a, float b)
{
    return make_float2(a.x * b, a.y * b);
}

__device__
float bessi0(float x) // only positive values expected
{
    float y, ax, ans;
    if (x < 3.75f)
    {
        y = x * x * 0.0711111f; // (x/3.75)^2
        ans = (y*(y*(y*(y*((0.0045813f*y + 0.0360768f)*y + 0.265973f) + 1.20675f) + 3.08994f) + 3.51562f) + 1.f);

    }
    else
    {
        y = 3.75f/x;
        ans = (expf(x)*(y*(y*(y*(y*(y*(y*((0.00392377f*y - 0.0164763f)*y + 0.0263554f) - 0.0205771f)
        		+ 0.009162f) - 0.001575f) + 0.002253f) + 0.013285f) + 0.398942f)) * rsqrtf(x);
    }
    return ans;
}


__device__ float bessi00(float x)
{
// -- See paper
// J.M. Blair, "Rational Chebyshev approximations for the modified Bessel functions I_0(x) and I_1(x)", Math. Comput., vol. 28, n. 126, pp. 581-583, Apr. 1974.

   float num, den, x2;

   x2 = x*x;

//   x=abs(x);
//
//   if (x > 15.0)
//   {
//        den = 1.0 / x;
//        num =              -4.4979236558557991E+006;
//        num = fma (num, den,  2.7472555659426521E+006);
//        num = fma (num, den, -6.4572046640793153E+005);
//        num = fma (num, den,  8.5476214845610564E+004);
//        num = fma (num, den, -7.1127665397362362E+003);
//        num = fma (num, den,  4.1710918140001479E+002);
//        num = fma (num, den, -1.3787683843558749E+001);
//        num = fma (num, den,  1.1452802345029696E+000);
//        num = fma (num, den,  2.1935487807470277E-001);
//        num = fma (num, den,  9.0727240339987830E-002);
//        num = fma (num, den,  4.4741066428061006E-002);
//        num = fma (num, den,  2.9219412078729436E-002);
//        num = fma (num, den,  2.8050629067165909E-002);
//        num = fma (num, den,  4.9867785050221047E-002);
//        num = fma (num, den,  3.9894228040143265E-001);
//        num = num * den;
//        den = sqrt (x);
//        num = num * den;
//        den = exp (0.5 * x);  /* prevent premature overflow */
//        num = num * den;
//        num = num * den;
//		return num;
//   }
//   else
//   {
	  num = -0.28840544803647313855232E-028f;
	  num = fmaf (num, x2, -0.72585406935875957424755E-025f);
	  num = fmaf (num, x2, -0.1247819710175804058844059E-021f);
	  num = fmaf (num, x2, -0.15795544211478823152992269E-018f);
	  num = fmaf (num, x2, -0.15587387207852991014838679E-015f);
	  num = fmaf (num, x2, -0.121992831543841162565677055E-012f);
	  num = fmaf (num, x2, -0.760147559624348256501094832E-010f);
	  num = fmaf (num, x2, -0.375114023744978945259642850E-007f);
	  num = fmaf (num, x2, -0.1447896113298369009581404138E-004f);
	  num = fmaf (num, x2, -0.4287350374762007105516581810E-002f);
	  num = fmaf (num, x2, -0.947449149975326604416967031E+000f);
	  num = fmaf (num, x2, -0.1503841142335444405893518061E+003f);
	  num = fmaf (num, x2, -0.1624100026427837007503320319E+005f);
	  num = fmaf (num, x2, -0.11016595146164611763171787004E+007f);
	  num = fmaf (num, x2, -0.4130296432630476829274339869E+008f);
	  num = fmaf (num, x2, -0.6768549084673824894340380223E+009f);
	  num = fmaf (num, x2, -0.27288446572737951578789523409E+010f);

	  den = 0.1E+001;
	  den = fmaf (den, x2, -0.38305191682802536272760E+004f);
	  den = fmaf (den, x2, 0.5356255851066290475987259E+007f);
	  den = fmaf (den, x2, -0.2728844657273795156746641315E+010f);

	  return num/den;
//   }
}

__device__
float bessi0_moje(float x) {
	// stable rational minimax approximations to the modified bessel functions, blair, edwards
	// from table 5
	float x2 = x*x;
	float num = -0.8436825781374849e-19f; // p11
	num = fmaf(num, x2, -0.93466495199548700e-17f); // -p4*2*22.567 (from inner bracket)
	num = fmaf(num, x2, -0.15716375332511895e-13f); // p4*22.567^2 + p2
	num = fmaf(num, x2, -0.42520971595532318e-11f); // p1
	num = fmaf(num, x2, -0.13704363824102120e-8f); //p0
	num = fmaf(num, x2, -0.28508770483148419e-6f); // -p4*2*22.567 (from inner bracket)
	num = fmaf(num, x2, -0.44322160233346062e-4f); // -p4*2*22.567 (from inner bracket)
	num = fmaf(num, x2, -0.46703811755736946e-2f); // -p4*2*22.567 (from inner bracket)
	num = fmaf(num, x2, -0.31112484643702141e-0f); // -p4*2*22.567 (from inner bracket)
	num = fmaf(num, x2, -0.11512633616429962e+2f); // -p4*2*22.567 (from inner bracket)
	num = fmaf(num, x2, -0.18720283332732112e+3f); // -p4*2*22.567 (from inner bracket)
	num = fmaf(num, x2, -0.75281108169006924e+3f); // -p4*2*22.567 (from inner bracket)


	float den = 1.f; //q1
	den = fmaf(den, x2, -0.75281109410939403e+3f); // q00

	return num/den;
}


__device__
float bessi1(float x)
{
    float ax, ans;
    float y;
    if ((ax = fabsf(x)) < 3.75f)
    {
        y = x / 3.75f;
        y *= y;
        ans = ax * (0.5f + y * (0.87890594f + y * (0.51498869f + y * (0.15084934f
                               + y * (0.2658733e-1f + y * (0.301532e-2f + y * 0.32411e-3f))))));
    }
    else
    {
        y = 3.75f / ax;
        ans = 0.2282967e-1f + y * (-0.2895312e-1f + y * (0.1787654e-1f
                                  - y * 0.420059e-2f));
        ans = 0.39894228f + y * (-0.3988024e-1f + y * (-0.362018e-2f
                                + y * (0.163801e-2f + y * (-0.1031555e-1f + y * ans))));
        ans *= (expf(ax) / sqrtf(ax));
    }
    return x < 0.0 ? -ans : ans;
}

__device__
float bessi2(float x)
{
    return (x == 0) ? 0 : bessi0(x) - ((2*1) / x) * bessi1(x);
}

__device__
float bessi3(float x)
{
    return (x == 0) ? 0 : bessi1(x) - ((2*2) / x) * bessi2(x);
}

__device__
float bessi4(float x)
{
    return (x == 0) ? 0 : bessi2(x) - ((2*3) / x) * bessi3(x);
}


__device__
float kaiserValue(float r)
{
    float rda, rdas, arg, w;

    rda = r / cBlobRadius;
    if (rda <= 1.f)
    {
        rdas = rda * rda;
        arg = cBlobAlpha * sqrtf(1.f - rdas);
        if (blobOrder == 0)
        {
            w = bessi0(arg) * 0.000002944f;
        }
        else if (blobOrder == 1)
        {
            w = sqrtf (1.f - rdas);
                w *= bessi1(arg) / bessi1(cBlobAlpha);
        }
        else if (blobOrder == 2)
        {
            w = sqrtf (1.f - rdas);
            w = w * w;
                w *= bessi2(arg) / bessi2(cBlobAlpha);
        }
        else if (blobOrder == 3)
        {
            w = sqrtf (1.f - rdas);
            w = w * w * w;
                w *= bessi3(arg) / bessi3(cBlobAlpha);
        }
        else if (blobOrder == 4)
        {
            w = sqrtf (1.f - rdas);
            w = w * w * w *w;
			w *= bessi4(arg) / bessi4(cBlobAlpha);
        }
        else {
        	printf("order (%d) out of range in kaiser_value(): %s, %d\n", blobOrder, __FILE__, __LINE__);
        }
    }
    else
        w = 0.f;

    return w;
}

__device__
float kaiserValue2(float dist) {
	float arg = cBlobAlpha * sqrtf(1.f - (dist * 0.27700831f)); // alpha * sqrt(1-(dist/blobRadius^2))
//	return bessi0(arg) * 0.000002944f;
//	 return cyl_bessel_i0f (arg ) *  0.000002944f;
	return bessi0_moje(arg) * 0.000002944f;
}




/** Returns true if x is in (min, max), i.e. opened, interval */
template <typename T>
__device__
static bool inRange(T x, T min, T max) {
	return (x > min) && (x < max);
}


/** Returns value within the range (included) */
template<typename T, typename U>
__device__
static U clamp(U val, T min, T max) {
	U res = (val > max) ? max : val;
	return (res < min) ? min : res;
}

__device__
const float* getNthItem(const float* array, int itemIndex, int fftSizeX, int fftSizeY, bool isFFT) {
//	if (array == FFTs)
		return array + (fftSizeX * fftSizeY * itemIndex * (isFFT ? 2 : 1)); // *2 since it's complex
//	if (array == CTFs) return array + (fftSizeX * fftSizeY * itemIndex);
//	if (array == modulators) return array + (fftSizeX * fftSizeY * itemIndex);
//	if (array == paddedImages) return array + (paddedImgSize * paddedImgSize * itemIndex);
//	return NULL; // undefined
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
 * Method will rotate box using transformation matrix around center of the
 * working space
 */
__device__
void rotate(Point3D<float>* box, const RecFourierProjectionTraverseSpace* tSpace) {
	for (int i = 0; i < 8; i++) {
		Point3D<float> imgPos;
		// transform current point to center
		imgPos.x = box[i].x - cMaxVolumeIndexX/2;
		imgPos.y = box[i].y - cMaxVolumeIndexYZ/2;
		imgPos.z = box[i].z - cMaxVolumeIndexYZ/2;
		// rotate around center
		multiply(tSpace, imgPos);
		// transform back just Y coordinate, since X now matches to picture and Z is irrelevant
		imgPos.y += cMaxVolumeIndexYZ / 2;

		box[i] = imgPos;
	}
}

/** Compute Axis Aligned Bounding Box of given cuboid */
__device__
void computeAABB(Point3D<float>* AABB, Point3D<float>* cuboid) {
	AABB[0].x = AABB[0].y = AABB[0].z = 0x7f800000;
	AABB[1].x = AABB[1].y = AABB[1].z = 0xff800000;
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

/**
 * Method will map one voxel from the temporal
 * spaces to the given projection and update temporal spaces
 * using the pixel value of the projection.
 */
__device__
void processVoxel(
	float2* tempVolumeGPU, float* tempWeightsGPU,
	const float* FFTs, const float* CTFs, const float* modulators,
	int fftSizeX, int fftSizeY,
	int x, int y, int z,
	const RecFourierProjectionTraverseSpace* const space)
{
	Point3D<float> imgPos;
	float wBlob = 1.f;
	float wCTF = 1.f;
	float wModulator = 1.f;
	const float* __restrict__ img = getNthItem(FFTs, space->projectionIndex, fftSizeX, fftSizeY, true);


	float dataWeight = space->weight;

	// transform current point to center
	imgPos.x = x - cMaxVolumeIndexX/2;
	imgPos.y = y - cMaxVolumeIndexYZ/2;
	imgPos.z = z - cMaxVolumeIndexYZ/2;
	if (imgPos.x*imgPos.x + imgPos.y*imgPos.y + imgPos.z*imgPos.z > space->maxDistanceSqr) {
		return; // discard iterations that would access pixel with too high frequency
	}
	// rotate around center
	multiply(space, imgPos);
	if (imgPos.x < 0.f) return; // reading outside of the image boundary. Z is always correct and Y is checked by the condition above
	// transform back and round
	// just Y coordinate needs adjusting, since X now matches to picture and Z is irrelevant
	int imgX = clamp((int)(imgPos.x + 0.5f), 0, fftSizeX - 1);
	int imgY = clamp((int)(imgPos.y + 0.5f + cMaxVolumeIndexYZ / 2), 0, fftSizeY - 1);

	int index3D = z * (cMaxVolumeIndexYZ+1) * (cMaxVolumeIndexX+1) + y * (cMaxVolumeIndexX+1) + x;
	int index2D = imgY * fftSizeX + imgX;

	if (false) {
		const float* __restrict__ CTF = getNthItem(CTFs, space->projectionIndex, fftSizeX, fftSizeY, false);
		const float* __restrict__ modulator = getNthItem(modulators, space->projectionIndex, fftSizeX, fftSizeY, false);
		wCTF = CTF[index2D];
		wModulator = modulator[index2D];
	}

	float weight = wBlob * wModulator * dataWeight;


	/* FIXME
	 // use atomic as two blocks can write to same voxel
	atomicAdd(&tempVolumeGPU[2*index3D], img[2*index2D] * weight * wCTF);
	atomicAdd(&tempVolumeGPU[2*index3D + 1], img[2*index2D + 1] * weight * wCTF);
	atomicAdd(&tempWeightsGPU[index3D], weight);
	*/
}

/**
 * Method will map one voxel from the temporal
 * spaces to the given projection and update temporal spaces
 * using the pixel values of the projection withing the blob distance.
 */
__device__
void processVoxelBlob(
	float2* tempVolumeGPU, float *tempWeightsGPU,
	const float* FFTs,
	const float* CTFs, const float* modulators,
	int fftSizeX, int fftSizeY,
	int x, int y, int z,
	const RecFourierProjectionTraverseSpace* const space,
	const float* blobTableSqrt,
	int imgCacheDim)
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

	float w;
	float2 vol;
	vol.x = vol.y = w = 0.f;

#if !SHARED_IMG
	const float2* __restrict__ img = (const float2*)getNthItem(FFTs, space->projectionIndex, fftSizeX, fftSizeY, true);
#endif
	float dataWeight = space->weight;

	// ugly spaghetti code, but improves performance by app. 10%
	if (false) {
		const float* __restrict__ CTF = getNthItem(CTFs, space->projectionIndex, fftSizeX, fftSizeY, false);
		const float* __restrict__ modulator = getNthItem(modulators, space->projectionIndex, fftSizeX, fftSizeY, false);

		// check which pixel in the vicinity should contribute
		for (int i = minY; i <= maxY; i++) {
			float ySqr = (imgPos.y - i) * (imgPos.y - i);
			float yzSqr = ySqr + zSqr;
			if (yzSqr > radiusSqr) continue;
			for (int j = minX; j <= maxX; j++) {
				float xD = imgPos.x - j;
				float distanceSqr = xD*xD + yzSqr;
				if (distanceSqr > radiusSqr) continue;

#if SHARED_IMG
				int index2D = (i - SHARED_AABB[0].y) * imgCacheDim + (j-SHARED_AABB[0].x); // position in img - offset of the AABB
#else
				int index2D = i * fftSizeX + j;
#endif

				float wCTF = CTF[index2D];
				float wModulator = modulator[index2D];
#if PRECOMPUTE_BLOB_VAL
				int aux = (int) ((distanceSqr * cIDeltaSqrt + 0.5f));
	#if SHARED_BLOB_TABLE
				float wBlob = BLOB_TABLE[aux];
	#else
				float wBlob = blobTableSqrt[aux];
	#endif
#else
				float wBlob = kaiserValue(sqrtf(distanceSqr));
#endif
				float weight = wBlob * wModulator * dataWeight;
				w += weight;

#if SHARED_IMG
				vol += IMG[index2D] * weight * wCTF;
#else
				vol += img[index2D] * weight * wCTF;
#endif
			}
		}
	} else {
		// check which pixel in the vicinity should contribute
		for (int i = minY; i <= maxY; i++) {
			float ySqr = (imgPos.y - i) * (imgPos.y - i);
			float yzSqr = ySqr + zSqr;
			if (yzSqr > radiusSqr) continue;
			for (int j = minX; j <= maxX; j++) {
				float xD = imgPos.x - j;
				float distanceSqr = xD*xD + yzSqr;
				if (distanceSqr > radiusSqr) continue;

#if SHARED_IMG
				int index2D = (i - SHARED_AABB[0].y) * imgCacheDim + (j-SHARED_AABB[0].x); // position in img - offset of the AABB
#else
				int index2D = i * fftSizeX + j;
#endif

#if PRECOMPUTE_BLOB_VAL
				int aux = (int) ((distanceSqr * cIDeltaSqrt + 0.5f));
#if SHARED_BLOB_TABLE
				float wBlob = BLOB_TABLE[aux];
#else
				float wBlob = blobTableSqrt[aux];
#endif
#else
				float wBlob = kaiserValue2(distanceSqr) * cIw0;
#endif
				float weight = wBlob * dataWeight;
				w += weight;
#if SHARED_IMG
				vol += IMG[index2D] * weight;
#else
				vol += img[index2D] * weight;
#endif
			}
		}
	}
#if USE_ATOMICS
	// use atomic as two blocks can write to same voxel
	atomicAdd(&tempVolumeGPU[index3D].x, vol.x);
	atomicAdd(&tempVolumeGPU[index3D].y, vol.y);
	atomicAdd(&tempWeightsGPU[index3D], w);
#else
	tempVolumeGPU[index3D] += vol;
	tempWeightsGPU[index3D] += w;
#endif
}

/**
  * Method will process one projection image and add result to temporal
  * spaces.
  */
__device__
void processProjection(
	float2* tempVolumeGPU, float *tempWeightsGPU,
	const float* FFTs,
	const float* CTFs, const float* modulators,
	int fftSizeX, int fftSizeY,
	const RecFourierProjectionTraverseSpace* const tSpace,
	const float* devBlobTableSqrt,
	int imgCacheDim)
{


#if TILE > 1
	int id = threadIdx.y * blockDim.x + threadIdx.x;
	int tidX = threadIdx.x % TILE + (id / (blockDim.y * TILE)) * TILE;
	int tidY = (id / TILE) % blockDim.y;
	int idx = blockIdx.x*blockDim.x + tidX;
	int idy = blockIdx.y*blockDim.y + tidY;
//	printf("%d,%d (%d) -> %d,%d\n", threadIdx.x, threadIdx.y, id, idx, idy);
#else
	// map thread to each (2D) voxel
	volatile int idx = blockIdx.x*blockDim.x + threadIdx.x;
	volatile int idy = blockIdx.y*blockDim.y + threadIdx.y;
#endif





	if (tSpace->XY == tSpace->dir) { // iterate XY plane
		if (idy >= tSpace->minY && idy <= tSpace->maxY) {
			if (idx >= tSpace->minX && idx <= tSpace->maxX) {
				if (false) {
					float hitZ = getZ(idx, idy, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ);
					int z = (int)(hitZ + 0.5f); // rounding
					processVoxel(tempVolumeGPU, tempWeightsGPU, FFTs, CTFs, modulators, fftSizeX, fftSizeY, idx, idy, z, tSpace);
				} else {
					float z1 = getZ(idx, idy, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ); // lower plane
					float z2 = getZ(idx, idy, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ); // upper plane
					z1 = clamp(z1, 0, cMaxVolumeIndexYZ);
					z2 = clamp(z2, 0, cMaxVolumeIndexYZ);
					int lower = floorf(fminf(z1, z2));
					int upper = ceilf(fmaxf(z1, z2));
					for (int z = lower; z <= upper; z++) {
						processVoxelBlob(tempVolumeGPU, tempWeightsGPU, FFTs, CTFs, modulators, fftSizeX, fftSizeY, idx, idy, z, tSpace, devBlobTableSqrt, imgCacheDim);
					}
				}
			}
		}
	} else if (tSpace->XZ == tSpace->dir) { // iterate XZ plane
		if (idy >= tSpace->minZ && idy <= tSpace->maxZ) { // map z -> y
			if (idx >= tSpace->minX && idx <= tSpace->maxX) {
				if (false) {
					float hitY =getY(idx, idy, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ);
					int y = (int)(hitY + 0.5f); // rounding
					processVoxel(tempVolumeGPU, tempWeightsGPU, FFTs, CTFs, modulators, fftSizeX, fftSizeY, idx, y, idy, tSpace);
				} else {
					float y1 = getY(idx, idy, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ); // lower plane
					float y2 = getY(idx, idy, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ); // upper plane
					y1 = clamp(y1, 0, cMaxVolumeIndexYZ);
					y2 = clamp(y2, 0, cMaxVolumeIndexYZ);
					int lower = floorf(fminf(y1, y2));
					int upper = ceilf(fmaxf(y1, y2));
					for (int y = lower; y <= upper; y++) {
						processVoxelBlob(tempVolumeGPU, tempWeightsGPU, FFTs, CTFs, modulators, fftSizeX, fftSizeY, idx, y, idy, tSpace, devBlobTableSqrt, imgCacheDim);
					}
				}
			}
		}
	} else { // iterate YZ plane
		if (idy >= tSpace->minZ && idy <= tSpace->maxZ) { // map z -> y
			if (idx >= tSpace->minY && idx <= tSpace->maxY) { // map y > x
				if (false) {
					float hitX = getX(idx, idy, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ);
					int x = (int)(hitX + 0.5f); // rounding
					processVoxel(tempVolumeGPU, tempWeightsGPU, FFTs, CTFs, modulators, fftSizeX, fftSizeY, x, idx, idy, tSpace);
				} else {
					float x1 = getX(idx, idy, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ); // lower plane
					float x2 = getX(idx, idy, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ); // upper plane
					x1 = clamp(x1, 0, cMaxVolumeIndexX);
					x2 = clamp(x2, 0, cMaxVolumeIndexX);
					int lower = floorf(fminf(x1, x2));
					int upper = ceilf(fmaxf(x1, x2));
					for (int x = lower; x <= upper; x++) {
						processVoxelBlob(tempVolumeGPU, tempWeightsGPU, FFTs, CTFs, modulators, fftSizeX, fftSizeY, x, idx, idy, tSpace, devBlobTableSqrt, imgCacheDim);
					}
				}
			}
		}
	}
}

/**
 * Method calculates an Axis Aligned Bounding Box in the image space.
 * AABB is guaranteed to be big enough that all threads in the block,
 * while processing the traverse space, will not read image data outside
 * of the AABB
 */
__device__
void calculateAABB(const RecFourierProjectionTraverseSpace* tSpace, Point3D<float>* dest) {
	Point3D<float> box[8];
	// calculate AABB for the whole working block
	if (tSpace->XY == tSpace->dir) { // iterate XY plane
		box[0].x = box[3].x = box[4].x = box[7].x = blockIdx.x*blockDim.x - cBlobRadius;
		box[1].x = box[2].x = box[5].x = box[6].x = (blockIdx.x+1)*blockDim.x + cBlobRadius - 1.f;

		box[2].y = box[3].y = box[6].y = box[7].y = (blockIdx.y+1)*blockDim.y + cBlobRadius - 1.f;
		box[0].y = box[1].y = box[4].y = box[5].y = blockIdx.y*blockDim.y- cBlobRadius;

		box[0].z = getZ(box[0].x, box[0].y, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ);
		box[4].z = getZ(box[4].x, box[4].y, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ);

		box[3].z = getZ(box[3].x, box[3].y, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ);
		box[7].z = getZ(box[7].x, box[7].y, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ);

		box[2].z = getZ(box[2].x, box[2].y, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ);
		box[6].z = getZ(box[6].x, box[6].y, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ);

		box[1].z = getZ(box[1].x, box[1].y, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ);
		box[5].z = getZ(box[5].x, box[5].y, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ);
	} else if (tSpace->XZ == tSpace->dir) { // iterate XZ plane
		box[0].x = box[3].x = box[4].x = box[7].x = blockIdx.x*blockDim.x - cBlobRadius;
		box[1].x = box[2].x = box[5].x = box[6].x = (blockIdx.x+1)*blockDim.x + cBlobRadius - 1.f;

		box[2].z = box[3].z = box[6].z = box[7].z = (blockIdx.y+1)*blockDim.y + cBlobRadius - 1.f;
		box[0].z = box[1].z = box[4].z = box[5].z = blockIdx.y*blockDim.y- cBlobRadius;

		box[0].y = getY(box[0].x, box[0].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ);
		box[4].y = getY(box[4].x, box[4].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ);

		box[3].y = getY(box[3].x, box[3].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ);
		box[7].y = getY(box[7].x, box[7].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ);

		box[2].y = getY(box[2].x, box[2].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ);
		box[6].y = getY(box[6].x, box[6].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ);

		box[1].y = getY(box[1].x, box[1].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ);
		box[5].y = getY(box[5].x, box[5].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ);
	} else { // iterate YZ plane
		box[0].y = box[3].y = box[4].y = box[7].y = blockIdx.x*blockDim.x - cBlobRadius;
		box[1].y = box[2].y = box[5].y = box[6].y = (blockIdx.x+1)*blockDim.x + cBlobRadius - 1.f;

		box[2].z = box[3].z = box[6].z = box[7].z = (blockIdx.y+1)*blockDim.y + cBlobRadius - 1.f;
		box[0].z = box[1].z = box[4].z = box[5].z = blockIdx.y*blockDim.y- cBlobRadius;

		box[0].x = getX(box[0].y, box[0].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ);
		box[4].x = getX(box[4].y, box[4].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ);

		box[3].x = getX(box[3].y, box[3].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ);
		box[7].x = getX(box[7].y, box[7].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ);

		box[2].x = getX(box[2].y, box[2].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ);
		box[6].x = getX(box[6].y, box[6].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ);

		box[1].x = getX(box[1].y, box[1].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->bottomOriginX, tSpace->bottomOriginY, tSpace->bottomOriginZ);
		box[5].x = getX(box[5].y, box[5].z, tSpace->unitNormalX, tSpace->unitNormalY, tSpace->unitNormalZ, tSpace->topOriginX, tSpace->topOriginY, tSpace->topOriginZ);
	}
	// transform AABB to the image domain
	rotate(box, tSpace);
	// AABB is projected on image. Create new AABB that will encompass all vertices
	computeAABB(dest, box);
}

/**
 * Method returns true if AABB lies within the image boundaries
 */
__device__
bool isWithin(Point3D<float>* AABB, int imgXSize, int imgYSize) {
	return (AABB[0].x < imgXSize)
			&& (AABB[1].x >= 0)
			&& (AABB[0].y < imgYSize)
			&& (AABB[1].y >= 0);
}

/**
 * Method will load data from image at position tXindex, tYindex
 * and return them.
 * In case the data lies outside of the image boundaries, zeros (0,0)
 * are returned
 */
__device__
void getImgData(Point3D<float>* AABB,
		int tXindex, int tYindex,
		const float* FFTs, int fftSizeX, int fftSizeY,
		int imgIndex,
		float& vReal, float& vImag) {
	int imgXindex = tXindex + AABB[0].x;
	int imgYindex = tYindex + AABB[0].y;
	if ((imgXindex >=0)
			&& (imgXindex < fftSizeX)
			&& (imgYindex >=0)
			&& (imgYindex < fftSizeY))	{
		int index = imgYindex * fftSizeX + imgXindex; // copy data from image
		vReal = getNthItem(FFTs, imgIndex, fftSizeX, fftSizeY, true)[2*index];
		vImag = getNthItem(FFTs, imgIndex, fftSizeX, fftSizeY, true)[2*index + 1];

	} else {
		vReal = vImag = 0.f; // out of image bound, so return zero
	}
}

/**
 * Method will copy imgIndex(th) data from buffer
 * to given destination (shared memory).
 * Only data within AABB will be copied.
 * Destination is expected to be continuous array of sufficient
 * size (imgCacheDim^2)
 */
__device__
void copyImgToCache(float2* dest, Point3D<float>* AABB,
		const float* FFTs, int fftSizeX, int fftSizeY,
		int imgIndex,
		int imgCacheDim) {
	for (int y = threadIdx.y; y < imgCacheDim; y += blockDim.y) {
		for (int x = threadIdx.x; x < imgCacheDim; x += blockDim.x) {
			int memIndex = y * imgCacheDim + x;
			getImgData(AABB, x, y, FFTs, fftSizeX, fftSizeY, imgIndex, dest[memIndex].x, dest[memIndex].y);
		}
	}
}

/**
 * Method will use data stored in the buffer and update temporal
 * storages appropriately.
 */
extern "C" __global__
void processBufferKernel(
		float2* tempVolumeGPU, float *tempWeightsGPU,
		RecFourierProjectionTraverseSpace* spaces, int noOfSpaces,
		const float* FFTs,
		int fftSizeX, int fftSizeY,
		float* devBlobTableSqrt,
		int imgCacheDim) {
#if SHARED_BLOB_TABLE
	if ( ! false) {
		// copy blob table to shared memory
		volatile int id = threadIdx.y*blockDim.x + threadIdx.x;
		volatile int blockSize = blockDim.x * blockDim.y;
		for (int i = id; i < BLOB_TABLE_SIZE_SQRT; i+= blockSize)
			BLOB_TABLE[i] = devBlobTableSqrt[i];
		__syncthreads();
	}
#endif

	for (int i = blockIdx.z; i < noOfSpaces; i += gridDim.z) {
		RecFourierProjectionTraverseSpace* space = &spaces[i];

#if SHARED_IMG
		if ( ! false) {
			// make sure that all threads start at the same time
			// as they can come from previous iteration
			__syncthreads();
			if ((threadIdx.x == 0) && (threadIdx.y == 0)) {
				// first thread calculates which part of the image should be shared
				calculateAABB(space, SHARED_AABB);
			}
			__syncthreads();
			// check if the block will have to copy data from image
			if (isWithin(SHARED_AABB, fftSizeX, fftSizeY)) {
				// all threads copy image data to shared memory
				copyImgToCache(IMG, SHARED_AABB,
						FFTs, fftSizeX, fftSizeY,
						space->projectionIndex, imgCacheDim);
				__syncthreads();
			} else {
				continue; // whole block can exit, as it's not reading from image
			}
		}
#endif

		processProjection(
			tempVolumeGPU, tempWeightsGPU,
			FFTs,
			NULL, NULL, fftSizeX, fftSizeY,
			space,
			devBlobTableSqrt,
			imgCacheDim);
		__syncthreads(); // sync threads to avoid write after read problems
	}
}
