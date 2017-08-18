
/** Struct representing all data regarding one projection */

#ifndef XMIPP_DATA_PROJECTION_DATA_H_
#define XMIPP_DATA_PROJECTION_DATA_H_

#include <data/array_2D.h>
#include <complex>
#include <data/matrix2d.h>

typedef float MATRIX[3][3];

struct ProjectionData
{
	Array2D<std::complex<float> >* img;
	Array2D<float>* CTF;
	Array2D<float>* modulator;
	int imgIndex;
	float weight;
	Matrix2D<double> localAInv;
	Matrix2D<double> localA;
	bool skip;
public:
	ProjectionData() {
		img = 0;
		CTF = modulator = 0;
		skip = true;
		weight = 0;
		imgIndex = -1;
	}
	~ProjectionData() {
		clean();
	}
	/** Remove stored data and set to skip */
	void clean() {
		delete img;
		delete CTF;
		delete modulator;
		img = 0;
		CTF = modulator = 0;
		skip = true;
	}
};

#endif /* XMIPP_DATA_PROJECTION_DATA_H_ */
