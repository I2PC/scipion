/*
 * fourier_reconstruction_traverse_space.h
 *
 *  Created on: Aug 21, 2017
 *      Author: david
 */

#ifndef XMIPP_LIBRARIES_DATA_FOURIER_RECONSTRUCTION_TRAVERSE_SPACE_H_
#define XMIPP_LIBRARIES_DATA_FOURIER_RECONSTRUCTION_TRAVERSE_SPACE_H_

#include <data/point3D.h>

struct TraverseSpace {

int minY, minX, minZ;
int maxY, maxX, maxZ;
float maxDistanceSqr;
enum Direction { XY, XZ, YZ } dir;
Point3D<float> u, v;
Point3D<float> topOrigin, bottomOrigin;
};

#endif /* XMIPP_LIBRARIES_DATA_FOURIER_RECONSTRUCTION_TRAVERSE_SPACE_H_ */
