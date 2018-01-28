/***************************************************************************
 *
 * Authors:     David Strelak (davidstrelak@gmail.com)
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

#ifndef XMIPP_LIBRARIES_DATA_RECONSTRUCT_FOURIER_PROJECTION_TRAVERSE_SPACE_H_
#define XMIPP_LIBRARIES_DATA_RECONSTRUCT_FOURIER_PROJECTION_TRAVERSE_SPACE_H_

#include "data/point3D.h"

/**
 * Struct describing a how to best traverse a single projection
 * during the Fourier Reconstruction.
 * It describes an Axis Aligned Bounding Box (AABB) that wraps
 * some projection plane (or cuboid, in case it has some thickness) oriented in the memory
 */
struct RecFourierProjectionTraverseSpace {

int minY, minX, minZ; // coordinates of bottom left front corner of the AABB
int maxY, maxX, maxZ; // coordinates of upper right back corner of the AABB
float maxDistanceSqr; // max squared distance from the center of the Fourier Space which should be processed (i.e. 'max frequency to process')
enum Direction { XY, XZ, YZ } dir; // optimal plane for traversing (i.e. you process this plane and iterate in last direction)
/**
 * Projection itself is a plane (with/without some thickness) somehow oriented in the AABB.
 * These variables hold normal to the plane
 */
Point3D<float> unitNormal;
/**
 * Projection can have some thickness due to the blob radius.
 * These variables hold the origin of the lower/upper plane.
 * Bear in mind that 'lower/upper' refers to initial orientation, before applying projection
 * specific rotation, so it can happen that 'topOrigin' is lower than 'bottomOrigin'.
 * In case the blob radius is zero, these variables hold the same point
 */
Point3D<float> topOrigin, bottomOrigin;
int projectionIndex; // index to array of projections, which holds appropriate (visual) data
float transformInv[3][3]; // rotation matrix, describing transformation to default position
float weight; // of the projection
};

#endif /* XMIPP_LIBRARIES_DATA_RECONSTRUCT_FOURIER_PROJECTION_TRAVERSE_SPACE_H_ */
