/*
 * metadata_extension.h
 *
 *  Created on: May 12, 2010
 *      Author: roberto
 */

#ifndef METADATA_EXTENSION_H_
#define METADATA_EXTENSION_H_

#include "filename.h"
#include "image.h"
#include "metadata.h"

void getStatistics(MetaData &MT, Image<double> & _ave, Image<double> & _sd, double& _min,
                    double& _max, bool apply_geo);

/** Get image size
 *
 */
void ImgSize(const MetaData &MD, int &Xdim, int &Ydim, int &Zdim, unsigned long &Ndim);

void ImgSize(const FileName &filename, int &Xdim, int &Ydim, int &Zdim, unsigned long &Ndim);

void getBlocksAvailableInMetaData(const FileName &inFile, StringVector& blockList);

int MaxFileNameLength(MetaData &MD);

void mpiSelectPart(MetaData &md, int rank, int size, int &num_img_tot);

/** Read a 1 or two column list of micrographs.
 *  Two column files are interpreted as Random Conical Tilt pairs.
 */
void readMetaDataWithTwoPossibleImages(const FileName &fn, MetaData &MD);

/** Substitute original images.
 *  This function substitutes the images in a given label of the metadata by
 *  a set of original images. The substituted images are supposed to be in a stack produced
 *  by processing the original images. The substitution is performed in all blocks.
 */
void substituteOriginalImages(const FileName &fn, const FileName &fnOrig, const FileName &fnOut,
		MDLabel label, bool skipFirstBlock);

#endif /* METADATA_EXTENSION_H_ */
