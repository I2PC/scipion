/*
 * metadata_extension.h
 *
 *  Created on: May 12, 2010
 *      Author: roberto
 */

#ifndef METADATA_EXTENSION_H_
#define METADATA_EXTENSION_H_

#include "xmipp_filename.h"
#include "xmipp_image.h"
#include "metadata.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#define BLOCKROW 1
#define BLOCKCOLUNM 2
#define BLOCKMIXED 3
#define BLOCKDIFFERENCE 30
/** @defgroup MetaDataExtension Extension to Metadata Stuff
 * @ingroup DataLibrary
 * @{
 */
/** Get the image statistics of a metadata.
 * Note that the mean and stddev are images, not values.*/
void getStatistics(MetaData md, Image<double> & _ave, Image<double> & _sd, bool apply_geo, bool wrap, MDLabel image_label=MDL_IMAGE);

/** Write images in MetaData into a stack */
void writeMdToStack(const MetaData &md, const FileName &fnStack, bool apply_geo, bool wrap, MDLabel image_label=MDL_IMAGE);

/** Get the average of a Metadata applying the header.
 * The md is not cleaned from disabled images (this option makes the call faster).
 */
void getAverageApplyGeo(MetaData md, MultidimArray<double> & _ave, MDLabel image_label=MDL_IMAGE);

/** Get the image statistics of a metadata.
 */
void getStatistics(MetaData md, double& _ave, double& _sd, double& _min,
        double& _max, bool apply_geo, MDLabel image_label=MDL_IMAGE);

/** Get Fourier statistics */
void getFourierStatistics(MetaData &MDin, double sam, MetaData &Mdout,
                          bool do_dpr, double max_sam, MDLabel image_label=MDL_IMAGE);

/** Get image size
 */
void getImageSize(const MetaData &md, size_t &Xdim, size_t &Ydim, size_t &Zdim, size_t &Ndim, MDLabel image_label=MDL_IMAGE);

/** Get image size and data type */
void getImageInfo(const MetaData &md, size_t &Xdim, size_t &Ydim, size_t &Zdim, size_t &Ndim, DataType &datatype, MDLabel image_label=MDL_IMAGE);

void getImageInfo(const MetaData &md, ImageInfo &imgInfo, MDLabel image_label=MDL_IMAGE);

/** Get image size and data type of a Metadata file */
void getImageSizeFromFilename(const FileName &filename, size_t &Xdim, size_t &Ydim, size_t &Zdim, size_t &Ndim, MDLabel image_label=MDL_IMAGE);

/// compare two image files
bool compareImage(const FileName &filename1, const FileName &filename2);

/// compare if same dimensions
bool compareImageSize(const FileName &filename1, const FileName &filename2);

/** Compare two metadata files */
bool compareTwoMetadataFiles(const FileName &fn1, const FileName &fn2);


/** Maximum length of the filenames inside */
int maxFileNameLength(const MetaData &md, MDLabel image_label=MDL_IMAGE);

/** Choose a part of the metadata for MPI */
void mpiSelectPart(MetaData &md, int rank, int size, int &num_img_tot);

/** Read a 1 or two column list of micrographs.
 *  Two column files are interpreted as Random Conical Tilt pairs.
 */
void readMetaDataWithTwoPossibleImages(const FileName &fn, MetaData &md);

/** Substitute original images.
 *  This function substitutes the images in a given label of the metadata by
 *  a set of original images. The substituted images are supposed to be in a stack produced
 *  by processing the original images. The substitution is performed in all blocks.
 */
void substituteOriginalImages(const FileName &fn, const FileName &fnOrig, const FileName &fnOut,
                              MDLabel label, bool skipFirstBlock);

/**convert bsoft metadas to xmipp style
 *
 */
void bsoftRemoveLoopBlock(const FileName &_inFile, const FileName &block);

/** undo bsoftRemoveLoopBlock actions
 *
 */
void bsoftRestoreLoopBlock(const FileName &_inFile, const FileName &block);

Matrix2D<double> getMatrix(char* matrix);

MDRow firstRow(const FileName &fnMetadata);
//@}

#endif /* METADATA_EXTENSION_H_ */
