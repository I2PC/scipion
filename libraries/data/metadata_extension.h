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

/** @defgroup MetaDataExtension Extension to Metadata Stuff
 * @ingroup DataLibrary
 * @{
 */
/** Get the image statistics of a metadata.
 * Note that the mean and stddev are images, not values.*/
void getStatistics(MetaData MD, Image<double> & _ave, Image<double> & _sd, double& _min,
                   double& _max, bool apply_geo, MDLabel image_label=MDL_IMAGE);

/** Get the average of a Metadata applying the header.
 * The MD is not cleaned from disabled images (this option makes the call faster).
 */
void getAverageApplyGeo(MetaData MD, MultidimArray<double> & _ave, MDLabel image_label=MDL_IMAGE);

/** Get the image statistics of a metadata.
 */
void getStatistics(MetaData MD, double& _ave, double& _sd, double& _min,
                   double& _max, bool apply_geo, MDLabel image_label=MDL_IMAGE);

/** Get Fourier statistics */
void getFourierStatistics(MetaData &MDin, double sam, MetaData &Mdout,
                          bool do_dpr, double max_sam, MDLabel image_label=MDL_IMAGE);

/** Get image size
 */
void getImageSize(const MetaData &MD, size_t &Xdim, size_t &Ydim, size_t &Zdim, size_t &Ndim, MDLabel image_label=MDL_IMAGE);

/** Get image size and data type */
void getImageInfo(const MetaData &MD, size_t &Xdim, size_t &Ydim, size_t &Zdim, size_t &Ndim, DataType &datatype, MDLabel image_label=MDL_IMAGE);

void getImageInfo(const MetaData &MD, ImageInfo &imgInfo, MDLabel image_label=MDL_IMAGE);

/** Get image size and data type of a Metadata file */
void getImageSizeFromFilename(const FileName &filename, size_t &Xdim, size_t &Ydim, size_t &Zdim, size_t &Ndim, MDLabel image_label=MDL_IMAGE);

/// compare two image files
bool compareImage(const FileName &filename1, const FileName &filename2);

/// compare if same dimensions
bool compareImageSize(const FileName &filename1, const FileName &filename2);

/** Compare two metadata files */
bool compareTwoMetadataFiles(const FileName &fn1, const FileName &fn2);

/**Copy images listed in metadata to some location using the
 * same logic as in xmipp_image_convert program
 * if independent is false, output will be treated as the output stack name
 * if is true, will be the prefix for the output images
*/
void copyImages(const MetaData &md, const char * output, bool independent, MDLabel image_label=MDL_IMAGE);

/** Maximum length of the filenames inside */
int maxFileNameLength(const MetaData &MD, MDLabel image_label=MDL_IMAGE);

/** Choose a part of the metadata for MPI */
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

//@}

#endif /* METADATA_EXTENSION_H_ */
