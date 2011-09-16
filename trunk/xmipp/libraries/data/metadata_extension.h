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
                   double& _max, bool apply_geo);

/** Get the average of a Metadata applying the header.
 * The MD is not cleaned from disabled images (this option makes the call faster).
 */
void getAverageApplyGeo(const MetaData &MD, MultidimArray<double> & _ave);

/** Get the image statistics of a metadata.
 */
void getStatistics(MetaData MD, double& _ave, double& _sd, double& _min,
                   double& _max, bool apply_geo);

/** Get Fourier statistics */
void getFourierStatistics(MetaData &MDin, double sam, MetaData &Mdout,
                          bool do_dpr, double max_sam);

/** Get image size
 */
void ImgSize(const MetaData &MD, int &Xdim, int &Ydim, int &Zdim, size_t &Ndim);

/** Get image size and data type */
void ImgSize(const MetaData &MD, int &Xdim, int &Ydim, int &Zdim, size_t &Ndim, DataType &datatype);

/** Get image size and data type of a Metadata file */
void ImgSize(const FileName &filename, int &Xdim, int &Ydim, int &Zdim, size_t &Ndim);

/// compare two image files
bool ImgCompare(const FileName &filename1, const FileName &filename2);

/// compare if same dimensions
bool ImgCompareSize(const FileName &filename1, const FileName &filename2);

/** Maximum length of the filenames inside */
int MaxFileNameLength(MetaData &MD);

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
