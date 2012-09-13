/*
 * metadata_extension.h
 *
 *  Created on: May 12, 2010
 *      Author: roberto
 */

#ifndef METADATA_EXTENSION_H_
#define METADATA_EXTENSION_H_

#include "image.h"
#include "metadata.h"

void get_statistics(MetaData MT, Image<double> & _ave, Image<double> & _sd, double& _min,
                    double& _max, bool apply_geo);

/** Get image size
 *
 */
static int null_object=-1;

void ImgSize(MetaData &MD, int &Xdim, int &Ydim=null_object,
             int &Zdim=null_object,
             int &Ndim=null_object);


#endif /* METADATA_EXTENSION_H_ */
