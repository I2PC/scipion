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

void getStatistics(MetaData &MT, Image<double> & _ave, Image<double> & _sd, double& _min,
                    double& _max, bool apply_geo);

/** Get image size
 *
 */

static int null_object=-1;

void SingleImgSize(const FileName &filename, int &Xdim, int &Ydim=null_object,
             int &Zdim=null_object,
             int &Ndim=null_object);

void ImgSize(const MetaData &MD, int &Xdim, int &Ydim=null_object,
             int &Zdim=null_object,
             int &Ndim=null_object);

void ImgSize(const FileName &filename, int &Xdim, int &Ydim=null_object,
             int &Zdim=null_object,
             int &Ndim=null_object);

int MaxFileNameLength(MetaData &MD);

void mpiSelectPart(MetaData &md, int rank, int size, int &num_img_tot);

///Swig interfaces
template<class T>
bool setValueSwig(MetaData &md,  MDLabel label,  T &valueIn, long int objectID = -1)
{
    md.setValue(label, valueIn, objectID);
}

template<class T>
bool getValueSwig(const MetaData &md, const MDLabel label, T &valueOut, long int objectID = -1)
{
    md.getValue(label, valueOut, objectID);
}

#endif /* METADATA_EXTENSION_H_ */
