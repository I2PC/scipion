#ifndef XMIPP_INTERNAL_DATA
#define XMIPP_INTERNAL_DATA
#include <jni.h>
#include <data/image.h>
#include <data/projection.h>

static jfieldID ImageDouble_peerId;
//static jfieldID ImageDouble_filenameId;

static jfieldID Projection_peerId;
//static jfieldID Projection_filenameId;

static jfieldID MetaData_peerId;
//static jfieldID MetaData_filenameId;

#define peerId ImageDouble_peerId
#define peerId Projection_peerId
#define peerId MetaData_peerId

#define GET_INTERNAL_IMAGE(obj) ((Image<double> *)(env->GetLongField(obj, peerId)))
#define GET_INTERNAL_PROJECTION(obj) ((Projection *)(env->GetLongField(obj, peerId)))
#define GET_INTERNAL_METADATA(obj) ((MetaData*)(env->GetLongField(obj, peerId)))

#endif
