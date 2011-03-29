#ifndef XMIPP_INTERNAL_DATA
#define XMIPP_INTERNAL_DATA
#include <jni.h>
#include <data/image.h>
#include <data/projection.h>

static jfieldID ImageDouble_peerId;

static jfieldID Projection_peerId;

static jfieldID MetaData_peerId;

#define GET_INTERNAL_IMAGE(obj) ((Image<double> *)(env->GetLongField(obj, ImageDouble_peerId)))
#define GET_INTERNAL_PROJECTION(obj) ((Projection *)(env->GetLongField(obj, Projection_peerId)))
#define GET_INTERNAL_METADATA(obj) ((MetaData*)(env->GetLongField(obj, MetaData_peerId)))

#endif
