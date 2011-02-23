#ifndef XMIPP_INTERNAL_DATA
#define XMIPP_INTERNAL_DATA
#include <jni.h>
#include <data/image.h>
#include <data/projection.h>

static jfieldID ImageDouble_peerId;
//static jfieldID ImageDouble_filenameId;

static jfieldID Projection_peerId;
//static jfieldID Projection_filenameId;

#define peerId ImageDouble_peerId
#define peerId Projection_peerId

#ifndef GET_INTERNAL_IMAGE
#define GET_INTERNAL_IMAGE(obj) ((Image<double> *)(env->GetLongField(obj, peerId)))
#endif
#ifndef GET_INTERNAL_PROJECTION
#define GET_INTERNAL_PROJECTION(obj) ((Projection *)(env->GetLongField(obj, peerId)))
#endif

#endif
