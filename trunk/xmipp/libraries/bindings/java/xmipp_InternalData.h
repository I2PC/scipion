#ifndef XMIPP_INTERNAL_DATA
#define XMIPP_INTERNAL_DATA
#include <jni.h>
#include <data/xmipp_image.h>
#include <data/xmipp_image_generic.h>
#include <data/projection.h>
#include <data/ctf.h>

static jfieldID ImageDouble_peerId;
static jfieldID ImageGeneric_peerId;
static jfieldID Projection_peerId;
static jfieldID MetaData_peerId;
static jfieldID CTFDescription_peerId;
static jfieldID TiltPairAligner_peerId;

#define GET_INTERNAL_IMAGE(obj) ((Image<double> *)(env->GetLongField(obj, ImageDouble_peerId)))
#define GET_INTERNAL_IMAGE_GENERIC(obj) ((ImageGeneric *)(env->GetLongField(obj, ImageGeneric_peerId)))
#define GET_INTERNAL_PROJECTION(obj) ((Projection *)(env->GetLongField(obj, Projection_peerId)))
#define GET_INTERNAL_METADATA(obj) ((MetaData *)(env->GetLongField(obj, MetaData_peerId)))
#define GET_INTERNAL_CTFDESCRIPTION(obj) ((CTFDescription *)(env->GetLongField(obj, CTFDescription_peerId)))
#define GET_INTERNAL_TPA(obj) ((TiltPairAligner *)(env->GetLongField(obj, TiltPairAligner_peerId)))

#endif
