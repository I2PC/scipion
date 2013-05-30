#ifndef XMIPP_INTERNAL_DATA
#define XMIPP_INTERNAL_DATA
#include <jni.h>
#include <data/xmipp_image.h>
#include <data/xmipp_image_generic.h>
#include <data/ctf.h>

#define IDFIELD "peer"

#define STORE_PEER_ID(jobj, id) (env->SetLongField(jobj, env->GetFieldID(env->GetObjectClass(jobj), IDFIELD, "J"), (long)id))
#define GET_INTERNAL_POINTER(jobj) (env->GetLongField(jobj, env->GetFieldID(env->GetObjectClass(jobj), IDFIELD, "J")))

#define GET_INTERNAL_IMAGE_GENERIC(jobj) ((ImageGeneric *)GET_INTERNAL_POINTER(jobj))
#define GET_INTERNAL_METADATA(jobj) ((MetaData *)GET_INTERNAL_POINTER(jobj))
#define GET_INTERNAL_AUTOPARTICLEPICKING2(jobj) ((AutoParticlePicking2 *)GET_INTERNAL_POINTER(jobj))
#define GET_INTERNAL_CTFDESCRIPTION(jobj) ((CTFDescription *)GET_INTERNAL_POINTER(jobj))

// @TODO Convert this to the way above.
#define GET_INTERNAL_TPA(obj) ((TiltPairAligner *)(env->GetLongField(obj, TiltPairAligner_peerId)))
#define GET_INTERNAL_PTA(obj) ((ProgTomographAlignment *)(env->GetLongField(obj, ProgTomographAlignment_peerId)))

#endif
