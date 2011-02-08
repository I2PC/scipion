#include <jni.h>
#include "xmipp_Projection.h"
#include "xmipp_ImageDouble.h"
#include "xmipp_ExceptionsHandler.h"
#include <data/image.h>
#include <data/projection.h>

static jfieldID Projection_peerId;
static jfieldID Projection_filenameId;

#define peerId Projection_peerId
#define GET_INTERNAL_PROJECTION() ((Projection *)(env->GetLongField(obj, peerId)))

JNIEXPORT void JNICALL Java_xmipp_Projection_reset
(JNIEnv *env, jobject obj, jint h, jint w) {
	Projection *projection = GET_INTERNAL_PROJECTION();

	if(projection != NULL) {
		projection->reset((int) h, (int) w);
	} else {
		handleXmippException(env, "projection is NULL");
	}
}

JNIEXPORT void JNICALL Java_xmipp_Projection_setXmippOrigin
(JNIEnv *env, jobject obj) {
	Projection *projection = GET_INTERNAL_PROJECTION();

	if(projection != NULL) {
		projection->data.setXmippOrigin();
	} else {
		handleXmippException(env, "projection is NULL");
	}
}

JNIEXPORT void JNICALL Java_xmipp_Projection_projectVolume
(JNIEnv *env, jclass class_, jobject jvolume, jobject jprojection, jint w, jint h, jdouble rot, jdouble tilt, jdouble pshi) {
	std::string msg = "";
	Image<double> *volume = *(Image<double> **)&jvolume;
//	Image<double> *volume = GET_INTERNAL_IMAGE();
	Projection *projection = *(Projection **)&jprojection;
//	Projection *projection = GET_INTERNAL_PROJECTION();

	if(volume != NULL) {
		if(projection != NULL) {
			try {
				//project_Volume(volume->data, &projection, (int) w, (int) h, (double) rot, (double) tilt, (double) pshi);
				//project_Volume(*arg1,*arg2,arg3,arg4,arg5,arg6,arg7);
			} catch (XmippError xe) {
				msg = xe.getDefaultMessage();
			} catch (std::exception& e) {
				msg = e.what();
			} catch (...) {
				msg = "Unhandled exception";
			}
		} else {
			msg = "projection is NULL";
		}
	} else {
		msg = "volume is NULL";
	}

	// If there was an exception, sends it to java environment.
	if(!msg.empty()) {
		handleXmippException(env, msg);
	}
}
