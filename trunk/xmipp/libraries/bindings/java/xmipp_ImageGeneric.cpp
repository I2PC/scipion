#include <jni.h>
#include "xmipp_ImageGeneric.h"
#include "xmipp_InternalData.h"
#include "xmipp_ExceptionsHandler.h"
#include <data/xmipp_image_generic.h>
#include <data/xmipp_fft.h>

JNIEXPORT void JNICALL Java_xmipp_ImageGeneric_storeIds
(JNIEnv *env, jclass cls) {
	ImageGeneric_peerId = env->GetFieldID(cls, "peer", "J");
}

JNIEXPORT void JNICALL Java_xmipp_ImageGeneric_create
(JNIEnv *env, jobject jobj) {
	ImageGeneric *image = new ImageGeneric();
	env->SetLongField(jobj, ImageGeneric_peerId, (long)image);
}

JNIEXPORT void JNICALL Java_xmipp_ImageGeneric_destroy
(JNIEnv *env, jobject jobj) {
	ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
	delete image;
	image = NULL;
	env->SetLongField(jobj, ImageGeneric_peerId, (long)image);
}

JNIEXPORT void JNICALL Java_xmipp_ImageGeneric_read
(JNIEnv *env, jobject jobj, jstring filename) {
	std::string msg = "";
	ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

	if (image != NULL) {
		try {
			const char * fnStr = env->GetStringUTFChars(filename, false);

			image->readMapped(fnStr);
		} catch (XmippError xe) {
			msg = xe.getDefaultMessage();
		} catch (std::exception& e) {
			msg = e.what();
		} catch (...) {
			msg = "Unhandled exception";
		}
	} else {
		msg = "Metadata is null";
	}

	// If there was an exception, sends it to java environment.
	if(!msg.empty()) {
		handleXmippException(env, msg);
	}
}

JNIEXPORT jdoubleArray JNICALL Java_xmipp_ImageGeneric_getStatistics(
		JNIEnv *env, jobject jobj) {
	std::string msg = "";
	ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

	if (image != NULL) {
		try {
			double avg, stddev, min, max;

			(*image->data).computeStats(avg, stddev, min, max);

			// Copies vector into array.
			double statistics[4];
			statistics[0] = min;
			statistics[1] = max;
			statistics[2] = avg;
			statistics[3] = stddev;

			// Sets array value
			jdoubleArray array = env->NewDoubleArray(4);
			env->SetDoubleArrayRegion(array, 0, 4, statistics);

			return array;
		} catch (XmippError xe) {
			msg = xe.getDefaultMessage();
		} catch (std::exception& e) {
			msg = e.what();
		} catch (...) {
			msg = "Unhandled exception";
		}
	} else {
		msg = "Metadata is null";
	}

	// If there was an exception, sends it to java environment.
	if (!msg.empty()) {
		handleXmippException(env, msg);
	}

	return NULL;
}
