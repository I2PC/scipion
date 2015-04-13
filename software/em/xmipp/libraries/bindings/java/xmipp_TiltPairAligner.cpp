/*
 * xmipp_TiltPairAligner.cpp
 *
 *  Created on: Sep 20, 2011
 *      Author: airen
 */
#
#include "xmipp_java_adapter.h"
#include <jni.h>
#include "xmipp_TiltPairAligner.h"
#include "data/micrograph.h"
#include "xmipp_ExceptionsHandler.h"
#include "xmipp_InternalData.h"
#include <iostream>

static jfieldID TiltPairAligner_peerId;

JNIEXPORT void JNICALL Java_xmipp_jni_TiltPairAligner_storeIds(JNIEnv *env,
		jclass cls) {
	TiltPairAligner_peerId = env->GetFieldID(cls, "peer", "J");
}

JNIEXPORT void JNICALL Java_xmipp_jni_TiltPairAligner_create(JNIEnv *env,
		jobject jobj) {
	TiltPairAligner * tpa = new TiltPairAligner();
	env->SetLongField(jobj, TiltPairAligner_peerId, (long) tpa);
}

JNIEXPORT void JNICALL Java_xmipp_jni_TiltPairAligner_destroy(JNIEnv *env,
		jobject jobj) {
	TiltPairAligner * tpa = GET_INTERNAL_TPA(jobj);
	delete tpa;
	tpa = NULL;
	env->SetLongField(jobj, TiltPairAligner_peerId, (long) tpa);
}

JNIEXPORT void JNICALL Java_xmipp_jni_TiltPairAligner_addParticleToAligner(
		JNIEnv * env, jobject jobj, jint jx1, jint jy1, jint jx2, jint jy2) {

	String msg;
	TiltPairAligner * tpa = GET_INTERNAL_TPA(jobj);
	if (tpa != NULL) {
		try {
			tpa->adjustPassingMatrix(jx1, jy1, jx2, jy2);
		} catch (XmippError &xe) {
			msg = xe.getDefaultMessage();
		} catch (std::exception& e) {
			msg = e.what();
		} catch (...) {
			msg = "Unhandled exception";
		}
	} else {
		msg = "TiltPairAligner is null";
	}

	// If there was an exception, sends it to java environment.
	if (!msg.empty()) {
		handleXmippException(env, msg);
	}

}

JNIEXPORT jobject JNICALL Java_xmipp_jni_TiltPairAligner_getTiltedParticle(
		JNIEnv *env, jobject jobj, jint jx1, jint jy1) {
	String msg;

	try {
		int x = 0, y = 0;
		TiltPairAligner * tpa = GET_INTERNAL_TPA(jobj);
		if (tpa != NULL) {
			tpa->passToTilted(jx1, jy1, x, y);
			jclass pclass = env->FindClass("xmipp/jni/Particle");
			jmethodID constructor = env->GetMethodID(pclass, "<init>", "(II)V");
			jobject particle = env->NewObject(pclass, constructor, x, y);

			return particle;
		} else {
			msg = "TiltPairAligner is null";
		}
	} catch (XmippError &xe) {
		msg = xe.getDefaultMessage();
	} catch (std::exception& e) {
		msg = e.what();
	} catch (...) {
		msg = "Unhandled exception";
	}

	// If there was an exception, sends it to java environment.
	if (!msg.empty()) {
		handleXmippException(env, msg);
	}
	return NULL;
}

JNIEXPORT void JNICALL Java_xmipp_jni_TiltPairAligner_clear(JNIEnv *env,
		jobject jobj) {
	String msg;

	try {
		TiltPairAligner * tpa = GET_INTERNAL_TPA(jobj);
		if (tpa != NULL) {
			tpa->clear();
		} else {
			msg = "TiltPairAligner is null";
		}
	} catch (XmippError &xe) {
		msg = xe.getDefaultMessage();
	} catch (std::exception& e) {
		msg = e.what();
	} catch (...) {
		msg = "Unhandled exception";
	}

	// If there was an exception, sends it to java environment.
	if (!msg.empty()) {
		handleXmippException(env, msg);
	}
}

JNIEXPORT jdoubleArray JNICALL Java_xmipp_jni_TiltPairAligner_computeAngles(JNIEnv *env,
		jobject jobj) {
	String msg;

	try {
		TiltPairAligner * tpa = GET_INTERNAL_TPA(jobj);
		if (tpa != NULL) {

			double alphas[3];
			tpa->computeGamma();
			alphas[2]=tpa->gamma;
			tpa->computeAngles(alphas[0], alphas[1], alphas[2]);
			jdoubleArray result = env->NewDoubleArray(3);
			env->SetDoubleArrayRegion(result, 0, 3, alphas);
			return result;
		} else {
			msg = "TiltPairAligner is null";
		}
	} catch (XmippError &xe) {
		msg = xe.getDefaultMessage();
	} catch (std::exception& e) {
		msg = e.what();
	} catch (...) {
		msg = "Unhandled exception";
	}

	// If there was an exception, sends it to java environment.
	if (!msg.empty()) {
		handleXmippException(env, msg);
	}
	return NULL;
}

