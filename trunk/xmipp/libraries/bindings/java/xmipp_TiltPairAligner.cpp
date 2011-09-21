/*
 * xmipp_TiltPairAligner.cpp
 *
 *  Created on: Sep 20, 2011
 *      Author: airen
 */
#include <jni.h>
#include "xmipp_TiltPairAligner.h"
#include "data/micrograph.h"
#include "xmipp_ExceptionsHandler.h"

JNIEXPORT void JNICALL Java_xmipp_TiltPairAligner_addParticleToAligner(JNIEnv * env,
		jobject jobj, jint jx1, jint jy1, jint jx2, jint jy2) {

	String msg;

	try {

		TiltPairAligner tpa;
		tpa.adjust_passing_matrix(jx1, jx2, jy1, jy2);
	} catch (XmippError xe) {
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

