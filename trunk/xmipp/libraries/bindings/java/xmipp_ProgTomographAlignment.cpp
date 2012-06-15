#include "xmipp_java_adapter.h"

#include <jni.h>
#include <iostream>
#include <unistd.h>
#include "xmipp_InternalData.h"
#include "xmipp_ExceptionsHandler.h"
#include "xmipp_ProgTomographAlignment.h"
#include <reconstruction/tomo_align_tilt_series.h>

#define MAXPATHLEN 1024

JNIEXPORT void JNICALL Java_xmipp_jni_ProgTomographAlignment_storeIds
  (JNIEnv *env, jclass cl){
  	ProgTomographAlignment_peerId = env->GetFieldID(cl, "peer", "J");
  }

JNIEXPORT void JNICALL Java_xmipp_jni_ProgTomographAlignment_create
  (JNIEnv *env, jobject jobj){
  	ProgTomographAlignment *program = new ProgTomographAlignment();
	env->SetLongField(jobj, ProgTomographAlignment_peerId, (long)program);
	// default parameters for running produceInfo successfully
	/*program->lastStep=1;
	program->numThreads = 1;*/
	// default parameters for running the whole program successfully
	program->lastStep = -1;
	program->numThreads = 1;
  }

JNIEXPORT void JNICALL Java_xmipp_jni_ProgTomographAlignment_setInputFilename
  (JNIEnv *env, jobject jobj, jstring inputFileName){
  ProgTomographAlignment *program = GET_INTERNAL_PTA(jobj);
  const char * fnStr = env->GetStringUTFChars(inputFileName, false);
  program->fnSel=fnStr;
  }

JNIEXPORT void JNICALL Java_xmipp_jni_ProgTomographAlignment_setRoot
  (JNIEnv *env, jobject jobj, jstring root){
	  ProgTomographAlignment *program = GET_INTERNAL_PTA(jobj);
	  const char * fnStr = env->GetStringUTFChars(root, false);
	  program->fnRoot=fnStr;
  }

std::string get_working_path()
{
   char temp[MAXPATHLEN];
   return ( getcwd(temp, MAXPATHLEN) ? std::string( temp ) : std::string("") );
}


JNIEXPORT void JNICALL Java_xmipp_jni_ProgTomographAlignment_produceSideInfo
  (JNIEnv *env, jobject jobj){
	std::string msg = "";
  	try{
  	  ProgTomographAlignment *program = GET_INTERNAL_PTA(jobj);
  	  // std::cout << "file: " << program->fnSel << ". Root: " << program->fnRoot << ". Cwd: " << get_working_path()<< std::endl;
  	  program->produceSideInfo();
  	  // std::cout << "produceSideInfo finished" << std::endl;
  		} catch (XmippError xe) {
			msg = xe.getDefaultMessage();
		} catch (std::exception& e) {
			msg = e.what();
		} catch (...) {
			msg = "Unhandled exception";
	}

	// If there was an exception, sends it to java environment.
	if(!msg.empty()) {
		handleXmippException(env, msg);
	}
  }

JNIEXPORT void JNICALL Java_xmipp_jni_ProgTomographAlignment_writeTransformations
  (JNIEnv *env, jobject jobj, jstring fileName){
	std::string msg = "";
  	try{
  	  const char * fnStr = env->GetStringUTFChars(fileName, false);
  	  ProgTomographAlignment *program = GET_INTERNAL_PTA(jobj);
  	  program->writeTransformations(fnStr);
  		} catch (XmippError xe) {
			msg = xe.getDefaultMessage();
		} catch (std::exception& e) {
			msg = e.what();
		} catch (...) {
			msg = "Unhandled exception";
	}

	// If there was an exception, sends it to java environment.
	if(!msg.empty()) {
		handleXmippException(env, msg);
	}
}

JNIEXPORT jint JNICALL Java_xmipp_jni_ProgTomographAlignment_getIteration
  (JNIEnv *env, jobject jobj){
	ProgTomographAlignment *program = GET_INTERNAL_PTA(jobj);
	return program->iteration;

}

JNIEXPORT void JNICALL Java_xmipp_jni_ProgTomographAlignment_destroy
  (JNIEnv *env, jobject jobj){
  	std::cerr<<" destroying"<<std::endl;
	ProgTomographAlignment *program = GET_INTERNAL_PTA(jobj);
delete program;
	program = NULL;
	env->SetLongField(jobj, ProgTomographAlignment_peerId, (long)program);
  }


JNIEXPORT void JNICALL Java_xmipp_jni_ProgTomographAlignment_run
 (JNIEnv * env, jobject jobj){
	std::string msg = "";
  	try{
  	  ProgTomographAlignment *program = GET_INTERNAL_PTA(jobj);
  	  program->run();
  		} catch (XmippError xe) {
			msg = xe.getDefaultMessage();
		} catch (std::exception& e) {
			msg = e.what();
		} catch (...) {
			msg = "Unhandled exception";
	}

	// If there was an exception, sends it to java environment.
	if(!msg.empty()) {
		handleXmippException(env, msg);
	}
}

