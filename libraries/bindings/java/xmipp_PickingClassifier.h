/* DO NOT EDIT THIS FILE - it is machine generated */
#include <jni.h>
/* Header for class xmipp_Filename */

#ifndef _Included_xmipp_jni_PickingClassifier
#define _Included_xmipp_jni_PickingClassifier
#ifdef __cplusplus
extern "C" {
#endif


JNIEXPORT void JNICALL
Java_xmipp_jni_PickingClassifier_create(JNIEnv *env, jobject jobj, jint, jstring);



JNIEXPORT void JNICALL
Java_xmipp_jni_PickingClassifier_destroy(JNIEnv *env, jobject jobj);


JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_autopick
  (JNIEnv *, jobject, jstring, jobject, jint percent);


JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_correct
  (JNIEnv *, jobject, jobject, jobject);

JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_train
  (JNIEnv *, jobject, jobject);


#ifdef __cplusplus
}
#endif
#endif
