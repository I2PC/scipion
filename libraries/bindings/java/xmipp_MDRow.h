
#include <jni.h>
/* Header for class xmipp_MDRow */

#ifndef _Included_xmipp_jni_MDRow
#define _Included_xmipp_jni_MDRow
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Class:     xmipp_MDRow
 * Method:    create
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_create
  (JNIEnv *, jobject);

/*
 * Class:     xmipp_MDRow
 * Method:    destroy
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_destroy
  (JNIEnv *, jobject);

/*
 * Class:     xmipp_MDRow
 * Method:    clear
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_clear
  (JNIEnv *, jobject);

/*
 * Class:     xmipp_MDRow
 * Method:    size
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_xmipp_jni_MDRow_size
  (JNIEnv *, jobject);

/*
 * Class:     xmipp_MDRow
 * Method:    containsLabel
 * Signature: (V)Z
 */
JNIEXPORT jboolean JNICALL Java_xmipp_jni_MDRow_emtpy
  (JNIEnv *, jobject);

/*
 * Class:     xmipp_MDRow
 * Method:    containsLabel
 * Signature: (I)Z
 */
JNIEXPORT jboolean JNICALL Java_xmipp_jni_MDRow_containsLabel
  (JNIEnv *, jobject, jint);

///*
// * Class:     xmipp_MDRow
// * Method:    containsLabel
// * Signature: (I)Z
// */
//JNIEXPORT jboolean JNICALL Java_xmipp_jni_MDRow_addLabel
//  (JNIEnv *, jobject, jint);

/*
 * Class:     xmipp_MDRow
 * Method:    getValueInt
 * Signature: (IJ)I
 */
JNIEXPORT jint JNICALL Java_xmipp_jni_MDRow_getValueInt
  (JNIEnv *, jobject, jint);

/*
 * Class:     xmipp_MDRow
 * Method:    getValueLong
 * Signature: (IJ)J
 */
JNIEXPORT jlong JNICALL Java_xmipp_jni_MDRow_getValueLong
  (JNIEnv *, jobject, jint);

/*
 * Class:     xmipp_MDRow
 * Method:    getValueDouble
 * Signature: (IJ)D
 */
JNIEXPORT jdouble JNICALL Java_xmipp_jni_MDRow_getValueDouble
  (JNIEnv *, jobject, jint);

/*
 * Class:     xmipp_MDRow
 * Method:    getValueString
 * Signature: (IJ)Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_xmipp_jni_MDRow_getValueString
  (JNIEnv *, jobject, jint);

/*
 * Class:     xmipp_MDRow
 * Method:    getValueBoolean
 * Signature: (IJ)Z
 */
JNIEXPORT jboolean JNICALL Java_xmipp_jni_MDRow_getValueBoolean
  (JNIEnv *, jobject, jint);

/*
 * Class:     xmipp_MDRow
 * Method:    setValueInt
 * Signature: (IIJ)Z
 */
JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_setValueInt
  (JNIEnv *, jobject, jint, jint);

/*
 * Class:     xmipp_MDRow
 * Method:    setValueLong
 * Signature: (IIJ)Z
 */
JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_setValueLong
  (JNIEnv *, jobject, jint, jlong);

/*
 * Class:     xmipp_MDRow
 * Method:    setValueDouble
 * Signature: (IDJ)Z
 */
JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_setValueDouble
  (JNIEnv *, jobject, jint, jdouble);

/*
 * Class:     xmipp_MDRow
 * Method:    setValueString
 * Signature: (ILjava/lang/String;J)Z
 */
JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_setValueString
  (JNIEnv *, jobject, jint, jstring);

/*
 * Class:     xmipp_MDRow
 * Method:    setValueBoolean
 * Signature: (IZJ)Z
 */
JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_setValueBoolean
  (JNIEnv *, jobject, jint, jboolean);


#ifdef __cplusplus
}
#endif
#endif
