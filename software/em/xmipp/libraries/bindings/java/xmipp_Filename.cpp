#include "xmipp_java_adapter.h"

#include <iostream>
#include "xmipp_Filename.h"
#include "xmipp_ExceptionsHandler.h"
#include "data/xmipp_filename.h"


JNIEXPORT jboolean JNICALL Java_xmipp_jni_Filename_hasStackExtension
(JNIEnv *env, jclass class_, jstring filename)
{
    XMIPP_JAVA_TRY
    {
        jboolean aux=false;
        const char *fnStr = env->GetStringUTFChars(filename, &aux);
        return FileName(fnStr).hasStackExtension();
    }
    XMIPP_JAVA_CATCH;

    return false;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_Filename_hasVolumeExtension
(JNIEnv *env, jclass class_, jstring filename)
{
    XMIPP_JAVA_TRY
    {
        jboolean aux=false;
        const char *fnStr = env->GetStringUTFChars(filename, &aux);
        return FileName(fnStr).hasVolumeExtension();
    }
    XMIPP_JAVA_CATCH;

    return false;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_Filename_isMetaDataFile
(JNIEnv *env, jclass class_, jstring filename)
{
  XMIPP_JAVA_TRY
  {
      jboolean aux=false;
      const char *fnStr = env->GetStringUTFChars(filename, &aux);
      return FileName(fnStr).isMetaData(false);
  }
  XMIPP_JAVA_CATCH;

  return false;
}

JNIEXPORT jstring JNICALL Java_xmipp_jni_Filename_compose
  (JNIEnv *env, jclass class_, jint slice, jstring path){
    XMIPP_JAVA_TRY
    {
        jboolean aux=false;
        const char *fnStr = env->GetStringUTFChars(path, &aux);
        FileName fnTemp;
        String filestring= String(fnStr);
        fnTemp.compose(slice,filestring);
        return env->NewStringUTF(fnTemp.c_str());;
    }
    XMIPP_JAVA_CATCH;
    return env->NewStringUTF("");
}

JNIEXPORT jstring JNICALL Java_xmipp_jni_Filename_getXmippPath
  (JNIEnv *env, jclass class_){
    XMIPP_JAVA_TRY
    {
        return env->NewStringUTF(getXmippPath());
    }
    XMIPP_JAVA_CATCH;
    return env->NewStringUTF("");
}
