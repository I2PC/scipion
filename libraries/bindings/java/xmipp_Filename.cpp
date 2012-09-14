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
        const char *fnStr = env->GetStringUTFChars(filename, false);
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
        const char *fnStr = env->GetStringUTFChars(filename, false);
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
      const char *fnStr = env->GetStringUTFChars(filename, false);
      return FileName(fnStr).isMetaData(false);
  }
  XMIPP_JAVA_CATCH;

  return false;
}

JNIEXPORT jstring JNICALL Java_xmipp_jni_Filename_compose
  (JNIEnv *env, jclass class_, jint slice, jstring path){
    XMIPP_JAVA_TRY
    {
        const char *fnStr = env->GetStringUTFChars(path, false);
        FileName fnTemp;
        String filestring= String(fnStr);
        fnTemp.compose(slice,filestring);
        return env->NewStringUTF(fnTemp.data());;
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT jstring JNICALL Java_xmipp_jni_Filename_getXmippPath
  (JNIEnv *env, jclass class_){
    XMIPP_JAVA_TRY
    {
        return env->NewStringUTF(getXmippPath());
    }
    XMIPP_JAVA_CATCH;
}
