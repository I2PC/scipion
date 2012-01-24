#include "xmipp_java_adapter.h"

#include <iostream>
#include "xmipp_Filename.h"
#include "xmipp_ExceptionsHandler.h"
#include "data/xmipp_filename.h"


JNIEXPORT jboolean JNICALL Java_xmipp_jni_Filename_hasStackExtension
(JNIEnv *env, jclass class_, jstring filename)
{
    XMIPP_TRY
    {
        const char *fnStr = env->GetStringUTFChars(filename, false);
        return FileName(fnStr).hasStackExtension();
    }
    XMIPP_CATCH;

    return false;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_Filename_hasVolumeExtension
(JNIEnv *env, jclass class_, jstring filename)
{
    XMIPP_TRY
    {
        const char *fnStr = env->GetStringUTFChars(filename, false);
        return FileName(fnStr).hasVolumeExtension();
    }
    XMIPP_CATCH;

    return false;
}

JNIEXPORT jstring JNICALL Java_xmipp_jni_Filename_compose
  (JNIEnv *env, jclass class_, jint slice, jstring path){
    XMIPP_TRY
    {
        const char *fnStr = env->GetStringUTFChars(path, false);
        FileName fnTemp;
        std::string filestring= std::string(fnStr);
        fnTemp.compose(slice,filestring);
        return env->NewStringUTF(fnTemp.data());;
    }
    XMIPP_CATCH;
}

JNIEXPORT jstring JNICALL Java_xmipp_jni_Filename_getXmippPath
  (JNIEnv *env, jclass class_){
    XMIPP_TRY
    {
        return env->NewStringUTF(FileName::getXmippPath().c_str());
    }
    XMIPP_CATCH;
}
