#include "xmipp_java_adapter.h"

#include <iostream>
#include "xmipp_Filename.h"
#include "xmipp_ExceptionsHandler.h"
#include "data/xmipp_filename.h"


JNIEXPORT jboolean JNICALL Java_xmipp_jni_PickingClassifier_autopick
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

JNIEXPORT jboolean JNICALL Java_xmipp_jni_PickingClassifier_correct
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


