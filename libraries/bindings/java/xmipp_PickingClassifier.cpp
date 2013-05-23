#include "xmipp_java_adapter.h"

#include <iostream>
#include "xmipp_PickingClassifier.h"
#include "xmipp_ExceptionsHandler.h"



JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_autopick
(JNIEnv *env, jclass class_, jstring filename)
{
    XMIPP_JAVA_TRY
    {
        std::cout << "autopick "<< filename << std::endl;
    }
    XMIPP_JAVA_CATCH;

}

JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_correct
(JNIEnv *env, jclass class_, jstring filename)
{
    XMIPP_JAVA_TRY
    {
    	std::cout << "correct "<< filename << std::endl;
    }
    XMIPP_JAVA_CATCH;

}

JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_train
(JNIEnv *env, jclass class_)
{
    XMIPP_JAVA_TRY
    {
    	std::cout << "train " << std::endl;
    }
    XMIPP_JAVA_CATCH;

}


