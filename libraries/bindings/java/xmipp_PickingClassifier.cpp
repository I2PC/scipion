#include "xmipp_java_adapter.h"
#include <iostream>
#include "xmipp_PickingClassifier.h"
#include "xmipp_ExceptionsHandler.h"



JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_autopick
(JNIEnv *env, jclass class_, jstring filename)
{
    XMIPP_JAVA_TRY
    {
    	jboolean iscopy = false;
        std::cout << "autopick "<< env->GetStringUTFChars(filename, &iscopy) << std::endl;
    }
    XMIPP_JAVA_CATCH;

}

JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_correct
(JNIEnv *env, jclass class_, jstring filename)
{
    XMIPP_JAVA_TRY
    {
    	jboolean iscopy = false;
    	std::cout << "correct "<< env->GetStringUTFChars(filename, &iscopy) << std::endl;
    }
    XMIPP_JAVA_CATCH;

}

JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_train
(JNIEnv *env, jclass class_, jobjectArray micrographs)
{
    XMIPP_JAVA_TRY
    {
    	std::cout << "train " << std::endl;
    	int size = env->GetArrayLength(micrographs);
    	const char* name;
    	jboolean iscopy = false;
    	jstring micrograph;
    	for (int i = 0; i< size; i++)
    	{
			micrograph = (jstring) env->GetObjectArrayElement(micrographs, i);
			name = env->GetStringUTFChars(micrograph, &iscopy);
			std::cout << name << std::endl;
    	}
    }
    XMIPP_JAVA_CATCH;

}


