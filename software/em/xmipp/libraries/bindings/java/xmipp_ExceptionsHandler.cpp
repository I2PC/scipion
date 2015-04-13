#include "xmipp_java_adapter.h"
#include "xmipp_ExceptionsHandler.h"

void handleXmippException(JNIEnv *env, std::string message) {
	jclass newExcCls = env->FindClass("java/lang/Exception");
	if (newExcCls != 0) {
		env->ThrowNew(newExcCls, message.c_str());
	}
}
