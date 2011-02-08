#include "xmipp_ExceptionsHandler.h"
//#include <jni.h>
//#include <string>

using namespace std;

void handleXmippException(JNIEnv *env, std::string message) {
	jclass newExcCls = env->FindClass("java/lang/Exception");
	if (newExcCls != 0) {
		env->ThrowNew(newExcCls, message.c_str());
	}
}
