#ifndef XMIPP_EXCEPTIONS_HANDLER
#define XMIPP_EXCEPTIONS_HANDLER
#include <jni.h>
#include <string>

void handleXmippException(JNIEnv *env, std::string message);

#endif
