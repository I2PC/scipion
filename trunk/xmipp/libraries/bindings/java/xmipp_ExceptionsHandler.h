#ifndef XMIPP_EXCEPTIONS_HANDLER
#define XMIPP_EXCEPTIONS_HANDLER
#include <jni.h>
#include "data/xmipp_error.h"
#include "data/xmipp_strings.h"

void handleXmippException(JNIEnv *env, std::string message);

//This is a macro for common error handling inside JNI implementations
#define XMIPP_TRY String __msg(""); try
#define XMIPP_CATCH catch (XmippError &xe){\
    __msg = xe.getMessage();}\
catch (std::exception& e){\
  __msg = e.what();}\
catch (...){\
  __msg = "Unhandled exception";}\
if(!__msg.empty()){\
    handleXmippException(env, __msg);}\

#endif
