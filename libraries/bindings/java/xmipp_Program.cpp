#include "xmipp_java_adapter.h"

#include <jni.h>
#include "xmipp_Program.h"
#include "xmipp_ExceptionsHandler.h"
#include <data/xmipp_filename.h>
#include <reconstruction/program_extension.h>

JNIEXPORT jint JNICALL Java_xmipp_jni_Program_runByName
(JNIEnv *env, jclass jc, jstring progName, jstring args)
{
     XMIPP_JAVA_TRY
     {
    	jboolean aux=false;
        String progNameStr = env->GetStringUTFChars(progName, &aux);
        aux=false;
        String argsStr = env->GetStringUTFChars(args, &aux);

        return runProgram(progNameStr,argsStr);
     }
     XMIPP_JAVA_CATCH;
     return 1;
}

