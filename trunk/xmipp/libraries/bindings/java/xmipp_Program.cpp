#include "xmipp_java_adapter.h"

#include <jni.h>
#include "xmipp_Program.h"
#include "xmipp_ExceptionsHandler.h"
#include <data/xmipp_filename.h>
#include <reconstruction/program_extension.h>

JNIEXPORT jint JNICALL Java_xmipp_jni_Program_runByName
(JNIEnv *env, jclass jc, jstring progName, jstring args)
{
    String msg = "";

    try
    {
        String progNameStr = env->GetStringUTFChars(progName, false);
        String argsStr = env->GetStringUTFChars(args, false);

        return runProgram(progNameStr,argsStr);
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if(!msg.empty())
    {
        handleXmippException(env, msg);
    }
}

JNIEXPORT jstring JNICALL Java_xmipp_jni_Program_getXmippPath
(JNIEnv *env, jclass jc)
{
    String msg = "";

    try
    {
        return env->NewStringUTF(FileName::getXmippPath().data());
    }
    catch (XmippError xe)
    {
        msg = xe.getDefaultMessage();
    }
    catch (std::exception& e)
    {
        msg = e.what();
    }
    catch (...)
    {
        msg = "Unhandled exception";
    }

    // If there was an exception, sends it to java environment.
    if(!msg.empty())
    {
        handleXmippException(env, msg);
    }
}
