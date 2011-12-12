#include <iostream>
#include "xmipp_Filename.h"
#include "xmipp_ExceptionsHandler.h"
#include <data/xmipp_filename.h>

JNIEXPORT jboolean JNICALL Java_xmipp_Filename_hasStackExtension
(JNIEnv *env, jclass class_, jstring filename)
{
    std::string msg = "";

    try
    {
        const char *fnStr = env->GetStringUTFChars(filename, false);

        return FileName(fnStr).hasStackExtension();
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

    return false;
}

JNIEXPORT jboolean JNICALL Java_xmipp_Filename_hasVolumeExtension
(JNIEnv *env, jclass class_, jstring filename)
{
    std::string msg = "";

    try
    {
        const char *fnStr = env->GetStringUTFChars(filename, false);

        return FileName(fnStr).hasVolumeExtension();
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

    return false;
}
