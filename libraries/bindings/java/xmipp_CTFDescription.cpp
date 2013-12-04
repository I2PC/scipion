#include "xmipp_java_adapter.h"
#include <iostream>
#include "xmipp_CTFDescription.h"
#include "xmipp_ExceptionsHandler.h"
#include "xmipp_InternalData.h"
#include <data/ctf.h>

/*JNIEXPORT void JNICALL Java_xmipp_jni_CTFDescription_storeIds
(JNIEnv *env, jclass cls) {
 CTFDescription_peerId = env->GetFieldID(cls, "peer", "J");
}*/

JNIEXPORT void JNICALL Java_xmipp_jni_CTFDescription_create
(JNIEnv *env, jobject jobj)
{
    CTFDescription *ctfDescription = new CTFDescription();
    //env->SetLongField(jobj, CTFDescription_peerId, (long)ctfDescription);
    STORE_PEER_ID(jobj, ctfDescription);
}

JNIEXPORT void JNICALL Java_xmipp_jni_CTFDescription_destroy
(JNIEnv *env, jobject jobj)
{
    CTFDescription *ctfDescription = GET_INTERNAL_CTFDESCRIPTION(jobj);
    delete ctfDescription;
    ctfDescription = NULL;
    //env->SetLongField(jobj, CTFDescription_peerId, (long)ctfDescription);
    STORE_PEER_ID(jobj, ctfDescription);
}

JNIEXPORT void JNICALL Java_xmipp_jni_CTFDescription_read_1
(JNIEnv *env, jobject jobj, jstring filename)
{
    std::string msg = "";
    CTFDescription *ctfDescription = GET_INTERNAL_CTFDESCRIPTION(jobj);

    if (ctfDescription != NULL)
    {
        try
        {
        	jboolean aux=false;
            const char *fnStr = env->GetStringUTFChars(filename, &aux);

            ctfDescription->enable_CTF = ctfDescription->enable_CTFnoise = true;
            ctfDescription->read(fnStr);
            ctfDescription->produceSideInfo();
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
    }
    else
    {
        msg = "CTFDescription is null";
    }

    // If there was an exception, sends it to java environment.
    if(!msg.empty())
    {
        handleXmippException(env, msg);
    }
}

JNIEXPORT jdouble JNICALL Java_xmipp_jni_CTFDescription_getFMAX(JNIEnv *env,
        jobject jobj)
{
    std::string msg = "";
    CTFDescription *ctfDescription = GET_INTERNAL_CTFDESCRIPTION(jobj);

    double FMAX=0.;

    if (ctfDescription != NULL)
    {
        try
        {
            FMAX = 1 / (2 * ctfDescription->Tm);
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
    }
    else
    {
        msg = "CTFDescription is null";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return FMAX;
}

JNIEXPORT jobjectArray JNICALL Java_xmipp_jni_CTFDescription_CTFProfile(
    JNIEnv *env, jobject jobj, jdouble angle, jdouble FMAX, jint samples)
{
    std::string msg = "";
    CTFDescription *ctfDescription = GET_INTERNAL_CTFDESCRIPTION(jobj);

    MultidimArray<double> profiles;

    if (ctfDescription != NULL)
    {
        try
        {
            ctfDescription->getProfile(angle, FMAX, samples, profiles);

            // Stores result to return data.
            jdoubleArray row = env->NewDoubleArray(1);
            size_t nprofiles = XSIZE(profiles);
            jobjectArray profilesArray = env->NewObjectArray(nprofiles,
                                         env->GetObjectClass(row), 0);

            double *aux=new double[samples];
            for (size_t i = 0; i < nprofiles; i++)
            {
                for (int j = 0; j < samples; j++)
                {
                    aux[j] = A2D_ELEM(profiles, j, i);
                }

                // Sets intermediate data.
                row = env->NewDoubleArray(samples);
                env->SetDoubleArrayRegion((jdoubleArray) row, 0, samples,
                                          (jdouble *) aux);

                // Assings the row to the result array object.

                env->SetObjectArrayElement(profilesArray, i, row);

            }
            delete []aux;
            return profilesArray;

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
    }
    else
    {
        msg = "CTFDescription is null";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return (jobjectArray) NULL;
}

JNIEXPORT jobjectArray JNICALL Java_xmipp_jni_CTFDescription_CTFAverageProfile(
    JNIEnv *env, jobject jobj, jdouble FMAX, jint samples)
{
    std::string msg = "";
    CTFDescription *ctfDescription = GET_INTERNAL_CTFDESCRIPTION(jobj);

    MultidimArray<double> profiles;

    if (ctfDescription != NULL)
    {
        try
        {
            ctfDescription->getAverageProfile(FMAX, samples, profiles);

            // Stores result to return data.
            jdoubleArray row = env->NewDoubleArray(1);
            size_t nprofiles = XSIZE(profiles);
            jobjectArray profilesArray = env->NewObjectArray(nprofiles,
                                         env->GetObjectClass(row), 0);

            double *aux=new double[samples];
            for (size_t i = 0; i < nprofiles; i++)
            {
                for (int j = 0; j < samples; j++)
                {
                    aux[j] = A2D_ELEM(profiles, j, i);
                }

                // Sets intermediate data.
                row = env->NewDoubleArray(samples);
                env->SetDoubleArrayRegion((jdoubleArray) row, 0, samples,
                                          (jdouble *) aux);

                // Assings the row to the result array object.
                env->SetObjectArrayElement(profilesArray, i, row);
            }
            delete []aux;

            return profilesArray;
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
    }
    else
    {
        msg = "CTFDescription is null";
    }

    // If there was an exception, sends it to java environment.
    if (!msg.empty())
    {
        handleXmippException(env, msg);
    }

    return (jobjectArray) NULL;
}
