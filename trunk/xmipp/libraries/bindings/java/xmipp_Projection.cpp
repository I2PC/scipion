#include <jni.h>
#include "xmipp_Projection.h"
#include "xmipp_ImageGeneric.h"
#include "xmipp_ExceptionsHandler.h"
#include "xmipp_InternalData.h"
#include <data/xmipp_image.h>
#include <data/projection.h>
#include <data/filters.h>
#include <data/geometry.h>
/*
JNIEXPORT void JNICALL Java_xmipp_Projection_projectVolume
(JNIEnv *env, jclass class_, jobject jvolume, jobject jprojection, jdouble rot, jdouble tilt, jdouble psi)
{
    std::string msg = "";
    ImageGeneric *volume = GET_INTERNAL_IMAGE_GENERIC(jvolume);
    ImageGeneric *projection=  GET_INTERNAL_IMAGE_GENERIC(jprojection);

    if(volume != NULL)
    {
        if(projection != NULL)
        {
            try
            {
                int w = XSIZE(* projection->data->im);
                int h = YSIZE(* projection->data->im);

                projection->convert2Datatype(Double);
                volume->convert2Datatype(Double);

                MultidimArray<double> *mdVolume;
                MULTIDIM_ARRAY_GENERIC(*volume).getMultidimArrayPointer(mdVolume);

                Projection auxProjection;
                MultidimArray<double> *mdProjection;

                MULTIDIM_ARRAY_GENERIC(*projection).getMultidimArrayPointer(mdProjection);
                MULTIDIM_ARRAY(auxProjection).alias(*mdProjection);

                projectVolume(*mdVolume, auxProjection, h, w, (double) rot, (double) tilt, (double) psi);
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
            msg = "projection is NULL";
        }
    }
    else
    {
        msg = "volume is NULL";
    }

    // If there was an exception, sends it to java environment.
    if(!msg.empty())
    {
        handleXmippException(env, msg);
    }
}*/

JNIEXPORT void JNICALL Java_xmipp_Projection_projectVolume
(JNIEnv *env, jclass class_, jobject jvolume, jobject jprojection, jdoubleArray matrix)
{
    std::string msg = "";
    ImageGeneric *volume = GET_INTERNAL_IMAGE_GENERIC(jvolume);
    ImageGeneric *projection=  GET_INTERNAL_IMAGE_GENERIC(jprojection);

    if(volume != NULL)
    {
        if(projection != NULL)
        {
            try
            {
                int w = XSIZE(* projection->data->im);
                int h = YSIZE(* projection->data->im);

                projection->convert2Datatype(Double);
                volume->convert2Datatype(Double);

                MultidimArray<double> *mdVolume;
                MULTIDIM_ARRAY_GENERIC(*volume).getMultidimArrayPointer(mdVolume);

                Projection auxProjection;
                MultidimArray<double> *mdProjection;

                MULTIDIM_ARRAY_GENERIC(*projection).getMultidimArrayPointer(mdProjection);
                MULTIDIM_ARRAY(auxProjection).alias(*mdProjection);

                MultidimArray<double> volumeRot;

                Matrix2D<double> m(4, 4);
                env->GetDoubleArrayRegion(matrix, 0, 16, (jdouble *)MATRIX2D_ARRAY(m));

                // Transforms...
                mdVolume->setXmippOrigin();
                applyGeometry(LINEAR, volumeRot, *mdVolume, m, false, false);

                // Projects...
                projectVolume(volumeRot, auxProjection, h, w, 0, 0, 0);

                Image<double> imageAux;
                imageAux().alias(volumeRot);
                imageAux.write("/home/jvega/volume_tmp.vol");
                //projectVolume(*mdVolume, auxProjection, h, w, (double) rot, (double) tilt, (double) psi);
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
            msg = "projection is NULL";
        }
    }
    else
    {
        msg = "volume is NULL";
    }

    // If there was an exception, sends it to java environment.
    if(!msg.empty())
    {
        handleXmippException(env, msg);
    }
}

JNIEXPORT jdouble JNICALL Java_xmipp_Projection_entropyOtsuSegmentation
(JNIEnv *env, jclass class_, jobject jvolume, jdouble percentile, jboolean binarize)
{
    std::string msg = "";
    ImageGeneric *volume = GET_INTERNAL_IMAGE_GENERIC(jvolume);

    try
    {
        MultidimArray<double> *mdarray;
        volume->convert2Datatype(Double);

        MULTIDIM_ARRAY_GENERIC(*volume).getMultidimArrayPointer(mdarray);

        return EntropyOtsuSegmentation(*mdarray, percentile, binarize);
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

    return -1;
}
/*
JNIEXPORT jdoubleArray JNICALL Java_xmipp_Projection_eulerDirection2Angles
(JNIEnv *env, jclass class_, jdoubleArray jvector)
{
    std::string msg = "";

    try
    {
        Matrix1D<double> v(3);

        env->GetDoubleArrayRegion(jvector, 0, 3, (jdouble *)MATRIX1D_ARRAY(v));

        double alpha, beta, gamma;

        Euler_direction2angles(v, alpha, beta, gamma);

        double angles[3];
        angles[0] = alpha;
        angles[1] = beta;
        angles[2] = gamma;

        // Sets array value
        jdoubleArray array = env->NewDoubleArray(3);
        env->SetDoubleArrayRegion(array, 0, 3, angles);

        return array;
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

    return NULL;
}
*/
JNIEXPORT jdoubleArray JNICALL Java_xmipp_Projection_eulerMatrix2Angles
(JNIEnv *env, jclass class_, jdoubleArray matrix)
{
    std::string msg = "";

    try
    {
        Matrix2D<double> m(4, 4);

        env->GetDoubleArrayRegion(matrix, 0, 16, (jdouble *)MATRIX2D_ARRAY(m));

        m.resize(3,3);

        MAT_ELEM(m,1,1) *= -1;
        MAT_ELEM(m,2,2) *= -1;
        m=m.inv();

        double alpha, beta, gamma;
        Euler_matrix2angles(m, alpha, beta, gamma);

        double angles[3];
        angles[0] = alpha;
        angles[1] = beta;
        angles[2] = gamma;

        // Sets array value
        jdoubleArray array = env->NewDoubleArray(3);
        env->SetDoubleArrayRegion(array, 0, 3, angles);

        return array;
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

    return NULL;
}
