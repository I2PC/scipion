#include <jni.h>
#include "xmipp_ImageDouble.h"
#include <data/image.h>

static jfieldID ImageDouble_peerId;
static jfieldID ImageDouble_filenameId;

#define peerId ImageDouble_peerId
#define GET_INTERNAL_IMAGE() ((Image<double>*)(env->GetLongField(obj, peerId)))

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_storeIds
(JNIEnv *env, jclass cls)
{
    peerId = env->GetFieldID(cls, "peer", "J");
}

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_create
(JNIEnv *env, jobject obj)
{
    Image<double> * image = new Image<double>();
    env->SetLongField(obj, peerId, (long)image);
}

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_destroy
(JNIEnv *env, jobject obj)
{
    Image<double> * image = GET_INTERNAL_IMAGE();
    delete image;
    image = NULL;
    env->SetLongField(obj, peerId, (long)image);
}

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_read
(JNIEnv *env, jobject obj, jstring filename)
{
    Image<double> * image = GET_INTERNAL_IMAGE();
    if (image != NULL)
    {
        const char * fnStr = env->GetStringUTFChars(filename, false);
        image->read(fnStr);
    }
}

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_write
(JNIEnv *env, jobject obj, jstring filename)
{
    Image<double> * image = GET_INTERNAL_IMAGE();

    if (image != NULL)
    {
        const char * fnStr = env->GetStringUTFChars(filename, false);
        image->write(fnStr);
    }
}

JNIEXPORT jdoubleArray JNICALL Java_xmipp_ImageDouble_getData
(JNIEnv *env, jobject obj)
{
    Image<double> * image = GET_INTERNAL_IMAGE();
    if (image != NULL)
    {
      unsigned long int size = image->getSize();
        jdoubleArray array = env->NewDoubleArray(size);
        env->SetDoubleArrayRegion(array, 0, size, MULTIDIM_ARRAY(image->data));
        return array;
    }
}

JNIEXPORT jintArray JNICALL Java_xmipp_ImageDouble_getSize
(JNIEnv *env, jobject obj)
{
    jintArray array = env->NewIntArray(4);
    Image<double> * image = GET_INTERNAL_IMAGE();
    if (image != NULL)
    {
        unsigned long n;
        int xyz[3];
        image->getDimensions(xyz[0], xyz[1], xyz[2], n);
        env->SetIntArrayRegion(array, 0, 3, xyz);
        return array;
    }
    return NULL;
}

/*
 * Class:     ImageDouble
 * Method:    printShape
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_xmipp_ImageDouble_printShape
(JNIEnv *env, jobject obj)
{
    Image<double> * image = GET_INTERNAL_IMAGE();
    std::cout << (*image) << std::endl;
}


