#include "xmipp_java_adapter.h"

#include <jni.h>
#include "xmipp_ImageGeneric.h"
#include "xmipp_InternalData.h"
#include "xmipp_ExceptionsHandler.h"
#include "data/xmipp_image_generic.h"
#include "data/xmipp_fft.h"
#include "data/ctf.h"
#include "reconstruction/transform_downsample.h"
#include "data/filters.h"
#include "data/geometry.h"

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_create(JNIEnv *env, jobject jobj)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = new ImageGeneric();
        STORE_PEER_ID(jobj, (long)image);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_destroy(JNIEnv *env, jobject jobj)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        delete image;
        image = NULL;
        STORE_PEER_ID(jobj, (long)image);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_resize(JNIEnv *env, jobject jobj, jint w, jint h,
                                   jint d, jlong n)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        image->resize(w, h, d, (size_t) n, false);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT jint JNICALL
Java_xmipp_jni_ImageGeneric_getXDim(JNIEnv *env, jobject jobj)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        return XSIZE(* image->data->im);
    }
    XMIPP_JAVA_CATCH;

    return -1;
}

JNIEXPORT jint JNICALL
Java_xmipp_jni_ImageGeneric_getYDim(JNIEnv *env, jobject jobj)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

        return YSIZE(* image->data->im);
    }
    XMIPP_JAVA_CATCH;

    return -1;
}

JNIEXPORT jint JNICALL
Java_xmipp_jni_ImageGeneric_getZDim(JNIEnv *env, jobject jobj)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

        return ZSIZE(* image->data->im);
    }
    XMIPP_JAVA_CATCH;

    return -1;
}

JNIEXPORT jlong JNICALL
Java_xmipp_jni_ImageGeneric_getNDim(JNIEnv *env, jobject jobj)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

        return NSIZE(* image->data->im);
    }
    XMIPP_JAVA_CATCH;

    return -1;
}

JNIEXPORT jint JNICALL
Java_xmipp_jni_ImageGeneric_getDataType(JNIEnv *env, jobject jobj)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

        return image->getDatatype();
    }
    XMIPP_JAVA_CATCH;

    return -1;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_readHeader(JNIEnv * env, jobject jobj,
                                       jstring filename)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        jboolean aux=false;
        const char *fnStr = env->GetStringUTFChars(filename, &aux);
        image->read(fnStr, HEADER);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_read(JNIEnv *env, jobject jobj, jstring filename,
                                 jint jx, jint jy, jint jz, jlong jn)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        jboolean aux=false;
        const char *fn = env->GetStringUTFChars(filename, &aux);
        image->readOrReadPreview(fn, jx, jy, jz, jn, true);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_readApplyGeo_1(JNIEnv *env, jobject jimage,
        jstring filename, jobject jmetadata, jlong id, jint w, jint h, jboolean wrap)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jimage);
        MetaData *metadata = GET_INTERNAL_METADATA(jmetadata);

        jboolean aux=false;
        const char *fnStr = env->GetStringUTFChars(filename, &aux);
        ApplyGeoParams params;
        params.wrap = wrap;
        image->readApplyGeo(fnStr, *metadata, (size_t) id, params);

#define SELF_SCALE_TO_SIZE(type) MultidimArray<type> *mdaFloat; \
 MULTIDIM_ARRAY_GENERIC(*image).getMultidimArrayPointer(mdaFloat); \
 selfScaleToSize(LINEAR, *mdaFloat, (int)w, (int)h);

        SWITCHDATATYPE(image->getDatatype(), SELF_SCALE_TO_SIZE);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT jbyteArray JNICALL
Java_xmipp_jni_ImageGeneric_getArrayByte(JNIEnv *env, jobject jobj, jlong select_image, jint nslice)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

        // Go to slice.
        image->movePointerTo(nslice, select_image);

        size_t size = image->getSize();
        jbyteArray array = env->NewByteArray(size);

        DataType dataType = image->getDatatype();

        switch (dataType)
    {
    case DT_UChar:
        {
            unsigned char *mdarray;

            // Get slice array.
            MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);

                env->SetByteArrayRegion(array, 0, size, (jbyte *) mdarray);
            }
            break;
        case DT_SChar:
            {
                char *mdarray;

                // Get slice array.
                MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);

                unsigned char * buffer = new unsigned char[rw_max_page_size];
                size_t page_size = rw_max_page_size;

                for (size_t written = 0; written < size; written += page_size)
                {
                    page_size = std::min(page_size, size - written);
                    unsigned char * iter = (unsigned char *) mdarray + written;

                    for (size_t k = 0; k < page_size; ++k)
                        buffer[k] = iter[k] - CHAR_MIN;

                    env->SetByteArrayRegion(array, written, page_size, (jbyte *) buffer);
                }
                delete [] buffer;
            }
            break;
        default:
                REPORT_ERROR(
                    ERR_ARG_INCORRECT,
                    (String)"Java_xmipp_jni_ImageGeneric_getArrayByte: reading invalid image datatype: " + datatype2Str((DataType)dataType));
        }

        // Resets slice pointer.
        image->movePointerTo();

        return array;
    }
    XMIPP_JAVA_CATCH;

    return NULL;
}

JNIEXPORT jshortArray JNICALL
Java_xmipp_jni_ImageGeneric_getArrayShort(JNIEnv *env, jobject jobj,
        jlong select_image, jint nslice)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

        // Go to slice.
        image->movePointerTo(nslice, select_image);

        size_t size = image->getSize();
        jshortArray array = env->NewShortArray(size);

        DataType dataType = image->getDatatype();

        switch (dataType)
    {
    case DT_UShort:
        {
            unsigned short *mdarray;

            // Get slice array.
            MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);

                env->SetShortArrayRegion(array, 0, size, (jshort *) mdarray);
            }
            break;
        case DT_Short:
            {
                short *mdarray;

                // Get slice array.
                MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);

                unsigned short * buffer = new unsigned short[rw_max_page_size];
                size_t page_size = rw_max_page_size;

                for (size_t written = 0; written < size; written += page_size)
                {
                    page_size = std::min(page_size, size - written);
                    unsigned short * iter = (unsigned short *) (mdarray + written);

                    for (size_t k = 0; k < page_size; ++k)
                        buffer[k] = iter[k] - SHRT_MIN;

                    env->SetShortArrayRegion(array, written, page_size,
                                             (jshort *) buffer);
                }
                delete [] buffer;
            }
            break;
        default:
                REPORT_ERROR(
                    ERR_ARG_INCORRECT,
                    (String)"Java_xmipp_jni_ImageGeneric_getArrayShort: reading invalid image datatype: " + datatype2Str((DataType)dataType));
        }

        // Resets slice pointer.
        image->movePointerTo();

        return array;
    }
    XMIPP_JAVA_CATCH;

    return NULL;
}

JNIEXPORT jfloatArray JNICALL
Java_xmipp_jni_ImageGeneric_getArrayFloat(JNIEnv *env, jobject jobj,
        jlong select_image, jint nslice)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

        // Go to slice.
        image->movePointerTo(nslice, select_image);

        size_t size = image->getSize();
        jfloatArray array = env->NewFloatArray(size);

        DataType dataType = image->getDatatype();

        switch (dataType)
    {
    case DT_Float:
        {
            float *mdarray;

            // Get slice array.
            MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);

                env->SetFloatArrayRegion(array, 0, size, mdarray);
            }
            break;
        default:
            {
#define CAST_PAGE(type) Image<type> *imageAux = (Image<type> *)image->image; \
                type *data = MULTIDIM_ARRAY(imageAux->data); \
                size_t page_size = rw_max_page_size; \
                char * buffer = new char[rw_max_page_size * gettypesize(DT_Float)]; \
                for (size_t written = 0; written < size; written += page_size) \
                { \
                    page_size = std::min(page_size, size - written); \
                    imageAux->castPage2Datatype(data + written, buffer, DT_Float, page_size); \
                    env->SetFloatArrayRegion(array, written, page_size, (jfloat*) buffer); \
                } \
                delete[] buffer;

                SWITCHDATATYPE(dataType, CAST_PAGE);

            }
            break;
        }

        // Resets slice pointer.
        image->movePointerTo();

        return array;
    }
    XMIPP_JAVA_CATCH;

    return NULL;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_setArrayByte(JNIEnv *env, jobject jobj,
        jbyteArray data, jlong select_image, jint nslice)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        DataType dataType = image->getDatatype();

        // Go to slice.
        image->movePointerTo(nslice, select_image);

        size_t size = image->getSize();

        switch (dataType)
    {
    case DT_UChar:
        {
            unsigned char *mdarray;

            // Get slice array.
            MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);

                env->GetByteArrayRegion(data, 0, size, (jbyte *) mdarray);
            }
            break;
        case DT_SChar:
            {
                char *mdarray;

                // Get slice array.
                MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);

                unsigned char * buffer = new unsigned char[rw_max_page_size];
                size_t page_size = rw_max_page_size;

                for (size_t written = 0; written < size; written += page_size)
                {
                    page_size = std::min(page_size, size - written);

                    env->GetByteArrayRegion(data, written, page_size, (jbyte *) buffer);

                    char * iter = mdarray + written;
                    for (size_t k = 0; k < page_size; ++k)
                        iter[k] = (char) (buffer[k] + CHAR_MIN);
                }
                delete [] buffer;
            }
            break;
        default:
                {
#undef CAST_PAGE
#define CAST_PAGE(type) \
  type *mdarray; \
                MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray); \
                Image<type> *imageAux = (Image<type> *) image->image; \
                char * buffer = new  char[rw_max_page_size * sizeof(char)]; \
                size_t page_size = rw_max_page_size; \
                for (size_t written = 0; written < size; written += page_size) \
                { \
                    page_size = std::min(page_size, size - written); \
                    env->GetByteArrayRegion(data, written, page_size, (jbyte *) buffer); \
                    imageAux->castPage2T(written, buffer, DT_UChar, page_size); \
                } \
                SWITCHDATATYPE(image->getDatatype(), CAST_PAGE);\
                  delete [] buffer;

                }
                break;
        }

        // Resets slice pointer.
        image->movePointerTo();
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_setArrayShort(JNIEnv *env, jobject jobj,
        jshortArray data, jlong select_image, jint nslice)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

        // Go to slice.
        image->movePointerTo(nslice, select_image);

        size_t size = image->getSize();

        switch (image->getDatatype())
    {
    case DT_UShort:
        {
            unsigned short *mdarray;

            // Get slice array.
            MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);

                env->GetShortArrayRegion(data, 0, size, (jshort *) mdarray);
            }
            break;
        case DT_Short:
            {
                short *mdarray;

                // Get slice array.
                MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);

                unsigned short * buffer = new unsigned short[rw_max_page_size];
                size_t page_size = rw_max_page_size;

                for (size_t written = 0; written < size; written += page_size)
                {
                    page_size = std::min(page_size, size - written);

                    env->GetShortArrayRegion(data, written, page_size, (jshort *) buffer);

                    short * iter = mdarray + written;
                    for (size_t k = 0; k < page_size; ++k)
                        iter[k] = (short) (buffer[k] + SHRT_MIN);
                }
                delete [] buffer;
            }
            break;
        case DT_Float:
                {
                    float *mdarray;

                    // Get slice array.
                    MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);
                    Image<float> *imageAux = (Image<float> *) image->image;

                    char * buffer = new char[rw_max_page_size * sizeof(short)];
                    size_t page_size = rw_max_page_size;

                    for (size_t written = 0; written < size; written += page_size)
                {
                    page_size = std::min(page_size, size - written);
                        env->GetShortArrayRegion(data, written, page_size, (jshort *) buffer);
                        imageAux->setPage2T(written, buffer, DT_UShort, page_size);
                    }
                delete [] buffer;
            }
            break;
        default:
            {
                REPORT_ERROR(
                    ERR_IO_NOWRITE,
                    (String)"Not supported conversion to dataType: " + datatype2Str(image->getDatatype()));
            }
            break;
        }

        // Resets slice pointer.
        image->movePointerTo();
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_setArrayFloat(JNIEnv *env, jobject jobj,
        jfloatArray data, jlong select_image, jint nslice)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

        // Go to slice.
        image->movePointerTo(nslice, select_image);
        size_t size = image->getSize();

        //jfloatArray array = env->NewFloatArray(size);
        switch (image->getDatatype())
    {
    case DT_Float:
        	{
				float *mdarray;

				// Get slice array.
				MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);
				env->GetFloatArrayRegion(data, 0, size, mdarray);
			}
			break;
    case DT_Double:
            {
                // Get slice array.
                Image<double> *imageAux = (Image<double> *) image->image;

                char * buffer = new char[rw_max_page_size * sizeof(float)];
                size_t page_size = rw_max_page_size;

                for (size_t written = 0; written < size; written += page_size)
                {
                    page_size = std::min(page_size, size - written);
                    env->GetFloatArrayRegion(data, written, page_size, (jfloat *) buffer);
                    imageAux->setPage2T(written, buffer, DT_Float, page_size);
                }
                delete [] buffer;
            }
            break;
    default:
            {
                REPORT_ERROR(
                    ERR_IO_NOWRITE,
                    (String)"Not supported conversion From jfloat to dataType: " + datatype2Str(image->getDatatype()));
            }
            break;
        }

        // Resets slice pointer.
        image->movePointerTo();
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_setDataType(JNIEnv *env, jobject jobj,
                                        jint dataType)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        image->setDatatype((DataType) dataType);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_convert2Datatype(JNIEnv *env, jobject jobj,
        jint dataType)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        image->convert2Datatype((DataType) dataType);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_mapFile2Write(JNIEnv *env, jobject jobj, jint w,
        jint h, jint z, jstring filename, jlong n)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        jboolean aux=false;
        const char *fnStr = env->GetStringUTFChars(filename, &aux);

        image->mapFile2Write(w, h, z, fnStr, false, (size_t) n);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_write(JNIEnv *env, jobject jobj, jstring filename)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        jboolean aux=false;
        const char *fnStr = env->GetStringUTFChars(filename, &aux);
        image->write(fnStr);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_printShape(JNIEnv *env, jobject jobj)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        image->print();
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT jboolean JNICALL
Java_xmipp_jni_ImageGeneric_equal(JNIEnv *env, jobject jobj1, jobject jobj2, jdouble accuracy)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image1 = GET_INTERNAL_IMAGE_GENERIC(jobj1);
        ImageGeneric *image2 = GET_INTERNAL_IMAGE_GENERIC(jobj2);
        return image1->equal(*image2,accuracy);
    }
    XMIPP_JAVA_CATCH;
    return false;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_subtract(JNIEnv *env, jobject jobj1, jobject jobj2, jobject jobj3)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image1 = GET_INTERNAL_IMAGE_GENERIC(jobj1);
        ImageGeneric *image2 = GET_INTERNAL_IMAGE_GENERIC(jobj2);
        ImageGeneric *image3 = GET_INTERNAL_IMAGE_GENERIC(jobj3);
        *image3 = *image1;
        image3->subtract(*image2);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_smooth(JNIEnv *env, jobject jobj1, jobject jobj2)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image1 = GET_INTERNAL_IMAGE_GENERIC(jobj1);
        ImageGeneric *image2 = GET_INTERNAL_IMAGE_GENERIC(jobj2);

        downsampleSmooth(*image1, *image2);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT jdoubleArray JNICALL
Java_xmipp_jni_ImageGeneric_getStatistics(JNIEnv *env, jobject jobj)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

        double avg, stddev, min, max;
        (*image)().computeStats(avg, stddev, min, max);

        // Copies vector into array.
        double statistics[4];
        statistics[0] = min;
        statistics[1] = max;
        statistics[2] = avg;
        statistics[3] = stddev;

        // Sets array value
        jdoubleArray array = env->NewDoubleArray(4);
        env->SetDoubleArrayRegion(array, 0, 4, statistics);

        return array;
    }
    XMIPP_JAVA_CATCH;

    return NULL;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_setXmippOrigin(JNIEnv *env, jobject jobj)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        (*image)().setXmippOrigin();
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_convertPSD(JNIEnv *env, jobject jobj,
                                       jboolean useLogarithm)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        image->convert2Datatype(DT_Double);
        MultidimArray<double> *in;
        MULTIDIM_ARRAY_GENERIC(*image).getMultidimArrayPointer(in);
        xmipp2PSD(*in, *in, useLogarithm);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_ImageGeneric_generatePSDCTF
  (JNIEnv *env, jobject jobj, jobject md)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        image->convert2Datatype(DT_Double);
        MultidimArray<double> *in;
        MULTIDIM_ARRAY_GENERIC(*image).getMultidimArrayPointer(in);
        MetaData * mdC = GET_INTERNAL_METADATA(md);
        generatePSDCTFImage(*in, *mdC);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_ImageGeneric_generateImageWithTwoCTFs
  (JNIEnv *env, jobject jobj, jobject md1, jobject md2, jint xdim)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        image->setDatatype(DT_Double);
        image->resize(xdim,xdim,1,1);
        MultidimArray<double> *in;
        MULTIDIM_ARRAY_GENERIC(*image).getMultidimArrayPointer(in);
        MetaData * mdC1 = GET_INTERNAL_METADATA(md1);
        MetaData * mdC2 = GET_INTERNAL_METADATA(md2);
        generateCTFImageWith2CTFs(*mdC1, *mdC2, xdim, *in);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_getReslice(JNIEnv *env, jobject jobj, jobject jimgOut,
                                       jint view)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        ImageGeneric *imageOut = GET_INTERNAL_IMAGE_GENERIC(jimgOut);
        image->reslice((AxisView)view, *imageOut);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_reslice(JNIEnv *env, jobject jobj, jint view)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        image->reslice((AxisView)view);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_ImageGeneric_getPreview
(JNIEnv * env, jobject jobj, jobject jimgOut, jint xdim, jint ydim,
 jint select_slice, jlong select_image)
{
    XMIPP_JAVA_TRY
    {
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
        ImageGeneric *imageOut = GET_INTERNAL_IMAGE_GENERIC(jimgOut);
        image->getPreview(*imageOut, xdim, ydim, select_slice, select_image);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT jobject JNICALL Java_xmipp_jni_ImageGeneric_bestShift
(JNIEnv * env, jobject jobj, jobject jimg)
{
    XMIPP_JAVA_TRY
    {
    	ImageGeneric *templates = GET_INTERNAL_IMAGE_GENERIC(jobj);
    	ImageGeneric *img = GET_INTERNAL_IMAGE_GENERIC(jimg);
    	img->convert2Datatype(DT_Double);
    	templates->convert2Datatype(DT_Double);

        MultidimArray<double> *I, tmpI, tmpI2, alignedI, *Tp, T;

        MULTIDIM_ARRAY_GENERIC(*img).getMultidimArrayPointer(I);
        MULTIDIM_ARRAY_GENERIC(*templates).getMultidimArrayPointer(Tp);

        Matrix2D<double> alignmentM, transformM, bestAlignmentM;
        alignmentM.initIdentity(3);

        CorrelationAux aux2;
        ArrayDim dim;

        templates->getDimensions(dim);
        double corr, max = -MAXDOUBLE;

        for (size_t i=0; i < dim.ndim; ++i)
        {
        	T.aliasImageInStack(*Tp, i);
            tmpI=*I;
            tmpI2=*I;
            T.setXmippOrigin();
            tmpI.setXmippOrigin();
            tmpI2.setXmippOrigin();
            bestNonwrappingShift(T, tmpI, MAT_ELEM(alignmentM, 0, 2), MAT_ELEM(alignmentM, 1, 2), aux2);
            applyGeometry(LINEAR, tmpI, tmpI2, alignmentM, IS_NOT_INV, true);
            corr = correlationIndex(T, tmpI2);
            if (corr > max)
            {
                max=corr;
                bestAlignmentM = alignmentM;
            }
        }



        int x, y;
        x = -MAT_ELEM(bestAlignmentM, 0, 2);
        y = -MAT_ELEM(bestAlignmentM, 1, 2);

        jclass pclass = env->FindClass("xmipp/jni/Particle");
        jmethodID constructor = env->GetMethodID(pclass, "<init>", "(II)V");
        jobject particle = env->NewObject(pclass, constructor, x, y);

        return particle;
    }

        XMIPP_JAVA_CATCH;
}

JNIEXPORT jdoubleArray JNICALL Java_xmipp_jni_ImageGeneric_alignImage
(JNIEnv * env, jobject jobj, jobject jimg, jboolean jupdate)
{
	double result[4];
	XMIPP_JAVA_TRY
	{

		ImageGeneric *templates = GET_INTERNAL_IMAGE_GENERIC(jobj);
		ImageGeneric *img = GET_INTERNAL_IMAGE_GENERIC(jimg);
		img->convert2Datatype(DT_Double);
		templates->convert2Datatype(DT_Double);

		MultidimArray<double> *I, tmpI, alignedI, *Tp, T;
		AlignmentAux aux;
		CorrelationAux aux2;
		RotationalCorrelationAux aux3;
		Matrix2D<double> transformM;
		ArrayDim dim;
		MULTIDIM_ARRAY_GENERIC(*img).getMultidimArrayPointer(I);
		MULTIDIM_ARRAY_GENERIC(*templates).getMultidimArrayPointer(Tp);

		templates->getDimensions(dim);
		double corr,max=-MAXDOUBLE;
		int maxIndex=0;

		for (size_t i=0;i<dim.ndim;++i)
		{
			T.aliasImageInStack(*Tp,i);
			tmpI=*I;
			T.setXmippOrigin();
			tmpI.setXmippOrigin();
			corr = alignImages(T, tmpI, transformM, true, aux, aux2, aux3);
			if (corr>max)
			{
				max=corr;
				maxIndex=i;
				alignedI=tmpI;
			}
		}
		bool update = jupdate;
		if(update)
		{
			T.aliasImageInStack(*Tp, maxIndex);
			T+=alignedI;
			centerImage(T, aux2, aux3, 3);
		}
		double rot, tilt, psi;
		Euler_matrix2angles(transformM, rot, tilt, psi);

		result[0] = maxIndex;//template aligned with particle
		result[1] = rot;
		result[2] = tilt;
		result[3] = psi;



	}
	XMIPP_JAVA_CATCH;

	// Sets array value
	jdoubleArray jresult = env->NewDoubleArray(4);
	env->SetDoubleArrayRegion(jresult, 0, 4, result);
	return jresult;

}

JNIEXPORT void JNICALL Java_xmipp_jni_ImageGeneric_removeAlignment
(JNIEnv * env, jobject jobj, jobject jimg, jint index, jdouble rot, jdouble tilt, jdouble psi)
{
   XMIPP_JAVA_TRY
   {

	   ImageGeneric *templates = GET_INTERNAL_IMAGE_GENERIC(jobj);
	   ImageGeneric *img = GET_INTERNAL_IMAGE_GENERIC(jimg);
	   img->convert2Datatype(DT_Double);
	   templates->convert2Datatype(DT_Double);

	   MultidimArray<double> *Tp, T, *Ip, I;
	   MULTIDIM_ARRAY_GENERIC(*img).getMultidimArrayPointer(Ip);
	   MULTIDIM_ARRAY_GENERIC(*templates).getMultidimArrayPointer(Tp);

	   AlignmentAux aux;
	   CorrelationAux aux2;
	   RotationalCorrelationAux aux3;
	   Matrix2D<double> transformM;
	   Euler_angles2matrix(rot, tilt, psi, transformM);

	   T.aliasImageInStack(*Tp, index);
	   I=*Ip;
	   T.setXmippOrigin();
	   I.setXmippOrigin();
	   alignImages(T, I, transformM, true, aux, aux2, aux3);
	   T-=I;


   }
   XMIPP_JAVA_CATCH;


}






