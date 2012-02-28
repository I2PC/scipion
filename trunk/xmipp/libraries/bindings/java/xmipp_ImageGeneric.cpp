#include "xmipp_java_adapter.h"

#include <jni.h>
#include "xmipp_ImageGeneric.h"
#include "xmipp_InternalData.h"
#include "xmipp_ExceptionsHandler.h"
#include <data/xmipp_image_generic.h>
#include <data/xmipp_fft.h>

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_create(JNIEnv *env, jobject jobj)
{
  ImageGeneric *image = new ImageGeneric();
  STORE_PEER_ID(jobj, (long)image);
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_destroy(JNIEnv *env, jobject jobj)
{
  ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
  delete image;
  image = NULL;
  STORE_PEER_ID(jobj, (long)image);
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_resize(JNIEnv *env, jobject jobj, jint w, jint h,
    jint d, jlong n)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

    image->resize(w, h, d, (size_t) n, false);
  }
  XMIPP_CATCH;
}

JNIEXPORT jint JNICALL
Java_xmipp_jni_ImageGeneric_getXDim(JNIEnv *env, jobject jobj)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

    return XSIZE(* image->data->im);
  }
  XMIPP_CATCH;

  return -1;
}

JNIEXPORT jint JNICALL
Java_xmipp_jni_ImageGeneric_getYDim(JNIEnv *env, jobject jobj)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

    return YSIZE(* image->data->im);
  }
  XMIPP_CATCH;

  return -1;
}

JNIEXPORT jint JNICALL
Java_xmipp_jni_ImageGeneric_getZDim(JNIEnv *env, jobject jobj)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

    return ZSIZE(* image->data->im);
  }
  XMIPP_CATCH;

  return -1;
}

JNIEXPORT jlong JNICALL
Java_xmipp_jni_ImageGeneric_getNDim(JNIEnv *env, jobject jobj)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

    return NSIZE(* image->data->im);
  }
  XMIPP_CATCH;

  return -1;
}

JNIEXPORT jint JNICALL
Java_xmipp_jni_ImageGeneric_getDataType(JNIEnv *env, jobject jobj)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

    return image->getDatatype();
  }
  XMIPP_CATCH;

  return -1;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_readHeader(JNIEnv * env, jobject jobj,
    jstring filename)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
    const char *fnStr = env->GetStringUTFChars(filename, false);
    image->read(fnStr, HEADER);
  }
  XMIPP_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_read(JNIEnv *env, jobject jobj, jstring filename,
    jint jx, jint jy, jint jz, jlong jn)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
    const char *fn = env->GetStringUTFChars(filename, false);
    image->readOrReadPreview(fn, jx, jy, jz, jn, true);
  }
  XMIPP_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_readApplyGeo_1(JNIEnv *env, jobject jimage,
    jstring filename, jobject jmetadata, jlong id, jint w, jint h)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jimage);
    MetaData *metadata = GET_INTERNAL_METADATA(jmetadata);

    const char *fnStr = env->GetStringUTFChars(filename, false);
    image->readApplyGeo(fnStr, *metadata, (size_t) id);

#define SELF_SCALE_TO_SIZE(type) MultidimArray<type> *mdaFloat; \
 MULTIDIM_ARRAY_GENERIC(*image).getMultidimArrayPointer(mdaFloat); \
 selfScaleToSize(LINEAR, *mdaFloat, (int)w, (int)h);

    SWITCHDATATYPE(image->getDatatype(), SELF_SCALE_TO_SIZE);
  }
  XMIPP_CATCH;
}

JNIEXPORT jbyteArray JNICALL
Java_xmipp_jni_ImageGeneric_getArrayByte(JNIEnv *env, jobject jobj, jint nslice)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

    // Go to slice.
    image->movePointerToSlice(nslice);

    size_t size = image->getSize();
    jbyteArray array = env->NewByteArray(size);

    DataType dataType = image->getDatatype();

    switch (dataType)
    {
      case UChar:
      {
        unsigned char *mdarray;

        // Get slice array.
        MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);

        env->SetByteArrayRegion(array, 0, size, (jbyte *) mdarray);
      }
        break;
      case SChar:
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
      }
        break;
      default:
        REPORT_ERROR(
            ERR_ARG_INCORRECT,
            (String)"Java_xmipp_jni_ImageGeneric_getArrayByte: reading invalid image datatype: " + datatype2Str((DataType)dataType));
    }

    // Resets slice pointer.
    image->movePointerToSlice();

    return array;
  }
  XMIPP_CATCH;

  return NULL;
}

JNIEXPORT jshortArray JNICALL
Java_xmipp_jni_ImageGeneric_getArrayShort(JNIEnv *env, jobject jobj,
    jint nslice)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

    // Go to slice.
    image->movePointerToSlice(nslice);

    size_t size = image->getSize();
    jshortArray array = env->NewShortArray(size);

    DataType dataType = image->getDatatype();

    switch (dataType)
    {
      case UShort:
      {
        unsigned short *mdarray;

        // Get slice array.
        MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);

        env->SetShortArrayRegion(array, 0, size, (jshort *) mdarray);
      }
        break;
      case Short:
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
      }
        break;
      default:
        REPORT_ERROR(
            ERR_ARG_INCORRECT,
            (String)"Java_xmipp_jni_ImageGeneric_getArrayShort: reading invalid image datatype: " + datatype2Str((DataType)dataType));
    }

    // Resets slice pointer.
    image->movePointerToSlice();

    return array;
  }
  XMIPP_CATCH;

  return NULL;
}

JNIEXPORT jfloatArray JNICALL
Java_xmipp_jni_ImageGeneric_getArrayFloat(JNIEnv *env, jobject jobj,
    jint nslice)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

    // Go to slice.
    image->movePointerToSlice(nslice);

    size_t size = image->getSize();
    jfloatArray array = env->NewFloatArray(size);

    DataType dataType = image->getDatatype();

    switch (dataType)
    {
      case Float:
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
                char * buffer = new char[rw_max_page_size * gettypesize(Float)]; \
                for (size_t written = 0; written < size; written += page_size) \
                { \
                    page_size = std::min(page_size, size - written); \
                    imageAux->castPage2Datatype(data + written, buffer, Float, page_size); \
                    env->SetFloatArrayRegion(array, written, page_size, (jfloat*) buffer); \
                } \

        SWITCHDATATYPE(dataType, CAST_PAGE);

      }
        break;
    }

    // Resets slice pointer.
    image->movePointerToSlice();

    return array;
  }
  XMIPP_CATCH;

  return NULL;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_setArrayByte(JNIEnv *env, jobject jobj,
    jbyteArray data, jint nslice)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

    // Go to slice.
    image->movePointerToSlice(nslice);

    size_t size = image->getSize();
    jbyteArray array = env->NewByteArray(size);

    DataType dataType = image->getDatatype();

    switch (dataType)
    {
      case UChar:
      {
        unsigned char *mdarray;

        // Get slice array.
        MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);

        env->GetByteArrayRegion(data, 0, size, (jbyte *) mdarray);
      }
        break;
      case SChar:
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
      }
        break;
      default:
      {
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
                    type * iter = mdarray + written; \
                    imageAux->castPage2T(buffer, iter, UChar, page_size); \
                } \
                SWITCHDATATYPE(image->getDatatype(), CAST_PAGE);

      }
    }

    // Resets slice pointer.
    image->movePointerToSlice();
  }
  XMIPP_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_setArrayShort(JNIEnv *env, jobject jobj,
    jshortArray data, jint nslice)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

    // Go to slice.
    image->movePointerToSlice(nslice);

    size_t size = image->getSize();

    switch (image->getDatatype())
    {
      case UShort:
      {
        unsigned short *mdarray;

        // Get slice array.
        MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);

        env->GetShortArrayRegion(data, 0, size, (jshort *) mdarray);
      }
        break;
      case Short:
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
      }
        break;
      case Float:
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

          float * iter = mdarray + written;
          imageAux->castPage2T(buffer, iter, UShort, page_size);
        }
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
    image->movePointerToSlice();
  }
  XMIPP_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_setArrayFloat(JNIEnv *env, jobject jobj,
    jfloatArray data, jint nslice)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

    // Go to slice.
    image->movePointerToSlice(nslice);

    size_t size = image->getSize();
    //jfloatArray array = env->NewFloatArray(size);

    switch (image->getDatatype())
    {
      case Float:
      {
        float *mdarray;

        // Get slice array.
        MULTIDIM_ARRAY_GENERIC(*image).getArrayPointer(mdarray);

        env->GetFloatArrayRegion(data, 0, size, mdarray);
      }
        break;
      default:
      {
#define CAST_PAGE(type) Image<type> *imageAux = (Image<type> *)image->image; \
             type *data = MULTIDIM_ARRAY(imageAux->data); \
             size_t page_size = rw_max_page_size; \
             char * buffer = new char[rw_max_page_size * gettypesize(Float)]; \
             for (size_t written = 0; written < size; written += page_size) \
             { \
              page_size = std::min(page_size, size - written); \
              imageAux->castPage2Datatype(data + written, buffer, Float, page_size); \
              env->SetFloatArrayRegion(array, written, page_size, (jfloat*) buffer); \
             } \
             SWITCHDATATYPE(dataType, CAST_PAGE);

      }
        break;
    }

    // Resets slice pointer.
    image->movePointerToSlice();
  }
  XMIPP_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_setDataType(JNIEnv *env, jobject jobj,
    jint dataType)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
    image->setDatatype((DataType) dataType);
  }
  XMIPP_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_convert2Datatype(JNIEnv *env, jobject jobj,
    jint dataType)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
    image->convert2Datatype((DataType) dataType);
  }
  XMIPP_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_mapFile2Write(JNIEnv *env, jobject jobj, jint w,
    jint h, jint z, jstring filename, jlong n)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
    const char *fnStr = env->GetStringUTFChars(filename, false);

    image->mapFile2Write(w, h, z, fnStr, false, (size_t) n);
  }
  XMIPP_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_write(JNIEnv *env, jobject jobj, jstring filename)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
    const char *fnStr = env->GetStringUTFChars(filename, false);
    image->write(fnStr);
  }
  XMIPP_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_printShape(JNIEnv *env, jobject jobj)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
    image->print();
  }
  XMIPP_CATCH;
}

JNIEXPORT jboolean JNICALL
Java_xmipp_jni_ImageGeneric_equal(JNIEnv *env, jobject jobj1, jobject jobj2, jdouble accuracy)
{
  XMIPP_TRY
  {
	   ImageGeneric *image1 = GET_INTERNAL_IMAGE_GENERIC(jobj1);
	   ImageGeneric *image2 = GET_INTERNAL_IMAGE_GENERIC(jobj2);
	   return image1->equal(*image2,accuracy);
  }
  XMIPP_CATCH;
}

JNIEXPORT jdoubleArray JNICALL
Java_xmipp_jni_ImageGeneric_getStatistics(JNIEnv *env, jobject jobj)
{
  XMIPP_TRY
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
  XMIPP_CATCH;

  return NULL;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_setXmippOrigin(JNIEnv *env, jobject jobj)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
    (*image)().setXmippOrigin();
  }
  XMIPP_CATCH;
}

JNIEXPORT void JNICALL
Java_xmipp_jni_ImageGeneric_convertPSD(JNIEnv *env, jobject jobj,
    jboolean useLogarithm)
{
  XMIPP_TRY
  {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

    MultidimArray<double> out(image->getSize());

    image->convert2Datatype(Double);
    MultidimArray<double> *in;
    MULTIDIM_ARRAY_GENERIC(*image).getMultidimArrayPointer(in);

    xmipp2PSD(*in, out, useLogarithm);

    image->clear(); // Frees memory.
    image->setDatatype(Float);
    MULTIDIM_ARRAY_GENERIC(*image).setImage(out);
  }
  XMIPP_CATCH;
}
