#include <jni.h>
#include "xmipp_ImageGeneric.h"
#include "xmipp_InternalData.h"
#include "xmipp_ExceptionsHandler.h"
#include <data/xmipp_image_generic.h>
#include <data/xmipp_fft.h>

JNIEXPORT void JNICALL Java_xmipp_ImageGeneric_storeIds
(JNIEnv *env, jclass cls)
{
    ImageGeneric_peerId = env->GetFieldID(cls, "peer", "J");
}

JNIEXPORT void JNICALL Java_xmipp_ImageGeneric_create
(JNIEnv *env, jobject jobj)
{
    ImageGeneric *image = new ImageGeneric();
    env->SetLongField(jobj, ImageGeneric_peerId, (long)image);
}

JNIEXPORT void JNICALL Java_xmipp_ImageGeneric_destroy
(JNIEnv *env, jobject jobj)
{
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
    delete image;
    image = NULL;
    env->SetLongField(jobj, ImageGeneric_peerId, (long)image);
}

/*
 * Class:     xmipp_ImageGeneric
 * Method:    readHeader
 * Signature: (Ljava/lang/String;)V
 */
JNIEXPORT void JNICALL Java_xmipp_ImageGeneric_read
(JNIEnv * env, jobject obj, jstring filename, jboolean onlyHeader)
{
    String msg = "";
    try
    {
        const char *fnStr = env->GetStringUTFChars(filename, false);
        ImageInfo imf;
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(obj);

        image->read(fnStr, onlyHeader ? HEADER : DATA);

        image->getInfo(imf);
		//Set object properties
		jclass class_ = env->GetObjectClass(obj);

        jfieldID fieldId = env->GetFieldID(class_, "xSize", "I");
        env->SetIntField(obj, fieldId, imf.adim.xdim);
        fieldId = env->GetFieldID(class_, "ySize", "I");
        env->SetIntField(obj, fieldId, imf.adim.ydim);
        fieldId = env->GetFieldID(class_, "zSize", "I");
        env->SetIntField(obj, fieldId, imf.adim.zdim);
        fieldId = env->GetFieldID(class_, "nSize", "J");
        env->SetIntField(obj, fieldId, imf.adim.ndim);
        fieldId = env->GetFieldID(class_, "dataType", "I");
        env->SetIntField(obj, fieldId, (int)imf.datatype);
    }
    catch (XmippError &xe)
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

JNIEXPORT void JNICALL Java_xmipp_ImageGeneric_write
(JNIEnv * env, jobject obj, jstring filename)
{
    String msg = "";
    try
    {
        const char *fnStr = env->GetStringUTFChars(filename, false);
        ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(obj);
        image->write(fnStr);
    }
    catch (XmippError &xe)
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

/*
 * Class:     xmipp_ImageGeneric
 * Method:    getArrayByte
 * Signature: (Ljava/lang/String;IIIJ)[F
 */
JNIEXPORT jbyteArray JNICALL Java_xmipp_ImageGeneric_getArrayByte(JNIEnv * env,
		jclass class_, jstring filename, jint jx, jint jy, jint jz, jlong jn,
		jint datatype) {
	String msg = "";
	try {
		const char *fnStr = env->GetStringUTFChars(filename, false);
		jbyteArray array;
		size_t size;
		//Image<char> image;
		switch (datatype) {
		case UChar: {
			Image<unsigned char> image;
			image.readOrReadPreview(fnStr, jx, jy, jz, jn, true);
			size = image.getSize();
			array = env->NewByteArray(size);
			env->SetByteArrayRegion(array, 0, size,
					(jbyte*) MULTIDIM_ARRAY(image.data));
		}
			break;
		case SChar: {
			Image<char> image;
			image.readOrReadPreview(fnStr, jx, jy, jz, jn, true);
			size = image.getSize();
			array = env->NewByteArray(size);
			unsigned char * buffer = new unsigned char[rw_max_page_size];
			char * imageArray = MULTIDIM_ARRAY(image.data);
			size_t page_size = rw_max_page_size;

			for (size_t written = 0; written < size; written += page_size) {
				page_size = std::min(page_size, size - written);
				char * iter = imageArray + written;

				for (size_t k = 0; k < page_size; ++k)
					buffer[k] = iter[k] - CHAR_MIN;

				env->SetByteArrayRegion(array, written, page_size,
						(jbyte*) buffer);
			}

		}
			break;
		default:
			REPORT_ERROR(
					ERR_ARG_INCORRECT,
					(String)"Java_xmipp_ImageGeneric_getArrayByte: reading invalid image datatype: " + datatype2Str((DataType)datatype));
		}
		return array;
	} catch (XmippError &xe) {
		msg = xe.getDefaultMessage();
	} catch (std::exception &e) {
		msg = e.what();
	} catch (...) {
		msg = "Unhandled exception";
	}

	// If there was an exception, sends it to java environment.
	if (!msg.empty()) {
		handleXmippException(env, msg);
	}
}

/*
 * Class:     xmipp_ImageGeneric
 * Method:    getArrayShort
 * Signature: (Ljava/lang/String;IIIJ)[S
 */
JNIEXPORT jshortArray JNICALL Java_xmipp_ImageGeneric_getArrayShort(
		JNIEnv * env, jclass class_, jstring filename, jint jx, jint jy,
		jint jz, jlong jn, jint datatype) {
	String msg = "";

	try {
		const char *fnStr = env->GetStringUTFChars(filename, false);
		Image<short> image;

		image.readOrReadPreview(fnStr, jx, jy, jz, jn, true);

		size_t size = image.getSize();

		jshortArray array = env->NewShortArray(size);
		env->SetShortArrayRegion(array, 0, size, MULTIDIM_ARRAY(image.data));

		return array;
	} catch (XmippError &xe) {
		msg = xe.getDefaultMessage();
	} catch (std::exception& e) {
		msg = e.what();
	} catch (...) {
		msg = "Unhandled exception";
	}

	// If there was an exception, sends it to java environment.
	if (!msg.empty()) {
		handleXmippException(env, msg);
	}
}

/*
 * Class:     xmipp_ImageGeneric
 * Method:    getArrayFloat
 * Signature: (Ljava/lang/String;IIIJ)[F
 */
JNIEXPORT jfloatArray JNICALL Java_xmipp_ImageGeneric_getArrayFloat(
		JNIEnv * env, jclass class_, jstring filename, jint jx, jint jy,
		jint jz, jlong jn, jint datatype) {
	String msg = "";
	try {
		const char *fnStr = env->GetStringUTFChars(filename, false);
		Image<float> image;

		image.readOrReadPreview(fnStr, jx, jy, jz, jn, true);
		//image.data.printShape(std::cerr);
		size_t size = image.getSize();
		jfloatArray array = env->NewFloatArray(size);
		env->SetFloatArrayRegion(array, 0, size, MULTIDIM_ARRAY(image.data));
		return array;
	} catch (XmippError &xe) {
		msg = xe.getDefaultMessage();
	} catch (std::exception& e) {
		msg = e.what();
	} catch (...) {
		msg = "Unhandled exception";
	}

	// If there was an exception, sends it to java environment.
	if (!msg.empty()) {
		handleXmippException(env, msg);
	}
}


JNIEXPORT jfloatArray JNICALL Java_xmipp_ImageGeneric_setArrayFloat
  (JNIEnv *env, jobject jobj, jint x, jint y, jint z, jlong N, jint datatype, jfloatArray data) {
    ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
    image->setDatatype((DataType)Float);//datatype);
    MultidimArrayGeneric & mag = MULTIDIM_ARRAY_GENERIC(*image);
    mag.resize(N, z, y, x, false);
    float * fdata;
    mag.getMultidimArrayPointer(fdata);
    std::cerr << "DEBUG_JM: N*z*y*x: " << N*z*y*x << std::endl;
    size_t n = env->GetArrayLength(data);
    std::cerr << "DEBUG_JM: n: " << n << std::endl;
    env->GetFloatArrayRegion(data, 0, x*y*z*N, fdata);
    std::cerr << "DEBUG_JM: data after GetFloatArrayRegion" << std::endl;
      for (int i = 0; i < N*z*y*x; ++i)
        std::cerr << " " << fdata[i];
      std::cerr << std::endl;
    image->print();
    image->write("tt.xmp");
}

//JNIEXPORT void JNICALL Java_xmipp_ImageGeneric_read
//(JNIEnv *env, jobject jobj, jstring filename) {
// std::string msg = "";
// ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);
//
// if (image != NULL) {
//  try {
//   const char * fnStr = env->GetStringUTFChars(filename, false);
//
//   image->readMapped(fnStr);
//  } catch (XmippError xe) {
//   msg = xe.getDefaultMessage();
//  } catch (std::exception& e) {
//   msg = e.what();
//  } catch (...) {
//   msg = "Unhandled exception";
//  }
// } else {
//  msg = "Metadata is null";
// }
//
// // If there was an exception, sends it to java environment.
// if(!msg.empty()) {
//  handleXmippException(env, msg);
// }
//}

JNIEXPORT jdoubleArray JNICALL Java_xmipp_ImageGeneric_getStatistics(
		JNIEnv *env, jobject jobj) {
	std::string msg = "";
	ImageGeneric *image = GET_INTERNAL_IMAGE_GENERIC(jobj);

	if (image != NULL) {
		try {
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
		} catch (XmippError &xe) {
			msg = xe.getDefaultMessage();
		} catch (std::exception& e) {
			msg = e.what();
		} catch (...) {
			msg = "Unhandled exception";
		}
	} else {
		msg = "Metadata is null";
	}

	// If there was an exception, sends it to java environment.
	if (!msg.empty()) {
		handleXmippException(env, msg);
	}

	return NULL;
}

