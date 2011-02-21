#include <jni.h>
#include "xmipp_ImageDouble.h"
#include "xmipp_ExceptionsHandler.h"
#include <data/image.h>
#include <data/fft.h>

static jfieldID ImageDouble_peerId;
static jfieldID ImageDouble_filenameId;

#define peerId ImageDouble_peerId
#define GET_INTERNAL_IMAGE() ((Image<double>*)(env->GetLongField(obj, peerId)))

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_storeIds
(JNIEnv *env, jclass cls) {
	peerId = env->GetFieldID(cls, "peer", "J");
}

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_create
(JNIEnv *env, jobject obj) {
	Image<double> * image = new Image<double>();
	env->SetLongField(obj, peerId, (long)image);
}

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_destroy
(JNIEnv *env, jobject obj) {
	Image<double> * image = GET_INTERNAL_IMAGE();
	delete image;
	image = NULL;
	env->SetLongField(obj, peerId, (long)image);
}

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_read
(JNIEnv *env, jobject obj, jstring filename, jboolean read_data, jlong nimage) {
	std::string msg = "";
	Image<double> * image = GET_INTERNAL_IMAGE();

	if (image != NULL) {
		const char * fnStr = env->GetStringUTFChars(filename, false);

		try {
			image->read(fnStr, (read_data)? DATA : HEADER, (size_t)nimage);
		} catch (XmippError xe) {
			msg = xe.getDefaultMessage();
		} catch (std::exception& e) {
			msg = e.what();
		} catch (...) {
			msg = "Unhandled exception";
		}
	} else {
		msg = "Image is null";
	}

	// If there was an exception, sends it to java environment.
	if(!msg.empty()) {
		handleXmippException(env, msg);
	}
}

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_readPreview
(JNIEnv *env, jobject obj, jstring filename, jint w, jint h, jint slice, jlong nimage) {
	std::string msg = "";
	Image<double> *image = GET_INTERNAL_IMAGE();

	if (image != NULL) {
		const char *fnStr = env->GetStringUTFChars(filename, false);

		try {
			image->readPreview(fnStr, (int)w, (int)h, (int)slice, (size_t)nimage);
		} catch (XmippError xe) {
			msg = xe.getDefaultMessage();
		} catch (std::exception& e) {
			msg = e.what();
		} catch (...) {
			msg = "Unhandled exception";
		}
	} else {
		msg = "Image is null";
	}

	// If there was an exception, sends it to java environment.
	if(!msg.empty()) {
		handleXmippException(env, msg);
	}
}

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_setData
(JNIEnv *env, jobject obj, jint w, jint h, jint d, jdoubleArray data) {
	std::string msg = "";
	Image<double> *image = GET_INTERNAL_IMAGE();

	if (image != NULL) {
		try {
			image->data.resizeNoCopy((int) d, (int) h, (int) w);
			env->GetDoubleArrayRegion(data, 0, w * h * d, MULTIDIM_ARRAY(image->data));
		} catch (XmippError xe) {
			msg = xe.getDefaultMessage();
		} catch (std::exception& e) {
			msg = e.what();
		} catch (...) {
			msg = "Unhandled exception";
		}
	} else {
		msg = "Image is null";
	}

	// If there was an exception, sends it to java environment.
	if(!msg.empty()) {
		handleXmippException(env, msg);
	}
}

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_write
(JNIEnv *env, jobject obj, jstring filename) {
	std::string msg = "";
	Image<double> * image = GET_INTERNAL_IMAGE();

	if (image != NULL) {
		const char * fnStr = env->GetStringUTFChars(filename, false);

		try {
			image->write(fnStr);
		} catch (XmippError xe) {
			msg = xe.getDefaultMessage();
		} catch (std::exception& e) {
			msg = e.what();
		} catch (...) {
			msg = "Unhandled exception";
		}
	} else {
		msg = "Image is null";
	}

	// If there was an exception, sends it to java environment.
	if(!msg.empty()) {
		handleXmippException(env, msg);
	}
}

JNIEXPORT jdoubleArray JNICALL Java_xmipp_ImageDouble_getData(JNIEnv *env,
		jobject obj) {
	Image<double> * image = GET_INTERNAL_IMAGE();
	if (image != NULL) {
		size_t size = image->getSize();
		jdoubleArray array = env->NewDoubleArray(size);
		env->SetDoubleArrayRegion(array, 0, size, MULTIDIM_ARRAY(image->data));
		return array;
	}

	return (jdoubleArray) NULL;
}

JNIEXPORT jdouble JNICALL Java_xmipp_ImageDouble_getPixel(JNIEnv *env,
		jobject obj, jint x, jint y) {
	double value = 0;
	Image<double> * image = GET_INTERNAL_IMAGE();

	if (image != NULL) {
		value = image->getPixel((int) x, (int) y);
	}

	return value;
}

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_setPixel
(JNIEnv *env, jobject obj, jint x, jint y, jdouble value) {
	Image<double> * image = GET_INTERNAL_IMAGE();

	if (image != NULL) {
		image->setPixel((int)x, (int)y, (double)value);
	}
}

JNIEXPORT jdouble JNICALL Java_xmipp_ImageDouble_getVoxel(JNIEnv *env,
		jobject obj, jint x, jint y, jint z) {
	double value = 0;
	Image<double> * image = GET_INTERNAL_IMAGE();

	if (image != NULL) {
		value = image->data.getVoxel((int) x, (int) y, (int) z);
	}

	return value;
}

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_setVoxel
(JNIEnv *env, jobject obj, jint x, jint y, jint z, jdouble value) {
	Image<double> * image = GET_INTERNAL_IMAGE();

	if (image != NULL) {
		image->data.setVoxel((int)x, (int)y, (int)z, (double)value);
	}
}

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_convertPSD(JNIEnv *env, jobject obj) {
	Image<double> * image = GET_INTERNAL_IMAGE();
	if (image != NULL) {
		MultidimArray<double> out(image->getSize());
		xmipp2PSD(image->data, out);
		image->data.clear(); // Frees memory.
		image->data = out;
	}
}

JNIEXPORT jint JNICALL Java_xmipp_ImageDouble_getXsize(JNIEnv * env,
		jobject obj) {
	Image<double> * image = GET_INTERNAL_IMAGE();

	if (image != NULL) {
		return XSIZE(image->data);
	}

	return 0;
}

JNIEXPORT jint JNICALL Java_xmipp_ImageDouble_getYsize(JNIEnv * env,
		jobject obj) {
	Image<double> * image = GET_INTERNAL_IMAGE();

	if (image != NULL) {
		return YSIZE(image->data);
	}

	return 0;
}

JNIEXPORT jint JNICALL Java_xmipp_ImageDouble_getZsize(JNIEnv * env,
		jobject obj) {
	Image<double> * image = GET_INTERNAL_IMAGE();

	if (image != NULL) {
		return ZSIZE(image->data);
	}

	return 0;
}

JNIEXPORT jlong JNICALL Java_xmipp_ImageDouble_getNsize(JNIEnv * env,
		jobject obj) {
	Image<double> * image = GET_INTERNAL_IMAGE();

	if (image != NULL) {
		return NSIZE(image->data);
	}

	return 0;
}

JNIEXPORT void JNICALL Java_xmipp_ImageDouble_printShape
(JNIEnv *env, jobject obj) {
	Image<double> * image = GET_INTERNAL_IMAGE();
	std::cout << (*image) << std::endl;
}
