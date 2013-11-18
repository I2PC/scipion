#include "xmipp_java_adapter.h"

#include <jni.h>
#include "xmipp_MDRow.h"
#include "xmipp_ExceptionsHandler.h"
#include "xmipp_InternalData.h"

/*
JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_storeIds
(JNIEnv *env, jclass cls) {
 MDRow_peerId = env->GetFieldID(cls, "peer", "J");
}
*/
JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_create
(JNIEnv *env, jobject jobj)
{
    MDRow * mdRow = new MDRow();
    STORE_PEER_ID(jobj, mdRow);
}

JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_destroy
(JNIEnv *env, jobject jobj)
{
    MDRow * mdRow = GET_INTERNAL_MDROW(jobj);
    delete mdRow;
    mdRow = NULL;
    STORE_PEER_ID(jobj, mdRow);
}

JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_clear
(JNIEnv *env, jobject jobj)
{
  XMIPP_JAVA_TRY
  {
      MDRow * mdRow = GET_INTERNAL_MDROW(jobj);
      mdRow->clear();
  }
  XMIPP_JAVA_CATCH;
}

JNIEXPORT jint JNICALL Java_xmipp_jni_MDRow_size(JNIEnv *env, jobject jobj)
{
    XMIPP_JAVA_TRY
    {
      MDRow * mdRow = GET_INTERNAL_MDROW(jobj);
        return mdRow->size();
    }
    XMIPP_JAVA_CATCH;

    return 0;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MDRow_containsLabel
(JNIEnv *env, jobject jobj, jint label)
{
    XMIPP_JAVA_TRY
    {
        MDRow * mdRow = GET_INTERNAL_MDROW(jobj);


        return mdRow->containsLabel((MDLabel) label);
    }
    XMIPP_JAVA_CATCH;

    return false;
}


JNIEXPORT jint JNICALL Java_xmipp_jni_MDRow_getValueInt(JNIEnv *env,
        jobject jobj, jint label)
{
    XMIPP_JAVA_TRY
    {
        MDRow * mdRow = GET_INTERNAL_MDROW(jobj);


        int value;

        if (mdRow->getValue((MDLabel) label, value))
            return value;
    }
    XMIPP_JAVA_CATCH;

    return 0;
}

JNIEXPORT jlong JNICALL Java_xmipp_jni_MDRow_getValueLong
(JNIEnv *env, jobject jobj, jint label)
{
    XMIPP_JAVA_TRY
    {
        MDRow * mdRow = GET_INTERNAL_MDROW(jobj);


        size_t value;
        if (mdRow->getValue((MDLabel) label, value))
            return value;
    }
    XMIPP_JAVA_CATCH;

    return 0;
}

JNIEXPORT jdouble JNICALL Java_xmipp_jni_MDRow_getValueDouble(JNIEnv *env,
        jobject jobj, jint label)
{
    XMIPP_JAVA_TRY
    {
        MDRow * mdRow = GET_INTERNAL_MDROW(jobj);

        double value;
        if (mdRow->getValue((MDLabel) label, value))
            return value;
    }
    XMIPP_JAVA_CATCH;

    return 0;
}

JNIEXPORT jstring JNICALL Java_xmipp_jni_MDRow_getValueString(JNIEnv *env,
        jobject jobj, jint label)
{
    XMIPP_JAVA_TRY
    {
        MDRow * mdRow = GET_INTERNAL_MDROW(jobj);

        MDObject obj((MDLabel) label);
        if (mdRow->getValue(obj))
        {
            jstring jstr = env->NewStringUTF(obj.toString().c_str());
            return jstr;
        }
    }
    XMIPP_JAVA_CATCH;

    return NULL;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MDRow_getValueBoolean(JNIEnv *env,
        jobject jobj, jint label)
{
    XMIPP_JAVA_TRY
    {
        MDRow * mdRow = GET_INTERNAL_MDROW(jobj);


        bool value;
        if (mdRow->getValue((MDLabel) label, value))
            return value;
    }
    XMIPP_JAVA_CATCH;

    return 0;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_setValueInt(JNIEnv *env,
        jobject jobj, jint label, jint value)
{
    XMIPP_JAVA_TRY
    {
        MDRow * mdRow = GET_INTERNAL_MDROW(jobj);
        mdRow->setValue((MDLabel) label, (int)value);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_setValueLong(JNIEnv *env,
        jobject jobj, jint label, jlong value)
{
    XMIPP_JAVA_TRY
    {
        MDRow * mdRow = GET_INTERNAL_MDROW(jobj);
        mdRow->setValue((MDLabel) label, (size_t)value);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_setValueDouble(JNIEnv *env,
        jobject jobj, jint label, jdouble value)
{
    XMIPP_JAVA_TRY
    {
        MDRow * mdRow = GET_INTERNAL_MDROW(jobj);
        mdRow->setValue((MDLabel) label, value);
    }
    XMIPP_JAVA_CATCH;
}
JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_setValueString(JNIEnv *env,
        jobject jobj, jint label, jstring value)
{
    XMIPP_JAVA_TRY
    {
        MDRow * mdRow = GET_INTERNAL_MDROW(jobj);

        jboolean aux = false;
        const char * strValue = env->GetStringUTFChars(value, &aux);
        mdRow->setValueFromStr((MDLabel) label, strValue);
        env->ReleaseStringUTFChars(value, strValue);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_setValueBoolean(JNIEnv *env,
        jobject jobj, jint label, jboolean value)
{
    XMIPP_JAVA_TRY
    {
        MDRow * mdRow = GET_INTERNAL_MDROW(jobj);
        mdRow->setValue((MDLabel) label, (bool) value);
    }
    XMIPP_JAVA_CATCH;
}


//JNIEXPORT void JNICALL Java_xmipp_jni_MDRow_addLabel(JNIEnv *env, jobject jobj, jint label)
//{
//    MDRow *mdRow = GET_INTERNAL_MDROW(jobj);
//
//    XMIPP_JAVA_TRY
//    {
//        mdRow->addLabel((MDLabel)label);
//    }
//    XMIPP_JAVA_CATCH;
//}

