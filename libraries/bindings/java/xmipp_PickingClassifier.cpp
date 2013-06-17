#include "xmipp_java_adapter.h"
#include <iostream>
#include "xmipp_PickingClassifier.h"
#include "xmipp_ExceptionsHandler.h"
#include "reconstruction/micrograph_automatic_picking2.h"
#include "xmipp_InternalData.h"
#include <data/metadata.h>

JNIEXPORT void JNICALL
Java_xmipp_jni_PickingClassifier_create(JNIEnv *env, jobject jobj, jint particle_size, jstring output)
{
    XMIPP_JAVA_TRY
    {
    	int size = particle_size, filter_num = 6, corr_num = 2, NPCA = 4;
    	const FileName &model_name = env->GetStringUTFChars(output, false);
    	AutoParticlePicking2 *picker = new AutoParticlePicking2(size, filter_num, corr_num, NPCA, model_name);
        STORE_PEER_ID(jobj, (long)picker);
    }
    XMIPP_JAVA_CATCH;
}


JNIEXPORT void JNICALL
Java_xmipp_jni_PickingClassifier_destroy(JNIEnv *env, jobject jobj)
{

    XMIPP_JAVA_TRY
    {
    	AutoParticlePicking2 *picker = GET_INTERNAL_AUTOPARTICLEPICKING2(jobj);
        delete picker;
        picker = NULL;
        STORE_PEER_ID(jobj, (long)picker);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_train
(JNIEnv *env, jobject jobj, jobject micrographs)
{
    XMIPP_JAVA_TRY
    {
    	AutoParticlePicking2 *picker = GET_INTERNAL_AUTOPARTICLEPICKING2(jobj);
    	MetaData * md = GET_INTERNAL_METADATA(micrographs);
    	picker->train(*md, false);

    }
    XMIPP_JAVA_CATCH;

}


JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_autopick
(JNIEnv *env, jobject jobj, jstring filename, jobject joutputmd, jint percent)
{
    XMIPP_JAVA_TRY
    {
    	AutoParticlePicking2 *picker = GET_INTERNAL_AUTOPARTICLEPICKING2(jobj);
    	MetaData * outputmd = GET_INTERNAL_METADATA(joutputmd);

    	const FileName micrograph = env->GetStringUTFChars(filename, false);
    	std::cout<<"autopick"<< std::endl;
    	picker->automaticallySelectParticles(micrograph, percent, *outputmd);
//    	outputmd->print();
//    	md.write(env->GetStringUTFChars(jautoposfile, false));
    }
    XMIPP_JAVA_CATCH;
}


JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_autopickOut
(JNIEnv *env, jobject jobj, jstring filename, jobject joutputmd, jint percent, jint x, jint y, jint width, jint height)
{
    XMIPP_JAVA_TRY
    {
    	AutoParticlePicking2 *picker = GET_INTERNAL_AUTOPARTICLEPICKING2(jobj);
    	MetaData * outputmd = GET_INTERNAL_METADATA(joutputmd);

    	const FileName micrograph = env->GetStringUTFChars(filename, false);
    	std::cout<<"autopick"<< std::endl;
    	picker->automaticallySelectParticles(micrograph, percent, *outputmd);
//    	outputmd->print();
//    	md.write(env->GetStringUTFChars(jautoposfile, false));
    }
    XMIPP_JAVA_CATCH;
}


JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_correct
(JNIEnv *env, jobject jobj, jobject jmanualmd, jobject jautomaticmd)
{
    XMIPP_JAVA_TRY
    {
    	MetaData * manualmd = GET_INTERNAL_METADATA(jmanualmd);
    	MetaData * automaticmd = GET_INTERNAL_METADATA(jautomaticmd);
    	AutoParticlePicking2 *picker = GET_INTERNAL_AUTOPARTICLEPICKING2(jobj);
    	std::cout<<"correct"<< std::endl;
    	manualmd->print();
    	automaticmd->print();
    	picker->correction(*manualmd, *automaticmd);
    }
    XMIPP_JAVA_CATCH;

}



JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_setSize
(JNIEnv *env, jobject jobj, jint psize)
{
    XMIPP_JAVA_TRY
    {
    	AutoParticlePicking2 *picker = GET_INTERNAL_AUTOPARTICLEPICKING2(jobj);
    	picker->setSize(psize);
    }
    XMIPP_JAVA_CATCH;
}

