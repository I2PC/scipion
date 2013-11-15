#include "xmipp_java_adapter.h"
#include <iostream>
#include "xmipp_PickingClassifier.h"
#include "xmipp_ExceptionsHandler.h"
#include "reconstruction/micrograph_automatic_picking2.h"
#include "xmipp_InternalData.h"
#include <data/metadata.h>

JNIEXPORT void JNICALL
Java_xmipp_jni_PickingClassifier_create(JNIEnv *env, jobject jobj, jint particle_size, jstring output, jstring micsFn)
{
    XMIPP_JAVA_TRY
    {
        int size = particle_size, filter_num = 6, corr_num = 2, NPCA = 4;
        jboolean aux=false;
        const FileName &model_name = env->GetStringUTFChars(output, &aux);
        const FileName &mics = env->GetStringUTFChars(micsFn, &aux);
        MetaData micList;
        std::vector<MDRow> vMicList;

        micList.read(mics);
        micList.metadataToVec(vMicList);

        AutoParticlePicking2 *picker = new AutoParticlePicking2(size, filter_num, corr_num, NPCA, model_name, vMicList);

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
(JNIEnv *env, jobject jobj, jobject micrographs, jint x, jint y, jint width, jint height)
{
    XMIPP_JAVA_TRY
    {
        AutoParticlePicking2 *picker = GET_INTERNAL_AUTOPARTICLEPICKING2(jobj);
        MetaData * md = GET_INTERNAL_METADATA(micrographs);
        MetaData md2;
        FileName micFile,posFile;
        std::vector<MDRow> vd;
        MDRow row;
        int cnt=0,x,y;
        FOR_ALL_OBJECTS_IN_METADATA((*md))
        {
            md->getValue(MDL_MICROGRAPH,micFile, __iter.objId);
            md->getValue(MDL_MICROGRAPH_PARTICLES,posFile, __iter.objId);
            md2.read("particles@"+posFile);

            FOR_ALL_OBJECTS_IN_METADATA(md2)
            {
                row.setValue(MDL_MICROGRAPH,micFile);
                md2.getValue(MDL_XCOOR,x, __iter.objId);
                row.setValue(MDL_XCOOR,x);
                md2.getValue(MDL_YCOOR,y, __iter.objId);
                row.setValue(MDL_YCOOR,y);
                vd.push_back(row);
            }
        }
//        md->metadataToVec(vd);
        picker->train(vd, false, x, y, width, height);
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
        std::vector<MDRow> vd;
        jboolean aux=false;
        const FileName micrograph = env->GetStringUTFChars(filename, &aux);
        picker->automaticallySelectParticles(micrograph, percent, vd);
        outputmd->vecToMetadata(vd);
        //     md.write(env->GetStringUTFChars(jautoposfile, false));
    }
    XMIPP_JAVA_CATCH;
}


JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_correct
(JNIEnv *env, jobject jobj, jobject jmanualmd, jobject jautomaticmd, jdouble threshold)
{
    XMIPP_JAVA_TRY
    {

        std::vector<MDRow> vManualMd, vAutomaticMd;
        MetaData * manualmd = GET_INTERNAL_METADATA(jmanualmd);
        MetaData * automaticmd = GET_INTERNAL_METADATA(jautomaticmd);
        AutoParticlePicking2 *picker = GET_INTERNAL_AUTOPARTICLEPICKING2(jobj);
        manualmd->metadataToVec(vManualMd);
        automaticmd->metadataToVec(vAutomaticMd);
        picker->correction(vManualMd, vAutomaticMd);
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

