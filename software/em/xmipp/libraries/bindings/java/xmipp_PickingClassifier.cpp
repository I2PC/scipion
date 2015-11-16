#include "xmipp_java_adapter.h"
#include <iostream>
#include "xmipp_PickingClassifier.h"
#include "xmipp_ExceptionsHandler.h"
#include "reconstruction/micrograph_automatic_picking2.h"
#include "xmipp_InternalData.h"
#include <data/metadata.h>

JNIEXPORT void JNICALL
Java_xmipp_jni_PickingClassifier_create(JNIEnv *env, jobject jobj, jint particle_size, jstring output, jobjectArray jmics)
{
    XMIPP_JAVA_TRY
    {
        int size = particle_size, filter_num = 6, corr_num = 2, NPCA = 4;
        jboolean aux=false;
        const FileName &model_name = env->GetStringUTFChars(output, &aux);
        int length = env->GetArrayLength(jmics);
	    MDRow row;
	    std::vector<MDRow> mics;
	    for(int i = 0; i < length; i++)
	    {
			 row = *GET_INTERNAL_MDROW(env->GetObjectArrayElement(jmics, i));
			 mics.push_back(row);
	    }
        AutoParticlePicking2 *picker = new AutoParticlePicking2(size, filter_num, corr_num, NPCA, model_name, mics);

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
(JNIEnv *env, jobject jobj, jobjectArray jtrainList, jint x, jint y, jint width, jint height)
{
    XMIPP_JAVA_TRY
    {
    	   AutoParticlePicking2 *picker = GET_INTERNAL_AUTOPARTICLEPICKING2(jobj);
    	   MetaData md;
    	   FileName micFile, posFile;
    	   std::vector<MDRow> vd;
    	   MDRow row;
    	   MDRow* trainRow;
    	   int size = env->GetArrayLength(jtrainList);
    	   int xcoor, ycoor;
    	   for(int i = 0; i < size; i++)
    	   {
				trainRow = GET_INTERNAL_MDROW(env->GetObjectArrayElement(jtrainList, i));
				trainRow->getValue(MDL_MICROGRAPH, micFile);
				trainRow->getValue(MDL_MICROGRAPH_PARTICLES, posFile);
				md.read("particles@" + posFile);
				FOR_ALL_OBJECTS_IN_METADATA(md)
				{
					row.setValue(MDL_MICROGRAPH, micFile);
					md.getValue(MDL_XCOOR, xcoor, __iter.objId);
					row.setValue(MDL_XCOOR, xcoor);
					md.getValue(MDL_YCOOR, ycoor, __iter.objId);
					row.setValue(MDL_YCOOR, ycoor);
					vd.push_back(row);
				}
    	   }

        picker->train(vd, false, x, y, width, height);
    }
    XMIPP_JAVA_CATCH;

}


JNIEXPORT jobjectArray JNICALL Java_xmipp_jni_PickingClassifier_autopick
(JNIEnv *env, jobject jobj, jstring filename, jint percent)
{
    XMIPP_JAVA_TRY
    {
        AutoParticlePicking2 *picker = GET_INTERNAL_AUTOPARTICLEPICKING2(jobj);
        std::vector<MDRow> vd;
        jboolean aux=false;
        const FileName micrograph = env->GetStringUTFChars(filename, &aux);
        picker->automaticallySelectParticles(micrograph, percent, vd);

        jclass pclass = env->FindClass("xmipp/jni/Particle");
        jmethodID constructor = env->GetMethodID(pclass, "<init>", "(IID)V");
        jobject particle;
		jobjectArray rows = env->NewObjectArray(vd.size(), pclass, 0);

		MDRow* row;
		int x, y;
		double cost;
		for (int i = 0; i < vd.size(); i++) {
			row = &vd[i];
			row->getValue(MDL_XCOOR, x);
			row->getValue(MDL_YCOOR, y);
			row->getValue(MDL_COST, cost);
			particle = env->NewObject(pclass, constructor, x, y, cost);

			env->SetObjectArrayElement(rows, i, particle);
		}
		return rows;
    }
    XMIPP_JAVA_CATCH;
}


JNIEXPORT void JNICALL Java_xmipp_jni_PickingClassifier_correct
(JNIEnv *env, jobject jobj, jobjectArray jmanualRows, jobjectArray jautoRows)
{
    XMIPP_JAVA_TRY
    {
    	  int size = env->GetArrayLength(jmanualRows);
    	  MDRow row;
    	  std::vector<MDRow> manualRows;
    	  for(int i = 0; i < size; i++)
    	  {
    		  	 row = *GET_INTERNAL_MDROW(env->GetObjectArrayElement(jmanualRows, i));
    	         manualRows.push_back(row);
    	  }

    	  size = env->GetArrayLength(jautoRows);
		  std::vector<MDRow> autoRows;
		  for(int i = 0; i < size; i++)
		  {
				 row = *GET_INTERNAL_MDROW(env->GetObjectArrayElement(jautoRows, i));
				 autoRows.push_back(row);
		  }


//        std::vector<MDRow> vManualMd, vAutomaticMd;
//        MetaData * manualmd = GET_INTERNAL_METADATA(jmanualmd);
//        MetaData * automaticmd = GET_INTERNAL_METADATA(jautomaticmd);
        AutoParticlePicking2 *picker = GET_INTERNAL_AUTOPARTICLEPICKING2(jobj);
        //manualmd->metadataToVec(vManualMd);
//        automaticmd->metadataToVec(vAutomaticMd);
        picker->correction(manualRows, autoRows);
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


JNIEXPORT jint JNICALL Java_xmipp_jni_PickingClassifier_getParticlesThreshold
(JNIEnv *env, jobject jobj)
{
    XMIPP_JAVA_TRY
    {
        AutoParticlePicking2 *picker = GET_INTERNAL_AUTOPARTICLEPICKING2(jobj);
        return picker->getParticlesThreshold();

    }
    XMIPP_JAVA_CATCH;
    return 30;
}
