#include "xmipp_java_adapter.h"

#include <jni.h>
#include "xmipp_MetaData.h"
#include "xmipp_InternalData.h"
#include "xmipp_ExceptionsHandler.h"

#include <data/metadata.h>
#include <data/metadata_label.h>
#include <data/metadata_extension.h>
#include <data/xmipp_program_sql.h>
#include <classification/analyze_cluster.h>

int debug = 0;
/*
JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_storeIds
(JNIEnv *env, jclass cls) {
 MetaData_peerId = env->GetFieldID(cls, "peer", "J");
}
*/
JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_create
(JNIEnv *env, jobject jobj)
{
    MetaData * md = new MetaData();
    //env->SetLongField(jobj, MetaData_peerId, (long)md);
    STORE_PEER_ID(jobj, md);
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_destroy
(JNIEnv *env, jobject jobj)
{
    MetaData * md = GET_INTERNAL_METADATA(jobj);
    delete md;
    md = NULL;
    //env->SetLongField(jobj, MetaData_peerId, (long)md);
    STORE_PEER_ID(jobj, md);
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_clear
(JNIEnv *env, jobject jobj)
{
  XMIPP_JAVA_TRY
  {
      MetaData * md = GET_INTERNAL_METADATA(jobj);
      md->clear();
  }
  XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_read_1
(JNIEnv *env, jobject jobj, jstring filename)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);

        jboolean aux=false;
        const char * fn = env->GetStringUTFChars(filename, &aux);
        md->read(fn);
        env->ReleaseStringUTFChars(filename, fn);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT jint JNICALL Java_xmipp_jni_MetaData_size(JNIEnv *env, jobject jobj)
{
    XMIPP_JAVA_TRY
    {
      MetaData * md = GET_INTERNAL_METADATA(jobj);
        return md->size();
    }
    XMIPP_JAVA_CATCH;

    return 0;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_setColumnFormat(JNIEnv *env, jobject jobj, jboolean format)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);

        md->setColumnFormat(format);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_write
(JNIEnv *env, jobject jobj, jstring filename)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);

        jboolean aux=false;
        const char * fn = env->GetStringUTFChars(filename, &aux);
        md->write(fn);
        env->ReleaseStringUTFChars(filename, fn);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_writeBlock
(JNIEnv *env, jobject jobj, jstring filename)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);

        jboolean aux=false;
        const char * fn = env->GetStringUTFChars(filename, &aux);
        md->write(fn, MD_APPEND);
        env->ReleaseStringUTFChars(filename, fn);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_print
(JNIEnv *env, jobject jobj)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);


        md->write(std::cout);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_isColumnFormat
(JNIEnv *env, jobject jobj)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);


        return md->isColumnFormat();
    }
    XMIPP_JAVA_CATCH;

    return false;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_containsLabel
(JNIEnv *env, jobject jobj, jint label)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);


        return md->containsLabel((MDLabel) label);
    }
    XMIPP_JAVA_CATCH;

    return false;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_removeLabel
(JNIEnv *env, jobject jobj, jint label)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);


        return md->removeLabel((MDLabel) label);
    }
    XMIPP_JAVA_CATCH;

    return false;
}

JNIEXPORT jstring JNICALL Java_xmipp_jni_MetaData_label2Str(JNIEnv *env,
        jclass class_, jint label)
{
    XMIPP_JAVA_TRY
    {
        String str = MDL::label2Str((MDLabel) label);
        return env->NewStringUTF(str.c_str());
    }
    XMIPP_JAVA_CATCH;

    return NULL;
}

JNIEXPORT jint JNICALL Java_xmipp_jni_MetaData_str2Label
  (JNIEnv * env, jclass class_, jstring jlabelName)
{
    XMIPP_JAVA_TRY
    {
        jboolean aux=false;
        const char * labelName = env->GetStringUTFChars(jlabelName, &aux);
        MDLabel label = MDL::str2Label(labelName);
        env->ReleaseStringUTFChars(jlabelName, labelName);
        return (jint) label;
    }
    XMIPP_JAVA_CATCH;

    return -1;
}

JNIEXPORT jobjectArray JNICALL Java_xmipp_jni_MetaData_getBlocksInMetaDataFile(
    JNIEnv *env, jclass class_, jstring filename)
{
    XMIPP_JAVA_TRY
    {
        StringVector blocks;

        jboolean aux=false;
        const char * fn = env->GetStringUTFChars(filename, &aux);

        getBlocksInMetaDataFile(fn, blocks);

        env->ReleaseStringUTFChars(filename, fn);

        // Sets array value
        jstring str;
        jobjectArray array = env->NewObjectArray(blocks.size(), env->FindClass(
                                 "java/lang/String"), NULL);

        for (size_t i = 0; i < blocks.size(); i++)
        {
            str = env->NewStringUTF(blocks[i].c_str());
            env->SetObjectArrayElement(array, i, str);
            env->DeleteLocalRef(str);
        }

        return array;
    }
    XMIPP_JAVA_CATCH;

    return NULL;
}

JNIEXPORT jintArray JNICALL Java_xmipp_jni_MetaData_getActiveLabels(JNIEnv *env,
        jobject jobj)
{
    MetaData *md = GET_INTERNAL_METADATA(jobj);
    XMIPP_JAVA_TRY
    {
        std::vector < MDLabel > labels = md->getActiveLabels();

        // Copies vector into array.
        size_t size = labels.size();
        jint *body = new jint[size];
        for (size_t i = 0; i < size; i++)
        {
            body[i] = labels[i];
        }

        // Sets array value
        jintArray array = env->NewIntArray(size);
        env->SetIntArrayRegion(array, 0, size, body);

        return array;
    }
    XMIPP_JAVA_CATCH;

    return NULL;
}

JNIEXPORT jint JNICALL Java_xmipp_jni_MetaData_getLabelType(JNIEnv *env,
        jclass jclass_, jint label)
{
    XMIPP_JAVA_TRY
    {
        return (MDL::labelType((MDLabel) label));
    }
    XMIPP_JAVA_CATCH;

    return LABEL_NOTYPE;
}

JNIEXPORT jstring JNICALL Java_xmipp_jni_MetaData_getLabelComment(JNIEnv *env,
        jclass jclass_, jint label)
{
    XMIPP_JAVA_TRY
    {
        ProgramDb db;
        String str = db.getLabelComment((MDLabel) label);
        return env->NewStringUTF(str.c_str());
    }
    XMIPP_JAVA_CATCH;

    return NULL;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_isTextFile(JNIEnv *env,
        jclass class_, jint label)
{
    bool result = false;

    XMIPP_JAVA_TRY
    {
        result = MDL::isTextFile((MDLabel) label);
    }
    XMIPP_JAVA_CATCH;

    return result;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_isMetadata(JNIEnv *env,
        jclass class_, jint label)
{
    bool result = false;

    XMIPP_JAVA_TRY
    {
        result = MDL::isMetadata((MDLabel) label);
    }
    XMIPP_JAVA_CATCH;

    return result;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_isCtfParam(JNIEnv *env,
        jclass class_, jint label)
{
    bool result = false;

    XMIPP_JAVA_TRY
    {
        result = MDL::isCtfParam((MDLabel) label);
    }
    XMIPP_JAVA_CATCH;

    return result;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_isImage(JNIEnv *env,
        jclass class_, jint label)
{
    bool result = false;

    XMIPP_JAVA_TRY
    {
        result = MDL::isImage((MDLabel) label);
    }
    XMIPP_JAVA_CATCH;

    return result;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_isStack(JNIEnv *env,
        jclass class_, jint label)
{
    bool result = false;

    XMIPP_JAVA_TRY
    {
        result = MDL::isStack((MDLabel) label);
    }
    XMIPP_JAVA_CATCH;

    return result;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_isMicrograph(JNIEnv *env,
        jclass class_, jint label)
{
    bool result = false;

    XMIPP_JAVA_TRY
    {
        result = MDL::isMicrograph((MDLabel) label);
    }
    XMIPP_JAVA_CATCH;

    return result;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_isPSD(JNIEnv *env,
        jclass class_, jint label)
{
    bool result = false;

    XMIPP_JAVA_TRY
    {
        result = MDL::isPSD((MDLabel) label);
    }
    XMIPP_JAVA_CATCH;

    return result;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_isMetadataFile(JNIEnv *env,
		jobject jobj)
{
    bool result = false;

    XMIPP_JAVA_TRY
    {
    	 MetaData * md = GET_INTERNAL_METADATA(jobj);
        result =  md->isMetadataFile;
    }
    XMIPP_JAVA_CATCH;

    return result;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_makeAbsPath
  (JNIEnv * env, jobject jobj, jint label)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);


        md->makeAbsPath((MDLabel) label);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT jint JNICALL Java_xmipp_jni_MetaData_getValueInt(JNIEnv *env,
        jobject jobj, jint label, jlong objId)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);


        int value;

        if (md->getValue((MDLabel) label, value, objId))
            return value;
    }
    XMIPP_JAVA_CATCH;

    return 0;
}

JNIEXPORT jlong JNICALL Java_xmipp_jni_MetaData_getValueLong
(JNIEnv *env, jobject jobj, jint label, jlong objId)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);


        size_t value;
        if (md->getValue((MDLabel) label, value, objId))
            return value;
    }
    XMIPP_JAVA_CATCH;

    return 0;
}

JNIEXPORT jdouble JNICALL Java_xmipp_jni_MetaData_getValueDouble(JNIEnv *env,
        jobject jobj, jint label, jlong objId)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);

        double value;
        if (md->getValue((MDLabel) label, value, objId))
            return value;
    }
    XMIPP_JAVA_CATCH;

    return 0;
}

JNIEXPORT jstring JNICALL Java_xmipp_jni_MetaData_getValueString(JNIEnv *env,
        jobject jobj, jint label, jlong objId)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);

        MDObject obj((MDLabel) label);
        if (md->getValue(obj, objId))
        {
          //String str("constant.kkkkk");
            jstring jstr = env->NewStringUTF(obj.toString().c_str());
            //env->DeleteLocalRef(jstr);
            return jstr;
        }
    }
    XMIPP_JAVA_CATCH;

    return NULL;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_getValueBoolean(JNIEnv *env,
        jobject jobj, jint label, jlong objId)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);


        bool value;
        if (md->getValue((MDLabel) label, value, objId))
            return value;
    }
    XMIPP_JAVA_CATCH;

    return 0;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_setValueInt(JNIEnv *env,
        jobject jobj, jint label, jint value, jlong objId)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);


        return md->setValue((MDLabel) label, (int)value, objId);
    }
    XMIPP_JAVA_CATCH;

    return false;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_setValueLong(JNIEnv *env,
        jobject jobj, jint label, jlong value, jlong objId)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);


        return md->setValue((MDLabel) label, (size_t)value, objId);
    }
    XMIPP_JAVA_CATCH;

    return false;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_setValueDouble(JNIEnv *env,
        jobject jobj, jint label, jdouble value, jlong objId)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);


        return md->setValue((MDLabel) label, value, objId);
    }
    XMIPP_JAVA_CATCH;

    return false;
}
JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_setValueString(JNIEnv *env,
        jobject jobj, jint label, jstring value, jlong objId)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);

//        std::cout <<  (MDLabel) label << ": " << MDL_IMAGE << " -> " << MDL::label2Str(MDL_IMAGE) << std::endl;

        jboolean aux=false;
        const char * strValue = env->GetStringUTFChars(value, &aux);
        bool result = md->setValueFromStr((MDLabel) label, strValue, objId);

        env->ReleaseStringUTFChars(value, strValue);

        return result;
    }
    XMIPP_JAVA_CATCH;

    return false;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_setValueBoolean(JNIEnv *env,
        jobject jobj, jint label, jboolean value, jlong objId)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);


        return md->setValue((MDLabel) label, (bool) value, objId);
    }
    XMIPP_JAVA_CATCH;

    return false;
}

JNIEXPORT jdoubleArray JNICALL Java_xmipp_jni_MetaData_getStatistics(JNIEnv *env,
        jobject jobj, jboolean applyGeo)
{
    MetaData *md = GET_INTERNAL_METADATA(jobj);

    XMIPP_JAVA_TRY
    {
        double avg, stddev, min, max;
        getStatistics((*md), avg, stddev, min, max, applyGeo);

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

JNIEXPORT jdoubleArray JNICALL Java_xmipp_jni_MetaData_getColumnValues(JNIEnv *env,
        jobject jobj, jint label)
{
    MetaData *md = GET_INTERNAL_METADATA(jobj);

    XMIPP_JAVA_TRY
    {
        std::vector<double> values;
        md->getColumnValues((MDLabel) label, values);

        // Copies vector into array.
        size_t size = values.size();
        jdouble *body = new jdouble[size];
        for (size_t i = 0; i < size; i++)
        {
            body[i] = values[i];
        }

        jdoubleArray array = env->NewDoubleArray(size);
        env->SetDoubleArrayRegion(array, 0, size, body);

        return array;
    }
    XMIPP_JAVA_CATCH;

    return NULL;
}

//Utility function to create jlongArray from std::vector<size_t>
jlongArray createLongArray(JNIEnv *env, const std::vector<size_t> & ids){
    // Copies vector into array.
    size_t size = ids.size();
    jlong *body = new jlong[size];
    for (size_t i = 0; i < size; i++)
    {
        body[i] = ids[i];
    }

    // Sets array value
    jlongArray array = env->NewLongArray(size);
    env->SetLongArrayRegion(array, 0, size, body);

    return array;
}

JNIEXPORT jlongArray JNICALL Java_xmipp_jni_MetaData_findObjects(JNIEnv *env,
        jobject jobj)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);


        std::vector < size_t > ids;
        md->findObjects(ids);
        return createLongArray(env, ids);
    }
    XMIPP_JAVA_CATCH;

    return NULL;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_sort
  (JNIEnv * env, jobject jobj, jint label, jboolean ascending){

  XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);


        MetaData mdSorted;
        mdSorted.sort(*md, (MDLabel)label, ascending);
        *md = mdSorted;
    }
    XMIPP_JAVA_CATCH;
}


JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_importObjects
(JNIEnv *env, jobject jobj, jobject from, jlongArray jids)
{
    MetaData * md = GET_INTERNAL_METADATA(jobj);
    MetaData * mdfrom = GET_INTERNAL_METADATA(from);

    XMIPP_JAVA_TRY
    {
        jlong *ids = env->GetLongArrayElements(jids, 0);
        int size = env->GetArrayLength(jids);

        std::vector<size_t> out_ids(size);
        for(int i=0; i<size; i++)
        {
            out_ids[i] = ids[i];
        }

        env->ReleaseLongArrayElements(jids, ids, 0);

        md->importObjects(*mdfrom, out_ids);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT jlong JNICALL Java_xmipp_jni_MetaData_firstObject(JNIEnv *env,
        jobject jobj)
{
    MetaData * md = GET_INTERNAL_METADATA(jobj);
    jlong id = 0;

    XMIPP_JAVA_TRY
    {
        id = md->firstObject();
    }
    XMIPP_JAVA_CATCH;

    return id;
}

JNIEXPORT jlong JNICALL Java_xmipp_jni_MetaData_addObject(JNIEnv *env, jobject jobj)
{
    jlong id = 0;
    MetaData *md = GET_INTERNAL_METADATA(jobj);

    XMIPP_JAVA_TRY
    {
        id = md->addObject();
    }
    XMIPP_JAVA_CATCH;

    return id;
}

JNIEXPORT jboolean JNICALL Java_xmipp_jni_MetaData_removeObject
  (JNIEnv * env, jobject jobj, jlong objId)
{
  MetaData *md = GET_INTERNAL_METADATA(jobj);

  XMIPP_JAVA_TRY
  {
    return md->removeObject((size_t)objId);
  }
  XMIPP_JAVA_CATCH;

  return false;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_removeDisabled
  (JNIEnv * env, jobject jobj)
{
  MetaData *md = GET_INTERNAL_METADATA(jobj);

  XMIPP_JAVA_TRY
  {
    return md->removeDisabled();
  }
  XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_addLabel(JNIEnv *env, jobject jobj, jint label)
{
    MetaData *md = GET_INTERNAL_METADATA(jobj);

    XMIPP_JAVA_TRY
    {
        md->addLabel((MDLabel)label);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_getStatsImages
(JNIEnv *env, jobject jmetadata,
		jobject jimageAvg, jobject jimageStd, jboolean applyGeo, jint label)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jmetadata);
        ImageGeneric *avg = GET_INTERNAL_IMAGE_GENERIC(jimageAvg);
        ImageGeneric *std = GET_INTERNAL_IMAGE_GENERIC(jimageStd);
        avg->setDatatype(DT_Double);
        std->setDatatype(DT_Double);
        Image<double> * imgAvg = (Image<double>*)avg->image;
        Image<double> * imgStd = (Image<double>*)std->image;
        double dum;
        getStatistics(*md, *imgAvg, *imgStd, dum, dum, applyGeo, (MDLabel)label);
    }
    XMIPP_JAVA_CATCH;
}


JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_getPCAbasis
(JNIEnv *env, jobject jmetadata, jobject jbasis, jint label)
{
    XMIPP_JAVA_TRY
    {
        MetaData * mdIn = GET_INTERNAL_METADATA(jmetadata);
        ImageGeneric *basis= GET_INTERNAL_IMAGE_GENERIC(jbasis);
        basis->setDatatype(DT_Double);
        MultidimArray<double> *mdArray;
        MULTIDIM_ARRAY_GENERIC(*basis).getMultidimArrayPointer(mdArray);

        ProgAnalyzeCluster program;
        program.verbose = 1;
        program.NPCA=4;
        program.Niter=10;
        program.dontMask=false;
        program.SFin=*mdIn;
        program.produceSideInfo((MDLabel)label);
        program.pcaAnalyzer.evaluateZScore(program.NPCA, program.Niter);
        program.produceBasis(*mdArray);
        //Notify the real dimensions to image generic
        ArrayDim adim;
        mdArray->getDimensions(adim);
        basis->image->setADimFile(adim);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_computeFourierStatistics
(JNIEnv *env, jobject jobj, jobject jmetadata, jint label)
{
    MetaData * mdOut = GET_INTERNAL_METADATA(jobj);

    XMIPP_JAVA_TRY
    {
        MetaData * mdIn = GET_INTERNAL_METADATA(jmetadata);
        getFourierStatistics(*mdIn, 1, *mdOut, true, 2, (MDLabel)label);
    }
    XMIPP_JAVA_CATCH;
}

/*
 * Class:     xmipp_MetaData
 * Method:    unionAll
 * Signature: (Ljava/lang/String;)V
 */
JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_unionAll
  (JNIEnv * env, jobject jobj, jobject jmdIn)
{
  XMIPP_JAVA_TRY
  {
      MetaData * md = GET_INTERNAL_METADATA(jobj);


      MetaData * mdIn = GET_INTERNAL_METADATA(jmdIn);
      md->unionAll(*mdIn);
  }
  XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_fillConstant
  (JNIEnv * env, jobject jobj, jint label, jstring value)
{
  XMIPP_JAVA_TRY
  {
      MetaData * md = GET_INTERNAL_METADATA(jobj);

      jboolean aux=false;
      const char * strValue = env->GetStringUTFChars(value, &aux);
      md->fillConstant((MDLabel)label, strValue);
  }
  XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_fillLinear
  (JNIEnv * env, jobject jobj, jint label, jdouble start, jdouble step)
{
  XMIPP_JAVA_TRY
  {
      MetaData * md = GET_INTERNAL_METADATA(jobj);


      md->fillLinear((MDLabel)label, start, step);
  }
  XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_fillRandom
  (JNIEnv * env, jobject jobj, jint label, jstring mode, jdouble op1, jdouble op2)
{
  XMIPP_JAVA_TRY
  {
      MetaData * md = GET_INTERNAL_METADATA(jobj);

      jboolean aux=false;
      const char * strMode = env->GetStringUTFChars(mode, &aux);
      md->fillRandom((MDLabel)label, strMode, op1, op2);
  }
  XMIPP_JAVA_CATCH;
}


JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_enableDebug
(JNIEnv *, jobject)
{
    extern int debug;
    debug = 1;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_readPlain
(JNIEnv * env, jobject jobj, jstring jfile, jstring jcolumns)
{
    XMIPP_JAVA_TRY
    {
        jboolean aux=false;
        const char * nfile = env->GetStringUTFChars(jfile, &aux);
        aux=false;
        const char * ncolumns = env->GetStringUTFChars(jcolumns, &aux);
        MetaData * metadata = GET_INTERNAL_METADATA(jobj);
        metadata->readPlain(nfile, ncolumns);
        env->ReleaseStringUTFChars(jfile, nfile);
        env->ReleaseStringUTFChars(jcolumns, ncolumns);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_writeImages
(JNIEnv * env, jobject jobj, jstring joutput, jboolean independent, jint image_label)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);
        jboolean aux=false;
        const char * output = env->GetStringUTFChars(joutput, &aux);
        copyImages(*md, output, independent, (MDLabel) image_label);
        env->ReleaseStringUTFChars(joutput, output);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_operate
(JNIEnv * env, jobject jobj, jstring operateString)
{
    XMIPP_JAVA_TRY
    {
        jboolean aux=false;
        const char * opStr = env->GetStringUTFChars(operateString, &aux);
        MetaData * md = GET_INTERNAL_METADATA(jobj);
        md->operate(opStr);
    }
    XMIPP_JAVA_CATCH;
}

JNIEXPORT jdouble JNICALL Java_xmipp_jni_MetaData_getColumnMax
(JNIEnv * env, jobject jobj, jint column)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);
        return md->getColumnMax((MDLabel)column);
    }
    XMIPP_JAVA_CATCH;
    return 0.;
}

JNIEXPORT jdouble JNICALL Java_xmipp_jni_MetaData_getColumnMin
(JNIEnv * env, jobject jobj, jint column)
{
    XMIPP_JAVA_TRY
    {
        MetaData * md = GET_INTERNAL_METADATA(jobj);
        return md->getColumnMin((MDLabel)column);
    }
    XMIPP_JAVA_CATCH;
    return 0.;
}

JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_getRow
(JNIEnv * env, jobject jobj, jobject jobjRow, jlong objId)
{
  XMIPP_JAVA_TRY
  {
      MetaData * md = GET_INTERNAL_METADATA(jobj);
      MDRow * mdRow = GET_INTERNAL_MDROW(jobjRow);
      md->getRow(*mdRow, (size_t)objId);
  }
  XMIPP_JAVA_CATCH;
}


JNIEXPORT void JNICALL Java_xmipp_jni_MetaData_setRow
(JNIEnv * env, jobject jobj, jobject jobjRow, jlong objId)
{
  XMIPP_JAVA_TRY
  {
      MetaData * md = GET_INTERNAL_METADATA(jobj);
      MDRow * mdRow = GET_INTERNAL_MDROW(jobjRow);
      md->setRow(*mdRow, (size_t)objId);
  }
  XMIPP_JAVA_CATCH;
}


