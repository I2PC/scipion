/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2011)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
#include <math.h>
#include "micrograph_automatic_picking2.h"
#include <data/filters.h>
#include <data/rotational_spectrum.h>
#include <reconstruction/denoise.h>
#include <data/xmipp_fft.h>
#include <data/xmipp_filename.h>
#include <algorithm>


AutoParticlePicking2::AutoParticlePicking2(const FileName &fn, Micrograph *_m,int size,int filterNum,int pcaNum,int corrNum)
{
    __m = _m;
    fn_micrograph = fn;
    microImage.read(fn_micrograph);
    scaleRate=std::max(0.25,50.0/size);
    particle_radius = (size*scaleRate)*0.5;
    particle_size = particle_radius*2;
    NRsteps=particle_size/2-3;
    filter_num = filterNum;
    corr_num = corrNum;
    NPCA = pcaNum;
    __numThreads = 1;
    piece_xsize = 512;
}

//Generate filter bank from the micrograph image
void filterBankGenerator(MultidimArray<double> &inputMicrograph, const FileName &fnFilterBankStack, int filter_num)
{
    fnFilterBankStack.deleteFile();

    Image<double> Iaux;
    Iaux()=inputMicrograph;

    FourierFilter filter;
    filter.raised_w = 0.02;
    filter.FilterShape = RAISED_COSINE;
    filter.FilterBand = BANDPASS;
    MultidimArray< std::complex<double> > micrographFourier;
    FourierTransformer transformer;
    transformer.FourierTransform(Iaux(), micrographFourier, true);

    for (int i=0; i<filter_num; i++)
    {
        filter.w1 =0.025*i;
        filter.w2 = (filter.w1)+0.025;

        transformer.setFourier(micrographFourier);
        filter.applyMaskFourierSpace(inputMicrograph,transformer.fFourier);
        transformer.inverseFourierTransform();

        Iaux.write(fnFilterBankStack,i+1,true,WRITE_APPEND);
    }
}
void correlationBetweenPolarChannels(int n1, int n2, int nF,
                                     MultidimArray<double> &mIpolar, MultidimArray<double> &mIpolarCorr,
                                     CorrelationAux &aux)
{
    MultidimArray<double> imgPolar1, imgPolar2, imgPolarCorr, corr2D;
    imgPolar1.aliasImageInStack(mIpolar,n1);
    imgPolar2.aliasImageInStack(mIpolar,n2);

    imgPolarCorr.aliasImageInStack(mIpolarCorr,nF);
    correlation_matrix(imgPolar1,imgPolar2,corr2D,aux);
    imgPolarCorr=corr2D;
}

void stackReadResize(int n,int z,int y,int x,ImageGeneric &imageStack,
                     MultidimArray<float> &vec)
{
    MultidimArray<float> * pTempArray;
    MultidimArray<float> tempArray;

    imageStack().getMultidimArrayPointer(pTempArray);
    imageStack().getImage(vec);
    typeCast(*pTempArray,vec);
    vec.resize(n,z,y,x);
}
SVMClassifier::SVMClassifier()
{
    param.svm_type = C_SVC;
    param.kernel_type = RBF;
    param.degree = 3;
    param.gamma = 0;
    param.coef0 = 0;
    param.nu = 0.5;
    param.cache_size = 100;
    param.C = 1;
    param.eps = 1e-3;
    param.p = 0.1;
    param.shrinking = 1;
    param.probability = 0;
    param.nr_weight = 0;
    param.weight_label = NULL;
    param.weight = NULL;
}
void SVMClassifier::SVMTrain(MultidimArray<double> &trainSet,MultidimArray<int> &lable)
{
    prob.l = YSIZE(trainSet);
    prob.y = new double[prob.l];
    prob.x = new svm_node *[prob.l];
    for (int i=0;i<YSIZE(trainSet);i++)
    {
        prob.x[i]=new svm_node[YSIZE(trainSet)];
        int cnt = 0;
        for (int j=0;j<XSIZE(trainSet);j++)
        {
            if (trainSet(i,j)==0)
                continue;
            else
            {
                prob.x[i][cnt].value=DIRECT_A2D_ELEM(trainSet,i,j);
                prob.x[i][cnt].index=j+1;
                cnt++;
            }
        }
        prob.x[i][cnt].index=-1;
        prob.y[i] = DIRECT_A1D_ELEM(lable,i);
    }
    model=svm_train(&prob, &param);
}
void  AutoParticlePicking2::polarCorrelation(MultidimArray<double> &Ipolar,MultidimArray<double> &IpolarCorr)
{
    int nF=NSIZE(Ipolar);
    CorrelationAux aux;

    for (int n=0; n<nF; ++n)
    {
        correlationBetweenPolarChannels(n,n,n,Ipolar,IpolarCorr,aux);
    }
    for (int i=0;i<(filter_num-corr_num);i++)
        for (int j=1;j<=corr_num;j++)
        {
            correlationBetweenPolarChannels(i,i+j,nF++,Ipolar,IpolarCorr,aux);
        }
}

AutoParticlePicking2::~AutoParticlePicking2()
{}

void AutoParticlePicking2::buildInvariant(MultidimArray<double> &invariantChannel,int x,int y,
        const FileName &fnFilterBank)
{
    ImageGeneric micrographStack;
    MultidimArray<double>  pieceImage;
    MultidimArray<double> Ipolar;
    MultidimArray<double> mIpolar;
    Ipolar.initZeros(filter_num,1,NangSteps,NRsteps);

    for (int j=0;j<filter_num;j++)
    {
        micrographStack.readMapped(fnFilterBank,j+1);
        extractParticle(x,y,micrographStack,pieceImage);
        mIpolar.aliasImageInStack(Ipolar,j);
        convert2Polar(pieceImage,mIpolar);
    }
    polarCorrelation(Ipolar,invariantChannel);
}
//void AutoParticlePicking2::PCAProject(MultidimArray<double>&inputVec,int firstIndex)
//{
//    inputVec-=pcaAnalyzer.avg;
//
//    for (int i=0; i<NPCA; i++)
//    {
//        double dotProduct=0;
//        for (int j=0;j<)
//        {
//            dotProduct += *ptrii * *ptrjj;
//            dotProduct += *(ptrii+1) * *(ptrjj+1);
//            dotProduct += *(ptrii+2) * *(ptrjj+2);
//            dotProduct += *(ptrii+3) * *(ptrjj+3);
//        }
//    }
//}
void AutoParticlePicking2::trainPCA(const FileName &fnPositiveFeat)
{
    ImageGeneric positiveInvariant;
    MultidimArray<float> pcaVec;
    ArrayDim aDim;

    FileName fnPositiveInvariatn=fnPositiveFeat + "_invariant_Positive.stk";
    positiveInvariant.read(fnPositiveInvariatn, HEADER);
    positiveInvariant.getDimensions(aDim);
    int num_correlation=filter_num+((filter_num-corr_num)*corr_num);
    int steps= aDim.ndim/num_correlation;

    for (int i=1;i<=num_correlation;i++)
    {
        for (int j=0;j<steps;j++)
        {
            positiveInvariant.readMapped(fnPositiveInvariatn,j*num_correlation+i);
            positiveInvariant().getImage(pcaVec);
            pcaVec.resize(1,1,1,XSIZE(pcaVec)*YSIZE(pcaVec));
            pcaAnalyzer.addVector(pcaVec);
        }
        pcaAnalyzer.subtractAvg();
        pcaAnalyzer.learnPCABasis(NPCA, 10);
        savePCAModel(fnPositiveFeat);
        pcaAnalyzer.clear();
    }
}
void AutoParticlePicking2::trainSVM(const FileName &fn_root)
{
    FileName fnPCAModel=fn_root + "_pca_model.stk";
    FileName fnPositiveInvariant=fn_root + "_invariant_Positive.stk";
    ImageGeneric positiveInvariant;
    ImageGeneric pcaBasis;
    MultidimArray<double> vec;
    MultidimArray<double> trainSet;
    ArrayDim aDim;

    positiveInvariant.read(fnPositiveInvariant, HEADER);
    int num_correlation=filter_num+((filter_num-corr_num)*corr_num);
    int steps= aDim.ndim/num_correlation;
    //trainSet.resize(1,1,steps,num_correlation*NPCA);

    for (int i=0;i<num_correlation;i++)
    {
        pcaAnalyzer.clear();
        pcaBasis.readMapped(fnPCAModel,i*(NPCA+1)+1);
        pcaBasis().getImage(pcaAnalyzer.avg);
        pcaAnalyzer.avg.resize(1,1,1,XSIZE(pcaAnalyzer.avg)*YSIZE(pcaAnalyzer.avg));
        for (int j=0;j<NPCA;j++)
        {

            pcaBasis.readMapped(fnPCAModel,(i*(NPCA+1)+1)+(j+1));
            MultidimArray<double> Ijj;

            pcaBasis().getImage(Ijj);
            typeCast(Ijj,pcaAnalyzer.PCAbasis[j]);
            pcaAnalyzer.PCAbasis[j].resize(1,1,1,XSIZE(pcaAnalyzer.PCAbasis[j]),YSIZE(pcaAnalyzer.PCAbasis[j]));
        }
        for (int k=0;k<steps;k++)
        {
            positiveInvariant.readMapped(fnPositiveInvariant,k*num_correlation+i+1);
            positiveInvariant().getImage(vec);
            vec.resize(1,1,1,XSIZE(vec)*YSIZE(vec));
        }
        //fnPositiveInvariant.readMapped(fnPositiveInvariatn,j*num_correlation+i);
    }
}

void AutoParticlePicking2::extractPositiveInvariant(const FileName &fnFilterBank, const FileName &fnInvariantFeat)
{

    MultidimArray<double> IpolarCorr;

    int num_part = __m->ParticleNo();
    int num_correlation=filter_num+((filter_num-corr_num)*corr_num);
    FileName fnPositiveInvariatn=fnInvariantFeat+"_Positive.stk";
    fnPositiveInvariatn.deleteFile();
    IpolarCorr.initZeros(num_correlation,1,NangSteps,NRsteps);
    for (int i=0;i<num_part;i++)
    {
        int x = (__m->coord(i).X)*scaleRate;
        int y = (__m->coord(i).Y)*scaleRate;
        buildInvariant(IpolarCorr,x,y,fnFilterBank);
        Image<double> II;
        II() = IpolarCorr;
        II.write(fnPositiveInvariatn,ALL_IMAGES,true,WRITE_APPEND);
    }
}
void AutoParticlePicking2::extractInvariant(const FileName &fnFilterBank,const FileName &fnInvariantFeat)
{
    extractPositiveInvariant(fnFilterBank,fnInvariantFeat);
}

void AutoParticlePicking2::extractParticle(const int x, const int y,ImageGeneric & filter,
        MultidimArray<double> &particleImage)
{
    int startX, startY, endX, endY;
    MultidimArray<double> Filter;

    startX=x-particle_radius;
    startY=y-particle_radius;
    endX=x+particle_radius;
    endY=y+particle_radius;

    filter().getImage(Filter);
    Filter.window(particleImage,startY,startX,endY,endX);
    normalize_OldXmipp(particleImage);
}
void AutoParticlePicking2::convert2Polar(MultidimArray<double> &particleImage, MultidimArray<double> &polar)
{
    Matrix1D<double> R;
    particleImage.setXmippOrigin();
    image_convertCartesianToPolar_ZoomAtCenter(particleImage,polar,
            R,1,3,XSIZE(particleImage)/2,NRsteps,0,2*PI,NangSteps);
}

void AutoParticlePicking2::buildVector(MultidimArray<double> &trainSet)
{}

void AutoParticlePicking2::savePCAModel(const FileName &fn_root) const
{
    FileName pcaModel=fn_root + "_pca_model.stk";
    Image<double> II;
    MultidimArray<double> avg;
    avg=pcaAnalyzer.avg;
    avg.resize(1,1,NangSteps,NRsteps);
    II()=avg;
    II.write(pcaModel,ALL_IMAGES,true,WRITE_APPEND);
    for (int i=0;i<NPCA;i++)
    {
        MultidimArray<double> pcaBasis;
        pcaBasis=pcaAnalyzer.PCAbasis[i];
        pcaBasis.resize(1,1,NangSteps,NRsteps);
        II()=pcaBasis;
        II.write(pcaModel,ALL_IMAGES,true,WRITE_APPEND);
    }
}



//void AutoParticlePicking2::loadPCAModel(const FileName &fn_root)
//{
//    FileName pcaModel=fn_root + "_pca_model.txt";
//    std::ifstream fh_pca;
//    fh_pca.open(pcaModel.c_str());
//    if (!fh_pca)
//        REPORT_ERROR(ERR_IO_NOWRITE,fn_root);
//    fh_pca >> NPCA;
//    pcaAnalyzer.avg.resize(NangSteps*NRsteps);
//    fh_pca >> (pcaAnalyzer.avg);
//    std::cerr<<pcaAnalyzer.avg;
//
//
//    fh_pca.close();
//}

// ==========================================================================
// Section: Program interface ===============================================
// ==========================================================================
void ProgMicrographAutomaticPicking2::readParams()
{
    fn_micrograph = getParam("-i");
    fn_model = getParam("--model");
    size = getIntParam("--particleSize");
    mode = getParam("--mode");
    if (mode == "buildinv")
    {
        fn_train = getParam("--mode", 1);
    }
    NPCA = getIntParam("--NPCA");
    filter_num = getIntParam("--filter_num");
    corr_num = getIntParam("--NCORR");
    Nthreads = getIntParam("--thr");
    fn_root = getParam("--outputRoot");
    fast = checkParam("--fast");
    incore = checkParam("--in_core");
}

void ProgMicrographAutomaticPicking2::defineParams()
{
    addUsageLine("Automatic particle picking for micrographs");
    addUsageLine("+The algorithm is designed to learn the particles from the user, as well as from its own errors.");
    addUsageLine("+The algorithm is fully described in [[http://www.ncbi.nlm.nih.gov/pubmed/19555764][this paper]].");
    addParamsLine("  -i <micrograph>               : Micrograph image");
    addParamsLine("  --outputRoot <rootname>       : Output rootname");
    addParamsLine("  --mode <mode>                 : Operation mode");
    addParamsLine("         where <mode>");
    addParamsLine("                    try              : Try to autoselect within the training phase.");
    addParamsLine("                    train <positivefile=\"\"> <negativefile=\"\"> : posfile contains the coordinates of manually picked particles");
    addParamsLine("                                     : <rootname>_auto_feature_vectors.txt contains the particle structure created by this program when used in automatic selection mode");
    addParamsLine("                                     : <rootname>_false_positives.xmd contains the list of false positives among the automatically picked particles");
    addParamsLine("                    autoselect  : Autoselect");
    addParamsLine("                    buildinv <posfile=\"\"> : posfile contains the coordinates of manually picked particles");
    addParamsLine("  --model <model_rootname>      : Bayesian model of the particles to pick");
    addParamsLine("  --particleSize <size>         : Particle size in pixels");
    addParamsLine("  [--thr <p=1>]                 : Number of threads for automatic picking");
    addParamsLine("  [--fast]                      : Perform a fast preprocessing of the micrograph (Fourier filter instead of Wavelet filter)");
    addParamsLine("  [--in_core]                   : Read the micrograph in memory");
    addParamsLine("  [--filter_num <n=7>]          : The number of filters in filter bank");
    addParamsLine("  [--NPCA <n=4>]               : The number of PCA components");
    addParamsLine("  [--NCORR <n=2>]               : The number of PCA components");
    addExampleLine("Automatically select particles during training:", false);
    addExampleLine("xmipp_micrograph_automatic_picking -i micrograph.tif --particleSize 100 --model model --thr 4 --outputRoot micrograph --mode try ");
    addExampleLine("Training:", false);
    addExampleLine("xmipp_micrograph_automatic_picking -i micrograph.tif --particleSize 100 --model model --thr 4 --outputRoot micrograph --mode train manual.pos");
    addExampleLine("Automatically select particles after training:", false);
    addExampleLine("xmipp_micrograph_automatic_picking -i micrograph.tif --particleSize 100 --model model --thr 4 --outputRoot micrograph --mode autoselect");
}

void ProgMicrographAutomaticPicking2::run()
{
    Micrograph m;
    m.open_micrograph(fn_micrograph);
    AutoParticlePicking2 *autoPicking = new AutoParticlePicking2(fn_micrograph,&m,size,filter_num,NPCA,corr_num);
    // Resize the Micrograph
    selfScaleToSizeFourier((m.Ydim)*autoPicking->scaleRate,(m.Xdim)*autoPicking->scaleRate,autoPicking->microImage(), 2);

    FileName familyName=fn_model.removeDirectories();
    FileName fnAutoParticles=familyName+"@"+fn_root+"_auto.pos";
    FileName fnVectors=fn_root + "_auto_feature_vectors_"+familyName+".txt";
    FileName fnFilterBank=fn_micrograph.removeLastExtension()+"_filterbank.stk";
    FileName fnInvariant=fn_root + "_invariant";
    // Generating the filter bank
    filterBankGenerator(autoPicking->microImage(), fnFilterBank, filter_num);

    if (mode == "buildinv")
    {
        MetaData MD;
        // Insert all true positives
        if (fn_train!="")
        {
            MD.read(fn_train);
            int x, y;
            FOR_ALL_OBJECTS_IN_METADATA(MD)
            {
                MD.getValue(MDL_XINT, x, __iter.objId);
                MD.getValue(MDL_YINT, y, __iter.objId);
                m.add_coord(x, y, 0, 1);
            }
        }
        autoPicking->extractInvariant(fnFilterBank,fnInvariant);
    }
    if (mode == "train")
    {
        autoPicking->trainPCA(fn_root);
        //autoPicking->trainSVM(fn_root);
    }
}


