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

//To calculate the euclidean distance between to points
double euclidean_distance(const Particle2 &p1, const Particle2 &p2)
{
    double dx = (p1.x - p2.x);
    double dy = (p1.y - p2.y);
    return sqrt(dx * dx + dy * dy);
}

bool AutoParticlePicking2::checkDist(Particle2 &p)
{
    int num_part = __m->ParticleNo();
    int dist=0,min;
    Particle2 posSample;

    posSample.x=(__m->coord(0).X)*scaleRate;
    posSample.y=(__m->coord(0).Y)*scaleRate;
    min=euclidean_distance(p,posSample);
    for (int i=1;i<num_part;i++)
    {
        posSample.x=(__m->coord(i).X)*scaleRate;
        posSample.y=(__m->coord(i).Y)*scaleRate;
        dist=euclidean_distance(p,posSample);
        if (dist<min)
            min=dist;
    }

    if (min>(0.25*particle_radius))
        return true;
    else
        return false;
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
    int startX,startY,endX,endY;

    for (int j=0;j<filter_num;j++)
    {
        startX=x-particle_radius;
        startY=y-particle_radius;
        endX=x+particle_radius;
        endY=y+particle_radius;

        micrographStack.readMapped(fnFilterBank,j+1);
        extractParticle(x,y,micrographStack,pieceImage);
        mIpolar.aliasImageInStack(Ipolar,j);
        convert2Polar(pieceImage,mIpolar);
    }
    polarCorrelation(Ipolar,invariantChannel);
}
double AutoParticlePicking2::PCAProject(MultidimArray<double> &pcaBasis,MultidimArray<double> &vec,
                                        MultidimArray<double> &avg)
{

    double dotProduct=0;
    vec/=avg;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(pcaBasis)
    dotProduct += DIRECT_A1D_ELEM(pcaBasis,i) * DIRECT_A1D_ELEM(vec,i);
    return dotProduct;
}
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
    FileName fnModel= fn_root + "_training.txt";
    classifier.SVMTrain(trainSet,classLabel);
    classifier.SaveModel(fnModel);
}

void AutoParticlePicking2::add2Dataset(const FileName &fn_Invariant,const FileName &fn_root,int lable)
{
    FileName fnPCAModel=fn_root + "_pca_model.stk";
    ImageGeneric positiveInvariant;
    ImageGeneric pcaBasisStack;
    MultidimArray<double> vec;
    MultidimArray<double> avg;
    MultidimArray<double> pcaBase;
    ArrayDim aDim;
    int yDataSet=YSIZE(trainSet);

    positiveInvariant.read(fn_Invariant, HEADER);
    int num_correlation=filter_num+((filter_num-corr_num)*corr_num);
    positiveInvariant.getDimensions(aDim);
    int steps= aDim.ndim/num_correlation;
    trainSet.resize(1,1,yDataSet+steps,num_correlation*NPCA);
    classLabel.resize(1,1,1,YSIZE(trainSet));
    for (int n=yDataSet;n<XSIZE(classLabel);n++)
        classLabel(n)=lable;
    for (int i=0;i<num_correlation;i++)
    {
        //pcaAnalyzer.clear();
        pcaBasisStack.readMapped(fnPCAModel,i*(NPCA+1)+1);
        pcaBasisStack().getImage(avg);
        avg.resize(1,1,1,XSIZE(avg)*YSIZE(avg));
        for (int j=0;j<NPCA;j++)
        {
            pcaBasisStack.readMapped(fnPCAModel,(i*(NPCA+1)+1)+(j+1));
            pcaBasisStack().getImage(pcaBase);
            pcaBase.resize(1,1,1,XSIZE(pcaBase)*YSIZE(pcaBase));
            for (int k=0;k<steps;k++)
            {
                positiveInvariant.readMapped(fn_Invariant,k*num_correlation+i+1);
                positiveInvariant().getImage(vec);
                vec.resize(1,1,1,XSIZE(vec)*YSIZE(vec));
                DIRECT_A2D_ELEM(trainSet,k+yDataSet,j+(i*NPCA))=PCAProject(pcaBase,vec,avg);
            }
        }
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

void AutoParticlePicking2::extractNegativeInvariant(const FileName &fnFilterBank,const FileName &fnInvariantFeat)
{
    MultidimArray<double> IpolarCorr;
    std::vector<Particle2> negativeSamples;

    int num_part = __m->ParticleNo();
    int num_correlation=filter_num+((filter_num-corr_num)*corr_num);
    FileName fnNegativeInvariatn=fnInvariantFeat+"_Negative.stk";
    fnNegativeInvariatn.deleteFile();
    IpolarCorr.initZeros(num_correlation,1,NangSteps,NRsteps);
    extractNonParticle(negativeSamples);

    for (int i=0;i<negativeSamples.size();i++)
    {
        int x = negativeSamples[i].x;
        int y = negativeSamples[i].y;
        buildInvariant(IpolarCorr,x,y,fnFilterBank);
        Image<double> II;
        II() = IpolarCorr;
        II.write(fnNegativeInvariatn,ALL_IMAGES,true,WRITE_APPEND);
    }
}

void AutoParticlePicking2::extractInvariant(const FileName &fnFilterBank,const FileName &fnInvariantFeat)
{
    extractPositiveInvariant(fnFilterBank,fnInvariantFeat);
    extractNegativeInvariant(fnFilterBank,fnInvariantFeat);
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

void AutoParticlePicking2::extractNonParticle(std::vector<Particle2> &negativePosition)
{
    int endX, endY;
    int gridStep=particle_radius/2;
    Particle2 negSample;

    endX=XSIZE(microImage())-particle_radius;
    endY=YSIZE(microImage())-particle_radius;

    for (int i=particle_radius;i<endY;i=i+gridStep)
        for (int j=particle_radius;j<endX;j=j+gridStep)
        {
            negSample.y=i;
            negSample.x=j;
            if (checkDist(negSample))
            {
                negativePosition.push_back(negSample);
            }
        }
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
    FileName fnFilterBank=fn_micrograph.removeLastExtension()+"_filterbank.stk";
    FileName fnInvariant=fn_root + "_invariant";
    AutoParticlePicking2 *autoPicking = new AutoParticlePicking2(fn_micrograph,&m,size,filter_num,NPCA,corr_num);

    // Resize the Micrograph
    selfScaleToSizeFourier((m.Ydim)*autoPicking->scaleRate,(m.Xdim)*autoPicking->scaleRate,autoPicking->microImage(), 2);
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
        // Train PCA using positive smaples
        autoPicking->trainPCA(fn_root);
        // Add positive samples to the dataset
        autoPicking->add2Dataset(fnInvariant+"_Positive.stk",fn_root,1);
        // Add negative samples to the dataset
        autoPicking->add2Dataset(fnInvariant+"_Negative.stk",fn_root,2);
        // Train the SVM with positive and negative samples
        autoPicking->trainSVM(fn_root);
        //autoPicking->trainSVM(fn_root);
    }
}


