/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2008)
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

#include "tomo_align_tilt_series.h"

#include <data/args.h>
#include <data/filters.h>
#include <data/xmipp_fftw.h>
#include <data/metadata.h>
#include <data/numerical_tools.h>
#include <data/morphology.h>
#include <fstream>
#include <queue>
#include <iostream>

#include "fourier_filter.h"

/* Generate mask ----------------------------------------------------------- */
//#define DEBUG
void generateMask(const MultidimArray<double> &I, MultidimArray<unsigned char> &mask,
                  int patchSize)
{
    MultidimArray<double> dmask;
    dmask.initZeros(I);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(dmask)
    if (I(i,j)!=0)
        dmask(i,j)=1;

    MultidimArray<double> maskEroded;
    maskEroded.initZeros(dmask);
    erode2D(dmask,maskEroded,8,0,patchSize);
    typeCast(maskEroded,mask);
    mask.setXmippOrigin();

#ifdef DEBUG

    ImageXmipp save;
    save()=I;
    save.write("PPP.xmp");
    save()=maskEroded;
    save.write("PPPmask.xmp");
    std::cout << "Masks saved. Press any key\n";
    char c;
    std::cin >> c;
#endif
}
#undef DEBUG

/* Compute affine matrix --------------------------------------------------- */
class AffineFitness
{
public:
    Matrix1D<double> minAllowed;
    Matrix1D<double> maxAllowed;
    std::vector< MultidimArray<double> *> I1;
    std::vector< MultidimArray<double> *> I2;
    std::vector< MultidimArray<double> *> Mask1;
    std::vector< MultidimArray<double> *> Mask2;
    bool showMode;
    bool checkRotation;

    AffineFitness()
    {
        showMode=false;
        checkRotation=true;
    }

    ~AffineFitness()
    {
        for (size_t i=0; i<I1.size();    i++)
            delete I1[i];
        for (size_t i=0; i<I2.size();    i++)
            delete I2[i];
        for (size_t i=0; i<Mask1.size(); i++)
            delete Mask1[i];
        for (size_t i=0; i<Mask2.size(); i++)
            delete Mask2[i];
    }

    double affine_fitness_individual(double *p)
    {
        // Check limits
        if (!showMode)
            FOR_ALL_ELEMENTS_IN_MATRIX1D(minAllowed)
            if (p[i]<minAllowed(i) || p[i]>maxAllowed(i))
                return 1e20;

        // Separate solution
        Matrix2D<double> A12, A21;
        A12.initIdentity(3);
        A12(0,0)=p[0];
        A12(0,1)=p[1];
        A12(0,2)=p[4];
        A12(1,0)=p[2];
        A12(1,1)=p[3];
        A12(1,2)=p[5];

        A21=A12.inv();

        // Check it is approximately a rotation
        if (!showMode && checkRotation)
        {
            // Check that the eigenvalues are close to 1
            std::complex<double> a=A12(0,0);
            std::complex<double> b=A12(0,1);
            std::complex<double> c=A12(1,0);
            std::complex<double> d=A12(1,1);
            std::complex<double> eig1=(a+d)/2.0+sqrt(4.0*b*c+(a-d)*(a-d))/2.0;
            std::complex<double> eig2=(a+d)/2.0-sqrt(4.0*b*c+(a-d)*(a-d))/2.0;
            double abs_eig1=abs(eig1);
            double abs_eig2=abs(eig2);

            if (abs_eig1<0.95 || abs_eig1>1.05)
                return 1e20*ABS(abs_eig1-1);
            if (abs_eig2<0.95 || abs_eig2>1.05)
                return 1e20*ABS(abs_eig2-1);

            // Compute the eigenvalues of A*A^t disregarding the shifts
            // A*A^t=[a b; c d]

            a=A12(0,0)*A12(0,0)+A12(0,1)*A12(0,1);
            b=A12(0,0)*A12(1,0)+A12(0,1)*A12(1,1);
            c=b;
            d=A12(1,0)*A12(1,0)+A12(1,1)*A12(1,1);
            eig1=(a+d)/2.0+sqrt(4.0*b*c+(a-d)*(a-d))/2.0;
            eig2=(a+d)/2.0-sqrt(4.0*b*c+(a-d)*(a-d))/2.0;
            abs_eig1=abs(eig1);
            abs_eig2=abs(eig2);

            if (abs_eig1<0.95 || abs_eig1>1.05)
                return 1e20*ABS(abs_eig1-1);
            if (abs_eig2<0.95 || abs_eig2>1.05)
                return 1e20*ABS(abs_eig2-1);
            if (abs(c)>0.05)
                return 1e20*ABS(abs(c));
        }

        // For each pyramid level
        double dist=0;
        Matrix2D<double> A12level=A12;
        Matrix2D<double> A21level=A21;
        for (size_t level=0; level<I1.size(); level++)
        {
            const MultidimArray<double> &Mask1_level=*(Mask1[level]);
            const MultidimArray<double> &Mask2_level=*(Mask2[level]);

            // Produce the transformed images
            MultidimArray<double> transformedI1, transformedI2;
            applyGeometry(LINEAR,transformedI1,*(I1[level]),A12level,IS_NOT_INV,DONT_WRAP);
            applyGeometry(LINEAR,transformedI2,*(I2[level]),A21level,IS_NOT_INV,DONT_WRAP);

            // Produce masks for the comparison
            MultidimArray<int> maskInTheSpaceOf1, maskInTheSpaceOf2;
            MultidimArray<double> maskAux;
            applyGeometry(LINEAR,maskAux,Mask1_level,A12level,IS_NOT_INV,DONT_WRAP);
            maskInTheSpaceOf2.initZeros(YSIZE(maskAux),XSIZE(maskAux));
            maskInTheSpaceOf2.setXmippOrigin();
            FOR_ALL_ELEMENTS_IN_ARRAY2D(maskAux)
            {
                if (A2D_ELEM(maskAux,i,j)>0.5)
                    A2D_ELEM(maskInTheSpaceOf2,i,j)=(int)A2D_ELEM(Mask2_level,i,j);
            }

            maskAux.initZeros();
            applyGeometry(LINEAR,maskAux,Mask2_level,A21level,IS_NOT_INV,DONT_WRAP);
            maskInTheSpaceOf1.initZeros(YSIZE(maskAux),XSIZE(maskAux));
            maskInTheSpaceOf1.setXmippOrigin();
            FOR_ALL_ELEMENTS_IN_ARRAY2D(maskAux)
            {
                if (A2D_ELEM(maskAux,i,j)>0.5)
                    A2D_ELEM(maskInTheSpaceOf1,i,j)=(int)A2D_ELEM(Mask1_level,i,j);
            }

            // Compare the two images
            double distLevel=0.5*(
                                 1-correlationIndex(transformedI1,*(I2[level]),&maskInTheSpaceOf2)+
                                 1-correlationIndex(transformedI2,*(I1[level]),&maskInTheSpaceOf1));

            if (showMode)
            {
                Image<double> save;
                save()=*(I1[level]);
                save.write("PPPimg1.xmp");
                save()=*(I2[level]);
                save.write("PPPimg2.xmp");
                save()=*(Mask1[level]);
                save.write("PPPmask1.xmp");
                save()=*(Mask2[level]);
                save.write("PPPmask2.xmp");
                save()=transformedI1;
                save.write("PPPTransformedImg1.xmp");
                save()=transformedI2;
                save.write("PPPTransformedImg2.xmp");
                typeCast(maskInTheSpaceOf1,save());
                save.write("PPPmaskInTheSpaceOf1.xmp");
                typeCast(maskInTheSpaceOf2,save());
                save.write("PPPmaskInTheSpaceOf2.xmp");
                save()=*(I1[level])-transformedI2;
                save.write("PPPDiffImg1.xmp");
                save()=*(I2[level])-transformedI1;
                save.write("PPPDiffImg2.xmp");
                std::cout << "Level=" << level << "\nA12level=\n" << A12level << "A21level=\n" << A21level << std::endl;
                std::cout << "distLevel=" << distLevel << std::endl;
                std::cout << "Do you like the alignment? (y/n)\n";
                char c;
                std::cin >> c;
                if (c=='n')
                {
                    dist=-1;
                    break;
                }
            }
            dist+=distLevel;

            A12level(0,2)/=2;
            A12level(1,2)/=2;
            A21level(0,2)/=2;
            A21level(1,2)/=2;
        }
        dist/=I1.size();

        return dist;
    }

    static double Powell_affine_fitness_individual(double *p, void *prm)
    {
        return ((AffineFitness*)prm)->affine_fitness_individual(p+1);
    }
};

void setupAffineFitness(AffineFitness &fitness, const MultidimArray<double> &I1,
                        const MultidimArray<double> &I2, int maxShift, bool isMirror,
                        bool checkRotation, int pyramidLevel)
{
    fitness.checkRotation=checkRotation;

    // Set images
    int level=0;
    MultidimArray<double> I1aux=I1;
    MultidimArray<double> I2aux=I2;

    // Remove the borders
    int borderY=CEIL(0.1*YSIZE(I1));
    int borderX=CEIL(0.1*XSIZE(I1));
    I1aux.selfWindow(STARTINGY(I1)+borderY,STARTINGX(I1)+borderX,
                     FINISHINGY(I1)-borderY,FINISHINGX(I1)-borderX);
    I2aux.selfWindow(STARTINGY(I1)+borderY,STARTINGX(I1)+borderX,
                     FINISHINGY(I1)-borderY,FINISHINGX(I1)-borderX);

    MultidimArray<double> Mask1;
    Mask1.initZeros(I1aux);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I1aux)
    if (A2D_ELEM(I1,i,j)!=0)
        A2D_ELEM(Mask1,i,j)=1;

    MultidimArray<double> Mask2;
    Mask2.initZeros(I2aux);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I2aux)
    if (A2D_ELEM(I2,i,j)!=0)
        A2D_ELEM(Mask2,i,j)=1;

    do
    {
        // Push back the current images
        I1aux.setXmippOrigin();
        I2aux.setXmippOrigin();
        Mask1.setXmippOrigin();
        Mask2.setXmippOrigin();

        MultidimArray<double> *dummy=NULL;
        dummy=new MultidimArray<double>;
        *dummy=I1aux;
        fitness.I1.push_back(dummy);
        dummy=new MultidimArray<double>;
        *dummy=I2aux;
        fitness.I2.push_back(dummy);
        dummy=new MultidimArray<double>;
        *dummy=Mask1;
        fitness.Mask1.push_back(dummy);
        dummy=new MultidimArray<double>;
        *dummy=Mask2;
        fitness.Mask2.push_back(dummy);

        // Prepare for next level
        level++;
        if (pyramidLevel>=level)
        {
            selfScaleToSize(LINEAR,I1aux,YSIZE(I1aux)/2,XSIZE(I1aux)/2);
            selfScaleToSize(LINEAR,I2aux,YSIZE(I2aux)/2,XSIZE(I2aux)/2);
            selfScaleToSize(LINEAR,Mask1,YSIZE(Mask1)/2,XSIZE(Mask1)/2);
            selfScaleToSize(LINEAR,Mask2,YSIZE(Mask2)/2,XSIZE(Mask2)/2);
            FOR_ALL_ELEMENTS_IN_ARRAY2D(Mask1)
            Mask1(i,j)=(Mask1(i,j)>0.5)? 1:0;
            FOR_ALL_ELEMENTS_IN_ARRAY2D(Mask2)
            Mask2(i,j)=(Mask2(i,j)>0.5)? 1:0;
        }
    }
    while (level<=pyramidLevel);

    // Set limits for the affine matrices
    // Order: 1->2: 4 affine params+2 translations
    // Order: 2->1: 4 affine params+2 translations
    fitness.minAllowed.resize(6);
    fitness.maxAllowed.resize(6);

    // Scale factors
    fitness.minAllowed(0)=fitness.minAllowed(3)=0.9;
    fitness.maxAllowed(0)=fitness.maxAllowed(3)=1.1;
    if (isMirror)
    {
        fitness.minAllowed(3)=-1.1;
        fitness.maxAllowed(3)=-0.9;
    }

    // Rotation factors
    fitness.minAllowed(1)=fitness.minAllowed(2)=-0.2;
    fitness.maxAllowed(1)=fitness.maxAllowed(2)= 0.2;

    // Shifts
    fitness.minAllowed(4)=fitness.minAllowed(5)=-maxShift;
    fitness.maxAllowed(4)=fitness.maxAllowed(5)= maxShift;
}

static pthread_mutex_t globalAffineMutex = PTHREAD_MUTEX_INITIALIZER;
double computeAffineTransformation(const MultidimArray<unsigned char> &I1,
                                   const MultidimArray<unsigned char> &I2, int maxShift, int maxIterDE,
                                   Matrix2D<double> &A12, Matrix2D<double> &A21, bool show,
                                   double thresholdAffine, bool globalAffine, bool isMirror,
                                   bool checkRotation, int pyramidLevel)
{
    try
    {
        AffineFitness affy;
        MultidimArray<double> I1d, I2d;
        typeCast(I1, I1d);
        typeCast(I2, I2d);
        setupAffineFitness(affy, I1d, I2d, maxShift, isMirror, checkRotation,
                           pyramidLevel);

        // Return result
        double cost;

        // Optimize with differential evolution
        Matrix1D<double> A(6);
        A(0)=A(3)=1;
        if (isMirror)
            A(3)*=-1;
        if (MAT_XSIZE(A12)==0)
        {
            if (globalAffine)
            {
                // Exhaustive search
                double bestCost=2;
                double stepX=(affy.maxAllowed(4)-affy.minAllowed(4))/40.0;
                double stepY=(affy.maxAllowed(5)-affy.minAllowed(5))/40.0;
                for (double shiftY=affy.minAllowed(5); shiftY<=affy.maxAllowed(5); shiftY+=stepY)
                    for (double shiftX=affy.minAllowed(4); shiftX<=affy.maxAllowed(4); shiftX+=stepX)
                    {
                        A(4)=shiftX;
                        A(5)=shiftY;
                        double cost=affy.affine_fitness_individual(MATRIX1D_ARRAY(A));
                        if (cost<bestCost)
                            bestCost=cost;
                    }
            }
            else
            {
                // Initialize with cross correlation
                double tx, ty;
                pthread_mutex_lock( &globalAffineMutex );
                CorrelationAux aux;
                if (!isMirror)
                    bestShift(I1d,I2d,tx,ty,aux);
                else
                {
                    MultidimArray<double> auxI2d=I2d;
                    auxI2d.selfReverseY();
                    STARTINGX(auxI2d)=STARTINGX(I2d);
                    STARTINGY(auxI2d)=STARTINGY(I2d);
                    bestShift(I1d,auxI2d,tx,ty,aux);
                    ty=-ty;
                }
                pthread_mutex_unlock( &globalAffineMutex );
                A(4)=-tx;
                A(5)=-ty;
            }

            // Optimize with Powell
            Matrix1D<double> steps(A);
            steps.initConstant(1);
            int iter;
            powellOptimizer(A, 1, VEC_XSIZE(A),
                            AffineFitness::Powell_affine_fitness_individual, &affy, 0.005,
                            cost, iter, steps, false);

            // Separate solution
            A12.initIdentity(3);
            A12(0,0)=A(0);
            A12(0,1)=A(1);
            A12(0,2)=A(4);
            A12(1,0)=A(2);
            A12(1,1)=A(3);
            A12(1,2)=A(5);

            A21=A12.inv();
        }

        if (show)
        {
            affy.showMode=true;
            Matrix1D<double> p(6);
            p(0)=A12(0,0);
            p(1)=A12(0,1);
            p(4)=A12(0,2);
            p(2)=A12(1,0);
            p(3)=A12(1,1);
            p(5)=A12(1,2);
            cost = affy.affine_fitness_individual(MATRIX1D_ARRAY(p));
            affy.showMode=false;
        }

        return cost;
    }
    catch (XmippError &XE)
    {
        std::cout << XE;
        exit(1);
    }
}

/* Parameters -------------------------------------------------------------- */
void ProgTomographAlignment::readParams()
{
    fnSel=getParam("-i");
    fnSelOrig=getParam("--iorig");
    fnRoot=getParam("--oroot");
    if (fnRoot=="")
        fnRoot=fnSel.withoutExtension();
    globalAffine=checkParam("--globalAffine");
    useCriticalPoints=checkParam("--useCriticalPoints");
    if (useCriticalPoints)
        Ncritical=getIntParam("--useCriticalPoints");
    else
        Ncritical=0;
    seqLength=getIntParam("--seqLength");
    blindSeqLength=getIntParam("--blindSeqLength");
    maxStep=getIntParam("--maxStep");
    gridSamples=getIntParam("--gridSamples");
    psiMax=getDoubleParam("--psiMax");
    deltaRot=getDoubleParam("--deltaRot");
    localSize=getDoubleParam("--localSize");
    optimizeTiltAngle=checkParam("--optimizeTiltAngle");
    isCapillar=checkParam("--isCapillar");
    dontNormalize=checkParam("--dontNormalize");
    difficult=checkParam("--difficult");
    corrThreshold=getDoubleParam("--threshold");
    maxShiftPercentage=getDoubleParam("--maxShiftPercentage");
    maxIterDE=getIntParam("--maxIterDE");
    showAffine=checkParam("--showAffine");
    thresholdAffine=getDoubleParam("--thresholdAffine");
    identifyOutliersZ = getDoubleParam("--identifyOutliers");
    doNotIdentifyOutliers = checkParam("--noOutliers");
    pyramidLevel = getIntParam("--pyramid");
    numThreads = getIntParam("--thr");
    if (numThreads<1)
        numThreads = 1;
    lastStep=getIntParam("--lastStep");
}

void ProgTomographAlignment::show()
{
    std::cout << "Input images:       " << fnSel              << std::endl
    << "Original images:    " << fnSelOrig          << std::endl
    << "Output rootname:    " << fnRoot             << std::endl
    << "Global affine:      " << globalAffine       << std::endl
    << "Use critical points:" << useCriticalPoints  << std::endl
    << "Num critical points:" << Ncritical          << std::endl
    << "SeqLength:          " << seqLength          << std::endl
    << "BlindSeqLength:     " << blindSeqLength     << std::endl
    << "MaxStep:            " << maxStep            << std::endl
    << "Grid samples:       " << gridSamples        << std::endl
    << "Maximum psi:        " << psiMax             << std::endl
    << "Delta rot:          " << deltaRot           << std::endl
    << "Local size:         " << localSize          << std::endl
    << "Optimize tilt angle:" << optimizeTiltAngle  << std::endl
    << "isCapillar:         " << isCapillar         << std::endl
    << "DontNormalize:      " << dontNormalize      << std::endl
    << "Difficult:          " << difficult          << std::endl
    << "Threshold:          " << corrThreshold      << std::endl
    << "MaxShift Percentage:" << maxShiftPercentage << std::endl
    << "MaxIterDE:          " << maxIterDE          << std::endl
    << "Show Affine:        " << showAffine         << std::endl
    << "Threshold Affine:   " << thresholdAffine    << std::endl
    << "Identify outliers Z:" << identifyOutliersZ     << std::endl
    << "No outliers:        " << doNotIdentifyOutliers << std::endl
    << "Pyramid level:      " << pyramidLevel       << std::endl
    << "Threads to use:     " << numThreads         << std::endl
    << "Last step:          " << lastStep           << std::endl
    ;
}

void ProgTomographAlignment::defineParams()
{
    addUsageLine("Align a single-axis tilt series without any marker.");
    addSeeAlsoLine("tomo_align_dual_tilt_series");
    addParamsLine(" == General Options == ");
    addParamsLine("   -i <metadatafile>              : Input images");
    addParamsLine("                                  : The selfile must contain the list of micrographs");
    addParamsLine("                                  : and its tilt angles");
    addParamsLine("  [--iorig <metadatafile=\"\">]   : Metadata with images at original scale");
    addParamsLine("  [--oroot <fn_out=\"\">]         : Output alignment");
    addParamsLine("                                  : If not given, the input selfile without extension");
    addParamsLine("  [--thr <num=1>]                 : Parallel processing using \"num\" threads");
    addParamsLine("  [--lastStep+ <step=-1>]         : Last step to perform");
    addParamsLine("                                  : Step -1 -> Perform all steps");
    addParamsLine("                                  : Step  0 -> Determination of affine transformations");
    addParamsLine("                                  : Step  1 -> Determination of landmark chains");
    addParamsLine("                                  : Step  2 -> Determination of alignment parameters");
    addParamsLine("                                  : Step  3 -> Writing aligned images");
    addParamsLine(" == Step 0 (Affine alignment) Options == ");
    addParamsLine("  [--maxShiftPercentage <p=0.2>]   : Maximum shift as percentage of image size");
    addParamsLine("  [--thresholdAffine <th=0.85>]    : Threshold affine");
    addParamsLine("  [--globalAffine]                 : Look globally for affine transformations");
    addParamsLine("  [--difficult+]                   : Apply some filters before affine alignment");
    addParamsLine("  [--maxIterDE+ <n=30>]            : Maximum number of iteration in Differential Evolution");
    addParamsLine("  [--showAffine+]                  : Show affine transformations as PPP*");
    addParamsLine("  [--identifyOutliers+ <z=5>]      : Z-score to be an outlier");
    addParamsLine("  [--noOutliers+]                  : Do not identify outliers");
    addParamsLine("  [--pyramid+ <level=1>]           : Multiresolution for affine transformations");
    addParamsLine(" == Step 1 (Landmark chain) Options == ");
    addParamsLine("  [--seqLength <n=5>]              : Sequence length");
    addParamsLine("  [--localSize <size=0.04>]        : In percentage");
    addParamsLine("  [--useCriticalPoints <n=0>]      : Use critical points instead of a grid");
    addParamsLine("                                   : n is the number of critical points to choose");
    addParamsLine("                                   : in each image");
    addParamsLine("  [--threshold+ <th=-1>]           : Correlation threshold");
    addParamsLine("  [--blindSeqLength+ <n=-1>]       : Blind sequence length, -1=No blind landmarks");
    addParamsLine("  [--maxStep+ <step=4>]            : Maximum step for chain refinement");
    addParamsLine("  [--gridSamples+ <n=40>]           : Total number of samples=n*n");
    addParamsLine("  [--isCapillar+]                  : Set this flag if the tilt series is of a capillar");
    addParamsLine(" == Step 2 (Determination of alignment parameters) Options == ");
    addParamsLine("  [--psiMax+ <psi=-1>]             : Maximum psi in absolute value (degrees)");
    addParamsLine("                                  : -1 -> do not optimize for psi");
    addParamsLine("  [--deltaRot+ <rot=5>]            : In degrees. For the first optimization stage");
    addParamsLine("  [--optimizeTiltAngle+]           : Optimize tilt angle");
    addParamsLine(" == Step 3 (Produce aligned images) Options == ");
    addParamsLine("  [--dontNormalize+]               : Don't normalize the output images");
    addExampleLine("Typical run",false);
    addExampleLine("xmipp_tomo_align_tilt_series -i tiltseries.sel --thr 8");
    addExampleLine("If there are image with large shifts",false);
    addExampleLine("xmipp_tomo_align_tilt_series -i tiltseries.sel --thr 8 --globalAffine");
    addExampleLine("If there are clear landmarks that can be tracked",false);
    addExampleLine("xmipp_tomo_align_tilt_series -i tiltseries.sel --thr 8 --criticalPoints");
}

/* Produce side info ------------------------------------------------------- */
static pthread_mutex_t printingMutex = PTHREAD_MUTEX_INITIALIZER;
struct ThreadComputeTransformParams
{
    int myThreadID;
    ProgTomographAlignment * parent;
};

void * threadComputeTransform( void * args )
{
    ThreadComputeTransformParams * master =
        (ThreadComputeTransformParams *) args;

    ProgTomographAlignment * parent = master->parent;
    int thread_id = master->myThreadID;
    int localnumThreads = parent->numThreads;
    bool isCapillar = parent->isCapillar;
    int Nimg = parent->Nimg;
    double maxShiftPercentage = parent->maxShiftPercentage;
    int maxIterDE = parent->maxIterDE;
    bool showAffine = parent->showAffine;
    double thresholdAffine = parent->thresholdAffine;
    std::vector < MultidimArray<unsigned char> *> & img = parent->img;
    std::vector< std::vector< Matrix2D<double> > > & affineTransformations = parent->affineTransformations;
    double globalAffine = parent->globalAffine;

    int maxShift=FLOOR(XSIZE(*img[0])*maxShiftPercentage);
    int initjj=1;
    if (isCapillar)
        initjj=0;

    initjj += thread_id;

    double cost;
    for (int jj=initjj; jj<Nimg; jj+= localnumThreads)
    {
        int jj_1;
        if (isCapillar)
            jj_1=intWRAP(jj-1,0,Nimg-1);
        else
            jj_1=jj-1;
        MultidimArray<unsigned char>& img_i=*img[jj_1];
        MultidimArray<unsigned char>& img_j=*img[jj];
        bool isMirror=(jj==0) && (jj_1==Nimg-1);

        if (MAT_XSIZE(affineTransformations[jj_1][jj])==0)
        {
            Matrix2D<double> Aij, Aji;
            cost = computeAffineTransformation(img_i, img_j, maxShift,
                                               maxIterDE, Aij, Aji,
                                               showAffine, thresholdAffine, globalAffine,
                                               isMirror,true, parent->pyramidLevel);
            parent->correlationList[jj]=1-cost;

            pthread_mutex_lock( &printingMutex );
            affineTransformations[jj_1][jj]=Aij;
            affineTransformations[jj][jj_1]=Aji;
            if (cost<1)
                parent->writeTransformations(
                    parent->fnRoot+"_transformations.txt");
            std::cout << "Cost for [" << jj_1 << "] - ["
            << jj << "] = " << cost << std::endl;
            parent->iteration++;
            pthread_mutex_unlock( &printingMutex );
        }
        else
        {
            Matrix2D<double> Aij;
            Aij=affineTransformations[jj_1][jj];

            MultidimArray<double> img_id, img_jd;
            typeCast(img_i, img_id);
            typeCast(img_j, img_jd);

            AffineFitness fitness;
            setupAffineFitness(fitness, img_id, img_jd, maxShift, isMirror,
                               false, parent->pyramidLevel);
            Matrix1D<double> p(6);
            p(0)=Aij(0,0);
            p(1)=Aij(0,1);
            p(4)=Aij(0,2);
            p(2)=Aij(1,0);
            p(3)=Aij(1,1);
            p(5)=Aij(1,2);
            cost = fitness.affine_fitness_individual(MATRIX1D_ARRAY(p));
            parent->correlationList[jj]=1-cost;

            pthread_mutex_lock( &printingMutex );
            std::cout << "Cost for [" << jj_1 << "] - ["
            << jj << "] = " << cost << std::endl;
            pthread_mutex_unlock( &printingMutex );
            parent->iteration++;
        }
    }

    return NULL;
}

void ProgTomographAlignment::computeAffineTransformations(
    bool globalAffineToUse)
{
    bool oldglobalAffine=globalAffine;
    globalAffine=globalAffineToUse;

    pthread_t * th_ids = new pthread_t[numThreads];
    ThreadComputeTransformParams * th_args = new ThreadComputeTransformParams[numThreads];

    for( int nt = 0 ; nt < numThreads ; nt ++ )
    {
        // Passing parameters to each thread
        th_args[nt].parent = this;
        th_args[nt].myThreadID = nt;
        pthread_create( (th_ids+nt) , NULL, threadComputeTransform, (void *)(th_args+nt) );
    }

    // Waiting for threads to finish
    for( int nt = 0 ; nt < numThreads ; nt ++ )
        pthread_join(*(th_ids+nt), NULL);

    // Threads structures are not needed any more
    delete[] th_ids;
    delete[] th_args;
    globalAffine=oldglobalAffine;
}

void ProgTomographAlignment::identifyOutliers(bool mark)
{
    isOutlier.initZeros(Nimg);
    MultidimArray<double> correlationListAux(Nimg);
    for (int i=0; i<Nimg; i++)
        correlationListAux(i)=correlationList[i];

    MultidimArray<double> diff;
    double medianCorr=correlationListAux.computeMedian();
    diff=correlationListAux-medianCorr;
    diff.selfABS();
    double madCorr=diff.computeMedian();

    std::cout << "Cost distribution= " << 1-medianCorr << " +- "
    << madCorr << std::endl;

    double thresholdCorr=medianCorr-identifyOutliersZ*madCorr;
    for (size_t i=1; i<VEC_XSIZE(isOutlier); i++)
    {
        bool potentialOutlier=(correlationListAux(i)<thresholdCorr);
        if (potentialOutlier)
        {
            affineTransformations[i-1][i].clear();
            if (mark)
            {
                isOutlier(i)=true;
                std::cout << name_list[i-1] << " [" << i
                << "] is considered as an outlier. "
                << "Its cost is " << 1-correlationListAux(i)
                << std::endl;
            }
            else
                std::cout << name_list[i-1] << " [" << i
                << "] might be an outlier. "
                << "Its cost is " << 1-correlationListAux(i)
                << std::endl;
        }
        iteration++;
    }
}

void ProgTomographAlignment::produceSideInfo()
{
    // Difficult images?
    if (difficult)
    {
        if (pyramidLevel==0)
            pyramidLevel=2;
        if (localSize<0.05)
            localSize=0.08;
        if (seqLength==5)
            seqLength=11;
    }

    bestPreviousAlignment=new Alignment(this);
    // Read input data
    SF.read(fnSel,NULL);
    if (SF.containsLabel(MDL_ENABLED))
        SF.removeObjects(MDValueEQ(MDL_ENABLED, -1));
    Nimg=SF.size();
    if (Nimg!=0)
    {
        // Clear the list of images if not empty
        if (!img.empty())
        {
            for (size_t i=0; i<img.size(); i++)
                delete img[i];
            img.clear();
        }

        // Clear the list of masks if not empty
        if (!maskImg.empty())
        {
            for (size_t i=0; i<maskImg.size(); i++)
                delete maskImg[i];
            maskImg.clear();
        }

        std::cerr << "Reading input data\n";
        init_progress_bar(Nimg);
        iteration=0;
        int n=0;

        iMinTilt=-1;
        double minTilt=1000;
        bool nonZeroTilt=false;
        Image<double> imgaux;
        FileName fn;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            SF.getValue( MDL_IMAGE, fn ,__iter.objId);
            imgaux.read(fn);
            if (difficult)
            {
                //ImageXmipp save;
                //save()=imgaux(); save.write("PPPoriginal.xmp");
                MultidimArray<double> Ifiltered;
                Ifiltered=imgaux();
                Ifiltered.setXmippOrigin();

                // Reject outliers
                reject_outliers(Ifiltered,0.5);

                // Substract background plane
                substractBackgroundPlane(Ifiltered);

                // Substract the background (rolling ball)
                substractBackgroundRollingBall(Ifiltered,
                                               XSIZE(Ifiltered)/10);

                // Bandpass the image
                FourierFilter FilterBP;
                FilterBP.FilterBand=BANDPASS;
                FilterBP.w1=1.0/XSIZE(Ifiltered);
                FilterBP.w2=100.0/XSIZE(Ifiltered);
                FilterBP.raised_w=1.0/XSIZE(Ifiltered);
                FilterBP.generateMask(Ifiltered);
                FilterBP.applyMaskSpace(Ifiltered);

                // Equalize histogram
                //histogram_equalization(Ifiltered,8);
                imgaux()=Ifiltered;
                //save()=Ifiltered; save.write("PPPpreprocessed.xmp");
                //std::cout << "Press\n";
                //char c; std::cin >> c;
            }

            if (!useCriticalPoints)
            {
                MultidimArray<unsigned char>* mask_i=new MultidimArray<unsigned char>;
                generateMask(imgaux(),*mask_i,
                             XMIPP_MAX(ROUND(localSize*XSIZE(imgaux()))/2,5));
                maskImg.push_back(mask_i);
            }

            MultidimArray<unsigned char>* img_i=new MultidimArray<unsigned char>;
            imgaux().rangeAdjust(0,255);
            typeCast(imgaux(),*img_i);
            img_i->setXmippOrigin();
            img.push_back(img_i);

            double tiltAngle;
            if (SF.containsLabel(MDL_ANGLE_TILT))
                SF.getValue(MDL_ANGLE_TILT,tiltAngle,__iter.objId);
            else
                tiltAngle=imgaux.tilt();
            tiltList.push_back(tiltAngle);
            if (tiltAngle!=0)
                nonZeroTilt=true;
            if (ABS(tiltAngle)<minTilt)
            {
                iMinTilt=n;
                minTilt=ABS(tiltAngle);
            }
            name_list.push_back(imgaux.name());

            progress_bar(n++);
            iteration++;
        }
        progress_bar(Nimg);
        if (!nonZeroTilt)
            REPORT_ERROR(ERR_VALUE_NOTSET,"Tilt angles have not been assigned to the input selfile");
    }

    // Read images at original scale
    if (!fnSelOrig.empty())
    {
        SForig.read(fnSelOrig,NULL);
        SForig.removeDisabled();

        if (SForig.size()!=SF.size())
            REPORT_ERROR(ERR_MD_OBJECTNUMBER,"The number of images in both selfiles (-i and -iorig) is different");
    }

    // Fill the affine transformations with empty matrices
    std::vector< Matrix2D<double> > emptyRow;
    for (int i=0; i<Nimg; i++)
    {
        Matrix2D<double> A;
        emptyRow.push_back(A);
    }
    for (int i=0; i<Nimg; i++)
    {
        affineTransformations.push_back(emptyRow);
        correlationList.push_back(-1);
    }

    // Check if there are transformations already calculated
    FileName fn_tmp = fnRoot+"_transformations.txt";
    if (fn_tmp.exists())
    {
        std::ifstream fhIn;
        fhIn.open(fn_tmp.c_str());
        if (!fhIn)
            REPORT_ERROR(ERR_IO_NOTEXIST,(String)"Cannot open "+ fn_tmp);
        size_t linesRead=0;
        while (!fhIn.eof())
        {
            std::string line;
            getline(fhIn, line);
            if (linesRead<affineTransformations.size()-1 && line!="")
            {
                std::vector< double > data;
                readFloatList(line.c_str(), 8, data);
                double tilt=data[0];
                int i=0;
                double bestDistance=ABS(tilt-tiltList[0]);
                for (int j=0; j<Nimg; j++)
                {
                    double distance=ABS(tilt-tiltList[j]);
                    if (bestDistance>distance)
                    {
                        bestDistance=distance;
                        i=j;
                    }
                }
                affineTransformations[i][i+1].resize(3,3);
                affineTransformations[i][i+1].initIdentity(3);
                affineTransformations[i][i+1](0,0)=data[2];
                affineTransformations[i][i+1](1,0)=data[3];
                affineTransformations[i][i+1](0,1)=data[4];
                affineTransformations[i][i+1](1,1)=data[5];
                affineTransformations[i][i+1](0,2)=data[6];
                affineTransformations[i][i+1](1,2)=data[7];

                affineTransformations[i+1][i]=
                    affineTransformations[i][i+1].inv();

                if (showAffine)
                {
                    MultidimArray<unsigned char>& img_i=*img[i];
                    MultidimArray<unsigned char>& img_j=*img[i+1];
                    int maxShift=FLOOR(XSIZE(*img[0])*maxShiftPercentage);
                    Matrix2D<double> Aij, Aji;
                    Aij=affineTransformations[i][i+1];
                    Aji=affineTransformations[i+1][i];
                    bool isMirror=false;
                    std::cout << "Transformation between " << i << " ("
                    << tiltList[i] << ") and " << i+1 << " ("
                    << tiltList[i+1] << std::endl;
                    double cost = computeAffineTransformation(img_i, img_j,
                                  maxShift, maxIterDE, Aij, Aji,
                                  showAffine, thresholdAffine, globalAffine,
                                  isMirror,true, pyramidLevel);
                    if (cost<0)
                    {
                        affineTransformations[i][i+1].clear();
                        affineTransformations[i+1][i].clear();
                        writeTransformations(fn_tmp);
                    }
                }
            }
            linesRead++;
        }
        fhIn.close();
    }

    // Compute the rest of transformations
    computeAffineTransformations(globalAffine);

    // Do not show refinement
    showRefinement = false;

    // Check which is the distribution of correlation
    if (!useCriticalPoints)
    {
        int X0=(int)(STARTINGX(*(img[0]))+2*localSize*XSIZE(*(img[0])));
        int XF=(int)(FINISHINGX(*(img[0]))-2*localSize*XSIZE(*(img[0])));
        int Y0=(int)(STARTINGY(*(img[0]))+2*localSize*YSIZE(*(img[0])));
        int YF=(int)(FINISHINGY(*(img[0]))-2*localSize*YSIZE(*(img[0])));
        avgForwardPatchCorr.initZeros(Nimg);
        avgBackwardPatchCorr.initZeros(Nimg);
        avgForwardPatchCorr.initConstant(1);
        avgBackwardPatchCorr.initConstant(1);

        for (int ii=0; ii<Nimg; ++ii)
        {
            // Compute average forward patch corr
            MultidimArray<double> corrList;
            if (ii<=Nimg-2)
            {
                corrList.initZeros(200);
                FOR_ALL_ELEMENTS_IN_ARRAY1D(corrList)
                {
                    Matrix1D<double> rii(3), rjj(3);
                    do
                    {
                        XX(rii)=rnd_unif(X0,XF);
                        YY(rii)=rnd_unif(Y0,YF);
                        rjj=affineTransformations[ii][ii+1]*rii;
                        refineLandmark(ii,ii+1,rii,rjj,corrList(i),false);
                    }
                    while (corrList(i)<-0.99);
                }
                avgForwardPatchCorr(ii)=corrList.computeAvg();
            }

            // Compute average backward patch corr
            if (ii>0)
            {
                corrList.initZeros(200);
                FOR_ALL_ELEMENTS_IN_ARRAY1D(corrList)
                {
                    Matrix1D<double> rii(3), rjj(3);
                    do
                    {
                        XX(rii)=rnd_unif(X0,XF);
                        YY(rii)=rnd_unif(Y0,YF);
                        rjj=affineTransformations[ii][ii-1]*rii;
                        refineLandmark(ii,ii-1,rii,rjj,corrList(i),false);
                    }
                    while (corrList(i)<-0.99);
                }
                avgBackwardPatchCorr(ii)=corrList.computeAvg();
            }

            std::cout << "Image " << ii << " Average correlation forward="
            << avgForwardPatchCorr(ii)
            << " backward=" << avgBackwardPatchCorr(ii) << std::endl;
            iteration++;
        }
    }

    // Identify outliers
    if (!doNotIdentifyOutliers)
    {
        identifyOutliers(false);
        computeAffineTransformations(false);
        identifyOutliers(true);
    }
    else
        isOutlier.initZeros(Nimg);
    if (lastStep==0)
        exit(0);
}

/* Generate landmark set --------------------------------------------------- */
//#define DEBUG
struct ThreadGenerateLandmarkSetParams
{
    int myThreadID;
    ProgTomographAlignment * parent;

    std::vector<LandmarkChain> *chainList;
};

void * threadgenerateLandmarkSetGrid( void * args )
{
    ThreadGenerateLandmarkSetParams * master =
        (ThreadGenerateLandmarkSetParams *) args;
    ProgTomographAlignment * parent = master->parent;
    int thread_id = master->myThreadID;
    int Nimg=parent->Nimg;
    int numThreads=parent->numThreads;
    const std::vector< std::vector< Matrix2D<double> > > &affineTransformations=
        parent->affineTransformations;
    int gridSamples=parent->gridSamples;

    int deltaShift=(int)floor(XSIZE(*(parent->img)[0])/gridSamples);
    master->chainList=new std::vector<LandmarkChain>;
    Matrix1D<double> rii(3), rjj(3);
    ZZ(rii)=1;
    ZZ(rjj)=1;
    if (thread_id==0)
        init_progress_bar(gridSamples);
    int includedPoints=0;
    Matrix1D<int> visited(Nimg);
    for (int nx=thread_id; nx<gridSamples; nx+=numThreads)
    {
        XX(rii)=STARTINGX(*(parent->img)[0])+ROUND(deltaShift*(0.5+nx));
        for (int ny=0; ny<gridSamples; ny+=1)
        {
            YY(rii)=STARTINGY(*(parent->img)[0])+ROUND(deltaShift*(0.5+ny));
            for (int ii=0; ii<=Nimg-1; ++ii)
            {
                // Check if inside the mask
                if (!(*(parent->maskImg[ii]))((int)YY(rii),(int)XX(rii)))
                    continue;
                if (parent->isOutlier(ii))
                    continue;

                LandmarkChain chain;
                chain.clear();
                Landmark l;
                l.x=XX(rii);
                l.y=YY(rii);
                l.imgIdx=ii;
                chain.push_back(l);
                visited.initZeros();
                visited(ii)=1;

                // Follow this landmark backwards
                bool acceptLandmark=true;
                int jj;
                if (parent->isCapillar)
                    jj=intWRAP(ii-1,0,Nimg-1);
                else
                    jj=ii-1;
                Matrix2D<double> Aij, Aji;
                Matrix1D<double> rcurrent=rii;
                while (jj>=0 && acceptLandmark && !visited(jj) &&
                       !parent->isOutlier(jj))
                {
                    // Compute the affine transformation between ii and jj
                    int jj_1;
                    if (parent->isCapillar)
                        jj_1=intWRAP(jj+1,0,Nimg-1);
                    else
                        jj_1=jj+1;
                    visited(jj)=1;
                    Aij=affineTransformations[jj][jj_1];
                    Aji=affineTransformations[jj_1][jj];
                    rjj=Aji*rcurrent;
                    double corr;
                    acceptLandmark=parent->refineLandmark(jj_1,jj,rcurrent,rjj,
                                                          corr,true);
                    if (acceptLandmark)
                    {
                        l.x=XX(rjj);
                        l.y=YY(rjj);
                        l.imgIdx=jj;
                        chain.push_back(l);
                        if (parent->isCapillar)
                            jj=intWRAP(jj-1,0,Nimg-1);
                        else
                            jj=jj-1;
                        rcurrent=rjj;
                    }
                }

                // Follow this landmark forward
                acceptLandmark=true;
                if (parent->isCapillar)
                    jj=intWRAP(ii+1,0,Nimg-1);
                else
                    jj=ii+1;
                rcurrent=rii;
                while (jj<Nimg && acceptLandmark && !visited(jj) &&
                       !parent->isOutlier(jj))
                {
                    // Compute the affine transformation between ii and jj
                    int jj_1;
                    if (parent->isCapillar)
                        jj_1=intWRAP(jj-1,0,Nimg-1);
                    else
                        jj_1=jj-1;
                    visited(jj)=1;
                    Aij=affineTransformations[jj_1][jj];
                    Aji=affineTransformations[jj][jj_1];
                    rjj=Aij*rcurrent;
                    double corr;
                    acceptLandmark=parent->refineLandmark(jj_1,jj,rcurrent,rjj,
                                                          corr,true);
                    if (acceptLandmark)
                    {
                        l.x=XX(rjj);
                        l.y=YY(rjj);
                        l.imgIdx=jj;
                        chain.push_back(l);
                        if (parent->isCapillar)
                            jj=intWRAP(jj+1,0,Nimg-1);
                        else
                            jj=jj+1;
                        rcurrent=rjj;
                    }
                }

#ifdef DEBUG
                std::cout << "img=" << ii << " chain length="
                << chain.size() << " [" << jjleft
                << " - " << jjright << "]";
#endif

                if (chain.size()>parent->seqLength)
                {
                    double corrChain;
                    bool accepted=parent->refineChain(chain,corrChain);
                    if (accepted)
                    {
#ifdef DEBUG
                        std::cout << " Accepted with length= "
                        << chain.size()
                        << " [" << chain[0].imgIdx << " - "
                        << chain[chain.size()-1].imgIdx << "]= ";
                        for (int i=0; i<chain.size(); i++)
                            std::cout << chain[i].imgIdx << " ";
#endif

                        master->chainList->push_back(chain);
                        includedPoints+=chain.size();
                    }
                }
#ifdef DEBUG
                std::cout << std::endl;
#endif

            }
#ifdef DEBUG
            std::cout << "Point nx=" << nx << " ny=" << ny
            << " Number of points="
            << includedPoints
            << " Number of chains=" << master->chainList->size()
            << " ( " << ((double) includedPoints)/
            master->chainList->size() << " )\n";
#endif

        }
        if (thread_id==0)
            progress_bar(nx);
    }
    if (thread_id==0)
        progress_bar(gridSamples);
    return NULL;
}

void * threadgenerateLandmarkSetBlind( void * args )
{
    ThreadGenerateLandmarkSetParams * master =
        (ThreadGenerateLandmarkSetParams *) args;
    ProgTomographAlignment * parent = master->parent;
    int thread_id = master->myThreadID;
    int Nimg=parent->Nimg;
    int numThreads=parent->numThreads;
    const std::vector< std::vector< Matrix2D<double> > > &affineTransformations=
        parent->affineTransformations;
    int gridSamples=parent->gridSamples;

    int deltaShift=(int)floor(XSIZE(*(parent->img)[0])/gridSamples);
    master->chainList=new std::vector<LandmarkChain>;
    Matrix1D<double> rii(3), rjj(3);
    ZZ(rii)=1;
    ZZ(rjj)=1;
    if (thread_id==0)
        init_progress_bar(gridSamples);
    int includedPoints=0;
    int maxSideLength=(parent->blindSeqLength-1)/2;
    for (int nx=thread_id; nx<gridSamples; nx+=numThreads)
    {
        XX(rii)=STARTINGX(*(parent->img)[0])+ROUND(deltaShift*(0.5+nx));
        for (int ny=0; ny<gridSamples; ny+=1)
        {
            YY(rii)=STARTINGY(*(parent->img)[0])+ROUND(deltaShift*(0.5+ny));
            for (int ii=0; ii<=Nimg-1; ++ii)
            {
                // Check if inside the mask
                if (!(*(parent->maskImg[ii]))((int)YY(rii),(int)XX(rii)))
                    continue;
                if (parent->isOutlier(ii))
                    continue;

                LandmarkChain chain;
                chain.clear();
                Landmark l;
                l.x=XX(rii);
                l.y=YY(rii);
                l.imgIdx=ii;
                chain.push_back(l);

                // Follow this landmark backwards
                bool acceptLandmark=true;
                int jj;
                int sideLength=0;
                if (parent->isCapillar)
                    jj=intWRAP(ii-1,0,Nimg-1);
                else
                    jj=ii-1;
                Matrix2D<double> Aij, Aji;
                Matrix1D<double> rcurrent=rii;
                while (jj>=0 && acceptLandmark && sideLength<maxSideLength &&
                       !parent->isOutlier(jj))
                {
                    // Compute the affine transformation between ii and jj
                    int jj_1;
                    if (parent->isCapillar)
                        jj_1=intWRAP(jj+1,0,Nimg-1);
                    else
                        jj_1=jj+1;
                    Aij=affineTransformations[jj][jj_1];
                    Aji=affineTransformations[jj_1][jj];
                    rjj=Aji*rcurrent;
                    int iYYrjj=(int)YY(rjj);
                    int iXXrjj=(int)XX(rjj);
                    if (!(*(parent->maskImg[jj])).outside(iYYrjj,iXXrjj))
                        acceptLandmark=(*(parent->maskImg[jj]))(iYYrjj,iXXrjj);
                    else
                        acceptLandmark=false;
                    if (acceptLandmark)
                    {
                        l.x=XX(rjj);
                        l.y=YY(rjj);
                        l.imgIdx=jj;
                        chain.push_back(l);
                        if (parent->isCapillar)
                            jj=intWRAP(jj-1,0,Nimg-1);
                        else
                            jj=jj-1;
                        rcurrent=rjj;
                        sideLength++;
                    }
                }

                // Follow this landmark forward
                acceptLandmark=true;
                if (parent->isCapillar)
                    jj=intWRAP(ii+1,0,Nimg-1);
                else
                    jj=ii+1;
                rcurrent=rii;
                sideLength=0;
                while (jj<Nimg && acceptLandmark && sideLength<maxSideLength &&
                       !parent->isOutlier(jj))
                {
                    // Compute the affine transformation between ii and jj
                    int jj_1;
                    if (parent->isCapillar)
                        jj_1=intWRAP(jj-1,0,Nimg-1);
                    else
                        jj_1=jj-1;
                    Aij=affineTransformations[jj_1][jj];
                    Aji=affineTransformations[jj][jj_1];
                    rjj=Aij*rcurrent;
                    int iYYrjj=(int)YY(rjj);
                    int iXXrjj=(int)XX(rjj);
                    if (!(*(parent->maskImg[jj])).outside(iYYrjj,iXXrjj))
                        acceptLandmark=(*(parent->maskImg[jj]))(iYYrjj,iXXrjj);
                    else
                        acceptLandmark=false;
                    if (acceptLandmark)
                    {
                        l.x=XX(rjj);
                        l.y=YY(rjj);
                        l.imgIdx=jj;
                        chain.push_back(l);
                        if (parent->isCapillar)
                            jj=intWRAP(jj+1,0,Nimg-1);
                        else
                            jj=jj+1;
                        rcurrent=rjj;
                        sideLength++;
                    }
                }

#ifdef DEBUG
                std::cout << "img=" << ii << " chain length="
                << chain.size() << " [" << jjleft
                << " - " << jjright << "]\n";
#endif

                master->chainList->push_back(chain);
                includedPoints+=chain.size();
            }
#ifdef DEBUG
            std::cout << "Point nx=" << nx << " ny=" << ny
            << " Number of points="
            << includedPoints
            << " Number of chains=" << master->chainList->size()
            << " ( " << ((double) includedPoints)/
            master->chainList->size() << " )\n";
#endif

        }
        if (thread_id==0)
            progress_bar(nx);
    }
    if (thread_id==0)
        progress_bar(gridSamples);
    return NULL;
}

//#define DEBUG
void * threadgenerateLandmarkSetCriticalPoints( void * args )
{
    ThreadGenerateLandmarkSetParams * master =
        (ThreadGenerateLandmarkSetParams *) args;
    ProgTomographAlignment * parent = master->parent;
    int thread_id = master->myThreadID;
    int Nimg=parent->Nimg;
    int numThreads=parent->numThreads;
    const std::vector< std::vector< Matrix2D<double> > > &affineTransformations=
        parent->affineTransformations;

    master->chainList=new std::vector<LandmarkChain>;
    std::vector<LandmarkChain> candidateChainList;
    if (thread_id==0)
        init_progress_bar(Nimg);
    int halfSeqLength=parent->seqLength/2;

    // Design a mask for the dilation
    int radius=4;
    MultidimArray<int> mask;
    mask.resize(2*radius+1,2*radius+1);
    mask.setXmippOrigin();
    BinaryCircularMask(mask,4,OUTSIDE_MASK);

    Image<double> I;
    for (int ii=thread_id; ii<=Nimg-1; ii+=numThreads)
    {
        if (parent->isOutlier(ii))
            continue;
        I.read(parent->name_list[ii]);

        // Generate mask
        MultidimArray<unsigned char> largeMask;
        generateMask(I(),largeMask,
                     XMIPP_MAX(ROUND(parent->localSize*XSIZE(I()))/2,5));

        // Filter the image
        MultidimArray<double> Ifiltered;
        Ifiltered=I();
        Ifiltered.setXmippOrigin();
        FourierFilter FilterBP;
        FilterBP.FilterBand=BANDPASS;
        FilterBP.w1=2.0/XSIZE(Ifiltered);
        FilterBP.w2=128.0/XSIZE(Ifiltered);
        FilterBP.raised_w=1.0/XSIZE(Ifiltered);
        ;
        FilterBP.generateMask(Ifiltered);
        FilterBP.applyMaskSpace(Ifiltered);

        // Identify low valued points and perform dilation
        MultidimArray<double> Iaux=Ifiltered;
        Iaux.selfWindow(
            -ROUND(0.45*YSIZE(Ifiltered)),-ROUND(0.45*XSIZE(Ifiltered)),
            ROUND(0.45*YSIZE(Ifiltered)), ROUND(0.45*XSIZE(Ifiltered)));
        Histogram1D hist;
        compute_hist(Iaux, hist, 400);
        double th=hist.percentil(2);
        std::vector< Matrix1D<double> > Q;
        FOR_ALL_ELEMENTS_IN_ARRAY2D(Iaux)
        if (Ifiltered(i,j)<th && largeMask(i,j))
        {
            // Check if it is a local minimum
            bool localMinimum=true;
            int x0=j;
            int y0=i;
            FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
            {
                int x=x0+j;
                int y=y0+i;
                if (x>=STARTINGX(Iaux) && x<=FINISHINGX(Iaux) &&
                    y>=STARTINGY(Iaux) && y<=FINISHINGY(Iaux))
                    if (A2D_ELEM(Iaux,y,x)<A2D_ELEM(Iaux,y0,x0))
                    {
                        localMinimum=false;
                        break;
                    }
            }
            if (localMinimum)
                Q.push_back(vectorR2(x0,y0));
        }

        // Check that the list is not too long
        if (Q.size()>10*parent->Ncritical)
        {
            size_t qmax=Q.size();
            MultidimArray<double> minDistance;
            minDistance.resize(qmax);
            minDistance.initConstant(1e20);
            for (size_t q1=0; q1<qmax; q1++)
                for (size_t q2=q1+1; q2< qmax; q2++)
                {
                    double diffX=XX(Q[q1])-XX(Q[q2]);
                    double diffY=YY(Q[q1])-YY(Q[q2]);
                    double d12=diffX*diffX+diffY*diffY;
                    if (d12<minDistance(q1))
                        minDistance(q1)=d12;
                    if (d12<minDistance(q2))
                        minDistance(q2)=d12;
                }
            MultidimArray<int> idxDistanceSort;
            minDistance.indexSort(idxDistanceSort);
            std::vector< Matrix1D<double> > Qaux;
            int qlimit=XMIPP_MIN(10*(parent->Ncritical),qmax);
            for (int q=0; q<qlimit; q++)
                Qaux.push_back(Q[idxDistanceSort(qmax-1-q)-1]);
            Q.clear();
            Q=Qaux;
        }

        int qmax=Q.size();
        Matrix1D<double> rii(3), rjj(3);
        ZZ(rii)=1;
        ZZ(rjj)=1;
        MultidimArray<double> corrQ;
        corrQ.initZeros(qmax);
        for (int q=0; q<qmax; q++)
        {
            XX(rii)=XX(Q[q]);
            YY(rii)=YY(Q[q]);

            // Initiate a chain here
            LandmarkChain chain;
            chain.clear();
            Landmark l;
            l.x=XX(rii);
            l.y=YY(rii);
            l.imgIdx=ii;
            chain.push_back(l);

            // Follow this landmark backwards
            int jjmin=XMIPP_MAX(0,ii-halfSeqLength);
            int jjmax=ii-1;
            Matrix2D<double> Aij, Aji;
            Matrix1D<double> rcurrent=rii;
            for (int jj=jjmax; jj>=jjmin; --jj)
            {
                if (parent->isOutlier(jj))
                    break;

                // Compute the affine transformation between jj and jj+1
                int jj_1=jj+1;
                Aij=affineTransformations[jj][jj_1];
                Aji=affineTransformations[jj_1][jj];
                rjj=Aji*rcurrent;
                double corr;
                parent->refineLandmark(jj_1,jj,rcurrent,rjj, corr,true);
                l.x=XX(rjj);
                l.y=YY(rjj);
                l.imgIdx=jj;
                chain.push_back(l);
                rcurrent=rjj;
            }

            // Follow this landmark forwards
            jjmin=ii+1;
            jjmax=XMIPP_MIN(Nimg-1,ii+halfSeqLength);
            rcurrent=rii;
            for (int jj=jjmin; jj<=jjmax; ++jj)
            {
                if (parent->isOutlier(jj))
                    break;

                // Compute the affine transformation between jj-1 and jj
                int jj_1=jj-1;
                Aij=affineTransformations[jj_1][jj];
                Aji=affineTransformations[jj][jj_1];
                rjj=Aij*rcurrent;
                double corr;
                parent->refineLandmark(jj_1,jj,rcurrent,rjj,corr,true);
                l.x=XX(rjj);
                l.y=YY(rjj);
                l.imgIdx=jj;
                chain.push_back(l);
                rcurrent=rjj;
            }

            // Refine chain
            parent->refineChain(chain,corrQ(q));
            candidateChainList.push_back(chain);
        }
        if (thread_id==0)
            progress_bar(ii);

        // Sort all chains according to its correlation
        MultidimArray<int> idx;
        corrQ.indexSort(idx);
        int imax=XMIPP_MIN(XSIZE(idx),parent->Ncritical);
        for (int iq=0; iq<imax; iq++)
        {
            int q=idx(XSIZE(idx)-1-iq)-1;
            if (corrQ(q)>0.5)
            {
                master->chainList->push_back(candidateChainList[q]);
#ifdef DEBUG

                std::cout << "Corr " << iq << ": " << corrQ(q) << ":";
                for (int i=0; i<candidateChainList[q].size(); i++)
                    std::cout << candidateChainList[q][i].imgIdx << " ";
                std::cout << std::endl;
#endif

            }
            else
                std::cout << "It's recommended to reduce the number of critical points\n";
        }
        candidateChainList.clear();

#ifdef DEBUG

        Image<double> save;
        typeCast(*(parent->img[ii]),save());
        save.write("PPPoriginal.xmp");
        save()=Ifiltered;
        save.write("PPPfiltered.xmp");
        double minval=Ifiltered.computeMin();
        FOR_ALL_ELEMENTS_IN_MATRIX2D(Ifiltered)
        if (Ifiltered(i,j)>=th || !largeMask(i,j))
            save(i,j)=th;
        for (int q=0; q<Q.size(); q++)
            FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
            if (YY(Q[q])+i>=STARTINGY(Ifiltered) && YY(Q[q])+i<=FINISHINGY(Ifiltered) &&
                XX(Q[q])+j>=STARTINGX(Ifiltered) && XX(Q[q])+j<=FINISHINGX(Ifiltered))
                save(YY(Q[q])+i,XX(Q[q])+j)=minval;
        imax=XMIPP_MIN(XSIZE(idx),parent->Ncritical);
        for (int iq=0; iq<imax; iq++)
        {
            int q=idx(XSIZE(idx)-1-iq)-1;
            FOR_ALL_ELEMENTS_IN_MATRIX2D(mask)
            if (YY(Q[q])+i>=STARTINGY(Ifiltered) && YY(Q[q])+i<=FINISHINGY(Ifiltered) &&
                XX(Q[q])+j>=STARTINGX(Ifiltered) && XX(Q[q])+j<=FINISHINGX(Ifiltered))
                save(YY(Q[q])+i,XX(Q[q])+j)=(minval+th)/2;
        }

        save.write("PPPcritical.xmp");
        std::cout << "Number of critical points=" << Q.size() << std::endl;
        std::cout << "CorrQ stats:";
        corrQ.printStats();
        std::cout << std::endl;
        std::cout << "Press any key\n";
        char c;
        std::cin >> c;
#endif

    }
    if (thread_id==0)
        progress_bar(Nimg);
    return NULL;
}
#undef DEBUG

ProgTomographAlignment::~ProgTomographAlignment()
{
    // Clear the list of images if not empty
    if (!img.empty())
    {
        for (size_t i=0; i<img.size(); i++)
            delete img[i];
        img.clear();
    }

    if (!maskImg.empty())
    {
        for (size_t i=0; i<maskImg.size(); i++)
            delete maskImg[i];
        maskImg.clear();
    }
}

void ProgTomographAlignment::generateLandmarkSet()
{
    FileName fn_tmp = fnRoot+"_landmarks.txt";
    if (!fn_tmp.exists())
    {
        pthread_t * th_ids = new pthread_t[numThreads];
        ThreadGenerateLandmarkSetParams * th_args=
            new ThreadGenerateLandmarkSetParams[numThreads];
        for( int nt = 0 ; nt < numThreads ; nt ++ )
        {
            th_args[nt].parent = this;
            th_args[nt].myThreadID = nt;
            if (useCriticalPoints)
                pthread_create( (th_ids+nt) , NULL, threadgenerateLandmarkSetCriticalPoints, (void *)(th_args+nt) );
            else
                pthread_create( (th_ids+nt) , NULL, threadgenerateLandmarkSetGrid, (void *)(th_args+nt) );
        }

        std::vector<LandmarkChain> chainList;
        int includedPoints=0;
        for( int nt = 0 ; nt < numThreads ; nt ++ )
        {
            pthread_join(*(th_ids+nt), NULL);
            int imax=th_args[nt].chainList->size();
            for (int i=0; i<imax; i++)
            {
                chainList.push_back((*th_args[nt].chainList)[i]);
                includedPoints+=(*th_args[nt].chainList)[i].size();
            }
            delete th_args[nt].chainList;
        }
        delete( th_ids );
        delete( th_args );

        // Add blind landmarks
        if (blindSeqLength>0)
        {
            th_ids  = new pthread_t[numThreads];
            th_args = new ThreadGenerateLandmarkSetParams[numThreads];
            for( int nt = 0 ; nt < numThreads ; nt ++ )
            {
                th_args[nt].parent = this;
                th_args[nt].myThreadID = nt;
                pthread_create( (th_ids+nt) , NULL, threadgenerateLandmarkSetBlind, (void *)(th_args+nt) );
            }

            for( int nt = 0 ; nt < numThreads ; nt ++ )
            {
                pthread_join(*(th_ids+nt), NULL);
                int imax=th_args[nt].chainList->size();
                for (int i=0; i<imax; i++)
                {
                    chainList.push_back((*th_args[nt].chainList)[i]);
                    includedPoints+=(*th_args[nt].chainList)[i].size();
                }
                delete th_args[nt].chainList;
            }
            delete( th_ids );
            delete( th_args );
        }

        // Generate the landmark "matrix"
        allLandmarksX.resize(chainList.size(),Nimg);
        allLandmarksY.resize(chainList.size(),Nimg);
        allLandmarksX.initConstant(XSIZE(*img[0]));
        allLandmarksY.initConstant(YSIZE(*img[0]));
        for (size_t i=0; i<chainList.size(); i++)
        {
            for (size_t j=0; j<chainList[i].size(); j++)
            {
                int idx=chainList[i][j].imgIdx;
                allLandmarksX(i,idx)=chainList[i][j].x;
                allLandmarksY(i,idx)=chainList[i][j].y;
            }
        }

        // Write landmarks
        if (includedPoints>0)
        {
            writeLandmarkSet(fnRoot+"_landmarks.txt");
            std::cout << " Number of points="
            << includedPoints
            << " Number of chains=" << chainList.size()
            << " ( " << ((double) includedPoints)/chainList.size() << " )\n";
        }
        else
            REPORT_ERROR(ERR_VALUE_INCORRECT,"There are no landmarks meeting this threshold. Try to lower it");
    }
    else
    {
        readLandmarkSet(fnRoot+"_landmarks.txt");
    }
}
#undef DEBUG

/* Refine landmark --------------------------------------------------------- */
bool ProgTomographAlignment::refineLandmark(int ii, int jj,
        const Matrix1D<double> &rii, Matrix1D<double> &rjj, double &maxCorr,
        bool tryFourier) const
{
    maxCorr=-1;
    int halfSize=XMIPP_MAX(ROUND(localSize*XSIZE(*img[ii]))/2,5);
    if (XX(rii)-halfSize<STARTINGX(*img[ii])  ||
        XX(rii)+halfSize>FINISHINGX(*img[ii]) ||
        YY(rii)-halfSize<STARTINGY(*img[ii])  ||
        YY(rii)+halfSize>FINISHINGY(*img[ii]) ||
        XX(rjj)-halfSize<STARTINGX(*img[jj])  ||
        XX(rjj)+halfSize>FINISHINGX(*img[jj]) ||
        YY(rjj)-halfSize<STARTINGY(*img[jj])  ||
        YY(rjj)+halfSize>FINISHINGY(*img[jj]))
        return 0;

    // Check if the two pieces are reversed
    bool reversed=isCapillar && ABS(ii-jj)>Nimg/2;

    // Select piece in image ii, compute its statistics and normalize
    MultidimArray<double> pieceii(2*halfSize+1,2*halfSize+1);
    pieceii.setXmippOrigin();
    const MultidimArray<unsigned char> &Iii=(*img[ii]);
    FOR_ALL_ELEMENTS_IN_ARRAY2D(pieceii)
    A2D_ELEM(pieceii,i,j)=A2D_ELEM(Iii,
                                   (int)(YY(rii)+i),(int)(XX(rii)+j));
    if (showRefinement)
    {
        Image<double> save;
        save()=pieceii;
        save.write("PPPpieceii.xmp");
        std::cout << "ii=" << ii << " jj=" << jj << std::endl;
        std::cout << "rii=" << rii.transpose() << std::endl;
    }

    // Choose threshold
    double actualCorrThreshold=corrThreshold;
    if (actualCorrThreshold<0 && VEC_XSIZE(avgForwardPatchCorr)>0)
    {
        if (ii<jj)
            actualCorrThreshold=XMIPP_MIN(avgForwardPatchCorr(ii),
                                          avgBackwardPatchCorr(jj));
        else
            actualCorrThreshold=XMIPP_MIN(avgBackwardPatchCorr(ii),
                                          avgForwardPatchCorr(jj));
    }
    if (showRefinement)
        std::cout << "actualCorrThreshold=" << actualCorrThreshold << std::endl;

    // Try Fourier
    if (tryFourier)
    {
        const MultidimArray<unsigned char> &Ijj=(*img[jj]);
        if (XX(rjj)-halfSize>=STARTINGX(Ijj) &&
            XX(rjj)+halfSize<=FINISHINGX(Ijj) &&
            YY(rjj)-halfSize>=STARTINGY(Ijj) &&
            YY(rjj)+halfSize<=FINISHINGY(Ijj))
        {
            // Take the piece at jj
            MultidimArray<double> piecejj(2*halfSize+1,2*halfSize+1);
            piecejj.setXmippOrigin();
            FOR_ALL_ELEMENTS_IN_ARRAY2D(piecejj)
            A2D_ELEM(piecejj,i,j)=A2D_ELEM(Ijj,
                                           (int)(YY(rjj)+i),(int)(XX(rjj)+j));

            // Compute its correlation with pieceii
            double corrOriginal=correlationIndex(pieceii,piecejj);
            if (showRefinement)
            {
                Image<double> save;
                save()=piecejj;
                save.write("PPPpiecejjOriginal.xmp");
                std::cout << "Corr original=" << corrOriginal << std::endl;
            }

            // Now try with the best shift
            double shiftX,shiftY;
            CorrelationAux aux;
            bestNonwrappingShift(pieceii,piecejj,shiftX,shiftY,aux);
            Matrix1D<double> fftShift(2);
            VECTOR_R2(fftShift,shiftX,shiftY);
            selfTranslate(LINEAR,piecejj,fftShift,WRAP);
            double corrFFT=correlationIndex(pieceii,piecejj);
            if (corrFFT>corrOriginal)
            {
                XX(rjj)-=shiftX;
                YY(rjj)-=shiftY;
            }

            if (showRefinement)
            {
                Image<double> save;
                save()=piecejj;
                save.write("PPPpiecejjFFT.xmp");
                std::cout << "FFT shift=" << fftShift.transpose() << std::endl;
                std::cout << "Corr FFT=" << corrFFT << std::endl;
                if (corrFFT>corrOriginal)
                {
                    if (XX(rjj)-halfSize>=STARTINGX(Ijj) &&
                        XX(rjj)+halfSize<=FINISHINGX(Ijj) &&
                        YY(rjj)-halfSize>=STARTINGY(Ijj) &&
                        YY(rjj)+halfSize<=FINISHINGY(Ijj))
                    {
                        FOR_ALL_ELEMENTS_IN_ARRAY2D(piecejj)
                        A2D_ELEM(piecejj,i,j)=A2D_ELEM(Ijj,
                                                       (int)(YY(rjj)+i),(int)(XX(rjj)+j));
                        save()=piecejj;
                        save.write("PPPpiecejjNew.xmp");
                    }
                }
            }
        }
    }

    bool retval=refineLandmark(pieceii,jj,rjj,actualCorrThreshold,
                               reversed,maxCorr);
    return retval;
}

bool ProgTomographAlignment::refineLandmark(const MultidimArray<double> &pieceii,
        int jj, Matrix1D<double> &rjj, double actualCorrThreshold,
        bool reversed, double &maxCorr) const
{
    int halfSize=XSIZE(pieceii)/2;

    double mean_ii=0, stddev_ii=0;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(pieceii)
    {
        mean_ii+=A2D_ELEM(pieceii,i,j);
        stddev_ii+=A2D_ELEM(pieceii,i,j)*A2D_ELEM(pieceii,i,j);
    }
    mean_ii/=MULTIDIM_SIZE(pieceii);
    stddev_ii = stddev_ii / MULTIDIM_SIZE(pieceii) - mean_ii * mean_ii;
    stddev_ii *= MULTIDIM_SIZE(pieceii) / (MULTIDIM_SIZE(pieceii) - 1);
    stddev_ii = sqrt(static_cast<double>((ABS(stddev_ii))));
    if (stddev_ii>XMIPP_EQUAL_ACCURACY)
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pieceii)
        DIRECT_MULTIDIM_ELEM(pieceii,n)=
            (DIRECT_MULTIDIM_ELEM(pieceii,n)-mean_ii)/stddev_ii;

    // Try all possible shifts
    MultidimArray<double> corr((int)(1.5*(2*halfSize+1)),(int)(1.5*(2*halfSize+1)));
    corr.setXmippOrigin();
    corr.initConstant(-1.1);
    bool accept=false;
    double maxval=-1;
    MultidimArray<double> piecejj(2*halfSize+1,2*halfSize+1);
    piecejj.setXmippOrigin();
    if (stddev_ii>XMIPP_EQUAL_ACCURACY)
    {
        int imax=0, jmax=0;
        std::queue< Matrix1D<double> > Q;
        Q.push(vectorR2(0,0));
        while (!Q.empty())
        {
            // Get the first position to evaluate
            int shifty=(int)YY(Q.front());
            int shiftx=(int)XX(Q.front());
            Q.pop();
            if (!corr(shifty,shiftx)<-1)
                continue;

            // Select piece in image jj, compute its statistics and normalize
            double mean_jj=0, stddev_jj=0;
            const MultidimArray<unsigned char> &Ijj=(*img[jj]);
            if (XX(rjj)+shiftx+STARTINGX(piecejj)<STARTINGX(Ijj) ||
                YY(rjj)+shifty+STARTINGY(piecejj)<STARTINGY(Ijj) ||
                XX(rjj)+shiftx+FINISHINGX(piecejj)>FINISHINGX(Ijj) ||
                YY(rjj)+shifty+FINISHINGY(piecejj)>FINISHINGY(Ijj))
                continue;
            FOR_ALL_ELEMENTS_IN_ARRAY2D(piecejj)
            {
                double pixval=A2D_ELEM(Ijj,
                                       (int)(YY(rjj)+shifty+i),
                                       (int)(XX(rjj)+shiftx+j));
                A2D_ELEM(piecejj,i,j)=pixval;
                mean_jj+=pixval;
                stddev_jj+=pixval*pixval;
            }
            mean_jj/=MULTIDIM_SIZE(piecejj);
            stddev_jj = stddev_jj / MULTIDIM_SIZE(piecejj) - mean_jj * mean_jj;
            stddev_jj *= MULTIDIM_SIZE(piecejj) / (MULTIDIM_SIZE(piecejj) - 1);
            stddev_jj = sqrt(static_cast<double>((ABS(stddev_jj))));
            if (reversed)
                piecejj.selfReverseY();

            // Compute the correlation
            corr(shifty,shiftx)=0;
            double &corrRef=corr(shifty,shiftx);
            if (stddev_jj>XMIPP_EQUAL_ACCURACY)
            {
                double istddev_jj=1.0/stddev_jj;
                FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(piecejj)
                corrRef+=
                    DIRECT_MULTIDIM_ELEM(pieceii,n)*
                    (DIRECT_MULTIDIM_ELEM(piecejj,n)-mean_jj)*istddev_jj;
                corrRef/=MULTIDIM_SIZE(piecejj);
            }

            if (corrRef>maxval)
            {
                maxval=corrRef;
                imax=shifty;
                jmax=shiftx;
                for (int step=1; step<=5; step+=2)
                    for (int stepy=-1; stepy<=1; stepy++)
                        for (int stepx=-1; stepx<=1; stepx++)
                        {
                            int newshifty=shifty+stepy*step;
                            int newshiftx=shiftx+stepx*step;
                            if (newshifty>=STARTINGY(corr) &&
                                newshifty<=FINISHINGY(corr) &&
                                newshiftx>=STARTINGX(corr) &&
                                newshiftx<=FINISHINGX(corr) &&
                                (XX(rjj)+newshiftx-halfSize)>=STARTINGX(Ijj) &&
                                (XX(rjj)+newshiftx+halfSize)<=FINISHINGX(Ijj) &&
                                (YY(rjj)+newshifty-halfSize)>=STARTINGY(Ijj) &&
                                (YY(rjj)+newshifty+halfSize)<=FINISHINGY(Ijj))
                                if (corr(newshifty,newshiftx)<-1)
                                    Q.push(vectorR2(newshiftx,newshifty));
                        }
            }
        }

        if (maxval>actualCorrThreshold)
        {
            XX(rjj)+=jmax;
            if (reversed)
                YY(rjj)-=imax;
            else
                YY(rjj)+=imax;
            accept=true;
        }
        if (showRefinement)
        {
            Image<double> save;
            FOR_ALL_ELEMENTS_IN_ARRAY2D(piecejj)
            piecejj(i,j)=(*img[jj])((int)(YY(rjj)+i),(int)(XX(rjj)+j));
            if (reversed)
                piecejj.selfReverseY();
            save()=piecejj;
            save.write("PPPpiecejj.xmp");
            save()=corr;
            save.write("PPPcorr.xmp");
            std::cout << "jj=" << jj << " rjj=" << rjj.transpose() << std::endl;
            std::cout << "imax=" << imax << " jmax=" << jmax << std::endl;
            std::cout << "maxval=" << maxval << std::endl;
        }
    }
    if (showRefinement)
    {
        std::cout << "Accepted=" << accept << std::endl;
        std::cout << "Press any key\n";
        char c;
        std::cin >> c;
    }
    maxCorr=maxval;
    return (accept);
}

/* Refine chain ------------------------------------------------------------ */
//#define DEBUG
bool ProgTomographAlignment::refineChain(LandmarkChain &chain,
        double &corrChain)
{
#ifdef DEBUG
    std::cout << "Chain for refinement: ";
    for (int i=0; i<chain.size(); i++)
        std::cout << chain[i].imgIdx << " ";
    std::cout << std::endl;
#endif

    for (int K=0; K<2 && (chain.size()>seqLength || useCriticalPoints); K++)
    {
        sort(chain.begin(), chain.end());
        Matrix1D<double> rii(2), rjj(2), newrjj(2);

        if (useCriticalPoints)
        {
            // compute average piece
            int halfSize=XMIPP_MAX(ROUND(localSize*XSIZE(*img[0]))/2,5);
            int chainLength=chain.size();
            MultidimArray<double> avgPiece(2*halfSize+1,2*halfSize+1), pieceAux;
            avgPiece.setXmippOrigin();
            pieceAux=avgPiece;
            for (int n=0; n<chainLength; n++)
            {
                int ii=chain[n].imgIdx;
                VECTOR_R2(rii,chain[n].x,chain[n].y);
                if (XX(rii)-halfSize<STARTINGX(*img[ii])  ||
                    XX(rii)+halfSize>FINISHINGX(*img[ii]) ||
                    YY(rii)-halfSize<STARTINGY(*img[ii])  ||
                    YY(rii)+halfSize>FINISHINGY(*img[ii]))
                {
                    corrChain=0;
                    return false;
                }
                FOR_ALL_ELEMENTS_IN_ARRAY2D(avgPiece)
                pieceAux(i,j)=(*img[ii])((int)(YY(rii)+i),(int)(XX(rii)+j));
                if (isCapillar && ABS(chain[n].imgIdx-chain[0].imgIdx)>Nimg/2)
                    pieceAux.selfReverseY();
                avgPiece+=pieceAux;
            }
            avgPiece/=chainLength;
#ifdef DEBUG

            Image<double> save;
            save()=avgPiece;
            save.write("PPPavg.xmp");
            std::cout << "Average " << K << " -------------------  \n";
            showRefinement=true;
#endif

            // Align all images with respect to this average
            corrChain=2;
            for (int j=0; j<chainLength; j++)
            {
                int jj=chain[j].imgIdx;
                VECTOR_R2(rjj,chain[j].x,chain[j].y);
                double corr;
                bool accepted=refineLandmark(avgPiece,jj,rjj,0,false,corr);
                if (accepted)
                {
                    chain[j].x=XX(rjj);
                    chain[j].y=YY(rjj);
                    corrChain=XMIPP_MIN(corrChain,corr);
                }
            }

#ifdef DEBUG
            showRefinement=false;
#endif

        }
        else
        {
            // Refine every step
            for (int step=2; step<=maxStep; step++)
            {
                int chainLength=chain.size();

                // Refine forwards every step
                int ileft=-1;
                for (int i=0; i<chainLength-step; i++)
                {
                    int ii=chain[i].imgIdx;
                    VECTOR_R2(rii,chain[i].x,chain[i].y);
                    int jj=chain[i+step].imgIdx;
                    VECTOR_R2(rjj,chain[i+step].x,chain[i+step].y);
                    newrjj=rjj;
                    double corr;
                    bool accepted=refineLandmark(ii,jj,rii,newrjj,corr,false);
                    if (((newrjj-rjj).module()<4 && accepted) || useCriticalPoints)
                    {
                        chain[i+step].x=XX(newrjj);
                        chain[i+step].y=YY(newrjj);
                    }
                    else
                        ileft=i;
                }

#ifdef DEBUG
                // COSS showRefinement=(step==maxStep);
#endif
                // Refine backwards all images
                int iright=chainLength;
                corrChain=2;
                for (int i=chainLength-1; i>=1; i--)
                {
                    int ii=chain[i].imgIdx;
                    VECTOR_R2(rii,chain[i].x,chain[i].y);
                    int jj=chain[i-1].imgIdx;
                    VECTOR_R2(rjj,chain[i-1].x,chain[i-1].y);
                    newrjj=rjj;
                    double corr;
                    bool accepted=refineLandmark(ii,jj,rii,newrjj,corr,false);
                    corrChain=XMIPP_MIN(corrChain,corr);
                    if (((newrjj-rjj).module()<4 && accepted) || useCriticalPoints)
                    {
                        chain[i-1].x=XX(newrjj);
                        chain[i-1].y=YY(newrjj);
                    }
                    else
                        iright=i;
                }
#ifdef DEBUG
                // COSS: showRefinement=false;
#endif

                LandmarkChain refinedChain;
                for (int ll=ileft+1; ll<=iright-1; ll++)
                    refinedChain.push_back(chain[ll]);
                chain=refinedChain;
                double tilt0=tiltList[chain[0].imgIdx];
                double tiltF=tiltList[chain[chain.size()-1].imgIdx];
                double lengthThreshold=XMIPP_MAX(3.,FLOOR(seqLength*cos(DEG2RAD(0.5*(tilt0+tiltF)))));
                if (chain.size()<lengthThreshold && !useCriticalPoints)
                    return false;
            }
        }
    }

    // Check that the chain is of the desired length
    double tilt0=tiltList[chain[0].imgIdx];
    double tiltF=tiltList[chain[chain.size()-1].imgIdx];
    double lengthThreshold=XMIPP_MAX(3,FLOOR(seqLength*cos(DEG2RAD(0.5*(tilt0+tiltF)))));
    return chain.size()>lengthThreshold;
}
#undef DEBUG

/* Read/Write landmark set ------------------------------------------------- */
void ProgTomographAlignment::writeLandmarkSet(const FileName &fnLandmark) const
{
    std::ofstream fhOut;
    fhOut.open(fnLandmark.c_str());
    if (!fhOut)
        REPORT_ERROR(ERR_IO_NOWRITE,(std::string)"Cannot open "+fnLandmark+" for output");
    fhOut << "Point     x       y       slice   color "
    << MAT_YSIZE(allLandmarksX) << " " << MAT_XSIZE(allLandmarksX) << std::endl;
    for (size_t i=0; i<MAT_XSIZE(allLandmarksX); i++)
    {
        int counter=0;
        for (size_t j=0; j<MAT_YSIZE(allLandmarksX); j++)
            if (allLandmarksX(j,i)!=XSIZE(*img[0]))
            {
                fhOut << counter << " \t"
                << ROUND(allLandmarksX(j,i)-STARTINGX(*img[0])) << " \t"
                << ROUND(allLandmarksY(j,i)-STARTINGY(*img[0])) << " \t"
                << i+1 << " \t" << j << std::endl;
                counter++;
            }
    }
    fhOut.close();
}

void ProgTomographAlignment::readLandmarkSet(const FileName &fnLandmark)
{
    std::ifstream fhIn;
    fhIn.open(fnLandmark.c_str());
    if (!fhIn)
        REPORT_ERROR(ERR_IO_NOTEXIST,(std::string)"Cannot open "+fnLandmark+" for input");
    std::string dummyStr;
    int Nlandmark;
    fhIn >> dummyStr >> dummyStr >> dummyStr >> dummyStr >> dummyStr
    >> Nlandmark >> Nimg;
    if (Nlandmark<=0)
        REPORT_ERROR(ERR_VALUE_INCORRECT,(std::string)"No landmarks are found in "+fnLandmark);
    allLandmarksX.resize(Nlandmark,Nimg);
    allLandmarksY.resize(Nlandmark,Nimg);
    allLandmarksX.initConstant(XSIZE(*img[0]));
    allLandmarksY.initConstant(XSIZE(*img[0]));
    fhIn.exceptions ( std::ifstream::eofbit |
                      std::ifstream::failbit |
                      std::ifstream::badbit );
    while (!fhIn.eof())
    {
        try
        {
            int dummyInt, x, y, i, j;
            fhIn >> dummyInt >> x >> y >> i >> j;
            i=i-1;
            allLandmarksX(j,i)=x+STARTINGX(*img[0]);
            allLandmarksY(j,i)=y+STARTINGY(*img[0]);
        }
        catch (std::ifstream::failure e)
        {
            // Do nothing with this line
        }
    }
    fhIn.close();
    std::cout << "The file " << fnLandmark << " has been read for the landmarks\n"
    << Nlandmark << " landmarks are read\n";
}

/* Write affine transformations -------------------------------------------- */
void ProgTomographAlignment::writeTransformations(
    const FileName &fnTransformations) const
{
    std::ofstream fhOut;
    fhOut.open(fnTransformations.c_str());
    if (!fhOut)
        REPORT_ERROR(ERR_IO_NOWRITE,(std::string)"Cannot open "+fnTransformations+" for output");
    int imax=affineTransformations.size();
    int counter=0;
    for (int i=0; i<imax; i++)
    {
        int jmax=affineTransformations[i].size();
        for (int j=i+1; j<jmax; j++)
            if (MAT_XSIZE(affineTransformations[i][j])!=0)
            {
                fhOut << tiltList[i] << "\t0.0\t"
                << affineTransformations[i][j](0,0) << "\t"
                << affineTransformations[i][j](1,0) << "\t"
                << affineTransformations[i][j](0,1) << "\t"
                << affineTransformations[i][j](1,1) << "\t"
                << affineTransformations[i][j](0,2) << "\t"
                << affineTransformations[i][j](1,2) << std::endl;
                counter++;
            }
    }
    if (counter==imax-1)
    {
        fhOut << tiltList[tiltList.size()-1] << "\t0.0\t1.0\t0.0\t0.0\t1.0\t0.0\t0.0\n";
        fhOut << "1 0 0 1 0 0\n";
    }
    fhOut.close();
}

/* Align images ------------------------------------------------------------ */
//#define DEBUG
void ProgTomographAlignment::alignImages(const Alignment &alignment)
{
    // Correct all landmarks
    Matrix1D<double> r(2);
    for (size_t i=0; i<MAT_XSIZE(allLandmarksX); i++)
    {
        Matrix2D<double> R;
        rotation2DMatrix(90-alignment.rot+alignment.psi(i),R,false);
        for (size_t j=0; j<MAT_YSIZE(allLandmarksX); j++)
            if (allLandmarksX(j,i)!=XSIZE(*img[0]))
            {
                VECTOR_R2(r,allLandmarksX(j,i),allLandmarksY(j,i));
                r=R*(r-alignment.di[i]+alignment.diaxis[i]);
                allLandmarksX(j,i)=XX(r);
                allLandmarksY(j,i)=YY(r);
            }
    }

    // Compute the average height of the 3D landmarks seen at 0 degrees
    Matrix1D<double> axis;
    Euler_direction(alignment.rot, alignment.tilt, 0, axis);
    double z0=0, z0N=0;
    for (size_t j=0; j<MAT_YSIZE(allLandmarksX); j++)
        if (allLandmarksX(j,iMinTilt)!=XSIZE(*img[0]))
        {
            Matrix2D<double> Raxismin, Rmin, RtiltYmin;
            rotation3DMatrix(tiltList[iMinTilt],axis,Raxismin,false);
            rotation2DMatrix(90-alignment.rot+alignment.psi(iMinTilt),Rmin);
            rotation3DMatrix(-tiltList[iMinTilt],'Y',RtiltYmin,false);
            Matrix1D<double> rjp=RtiltYmin*Rmin*Raxismin*alignment.rj[j];
            z0+=ZZ(rjp);
            z0N++;
        }
    if (z0N==0)
        REPORT_ERROR(ERR_VALUE_INCORRECT,"There is no landmark at 0 degrees");
    z0/=z0N;
    std::cout << "Average height of the landmarks at 0 degrees=" << z0 << std::endl;
    MetaData DF;

    MDIterator * iter = NULL;

    if (!fnSelOrig.empty())
    {
        iter = new MDIterator(SForig);
    }
    DF.setComment("First shift by -(shiftX,shiftY), then rotate by psi");

    MultidimArray<double> mask;
    MultidimArray<int> iMask;
    Image<double> I;
    Matrix2D<double> M, M1, M2, M3;
    FileName fn_corrected;

    for (int n=0;n<Nimg; n++)
    {
        // Align the normal image
        I.read(name_list[n]);
        mask.initZeros(I());
        FOR_ALL_ELEMENTS_IN_ARRAY2D(I())
        if (I(i,j)!=0)
            mask(i,j)=1;
        translation2DMatrix(vectorR2(-z0*sin(DEG2RAD(tiltList[n])),0),M1);
        rotation2DMatrix(90-alignment.rot+alignment.psi(n),M2);
        translation2DMatrix(-(alignment.di[n]+alignment.diaxis[n]),M3);
        M=M1*M2*M3;
        selfApplyGeometry(BSPLINE3,I(),M,IS_NOT_INV,DONT_WRAP);
        selfApplyGeometry(LINEAR,mask,M,IS_NOT_INV,DONT_WRAP);
        mask.binarize(0.5);
        typeCast(mask,iMask);
        double minval, maxval, avg, stddev;
        computeStats_within_binary_mask(iMask,I(),minval, maxval, avg, stddev);
        FOR_ALL_ELEMENTS_IN_ARRAY2D(iMask)
        if (iMask(i,j)==0)
            I(i,j)=0;
        else if (!dontNormalize)
            I(i,j)-=avg;
        double rot=0;
        double tilt=tiltList[n];
        double psi=0;
        I.setEulerAngles(rot, tilt, psi);
        fn_corrected.compose(n+1, fnRoot+"_corrected_", "stk");
        I.write(fn_corrected);

        // Align the original image
        if (fnSelOrig!="")
        {
            FileName auxFn;
            SForig.getValue( MDL_IMAGE, auxFn, iter->objId);
            Image<double> Iorig;
            Iorig.read( auxFn );
            //SForig.nextObject();
            iter->moveNext();
            mask.initZeros(Iorig());
            FOR_ALL_ELEMENTS_IN_ARRAY2D(Iorig())
            if (Iorig(i,j)!=0)
                mask(i,j)=1;
            translation2DMatrix(vectorR2(-z0*sin(DEG2RAD(tiltList[n]))*
                                         (((double)XSIZE(Iorig()))/XSIZE(I())),0),M1);
            rotation2DMatrix(90-alignment.rot+alignment.psi(n),M2);
            translation2DMatrix(-(alignment.di[n]+alignment.diaxis[n])*
                                (((double)XSIZE(Iorig()))/XSIZE(I())),M3);
            M=M1*M2*M3;
            selfApplyGeometry(BSPLINE3,Iorig(),M,IS_NOT_INV,DONT_WRAP);
            selfApplyGeometry(LINEAR,mask,M,IS_NOT_INV,DONT_WRAP);
            mask.binarize(0.5);
            typeCast(mask,iMask);
            computeStats_within_binary_mask(iMask,Iorig(),minval, maxval,
                                            avg, stddev);
            FOR_ALL_ELEMENTS_IN_ARRAY2D(iMask)
            if (iMask(i,j)==0)
                Iorig(i,j)=0;
            else if (!dontNormalize)
                Iorig(i,j)-=avg;
            Iorig.setEulerAngles(rot, tilt, psi);
            fn_corrected.compose(n+1, fnRoot+"_corrected_originalsize_", "stk");
            Iorig.write(fn_corrected);
        }

        // Prepare data for the docfile
        size_t id = DF.addObject();
        DF.setValue(MDL_IMAGE, fn_corrected, id);
        DF.setValue(MDL_ANGLE_PSI, 90.-alignment.rot+alignment.psi(n), id);
        DF.setValue(MDL_SHIFT_X, XX(alignment.di[n]+alignment.diaxis[n]), id);
        DF.setValue(MDL_SHIFT_Y, YY(alignment.di[n]+alignment.diaxis[n]), id);
    }
    DF.write(fnRoot+"_correction_parameters.txt");
    delete iter;
#ifdef DEBUG

    Image<double> save;
    save().initZeros(*img[0]);
    save().setXmippOrigin();
    for (int j=0; j<YSIZE(allLandmarksX); j++)
        if (allLandmarksX(j,iMinTilt)!=XSIZE(*img[0]))
        {
            Matrix2D<double> Raxismin=rotation3DMatrix(tiltList[iMinTilt],axis);
            Raxismin.resize(3,3);
            Matrix2D<double> Rmin=rotation2DMatrix(90-alignment.rot+alignment.psi(iMinTilt));
            Matrix2D<double> RtiltYmin=rotation3DMatrix(-tiltList[iMinTilt],'Y');
            RtiltYmin.resize(3,3);
            Matrix1D<double> rjp=RtiltYmin*Rmin*Raxismin*alignment.rj[j];
            std::cout << rjp.transpose() << std::endl;
            for (int i=0; i<XSIZE(allLandmarksX); i++)
                if (allLandmarksX(j,i)!=XSIZE(*img[0]))
                {
                    save(allLandmarksY(j,i),allLandmarksX(j,i))=ABS(i-iMinTilt);

                    /*
                    Matrix2D<double> Raxis=rotation3DMatrix(tiltList[i],axis);
                    Raxis.resize(3,3);
                    Matrix2D<double> R=rotation2DMatrix(90-alignment.rot+alignment.psi(i));
                    Matrix1D<double> p=R*Raxis*alignment.rj[j];
                    if (YY(p)>=STARTINGY(save()) && YY(p)<=FINISHINGY(save()) &&
                        XX(p)>=STARTINGX(save()) && XX(p)<=FINISHINGX(save()))
                       save(YY(p),XX(p))=-ABS(i-iMinTilt);
                    */
                    Matrix2D<double> RtiltY=rotation3DMatrix(-tiltList[i],'Y');
                    RtiltY.resize(3,3);
                    Matrix1D<double> p=RtiltY*rjp;
                    if (YY(p)>=STARTINGY(save()) && YY(p)<=FINISHINGY(save()) &&
                        XX(p)>=STARTINGX(save()) && XX(p)<=FINISHINGX(save()))
                        save(YY(p),XX(p))=-ABS(i-iMinTilt);
                }
        }
    save.write("PPPmovementsLandmarks0.xmp");
#endif
}
#undef DEBUG

/* Produce information from landmarks -------------------------------------- */
void ProgTomographAlignment::produceInformationFromLandmarks()
{
    // Produce V sets
    std::vector<int> emptyVector;

    Vseti.clear();
    for (size_t i=0; i<MAT_XSIZE(allLandmarksX); i++)
        Vseti.push_back(emptyVector);

    Vsetj.clear();
    for (size_t j=0; j<MAT_YSIZE(allLandmarksX); j++)
        Vsetj.push_back(emptyVector);

    for (size_t j=0; j<MAT_YSIZE(allLandmarksX); j++)
        for (size_t i=0; i<MAT_XSIZE(allLandmarksX); i++)
            if (allLandmarksX(j,i)!=XSIZE(*img[0]))
            {
                Vseti[i].push_back(j);
                Vsetj[j].push_back(i);
            }

    // Count the number of landmarks per image and average projection
    // of all landmarks in a given image
    ni.initZeros(Nimg);
    barpi.clear();
    for (int i=0; i<Nimg; i++)
    {
        ni(i)=static_cast<int>(Vseti[i].size());
        Matrix1D<double> pi(2);
        for (int jj=0; jj<ni(i); jj++)
        {
            int j=Vseti[i][jj];
            XX(pi)+=allLandmarksX(j,i);
            YY(pi)+=allLandmarksY(j,i);
        }
        pi/=ni(i);
        barpi.push_back(pi);
    }
}

/* Remove outliers ---------------------------------------------------------*/
void ProgTomographAlignment::removeOutlierLandmarks(
    const Alignment &alignment)
{
    std::cout << "Removing outliers ...\n";

    // Compute threshold for outliers
    Histogram1D hist;
    compute_hist(alignment.errorLandmark, hist, 100);
    // double threshold0=hist.percentil(10);
    double thresholdF=hist.percentil(90);

    // Identify outliers
    std::vector<bool> validLandmark;
    int invalidCounter=0;
    int Nlandmark=MAT_YSIZE(allLandmarksX);
    for (int j=0; j<Nlandmark; j++)
        if (//COSS alignment.errorLandmark(j)<threshold0 ||
            alignment.errorLandmark(j)>thresholdF)
        {
            validLandmark.push_back(false);
            invalidCounter++;
        }
        else
            validLandmark.push_back(true);

    // Remove outliers
    Matrix2D<double> newAllLandmarksX(Nlandmark-invalidCounter,Nimg);
    Matrix2D<double> newAllLandmarksY(Nlandmark-invalidCounter,Nimg);
    int jj=0;
    for (int j=0; j<Nlandmark; j++)
    {
        if (!validLandmark[j])
            continue;
        for (int i=0; i<Nimg; i++)
        {
            newAllLandmarksX(jj,i)=allLandmarksX(j,i);
            newAllLandmarksY(jj,i)=allLandmarksY(j,i);
        }
        jj++;
    }
    allLandmarksX=newAllLandmarksX;
    allLandmarksY=newAllLandmarksY;

    std::cout << invalidCounter << " out of " << Nlandmark
    << " landmarks have been removed\n";
}

/* Run --------------------------------------------------------------------- */
namespace TomographAlignment
{
const ProgTomographAlignment* global_prm;
}

double wrapperError(double *p, void *prm)
{
    Alignment alignment(TomographAlignment::global_prm);
    alignment=*(TomographAlignment::global_prm->bestPreviousAlignment);
    alignment.rot=p[1];
    alignment.tilt=p[2];
    return alignment.optimizeGivenAxisDirection();
}

#define DEBUG
void ProgTomographAlignment::run()
{
    produceSideInfo();
    std::cerr << "generateLandmarkSet"<< std::endl;
    generateLandmarkSet();
    std::cerr << "produceInformationFromLandmarks" << std::endl;
    produceInformationFromLandmarks();
    std::cerr << "alignment" << std::endl;
    Alignment *alignment=new Alignment(this);

    // Exhaustive search for rot
    double bestError=0, bestRot=-1;
    for (double rot=0; rot<=180-deltaRot; rot+=deltaRot)
    {
        alignment->clear();
        alignment->rot=rot;
        double error=alignment->optimizeGivenAxisDirection();
#ifdef DEBUG

        std::cout << "rot= " << rot
        << " error= " << error << std::endl;
#endif

        if (bestRot<0 || bestError>error)
        {
            bestRot=rot;
            bestError=error;
            *bestPreviousAlignment=*alignment;
        }
    }
    delete alignment;
    std::cout << "Best rot=" << bestRot
    << " Best error=" << bestError << std::endl;

    // Continuous optimization for the axis direction
    Matrix1D<double> axisAngles(2), steps(2);
    axisAngles(0)=bestRot;
    axisAngles(1)=90;
    steps.initConstant(1);
    if (!optimizeTiltAngle)
        steps(1)=0;
    double fitness;
    int iter;
    TomographAlignment::global_prm=this;
    powellOptimizer(axisAngles,1,2,&wrapperError, NULL,
                    0.01,fitness,iter,steps,true);

    // Outlier removal
    for (int i=0; i<3; i++)
    {
        // Compute the best alignment
        bestPreviousAlignment->rot=axisAngles(0);
        bestPreviousAlignment->tilt=axisAngles(1);
        fitness=bestPreviousAlignment->optimizeGivenAxisDirection();
        bestPreviousAlignment->computeErrorForLandmarks();

        // Remove landmarks that are outliers in the current model
        removeOutlierLandmarks(*bestPreviousAlignment);
        produceInformationFromLandmarks();
        delete bestPreviousAlignment;
        bestPreviousAlignment=new Alignment(this);
        bestPreviousAlignment->rot=axisAngles(0);
        bestPreviousAlignment->tilt=axisAngles(1);
        fitness=bestPreviousAlignment->optimizeGivenAxisDirection();

        // Optimize again
        powellOptimizer(axisAngles,1,2,&wrapperError,NULL,
                        0.01,fitness,iter,steps,true);
    }
    bestPreviousAlignment->rot=axisAngles(0);
    bestPreviousAlignment->tilt=axisAngles(1);
    fitness=bestPreviousAlignment->optimizeGivenAxisDirection();

    // Save the alignment
    writeLandmarkSet(fnRoot+"_good_landmarks.txt");
    std::ofstream fh_out;
    fh_out.open((fnRoot+"_alignment.txt").c_str());
    if (!fh_out)
        REPORT_ERROR(ERR_IO_NOWRITE,
                     (std::string)"Cannot open "+fnRoot+"_alignment.txt for output");
    fh_out << *bestPreviousAlignment;
    fh_out.close();

    // Correct the input images
    alignImages(*bestPreviousAlignment);
    writeLandmarkSet(fnRoot+"_good_landmarks_corrected.txt");
    delete bestPreviousAlignment;
}
#undef DEBUG

/* Optimize for rot -------------------------------------------------------- */
//#define DEBUG
double Alignment::optimizeGivenAxisDirection()
{
    double bestError;
    bool firstIteration=true, finish=false;
    int Niterations=0;
    computeGeometryDependentOfAxis();
    do
    {
        computeGeometryDependentOfRotation();
        double error=computeError();
#ifdef DEBUG

        std::cout << "it=" << Niterations << " Error= " << error
        << " raxis=" << raxis.transpose() << std::endl;
#endif

        updateModel();
        Niterations++;
        if (firstIteration)
        {
            bestError=error;
            firstIteration=false;
        }
        else
        {
            finish=((error>bestError) || (Niterations>1000) ||
                    (ABS(error-bestError)/bestError<0.001)) && Niterations>20;
            if (error<bestError)
                bestError=error;
        }
    }
    while (!finish);
    return bestError;
}
#undef DEBUG

/* Compute the geometry part corresponding to axis ------------------------- */
void Alignment::computeGeometryDependentOfAxis()
{
    Matrix1D<double> axis;
    Euler_direction(rot, tilt, 0, axis);

    // Compute Aip, Aipt
    Matrix2D<double> Raxis, I, Ip, B;
    Binvraxis.initZeros(3,3);
    I.initIdentity(3);
    Ip=I;
    Ip(2,2)=0;
    for (int i=0; i<Nimg; i++)
    {
        rotation3DMatrix(prm->tiltList[i],axis,Raxis,false);
        Aipt[i](0,0)=Aip[i](0,0)=Raxis(0,0);
        Aipt[i](1,0)=Aip[i](0,1)=Raxis(0,1);
        Aipt[i](2,0)=Aip[i](0,2)=Raxis(0,2);
        Aipt[i](0,1)=Aip[i](1,0)=Raxis(1,0);
        Aipt[i](1,1)=Aip[i](1,1)=Raxis(1,1);
        Aipt[i](2,1)=Aip[i](1,2)=Raxis(1,2);

        B=I-Raxis;
        B=B.transpose()*Ip*B;
        Binvraxis+=((double)prm->ni(i))*(B+B.transpose());

        B1i[i]=I-Raxis.transpose();
        B2i[i]=B1i[i]*Ip*Raxis;
        /*COSS:
                std::cout << "Dependent on axis i=" << i << std::endl
                          << "tilt[i]=" << prm->tiltList[i] << std::endl
                          << "Raxis=\n" << Raxis
                          << "Aip[i]=\n" << Aip[i]
                          << "B1i[i]=\n" << B1i[i] << std::endl
                          << "B2i[i]=\n" << B2i[i] << std::endl;
        */
    }
    Binvraxis=Binvraxis.inv();
}

/* Compute the geometry part corresponding to axis ------------------------- */
void Alignment::computeGeometryDependentOfRotation()
{
    // Compute Ai, Ait
    Matrix2D<double> Rinplane, I(2,3);
    I(0,0)=I(1,1)=1;
    for (int i=0; i<Nimg; i++)
    {
        rotation3DMatrix(-psi(i),'Z',Rinplane,false);
        B1i[i]=B1i[i]*Rinplane.transpose();

        Rinplane.resize(2,2);
        Ai[i]=Rinplane*Aip[i];
        Ait[i]=Ai[i].transpose();
        diaxis[i]=Rinplane*(I-Aip[i])*raxis;
        /* COSS:
                std::cout << "Dependent on rotation i=" << i << std::endl
                          << "B1i[i]=\n" << B1i[i] << std::endl
                          << "I-Raxis=\n" << (I-Aip[i])
                          << "raxis=" << raxis.transpose() << std::endl;
                std::cout << "diaxis[" << i << "]=" << diaxis[i].transpose() << std::endl;
        */
    }
}

/* Compute error ----------------------------------------------------------- */
double Alignment::computeError() const
{
    Matrix1D<double> pijp;
    double error=0;
    double N=0;
    for (int i=0; i<Nimg; i++)
    {
        int jjmax=prm->Vseti[i].size();
        for (int jj=0; jj<jjmax; jj++)
        {
            int j=prm->Vseti[i][jj];
            pijp=Ai[i]*rj[j]+di[i]+diaxis[i];
            MAT_ELEM(allLandmarksPredictedX,j,i)=XX(pijp);
            MAT_ELEM(allLandmarksPredictedY,j,i)=YY(pijp);
            double diffx=MAT_ELEM(prm->allLandmarksX,j,i)-XX(pijp);
            double diffy=MAT_ELEM(prm->allLandmarksY,j,i)-YY(pijp);
            error+=diffx*diffx+diffy*diffy;
            N++;
        }
    }
    return sqrt(error/N);
}

void Alignment::computeErrorForLandmarks()
{
    Nlandmark=MAT_YSIZE(prm->allLandmarksX);
    errorLandmark.initZeros(Nlandmark);
    for (int j=0; j<Nlandmark; j++)
    {
        int counterj=0;
        for (int i=0; i<Nimg; i++)
        {
            if (prm->allLandmarksX(j,i)!=XSIZE(*(prm->img[0])))
            {
                double diffx=MAT_ELEM(prm->allLandmarksX,j,i)-
                             MAT_ELEM(allLandmarksPredictedX,j,i);
                double diffy=MAT_ELEM(prm->allLandmarksY,j,i)-
                             MAT_ELEM(allLandmarksPredictedY,j,i);
                errorLandmark(j)+=sqrt(diffx*diffx+diffy*diffy);
                counterj++;
            }
        }
        errorLandmark(j)/=counterj;
    }
}

/* Update model ------------------------------------------------------------ */
//#define DEBUG
void Alignment::updateModel()
{
    Matrix1D<double> pij(2), piN(2);
    Matrix2D<double> A(3,3), AitAi;
    Matrix1D<double> b(3);

#ifdef DEBUG

    std::cout << "Step 0: error=" << computeError() << std::endl;
#endif

    // Update the 3D positions of the landmarks
    for (int j=0; j<Nlandmark; j++)
    {
        A.initZeros();
        b.initZeros();
        // Compute the part correponding to Vj
        int iimax=prm->Vsetj[j].size();
        for (int ii=0; ii<iimax; ii++)
        {
            int i=prm->Vsetj[j][ii];
            XX(pij)=MAT_ELEM(prm->allLandmarksX,j,i);
            YY(pij)=MAT_ELEM(prm->allLandmarksY,j,i);
            A+=Ait[i]*Ai[i];
            b+=Ait[i]*(pij-(di[i]+diaxis[i]));
        }

        // Update rj[j]
        rj[j]=A.inv()*b;
    }
#ifdef DEBUG
    std::cout << "Step 1: error=" << computeError() << std::endl;
#endif

    // Compute the average landmarks seen in each image
    for (int i=0; i<Nimg; i++)
    {
        barri[i].initZeros();
        for (int jj=0; jj<prm->ni(i); jj++)
        {
            int j=prm->Vseti[i][jj];
            barri[i]+=rj[j];
        }
        barri[i]/=prm->ni(i);
    }

    // Update shifts
    if (prm->isCapillar)
    {
        // Update the raxis
        Matrix1D<double> tmp2(3), Pij(3);
        for (int i=0; i<Nimg; i++)
        {
            for (int jj=0; jj<prm->ni(i); jj++)
            {
                int j=prm->Vseti[i][jj];
                XX(Pij)=MAT_ELEM(prm->allLandmarksX,j,i);
                YY(Pij)=MAT_ELEM(prm->allLandmarksY,j,i);
                tmp2+=B1i[i]*Pij-B2i[i]*rj[j];
            }
        }
        tmp2*=2.0;
        tmp2=Binvraxis*tmp2;
        XX(raxis)=XX(tmp2);
        YY(raxis)=YY(tmp2);
        computeGeometryDependentOfRotation();
    }
    else
    {
        // Update the individual di
        for (int i=0; i<Nimg; i++)
            if (i!=prm->iMinTilt)
            {
                di[i] = prm->barpi[i]-Ai[i]*barri[i]-diaxis[i];

                if (di[i].isAnyNaN())
                    di[i].initZeros();
            }
    }
#ifdef DEBUG
    std::cout << "Step 2: error=" << computeError() << std::endl;
#endif

    // Update rotations
    if (prm->psiMax>0)
    {
        Matrix2D<double> Ri(2,2);
        Matrix2D<double> Aiprj(2,1), Aiprjt(1,2), dim(2,1), Pij(2,1);
        for (int i=0; i<Nimg; i++)
        {
            Ri.initZeros();
            dim.fromVector(di[i]+diaxis[i]);
            for (int jj=0; jj<prm->ni(i); jj++)
            {
                int j=prm->Vseti[i][jj];
                Aiprj.fromVector(Aip[i]*rj[j]);
                Aiprjt=Aiprj.transpose();

                MAT_ELEM(Pij,0,0)=MAT_ELEM(prm->allLandmarksX,j,i);
                MAT_ELEM(Pij,1,0)=MAT_ELEM(prm->allLandmarksY,j,i);
                Ri+=(dim-Pij)*Aiprjt;
            }
            psi(i)=CLIP(RAD2DEG(atan(((Ri(0,1)-Ri(1,0))/(Ri(0,0)+Ri(1,1))))),
                        -(prm->psiMax),prm->psiMax);
            Matrix2D<double> Rinplane;
            rotation3DMatrix(-psi(i),'Z',Rinplane,false);
        }
    }
#ifdef DEBUG
    std::cout << "Step 3: error=" << computeError() << std::endl;
#endif

    // Update the rotation dependent part
    computeGeometryDependentOfRotation();
}
#undef DEBUG

/* Print ------------------------------------------------------------------- */
std::ostream& operator << (std::ostream &out, Alignment &alignment)
{
    out << "Alignment parameters ===========================================\n"
    << "rot=" << alignment.rot << " tilt=" << alignment.tilt << std::endl
    << "raxis=" << alignment.raxis.transpose() << std::endl;
    out << "Images ---------------------------------------------------------\n";
    for (int i=0; i<alignment.Nimg; i++)
        out << "Image " << i << " psi= " << alignment.psi(i)
        << " di= " << alignment.di[i].transpose()
        << " diaxis= " << alignment.diaxis[i].transpose()
        << " (" << alignment.prm->ni(i) << ")\n";
    out << "Landmarks ------------------------------------------------------\n";
    alignment.computeErrorForLandmarks();
    for (int j=0; j<alignment.Nlandmark; j++)
        out << "Landmark " << j << " rj= " << alignment.rj[j].transpose()
        << " " << alignment.errorLandmark(j) << std::endl;
    return out;
}
