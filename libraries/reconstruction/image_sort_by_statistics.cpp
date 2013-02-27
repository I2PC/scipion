/***************************************************************************
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
 *    Javier Vargas (jvargas@cnb.csic.es)
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

#include "image_sort_by_statistics.h"

void ProgSortByStatistics::clear()
{
    pcaAnalyzer.clear();
}

void ProgSortByStatistics::readParams()
{
    fn = getParam("-i");
    fn_out = getParam("-o");
    addToInput = checkParam("--addToInput");
    fn_train = getParam("--train");
    cutoff = getDoubleParam("--zcut");
    per = getDoubleParam("--percent");
}

void ProgSortByStatistics::defineParams()
{
    addUsageLine("Sorts the input images for identifying junk particles");
    addUsageLine("+The program associates to each image four vectors.");
    addUsageLine("+One vector is composed by descriptors that encodes the particle shape.");
    addUsageLine("+Another two vectors give information about the SNR of the objects.");
    addUsageLine("+Finally, the last vector provides information of the particle histogram");
    addUsageLine("+");
    addUsageLine("+These vector are then scored according to a multivariate Gaussian distribution");
    addUsageLine("+");
    addUsageLine("+You can reject erroneous particles choosing a threshold, you must take into account");
    addUsageLine("+that it is a zscore.");
    addUsageLine("+For univariate and mulivariate Gaussian distributions, 99% of the individuals");
    addUsageLine("+have a Z-score below 3.");
    addUsageLine("+Additionally, you can discard bad particles selecting an inaccurate particle percentage");
    addUsageLine("+typically around 10-20%.");
    addParamsLine(" -i <selfile>            : Selfile with input images");
    addParamsLine(" [-o <rootname=\"\">]    : Output rootname");
    addParamsLine("                         : rootname.xmd contains the list of sorted images with their Z-score");
    addParamsLine("                         : rootname_vectors.xmd (if verbose>=2) contains the vectors associated to each image");
    addParamsLine("                         : If no rootname is given, these two files are not created");
    addParamsLine(" [--train <selfile=\"\">]: Train on selfile with good particles");
    addParamsLine(" [--zcut <float=-1>]     : Cut-off for Z-scores (negative for no cut-off) ");
    addParamsLine("                         : Images whose Z-score is larger than the cutoff are disabled");
    addParamsLine(" [--percent <float=0>]   : Cut-off for particles (zero for no cut-off) ");
    addParamsLine("                         : Percentage of images with larger Z-scores are disabled");
    addParamsLine(" [--addToInput]          : Add columns also to input MetaData");
}


//majorAxis and minorAxis is the estimated particle size in px
void ProgSortByStatistics::processInprocessInputPrepareSPTH(MetaData &SF)
{
    //#define DEBUG
    pcaAnalyzer[5];
    PCAMahalanobisAnalyzer tempPcaAnalyzer0;
    PCAMahalanobisAnalyzer tempPcaAnalyzer1;
    PCAMahalanobisAnalyzer tempPcaAnalyzer2;
    PCAMahalanobisAnalyzer tempPcaAnalyzer3;
    PCAMahalanobisAnalyzer tempPcaAnalyzer4;

    //Morphology
    tempPcaAnalyzer0.clear();
    //Signal to noise ratio
    tempPcaAnalyzer1.clear();
    tempPcaAnalyzer2.clear();
    tempPcaAnalyzer3.clear();
    //Histogram analysis, to detect black points and saturated parts
    tempPcaAnalyzer4.clear();

    double sign = 1;//;-1;
    int numNorm = 3;
    int numDescriptors0=numNorm;
    int numDescriptors2=4;
    int numDescriptors3=11;
    int numDescriptors4 = 10;

    MultidimArray<float> v0(numDescriptors0);
    MultidimArray<float> v2(numDescriptors2);
    MultidimArray<float> v3(numDescriptors3);
    MultidimArray<float> v4(numDescriptors4);

    if (verbose>0)
    {
        std::cout << " Sorting particle set by new xmipp method..." << std::endl;
    }

    int nr_imgs = SF.size();
    if (verbose>0)
        init_progress_bar(nr_imgs);

    int c = XMIPP_MAX(1, nr_imgs / 60);
    int imgno = 0, imgnoPCA=0;

    bool thereIsEnable=SF.containsLabel(MDL_ENABLED);
    bool first=true;

    // We assume that at least there is one particle
    size_t Xdim, Ydim, Zdim, Ndim;
    getImageSize(SF,Xdim,Ydim,Zdim,Ndim);

    //Initialization:
    MultidimArray<double> nI, modI, tempI, tempM, ROI;
    MultidimArray<bool> mask;
    nI.resizeNoCopy(Ydim,Xdim);
    modI.resizeNoCopy(Ydim,Xdim);
    tempI.resizeNoCopy(Ydim,Xdim);
    tempM.resizeNoCopy(Ydim,Xdim);
    mask.resizeNoCopy(Ydim,Xdim);
    mask.initConstant(true);

    MultidimArray<double> autoCorr(2*Ydim,2*Xdim);
    MultidimArray<double> smallAutoCorr;

    Histogram1D hist;
    Matrix2D<double> U,V,temp;
    Matrix1D<double> D;

    MultidimArray<int> radial_count;
    MultidimArray<double> radial_avg;
    Matrix1D<int> center(2);
    MultidimArray<int> distance;
    int dim;
    center.initZeros();

    v0.initZeros(numDescriptors0);
    v2.initZeros(numDescriptors2);
    v3.initZeros(numDescriptors3);
    v4.initZeros(numDescriptors4);

    ROI.resizeNoCopy(Ydim,Xdim);
    ROI.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY2D(ROI)
    {
        double temp = std::sqrt(i*i+j*j);
        if ( temp < (Xdim/2))
            A2D_ELEM(ROI,i,j)= 1;
        else
            A2D_ELEM(ROI,i,j)= 0;
    }

    Image<double> img;
    FourierTransformer transformer(FFTW_BACKWARD);

    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        if (thereIsEnable)
        {
            int enabled;
            SF.getValue(MDL_ENABLED,enabled,__iter.objId);
            if ( (enabled==-1)  )
            {
                imgno++;
                continue;
            }

            img.readApplyGeo(SF,__iter.objId);

            MultidimArray<double> &mI=img();
            mI.setXmippOrigin();
            mI.statisticsAdjust(0,1);
            mask.setXmippOrigin();
            //The size of v1 depends on the image size and must be declared here
            int numDescriptors1 = XSIZE(mI)/2; //=100;
            MultidimArray<float> v1(numDescriptors1);
            v1.initZeros(numDescriptors1);

            double var = 1;
            normalize(transformer,mI,tempI,modI,0,var,mask);
            modI.setXmippOrigin();
            tempI.setXmippOrigin();
            nI = sign*tempI*(modI*modI);
            tempM = (modI*modI);

            A1D_ELEM(v0,0) = (tempM*ROI).sum();
            int index = 1;
            var+=2;
            while (index < numNorm)
            {
                normalize(transformer,mI,tempI,modI,0,var,mask);
                modI.setXmippOrigin();
                tempI.setXmippOrigin();
                nI += sign*tempI*(modI*modI);
                tempM += (modI*modI);
                A1D_ELEM(v0,index) = (tempM*ROI).sum();
                index++;
                var+=2;
            }

            nI /= tempM;
            tempPcaAnalyzer0.addVector(v0);
            nI=(nI*ROI);

            auto_correlation_matrix(mI,autoCorr);
            if (first)
            {
                radialAveragePrecomputeDistance(autoCorr, center, distance, dim);
                first=false;
            }
            fastRadialAverage(autoCorr, distance, dim, radial_avg, radial_count);

            for (int n = 0; n < numDescriptors1; ++n)
            	A1D_ELEM(v1,n)=(float)DIRECT_A1D_ELEM(radial_avg,n);

            tempPcaAnalyzer1.addVector(v1);

#ifdef DEBUG

            //String name = "000005@Images/Extracted/run_002/extra/BPV_1386.stk";
            String name = "000010@Images/Extracted/run_001/extra/KLH_Dataset_I_Training_0028.stk";
            //String name = "001160@Images/Extracted/run_001/DefaultFamily5";

            std::cout << img.name() << std::endl;

            if (img.name()==name2)
            {
                FileName fpName    = "test_1.txt";
                mI.write(fpName);
                fpName    = "test_2.txt";
                nI.write(fpName);
                fpName    = "test_3.txt";
                tempM.write(fpName);
                fpName    = "test_4.txt";
                ROI.write(fpName);
                //exit(1);
            }
#endif
            nI.binarize(0);
            int im = labelImage2D(nI,nI,8);
            compute_hist(nI, hist, 0, im, im+1);
            int l,k,i,j;
            hist.maxIndex(l,k,i,j);
            A1D_ELEM(hist,j)=0;
            hist.maxIndex(l,k,i,j);
            nI.binarizeRange(j-1,j+1);

            double x0=0,y0=0,majorAxis=0,minorAxis=0,ellipAng=0;
            size_t area=0;
            fitEllipse(nI,x0,y0,majorAxis,minorAxis,ellipAng,area);

            A1D_ELEM(v2,0)=majorAxis/((img().xdim) );
            A1D_ELEM(v2,1)=minorAxis/((img().xdim) );
            A1D_ELEM(v2,2)= (fabs((img().xdim)/2-x0)+fabs((img().ydim)/2-y0))/((img().xdim)/2);
            A1D_ELEM(v2,3)=area/( ((img().xdim)/2)*((img().ydim)/2) );

            for (int n=0 ; n < numDescriptors2 ; n++)
            {
                if ( std::isnan(A1D_ELEM(v2,n)) )
                    A1D_ELEM(v2,n)=0;
            }
            tempPcaAnalyzer2.addVector(v2);

            //mI.setXmippOrigin();
            //auto_correlation_matrix(mI*ROI,autoCorr);
            //auto_correlation_matrix(nI,autoCorr);
            autoCorr.window(smallAutoCorr,-5,-5, 5, 5);
            smallAutoCorr.copy(temp);
            svdcmp(temp,U,D,V);

            for (int n = 0; n < numDescriptors3; ++n)
                A1D_ELEM(v3,n)=(float)VEC_ELEM(D,n); //A1D_ELEM(v3,n)=(float)VEC_ELEM(D,n)/VEC_ELEM(D,0);

            tempPcaAnalyzer3.addVector(v3);


            double minVal;
            double maxVal;
            mI.computeDoubleMinMax(minVal,maxVal);
            compute_hist(mI, hist, minVal, maxVal, 100);

            for (int n=0 ; n <= numDescriptors4-1 ; n++)
            {
                A1D_ELEM(v4,n)= (hist.percentil((n+1)*10));
            }
            tempPcaAnalyzer4.addVector(v4);

#ifdef DEBUG

            if (img.name()==name1)
            {
                FileName fpName    = "test.txt";
                mI.write(fpName);
                fpName    = "test3.txt";
                nI.write(fpName);
            }
#endif
            imgno++;
            imgnoPCA++;

            if (imgno % c == 0 && verbose>0)
                progress_bar(imgno);
        }
    }

    tempPcaAnalyzer0.evaluateZScore(2,20);
    tempPcaAnalyzer1.evaluateZScore(2,20);
    tempPcaAnalyzer2.evaluateZScore(2,20);
    tempPcaAnalyzer3.evaluateZScore(2,20);
    tempPcaAnalyzer4.evaluateZScore(2,20);

    pcaAnalyzer.push_back(tempPcaAnalyzer0);
    pcaAnalyzer.push_back(tempPcaAnalyzer1);
    pcaAnalyzer.push_back(tempPcaAnalyzer2);
    pcaAnalyzer.push_back(tempPcaAnalyzer3);
    pcaAnalyzer.push_back(tempPcaAnalyzer4);

}

void ProgSortByStatistics::processInputPrepare(MetaData &SF)
{

    pcaAnalyzer[1];
    PCAMahalanobisAnalyzer tempPcaAnalyzer;
    tempPcaAnalyzer.clear();

    Image<double> img;
    MultidimArray<double> img2;
    MultidimArray<int> radial_count;
    MultidimArray<double> radial_avg;
    Matrix1D<int> center(2);
    center.initZeros();

    if (verbose>0)
        std::cout << " Processing training set ..." << std::endl;

    int nr_imgs = SF.size();
    if (verbose>0)
        init_progress_bar(nr_imgs);
    int c = XMIPP_MAX(1, nr_imgs / 60);
    int imgno = 0, imgnoPCA=0;
    MultidimArray<float> v;
    MultidimArray<int> distance;
    int dim;

    bool thereIsEnable=SF.containsLabel(MDL_ENABLED);
    bool first=true;
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        if (thereIsEnable)
        {
            int enabled;
            SF.getValue(MDL_ENABLED,enabled,__iter.objId);
            if (enabled==-1)
                continue;
        }
        img.readApplyGeo(SF,__iter.objId);
        MultidimArray<double> &mI=img();
        mI.setXmippOrigin();
        mI.statisticsAdjust(0,1);

        // Overall statistics
        Histogram1D hist;
        compute_hist(mI,hist,-4,4,31);

        // Radial profile
        img2.resizeNoCopy(mI);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img2)
        {
            double val=DIRECT_MULTIDIM_ELEM(mI,n);
            DIRECT_MULTIDIM_ELEM(img2,n)=val*val;
        }
        if (first)
        {
            radialAveragePrecomputeDistance(img2, center, distance, dim);
            first=false;
        }
        fastRadialAverage(img2, distance, dim, radial_avg, radial_count);

        // Build vector
        v.initZeros(XSIZE(hist)+XSIZE(img2)/2);
        int idx=0;
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(hist)
        v(idx++)=(float)DIRECT_A1D_ELEM(hist,i);
        for (size_t i=0; i<XSIZE(img2)/2; i++)
            v(idx++)=(float)DIRECT_A1D_ELEM(radial_avg,i);

        tempPcaAnalyzer.addVector(v);

        if (imgno % c == 0 && verbose>0)
            progress_bar(imgno);
        imgno++;
        imgnoPCA++;
    }
    if (verbose>0)
        progress_bar(nr_imgs);

    MultidimArray<double> vavg,vstddev;
    tempPcaAnalyzer.computeStatistics(vavg,vstddev);
    tempPcaAnalyzer.evaluateZScore(2,20);
    pcaAnalyzer.insert(pcaAnalyzer.begin(), tempPcaAnalyzer);
}

void ProgSortByStatistics::run()
{
    // Process input selfile ..............................................
    SF.read(fn);
    SF.removeDisabled();
    MetaData SF2 = SF;
    SF = SF2;

    if (fn_train != "")
        SFtrain.read(fn_train);
    else
        processInprocessInputPrepareSPTH(SF);

    int imgno = 0;
    int numPCAs = pcaAnalyzer.size();

    MultidimArray<double> finalZscore(SF.size());
    MultidimArray<double> ZscoreShape1(SF.size()), sortedZscoreShape1;
    MultidimArray<double> ZscoreShape2(SF.size()), sortedZscoreShape2;
    MultidimArray<double> ZscoreSNR1(SF.size()), sortedZscoreSNR1;
    MultidimArray<double> ZscoreSNR2(SF.size()), sortedZscoreSNR2;
    MultidimArray<double> ZscoreHist(SF.size()), sortedZscoreHist;

    finalZscore.initConstant(0);
    ZscoreShape1.resizeNoCopy(finalZscore);
    ZscoreShape2.resizeNoCopy(finalZscore);
    ZscoreSNR1.resizeNoCopy(finalZscore);
    ZscoreSNR2.resizeNoCopy(finalZscore);
    ZscoreHist.resizeNoCopy(finalZscore);
    sortedZscoreShape1.resizeNoCopy(finalZscore);
    sortedZscoreShape2.resizeNoCopy(finalZscore);
    sortedZscoreSNR1.resizeNoCopy(finalZscore);
    sortedZscoreSNR2.resizeNoCopy(finalZscore);
    sortedZscoreHist.resizeNoCopy(finalZscore);

    double zScore=0;
    int enabled;

    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        SF.getValue(MDL_ENABLED,enabled,__iter.objId);
        if ( (enabled==-1)  )
        {
            A1D_ELEM(finalZscore,imgno) = 1e3;
            A1D_ELEM(ZscoreShape1,imgno) = 1e3;
            A1D_ELEM(ZscoreShape2,imgno) = 1e3;
            A1D_ELEM(ZscoreSNR1,imgno) = 1e3;
            A1D_ELEM(ZscoreSNR2,imgno) = 1e3;
            A1D_ELEM(ZscoreHist,imgno) = 1e3;
            imgno++;
            enabled = 0;
        }
        else
        {
            for (int num = 0; num < numPCAs; ++num)
            {
                if (num == 0)
                {
                    A1D_ELEM(ZscoreSNR1,imgno) = pcaAnalyzer[num].getZscore(imgno);
                }
                else if (num == 1)
                {
                    A1D_ELEM(ZscoreShape2,imgno) = pcaAnalyzer[num].getZscore(imgno);
                }
                else if (num == 2)
                {
                    A1D_ELEM(ZscoreShape1,imgno) = pcaAnalyzer[num].getZscore(imgno);
                }
                else if (num == 3)
                {
                    A1D_ELEM(ZscoreSNR2,imgno) = pcaAnalyzer[num].getZscore(imgno);
                }
                else
                {
                	A1D_ELEM(ZscoreHist,imgno) = pcaAnalyzer[num].getZscore(imgno);
                }

                if(zScore < pcaAnalyzer[num].getZscore(imgno))
                    zScore = pcaAnalyzer[num].getZscore(imgno);
            }

            A1D_ELEM(finalZscore,imgno)=zScore;
            imgno++;
            zScore = 0;
        }
    }

    // Produce output .....................................................
    MetaData SFout;
    std::ofstream fh_zind;

    if (verbose==2 && !fn_out.empty())
        fh_zind.open((fn_out.withoutExtension() + "_vectors.xmd").c_str(), std::ios::out);

    MultidimArray<int> sorted;
    finalZscore.indexSort(sorted);

    int nr_imgs = SF.size();
    bool thereIsEnable=SF.containsLabel(MDL_ENABLED);
    MDRow row;

    for (int imgno = 0; imgno < nr_imgs; imgno++)
    {
        int isort_1 = DIRECT_A1D_ELEM(sorted,imgno);
        int isort = isort_1 - 1;
        SF.getRow(row, isort_1);

        if (thereIsEnable)
        {
            int enabled;
            row.getValue(MDL_ENABLED, enabled);
            if (enabled==-1)
                continue;
        }

        double zscore=DIRECT_A1D_ELEM(finalZscore,isort);
        double zscoreShape1=DIRECT_A1D_ELEM(ZscoreShape1,isort);
        double zscoreShape2=DIRECT_A1D_ELEM(ZscoreShape2,isort);
        double zscoreSNR1=DIRECT_A1D_ELEM(ZscoreSNR1,isort);
        double zscoreSNR2=DIRECT_A1D_ELEM(ZscoreSNR2,isort);
        double zscoreHist=DIRECT_A1D_ELEM(ZscoreHist,isort);

        DIRECT_A1D_ELEM(sortedZscoreShape1,imgno)=DIRECT_A1D_ELEM(ZscoreShape1,isort);
        DIRECT_A1D_ELEM(sortedZscoreShape2,imgno)=DIRECT_A1D_ELEM(ZscoreShape2,isort);
        DIRECT_A1D_ELEM(sortedZscoreSNR1,imgno)=DIRECT_A1D_ELEM(ZscoreSNR1,isort);
        DIRECT_A1D_ELEM(sortedZscoreSNR2,imgno)=DIRECT_A1D_ELEM(ZscoreSNR2,isort);
        DIRECT_A1D_ELEM(sortedZscoreHist,imgno)=DIRECT_A1D_ELEM(ZscoreHist,isort);

        if (zscore>cutoff && cutoff>0)
        {
            row.setValue(MDL_ENABLED,-1);
            if (addToInput)
                SF.setValue(MDL_ENABLED,-1,isort_1);
        }
        else
        {
            row.setValue(MDL_ENABLED,1);
            if (addToInput)
                SF.setValue(MDL_ENABLED,1,isort_1);
        }

        row.setValue(MDL_ZSCORE,zscore);
        row.setValue(MDL_ZSCORE_SHAPE1,zscoreShape1);
        row.setValue(MDL_ZSCORE_SHAPE2,zscoreShape2);
        row.setValue(MDL_ZSCORE_SNR1,zscoreSNR1);
        row.setValue(MDL_ZSCORE_SNR2,zscoreSNR2);
        row.setValue(MDL_ZSCORE_HISTOGRAM,zscoreHist);

        if (addToInput)
        {
            SF.setValue(MDL_ZSCORE,zscore,isort_1);
            SF.setValue(MDL_ZSCORE_SHAPE1,zscoreShape1,isort_1);
            SF.setValue(MDL_ZSCORE_SHAPE2,zscoreShape2,isort_1);
            SF.setValue(MDL_ZSCORE_SNR1,zscoreSNR1,isort_1);
            SF.setValue(MDL_ZSCORE_SNR2,zscoreSNR2,isort_1);
            SF.setValue(MDL_ZSCORE_HISTOGRAM,zscoreHist,isort_1);
        }

        SFout.addRow(row);
    }

    //Sorting taking into account a given percentage
    if (per > 0)
    {
        MultidimArray<int> sortedShape1,sortedShape2,sortedSNR1,sortedSNR2,sortedHist,
        					sortedShapeSF1,sortedShapeSF2,sortedSNR1SF,sortedSNR2SF,sortedHistSF;

        sortedZscoreShape1.indexSort(sortedShape1);
        sortedZscoreShape2.indexSort(sortedShape2);
        sortedZscoreSNR1.indexSort(sortedSNR1);
        sortedZscoreSNR2.indexSort(sortedSNR2);
        sortedZscoreHist.indexSort(sortedHist);
        size_t numPartReject = (size_t)std::floor((per/100)*SF.size());

        for (size_t numPar = SF.size()-1; numPar > (SF.size()-numPartReject); --numPar)
        {
            int isort_1 = DIRECT_A1D_ELEM(sortedShape1,numPar);
            SFout.getRow(row, isort_1);
            row.setValue(MDL_ENABLED,-1);
            SFout.setRow(row,isort_1);

            isort_1 = DIRECT_A1D_ELEM(sortedShape2,numPar);
            SFout.getRow(row, isort_1);
            row.setValue(MDL_ENABLED,-1);
            SFout.setRow(row,isort_1);

            isort_1 = DIRECT_A1D_ELEM(sortedSNR1,numPar);
            SFout.getRow(row, isort_1);
            row.setValue(MDL_ENABLED,-1);
            SFout.setRow(row,isort_1);

            isort_1 = DIRECT_A1D_ELEM(sortedSNR2,numPar);
            SFout.getRow(row, isort_1);
            row.setValue(MDL_ENABLED,-1);
            SFout.setRow(row,isort_1);

            isort_1 = DIRECT_A1D_ELEM(sortedHist,numPar);
            SFout.getRow(row, isort_1);
            row.setValue(MDL_ENABLED,-1);
            SFout.setRow(row,isort_1);

            if (addToInput)
            {
                ZscoreShape1.indexSort(sortedShapeSF1);
                ZscoreShape2.indexSort(sortedShapeSF2);
                ZscoreSNR1.indexSort(sortedSNR1SF);
                ZscoreSNR2.indexSort(sortedSNR2SF);
                ZscoreHist.indexSort(sortedHistSF);

                isort_1 = DIRECT_A1D_ELEM(sortedShapeSF1,numPar);
                SF.getRow(row, isort_1);
                row.setValue(MDL_ENABLED,-1);
                SF.setRow(row,isort_1);

                isort_1 = DIRECT_A1D_ELEM(sortedShapeSF2,numPar);
                SF.getRow(row, isort_1);
                row.setValue(MDL_ENABLED,-1);
                SF.setRow(row,isort_1);

                isort_1 = DIRECT_A1D_ELEM(sortedSNR1SF,numPar);
                SF.getRow(row, isort_1);
                row.setValue(MDL_ENABLED,-1);
                SF.setRow(row,isort_1);

                isort_1 = DIRECT_A1D_ELEM(sortedSNR2SF,numPar);
                SF.getRow(row, isort_1);
                row.setValue(MDL_ENABLED,-1);
                SF.setRow(row,isort_1);

                isort_1 = DIRECT_A1D_ELEM(sortedHistSF,numPar);
                SF.getRow(row, isort_1);
                row.setValue(MDL_ENABLED,-1);
                SF.setRow(row,isort_1);
            }
        }
    }

    if (verbose==2)
        fh_zind.close();
    if (!fn_out.empty())
    {
    	MetaData SFsorted;
    	SFsorted.sort(SFout,MDL_ZSCORE);
        SFout.write(fn_out,MD_OVERWRITE);
    }
    if (addToInput)
    {
    	MetaData SFsorted;
    	SFsorted.sort(SF,MDL_ZSCORE);
        SFsorted.write(fn,MD_APPEND);
    }
}
