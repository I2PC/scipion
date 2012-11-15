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

void ProgSortByStatistics::readParams()
{
    fn = getParam("-i");
    fn_out = getParam("-o");
    addToInput = checkParam("--addToInput");
    fn_train = getParam("--train");
    cutoff = getDoubleParam("--zcut");
}

void ProgSortByStatistics::defineParams()
{
    addUsageLine("Sorts the input images for identifying junk particles");
    addUsageLine("+The program associates to each image a vector composed by");
    addUsageLine("+the histogram of the image (this accounts for factors such as");
    addUsageLine("+min, max, avg, and standard deviation, plus a more complete");
    addUsageLine("+description of the image gray levels) and the radial profile of");
    addUsageLine("+the image squared.");
    addUsageLine("+");
    addUsageLine("+These vectors are then scored according to a Gaussian distribution");
    //addUsageLine("+that can be chosen to be univariate or multivariate. The multivariate");
    addUsageLine("+gaussian is more powerful in the sense that it can capture relationships");
    addUsageLine("+among variables.");
    addUsageLine("+");
    addUsageLine("+If you choose a threshold, you must take into account that it is a zscore.");
    addUsageLine("+For univariate and mulivariate Gaussian distributions, 99% of the individuals");
    addUsageLine("+have a Z-score below 3");
    addParamsLine(" -i <selfile>            : Selfile with input images");
    addParamsLine(" [-o <rootname=\"\">]    : Output rootname");
    addParamsLine("                         : rootname.xmd contains the list of sorted images with their Z-score");
    addParamsLine("                         : rootname_vectors.xmd (if verbose>=2) contains the vectors associated to each image");
    addParamsLine("                         : If no rootname is given, these two files are not created");
    addParamsLine(" [--train <selfile=\"\">]: Train on selfile with good particles");
    addParamsLine(" [--zcut <float=-1>]     : Cut-off for Z-scores (negative for no cut-off) ");
    addParamsLine("                         : Images whose Z-score is larger than the cutoff are disabled");
    addParamsLine(" [--addToInput]          : Add columns also to input MetaData");
}


//majorAxis and minorAxis is the estimated particle size in px
void ProgSortByStatistics::processInprocessInputPrepareSPTH(MetaData &SF)
{
#define DEBUG

    pcaAnalyzer[2];
    PCAMahalanobisAnalyzer tempPcaAnalyzer0;
    PCAMahalanobisAnalyzer tempPcaAnalyzer1;
    PCAMahalanobisAnalyzer tempPcaAnalyzer2;

    //Morphology
    tempPcaAnalyzer0.clear();
    //Signal to noise ratio
    tempPcaAnalyzer1.clear();
    //Histogram analysis, to detect black points and saturated parts
    tempPcaAnalyzer2.clear();

    Image<double> img;
    MultidimArray<double> img2;
    Matrix1D<double> center(2);
    center.initZeros();
    FringeProcessing fp;

    int numHistElem = 21;
    int numDescriptors=4;
    int numDescriptors2=3;

    MultidimArray<float> v(numDescriptors);
    MultidimArray<float> v2(numHistElem);
    MultidimArray<float> v3(numDescriptors2);

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
    img.readApplyGeo(SF,1);

    MultidimArray<double> nI, modI, tempI, tempM, Imean;
    MultidimArray<bool> mask;
    nI.resizeNoCopy(img());
    modI.resizeNoCopy(img());
    Imean.resizeNoCopy(img());
    tempI.resizeNoCopy(img());
    tempM.resizeNoCopy(img());
    mask.resizeNoCopy(img());

    Imean.initConstant(0);
    mask.initConstant(true);

    MultidimArray<double> autoCorr(2*img().ydim,2*img().xdim);
    MultidimArray<double> smallAutoCorr;

    Histogram1D hist;
    Matrix2D<double> U,V,temp;
    Matrix1D<double> D;
    v.initZeros(numDescriptors);
    v2.initZeros(numHistElem);
    v3.initZeros(numDescriptors2);

    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        if (thereIsEnable)
        {
            int enabled;
            SF.getValue(MDL_ENABLED,enabled,__iter.objId);
            if ( (enabled==-1)  )
            {
                continue;
            }
            img.readApplyGeo(SF,__iter.objId);
            //
            img().statisticsAdjust(0,1);
            Imean += img();
        }
    }

    Imean /= SF.size();
    Imean.setXmippOrigin();

#ifdef DEBUG

    FileName fpName    = "test0.txt";
    Imean.write(fpName);
#endif

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
            bool debug=false;

            MultidimArray<double> &mI=img();
            mI.setXmippOrigin();
            mI.statisticsAdjust(0,1);
            mask.setXmippOrigin();

            fp.normalize(mI,tempI,modI,1,5,mask);
            modI.setXmippOrigin();
            tempI.setXmippOrigin();
            nI = tempI*(modI*modI);
            tempM = (modI*modI);

            //fp.normalize(mI,tempI,modI,1,4,mask);
            //tempI.setXmippOrigin();
            //modI.setXmippOrigin();
            //nI += (tempI*(modI*modI));
            //tempM += (modI*modI);

            //fp.normalize(mI,tempI,modI,1,6,mask);
            //tempI.setXmippOrigin();
            //modI.setXmippOrigin();
            //nI += (tempI*(modI*modI));
            //tempM += (modI*modI);

            nI /= tempM;

#ifdef DEBUG
            //fn = img.name().removeLastExtension().addExtension("spi");
            //nI.write(fn);
            //SF.setValue(MDL_IMAGE,fn,__iter.objId);
            if (img.name()=="002870@Images/Extracted/run_001/DefaultFamily5")
            {
                FileName fpName    = "test2.txt";
                nI.write(fpName);
                fpName    = "test4.txt";
                tempM.write(fpName);
            }
#endif
            nI.binarize();
            int im = labelImage2D(nI,nI,8);
            compute_hist(nI, hist, 0, im, im+1);
            int l,k,i,j;
            hist.maxIndex(l,k,i,j);
            A1D_ELEM(hist,j)=0;
            hist.maxIndex(l,k,i,j);
            nI.binarizeRange(j-1,j+1);
            FringeProcessing fp;
            double x0=0,y0=0,majorAxis=0,minorAxis=0,ellipAng=0,area=0;
            fp.fitEllipse(nI,x0,y0,majorAxis,minorAxis,ellipAng,area);
            // Build vector
            v(0)=majorAxis;
            v(1)=minorAxis;
            v(2)= fabs((img().xdim)/2-x0)+fabs((img().ydim)/2-y0) ;
            v(3)=area;

            for (int n=0 ; n < numDescriptors ; n++)
            {
                if ( std::isnan(v(n)) )
                    v(n)=0;
            }
            tempPcaAnalyzer0.addVector(v);


            //tempI = (Imean*mI)/mI.sum();
            //double val = tempI.sum();

            mI-=Imean;
            mI.setXmippOrigin();
            auto_correlation_matrix(mI,autoCorr);
            autoCorr.window(smallAutoCorr,-10,-10, 10, 10);
            smallAutoCorr.copy(temp);
            svdcmp(temp,U,D,V);

            for (int n = 0; n < numHistElem; ++n)
            	v2(n)=(float)VEC_ELEM(D,n);
            tempPcaAnalyzer1.addVector(v2);


            double minVal;
            double maxVal;
            mI.computeDoubleMinMax(minVal,maxVal);
            compute_hist(mI, hist, minVal, maxVal, 100);

            v3(0)= fabs(hist.percentil(5));
            v3(1)= fabs(hist.percentil(95));
            v3(2) = fabs(hist.percentil(50));
            tempPcaAnalyzer2.addVector(v3);

#ifdef DEBUG

            if (img.name()=="002870@Images/Extracted/run_001/DefaultFamily5")
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


    MultidimArray<double> vavg,vstddev;
    tempPcaAnalyzer0.computeStatistics(vavg,vstddev);
    tempPcaAnalyzer1.computeStatistics(vavg,vstddev);
    tempPcaAnalyzer2.computeStatistics(vavg,vstddev);

    tempPcaAnalyzer0.evaluateZScore(2,20);
    tempPcaAnalyzer1.evaluateZScore(2,20);
    tempPcaAnalyzer2.evaluateZScore(2,20);

    pcaAnalyzer.insert(pcaAnalyzer.begin(), tempPcaAnalyzer2);
    pcaAnalyzer.insert(pcaAnalyzer.begin(), tempPcaAnalyzer1);
    pcaAnalyzer.insert(pcaAnalyzer.begin(), tempPcaAnalyzer0);

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
        for (int i=0; i<XSIZE(img2)/2; i++)
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

    if (fn_train != "")
        SFtrain.read(fn_train);

    else
        processInprocessInputPrepareSPTH(SF);
    //processInputPrepare(SF);

    std::cout << " OUT " << std::endl;
    int imgno = 0;
    int numPCAs = pcaAnalyzer.size();

    MultidimArray<double> ZscoreMultivariate(SF.size());
    MultidimArray<double> weights(numPCAs);
    ZscoreMultivariate.initConstant(0);
    weights.initConstant(1);
    weights(0) = 0.3;
    weights(1) = 0.5;
    weights(2) = 0.2;

    int enabled;

    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        for (int num = 0; num < numPCAs; ++num)
        {
            SF.getValue(MDL_ENABLED,enabled,__iter.objId);
            ZscoreMultivariate(imgno)+=(pcaAnalyzer[num].getZscore(imgno))*weights(num)*(pcaAnalyzer[num].getZscore(imgno));

        }
        imgno++;
    }

    // Produce output .....................................................
    MetaData SFout;
    std::ofstream fh_zind;

    if (verbose==2 && !fn_out.empty())
        fh_zind.open((fn_out.withoutExtension() + "_vectors.xmd").c_str(), std::ios::out);
    MultidimArray<double> finalZscore=ZscoreMultivariate;
    double L=1;

    finalZscore*=ZscoreMultivariate;
    L++;

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(finalZscore)
    DIRECT_A1D_ELEM(finalZscore,i)=pow(DIRECT_A1D_ELEM(finalZscore,i),1.0/L);

    MultidimArray<int> sorted;
    finalZscore.indexSort(sorted);
    int nr_imgs = SF.size();
    bool thereIsEnable=SF.containsLabel(MDL_ENABLED);
    FileName fnImg;
    MDRow row;

    for (int imgno = 0; imgno < nr_imgs; imgno++)
    {
        int isort_1 = DIRECT_A1D_ELEM(sorted,imgno);
        int isort = isort_1 - 1;
        //SF.getValue(MDL_IMAGE,fnImg,isort+1);
        SF.getRow(row, isort_1);
        row.getValue(MDL_IMAGE, fnImg);

        if (thereIsEnable)
        {
            int enabled;
            row.getValue(MDL_ENABLED, enabled);
            if (enabled==-1)
                continue;
        }
        //size_t objId = SFout.addObject();
        row.setValue(MDL_IMAGE,fnImg);
        double zscore=DIRECT_A1D_ELEM(finalZscore,isort);
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
        if (addToInput)
            SF.setValue(MDL_ZSCORE,zscore,isort_1);
        if (verbose==2)
        {
            fh_zind << fnImg << " : ";
            //FOR_ALL_ELEMENTS_IN_ARRAY1D(tempPcaAnalyzer.v[isort])
            //fh_zind << tempPcaAnalyzer.v[isort](i) << "\n";
        }
        SFout.addRow(row);
    }
    if (verbose==2)
        fh_zind.close();
    if (!fn_out.empty())
        SFout.write(fn_out,MD_OVERWRITE);
    if (addToInput)
        SF.write(fn,MD_APPEND);

}
