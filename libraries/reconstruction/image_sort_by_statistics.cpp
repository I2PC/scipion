/***************************************************************************
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
 *				Javier Vargas (jvargas@cnb.csic.es)
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
    pcaAnalyzer[2];
    PCAMahalanobisAnalyzer tempPcaAnalyzer0;
    PCAMahalanobisAnalyzer tempPcaAnalyzer1;

    tempPcaAnalyzer0.clear();
    tempPcaAnalyzer1.clear();

    Image<double> img;
    MultidimArray<double> img2;
    Matrix1D<double> center(2);
    center.initZeros();
    FringeProcessing fp;

    int numHistElem = 31;
    int numDescriptors =5;

    MultidimArray<float> v(5);

    if (verbose>0)
    {
        std::cout << " Sorting particle set by new xmipp method..." << std::endl;
    }

    int nr_imgs = SF.size();
    if (verbose>0)
        init_progress_bar(nr_imgs);

    int c = XMIPP_MAX(1, nr_imgs / 60);
    int imgno = 0, imgnoPCA=0;

    //Zscore.initZeros(SF.size());
    bool thereIsEnable=SF.containsLabel(MDL_ENABLED);
    bool first=true;

    // We assume that at least there is one particle
    img.readApplyGeo(SF,1);
    MultidimArray<double> nI, modI;
    MultidimArray<bool> mask;
    nI.resizeNoCopy(img());
    modI.resizeNoCopy(img());
    mask.resizeNoCopy(img());
    mask.initConstant(true);
    MultidimArray<double> autoCorr(2*img().ydim,2*img().xdim);
    MultidimArray<double> smallAutoCorr;
    Histogram1D hist;
    Matrix2D<double> U,V,temp;
    Matrix1D<double> D;
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        if (thereIsEnable)
        {
            int enabled;
            SF.getValue(MDL_ENABLED,enabled,__iter.objId);
            if (enabled==-1)
            {
                //Zscore(imgno)=1000;
                imgno++;
                continue;
            }

            img.readApplyGeo(SF,__iter.objId);
            MultidimArray<double> &mI=img();
            mI.setXmippOrigin();
            mI.statisticsAdjust(0,1);
            mask.setXmippOrigin();

            auto_correlation_matrix(mI,autoCorr);
            mI.setXmippOrigin();
            autoCorr.window(smallAutoCorr,-15,-15, 15, 15);
            smallAutoCorr.copy(temp);
            svdcmp(temp,U,D,V);


            //Here is done all the processing. Here we will need probably
            //some input arguments to tune the desired frequencies and tune
            //the filer
            //ftrans.auto_correlation_matrix
            fp.normalize(mI,nI,modI,0.8,2,mask);
            nI.binarize();
            int im = labelImage2D(nI,nI,8);
            compute_hist(nI, hist, 0, im, im+1);

            int l,k,i,j;

            //We supose that the biggest part if the background!
            //This can be problematic
            hist.maxIndex(l,k,i,j);
            A1D_ELEM(hist,j)=0;
            hist.maxIndex(l,k,i,j);
            nI.binarizeRange(j-1,j+1);

            FringeProcessing fp;
            double x0=0,y0=0,majorAxis=0,minorAxis=0,ellipAng=0,area=0;
            fp.fitEllipse(nI,x0,y0,majorAxis,minorAxis,ellipAng,area);

            // Build vector
            v.initZeros(numDescriptors);
            v(0)=x0;
            v(1)=y0;
            v(2)=majorAxis;
            v(3)=minorAxis;
            v(4)=area;

            //auto_correlation_vector(mI,mI);
            /*          int idx = 5;

                        double minI=0, maxI=0;
                        mI.computeDoubleMinMax(minI,maxI);
                        compute_hist(mI,hist,minI,maxI,numHistElem);
                        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(hist)
                        v(idx++)=(float)DIRECT_A1D_ELEM(hist,i);
            */
            tempPcaAnalyzer0.addVector(v);
            //pcaAnalyzer.addVector(v);

            if (imgno % c == 0 && verbose>0)
                progress_bar(imgno);

            if (imgno == 248)
            {
                std::cout << im << std::endl;
                std::cout << "Max hist " << hist.hmax << std::endl;
                std::cout << "Min hist " << hist.hmin << std::endl;

                std::cout << "x0 " << x0 << std::endl;
                std::cout << "y0 " << y0 << std::endl;
                std::cout << "majorAxis " << majorAxis << std::endl;
                std::cout << "minorAxis " << minorAxis << std::endl;
                std::cout << "ellipAng " << ellipAng << std::endl;

                FileName fpName    = "test.txt";
                mI.write(fpName);
                fpName    = "test3.txt";
                smallAutoCorr.write(fpName);
            }

            imgno++;
            imgnoPCA++;
        }
    }

    MultidimArray<double> vavg,vstddev;
    tempPcaAnalyzer0.computeStatistics(vavg,vstddev);
    tempPcaAnalyzer1.computeStatistics(vavg,vstddev);

    tempPcaAnalyzer0.evaluateZScore(2,20);

    pcaAnalyzer.insert(pcaAnalyzer.begin(), tempPcaAnalyzer0);
    pcaAnalyzer.insert(pcaAnalyzer.begin(), tempPcaAnalyzer1);

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

    /*
    //Process input selfile ..............................................
    SF.read(fn);
    SF.removeDisabled();
    pcaAnalyzer.clear();
    processInput2(SF);
    processInput(SF, false, multivariate);
    */

    // Process input selfile ..............................................
    SF.read(fn);
    SF.removeDisabled();

    if (fn_train != "")
        SFtrain.read(fn_train);

    else
    	processInputPrepare(SF);

    int imgno = 0;
    int numPCAs = pcaAnalyzer.size();
    MultidimArray<double> ZscoreMultivariate(SF.size());
    MultidimArray<double> weights(numPCAs);
    ZscoreMultivariate.initConstant(1);
    weights.initConstant(1);

    for (int num = 0; num < numPCAs; ++num)
		FOR_ALL_OBJECTS_IN_METADATA(SF)
		{
			ZscoreMultivariate(imgno)*=(pcaAnalyzer[num].getZscore(imgno))*weights(numPCAs);
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
