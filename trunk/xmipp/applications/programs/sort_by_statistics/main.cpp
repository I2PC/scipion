/***************************************************************************
 *
 * Authors: Sjors Scheres (scheres@cnb.csic.es)
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

#include <data/image.h>
#include <data/fftw.h>
#include <data/args.h>
#include <data/metadata.h>
#include <data/metadata_extension.h>
#include <classification/pca.h>

class Sort_junk_parameters
{

public:
    double cutoff;
    MultidimArray<double> avg, stddev, Zscore, ZscoreMultivariate;
    PCAMahalanobisAnalyzer pcaAnalyzer;
    FileName fn_out;

public:
    Sort_junk_parameters()
    {}

    void usage()
    {
        std::cout  << " A sorting program for identifying junk particles \n"
        << " Parameters:\n"
        << " -i <selfile>            : Selfile with (normalized) input images\n"
        << " [-o <root=\"sort_junk\">] : Output rootname\n"
        << " [-train <selfile>]      : Train on selfile with good particles\n"
        << " [-zcut <float=-1>]      : Cut-off for Z-scores (negative for no cut-off) \n"
        << " [-multivariate]         : Identify also multivariate outliers\n"
        ;
    }

    void process_selfile(MetaData &SF, bool do_prepare, bool multivariate)
    {
        Image<double> img;
        MultidimArray<std::complex<double> > IMG;
        MultidimArray<int> radial_count;
        Matrix1D<int> center(2);
        center.initZeros();

        if (do_prepare)
            std::cerr << " Processing training set ..." << std::endl;
        else
            std::cerr << " Sorting particle set ..." << std::endl;

        MultidimArray<double> Mrad, ampl2, rmean_ampl2;
        int dim;
        if (do_prepare)
        {
            ImgSize(SF, dim, dim);
            Mrad.resize(dim, dim);
            Mrad.setXmippOrigin();
            FOR_ALL_ELEMENTS_IN_ARRAY2D(Mrad)
            Mrad(i, j) = sqrt((double)(i * i + j * j));
        }

        int nr_imgs = SF.size();
        init_progress_bar(nr_imgs);
        int c = XMIPP_MAX(1, nr_imgs / 60);
        int imgno = 0;
        MultidimArray<float> v;
        v.initZeros(16);
        if (do_prepare)
        {
            Zscore.initZeros(SF.size());
            ZscoreMultivariate=Zscore;
        }
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            if (do_prepare)
            {
                FileName fn;
                SF.getValue(MDL_IMAGE,fn);
                img.read(fn);
                img().setXmippOrigin();

                // Overall statistics
                double mean, stddev, minval, maxval;
                img().computeStats(mean, stddev, minval, maxval);

                // Number of low or high-valued pixels
                double nlowpix, nhighpix, nradlow, nradhigh;
                nhighpix = nlowpix = 0.;
                nradhigh = nradlow = 0.;
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(img())
                {
                    if (dAij(img(), i, j) > stddev)
                        nhighpix += 1.;
                    if (dAij(img(), i, j) < -stddev)
                        nlowpix += 1.;
                    if (dAij(img(), i, j) > stddev)
                        nradhigh += dAij(Mrad, i, j);
                    if (dAij(img(), i, j) < -stddev)
                        nradlow += dAij(Mrad, i, j);
                }

                // Quadrant statistics
                double sum_quadsig, sum2_quadsig, sum_quadmean, sum2_quadmean;
                quadrant_stats(img, sum_quadsig, sum2_quadsig, sum_quadmean, sum2_quadmean);

                // Fourier-space stats
                XmippFftw transformer;
                transformer.FourierTransform(img(),IMG,false);
                FFT_magnitude(IMG, ampl2);
                CenterFFT(ampl2, true);
                ampl2 *= ampl2;
                rmean_ampl2.initZeros();
                radialAverage(ampl2, center, rmean_ampl2, radial_count);
                double ampl2_1, ampl2_2, ampl2_3, ampl2_4;
                ampl2_1 = dAi(rmean_ampl2, 1);
                ampl2_2 = dAi(rmean_ampl2, 2);
                ampl2_3 = dAi(rmean_ampl2, 3);
                ampl2_4 = dAi(rmean_ampl2, 4);

                v( 0)=(float)mean;
                v( 1)=(float)stddev;
                v( 2)=(float)minval;
                v( 3)=(float)maxval;
                v( 4)=(float)nhighpix;
                v( 5)=(float)nlowpix;
                v( 6)=(float)nradhigh;
                v( 7)=(float)nradlow;
                v( 8)=(float)sum_quadsig;
                v( 9)=(float)sum_quadmean;
                v(10)=(float)sum2_quadsig;
                v(11)=(float)sum2_quadmean;
                v(12)=(float)ampl2_1;
                v(13)=(float)ampl2_2;
                v(14)=(float)ampl2_3;
                v(15)=(float)ampl2_4;
                pcaAnalyzer.addVector(v);
            }
            else
            {
                v=pcaAnalyzer.v[imgno];
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(v)
                {
                    if (DIRECT_A1D_ELEM(stddev,i)>0)
                    {
                        DIRECT_A1D_ELEM(v,i)=(DIRECT_A1D_ELEM(v,i)-DIRECT_A1D_ELEM(avg,i))/
                                             DIRECT_A1D_ELEM(stddev,i);
                        DIRECT_A1D_ELEM(v,i)=ABS(DIRECT_A1D_ELEM(v,i));
                    }
                    else
                        DIRECT_A1D_ELEM(v,i)=0;
                }
                Zscore(imgno)=v.computeAvg();
                if (multivariate)
                    ZscoreMultivariate(imgno)=pcaAnalyzer.getZscore(imgno);
            }

            if (imgno % c == 0)
                progress_bar(imgno);
            imgno++;
        }
        progress_bar(nr_imgs);

        if (do_prepare)
        {
            pcaAnalyzer.computeStatistics(avg,stddev);
            if (multivariate)
                pcaAnalyzer.evaluateZScore(2,10);
        }
        else
        	Zscore.rangeAdjust(0,1);
    }

    void quadrant_stats(Image<double> &img,
                        double &sum_quadsig, double &sum2_quadsig,
                        double &sum_quadmean, double &sum2_quadmean)
    {
        double mean, stddev, minval, maxval;
        Matrix1D<int> corner1(2), corner2(2);

        sum_quadsig = sum2_quadsig = 0.;
        sum_quadmean = sum2_quadmean = 0.;
        XX(corner1) = STARTINGX(img());
        YY(corner1) = STARTINGY(img());
        XX(corner2) = 0;
        YY(corner2) = 0;
        img().computeStats(mean, stddev, minval, maxval, corner1, corner2);
        sum2_quadmean += mean * mean;
        sum_quadmean += mean;
        sum2_quadsig += stddev * stddev;
        sum_quadsig += stddev;
        XX(corner1) = STARTINGX(img());
        YY(corner1) = 0;
        XX(corner2) = 0;
        YY(corner2) = FINISHINGY(img());
        img().computeStats(mean, stddev, minval, maxval, corner1, corner2);
        sum2_quadmean += mean * mean;
        sum_quadmean += mean;
        sum2_quadsig += stddev * stddev;
        sum_quadsig += stddev;
        XX(corner1) = 0;
        YY(corner1) = STARTINGY(img());
        XX(corner2) = FINISHINGX(img());
        YY(corner2) = 0;
        img().computeStats(mean, stddev, minval, maxval, corner1, corner2);
        sum2_quadmean += mean * mean;
        sum_quadmean += mean;
        sum2_quadsig += stddev * stddev;
        sum_quadsig += stddev;
        XX(corner1) = 0;
        YY(corner1) = 0;
        XX(corner2) = FINISHINGX(img());
        YY(corner2) = FINISHINGX(img());
        img().computeStats(mean, stddev, minval, maxval, corner1, corner2);
        sum2_quadmean += mean * mean;
        sum_quadmean += mean;
        sum2_quadsig += stddev * stddev;
        sum_quadsig += stddev;

        sum_quadsig /= 4.;
        sum2_quadsig = sqrt(sum2_quadsig / 4. - sum_quadsig * sum_quadsig);
        sum_quadmean /= 4.;
        sum2_quadmean = sqrt(sum2_quadmean / 4. - sum_quadmean * sum_quadmean);
    }
};

/**************************************************************************
        Main
/**************************************************************************/
int main(int argc, char **argv)
{
    MetaData SF, SFtrain;
    bool multivariate;

    // Read input parameters ............................................
    FileName fn, fn_train;
    Sort_junk_parameters prm;
    try
    {
        fn = getParameter(argc, argv, "-i");
        SF.read(fn);
        prm.fn_out = getParameter(argc, argv, "-o", "sort_junk");
        fn_train = getParameter(argc, argv, "-train", "");
        multivariate = checkParameter(argc, argv, "-multivariate");
        if (fn_train != "")
            SFtrain.read(fn_train);
        prm.cutoff = textToFloat(getParameter(argc, argv, "-zcut", "-1"));
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        prm.usage();
        exit(1);
    }

    try
    {
        // Process input selfile ..............................................
        if (fn_train != "")
            prm.process_selfile(SFtrain, true, multivariate);
        else
            prm.process_selfile(SF, true, multivariate);
        prm.process_selfile(SF, false, multivariate);

        // Produce output .....................................................
        MetaData SFout, SFoutGood;
        std::ofstream fh_zind;
        fh_zind.open((prm.fn_out + ".indZ").c_str(), std::ios::out);
        MultidimArray<double> finalZscore=prm.Zscore;
        if (multivariate)
        	finalZscore+=prm.ZscoreMultivariate;

        MultidimArray<int> sorted = finalZscore.indexSort();
        fh_zind << "image : avg, stddev, min, max, nhighpix, nlowpix, nradhigh, nradlow, quadsig, quadmean, ampl2_1, ampl2_2, ampl2_3, ampl2_4 ";
        if (multivariate)
        	fh_zind << "MultivariateZscore ";
        fh_zind << std::endl;
        int nr_imgs = SF.size();
        for (int imgno = 0; imgno < nr_imgs; imgno++)
        {
            int isort = sorted(imgno) - 1;
            FileName fnImg;
            SF.getValue(MDL_IMAGE,fnImg,isort+1);
            SFout.addObject();
            SFout.setValue(MDL_IMAGE,fnImg);
            SFout.setValue(MDL_ENABLED,1);
            SFout.setValue(MDL_ZSCORE,finalZscore(isort));
            if (prm.Zscore(isort)<prm.cutoff && prm.cutoff>0)
            {
                SFoutGood.addObject();
                SFoutGood.setValue(MDL_IMAGE,fnImg);
                SFoutGood.setValue(MDL_ENABLED,1);
                SFoutGood.setValue(MDL_ZSCORE,finalZscore(isort));
            }
            fh_zind << fnImg << " : ";
            FOR_ALL_ELEMENTS_IN_ARRAY1D(prm.pcaAnalyzer.v[isort])
            fh_zind << prm.pcaAnalyzer.v[isort](i) << "\t";
            if (multivariate)
            	fh_zind << prm.pcaAnalyzer.getSortedZscore(isort);
            fh_zind << std::endl;
        }
        fh_zind.close();
        SFout.write(prm.fn_out + ".sel");
        if (prm.cutoff>0)
            SFoutGood.write(prm.fn_out + "_good.sel");
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        exit(1);
    }
    return 0;
}
