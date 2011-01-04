/***************************************************************************
 *
 * Authors: Carlos Oscar Sorzano (coss@cnb.csic.es)
 *          Sjors Scheres (scheres@cnb.csic.es)
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

#include <data/basic_pca.h>
#include <data/histogram.h>
#include <data/image.h>
#include <data/fftw.h>
#include <data/args.h>
#include <data/metadata.h>
#include <data/metadata_extension.h>

class Sort_junk_parameters
{

public:
    double cutoff;
    MultidimArray<double> vavg, vstddev, Zscore, ZscoreMultivariate;
    PCAMahalanobisAnalyzer pcaAnalyzer;
    FileName fn_out;

public:
    Sort_junk_parameters()
    {}

    void usage()
    {
        std::cout  << " A sorting program for identifying junk particles \n"
        << " Parameters:\n"
        << " -i <selfile>            : Selfile with input images\n"
        << " [--block <blockname>]   : Block within the input selfile\n"
        << " [-o <root=\"sort_junk\">] : Output rootname\n"
        << " [--train <selfile>]     : Train on selfile with good particles\n"
        << " [--zcut <float=-1>]     : Cut-off for Z-scores (negative for no cut-off) \n"
        << " [--multivariate]        : Identify also multivariate outliers\n"
        << " [--verbose]             : Save the vectors associated to each image\n"
        << " [--quiet]               : Do not show anything on screen\n"
        ;
    }

    void process_selfile(MetaData &SF, bool do_prepare, bool multivariate, bool quiet)
    {
        Image<double> img;
        MultidimArray<double> img2;
        MultidimArray<int> radial_count;
        MultidimArray<double> radial_avg;
        Matrix1D<int> center(2);
        center.initZeros();

        if (!quiet)
        {
            if (do_prepare)
                std::cerr << " Processing training set ..." << std::endl;
            else
                std::cerr << " Sorting particle set ..." << std::endl;
        }

        int nr_imgs = SF.size();
        if (!quiet)
            init_progress_bar(nr_imgs);
        int c = XMIPP_MAX(1, nr_imgs / 60);
        int imgno = 0, imgnoPCA=0;
        MultidimArray<float> v;
        MultidimArray<int> distance;
        int dim;
        if (do_prepare)
        {
            Zscore.initZeros(SF.size());
            ZscoreMultivariate=Zscore;
        }
        bool thereIsEnable=SF.containsLabel(MDL_ENABLED);
        bool first=true;
        FOR_ALL_OBJECTS_IN_METADATA(SF)
        {
            if (do_prepare)
            {
                FileName fn;
                SF.getValue(MDL_IMAGE,fn);
                if (thereIsEnable)
                {
                    int enabled;
                    SF.getValue(MDL_ENABLED,enabled);
                    if (enabled==-1)
                        continue;
                }
                img.read(fn);
                img().setXmippOrigin();
                img().statisticsAdjust(0,1);

                // Overall statistics
                Histogram1D hist;
                compute_hist(img(),hist,-4,4,31);

                // Radial profile
                img2.resizeNoCopy(img());
                FOR_ALL_ELEMENTS_IN_ARRAY2D(img2)
                {
                    double val=IMGPIXEL(img,i,j);
                    A2D_ELEM(img2,i,j)=val*val;
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
                pcaAnalyzer.addVector(v);
            }
            else
            {
                if (thereIsEnable)
                {
                    int enabled;
                    SF.getValue(MDL_ENABLED,enabled);
                    if (enabled==-1)
                    {
                        Zscore(imgno)=1000;
                        if (multivariate)
                            ZscoreMultivariate(imgno)=1000;
                        imgno++;
                        continue;
                    }
                }
                v=pcaAnalyzer.v[imgnoPCA];
                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(v)
                {
                    if (DIRECT_A1D_ELEM(vstddev,i)>0)
                    {
                        DIRECT_A1D_ELEM(v,i)=(DIRECT_A1D_ELEM(v,i)-DIRECT_A1D_ELEM(vavg,i))/
                                             DIRECT_A1D_ELEM(vstddev,i);
                        DIRECT_A1D_ELEM(v,i)=ABS(DIRECT_A1D_ELEM(v,i));
                    }
                    else
                        DIRECT_A1D_ELEM(v,i)=0;
                }
                Zscore(imgno)=v.computeAvg();
                if (multivariate)
                    ZscoreMultivariate(imgno)=pcaAnalyzer.getZscore(imgno);
            }

            if (imgno % c == 0 && !quiet)
                progress_bar(imgno);
            imgno++;
            imgnoPCA++;
        }
        if (!quiet)
            progress_bar(nr_imgs);

        if (do_prepare)
        {
            pcaAnalyzer.computeStatistics(vavg,vstddev);
            if (multivariate)
                pcaAnalyzer.evaluateZScore(2,20);
        }
    }
};

/**************************************************************************
        Main
/**************************************************************************/
int main(int argc, char **argv)
{
    MetaData SF, SFtrain;
    bool multivariate, verbose, quiet;

    // Read input parameters ............................................
    FileName fn, fn_train, fn_block;
    Sort_junk_parameters prm;
    try
    {
        fn = getParameter(argc, argv, "-i");
        if (checkParameter(argc,argv,"--block"))
            fn_block=getParameter(argc, argv, "--block");
        SF.read(fn,NULL,fn_block);
        prm.fn_out = getParameter(argc, argv, "-o", "sort_junk");
        fn_train = getParameter(argc, argv, "--train", "");
        multivariate = checkParameter(argc, argv, "--multivariate");
        verbose = checkParameter(argc, argv, "--verbose");
        quiet = checkParameter(argc, argv, "--quiet");
        if (fn_train != "")
            SFtrain.read(fn_train);
        prm.cutoff = textToFloat(getParameter(argc, argv, "--zcut", "-1"));
    }
    catch (XmippError XE)
    {
        std::cout << XE;
        prm.usage();
        exit(1);
    }

    try
    {
        // Process input selfile ..............................................
        if (fn_train != "")
            prm.process_selfile(SFtrain, true, multivariate, quiet);
        else
            prm.process_selfile(SF, true, multivariate, quiet);
        prm.process_selfile(SF, false, multivariate, quiet);

        // Produce output .....................................................
        MetaData SFout, SFoutGood;
        std::ofstream fh_zind;
        if (verbose)
            fh_zind.open((prm.fn_out + ".indZ").c_str(), std::ios::out);
        MultidimArray<double> finalZscore=prm.Zscore;
        double L=1;
        if (multivariate)
        {
            finalZscore*=prm.ZscoreMultivariate;
            L++;
        }

        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(finalZscore)
        DIRECT_A1D_ELEM(finalZscore,i)=pow(DIRECT_A1D_ELEM(finalZscore,i),1.0/L);

        MultidimArray<int> sorted = finalZscore.indexSort();
        int nr_imgs = SF.size();
        bool thereIsEnable=SF.containsLabel(MDL_ENABLED);
        for (int imgno = 0; imgno < nr_imgs; imgno++)
        {
            int isort = sorted(imgno) - 1;
            FileName fnImg;
            SF.getValue(MDL_IMAGE,fnImg,isort+1);
            if (thereIsEnable)
            {
                int enabled;
                SF.getValue(MDL_ENABLED,enabled,isort+1);
                if (enabled==-1)
                    continue;
            }
            SFout.addObject();
            SFout.setValue(MDL_IMAGE,fnImg);
            SFout.setValue(MDL_ENABLED,1);
            SFout.setValue(MDL_ZSCORE,finalZscore(isort));
            if (finalZscore(isort)<prm.cutoff && prm.cutoff>0)
            {
                SFoutGood.addObject();
                SFoutGood.setValue(MDL_IMAGE,fnImg);
                SFoutGood.setValue(MDL_ENABLED,1);
                SFoutGood.setValue(MDL_ZSCORE,finalZscore(isort));
            }
            if (verbose)
            {
                fh_zind << fnImg << " : ";
                FOR_ALL_ELEMENTS_IN_ARRAY1D(prm.pcaAnalyzer.v[isort])
                fh_zind << prm.pcaAnalyzer.v[isort](i) << "\t";
            }
        }
        if (verbose)
            fh_zind.close();
        SFout.write(prm.fn_out + ".sel");
        if (prm.cutoff>0)
            SFoutGood.write(prm.fn_out + "_good.sel");
    }
    catch (XmippError XE)
    {
        std::cout << XE;
        exit(1);
    }
    return 0;
}
