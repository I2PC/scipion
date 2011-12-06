/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.csic.es)
 *              Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es)
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

#include <data/args.h>
#include <data/xmipp_program.h>
#include <graphics/show_2d.h>
#include <graphics/show_selfile.h>
#include <graphics/show_vol.h>
#include <graphics/show_spectra.h>
#include <graphics/show_som.h>
#include <graphics/show_spectra_som.h>
#include <graphics/show_cl2d.h>

#include <qapplication.h>


enum ShowMode { MODE_INPUT, MODE_SPECT, MODE_SOM, MODE_SPECTSOM,
                MODE_PSD, MODE_CL2D };

class ProgShow: public XmippProgram
{
    ShowMode mode;
    int numCols, numRows, ifirst;
    FileName fn_dat;
    bool poll, apply_geo, common_normalization;
    String filterSuffix;
    StringVector files;

protected:
    void defineParams()
    {
        addUsageLine("Provides a Graphical User Interface for visualization and manipulation of");
        addUsageLine("Electron Microscopy Single Particle Images. The program can show images,");
        addUsageLine("volumes, stacks or selfiles. Also are provided some specific visualizations");
        addUsageLine("(like som, cl2d, spectra...etc)");
        addParamsLine("    -i <...>                   : Input files: accept images, volumes, stacks or metadatas");
        addParamsLine("or  --psd <...>                : Input PSDs (in image format)");
        addParamsLine("or  --spect <datafile>         : Spectra .dat file");
        addParamsLine("or  --som <SOM_rootname>       : SOM images");
        addParamsLine("or  --cl2d <CL2D_rootname>     : CL2D images");
        addParamsLine("or  --spectsom <SOM_rootname>  : SOM spectra");
        addParamsLine("[--spectsom_din <filename>]    : Original data");
        addParamsLine("      requires --spectsom;");
        addParamsLine("[--cl2d_filter_suffix <s>]     : Filter suffix for the CL2D");
        addParamsLine("      requires --cl2d;");
        addParamsLine("   [--dim <w=10> <h=10>]       : Dimensions of the table");
        addParamsLine("   [--common_norm]             : Normalize all the volumes or images with the same factors");
        addParamsLine("   [--showall]                 : Only for sel mode, show all images even if sel = -1");
        addParamsLine("   [--poll]                    : Check file change, NOT for sel files");
        addParamsLine("   [--dont_apply_geo]          : Do not apply transformations of 2D-images stored in metadata");
    }

    void readParams()
    {
        if (checkParam("-i"))
        {
            mode = MODE_INPUT;
            getListParam("-i", files);
        }
        else if (checkParam("--spect"))
        {
            mode = MODE_SPECT;
            getListParam("--spect", files);
        }
        else if (checkParam("--som"))
        {
            mode = MODE_SOM;
            getListParam("--som", files);
        }
        else if (checkParam("--cl2d"))
        {
            mode = MODE_CL2D;
            getListParam("--cl2d", files);
            filterSuffix = getParam("--cl2d_filter_suffix");
        }
        else if (checkParam("--psd"))
        {
            mode = MODE_PSD;
            getListParam("--psd", files);
        }
        else if (checkParam("--spectsom"))
        {
            mode = MODE_SPECTSOM;
            //getListParam("--spectsom", files);
            //ifirst = paremeterPosition(argc, argv, "-spectsom");
            fn_dat = getParam("--spectsom_din");
        }
        else
            REPORT_ERROR(ERR_ARG_MISSING, "No mode (img/sel/vol) supplied");
        numCols = getIntParam("--dim", 0);
        numRows = getIntParam("--dim", 1);
        apply_geo = !checkParam("--dont_apply_geo");
        poll = checkParam("--poll");
        common_normalization = checkParam("--common_norm");
    }

    void run()
    {
        QApplication::setFont(QFont("Helvetica", 12));
        QApplication a(argc, argv);

        // Get common normalization
        double m = 1e38, M = -m;
        Image<double> I;
        FileName fn;
        MetaData md;

        double maux, Maux;

        if (common_normalization)
        {
            //            for (StringVector::const_iterator it = files.begin(); it < files.end(); ++it)
            //            {
            //                fn = *it;
            //                if (exists(fn))
            //                {
            //                    if (mode == MODE_INPUT)
            //                    {
            //                        I.read(fn);
            //                        I().computeDoubleMinMax(maux, Maux);
            //                        m = XMIPP_MIN(m, maux);
            //                        M = XMIPP_MAX(M, Maux);
            //                    }
            //                    else if (mode == MODE_SEL)
            //                    {
            //                        // Not implemented
            //                    }
            //                    else if (mode == MODE_SPECT)
            //                    {
            //                        // Not implemented
            //                    }
            //                    else if (mode == MODE_SOM)
            //                    {
            //                        // Not implemented
            //                    }
            //                    else if (mode == MODE_CL2D)
            //                    {
            //                        // Not implemented
            //                    }
            //                    else if (mode == MODE_PSD)
            //                    {
            //                        // Not implemented
            //                    }
            //                    else if (mode == MODE_CTF)
            //                    {
            //                        // Not implemented
            //                    }
            //                    else if (mode == MODE_SPECTSOM)
            //                    {
            //                        // Not implemented
            //                    }
            //                }
            //            }
        }

        // Show
        int shown = 0;
        for (StringVector::const_iterator it = files.begin(); it < files.end(); ++it)
        {
            fn = *it;
            if (!fn.existsTrim())
            {
                switch (mode)
                {
                case MODE_INPUT:
                case MODE_SPECT:
                    continue;
                case MODE_SOM:
                    break;
                case MODE_CL2D:
                    break;
                case MODE_SPECTSOM:
                    break;
                }
            }
            if (mode == MODE_INPUT)
            {
                try
                {
                    md.read(fn);
                    if (md.size() > 1)//more than one object
                    {
                        ShowSel *showsel = new ShowSel;
                        showsel->apply_geo = apply_geo;
                        showsel->showonlyactive = !checkParam("--showall");
                        md.read(fn);
                        showsel->initWithObject(numRows, numCols, md, fn.c_str());
                        showsel->show();
                    }
                    else
                    {
                        md.getValue(MDL_IMAGE, fn, md.firstObject());
                        Image<char> img;
                        img.read(fn, HEADER);
                        if (img().zdim == 1)
                        {
                            ImageViewer *showimg = new ImageViewer(fn.c_str(), poll);
                            showimg->apply_geo = apply_geo;
                            showimg->loadImage(fn.c_str());
                            showimg->show();
                        }
                        else
                        {
                            ShowVol *showvol = new ShowVol;
                            if (poll)
                                showvol->setPoll();
                            showvol->initWithFile(numRows, numCols, fn);
                            showvol->show();
                        }
                    }
                }
                catch (XmippError xe)
                {
                  std::cerr << xe << std::endl;
                }
                shown++;
            }
            else if (mode == MODE_PSD)
            {
                ImageViewer *showimg = new ImageViewer(fn.c_str(), poll);
                showimg->apply_geo = false;
                showimg->loadImage(fn.c_str(), 0, 0, ImageViewer::PSD_mode);
                showimg->show();
                shown++;
            }
            else if (mode == MODE_SPECT)
            {
                ShowSpectra *showspectra = new ShowSpectra;
                showspectra->initWithFile(numRows, numCols, fn);
                showspectra->show();
                shown++;
            }
            else if (mode == MODE_SOM)
            {
                ShowSOM *showsom = new ShowSOM;
                showsom->apply_geo = apply_geo;
                showsom->initWithFile(fn);
                showsom->show();
                shown++;
            }
            else if (mode == MODE_CL2D)
            {
                ShowCL2D *showcl2d = new ShowCL2D;
                showcl2d->apply_geo = apply_geo;
                showcl2d->filterSuffix = filterSuffix;
                showcl2d->initWithFile(fn);
                showcl2d->show();
                shown++;
            }
            else if (mode == MODE_SPECTSOM)
            {
                ShowSpectraSOM *showspectrasom = new ShowSpectraSOM;
                showspectrasom->apply_geo = apply_geo;
                showspectrasom->initWithFile(fn, fn_dat);
                showspectrasom->show();
                shown++;
            }
        }

        if (!shown)
            return;

        QObject::connect(qApp, SIGNAL(lastWindowClosed()), qApp, SLOT(quit()));
        a.exec();

    }
}
;//end of class ProgShow

int main(int argc, char **argv)
{
    ProgShow prm;
    prm.read(argc, argv);
    prm.tryRun();
}

