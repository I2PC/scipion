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
#include <graphics/show_2d.h>
#include <graphics/show_selfile.h>
#include <graphics/show_vol.h>
#include <graphics/show_spectra.h>
#include <graphics/show_som.h>
#include <graphics/show_spectra_som.h>
#include <graphics/show_cl2d.h>

#include <qapplication.h>

void Usage();

#define MODE_IMG      0
#define MODE_SEL      1
#define MODE_VOL      2
#define MODE_SPECT    3
#define MODE_SOM      4
#define MODE_SPECTSOM 5
#define MODE_PSD      6
#define MODE_CTF      7
#define MODE_PSDSEL   8
#define MODE_CTFSEL   9
#define MODE_CL2D    10

int main(int argc, char **argv)
{
    int numCols, numRows, mode, ifirst;
    FileName fn_dat, fn_assign, fn_assignsel;
    bool poll, apply_geo = true, common_normalization = false;
    std::string filterSuffix;

    try
    {
        if (checkParameter(argc, argv, "-img"))
        {
            mode = MODE_IMG;
            ifirst = paremeterPosition(argc, argv, "-img");
        }
        else if (checkParameter(argc, argv, "-psdsel"))
        {
            mode = MODE_PSDSEL;
            ifirst = paremeterPosition(argc, argv, "-psdsel");
        }
        else if (checkParameter(argc, argv, "-ctfsel"))
        {
            mode = MODE_CTFSEL;
            ifirst = paremeterPosition(argc, argv, "-ctfsel");
            fn_assign = getParameter(argc, argv, "-assign", "");
            fn_assignsel = getParameter(argc, argv, "-assignsel", "");
        }
        else if (checkParameter(argc, argv, "-sel"))
        {
            mode = MODE_SEL;
            ifirst = paremeterPosition(argc, argv, "-sel");
        }
        else if (checkParameter(argc, argv, "-vol"))
        {
            mode = MODE_VOL;
            ifirst = paremeterPosition(argc, argv, "-vol");
        }
        else if (checkParameter(argc, argv, "-spect"))
        {
            mode = MODE_SPECT;
            ifirst = paremeterPosition(argc, argv, "-spect");
        }
        else if (checkParameter(argc, argv, "-som"))
        {
            mode = MODE_SOM;
            ifirst = paremeterPosition(argc, argv, "-som");
        }
        else if (checkParameter(argc, argv, "-cl2d"))
        {
            mode = MODE_CL2D;
            ifirst = paremeterPosition(argc, argv, "-cl2d");
            filterSuffix = getParameter(argc,argv,"-filterSuffix","");
        }
        else if (checkParameter(argc, argv, "-psd"))
        {
            mode = MODE_PSD;
            ifirst = paremeterPosition(argc, argv, "-psd");
        }
        else if (checkParameter(argc, argv, "-ctf"))
        {
            mode = MODE_CTF;
            ifirst = paremeterPosition(argc, argv, "-ctf");
            fn_assign = getParameter(argc, argv, "-assign", "");
        }
        else if (checkParameter(argc, argv, "-spectsom"))
        {
            mode = MODE_SPECTSOM;
            ifirst = paremeterPosition(argc, argv, "-spectsom");
            fn_dat = getParameter(argc, argv, "-din");
        }
        else
            REPORT_ERROR(1, "No mode (img/sel/vol) supplied");
        numCols = textToInteger(getParameter(argc, argv, "-w", "-1"));
        numRows = textToInteger(getParameter(argc, argv, "-h", "-1"));
        apply_geo = !checkParameter(argc, argv, "-dont_apply_geo");
        poll = checkParameter(argc, argv, "-poll");
        common_normalization = checkParameter(argc, argv, "-common_norm");
    }
    catch (Xmipp_error)
    {
        Usage();
        exit(1);
    }

    try
    {
        QApplication::setFont(QFont("Helvetica", 12));
        QApplication a(argc, argv);

        // Get common normalization
        double m = 0, M = 0;
        if (common_normalization)
        {
            for (int i = ifirst + 1; i < argc; i++)
            {
                if (!exists(argv[i]))
                {
                    if (argv[i][0] == '-') break; // There is nothing else to show
                    if (mode == MODE_VOL)
                    {
                        FileName fn = argv[i];
                        if (fn[fn.length()-1] == 'x' || fn[fn.length()-1] == 'y')
                            fn = fn.substr(0, fn.length() - 1);
                        if (exists(fn.c_str()))
                        {
                            VolumeXmipp V(fn);
                            double maux, Maux;
                            V().computeDoubleMinMax(maux, Maux);
                            if (i == ifirst + 1)
                            {
                                m = maux;
                                M = Maux;
                            }
                            else
                            {
                                m = XMIPP_MIN(m, maux);
                                M = XMIPP_MAX(M, Maux);
                            }
                        }
                    }
                }
                else
                {
                    if (mode == MODE_IMG)
                    {
                        ImageXmipp I(argv[i]);
                        double maux, Maux;
                        I().computeDoubleMinMax(maux, Maux);
                        if (i == ifirst + 1)
                        {
                            m = maux;
                            M = Maux;
                        }
                        else
                        {
                            m = XMIPP_MIN(m, maux);
                            M = XMIPP_MAX(M, Maux);
                        }
                    }
                    else if (mode == MODE_SEL)
                    {
                        // Not implemented
                    }
                    else if (mode == MODE_VOL)
                    {
                        VolumeXmipp V(argv[i]);
                        double maux, Maux;
                        V().computeDoubleMinMax(maux, Maux);
                        if (i == ifirst + 1)
                        {
                            m = maux;
                            M = Maux;
                        }
                        else
                        {
                            m = XMIPP_MIN(m, maux);
                            M = XMIPP_MAX(M, Maux);
                        }
                    }
                    else if (mode == MODE_SPECT)
                    {
                        // Not implemented
                    }
                    else if (mode == MODE_SOM)
                    {
                        // Not implemented
                    }
                    else if (mode == MODE_CL2D)
                    {
                        // Not implemented
                    }
                    else if (mode == MODE_PSD)
                    {
                        // Not implemented
                    }
                    else if (mode == MODE_CTF)
                    {
                        // Not implemented
                    }
                    else if (mode == MODE_SPECTSOM)
                    {
                        // Not implemented
                    }
                }
            }
        }

        // Show
        int shown = 0;
        for (int i = ifirst + 1; i < argc; i++)
        {
            if (!exists(argv[i]))
            {
                if (argv[i][0] == '-') break; // There is nothing else to show
                FileName fn;
                switch (mode)
                {
                case MODE_IMG:
                    fn = argv[i];
                    if (fn.find("imagic:") != -1) break;
                    std::cerr << argv[i] << " is not a valid filename\n";
                    continue;
                case MODE_PSDSEL:
                case MODE_CTFSEL:
                case MODE_SEL:
                    std::cerr << argv[i] << " is not a valid filename\n";
                    continue;
                case MODE_VOL:
                    fn = argv[i];
                    if (fn[fn.length()-1] == 'x' || fn[fn.length()-1] == 'y')
                    {
                        fn = fn.substr(0, fn.length() - 1);
                        if (!exists(fn.c_str()))
                        {
                            std::cerr << fn << " is not a valid filename\n";
                            continue;
                        }
                    }
                    else
                    {
                        continue;
                        std::cerr << fn << " is not a valid filename\n";
                    }
                    break;
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
            if (mode == MODE_IMG)
            {
                ImageViewer *showimg = new ImageViewer(argv[i], poll);
                showimg->apply_geo = apply_geo;
                showimg->loadImage(argv[i]);
                showimg->show();
                shown++;
            }
            else if (mode == MODE_PSD)
            {
                ImageViewer *showimg = new ImageViewer(argv[i], poll);
                showimg->apply_geo = false;
                showimg->loadImage(argv[i], 0, 0, ImageViewer::PSD_mode);
                showimg->show();
                shown++;
            }
            else if (mode == MODE_CTF)
            {
                ImageViewer *showimg = new ImageViewer(argv[i], poll);
                showimg->apply_geo = false;
                showimg->setAssignCTFfile(fn_assign);
                showimg->loadImage(argv[i], 0, 0, ImageViewer::CTF_mode);
                showimg->show();
                shown++;
            }
            else if (mode == MODE_SEL)
            {
                ShowSel *showsel = new ShowSel;
                showsel->apply_geo = apply_geo;
                showsel->showonlyactive = !checkParameter(argc, argv, "-showall");
                showsel->initWithFile(numRows, numCols, argv[i]);
                showsel->show();
                shown++;
            }
            else if (mode == MODE_PSDSEL)
            {
                ShowSel *showsel = new ShowSel;
                showsel->apply_geo = false;
                showsel->showonlyactive = false;
                showsel->initWithFile(numRows, numCols, argv[i], ShowSel::PSD_mode);
                showsel->show();
                shown++;
            }
            else if (mode == MODE_CTFSEL)
            {
                ShowSel *showsel = new ShowSel;
                showsel->apply_geo = false;
                showsel->showonlyactive = false;
                showsel->initWithFile(numRows, numCols, argv[i], ShowSel::CTF_mode);
                showsel->setAssignCTFfile(fn_assign);
                showsel->setAssignCTFselfile(fn_assignsel);
                showsel->show();
                shown++;
            }
            else if (mode == MODE_VOL)
            {
                ShowVol *showvol = new ShowVol;
                if (poll) showvol->setPoll();
                showvol->initWithFile(numRows, numCols, argv[i]);
                showvol->show();
                shown++;
            }
            else if (mode == MODE_SPECT)
            {
                ShowSpectra *showspectra = new ShowSpectra;
                showspectra->initWithFile(numRows, numCols, argv[i]);
                showspectra->show();
                shown++;
            }
            else if (mode == MODE_SOM)
            {
                ShowSOM *showsom = new ShowSOM;
                showsom->apply_geo = apply_geo;
                showsom->initWithFile(argv[i]);
                showsom->show();
                shown++;
            }
            else if (mode == MODE_CL2D)
            {
                ShowCL2D *showcl2d = new ShowCL2D;
                showcl2d->apply_geo = apply_geo;
                showcl2d->filterSuffix = filterSuffix;
                showcl2d->initWithFile(argv[i]);
                showcl2d->show();
                shown++;
            }
            else if (mode == MODE_SPECTSOM)
            {
                ShowSpectraSOM *showspectrasom = new ShowSpectraSOM;
                showspectrasom->apply_geo = apply_geo;
                showspectrasom->initWithFile(argv[i], fn_dat);
                showspectrasom->show();
                shown++;
            }
        }

        if (!shown) return 0;

        QObject::connect(qApp, SIGNAL(lastWindowClosed()), qApp, SLOT(quit()));
        return a.exec();
    }
    catch (Xmipp_error XE)
    {
        std::cerr << XE;
    }
}

void Usage()
{
    std::cout << "Usage: show [options]\n"
    << "    -img <images> |       : Input images\n"
    << "    -ctf <images> |       : Input CTFs (in image format)\n"
    << "    -psd <images> |       : Input PSDs (in image format)\n"
    << "    -sel <selfiles> |     : Input selfiles\n"
    << "    -psdsel <selfiles> |  : Input PSD Selfiles (in image format)\n"
    << "    -ctfsel <selfiles> |  : Input CTF Selfiles (in image format)\n"
    << "      [-assign <filename>]: Input Assign CTF Parameters file\n"
    << "      [-assignsel <filename>]: Selfile of Input Assign CTF Parameters files\n"
    << "    -vol <XmippVolumes>   : Add x or y to the filename\n"
    << "                            to see slices in that direction\n"
    << "    -spect <datafile>     : Spectra .dat file\n"
    << "    -som <SOM rootname>   : SOM images\n"
    << "    -cl2d <CL2D rootname> : CL2D images\n"
    << "      [-filterSuffix <s>] : Filter suffix for the CL2D\n"
    << "    -spectsom <SOM root>  : SOM spectra\n"
    << "       -din <Original.dat>: Original data\n"
    << "   [-w]                   : width (default: 10)\n"
    << "   [-h]                   : height (default: 10)\n"
    << "   [-showall]             : only for sel mode, show all images even\n"
    << "                            if sel = -1\n"
    << "   [-poll]                : check file change, NOT for sel files\n"
    << "   [-dont_apply_geo]      : Do not apply transformation stored in the header of 2D-images\n"
    ;
}
