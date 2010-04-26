/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Carlos Manzanares       (cmanzana@cnb.csic.es)
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

#include <qapplication.h>

#include "main_widget_mark.h"
#include "widget_psd.h"
#include "widget_micrograph.h"

void Usage();

int main(int argc, char **argv)
{
    FileName fnRaw;
    FileName fnRawTilted;
    FileName fnAutomaticModel;
    bool     reversed;
    FileName fn_assign_CTF;
    bool     ctf_mode = false;
    bool     autoSelect = false;
    int      numThreads = 1;

    // Get input parameters .................................................
    try
    {
        fnRaw         = getParameter(argc, argv, "-i");
        fnRawTilted   = getParameter(argc, argv, "-tilted", "");
        reversed      = checkParameter(argc, argv, "-reverse_endian");
        fn_assign_CTF = getParameter(argc, argv, "-psd", "");
        if (checkParameter(argc, argv, "-ctf"))
        {
            ctf_mode = true;
            fn_assign_CTF = getParameter(argc, argv, "-ctf");
        }
        fnAutomaticModel = getParameter(argc, argv, "-auto", "");
        autoSelect = checkParameter(argc, argv, "-autoSelect");
        if (fnRawTilted!="" && autoSelect)
            REPORT_ERROR(1,"Automatic particle picking cannot be performed on tilt pairs");
        numThreads = textToInteger(getParameter(argc, argv, "-thr","1"));
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage();
        exit(1);
    }

    try
    {
        Micrograph m, mTilted;
        FileName fn8bits="", fn8bitsTilted="";

        m.open_micrograph(fnRaw, reversed);
	// Sjors & Roberto: 25jan08
	// The following is a "chapuza" because the 8-bit
	// visualization routine from XVsmooth is much nicer than our
	// own routine. We save the micrograph with name fn8bits,
	// re-read it in 8bit-format, and then set the name explicitly
	// to the original fnRaw to get the correct .pos, .ang etc
	// files.
	// Note that this also requires that the name of the original
	// micrograph in the dialog window under "generate images"
	// should be the original fnRaw (and not be left blank as before!)
        if (m.depth()!=8)
        {
            fn8bits=fnRaw+".8bits";
            m.write_as_8_bits(fn8bits);
            m.close_micrograph();
            m.open_micrograph(fn8bits,false);
	    m.set_micrograph_name(fnRaw);
            m.compute_8_bit_scaling();
	    system(((std::string)"rm -rf "+fn8bits+"*").c_str());
        }
	else
	{
	    m.compute_8_bit_scaling();
	}

        if (fnRawTilted != "")
        {
            mTilted.open_micrograph(fnRawTilted, reversed);
            if (mTilted.depth()!=8)
            {
                fn8bitsTilted=fnRawTilted+".8bits";
                mTilted.write_as_8_bits(fn8bitsTilted);
                mTilted.close_micrograph();
                mTilted.open_micrograph(fn8bitsTilted,false);
		mTilted.set_micrograph_name(fnRawTilted);
                mTilted.compute_8_bit_scaling();
		system(((std::string)"rm -rf "+fn8bitsTilted+"*").c_str());
            }
	    else
	    {
		mTilted.compute_8_bit_scaling();
	    }
        }

        // Configure application .............................................
        AutoParticlePicking *autoPicking=NULL;
        autoPicking=new AutoParticlePicking(&m);

        QApplication *app=NULL;
        QtMainWidgetMark *mainWidget=NULL;
        if (!autoSelect)
        {
            app=new QApplication(argc, argv);
            if (fnRawTilted == "")
                mainWidget = new QtMainWidgetMark(&m);
            else
                mainWidget = new QtMainWidgetMark(&m, &mTilted);
            // Sjors 26apr2010: always set autoPicking to prevent problems with tilt pairs...
            mainWidget->untilted_widget()->
                setAutoParticlePicking(autoPicking);
        }

        // Check if the PSDs must be shown ...................................
        if (fn_assign_CTF != "")
        {
            QtWidgetPSD PSDshow;
            if (ctf_mode) PSDshow.set_CTF_mode();
            PSDshow.set_assign_CTF_file(m, fn_assign_CTF);
            PSDshow.show();
        }

        // Check if a model has been provided ................................
        autoPicking->setNumThreads(numThreads);
        if (fnAutomaticModel!="")
            autoPicking->loadModels(fnAutomaticModel);

        // Run application ...................................................
        if (!autoSelect)
        {
            app->setMainWidget(mainWidget);
            mainWidget->openAllWindows();
            app->exec();
        }
        else
        {
            autoPicking->automaticallySelectParticles();
            autoPicking->saveAutoParticles();
        }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
    return 0;
}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    std::cerr << "Purpose: Mark particles in a Raw image\n"
              << "         There must exist the image and the corresponding .inf file\n"
              << "\n"
              << "Usage: mark [options]\n"
              << "   -i <input raw file>                : File with the image\n"
              << "  [-tilted <tilted raw file>]         : Image with the tilted pair\n"
              << "  [-reverse_endian]                   : Raw 16-bit file with reversed endian\n"
              << "  [-psd <assign_CTF_prm_file>]        : Show the PSDs\n"
              << "  [-ctf <assign_CTF_prm_file>]        : Show the CTF models\n"
              << "  [-auto <model rootname>]            : For autoselection\n"
              << "  [-autoSelect]                       : Autoselect without user interaction\n"
              << "  [-thr <p=1>]                        : Number of threads for automatic picking\n"
    ;
}
