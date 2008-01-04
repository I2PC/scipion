/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Carlos Manzanares       (cmanzana@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include <data/args.h>

#include <qapplication.h>

#include "main_widget_mark.h"
#include "widget_psd.h"

void Usage();

int main(int argc, char **argv)
{
    FileName fnRaw;
    FileName fnRawTilted;
    bool     reversed;
    FileName fn_assign_CTF;
    bool     ctf_mode = false;

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
        FileName fn8bits, fn8bitsTilted;

        m.open_micrograph(fnRaw, reversed);
        m.compute_8_bit_scaling();
        if (m.depth()!=8)
        {
            fn8bits=fnRaw+".8bits";
            m.write_as_8_bits(fn8bits);
            m.close_micrograph();
            m.open_micrograph(fn8bits,false);
            m.compute_8_bit_scaling();
        }
        
        if (fnRawTilted != "")
        {
            mTilted.open_micrograph(fnRawTilted, reversed);
            mTilted.compute_8_bit_scaling();
            if (mTilted.depth()!=8)
            {
                fn8bitsTilted=fnRawTilted+".8bits";
                mTilted.write_as_8_bits(fn8bitsTilted);
                mTilted.close_micrograph();
                mTilted.open_micrograph(fn8bitsTilted,false);
                mTilted.compute_8_bit_scaling();
            }
        }

        // Configure application .............................................
        QApplication app(argc, argv);
        QtMainWidgetMark *mainWidget;

        if (fnRawTilted == "") mainWidget = new QtMainWidgetMark(&m);
        else mainWidget = new QtMainWidgetMark(&m, &mTilted);

        // Check if the PSDs must be shown ...................................
        if (fn_assign_CTF != "")
        {
            QtWidgetPSD PSDshow;
            if (ctf_mode) PSDshow.set_CTF_mode();
            PSDshow.set_assign_CTF_file(m, fn_assign_CTF);
            PSDshow.show();
        }

        // Run application ...................................................
        app.setMainWidget(mainWidget);
        mainWidget->show();

        app.exec();

        // Finish ............................................................
        m.close_micrograph();
        if (fnRawTilted != "") mTilted.close_micrograph();
        delete mainWidget;
        if (fn8bits!="") system(((std::string)"rm -rf "+fn8bits+"*").c_str());
        if (fn8bitsTilted!="") system(((std::string)"rm -rf "+fn8bitsTilted+"*").c_str());
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
    ;
}
