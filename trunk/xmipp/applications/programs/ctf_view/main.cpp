/***************************************************************************
 *
 * Authors:     Javier Rodr�guez Fern�ndez (javrodri@gmail.com)
 *              Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Universidad San Pablo CEU (Monteprincipe, Madrid)
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

#include <qapplication.h>

#include <data/args.h>

#include "module.h"

int main(int argc, char *argv[])
{
    FileName fn_ctf;
    try
    {
        fn_ctf = getParameter(argc, argv, "-i", "");
    }
    catch (Xmipp_error XE)
    {
        std::cerr << XE;
        std::cerr << "Usage: CTFViewer\n"
        << "   [-i <CTF file>] : The file is assumed to be of kind .ctfparam\n";
        exit(1);
    }
    try
    {
        QApplication app(argc, argv);
        CTFViewer * Viewer = new CTFViewer(0, 0, fn_ctf);
        app.setMainWidget(Viewer);
        return app.exec();
    }
    catch (Xmipp_error XE)
    {
        std::cerr << XE;
        exit(1);
    }
    return 0;
}
