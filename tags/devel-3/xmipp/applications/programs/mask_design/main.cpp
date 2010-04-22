/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.csic.es)
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

#include "maskimg.h"

#include <qapplication.h>
#include <qimage.h>
#include <qwindowdefs.h>
#include <qnamespace.h>

#ifdef QIMGIO
#include <qimageio.h>
#endif

#include <data/funcs.h>
#include <data/args.h>

int main(int argc, char **argv)
{
    QApplication::setFont(QFont("Helvetica", 12));
    QApplication a(argc, argv);

#ifdef QIMGIO
    qInitImageIO();
#endif

    std::string selname = "", imgname = "", saveasname = "";
    bool sdflag = false;
    bool apply_geo;

    try
    {
        if (checkParameter(argc, argv, "-sel"))
            selname = getParameter(argc, argv, "-sel");
        else
            imgname = getParameter(argc, argv, "-img");
        saveasname = getParameter(argc, argv, "-save_as", "");
        if (checkParameter(argc, argv, "-sd"))
            sdflag = true;
        apply_geo = !checkParameter(argc, argv, "-dont_apply_geo");
    }
    catch (Xmipp_error)
    {
        std::cout << "Xmask: Creates a mask using a Graphical User Interface" << std::endl;
        std::cout << "Usage:" << std::endl;
        std::cout << "-img              : Image to visualize" << std::endl;
        std::cout << "-sel              : Use Sel file average or SD Image" << std::endl;
        std::cout << "[-sd]             : Uses SD image instead of Average image (default: false)" << std::endl;
        std::cout << "[-save_as <name>] : Always save mask with this name" << std::endl;
        std::cout << "[-dont_apply_geo] : Do not apply transformation stored in the header" << std::endl;

        exit(1);
    }

    if (imgname != "")
    {
        maskImg *w = new maskImg(0, imgname.c_str(), Qt::WDestructiveClose);
        w->apply_geo = apply_geo;
        w->saveasname = saveasname;
        w->loadImage(imgname.c_str());
        w->show();
    }
    else
    {
        std::cout << "Calculating average and SD images from sel file......" << std::endl;
        SelFile SF((FileName) selname);
        Image<double> ave, sd;
        double min, max;
        SF.get_statistics(ave, sd, min, max, apply_geo);
        maskImg* w;
        if (sdflag)
            w = new maskImg(NULL, &sd, CIRCLE, selname.c_str(), Qt::WDestructiveClose);
        else
            w = new maskImg(NULL, &ave, CIRCLE, selname.c_str(), Qt::WDestructiveClose);
        w->apply_geo = apply_geo;
        w->saveasname = saveasname;
        w->show();
    }
    QObject::connect(qApp, SIGNAL(lastWindowClosed()), qApp, SLOT(quit()));
    return a.exec();
}
