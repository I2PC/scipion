/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
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

#include "maskimg.hh"
#include <qapplication.h>
#include <qimage.h>
#include <qwindowdefs.h>
#include <qnamespace.h>

#ifdef QIMGIO
#include <qimageio.h>
#endif

#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippArgs.hh>


int main( int argc, char **argv )
{
    QApplication::setFont( QFont("Helvetica", 12) );
    QApplication a( argc, argv );

#ifdef QIMGIO
    qInitImageIO();
#endif

  string selname = "", imgname = "";
  bool sdflag = false;

  try {     
    if (check_param(argc, argv, "-sel"))
    	selname = get_param(argc, argv, "-sel");
    else 
    	imgname = get_param(argc, argv, "-img");
    if (check_param(argc, argv, "-sd"))
     sdflag = true; 
  } 
  catch (Xmipp_error) {
    cout << "Xmask: Creates a mask using a Graphical User Interface" << endl;
    cout << "Usage:" << endl;
    cout << "-img           : Image to visualize" << endl;
    cout << "-sel           : Use Sel file average or SD Image" << endl;
    cout << "[-sd]          : Uses SD image instead of Average image (default: false)" << endl;
    exit(1);
   }

   if (imgname != "")  {
	    maskImg *w = new maskImg(0, imgname.c_str(), QWidget::WDestructiveClose);
	    w->loadImage( imgname.c_str() );
	    w->show();
   } else {
        cout << "Calculating average and SD images from sel file......" << endl;
  	SelFile SF((FileName) selname); 
	Image ave, sd;
	double min, max;
	SF.get_statistics(ave, sd, min, max); 
	maskImg* w;
	if (sdflag) 
	   w = new maskImg(NULL, &sd, CIRCLE, selname.c_str(), QWidget::WDestructiveClose);
	else    
	   w = new maskImg(NULL, &ave, CIRCLE, selname.c_str(), QWidget::WDestructiveClose);
	w->show();
   }
   QObject::connect(qApp, SIGNAL(lastWindowClosed()), qApp, SLOT(quit()));
   return a.exec();
}
