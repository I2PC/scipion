/***************************************************************************
 *
 * Authors:      Alberto Pascual Montano (pascual@cnb.uam.es)
 *               Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#ifndef _SHOWTOOLS_H
#define _SHOWTOOLS_H
#include <qwidget.h>
#include <qlabel.h>
#include <qimage.h>
#include <qpixmap.h>
#include <XmippData/xmippImages.hh>

/**@name ShowTools */
//@{

/**Scroll param class.
    This class opens a window for asking for the Scroll parameter.
    It emits a signal called new_value(float). (precision is the number of 
    digits used by the mantise)
    
    An example of use of this class is
    \begin{verbatim}
      // Create window
      ScrollParam* param_window=
           new ScrollParam(min,max, spacing, "Set spacing", 
             0, "new window", WDestructiveClose, precision);

      // Connect its output to my input (set_spacing)
      connect( param_window, SIGNAL(new_value(float)), 
               this,         SLOT(set_spacing(float)));

      // Show
      param_window->setFixedSize(250,200);
      param_window->show();
    \end{verbatim}
*/
class ScrollParam : public QWidget
{ 
  Q_OBJECT
public:
    /** Constructor.
        Provide the min_value, max_value, caption and initial_value.*/
    ScrollParam( float min, float max, float initial_value, 
       char *caption, QWidget *parent = 0,
       const char *name = 0, int wFlags=0, int precision=2  );
    ~ScrollParam();

private:
   float     value;
   QLabel   *value_lab;   //test label for the current value of the slider
   int       my_precision;
private slots:
   void scrollValueChanged(int);
   void but_close_clicked();
   void but_ok_clicked();   
signals:
   /** Signal emitted when the value is changed*/
   void new_value(float);
};


/**Scroll param2 class.
    This class opens a window and asks for TWO Scroll parameters.
    It emits a signal called new_value(float,float). (precision is the number of 
    digits used by the mantise)
    
    An example of use of this class is
    \begin{verbatim}
      // Create window
      ScrollParam* param_window=
           new ScrollParam(min,max, first init parameter, 
	      second_init_parameter,"Set spac/x_ticks", 
             0, "new window", WDestructiveClose, precision);

      // Connect its output to my input (set_spacing)
      connect( param_window, SIGNAL(new_value(float,float)), 
               this,         SLOT(set_spacing(float,float)));

      // Show
      param_window->setFixedSize(250,200);
      param_window->show();
    \end{verbatim}
*/
class ScrollParam2 : public QWidget
{ 
  Q_OBJECT
public:
    /** Constructor.
        Provide the min_value, max_value, caption and initial_value.*/
    ScrollParam2( float min, float max, float initial_value, 
       float initial_value2,
       char *caption, QWidget *parent = 0,
       const char *name = 0, int wFlags=0, int precision=2  );
    ~ScrollParam2();

private:
   float     value;
   float     value2;
   QLabel   *value_lab;   //test label for the current value of the slider
   QLabel   *value_lab2;   //test label for the current value of the slider
   int       my_precision;
private slots:
   void scrollValueChanged(int);
   void scrollValueChanged2(int);
   void but_close_clicked();
   void but_ok_clicked();   
signals:
   /** Signal emitted when the value is changed*/
   void new_value(float,float);
};

/**@name Image conversions */
//@{
/** QImage -> Xmipp.*/
void xmipp2Qt(Image& _ximage, QImage &_qimage,
   int _minScale = 0, int _maxScale = 255);

/** Xmipp -> QImage.*/
void Qt2xmipp(QImage &_qimage, Image &_ximage);

/** Xmipp -> PixMap */
void xmipp2Pixmap(Image &xmippImage, QPixmap* pixmap,
   int _minScale = 0, int _maxScale = 255);
//@}

//@}
#endif
