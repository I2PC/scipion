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
#include <qscrollbar.h>
#include <qimage.h>
#include <qpixmap.h>
#include <qradiobutton.h>
#include <qtextedit.h>

#include <data/image.h>

#include <vector>
#include <string>

/**@defgroup ShowTools ShowTools
   @ingroup GraphicsLibrary */
//@{

/**Scroll param class.
    This class opens a window for asking for the Scroll parameter.
    It emits a signal called new_value(float). (precision is the number of
    digits used by the mantise)

    An example of use of this class with a single parameter is
    @code
      // Create window
      ScrollParam* param_window;
      param_window = new ScrollParam(min, max, spacing, "Set spacing", "Spacing",
         0, "new window", WDestructiveClose);

      // Connect its output to my input (set_spacing)
      connect( param_window, SIGNAL(new_value(float)),
               this,         SLOT(set_spacing(float)));

      // Show
      param_window->setFixedSize(250,150);
      param_window->show();
    @endcode

    With two parameters
    @code
      // Create window
      ScrollParam* param_window;
      vector<float> min; min.push_back(1); min.push_back(1);
      vector<float> max; max.push_back(N); max.push_back(N);
      vector<float> initial_value;
         initial_value.push_back(spacing);
         initial_value.push_back(x_tick_off);
      vector<char *> prm_name;
         prm_name.push_back("Spacing");
         prm_name.push_back("Tick offset");
      param_window = new ScrollParam(min,max,initial_value,prm_name,
  "Set spacing", 0, "new window", WDestructiveClose,0);

      // Connect its output to my input (set_spacing)
      connect( param_window, SIGNAL(new_value(vector<float>)),
              this,          SLOT(set_spacing(vector<float>)) );

      // Show
      param_window->setFixedSize(200,175);
      param_window->show();
    @endcode
*/
class ScrollParam : public QWidget
{
    Q_OBJECT
public:
    /** Constructor for a single scroll.
        Provide the min_value, max_value, caption and initial_value.*/
    ScrollParam(float min, float max, float initial_value, char *prm_name,
                char *caption, QWidget *parent = 0,
                const char *name = 0, int wFlags = 0, int precision = 2);

    /** Constructor for several scrolls.
        Provide the min_value, max_value, caption and initial_value.*/
    ScrollParam(vector<float> &min, vector<float> &max,
                vector<float> &initial_value, vector<char *> &prm_name,
                char *caption, QWidget *parent = 0,
                const char *name = 0, int wFlags = 0, int precision = 2);

    /** Destructor. */
    ~ScrollParam();

    /** Init. */
    void init(vector<float> &min, vector<float> &max,
              vector<float> &initial_value, vector<char *> &prm_name,
              char *caption, int precision = 2);

    /** Get current values. */
    vector<float> getCurrentValues();
private:
    vector<float>     value;
    vector<QLabel *>  value_lab;   // label for the current value of the slider
    vector<QScrollBar *> scroll;   // sliders
    int       my_precision;
private slots:
    void scrollValueChanged(int);
    void slot_close_clicked();
    void slot_ok_clicked();
signals:
    /** Signal emitted when the value is changed*/
    void new_value(float);
    /** Signal emitted when the value is changed*/
    void new_value(vector<float>);
    /** Signal emitted when the close button is clicked */
    void signal_close_clicked();
    /** Signal emitted when the ok button is clicked */
    void signal_ok_clicked();
};

/**Exclusive param class.
    This class opens a window for asking for a exclusive parameter.
    It emits a signal called new_value(int) with the selected value

    An example of use of this class is
    @code
       vector<string> list_values;
       list_values.push_back("Option 1");
       list_values.push_back("Option 2");
      // Create window
      ExclusiveParam* param_window=
           new ExclusiveParam(list_values, parameter, "Set this exclusive parameter",
             0, "new window", WDestructiveClose);

      // Connect its output to my input (set_spacing)
      connect( param_window, SIGNAL(new_value(int)),
               this,         SLOT(set_spacing(int)));

      // Show
      param_window->setFixedSize(250,200);
      param_window->show();
    @endcode
*/
class ExclusiveParam : public QWidget
{
    Q_OBJECT
public:
    /** Constructor.
        Provide the min_value, max_value, caption and initial_value.*/
    ExclusiveParam(vector<string> &list_values, int initial_value,
                   char *caption, QWidget *parent = 0,
                   const char *name = 0, int wFlags = 0);
    ~ExclusiveParam();

private:
    int       value;
    vector< QRadioButton *> button;
private slots:
    void but_close_clicked();
    void but_ok_clicked();
    void exclusiveValueChanged();
signals:
    /** Signal emitted when the value is changed*/
    void new_value(int);
};

/** Image conversions */
//@{
/** Xmipp -> QImage.*/
void xmipp2Qt(Image& _ximage, QImage &_qimage,
              int _minScale = 0, int _maxScale = 255, double _m = 0, double _M = 0);

/** Qimage -> Xmipp.*/
void Qt2xmipp(QImage &_qimage, Image &_ximage);

/** Xmipp -> PixMap */
void xmipp2Pixmap(Image &xmippImage, QPixmap* pixmap,
                  int _minScale = 0, int _maxScale = 255, double _m = 0, double _M = 0);

/** Xmipp image -> Xmipp PSD.
    The log10 is taken, outliers rejected and the image is reorganized. */
void xmipp2PSD(const Matrix2D<double> &input, Matrix2D<double> &output);

/** Xmipp image -> Xmipp CTF.
    The log10 is taken, outliers rejected and the image is reorganized. */
void xmipp2CTF(const Matrix2D<double> &input, Matrix2D<double> &output);
//@}

/** Miscellanea */
//@{
/** Return a QPixmap with the file provided.
    Function taken from qt/src/kernel/qpixmap.cpp for compatibility reasons. */
QPixmap xmipp_qPixmapFromMimeSource(const QString &abs_name);
//@}

//@}
#endif
