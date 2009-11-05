/***************************************************************************
 *
 * Authors:     Javier Rodríguez Fernández (javrodri@gmail.com)
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
/* ------------------------------------------------------------------------- */
/* CTFVIEWER                                                                  */
/* ------------------------------------------------------------------------- */

#ifndef CTFVIEWER_H
#define CTFVIEWER_H

/* ************************************************************************* */
/* INCLUDES                                                                  */
/* ************************************************************************* */

#include <graphics/show_plotter.h>
#include <graphics/show_2d.h>
#include <data/ctf.h>

#include <qcombobox.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qpushbutton.h>
#include <qradiobutton.h>
#include <qspinbox.h>

#ifdef QT3_SUPPORT
//Added by qt3to4:
#include <QKeyEvent>
#endif

/* ************************************************************************* */
/* CLASS DEFINITION AND PROTOTYPES                                           */
/* ************************************************************************* */
/* Class to handle Plotter to show CTF curves (CTF,Damping,Noise) and
    Experimental Curve, it also handles the Plotter settings as MaxX, MinX,
    numXTicks,etc and allows the user to load files, make photos of the plotter,
    change the variables that define the CTF and select  the angle that defines
    the section it has represented in the plotter. */
#ifdef QT3_SUPPORT
// MOC_SKIP_BEGIN
class CTFViewer: public Q3MainWindow
// MOC_SKIP_END
#else
class CTFViewer : public QMainWindow
#endif
{
    // The Q_OBJECT macro is necessary for all clases that defines signals or
    // slots.
    Q_OBJECT
    // External interface
public:
    // Constructor for CTFViewer. The fn_ctf filename may be empty
    CTFViewer(QWidget *parent, const char *name,
              const FileName &fn_ctfparam = "");

    // Internal data
public:
    //Variable to handle the plotter
    Plotter plotter;
    // CTF
    XmippCTF ctf;
    // CTF valid
    bool ctf_valid;
    // Image model shown in the imageViewer
    ImageXmipp I;
    // Widget to show the Image
    ImageViewer *imageViewer;
    // Filename root of the model being visualized
    FileName fn_root;
    // PSD present
    bool psdPresent;

    // First settings setted when loaded the File,
    // used when button reset is pressed
    PlotSettings firstSettings;

    // Radio Buttons of show options
    QRadioButton *ShowCTF;
    QRadioButton *ShowBackground;
    QRadioButton *ShowDamping;
    QRadioButton *TenLog;
    QRadioButton *ShowExperimental;
    QLabel *LabelTenLog;

    // LineEdit for plotsettings
    QLineEdit *lineEditMaxX;
    QLineEdit *lineEditMaxY;
    QLineEdit *lineEditXTicks;
    QLineEdit *lineEditYTicks;
    QLineEdit *lineEditMinX;
    QLineEdit *lineEditMinY;

    // LineEdits for main CTF values
    QLineEdit *lineEditDfU;
    QLineEdit *lineEditDfV;
    QLineEdit *lineEditAzi;

    // Combo and lineEdit to change CTF values
    QComboBox *SelecVar;
    QLineEdit *lineEditVar;

    // Different buttons
    QPushButton *PBreset;
    QPushButton *quitButton;
    QPushButton *photoButton;

    //Spinbox to select section angle
    QSpinBox *spinBoxAngle;
    int old_angle;

    //Objects to select and open a file
    QPushButton *getFile;
    QLineEdit *selectedFile;

    // Internal functions
private:
    // Function to create a Curve for CTF
    // Valid types are "pure", "damping", "noise"
    void getCTFcurve(const std::string &type, int angle,
                     Matrix2D<double> &curve, int Nsamples = 200);
    // Function to create the curve for the experimental curve
    void getExperimentalCurve(int angle, Matrix2D<double> &curve,
                              int Nsamples = 200);
    // Function to make the CTFModel of the image,
    // used when any Variable change. It regenerates the curves
    void generate_ctfmodel();

    // Function used to open a Widget that sets the image corresponding
    // to the file loaded in CTFViewer
    void setImageViewer();

public:
    // Key handler
    void keyPressEvent(QKeyEvent *event);

    // Slots
private slots:
    // Slots to change Plotsettings, connected to the plotter via plotter
    // pointer
    void setPlotSettings();
    // Slot to open a file
    void openFile();
    // Slot to select the file to open
    void selectFile();
    // Refresh curves
    void refreshCurves();
    // Recompute curves and show them
    void recomputeCurves();
    // Draw the angle on the CTFmodel
    void drawAngle();
    // Slot to change not main CTF variables
    void setVar(int);
    // Slot to set the var selected in the LineEdit
    void changeVar();
    // Slots to change main CTF variables
    void changeMainCTFVariables();
    // Slot to reset the plotter
    void reset();
    // Slot to take a photo of the plotter
    void photo();
};
#endif
