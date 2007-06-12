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

#include "module.h"

#include <qgroupbox.h>
#include <qhbox.h>
#include <qfiledialog.h>

/* Empty constructor ------------------------------------------------------- */
CTFViewer::CTFViewer(QWidget *parent, const char *name,
                     const FileName &fn_ctfparam): QMainWindow(parent, name)
{
    // Initialize variables
    firstSettings = plotter.zoomStack[plotter.curZoom];
    ctf.clear();
    ctf_valid = false;
    psdPresent = false;
    
    // Connect plotter resize to recompute curves
    connect(&plotter, SIGNAL(resizeDone()),
            this, SLOT(recomputeCurves()));

    // Initialize geometry
    this->setGeometry(0, 0, 280, 500);
    this->setMinimumSize(280, 500);
    this->setMaximumSize(790, 500);

    // Set parameter window layout
    // Group A
    QGroupBox *GroupA = new QGroupBox("Delta fU", this);
    GroupA->setGeometry(10, 10, 130, 50);
    GroupA->show();

    QHBox *hboxA = new QHBox(GroupA);
    hboxA->setCaption("Enter value for DfU");
    hboxA->setMargin(6);
    hboxA->setSpacing(6);
    hboxA->setGeometry(5, 15, 120, 30);
    hboxA->show();

    lineEditDfU = new QLineEdit(hboxA);
    lineEditDfU->setGeometry(1, 1, 44, 30);
    lineEditDfU->setAlignment(Qt::AlignRight);
    lineEditDfU->setText("");
    lineEditDfU->setEnabled(false);

    connect(lineEditDfU, SIGNAL(returnPressed()),
            this, SLOT(changeMainCTFVariables()));

    // Group B
    QGroupBox *GroupB = new QGroupBox("Delta fV", this);
    GroupB->setGeometry(140, 10, 130, 50);
    GroupB->show();

    QHBox *hboxB = new QHBox(GroupB);
    hboxB->setCaption("Enter value for DfV");
    hboxB->setMargin(6);
    hboxB->setSpacing(6);
    hboxB->setGeometry(5, 15, 120, 30);
    hboxB->show();

    lineEditDfV = new QLineEdit(hboxB);
    lineEditDfV->setAlignment(Qt::AlignRight);
    lineEditDfV->setEnabled(false);
    lineEditDfV->setText("");
    connect(lineEditDfV, SIGNAL(returnPressed()),
            this, SLOT(changeMainCTFVariables()));

    // Group C
    QGroupBox *GroupC = new QGroupBox("Azimuthal Angle", this);
    GroupC->setGeometry(10, 70, 130, 50);
    GroupC->show();

    QHBox *hboxC = new QHBox(GroupC);
    hboxC->setCaption("Enter value for Azimuthal Angle");
    hboxC->setMargin(6);
    hboxC->setSpacing(6);
    hboxC->setGeometry(5, 15, 120, 30);
    hboxC->show();

    lineEditAzi = new QLineEdit(hboxC);
    lineEditAzi->setAlignment(Qt::AlignRight);

    lineEditAzi->setText("");
    lineEditAzi->setEnabled(false);
    connect(lineEditAzi, SIGNAL(returnPressed()),
            this, SLOT(changeMainCTFVariables()));

    // Combobox to select variable
    SelecVar = new QComboBox(this);
    SelecVar->setGeometry(140, 65, 140, 25);
    QFont serifFont("Times", 10);
    SelecVar->setFont(serifFont);
    SelecVar->show();
    SelecVar->insertItem("Sampling Rate", -1);                //Tm
    SelecVar->insertItem("Aceleration voltage", -2);          //kV
    SelecVar->insertItem("Spherical aberration", -3);         //Cs
    SelecVar->insertItem("Chromatic aberration", -4);         //Ca
    SelecVar->insertItem("Energy Loss", -5);                  //espr
    SelecVar->insertItem("Lens Stability", -6);               //ispr
    SelecVar->insertItem("Convergence Cone", -7);             //alpha
    SelecVar->insertItem("Longitudinal Displace", -8);        //DeltaF
    SelecVar->insertItem("Transversal Displace", -9);         //DeltaR
    SelecVar->insertItem("Amplitud Contrast", -10);           //Q0
    SelecVar->insertItem("Global Gain", -11);                 //K
    SelecVar->insertItem("Gain for the gaussian term", -12);  //gaussian_K
    SelecVar->insertItem("Gaussian Width U", -13);            //sigmaU
    SelecVar->insertItem("Gaussian Width V", -14);            //sigmaV
    SelecVar->insertItem("Gaussian Center U", -15);           //cU
    SelecVar->insertItem("Gaussian Center V", -16);           //cV
    SelecVar->insertItem("Gaussian Angle", -17);              //rad_gaussian
    SelecVar->insertItem("Gain for the square root", -18);    //sqrt_K
    SelecVar->insertItem("sqrt Witdh U", -19);                //sqU
    SelecVar->insertItem("sqrt Witdh V", -20);                //sqV
    SelecVar->insertItem("sqrt angle", -21);                  //sqrt_angle
    SelecVar->insertItem("Global Base Line", -22);            //base_line

    connect(SelecVar, SIGNAL(activated(int)), this, SLOT(setVar(int)));

    lineEditVar = new QLineEdit(this);
    lineEditVar->setGeometry(145, 90, 130, 20);

    lineEditVar->setAlignment(Qt::AlignRight);
    lineEditVar->show();
    lineEditVar->setEnabled(false);
    connect(lineEditVar, SIGNAL(returnPressed()), this, SLOT(changeVar()));
    connect(lineEditVar, SIGNAL(lostFocus()), this, SLOT(changeVar()));

    // Plot settings
    QGroupBox *GroupSet = new QGroupBox("Plot Settings", this);
    GroupSet->setGeometry(10, 130, 260, 130);
    GroupSet->show();

    QLabel *ValueMaxX = new QLabel("Max X", GroupSet);
    ValueMaxX->setGeometry(15, 20, 40, 20);
    ValueMaxX->show();

    lineEditMaxX = new QLineEdit(GroupSet);
    lineEditMaxX->setAlignment(Qt::AlignRight);

    lineEditMaxX->setText("");
    lineEditMaxX->setGeometry(55, 20, 60, 20);
    lineEditMaxX->setEnabled(false);
    QObject::connect(lineEditMaxX, SIGNAL(returnPressed()),
                     this, SLOT(setPlotSettings()));
    QObject::connect(lineEditMaxX, SIGNAL(lostFocus()),
                     this, SLOT(setPlotSettings()));

    QLabel *ValueMaxY = new QLabel("Max Y", GroupSet);
    ValueMaxY->setGeometry(125, 20, 40, 20);
    ValueMaxY->show();

    lineEditMaxY = new QLineEdit(GroupSet);
    lineEditMaxY->setAlignment(Qt::AlignRight);
    lineEditMaxY->setText("");
    lineEditMaxY->setGeometry(170, 25, 60, 20);
    lineEditMaxY->setEnabled(false);
    QObject::connect(lineEditMaxY, SIGNAL(returnPressed()),
                     this, SLOT(setPlotSettings()));
    QObject::connect(lineEditMaxY, SIGNAL(lostFocus()),
                     this, SLOT(setPlotSettings()));

    QLabel *ValueXTicks = new QLabel("XTicks", GroupSet);
    ValueXTicks->setGeometry(15, 60, 40, 20);
    ValueXTicks->show();

    lineEditXTicks = new QLineEdit(GroupSet);
    lineEditXTicks->setAlignment(Qt::AlignRight);
    lineEditXTicks->setText("");
    lineEditXTicks->setGeometry(55, 60, 60, 20);
    lineEditXTicks->setEnabled(false);
    QObject::connect(lineEditXTicks, SIGNAL(returnPressed()),
                     this, SLOT(setPlotSettings()));
    QObject::connect(lineEditXTicks, SIGNAL(lostFocus()),
                     this, SLOT(setPlotSettings()));

    QLabel *ValueYTicks = new QLabel("YTicks", GroupSet);
    ValueYTicks->setGeometry(125, 60, 40, 20);
    ValueYTicks->show();

    lineEditYTicks = new QLineEdit(GroupSet);
    lineEditYTicks->setAlignment(Qt::AlignRight);
    lineEditYTicks->setText("");
    lineEditYTicks->setGeometry(170, 60, 60, 20);
    lineEditYTicks->setEnabled(false);
    QObject::connect(lineEditYTicks, SIGNAL(returnPressed()),
                     this, SLOT(setPlotSettings()));
    QObject::connect(lineEditYTicks, SIGNAL(lostFocus()),
                     this, SLOT(setPlotSettings()));

    QLabel *ValueMinX = new QLabel("MinX", GroupSet);
    ValueMinX->setGeometry(15, 100, 40, 20);
    ValueMinX->show();

    lineEditMinX = new QLineEdit(GroupSet);
    lineEditMinX->setAlignment(Qt::AlignRight);
    lineEditMinX->setText("");
    lineEditMinX->setGeometry(55, 100, 60, 20);
    lineEditMinX->setEnabled(false);
    QObject::connect(lineEditMinX, SIGNAL(returnPressed()),
                     this, SLOT(setPlotSettings()));
    QObject::connect(lineEditMinX, SIGNAL(lostFocus()),
                     this, SLOT(setPlotSettings()));

    QLabel *ValueMinY = new QLabel("MinY", GroupSet);
    ValueMinY->setGeometry(125, 100, 40, 20);
    ValueMinY->show();

    lineEditMinY = new QLineEdit(GroupSet);
    lineEditMinY->setAlignment(Qt::AlignRight);
    lineEditMinY->setText("");
    lineEditMinY->setGeometry(170, 100, 60, 20);
    lineEditMinY->setEnabled(false);
    QObject::connect(lineEditMinY, SIGNAL(returnPressed()),
                     this, SLOT(setPlotSettings()));
    QObject::connect(lineEditMinY, SIGNAL(lostFocus()),
                     this, SLOT(setPlotSettings()));

    // Radio buttons
    QGroupBox *GroupRadio = new QGroupBox("View Options", this);
    GroupRadio->setGeometry(10, 270, 260, 110);
    GroupRadio->show();

    ShowDamping = new QRadioButton("Show Damping", GroupRadio);
    ShowDamping->setGeometry(10, 15, 200, 20);
    ShowDamping->show();
    ShowDamping->setChecked(true);
    ShowDamping->setEnabled(false);
    connect(ShowDamping, SIGNAL(clicked()), this, SLOT(refreshCurves()));

    ShowBackground = new QRadioButton("Show Background (Noise)", GroupRadio);
    ShowBackground->setGeometry(10, 35, 200, 20);
    ShowBackground->show();
    ShowBackground->setChecked(true);
    ShowBackground->setEnabled(false);
    connect(ShowBackground, SIGNAL(clicked()), this, SLOT(refreshCurves()));

    ShowCTF = new QRadioButton("Show CTF", GroupRadio);
    ShowCTF->setGeometry(10, 55, 100, 20);
    ShowCTF->show();
    ShowCTF->setChecked(true);
    ShowCTF->setEnabled(false);
    connect(ShowCTF, SIGNAL(clicked()), this, SLOT(refreshCurves()));

    QLabel *ValueAngle = new QLabel("Angle (In Degrees)", GroupRadio);
    ValueAngle->setGeometry(60, 80, 110, 20);
    ValueAngle->show();

    spinBoxAngle = new QSpinBox(GroupRadio);
    spinBoxAngle->setRange(-90, 90);
    spinBoxAngle->setLineStep(1);
    spinBoxAngle->setValue(0);
    spinBoxAngle->setGeometry(10, 80, 40, 18);
    spinBoxAngle->show();
    spinBoxAngle->setEnabled(false);
    connect(spinBoxAngle, SIGNAL(valueChanged(int)), this, SLOT(drawAngle()));

    TenLog = new QRadioButton("10 Log (x)", GroupRadio);
    TenLog->setGeometry(155, 15, 20, 20);
    TenLog->show();
    TenLog->setChecked(false);
    TenLog->setEnabled(false);
    connect(TenLog, SIGNAL(clicked()), this, SLOT(recomputeCurves()));

    QLabel *LabTenLog = new QLabel("10Log10(x)", GroupRadio);
    LabTenLog->setGeometry(175, 10, 70, 30);
    LabTenLog->show();

    ShowExperimental = new QRadioButton("Show Experimental", GroupRadio);
    ShowExperimental->setGeometry(110, 55, 125, 20);
    ShowExperimental->show();
    ShowExperimental->setChecked(true);
    ShowExperimental->setEnabled(false);
    connect(ShowExperimental, SIGNAL(clicked()), this, SLOT(refreshCurves()));

    // Select file
    QGroupBox *InputFile = new QGroupBox("Input file", this);
    InputFile->setGeometry(10, 390, 260, 70);
    InputFile->show();

    getFile = new QPushButton("&Browse", InputFile);
    getFile->setGeometry(195, 14, 50, 25);
    getFile->show();
    connect(getFile, SIGNAL(clicked()), this, SLOT(selectFile()));

    selectedFile = new QLineEdit("", InputFile);
    selectedFile->setGeometry(10, 15, 180, 20);
    selectedFile->show();

    // Quit
    quitButton = new QPushButton("&Quit", this);
    connect(quitButton, SIGNAL(clicked()), qApp, SLOT(quit()));
    quitButton->setGeometry(185, 465, 90, 31);
    quitButton->show();

    // Photo
    photoButton = new QPushButton("&Photo", this);
    connect(photoButton, SIGNAL(clicked()), this, SLOT(photo()));
    photoButton->setGeometry(95, 465, 90, 31);
    photoButton->show();
    photoButton->setEnabled(false);

    // Reset
    PBreset = new QPushButton("&Reset", this);
    connect(PBreset, SIGNAL(clicked()), this, SLOT(reset()));
    PBreset->setGeometry(5, 465, 90, 31);
    PBreset->show();
    PBreset->setEnabled(false);
    show();
    setCaption("Options");
    if (fn_ctfparam != "")
    {
        selectedFile->setText(fn_ctfparam.c_str());
        fn_root = fn_ctfparam.without_extension();
        openFile();
    }
}

// Take a photo ------------------------------------------------------------
void CTFViewer::photo()
{
    Plotter *foto = new Plotter();
    foto->copy(plotter);
    foto->setCaption("Photo");
    foto->show();
}

// Reset plot settings -----------------------------------------------------
void CTFViewer::reset()
{
    plotter.setPlotSettings(firstSettings);
    lineEditMaxX->setText(QString::number(firstSettings.maxX, 'f', 3));
    lineEditMaxY->setText(QString::number(firstSettings.maxY, 'f', 3));
    lineEditMinX->setText(QString::number(firstSettings.minX, 'f', 3));
    lineEditMinY->setText(QString::number(firstSettings.minY, 'f', 3));
    lineEditYTicks->setText(QString::number(firstSettings.numYTicks));
    lineEditXTicks->setText(QString::number(firstSettings.numXTicks));
    refreshCurves();
}

// Select a file from disk -------------------------------------------------
void CTFViewer::selectFile()
{
    int lastSlashPos;
    int length;
    QString fn_ctfparam = QFileDialog::getOpenFileName(QDir::currentDirPath(),
                          "Parameter files (*.ctfparam)", this);
    if (!fn_ctfparam.isEmpty())
    {
        selectedFile->setText(fn_ctfparam);
        fn_root = (static_cast<FileName>(fn_ctfparam.ascii())).without_extension();
        openFile();
    }
}

// Open file ---------------------------------------------------------------
void CTFViewer::openFile()
{
    // Show CTF-Experimental Image
    setImageViewer();

    // Prepare XmippCTF
    ctf.enable_CTFnoise = true;
    ctf.enable_CTF = true;
    ctf.read(fn_root + ".ctfparam");
    ctf.Produce_Side_Info();
    ctf_valid = true;

    // Recompute the curves shown
    recomputeCurves();

    // Initialize firstSettings
    firstSettings = plotter.zoomStack[plotter.curZoom];

    // Enable all LineEdits and other widgets
    lineEditDfU->setText(QString::number(ctf.DeltafU, 'f', 6));
    lineEditDfU->setEnabled(true);
    lineEditDfV->setText(QString::number(ctf.DeltafV, 'f', 6));
    lineEditDfV->setEnabled(true);
    lineEditAzi->setText(QString::number(ctf.azimuthal_angle, 'f', 6));
    lineEditAzi->setEnabled(true);
    lineEditVar->setText(QString::number(ctf.Tm, 'f', 6));
    lineEditVar->setEnabled(true);
    lineEditMaxX->setText(QString::number(firstSettings.maxX, 'f', 3));
    lineEditMaxX->setEnabled(true);
    lineEditMaxY->setText(QString::number(firstSettings.maxY, 'f', 3));
    lineEditMaxY->setEnabled(true);
    lineEditMinX->setText(QString::number(firstSettings.minX, 'f', 3));
    lineEditMinX->setEnabled(true);
    lineEditMinY->setText(QString::number(firstSettings.minY, 'f', 3));
    lineEditMinY->setEnabled(true);
    lineEditYTicks->setText(QString::number(firstSettings.numYTicks));
    lineEditYTicks->setEnabled(true);
    lineEditXTicks->setText(QString::number(firstSettings.numXTicks));
    lineEditXTicks->setEnabled(true);
    spinBoxAngle->setEnabled(true);
    ShowCTF->setEnabled(true);
    ShowBackground->setEnabled(true);
    ShowDamping->setEnabled(true);
    TenLog->setEnabled(true);
    ShowExperimental->setEnabled(true);
    photoButton->setEnabled(true);
    PBreset->setEnabled(true);
}

// Change one variable of the combobox -------------------------------------
void CTFViewer::changeVar()
{
    bool Ok;
    double val = lineEditVar->text().toDouble(&Ok);
    if (!Ok) return;
    switch (SelecVar->currentItem())
    {
    case  0:
        ctf.Tm = val;
        break;
    case  1:
        ctf.kV = val;
        break;
    case  2:
        ctf.Cs = val;
        break;
    case  3:
        ctf.Ca = val;
        break;
    case  4:
        ctf.espr = val;
        break;
    case  5:
        ctf.ispr = val;
        break;
    case  6:
        ctf.alpha = val;
        break;
    case  7:
        ctf.DeltaF = val;
        break;
    case  8:
        ctf.DeltaR = val;
        break;
    case  9:
        ctf.Q0 = val;
        break;
    case 10:
        ctf.K = val;
        break;
    case 11:
        ctf.gaussian_K = val;
        break;
    case 12:
        ctf.sigmaU = val;
        break;
    case 13:
        ctf.sigmaV = val;
        break;
    case 14:
        ctf.cU = val;
        break;
    case 15:
        ctf.cV = val;
        break;
    case 16:
        ctf.rad_gaussian = val;
        break;
    case 17:
        ctf.sqrt_K = val;
        break;
    case 18:
        ctf.sqU = val;
        break;
    case 19:
        ctf.sqV = val;
        break;
    case 20:
        ctf.sqrt_angle = val;
        break;
    case 21:
        ctf.base_line = val;
        break;
    }
    generate_ctfmodel();
}

void CTFViewer::setVar(int e)
{
    if (!ctf_valid) return;
    switch (e)
    {
    case  0:
        lineEditVar->setText(QString::number(ctf.Tm          , 'f', 6));
        break;
    case  1:
        lineEditVar->setText(QString::number(ctf.kV          , 'f', 6));
        break;
    case  2:
        lineEditVar->setText(QString::number(ctf.Cs          , 'f', 6));
        break;
    case  3:
        lineEditVar->setText(QString::number(ctf.Ca          , 'f', 6));
        break;
    case  4:
        lineEditVar->setText(QString::number(ctf.espr        , 'f', 6));
        break;
    case  5:
        lineEditVar->setText(QString::number(ctf.ispr        , 'f', 6));
        break;
    case  6:
        lineEditVar->setText(QString::number(ctf.alpha       , 'f', 6));
        break;
    case  7:
        lineEditVar->setText(QString::number(ctf.DeltaF      , 'f', 6));
        break;
    case  8:
        lineEditVar->setText(QString::number(ctf.DeltaR      , 'f', 6));
        break;
    case  9:
        lineEditVar->setText(QString::number(ctf.Q0          , 'f', 6));
        break;
    case 10:
        lineEditVar->setText(QString::number(ctf.K           , 'f', 6));
        break;
    case 11:
        lineEditVar->setText(QString::number(ctf.gaussian_K  , 'f', 6));
        break;
    case 12:
        lineEditVar->setText(QString::number(ctf.sigmaU      , 'f', 6));
        break;
    case 13:
        lineEditVar->setText(QString::number(ctf.sigmaV      , 'f', 6));
        break;
    case 14:
        lineEditVar->setText(QString::number(ctf.cU          , 'f', 6));
        break;
    case 15:
        lineEditVar->setText(QString::number(ctf.cV          , 'f', 6));
        break;
    case 16:
        lineEditVar->setText(QString::number(ctf.rad_gaussian, 'f', 6));
        break;
    case 17:
        lineEditVar->setText(QString::number(ctf.sqrt_K      , 'f', 6));
        break;
    case 18:
        lineEditVar->setText(QString::number(ctf.sqU         , 'f', 6));
        break;
    case 19:
        lineEditVar->setText(QString::number(ctf.sqV         , 'f', 6));
        break;
    case 20:
        lineEditVar->setText(QString::number(ctf.sqrt_angle  , 'f', 6));
        break;
    case 21:
        lineEditVar->setText(QString::number(ctf.base_line   , 'f', 6));
        break;
    }
    generate_ctfmodel();
}

// Determine active curves -------------------------------------------------
void CTFViewer::refreshCurves()
{
    plotter.setActiveState(1, ShowCTF->isOn());
    plotter.setActiveState(2, ShowDamping->isOn());
    plotter.setActiveState(3, ShowBackground->isOn());
    if (psdPresent)
        plotter.setActiveState(4, ShowExperimental->isOn());
    plotter.refreshPixmap();
}

// Recompute curves --------------------------------------------------------
void CTFViewer::recomputeCurves()
{
    int angle = spinBoxAngle->value();
    Matrix2D<double> data;
    getCTFcurve("pure", angle, data);
    plotter.setCurveData(1, data);
    getCTFcurve("damping", angle, data);
    plotter.setCurveData(2, data);
    getCTFcurve("background", angle, data);
    plotter.setCurveData(3, data);
    if (psdPresent)
    {
        getExperimentalCurve(angle, data);
        plotter.setCurveData(4, data);
    }
    setPlotSettings();
    refreshCurves();
}

// Draw angle ------------------------------------------------------------
void CTFViewer::drawAngle()
{
    static bool first_time_in_this_function = true;
    if (!ctf_valid) return;
    int angle = spinBoxAngle->value();
    if (!first_time_in_this_function) imageViewer->drawAngle(old_angle);
    imageViewer->drawAngle(angle);
    first_time_in_this_function = false;
    old_angle = angle;

    recomputeCurves();
}

// Parameter change --------------------------------------------------------
void CTFViewer::changeMainCTFVariables()
{
    if (!ctf_valid) return;
    bool Ok;
    double val;
    val = lineEditDfU->text().toDouble(&Ok);
    if (Ok) ctf.DeltafU = val;
    else return;
    val = lineEditDfV->text().toDouble(&Ok);
    if (Ok) ctf.DeltafV = val;
    else return;
    val = lineEditAzi->text().toDouble(&Ok);
    if (Ok) ctf.azimuthal_angle = val;
    else return;
    generate_ctfmodel();
}

// Plot settings -----------------------------------------------------------
void CTFViewer::setPlotSettings()
{
    bool Ok;
    double dval;
    int ival;
    dval = lineEditMinX->text().toDouble(&Ok);
    if (Ok) plotter.zoomStack[plotter.curZoom].minX = dval;
    dval = lineEditMaxX->text().toDouble(&Ok);
    if (Ok) plotter.zoomStack[plotter.curZoom].maxX = dval;
    dval = lineEditMinY->text().toDouble(&Ok);
    if (Ok) plotter.zoomStack[plotter.curZoom].minY = dval;
    dval = lineEditMaxY->text().toDouble(&Ok);
    if (Ok) plotter.zoomStack[plotter.curZoom].maxY = dval;
    ival = lineEditXTicks->text().toInt(&Ok);
    if (Ok) plotter.zoomStack[plotter.curZoom].numXTicks = ival;
    ival = lineEditYTicks->text().toInt(&Ok);
    if (Ok) plotter.zoomStack[plotter.curZoom].numYTicks = ival;
    refreshCurves();
}

// Get CTF curve -----------------------------------------------------------
void CTFViewer::getCTFcurve(const string &type, int angle,
                            Matrix2D<double> &data, int Nsamples)
{
    double sampling_rate = ctf.Tm; // in Angstroms/pixel
    double fmax = 1.0 / (2.0 * sampling_rate);

    data.resize(Nsamples, 2);
    int i = 0;
    double angle_rad = DEG2RAD(angle);
    for (double f = 0; f < fmax; f += fmax / Nsamples)
    {
        data(i, 0) = f;
        double ctf_bg  = ctf.CTFnoise_at(f * cos(angle_rad), f * sin(angle_rad));
        double ctf_pure = ctf.CTFpure_at(f * cos(angle_rad), f * sin(angle_rad));
        double ctf_E   = ctf.CTFdamping_at(f * cos(angle_rad), f * sin(angle_rad));
        if (type == "pure")
        {
            data(i, 1) = ctf_bg + ctf_pure * ctf_pure;
        }
        else if (type == "damping")
        {
            data(i, 1) = ctf_bg + ctf_E * ctf_E;
        }
        else if (type == "background")
        {
            data(i, 1) = ctf_bg;
        }

        // Check if logarithm
        if (TenLog->isOn()) data(i, 1) = 10 * log10(data(i, 1));
        i++;
        if (i == Nsamples) break;
    }
}

// Set Experimental image --------------------------------------------------
void CTFViewer::setImageViewer()
{
    // Show the CTFmodel
    imageViewer = new ImageViewer((fn_root + ".ctfmodel_halfplane").c_str() , false);
    imageViewer->apply_geo = false;
    imageViewer->loadImage((fn_root + ".ctfmodel_halfplane").c_str(), 0, 0,
                           ImageViewer::CTF_mode);
    imageViewer->show();

    // Load the CTF in image I
    try {
       I.read(fn_root + ".psd");
       CenterFFT(I(), true);
       I().setXmippOrigin();
       psdPresent = true;

       // Recompute curves since the experimental curve has changed
       recomputeCurves();
    } catch (Xmipp_error XE) {
       I.clear();
       psdPresent = false;
       ShowExperimental->setChecked(true);
       ShowExperimental->setEnabled(false);
    }
}

// Get experimental curve --------------------------------------------------
void CTFViewer::getExperimentalCurve(int angle, Matrix2D<double> &data,
                                     int Nsamples)
{
    if (!psdPresent) return;

    double angle_rad = DEG2RAD(angle + 180.0);
    double cos_ang = cos(angle_rad);
    double sin_ang = sin(angle_rad);
    double t_max = XSIZE(I()) / 2;
    double Ts = t_max / Nsamples;
    double sampling_rate = ctf.Tm;

    int Nmax = CEIL(t_max / Ts);
    data.resize(Nmax, 2);
     int i = 0;
     for (double t = 0; t < t_max; t += Ts, i++)
     {
         data(i, 0) = t * 1.0 / (XSIZE(I()) * sampling_rate);
         data(i, 1) = I().interpolatedElement(t * cos_ang, t * sin_ang);
         if (TenLog->isOn()) data(i, 1) = 10 * log10(data(i, 1));
         if (i == Nmax) break;
     }
}

// Generate CTF model ------------------------------------------------------
void CTFViewer::generate_ctfmodel()
{
    ctf.Produce_Side_Info();

    // Get the image size from the .ctfmodel
    ImageXmipp model;
    model.read(fn_root + ".ctfmodel_halfplane");
    int Ydim = YSIZE(model());
    int Xdim = XSIZE(model());

    // Substitute the left part by the CTF model
    Matrix1D<double> freq(2); // Frequencies for Fourier plane
    double minval=1e38, maxval=-1e38;
    FOR_ALL_ELEMENTS_IN_MATRIX2D(model())
    {
        if (j < Xdim / 2) continue;
        XX(freq) = ((double)j-Xdim/2.0)/Xdim;
        YY(freq) = ((double)i-Ydim/2.0)/Xdim;
        digfreq2contfreq(freq, freq, ctf.Tm);
        model(i, j) = ctf.CTFpure_at(XX(freq), YY(freq));
	model(i,j)*=model(i,j);
    	minval=MIN(minval,model(i,j));
    	maxval=MAX(maxval,model(i,j));
    }
    FOR_ALL_ELEMENTS_IN_MATRIX2D(model())
    {
        if (j < Xdim / 2) continue;
	model(i,j)=(model(i,j)-minval)/(maxval-minval);
    }

    // Recompute curves
    recomputeCurves();

    // Reload and redraw image
    imageViewer->loadMatrix(model(), 0, 0, ImageViewer::CTF_mode);
}

// Key press events --------------------------------------------------------
void CTFViewer::keyPressEvent(QKeyEvent *event)
{
    switch (event->key())
    {
    case Key_Q:
        if (event->state() == ControlButton) // If 'Ctrol Q' key,
            exit(0);
        break;
    }
}
