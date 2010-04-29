/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Alberto Pascual (pascual@cnb.csic.es)
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

#include "show_spectra.h"
#include "show_tools.h"

#include <qcolordialog.h>
#include <qfontdialog.h>
#include <qmessagebox.h>
#include <qlayout.h>
#include <qpushbutton.h>
#include <qtooltip.h>

#ifdef QT3_SUPPORT
//Added by qt3to4:
#include <QLabel>
#include <Q3GridLayout>
#include <QPixmap>
#include <Q3PopupMenu>
#include <Q3VBoxLayout>
#endif

/* Init/Clear data --------------------------------------------------------- */
void ShowSpectra::init()
{
    V = NULL;
    offX = 10;
    offY = 10;
    spacing = 3;
    x_tick_off = 1;
    backColor = Qt::black;
    axisColor = Qt::white;
    curveColor = Qt::white;
    QFont tmpFont("Clean", 6);
    axisFont = tmpFont;
    ShowSel::init();
}

void ShowSpectra::clear()
{
    if (V != NULL) delete V;
    ShowSel::clear();
}

/* Init with vectors ------------------------------------------------------- */
void ShowSpectra::initWithVectors(int _numRows, int _numCols,
                                  xmippCTVectors *_V, const char *_title)
{
    init();
    V = _V;
    fn = "";
    initFromVectors();
    setCaption(_title);
    NumRows = _numRows;
    NumCols = _numCols;
    initTable();
    initRightclickMenubar();
    repaint();
}

/* Read a Spectra ---------------------------------------------------------- */
void ShowSpectra::readFile(const FileName &_fn, double _minGray, double _maxGray)
{
    clear();
    fn = _fn;
    setCaption(fn.c_str());
    readDatFile(_fn);
}

void ShowSpectra::readDatFile(const FileName &_fn)
{
    std::ifstream fh_in(_fn.c_str());
    if (!fh_in)
        REPORT_ERROR(1, (std::string)"ShowSpectra::readFile: Cannot open" + _fn);
    V = new xmippCTVectors(fh_in);
    fh_in.close();

    annotateTime(_fn);
    initFromVectors();
}

/* Read vectors ------------------------------------------------------------ */
void ShowSpectra::initFromVectors()
{
    listSize = V->size();
    if (listSize == 0)
        REPORT_ERROR(1, "ShowSpectra::readFile: Input file is empty");
    imgnames        = new FileName[listSize];
    selstatus       = new bool[listSize];
    initContents();

    // Determine min, max and average
    minPixel = MAXFLOAT;
    maxPixel = MINFLOAT;
    for (long i = 0; i < V->size(); i++)
        for (int j = 0; j < V->theItems[0].size(); j++)
        {
            double val = V->theItems[i][j];
            if (minPixel > val) minPixel = val;
            if (maxPixel < val) maxPixel = val;
        }
    projXdim = projYdim = 100;
    for (long i = 0; i < V->size(); i++)
    {
        imgnames[i] = V->theTargets[i];
        selstatus[i] = SelLine::ACTIVE;
    }
}

/* Init Rightclick menubar ------------------------------------------------- */
void ShowSpectra::initRightclickMenubar()
{
#ifdef QT3_SUPPORT
    menubar = new Q3PopupMenu();
#else
    menubar = new QPopupMenu();
#endif
    setFileRightclickMenubar();

    // Options .............................................................
#ifdef QT3_SUPPORT
    options =  new Q3PopupMenu();
#else
    options = new QPopupMenu();
#endif
    setCommonOptionsRightclickMenubar();
    labeltype = Filename_LABEL;
    // spectra common options (pased to spectraSOM too)
    setCommonSpectraOptionsRightclickMenubar();

    // Statistics
    options->insertItem("View average and SD Images",  this,  SLOT(showStats()));
    options->insertItem("Show average and SD Spectra", this,  SLOT(showSpectraStats()));
    options->insertSeparator();

    // Insert options the menu
    menubar->insertItem("&Options", options);
    menubar->insertSeparator();

    // Inser Help and Quit
    insertGeneralItemsInRightclickMenubar();
}

/* Produce pixmap ---------------------------------------------------------- */
void ShowSpectra::paintCell(QPainter *p, int row, int col, const QRect & cr,
                            bool selected, const QColorGroup & cg)
{
    int scprojXdim = (int)((currScale * projXdim) / 100.0);
    int scprojYdim = (int)((currScale * projYdim) / 100.0);
    int i = indexOf(row, col);
    if (i >= listSize) return;
    int N = V->theItems[i].size();

    if (indexOf(row, col) >= listSize) return;
    QPixmap background(columnWidth(col), rowHeight(row));
    background.fill(backColor);
    p->drawPixmap(0, 0, background);
    p->setFont(axisFont);

    // Get minimum and Maximum of this spectrum
    double myMinValue, myMaxValue;
    if (!options->isItemEnabled(mi_Individual_norm))
    {
        // Individual normalization
        myMaxValue = MINFLOAT;
        myMinValue = MAXFLOAT;
        for (int l = 0; l < N ; l++)
        {
            double val = V->theItems[i][l];
            if (myMinValue > val) myMinValue = val;
            if (myMaxValue < val) myMaxValue = val;
        }
    }
    else
    {
        // Global normalization
        myMaxValue = maxPixel;
        myMinValue = minPixel;
    }

    // Calculate slope and intersection of linear transformation
    double slope = 0;
    if (myMinValue != myMaxValue)
        slope = (scprojXdim - 2 * offX) / (myMaxValue - myMinValue);

    // Paint X axis .........................................................
    QPen pen(axisColor);
    p->setPen(pen);
    p->drawLine(offX, scprojYdim - offY, scprojXdim - offX, scprojYdim - offY);
    double scaleX = (scprojXdim - 2 * offX) / N;
    int    istep = spacing;
    for (int l = x_tick_off - 1; l <= N; l += istep)
    {
        int x = offX + (int)(l * scaleX);
        // Draw legend
        if (!options->isItemEnabled(mi_showXlegend))
        {
            QString tmpS;
            tmpS.setNum(l + 1);
            p->drawText(x, (int)(scprojYdim - offY + p->font().pixelSize()), tmpS);
        }
        // Draw X grid lines
        if (!options->isItemEnabled(mi_showgrid))
            p->drawLine(x, offY, x, scprojYdim - offY);
        else
            p->drawLine(x, scprojYdim - offY, x, scprojYdim - offY - 3);
    }

    // Paint Y axis .........................................................
    p->drawLine(offX, offY, offX, scprojYdim - offY);
    for (int l = 0; l <= 3; l++)
    {
        // Draw legend
        if (!options->isItemEnabled(mi_showYlegend))
        {
            QString tmpS;
            if (l == 0) tmpS.setNum(myMinValue, 'f', 2);
            else      tmpS.setNum((myMinValue + l*((myMaxValue - myMinValue) / 3.0)), 'f', 2);
            p->drawText(2, (int)(scprojYdim - offY - (l*((scprojYdim - 2*offY) / 3.0))), tmpS);
        }

        // Draw Y grid lines
        if (!options->isItemEnabled(mi_showgrid))
            p->drawLine(offX, (int)(scprojYdim - offY - (l*((scprojYdim - 2*offY) / 3.0))),
                        scprojXdim - offX, (int)(scprojYdim - offY - (l*((scprojYdim - 2*offY) / 3.0))));
        else
            p->drawLine(offX, (int)(scprojYdim - offY - (l*((scprojYdim - 2*offY) / 3.0))),
                        2 + offX, (int)(scprojYdim - offY - (l*((scprojYdim - 2*offY) / 3.0))));
    }

    // Paint curves .........................................................
    pen.setColor(curveColor);
    pen.setWidth(3);
    p->setPen(pen);

    int      x0 = offX;
    double myY0 = offY + slope * (V->theItems[i][0] - myMinValue);
    int      y0 = scprojYdim - (int) myY0;

    for (int l = 1; l < N; l++)
    {
        int xF = offX + (int)(l * scaleX);
        double myYF = offY + slope * (V->theItems[i][l] - myMinValue);
        int yF = scprojYdim - (int) myYF;
        p->drawLine(x0, y0, xF, yF);
        x0=xF;
        y0=yF;
    }

    // Draw Frame and label
    drawFrameAndLabel(p, row, col, i, 1);
}

/* Open new file ----------------------------------------------------------- */
void ShowSpectra::openNewFile(const FileName &_fn)
{
    init();
    readFile(_fn);
    initTable();
    repaint();
}

/* Select by value ---------------------------------------------------------- */
void ShowSpectra::selectByValues()
{
    SpectraFilter* filter_window = new SpectraFilter(
        FLOOR(minPixel), CEIL(maxPixel), V->theItems[0], this,
        0, "new window", Qt::WDestructiveClose);
    filter_window->show();
}

/* Apply filter ------------------------------------------------------------ */
void ShowSpectra::applyFilter(const std::vector<int> &min, const std::vector<int> &max)
{
    int N = V->theItems[0].size();

    for (int i = 0; i < listSize; i++)
    {
        bool current_status = cellMarks[i];
        bool new_status = true;
        for (int j = 0; j < N; j++)
        {
            float v = V->theItems[i][j];
            if (v < min[j] || v > max[j])
            {
                new_status = false;
                break;
            }
        }
        if (current_status != new_status)
        {
            cellMarks[i] = new_status;
            updateCellIdx(i);
        }
    }
}

/* Change options ---------------------------------------------------------- */
void ShowSpectra::changeGrid()
{
    changeBoolOption(mi_showgrid, mi_hidegrid);
}
void ShowSpectra::changeXlegend()
{
    changeBoolOption(mi_showXlegend, mi_hideXlegend);
}
void ShowSpectra::changeYlegend()
{
    changeBoolOption(mi_showYlegend, mi_hideYlegend);
}

// Show Spectra Stats ------------------------------------------------------
void ShowSpectra::showSpectraStats()
{
    long counter = 0;
    // get number of marked spectra
    for (long i = 0; i < listSize; i++) if (cellMarks[i]) counter++;
    if (counter < 2)
        QMessageBox::about(this, "Error!", "No enough spectra selected\n");
    else
    {
        xmippCTVectors mySpectra(0, true);
        mySpectra.theItems.resize(counter);
        mySpectra.theTargets.resize(counter);
        long myIndex = 0;
        for (long i = 0; i < listSize; i++)
            if (cellMarks[i])
            {
                mySpectra.theItems[myIndex] = V->theItems[i];
                mySpectra.theTargets[myIndex] = V->theTargets[i];
                myIndex++;
            }

        xmippCTVectors *myVect = new xmippCTVectors(0, true);
        *myVect = mySpectra.getStatVector();
        ShowSpectra *myST = new ShowSpectra;
        myST->initWithVectors(1, 2, myVect, "Average and SD");
        myST->show();
    }
}

// Change colors -----------------------------------------------------------
void ShowSpectra::GUIchangeColor(QColor &_color, const char *_color_title)
{
    QColor tmpColor = QColorDialog::getColor(axisColor, this, _color_title);
    if (tmpColor.isValid())
    {
        _color = tmpColor;
        repaintContents();
    }
}

void ShowSpectra::changeBackColor()
{
    GUIchangeColor(backColor, "Back color");
}
void ShowSpectra::changeAxisColor()
{
    GUIchangeColor(axisColor, "Axis color");
}
void ShowSpectra::changeCurveColor()
{
    GUIchangeColor(curveColor, "Curve color");
}

// Change Font -------------------------------------------------------------
void ShowSpectra::changeFont()
{
    QFont tmpFont;
    bool ok;
    tmpFont = QFontDialog::getFont(&ok, axisFont, this, "Font type");
    if (ok)
    {
        axisFont = tmpFont;
        repaintContents();
    }
}
// Change Grid spacing in X axis -------------------------------------
void ShowSpectra::changeXstep()
{
    int N = V->theItems[0].size();
    ScrollParam* param_window;
    std::vector<float> min;
    min.push_back(1);
    min.push_back(1);
    std::vector<float> max;
    max.push_back(N);
    max.push_back(N);
    std::vector<float> initial_value;
    initial_value.push_back(spacing);
    initial_value.push_back(x_tick_off);
    std::vector<char *> prm_name;
    prm_name.push_back("spacing");
    prm_name.push_back("Tick offset");
    param_window = new ScrollParam(min, max, initial_value, prm_name,
                                   "Set spacing", 0, "new window", Qt::WDestructiveClose, 0);
    connect(param_window, SIGNAL(new_value(std::vector<float>)), this,
            SLOT(set_spacing(std::vector<float>)));
    param_window->setFixedSize(200, 175);
    param_window->show();
}

/****************************************************/
void ShowSpectra::set_spacing(std::vector<float> prm)
{
    spacing = (int) prm[0];
    x_tick_off = (int) prm[1];
    clearContents();
    repaintContents();
}

void ShowSpectra::setCommonSpectraOptionsRightclickMenubar()
{
    // Select by value
    options->insertItem("Select by spectral values", this, SLOT(selectByValues()));
    options->insertSeparator();

    // Show grid
    mi_showgrid = options->insertItem("Show Grid", this,  SLOT(changeGrid()));
    mi_hidegrid = options->insertItem("Hide Grid", this,  SLOT(changeGrid()));
    options->setItemEnabled(mi_showgrid, true);
    options->setItemEnabled(mi_hidegrid, false);
    options->insertSeparator();

    // Show X legend
    mi_showXlegend = options->insertItem("Show X legend", this,  SLOT(changeXlegend()));
    mi_hideXlegend = options->insertItem("Hide X legend", this,  SLOT(changeXlegend()));
    options->insertItem("&Change X-axis step...", this, SLOT(changeXstep()));
    options->setItemEnabled(mi_showXlegend, false);
    options->setItemEnabled(mi_hideYlegend, true);
    options->insertSeparator();

    // Show Y legend
    mi_showYlegend = options->insertItem("Show Y legend", this,  SLOT(changeYlegend()));
    mi_hideYlegend = options->insertItem("Hide Y legend", this,  SLOT(changeYlegend()));
    options->setItemEnabled(mi_showYlegend, false);
    options->setItemEnabled(mi_hideYlegend, true);
    options->insertSeparator();

    // Colors and change font
#ifdef QT3_SUPPORT
    Q3PopupMenu* colorMenu = new Q3PopupMenu();
#else
    QPopupMenu* colorMenu = new QPopupMenu();
#endif
    colorMenu->insertItem("&Background", this, SLOT(changeBackColor()));
    colorMenu->insertItem("&Curve", this, SLOT(changeCurveColor()));
    colorMenu->insertItem("&Axis", this, SLOT(changeAxisColor()));
    menubar->insertItem("&Colors", colorMenu);
    menubar->insertSeparator();
}

/* ------------------------------------------------------------------------- */
// Spectra filter constructor
SpectraFilter::SpectraFilter(int min, int max,
                             const std::vector<float> &_x, ShowSpectra *_show_spectra,
                             QWidget *parent, const char *name, int wflags):
#ifdef QT3_SUPPORT
    QWidget(parent, name, (Qt::WindowFlags) wflags)
#else
    QWidget(parent, name, wflags)
#endif
{
    __N = _x.size();
    __show_spectra = _show_spectra;

    // Allocate memory
    __current_values = new float[__N];
    if (!__current_values)
        REPORT_ERROR(1, "SpectraFilter: Cannot allocate memory");
    __scroll_min = new QScrollBar *[__N];
    if (!__scroll_min)
        REPORT_ERROR(1, "SpectraFilter: Cannot allocate memory");
    __scroll_max = new QScrollBar *[__N];
    if (!__scroll_max)
        REPORT_ERROR(1, "SpectraFilter: Cannot allocate memory");
    __label_min = new QLabel *[__N];
    if (!__label_min)
        REPORT_ERROR(1, "SpectraFilter: Cannot allocate memory");
    __label_max = new QLabel *[__N];
    if (!__label_max)
        REPORT_ERROR(1, "SpectraFilter: Cannot allocate memory");

    // Set this window caption
    setCaption("Spectra Filter");

    // Create a layout to position the widgets
#ifdef QT3_SUPPORT
    Q3BoxLayout *Layout = new Q3VBoxLayout(this, 10);
#else
    QBoxLayout* Layout = new QVBoxLayout(this, 10);
#endif

    // Create a grid layout to hold most of the widgets
#ifdef QT3_SUPPORT
    Q3GridLayout *grid = new Q3GridLayout(__N + 2, 5);
#else
    QGridLayout* grid = new QGridLayout(__N + 2, 5);
#endif
    Layout->addLayout(grid, 5);

    // Label
    QLabel     *label_min = new QLabel(this, "label");
    label_min->setFont(QFont("times", 12, QFont::Bold));
    label_min->setText("Minimum value");
    label_min->setFixedSize(label_min->sizeHint());
    grid->addWidget(label_min, 0, 1, Qt::AlignCenter);

    // Label
    QLabel     *label_max = new QLabel(this, "label");
    label_max->setFont(QFont("times", 12, QFont::Bold));
    label_max->setText("Maximum value");
    label_max->setFixedSize(label_max->sizeHint());
    grid->addWidget(label_max, 0, 3, Qt::AlignCenter);

    // Create all sliders
    for (int i = 0; i < __N; i++)
    {

        // Label
        __label_min[i] = new QLabel(this, "label");
        __label_min[i]->setFont(QFont("courier", 14));
        __label_min[i]->setText(integerToString(min, 3).c_str());
        __label_min[i]->setFixedSize(__label_min[i]->sizeHint());
        grid->addWidget(__label_min[i], i + 1, 0, Qt::AlignCenter);

        //Scroll Bar
        __scroll_min[i] = new QScrollBar(min, max, 1, 1, min,
                                         Qt::Horizontal, this, "scroll");
        __scroll_min[i]->setFixedWidth(100);
        __scroll_min[i]->setFixedHeight(15);
        grid->addWidget(__scroll_min[i], i + 1, 1, Qt::AlignCenter);
        connect(__scroll_min[i], SIGNAL(valueChanged(int)),
                SLOT(scrollValueChanged(int)));

        // Label
        QLabel     *label = new QLabel(this, "label");
        label->setFont(QFont("times", 12));
        label->setText(((std::string)"Harmonic " + integerToString(i + 1, 2)).c_str());
        label->setFixedSize(label->sizeHint());
        grid->addWidget(label, i + 1, 2, Qt::AlignCenter);

        //Scroll Bar
        __scroll_max[i] = new QScrollBar(min, max, 1, 1, max,
                                         Qt::Horizontal, this, "scroll");
        __scroll_max[i]->setFixedWidth(100);
        __scroll_max[i]->setFixedHeight(15);
        grid->addWidget(__scroll_max[i], i + 1, 3, Qt::AlignCenter);
        connect(__scroll_max[i], SIGNAL(valueChanged(int)),
                SLOT(scrollValueChanged(int)));

        // Label
        __label_max[i] = new QLabel(this, "label");
        __label_max[i]->setFont(QFont("courier", 14));
        __label_max[i]->setText(integerToString(max, 3).c_str());
        __label_max[i]->setFixedSize(__label_max[i]->sizeHint());
        grid->addWidget(__label_max[i], i + 1, 4, Qt::AlignCenter);

    }

    // Filter Button
    QPushButton *do_it;
    do_it = new QPushButton(this, "do_it");   // create button 3
    do_it->setFont(QFont("times", 12, QFont::Bold));
    do_it->setText("Filter");
    do_it->setFixedSize(do_it->sizeHint());
    grid->addWidget(do_it, __N + 1, 0, Qt::AlignCenter);
    QToolTip::add(do_it, "Select Spectra according to filter");
    connect(do_it, SIGNAL(clicked()), SLOT(but_ok_clicked()));
}

// Spectra filter destructor
SpectraFilter::~SpectraFilter()
{
    delete __current_values;
}

// One of the sliders changed ----------------------------------------------
void SpectraFilter::scrollValueChanged(int new_val)
{
    int dummy = new_val;  // This is to cheat some compilers
    if (dummy);         // who complain that the values
    // are never used

    // Read all sliders
    for (int i = 0; i < __N; i++)
    {
        __label_min[i]->setText((integerToString(__scroll_min[i]->value(), 3)).c_str());
        __label_max[i]->setText((integerToString(__scroll_max[i]->value(), 3)).c_str());
    }
}

/* Filter ------------------------------------------------------------------ */
void SpectraFilter::but_ok_clicked()
{
    std::vector<int> min(__N);
    std::vector<int> max(__N);
    for (int i = 0; i < __N; i++)
    {
        min[i] = __scroll_min[i]->value();
        max[i] = __scroll_max[i]->value();
    }
    __show_spectra->applyFilter(min, max);
}
