/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Carlos Manzanares       (cmanzana@cnb.csic.es)
 *              Arun Kulshreshth        (arun_2000_iitd@yahoo.com)
 *              Enrique Recarte Llorens (erecallo@hotmail.com)
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

#include "widget_micrograph.h"
#include "filter_menu.h"
#include "auto_menu.h"
#include "image_micrograph.h"
#include "image_overview_micrograph.h"

#include <data/micrograph.h>
#include <data/args.h>
#include <graphics/show_tools.h>

#include <sys/time.h>

#include <qinputdialog.h>
#include <qmessagebox.h>
#include <qlabel.h>
#include <qlayout.h>
#include <qpushbutton.h>

#ifdef QT3_SUPPORT
#include <q3grid.h>
//Added by qt3to4:
#include <Q3GridLayout>
#include <Q3VBoxLayout>
#else
#include <qgrid.h>
#endif

/* Constructor ------------------------------------------------------------ */
QtWidgetMicrograph::QtWidgetMicrograph(QtMainWidgetMark *_mainWidget,
                                       QtFiltersController *_f,
                                       Micrograph *_m) :
        QWidget((QWidget*) _mainWidget)
{
    __filtersController              = _f;
    __m                              = NULL;
    __activeFamily                   = -1;
    __tilted                         = false;
#ifdef QT3_SUPPORT

    __gridLayout                     = new Q3VBoxLayout(this);
#else

    __gridLayout                     = new QVBoxLayout(this);
#endif

    __menuBar                        = new QMenuBar(this);
    __menuBar->setSeparator(QMenuBar::InWindowsStyle);
    __mImage                         = new QtImageMicrograph(0);
    __mImage->setWidgetMicrograph(this);
    __mImageOverview                 = new QtImageOverviewMicrograph(this);
    __mImageOverview->setWidgetMicrograph(this);
    __file_menu                      = NULL;
    __ellipse_radius                 = 5;
    __ellipse_type                   = MARK_CIRCLE;

    __mImage->setFiltersController(_f);
    __mImageOverview->setFiltersController(_f);

    connect(__mImageOverview, SIGNAL(signalSetCoords(int, int)),
            __mImage, SLOT(slotSetCoords(int, int)));
    connect(__mImage, SIGNAL(signalSetCoords(int, int)),
            __mImageOverview, SLOT(slotSetCoords(int, int)));
    connect(__mImage, SIGNAL(signalSetWidthHeight(int, int)),
            __mImageOverview, SLOT(slotSetWidthHeight(int, int)));
    // COSS: Make this connection to update the overview when
    //       a particle is deleted or moved
    //    connect(__mImage, SIGNAL(signalRepaint(void)),
    //            __mImageOverview, SLOT(repaint(void)));
    connect(__mImageOverview, SIGNAL(signalRepaint(void)),
            __mImage, SLOT(repaint(void)));
    connect(__mImage, SIGNAL(signalRepaint(void)),
            __mImage, SLOT(repaint(void)));
    connect(__mImageOverview, SIGNAL(signalRepaint(void)),
            __mImageOverview, SLOT(repaint(void)));
    connect(__mImage, SIGNAL(signalAddCoordOther(int, int, int)),
            this, SLOT(slotDrawEllipse(int, int, int)));

    connect(this, SIGNAL(signalActiveFamily(int)),
            __mImage, SLOT(slotActiveFamily(int)));
    connect(this, SIGNAL(signalActiveFamily(int)),
            __mImageOverview, SLOT(slotActiveFamily(int)));
    connect(this, SIGNAL(signalRepaint()),
            this, SLOT(slotRepaint()));

#ifdef QT3_SUPPORT

    Q3Accel *ctrl = new Q3Accel(this);
#else

    QAccel *ctrl  = new QAccel(this);
#endif

    ctrl->connectItem(ctrl->insertItem(Qt::Key_G + Qt::CTRL, 200),
                      this, SLOT(slotChangeContrast(void)));
#ifdef QT3_SUPPORT

    Q3Accel *ctrl2 = new Q3Accel(this);
#else

    QAccel *ctrl2  = new QAccel(this);
#endif

    ctrl2->connectItem(ctrl2->insertItem(Qt::Key_R + Qt::CTRL, 200),
                       this, SLOT(slotChangeCircleRadius(void)));

    setMicrograph(_m);

    __gridLayout->setMenuBar(__menuBar);
    __gridLayout->addWidget(__mImageOverview);

    openMenus();
}

/* Destructor--------------------------------------------------------------- */
QtWidgetMicrograph::~QtWidgetMicrograph()
{
    delete __mImage;
    delete __mImageOverview;
    delete __menuBar;
    delete __gridLayout;
}

/* Open all windows -------------------------------------------------------- */
void QtWidgetMicrograph::openAllWindows()
{
    __mImage->show();
    __mImageOverview->show();
}

/* Set Micrograph ---------------------------------------------------------- */
void QtWidgetMicrograph::setMicrograph(Micrograph *_m)
{
    if (_m != NULL)
    {
        __m = _m;
        __mImage->setMicrograph(_m);
        __mImageOverview->setMicrograph(_m);
    }
}

/* Set Output directory ---------------------------------------------------- */
void QtWidgetMicrograph::setOutputRoot(const FileName& outputRoot)
{
    __outputRoot=outputRoot;
    FileName fn1(__outputRoot+".Common.pos"), fn2(__outputRoot+".Common.auto.pos");
    if (fn1.exists())
    	__file_menu->loadCoords(__outputRoot+".Common.pos");
    if (fn2.exists())
    	__file_menu->loadCoords(__outputRoot+".Common.auto.pos");
}

/* Set the class for automatic particle picking ---------------------------- */
void QtWidgetMicrograph::setAutoParticlePicking(
    AutoParticlePickingQt *_autoPicking)
{
    __autoPicking=_autoPicking;
    if (_autoPicking!=NULL)
    {
    	QtAutoMenu *autoMenu = new QtAutoMenu(this);
    	addMenuItem("&AutoSelection", (QtPopupMenuMark *)(autoMenu));
    }
}

AutoParticlePickingQt * QtWidgetMicrograph::getAutoParticlePicking() const
{
    return __autoPicking;
}

/* Open menus -------------------------------------------------------------- */
void QtWidgetMicrograph::openMenus()
{
    __file_menu = new QtFileMenu(this);
    connect(__mImage, SIGNAL(signalAddCoordOther(int, int, int)),
            __file_menu, SLOT(slotCoordChange()));

    QtFilterMenu *filterMenu = new QtFilterMenu(this);


    addMenuItem("&File", (QtPopupMenuMark *)(__file_menu));
    addMenuItem("F&ilters", (QtPopupMenuMark *)(filterMenu));

    connect(__file_menu, SIGNAL(signalAddFamily(const char *)),
            this, SLOT(slotAddFamily(const char*)));
    connect((QObject*)filterMenu, SIGNAL(signalAdjustContrast()),
            this, SLOT(slotChangeContrast(void)));
    connect((QObject*)filterMenu, SIGNAL(signalCrop()),
            this, SLOT(slotChangeCrop(void)));
    connect((QObject*)filterMenu, SIGNAL(signalAddFilter()),
            (QObject*)__filtersController, SLOT(slotAddFilter()));
    connect((QObject*)filterMenu, SIGNAL(signalCleanFilters()),
            (QObject*)__filtersController, SLOT(slotCleanFilters()));
    connect((QObject*)filterMenu, SIGNAL(signalCleanFilters()),
            this, SLOT(repaint()));

    // *** Add your own menus
}

/* Add family -------------------------------------------------------------- */
int QtWidgetMicrograph::add_family(std::vector<ParticleQt> &_list,
                                   const std::string &_family_name)
{
    int ilabel = __m->add_label(_family_name);
    int imax = _list.size();
    for (int i = 0; i < imax; i++)
    {
        int idx = __m->add_coord(_list.at(i).x, _list.at(i).y, ilabel,
                                 _list.at(i).cost);
        _list.at(i).idx = idx;
    }
    return ilabel;
}

/* Load models ------------------------------------------------------------- */
void QtWidgetMicrograph::changeContrast(int _mingray, int _maxgray, float _gamma)
{
    __mingray = _mingray;
    __maxgray = _maxgray;
    __gamma  = _gamma;
    __mImage->changeContrast(_mingray, _maxgray, _gamma);
    __mImageOverview->changeContrast(_mingray, _maxgray, _gamma);
}

void QtWidgetMicrograph::changeMarkType(int _type)
{
    __ellipse_type = _type;
    __mImage->__ellipse_type = __ellipse_type;
    __mImage->repaint(true);
}

void QtWidgetMicrograph::changeCircleRadius(float _circle_radius)
{
    __ellipse_radius = _circle_radius;
    __mImage->__ellipse_radius = __ellipse_radius;
    __mImage->repaint(true);
}

void QtWidgetMicrograph::repaint()
{
    __mImage->repaint(true);
    __mImageOverview->repaint(true);
}

void QtWidgetMicrograph::slotDrawEllipse(int _x, int _y, int _f)
{
    __mImage->drawEllipse(_x, _y, _f, __ellipse_radius, __ellipse_type);
    __mImageOverview->drawEllipse(_x, _y, _f);
}

void QtWidgetMicrograph::slotDrawLastEllipse(int _x, int _y, int _f)
{
    __mImage->drawEllipse(_x, _y, _f, __ellipse_radius, __ellipse_type);
    __mImageOverview->drawEllipse(_x, _y, _f);
    __mImage->drawLastEllipse(_x, _y, _f, __ellipse_radius, __ellipse_type);
}

/* Active family ----------------------------------------------------------- */
void QtWidgetMicrograph::slotActiveFamily(int _f)
{
    __activeFamily = _f;
    emit signalActiveFamily(_f);
}

void QtWidgetMicrograph::slotAddFamily(const char *_familyName)
{
    emit signalAddFamily(_familyName);
}

void QtWidgetMicrograph::slotDeleteMarkOther(int _coord)
{
    __m->coord(_coord).valid = false;
    emit signalRepaint();
}

void QtWidgetMicrograph::slotQuit()
{
    __file_menu->slotQuit();
}

void QtWidgetMicrograph::slotChangeContrast()
{
    AdjustContrastWidget *adjustContrast = new
                                           AdjustContrastWidget(0, 255, 1.0F, this,
                                                                0, "new window", Qt::WDestructiveClose);
    adjustContrast->show();
}

void QtWidgetMicrograph::slotChangeCrop()
{
    CropWidget *crop = new CropWidget(this, 0, "new window", Qt::WDestructiveClose);
    connect(crop, SIGNAL(new_value(std::vector<int>)),
            __mImageOverview, SLOT(slotDrawCropArea(std::vector<int>)));
    crop->show();
}

void QtWidgetMicrograph::slotChangeCircleRadius()
{
    AdjustCircleRadiustWidget *adjustCircleRadius = new
            AdjustCircleRadiustWidget(0, 255, 10, this,
                                      0, "new window", Qt::WDestructiveClose);
    adjustCircleRadius->show();
}

// Learn particles
bool QtWidgetMicrograph::learnParticles()
{
    if (__ellipse_radius<12 && !__autoPicking->__is_model_loaded)
    {
        if (!QMessageBox::information(this, "Mark",
                                      (std::string)"The particle radius is "+
                                      floatToString(__ellipse_radius,3,2)+
                                      "\nAre you sure to continue?",
                                      "&Yes", "&No",
                                      0,      // Enter == button 0
                                      1))
        {
            __autoPicking->learnParticles(__ellipse_radius);
            return true;
        }
        else
        	return false;
    }
    else
    {
        __autoPicking->learnParticles(__ellipse_radius);
        return true;
    }
}

// Automatically Select Particles
void QtWidgetMicrograph::automaticallySelectParticles()
{
    // Add a family to the micrograph is necessary
    bool addFamily=true;
    for (int i = 0; i < __m->LabelNo(); i++)
        if (__m->get_label(i)=="auto")
        {
            addFamily=false;
            break;
        }
    if (addFamily)
    {
        __m->add_label("auto");
        emit signalAddFamily("auto");
    }

    int Nalive=__autoPicking->automaticallySelectParticles();
    emit signalRepaint();
    if (Nalive>0)
    {
        ScrollParam* param_window;
        param_window = new ScrollParam(0, 1, 0.01, "Quality",
                                       "Restrict selection", 0, "new window", Qt::WDestructiveClose);

        // Connect its output to my input (set_spacing)
        connect( param_window, SIGNAL(new_value(float)), this,
                 SLOT(slotRestrictSelection(float)));

        // Show
        param_window->setFixedSize(300,150);
        param_window->show();
    }
    slotActiveFamily(0);
}

// Save models
void QtWidgetMicrograph::saveModels(bool askFilename)
{
    // Get the rootname
    std::string fn_root="";
    if (askFilename)
    {
        bool ok;
        QString qfn_root = QInputDialog::getText("Saving model",
                           "Model", QLineEdit::Normal,
                           "Model", &ok);
        if (!ok || qfn_root.isEmpty())
            return;
        fn_root = qfn_root.ascii();
    }
    __autoPicking->saveModels(fn_root);
}

/* Restrict Selection ------------------------------------------------------ */
void QtWidgetMicrograph::slotRestrictSelection(float _cost)
{
    __mImage->__minCost=__mImageOverview->__minCost=_cost;
    __autoPicking->restrictSelection(_cost);
    emit signalRepaint();
}

/* ------------------------------------------------------------------------- */
/* Adjust contrast widget                                                    */
/* ------------------------------------------------------------------------- */
/* AdjustContrastWidget ---------------------------------------------------- */
// Constructor
AdjustContrastWidget::AdjustContrastWidget(int min, int max, float gamma,
        QtWidgetMicrograph *_qtwidgetmicrograph,
        QWidget *parent, const char *name, int wflags):
#ifdef QT3_SUPPORT
        QWidget(parent, name, (Qt::WindowFlags) wflags)
#else
        QWidget(parent, name, wflags)
#endif
{
    __qtwidgetmicrograph = _qtwidgetmicrograph;

    // Set this window caption
    setCaption("Adjust Contrast");

    // Create a layout to position the widgets
#ifdef QT3_SUPPORT

    Q3BoxLayout *Layout = new Q3VBoxLayout(this, 10);
#else

    QBoxLayout *Layout = new QVBoxLayout(this, 10);
#endif

    // Create a grid layout to hold most of the widgets
#ifdef QT3_SUPPORT

    Q3GridLayout *grid = new Q3GridLayout(3, 3);
#else

    QGridLayout *grid = new QGridLayout(3, 3);
#endif

    Layout->addLayout(grid, 5);

    // Minimum
    QLabel     *label_min = new QLabel(this, "label");
    label_min->setFont(QFont("times", 12, QFont::Bold));
    label_min->setText("Minimum");
    label_min->setFixedSize(label_min->sizeHint());
    grid->addWidget(label_min, 0, 0, Qt::AlignCenter);

    __scroll_min = new QScrollBar(0, 255, 1, 1, min,
                                  Qt::Horizontal, this, "scroll");
    __scroll_min->setFixedWidth(100);
    __scroll_min->setFixedHeight(15);
    grid->addWidget(__scroll_min, 0, 1, Qt::AlignCenter);
    connect(__scroll_min, SIGNAL(valueChanged(int)),
            SLOT(scrollValueChanged(int)));

    __label_min = new QLabel(this, "label");
    __label_min->setFont(QFont("courier", 14));
    __label_min->setText(integerToString(min, 3).c_str());
    __label_min->setFixedSize(__label_min->sizeHint());
    grid->addWidget(__label_min, 0, 2, Qt::AlignCenter);

    // Maximum
    QLabel     *label_max = new QLabel(this, "label");
    label_max->setFont(QFont("times", 12, QFont::Bold));
    label_max->setText("Maximum");
    label_max->setFixedSize(label_max->sizeHint());
    grid->addWidget(label_max, 1, 0, Qt::AlignCenter);

    __scroll_max = new QScrollBar(0, 255, 1, 1, max,
                                  Qt::Horizontal, this, "scroll");
    __scroll_max->setFixedWidth(100);
    __scroll_max->setFixedHeight(15);
    grid->addWidget(__scroll_max, 1, 1, Qt::AlignCenter);
    connect(__scroll_max, SIGNAL(valueChanged(int)),
            SLOT(scrollValueChanged(int)));

    __label_max = new QLabel(this, "label");
    __label_max->setFont(QFont("courier", 14));
    __label_max->setText(integerToString(max, 3).c_str());
    __label_max->setFixedSize(__label_max->sizeHint());
    grid->addWidget(__label_max, 1, 2, Qt::AlignCenter);

    // Gamma
    QLabel     *label_gamma = new QLabel(this, "label");
    label_gamma->setFont(QFont("times", 12, QFont::Bold));
    label_gamma->setText("Gamma");
    label_gamma->setFixedSize(label_gamma->sizeHint());
    grid->addWidget(label_gamma, 2, 0, Qt::AlignCenter);

    __scroll_gamma = new QScrollBar(0, 40, 1, 1, (int)(10*gamma),
                                    Qt::Horizontal, this, "scroll");
    __scroll_gamma->setFixedWidth(100);
    __scroll_gamma->setFixedHeight(15);
    grid->addWidget(__scroll_gamma, 2, 1, Qt::AlignCenter);
    connect(__scroll_gamma, SIGNAL(valueChanged(int)),
            SLOT(scrollValueChanged(int)));

    __label_gamma = new QLabel(this, "label");
    __label_gamma->setFont(QFont("courier", 14));
    __label_gamma->setText(floatToString(gamma, 3, 2).c_str());
    __label_gamma->setFixedSize(__label_gamma->sizeHint());
    grid->addWidget(__label_gamma, 2, 2, Qt::AlignCenter);

}

// One of the sliders changed ----------------------------------------------
void AdjustContrastWidget::scrollValueChanged(int new_val)
{
    __label_min  ->setText(integerToString(__scroll_min  ->value(), 3).c_str());
    __label_max  ->setText(integerToString(__scroll_max  ->value(), 3).c_str());
    __label_gamma->setText(floatToString((__scroll_gamma->value()) / 10.0, 3, 2).c_str());
    __qtwidgetmicrograph->changeContrast(__scroll_min->value(),
                                         __scroll_max->value(), __scroll_gamma->value() / 10.0);
}

/* ------------------------------------------------------------------------- */
/* Crop widget                                                               */
/* ------------------------------------------------------------------------- */
/* CropWidget -------------------------------------------------------------- */
// Constructor
CropWidget::CropWidget(QtWidgetMicrograph *_qtwidgetmicrograph,
                       QWidget *parent, const char *name, int wflags):
#ifdef QT3_SUPPORT
        QWidget(parent, name, (Qt::WindowFlags) wflags)
#else
        QWidget(parent, name, wflags)
#endif
{
    __qtwidgetmicrograph = _qtwidgetmicrograph;
    int Xdim, Ydim;
    __qtwidgetmicrograph->getMicrograph()->size(Xdim, Ydim);

    // Set this window caption
    setCaption("Crop micrograph");

    // Create a layout to position the widgets
#ifdef QT3_SUPPORT

    Q3BoxLayout *Layout = new Q3VBoxLayout(this, 10);
#else

    QBoxLayout *Layout = new QVBoxLayout(this, 10);
#endif

    // Create a grid layout to hold most of the widgets
#ifdef QT3_SUPPORT

    Q3GridLayout *grid = new Q3GridLayout(6, 3);
#else

    QGridLayout *grid = new QGridLayout(6, 3);
#endif

    Layout->addLayout(grid, 5);

    // Layout the four bars
    std::vector<int> min, max, init_value;
    std::vector<char *> prm_name;
    prm_name.push_back("x0");
    min.push_back(0);
    max.push_back(Xdim);
    init_value.push_back(ROUND(0.25*Xdim));
    prm_name.push_back("y0");
    min.push_back(0);
    max.push_back(Ydim);
    init_value.push_back(ROUND(0.25*Ydim));
    prm_name.push_back("xF");
    min.push_back(0);
    max.push_back(Xdim);
    init_value.push_back(ROUND(0.75*Xdim));
    prm_name.push_back("yF");
    min.push_back(0);
    max.push_back(Ydim);
    init_value.push_back(ROUND(0.75*Ydim));
    for (int i = 0; i < min.size(); i++)
    {

        // Add Parameter name
        QLabel     *lab1 = new QLabel(this, "lab1");
        lab1->setFont(QFont("times", 12, QFont::Bold));
        lab1->setText(prm_name[i]);
        lab1->setFixedSize(lab1->sizeHint());
        grid->addWidget(lab1, i, 0, Qt::AlignLeft);

        // Add Scroll Bar
        QScrollBar  *scroll_aux = new QScrollBar(min[i], max[i], 1, 50, (int)init_value[i], Qt::Horizontal, this, "scroll");
        scroll_aux->setFixedWidth(100);
        scroll_aux->setFixedHeight(15);
        grid->addWidget(scroll_aux, i, 1, Qt::AlignCenter);
        __scroll.push_back(scroll_aux);

        // Label for the current value
        QLabel * value_lab_aux;
        value_lab_aux = new QLabel(this, "value_lab");
        value_lab_aux->setFont(QFont("times", 12));
        value_lab_aux->setNum(init_value[i]);
        grid->addWidget(value_lab_aux, i, 2, Qt::AlignLeft);
        __label.push_back(value_lab_aux);

        connect(scroll_aux, SIGNAL(valueChanged(int)), SLOT(scrollValueChanged(int)));
    }

    // Layout the output name
    QLabel     *lab2 = new QLabel(this, "lab2");
    lab2->setFont(QFont("times", 12, QFont::Bold));
    lab2->setText("Output image");
    lab2->setFixedSize(lab2->sizeHint());
    grid->addWidget(lab2, min.size() + 1, 0, Qt::AlignLeft);
    __outputNameLineEdit = new QLineEdit(this, "output name");
    grid->addWidget(__outputNameLineEdit, min.size() + 1, 1, Qt::AlignLeft);

    // Cancel Button
    QPushButton *cancel;
    cancel = new QPushButton(this, "cancel");     // create button 1
    cancel->setFont(QFont("times", 12, QFont::Bold));
    cancel->setText("Cancel");
    cancel->setFixedSize(cancel->sizeHint());
    grid->addWidget(cancel, min.size() + 2, 0, Qt::AlignVCenter);
    connect(cancel, SIGNAL(clicked()), this, SLOT(cancel()));

    // OK button
    QPushButton *do_it;
    do_it = new QPushButton(this, "do_it");     // create button 3
    do_it->setFont(QFont("times", 12, QFont::Bold));
    do_it->setText("Ok");
    do_it->setFixedHeight(do_it->sizeHint().height());
    do_it->setFixedWidth(80);
    grid->addWidget(do_it, min.size() + 2, 2, Qt::AlignVCenter);
    connect(do_it, SIGNAL(clicked()), this, SLOT(accept()));

    __qtwidgetmicrograph->overview()->init_crop_area();
}

// Destructor --------------------------------------------------------------
CropWidget::~CropWidget()
{
    for (int i = 0; i < __label.size(); i++)
    {
        delete __label[i];
        delete __scroll[i];
    }
}

// One of the sliders changed ----------------------------------------------
void CropWidget::scrollValueChanged(int new_val)
{
    std::vector<int> value;
    // Get values
    for (int i = 0; i < __label.size(); i++)
    {
        int v = __scroll[i]->value();
        value.push_back(ROUND((float)v));
    }

    // Check value validity
    value[0] = XMIPP_MIN(value[0], value[2]);
    value[2] = XMIPP_MAX(value[0], value[2]);
    value[1] = XMIPP_MIN(value[1], value[3]);
    value[3] = XMIPP_MAX(value[1], value[3]);

    // Set these values
    for (int i = 0; i < __label.size(); i++)
    {
        __label[i]->setNum(value[i]);
        __scroll[i]->setValue(value[i]);
    }

    emit new_value(value);
}

void CropWidget::accept()
{
    __qtwidgetmicrograph->overview()->finish_crop_area();
    // Get values
    std::vector<int> value;
    for (int i = 0; i < __label.size(); i++)
        value.push_back(__scroll[i]->value());

    // Get output image
    std::string fn_out = __outputNameLineEdit->text().ascii();
    if (fn_out == "")
    {
        QMessageBox::information(this, "Mark",
                                 "The output image is empty\n Cropping is not carried out\n");
        close();
        return;
    }

    // Do the cropping
    int w = value[2] - value[0];
    int h = value[3] - value[1];
    std::string command = (std::string)"xmipp_window_micrograph " +
                          "-i " + __qtwidgetmicrograph->getMicrograph()->micrograph_name() +
                          " -o " + fn_out +
                          " -size " + integerToString(w, 0) + " " + integerToString(h, 0) +
                          " -top_left_corner " + integerToString(value[0], 0) + " " + integerToString(value[1], 0);
    std::cout << "Executing:\n" << command << std::endl;
    system(command.c_str());

    // Close the parameters window
    close();
}

void CropWidget::cancel()
{
    __qtwidgetmicrograph->overview()->finish_crop_area();
    close();
}

/* ------------------------------------------------------------------------- */
/* Adjust circle widget                                                      */
/* ------------------------------------------------------------------------- */
/* AdjustCircleWidget ------------------------------------------------------ */
// Constructor
AdjustCircleRadiustWidget::AdjustCircleRadiustWidget(int min, int max,
        int start_with, QtWidgetMicrograph *_qtwidgetmicrograph,
        QWidget *parent, const char *name, int wflags):
#ifdef QT3_SUPPORT
        QWidget(parent, name, (Qt::WindowFlags) wflags)
#else
        QWidget(parent, name, wflags)
#endif
{
    __qtwidgetmicrograph = _qtwidgetmicrograph;

    // Set this window caption
    setCaption("Change Circle Radius");

    // Create a layout to position the widgets
#ifdef QT3_SUPPORT

    Q3BoxLayout *Layout = new Q3VBoxLayout(this, 3);
#else

    QBoxLayout *Layout = new QVBoxLayout(this, 3);
#endif

    // Create a grid layout to hold most of the widgets
#ifdef QT3_SUPPORT

    Q3GridLayout *grid = new Q3GridLayout(1, 3);
#else

    QGridLayout *grid = new QGridLayout(1, 3);
#endif

    Layout->addLayout(grid);

    // Radius
    QLabel     *label_radius = new QLabel(this, "label");
    label_radius->setFont(QFont("times", 12, QFont::Bold));
    label_radius->setText("Radius");
    label_radius->setFixedSize(label_radius->sizeHint());
    grid->addWidget(label_radius, 0, 0, Qt::AlignCenter);

    __scroll_radius = new QScrollBar(0, 255, 1, 10, start_with,
                                     Qt::Horizontal, this, "scroll");
    __scroll_radius->setFixedWidth(100);
    __scroll_radius->setFixedHeight(15);
    grid->addWidget(__scroll_radius, 0, 1, Qt::AlignCenter);
    connect(__scroll_radius, SIGNAL(valueChanged(int)),
            SLOT(scrollValueChanged(int)));

    __label_radius = new QLabel(this, "label");
    __label_radius->setFont(QFont("courier", 14));
    __label_radius->setText(integerToString(start_with, 3).c_str());
    __label_radius->setFixedSize(__label_radius->sizeHint());
    grid->addWidget(__label_radius, 0, 2, Qt::AlignCenter);
}

void AdjustCircleRadiustWidget::scrollValueChanged(int new_val)
{
    __label_radius  ->setText(integerToString(__scroll_radius  ->value(), 3).c_str());
    __qtwidgetmicrograph->changeCircleRadius(__scroll_radius->value());
}
