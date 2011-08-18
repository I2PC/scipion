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

#ifndef __QT_WIDGET_MICROGRAPH_HH__
#define __QT_WIDGET_MICROGRAPH_HH__

#include <qwidget.h>
#include <qpainter.h>
#include <qlayout.h>
#include <qmenubar.h>
#include <qscrollbar.h>
#include <qlabel.h>
#include <qlineedit.h>

#ifdef QT3_SUPPORT
#include <q3accel.h>
//Added by qt3to4:
#include <Q3PopupMenu>
#include <Q3VBoxLayout>
#else
#include <qaccel.h>
#endif

#include "image_micrograph.h"
#include "image_overview_micrograph.h"
#include "file_menu.h"

#include <reconstruction/micrograph_automatic_picking_for_qt.h>

#include <vector>

/* Forward declarations ---------------------------------------------------- */
class QtMainWidgetMark;
class QtImageMicrograph;
class QtPopupMenuMark;
class QtFiltersController;
class Micrograph;

/* Widget for the micrograph ----------------------------------------------- */
class QtWidgetMicrograph : public QWidget
{
    Q_OBJECT

public:
    Micrograph                *__m;
    QtFiltersController       *__filtersController;
    int                        __activeFamily;
    QMenuBar                  *__menuBar;
    QtImageMicrograph         *__mImage;
    QtImageOverviewMicrograph *__mImageOverview;
#ifdef QT3_SUPPORT

    Q3VBoxLayout               *__gridLayout;
#else

    QVBoxLayout               *__gridLayout;
#endif

    FileName                   __outputRoot;
    QtFileMenu                *__file_menu;
    bool                       __tilted;
    int                        __mingray;
    int                        __maxgray;
    float                      __gamma;
    float                      __ellipse_radius;
    int                        __ellipse_type;
    AutoParticlePickingQt       *__autoPicking;
public:
    // Constructor
    QtWidgetMicrograph(QtMainWidgetMark *_mainWidget,
                       QtFiltersController *_f,
                       Micrograph *_m = NULL);
    ~QtWidgetMicrograph();

    // Open all windows
    void openAllWindows();

    // Set Micrograph
    void setMicrograph(Micrograph *_m);

    // Get Micrograph
    Micrograph *getMicrograph()
    {
        return(__m);
    }

    // Set Output Directory
    void setOutputRoot(const FileName& outputRoot);

    // Get Output Directory
    const FileName& getOutputRoot() const
    {
    	return __outputRoot;
    }

    // Set this as tilted micrograph
    void setTilted()
    {
        __tilted = TRUE;
        __mImage->setTilted();
    }

    // Is tilted?
    bool isTilted()
    {
        return __tilted;
    }

    // Set Automatic Particle Picking
    void setAutoParticlePicking(AutoParticlePickingQt *_autoPicking);

    // Set Automatic Particle Picking
    AutoParticlePickingQt * getAutoParticlePicking() const;

    // Get filters controller
    QtFiltersController *getFiltersController()
    {
        return(__filtersController);
    }

    // Get active family
    int activeFamily()
    {
        return(__activeFamily);
    }

    // Get overview
    QtImageOverviewMicrograph *overview()
    {
        return(__mImageOverview);
    }

    // Get Image
    QtImageMicrograph *image()
    {
        return(__mImage);
    }

    // Get Filemenu
    QtFileMenu *file_menu()
    {
        return __file_menu;
    }

    // Add menu item
    void addMenuItem(const char *_msg, const QtPopupMenuMark *_item)
    {
#ifdef QT3_SUPPORT
        __menuBar->insertItem(_msg, (Q3PopupMenu*)_item);
#else

        __menuBar->insertItem(_msg, (QPopupMenu*) _item);
#endif

    }

    // Draw axis
    void draw_axis(double _ang)
    {
        __mImageOverview->enableAxis();
        __mImageOverview->draw_axis(_ang);
    }

    // Open menu.
    // Add your menus to this function
    void openMenus();

    // Change contrast
    void changeContrast(int _mingray, int _maxgray, float _gamma);

    // Change mark type
    void changeMarkType(int _type);

    // Change circle radius
    void changeCircleRadius(float _circle_radius);

    // Repaint
    void repaint();

    // Add family.
    // The family label is returned
    int add_family(std::vector<ParticleQt> &_list,
                   const std::string &_family_name);

    // Learn particles
    bool learnParticles();

    // Automatically Select Particles
    void automaticallySelectParticles();

    // Save models
    void saveModels(bool askFilename);

public slots:
    void slotActiveFamily(int _f);
    void slotAddFamily(const char *_familyName);
    void slotDeleteMarkOther(int _coord);
    void slotRepaint()
    {
        repaint();
    }
    void slotDrawEllipse(int _x, int _y, int _f);
    void slotDrawLastEllipse(int _x, int _y, int _f);
    void slotQuit();
    void slotChangeContrast();
    void slotChangeCrop();
    void slotChangeCircleRadius();
    void slotRestrictSelection(float _cost);
signals:
    void signalActiveFamily(int _f);
    void signalAddFamily(const char *_familyName);
    void signalRepaint();
};

/** Class to adjust contrast
*/
class AdjustContrastWidget : public QWidget
{
    Q_OBJECT
public:
    /** Constructor */
    AdjustContrastWidget(int min, int max, float gamma,
                         QtWidgetMicrograph *_qtwidgetmicrograph,
                         QWidget *parent = 0, const char *name = 0, int wflags = 0);
private:
    QtWidgetMicrograph *__qtwidgetmicrograph;
    QScrollBar        *__scroll_min;
    QScrollBar        *__scroll_max;
    QScrollBar        *__scroll_gamma;
    QLabel            *__label_min;
    QLabel            *__label_max;
    QLabel            *__label_gamma;
private slots:
    void scrollValueChanged(int);
};

/** Class to adjust contrast
*/
class CropWidget : public QWidget
{
    Q_OBJECT
public:
    /** Constructor */
    CropWidget(QtWidgetMicrograph *_qtwidgetmicrograph,
               QWidget *parent = 0, const char *name = 0, int wflags = 0);

    /** Destructor */
    ~CropWidget();
private:
    QtWidgetMicrograph      *__qtwidgetmicrograph;
    std::vector < QScrollBar * >  __scroll;
    std::vector < QLabel * >      __label;
    QLineEdit               *__outputNameLineEdit;
private slots:
    void scrollValueChanged(int);
    void accept();
    void cancel();
signals:
    void new_value(std::vector<int>);
};

/** Class to adjust circle radius
*/
class AdjustCircleRadiustWidget : public QWidget
{
    Q_OBJECT
public:
    /** Constructor */
    AdjustCircleRadiustWidget(int min, int max, int start_with,
                              QtWidgetMicrograph *_qtwidgetmicrograph,
                              QWidget *parent = 0, const char *name = 0, int wflags = 0);
private:
    QtWidgetMicrograph *__qtwidgetmicrograph;
    QScrollBar        *__scroll_radius;
    QLabel            *__label_radius;

private slots:
    void scrollValueChanged(int);
};
#endif
