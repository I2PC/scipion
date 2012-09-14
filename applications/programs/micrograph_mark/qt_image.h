/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Carlos Manzanares       (cmanzana@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Copyright (c) 2000 , CSIC .
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 *
 *****************************************************************************/

#ifndef __QT_IMAGE_HH__
#define __QT_IMAGE_HH__

#include <qwidget.h>
#include <qimage.h>

#ifdef QT3_SUPPORT
//Added by qt3to4:
#include <QResizeEvent>
#include <QPaintEvent>
#endif

#include "color_label.h"
#include "filters_controller.h"

/* Forward declarations ---------------------------------------------------- */
class QImage;
class Micrograph;
class QtWidgetMicrograph;

/* Widget for the micrograph ----------------------------------------------- */
class QtImage : public QWidget
{
    Q_OBJECT

private:
    QImage              *__img;
    Micrograph          *__m;
    QtFiltersController *__filtersController;
    QtWidgetMicrograph  *__wm;

    int                 __mingray;
    int                 __maxgray;
    float               __gamma;
public:
    // Constructor
    QtImage(QWidget *_parent = 0, const char *_name = 0, Qt::WFlags _f = 0);
    ~QtImage();

    // Set the micrograph
    void setMicrograph(Micrograph *_m);

    // Set the widget micrograph
    void setWidgetMicrograph(QtWidgetMicrograph *_wm);

    // Set the filters controller
    void setFiltersController(QtFiltersController *_f);

    // Get coordinates
    void getCoordinates(int &_x, int &_y)
    {
        _x = (int)(__x * __zoom);
        _y = (int)(__y * __zoom);
    }
    void setCoordinates(int _x, int _y)
    {
        emit signalSetCoords(_x, _y);
    }

    // Tilt axis
    void enableAxis()
    {
        __enable_axis = TRUE;
    }
    virtual void draw_axis(double _ang)
    {};

    // Change Image contrast
    void changeContrast(int _mingray, int _maxgray, float _gamma);
protected:
    int               __x, __y;
    static double     __zoom;
    int               __activeFamily;
    QtColorLabel      __col;
    QPainter         *__paint;
    double            __axis_ang;
    bool              __enable_axis;

    // Coordinates transformations
    void micrographToImage(const int _x, const int _y, int &_rx, int &_ry);
    void imageToMicrograph(const int _x, const int _y, int &_rx, int &_ry);
    void exact_micrographToImage(const int _x, const int _y, double &_rx, double &_ry);
    void exact_imageToMicrograph(const int _x, const int _y, double &_rx, double &_ry);
    virtual void loadImage() = 0;
    virtual void loadSymbols() = 0;
    virtual void applyFilters(QtFiltersController *_f, QImage *_img) = 0;

    // Get the micrograph
    Micrograph *getMicrograph()
    {
        return(__m);
    }
    // Get the widget micrograph
    QtWidgetMicrograph *getWidgetMicrograph()
    {
        return(__wm);
    }
    // Get the image handler
    QImage     *image()
    {
        return(__img);
    }

    // Set pixel
    void setPixel(int _x, int _y, int _value);

    void resizeEvent(QResizeEvent *);
    void paintEvent(QPaintEvent *);
public slots:
    void slotRepaint()
    {
        repaint(false);
    }
    void slotActiveFamily(int _f);
    void slotSetCoords(int _x, int _y);

signals:
    void signalRepaint();
    void signalSetCoords(int _x, int _y);
};

#endif
