/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Carlos Manzanares       (cmanzana@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 *
 *****************************************************************************/

#ifndef __QT_IMAGE_OVERVIEW_MICROGRAPH_HH__
#define __QT_IMAGE_OVERVIEW_MICROGRAPH_HH__

#include "image.h"

#ifdef QT3_SUPPORT
//Added by qt3to4:
#include <QMouseEvent>
#include <QResizeEvent>
#endif

/* Forward declarations ---------------------------------------------------- */
class Micrograph;

/* Widget for the micrograph ----------------------------------------------- */
class QtImageOverviewMicrograph : public QtImage
{
    Q_OBJECT

private:
    double  __w, __h;
    int     __x0_crop, __y0_crop, __xF_crop, __yF_crop;

public:
    // Minimum cost for drawing
    float __minCost;

    // Constructor
    QtImageOverviewMicrograph(QWidget *_parent = 0, const char *_name = 0, Qt::WFlags _f = 0);

    // Set the micrograph
    void setMicrograph(Micrograph *_m);
    void drawEllipse(int _x, int _y, int _color);
    void draw_axis(double _ang);

    // Crop area
    void init_crop_area();
    void finish_crop_area();

protected:
    // Coordinate transformations
    void micrographToOverview(const int _x, const int _y,
                              int &_rx, int &_ry);
    void overviewToMicrograph(const int _x, const int _y,
                              int &_rx, int &_ry);
    void exact_micrographToOverview(const int _x, const int _y,
                                    double &_rx, double &_ry);
    void exact_overviewToMicrograph(const int _x, const int _y,
                                    double &_rx, double &_ry);
    void loadImage();
    void loadSymbols();
    void applyFilters(QtFiltersController *_f, QImage *_img)
    {}

public slots:
    void slotSetWidthHeight(int _w, int _h);
    void slotActualizeOtherOverview(int _x, int _y);
    void slotDrawCropArea(std::vector<int> value);

signals:
    void signalActualizeOtherOverview(int _x, int _y);

protected:
    void mouseMoveEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);
    void resizeEvent(QResizeEvent *e);
};

#endif
