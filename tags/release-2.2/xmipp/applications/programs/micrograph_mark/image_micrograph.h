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

#ifndef __QT_IMAGE_MICROGRAPH_HH__
#define __QT_IMAGE_MICROGRAPH_HH__

#define MARK_CIRCLE 0
#define MARK_SQUARE 1

#include "image.h"

#ifdef QT3_SUPPORT
//Added by qt3to4:
#include <QResizeEvent>
#include <QMouseEvent>
#endif

/* Forward declarations ---------------------------------------------------- */
class Micrograph;

/* Widget for the micrograph ----------------------------------------------- */
class QtImageMicrograph : public QtImage
{
    Q_OBJECT
public:
    float   __ellipse_radius;
    int     __ellipse_type;
private:
    bool    __pressed;
    int     __movingMark;
    bool    __tilted;

public:
    // Constructor
    QtImageMicrograph(QWidget *_parent = 0, const char *_name = 0, Qt::WFlags _f = 0);

    // Set the micrograph
    void setMicrograph(Micrograph *_m);
    void setTilted()
    {
        __tilted = TRUE;
    }
    bool isTilted()
    {
        return __tilted;
    }
    void movingMark(int _coord)
    {
        __movingMark = _coord;
    }
    void drawEllipse(int _x, int _y, int _color, float _ellipse_radius = 5.0, int _type = MARK_CIRCLE);
    void drawLastEllipse(int _x, int _y, int _color, float _ellipse_radius = 5.0, int _type = MARK_CIRCLE);

protected:
    void loadImage();
    void loadSymbols();
    void applyFilters(QtFiltersController *_f, QImage *_img)
    {
        _f->applyFilters(_img);
    }

public slots:
    void slotDeleteMarkOther(int _coord);
    void slotChangeFamilyOther(int _coord, int _f);
    void slotZoomIn();
    void slotZoomOut();

signals:
    void signalSetWidthHeight(int _w, int _h);
    void signalAddCoordOther(int _mX, int _mY, int _f);
    void signalDeleteMarkOther(int _coord);
    void signalChangeFamilyOther(int _coord, int _f);
    void signalRecalculateTiltMatrix();

protected:
    void resizeEvent(QResizeEvent *e);
    void mousePressEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);
    void mouseMoveEvent(QMouseEvent *e);
    void changeProperties(int mX, int mY);
};

#endif
