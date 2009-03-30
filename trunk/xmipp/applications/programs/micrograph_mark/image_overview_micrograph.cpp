/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Carlos Manzanares       (cmanzana@cnb.uam.es)
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

#include "image_overview_micrograph.h"
#include "color_label.h"

#include <data/micrograph.h>
#include <data/xvsmooth.h>

#include <qpainter.h>

#ifdef QT3_SUPPORT
//Added by qt3to4:
#include <QMouseEvent>
#endif

/* Coordiantes transformations --------------------------------------------- */
void QtImageOverviewMicrograph::micrographToOverview(const int _x, const int _y,
        int &_rx, int &_ry)
{
    double rx, ry;
    exact_micrographToOverview(_x, _y, rx, ry);
    _rx = (int)(rx);
    _ry = (int)(ry);
}

void QtImageOverviewMicrograph::overviewToMicrograph(int _x, int _y,
        int &_rx, int &_ry)
{
    double rx, ry;
    exact_overviewToMicrograph(_x, _y, rx, ry);
    _rx = (int)(rx);
    _ry = (int)(ry);
}

void QtImageOverviewMicrograph::exact_micrographToOverview(const int _x, const int _y,
        double &_rx, double &_ry)
{
    int mMaxX, mMaxY;

    if (getMicrograph() != NULL) getMicrograph()->size(mMaxX, mMaxY);
    else return;

    double ratioX = (double)mMaxX / image()->width(),
                    ratioY = (double)mMaxY / image()->height();

    _rx = (_x / ratioX);
    _ry = (_y / ratioY);
}

void QtImageOverviewMicrograph::exact_overviewToMicrograph(int _x, int _y,
        double &_rx, double &_ry)
{
    int mMaxX, mMaxY;

    if (getMicrograph() != NULL) getMicrograph()->size(mMaxX, mMaxY);
    else return;

    double ratioX = (double)mMaxX / image()->width(),
                    ratioY = (double)mMaxY / image()->height();

    _rx = (_x * ratioX);
    _ry = (_y * ratioY);
}

/* Constructor ------------------------------------------------------------- */
QtImageOverviewMicrograph::QtImageOverviewMicrograph(QWidget *_parent,
        const char *_name,
        Qt::WFlags _f) :
        QtImage(_parent, _name, _f)
{
    __x = 0;
    __y = 0;
    __w = 0;
    __h = 0;
    __minCost = 0;
    __axis_ang = 0;
    __enable_axis = FALSE;
    __x0_crop = __y0_crop = __xF_crop = __yF_crop = -1;
}

/* Load Image -------------------------------------------------------------- */
void QtImageOverviewMicrograph::loadImage()
{
    int mMaxX, mMaxY, mX, mY;
    double emX, emY;

    if (getMicrograph() != NULL) getMicrograph()->size(mMaxX, mMaxY);
    else return;

    if (mMaxX >= image()->width() && mMaxY >= image()->height() &&
        getMicrograph()->depth() == 8)
    {
        byte rgb[256];
        for (int i = 0; i < 256; i++) rgb[i] = i;
        byte *result = SmoothResize((byte *)(getMicrograph()->array8()),
                                    mMaxX, mMaxY, image()->width(), image()->height(),
                                    rgb, rgb, rgb, rgb, rgb, rgb, 256);
        byte *ptr = result;
        double a,b;
        getMicrograph()->getLinearTransformatioVal8(a,b);
        for (int y = 0; y < image()->height(); y++)
            for (int x = 0; x < image()->width(); x++)
            {
                byte proposedVal=(*ptr++);
                setPixel(x, y, CLIP(a*proposedVal+b,0,255));
            }
        free(result);
    }
    else
    {
        // Apply bilinear interpolation
        for (int y = 0; y < image()->height(); y++)
            for (int x = 0; x < image()->width(); x++)
                if (getMicrograph() != NULL)
                {
                    exact_overviewToMicrograph(x, y, emX, emY);
                    double val = 0;
                    if (emX >= 0 && emX < mMaxX - 1 && emY >= 0 && emY < mMaxY - 1)
                    {
                        double wx = emX - (int)emX;
                        double wy = emY - (int)emY;
                        int    mX1 = (int)emX, mX2 = mX1 + 1;
                        int    mY1 = (int)emY, mY2 = mY1 + 1;
                        val += (1 - wy) * (1 - wx) * getMicrograph()->val8(mX1, mY1) +
                               (1 - wy) *   wx * getMicrograph()->val8(mX2, mY1) +
                               wy * (1 - wx) * getMicrograph()->val8(mX1, mY2) +
                               wy *   wx * getMicrograph()->val8(mX2, mY2);
                    }
                    setPixel(x, y, (unsigned int)val);
                }
                else setPixel(x, y, 0);
    }
}

/* Draw ellipse ------------------------------------------------------------ */
void QtImageOverviewMicrograph::drawEllipse(int _x, int _y, int _color)
{
    int mX, mY;
    micrographToOverview(_x, _y, mX, mY);
    if ((mX > 0 && mX < image()->width()) &&
        (mY > 0 && mY < image()->height()))
    {
        __paint->setPen(__col.col(_color));
        __paint->drawRect(mX - 2, mY - 2, 4, 4);
    }
}

/* Draw axis --------------------------------------------------------------- */
void QtImageOverviewMicrograph::draw_axis(double _ang)
{
    if (!__enable_axis) return;
    int h = image()->height();
    int w = image()->width();
    double limit_angle = RAD2DEG(atan2((double)w, (double)h));
    __axis_ang = _ang;

    // Set angle between -90 and 90
    if (_ang > 90 && _ang <= 270)      _ang -= 180;
    else if (_ang > 270 && _ang < 360) _ang -= 360;

    // Vertical or horizontal lines?
    int x1, y1, x2, y2;
    if (ABS(_ang) < limit_angle)
    {
        // Vertical
        double tg = tan(DEG2RAD(_ang));
        y1 = 0;
        x1 = (int)(w / 2 + tg * h / 2);
        y2 = (int)(h - 1);
        x2 = (int)(w / 2 - tg * h / 2);
    }
    else
    {
        // Horizontal
        double cotg = cos(DEG2RAD(_ang)) / sin(DEG2RAD(_ang));
        x1 = 0;
        y1 = (int)(h / 2 + cotg * w / 2);
        x2 = (int)(w - 1);
        y2 = (int)(h / 2 - cotg * w / 2);
    }

    __paint->setPen(__col.col(0));
    __paint->drawLine(x1, y1, x2, y2);
}

/* Load Symbols ------------------------------------------------------------ */
void QtImageOverviewMicrograph::loadSymbols()
{
    int x, y, w, h;

    if (getMicrograph() == NULL) return;

    micrographToOverview(__x, __y, x, y);
    micrographToOverview((int)(__w * __zoom), (int)(__h * __zoom), w, h);
    if ((x + w) > image()->width())  w = image()->width() - x;
    if ((y + h) > image()->height()) h = image()->height() - y;

    __paint->setPen(Qt::yellow);
    __paint->drawRect(x , y, w, h);

    for (int i = 0; i < getMicrograph()->ParticleNo(); i++)
    {
        if (!getMicrograph()->coord(i).valid ||
            getMicrograph()->coord(i).cost<__minCost) continue;
        drawEllipse(getMicrograph()->coord(i).X,
                    getMicrograph()->coord(i).Y, getMicrograph()->coord(i).label);
    }

    // Draw also the crop area if necessary
    if (__x0_crop != -1 && __y0_crop != -1 && __xF_crop != -1 && __yF_crop != -1)
    {
        int x0_crop, y0_crop, xF_crop, yF_crop;
        micrographToOverview(__x0_crop, __y0_crop, x0_crop, y0_crop);
        micrographToOverview(__xF_crop, __yF_crop, xF_crop, yF_crop);
        __paint->setPen(Qt::blue);
        __paint->drawRect(x0_crop , y0_crop, xF_crop - x0_crop, yF_crop - y0_crop);
    }
}

void QtImageOverviewMicrograph::mouseMoveEvent(QMouseEvent *e)
{
    int mX, mY, x, y;

    if (getMicrograph() == NULL) return;
    if (e->button() != Qt::LeftButton) return;

    x = e->pos().x();
    y = e->pos().y();

    overviewToMicrograph(x, y, __x, __y);

    // But this now must be the center
    int Xdim, Ydim;
    getMicrograph()->size(Xdim, Ydim);
    __x = XMIPP_MAX(0, (int)(__x - __w / 2.0));
    __y = XMIPP_MAX(0, (int)(__y - __h / 2.0));
    __x = XMIPP_MIN((int)(Xdim - 1 - __w * __zoom), __x);
    __y = XMIPP_MIN((int)(Ydim - 1 - __h * __zoom), __y);

    emit signalSetCoords(__x, __y);
    emit signalRepaint();
    emit signalActualizeOtherOverview(__x, __y);
}

void QtImageOverviewMicrograph::mouseReleaseEvent(QMouseEvent *e)
{
    if (e->button() != Qt::LeftButton) return;
    mouseMoveEvent(e);
}

void QtImageOverviewMicrograph::resizeEvent(QResizeEvent *e)
{
    QtImage::resizeEvent(e);
    if (getMicrograph() == NULL) return;
    emit signalRepaint();
}

void QtImageOverviewMicrograph::slotSetWidthHeight(int _w, int _h)
{
    __w = _w;
    __h = _h;
}

void QtImageOverviewMicrograph::slotActualizeOtherOverview(int _x, int _y)
{
    bool inside_current_window = true;
    if (_x < __x) inside_current_window = false;
    else if (_x > __x + __w*__zoom) inside_current_window = false;
    else if (_y < __y) inside_current_window = false;
    else if (_y > __y + __h*__zoom) inside_current_window = false;

    if (!inside_current_window)
    {
        __x = _x;
        __y = _y;

        // This now should be the center
        int Xdim, Ydim;
        getMicrograph()->size(Xdim, Ydim);
        __x = XMIPP_MAX(0, (int)(__x - __w / 2.0));
        __y = XMIPP_MAX(0, (int)(__y - __h / 2.0));
        __x = XMIPP_MIN((int)(Xdim - 1 - __w * __zoom), __x);
        __y = XMIPP_MIN((int)(Ydim - 1 - __h * __zoom), __y);

        emit signalSetCoords(__x, __y);
    }
    emit signalRepaint();
}

void QtImageOverviewMicrograph::setMicrograph(Micrograph *_m)
{
    QtImage::setMicrograph(_m);

    if (getMicrograph() != NULL)
    {
        int mMaxX, mMaxY;
        getMicrograph()->size(mMaxX, mMaxY);
        double aspect = (double)mMaxX / mMaxY;
        setGeometry(geometry().x(), geometry().y(),
                    geometry().width(), (int)(geometry().width() * aspect));
    }
}

/* Crop area --------------------------------------------------------------- */
void QtImageOverviewMicrograph::init_crop_area()
{
    int Xdim, Ydim;
    getMicrograph()->size(Xdim, Ydim);
    __x0_crop = ROUND(0.25 * Xdim);
    __y0_crop = ROUND(0.25 * Ydim);
    __xF_crop = ROUND(0.75 * Xdim);
    __yF_crop = ROUND(0.75 * Ydim);
    emit signalRepaint();
}

void QtImageOverviewMicrograph::finish_crop_area()
{
    __x0_crop = __y0_crop = __xF_crop = __yF_crop = -1;
    emit signalRepaint();
}

void QtImageOverviewMicrograph::slotDrawCropArea(std::vector<int> value)
{
    __x0_crop = value[0];
    __y0_crop = value[1];
    __xF_crop = value[2];
    __yF_crop = value[3];
    emit signalRepaint();
}
