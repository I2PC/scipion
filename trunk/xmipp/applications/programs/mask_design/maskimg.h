/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Copyright (c) 2001 , CSIC .
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

#ifndef MASKIMG_H
#define MASKIMG_H

#include <qwidget.h>
#include <qimage.h>
#include <qpainter.h>

#ifdef QT3_SUPPORT
//Added by qt3to4:
#include <QPaintEvent>
#include <QResizeEvent>
#include <QPixmap>
#include <QLabel>
#include <QMouseEvent>
#include <Q3PopupMenu>
#include <QKeyEvent>
#endif

#include <data/selfile.h>
#include <data/image.h>

///  Ways the training set can be used
typedef enum {CIRCLE, ELLIPSE,  RECTANGLE, SQUARE, RING, ELLIPRING,
              RECTFRAME, SQUFRAME, POLYGON
             } maskType;

/**@name maskFigure*/
//@{
/**
 * This class implements an abstract mask figure.
 * It is used for generating geometric masks.
 */
class maskFigure: public QPainter
{
public:
    /** maskFigure default constructor
    */

    /** maskFigure constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    */
    maskFigure(QPaintDevice* _paintDevice, QPixmap* _pixmap, int _width, int _height):
        QPainter(_paintDevice), pixmap(_pixmap)
    {
        width = _width;
        height = _height;
        rcx = width / 2;
        ocx = rcx;
        rcy = height / 2;
        ocy = rcy;
    }

    maskFigure(int _cx = 5, int _cy = 5): rcx(_cx), rcy(_cy)
    {};

    /** maskFigure destructor
    */
    ~maskFigure()
    {};

    /** Gets a constant reference to the X-coordinate of the figure's center
    */
    const int& cx() const
    {
        return rcx;
    };

    /** Gets a constant reference to the Y-coordinate of the figure's center
    */
    const int& cy() const
    {
        return rcy;
    };

    /** Gets a reference to the X-coordinate of the figure's center
    */
    int& cx()
    {
        return rcx;
    };

    /** Gets a reference to the Y-coordinate of the figure's center
    */
    int& cy()
    {
        return rcy;
    };

    /** Increases the area in the case of the ring
    */
    virtual bool increaseArea()
    {
        return true;
    }

    /** Decreases the area in the case of the ring
    */
    virtual bool decreaseArea()
    {
        return true;
    }

    /** Moves the figure
    *   Parameter: _cx: new X-coordinate
    *   Parameter: _cy: new Y-coordinate
    */
    virtual bool move(int _cx, int _cy) = 0;

    /** Moves the figure to the left one pixel
    */
    virtual bool moveLeft() = 0;

    /** Moves the figure up one pixel
    */
    virtual bool moveUp() = 0;

    /** Moves the figure to the right one pixel
    */
    virtual bool moveRight() = 0;

    /** Moves the figure down one pixel
    */
    virtual bool moveDown() = 0;

    /** Decreases the width one pixel
    */
    virtual bool decreaseWidth() = 0;

    /** Decreases the height one pixel
    */
    virtual bool decreaseHeight() = 0;

    /** Increases the width one pixel
    */
    virtual bool increaseWidth() = 0;

    /** Increases the height one pixel
    */
    virtual bool increaseHeight() = 0;


    /** Resizes the figure
    *   Parameter: _width: new width
    *   Parameter: _height: new height
    */
    virtual void resize(int _width, int _height)
    {
        rcx = (int)((float)(rcx * _width) / (float)width);
        rcy = (int)((float)(rcy * _width) / (float)width);
        ocx = rcx;
        ocy = rcy;
    };


    /** Paints the figure (virtual method)
    *
    */
    virtual void paint()
    {
        QBrush brush(Qt::NoBrush);    // void brush
        QPen myPen(Qt::red);
        setPen(myPen);
        setBrush(brush);
    };

    /** this method returns true if the point passed as parameter is inside teh figure (virtual method)
    *   Parameter: _cx: X xoordinate of the point to test
    *   Parameter: _cy: X xoordinate of the point to test
    */
    virtual bool isIn(int _cx, int _cy) = 0;

    /** This method prints the class parameters.
    */
    virtual void print(float scalex = 1.0, float scaley = 1.0) = 0;


protected:
    int rcx, rcy;     // actual center of the figure
    int ocx, ocy;   // old center of the figure
    int height, width;   // Height and Width of the figure
    QPaintDevice* paintDevice;         // Painter device
    QPixmap* pixmap;   // Pixmap to paint on
};



/**@name maskEllipse*/
//@{
/**
 * This class implements circle figure to b used as mask.
 * It is used for generating geometric masks.
 */
class maskEllipse: public maskFigure
{
public:

    /** maskEllipse constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    */
    maskEllipse(QPaintDevice* _paintDevice, QPixmap* _pixmap, int _width, int _height): maskFigure(_paintDevice, _pixmap, _width, _height)
    {
        rxradius = height / 4;
        ryradius = height / 4;
        if (rxradius < 1) rxradius = 20;
        if (ryradius < 1) ryradius = 20;
        oxradius = rxradius;
        oyradius = ryradius;
    };


    /** maskEllipse constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    *   Parameter: _xradius: radius of the ellipse
    *   Parameter: _yradius: radius of the ellipse
    */
    maskEllipse(int _cx = 0, int _cy = 0, int _xradius = 1, int _yradius = 1): maskFigure(_cx, _cy), rxradius(_xradius), ryradius(_yradius)
    {};

    /** maskCircle destructor
    */
    ~maskEllipse()
    {};

    /** Gets a constant reference to the x radius of the ellipse
    */
    const int& xradius() const
    {
        return rxradius;
    };

    /** Gets a constant reference to the y radius of the ellipse
    */
    const int& yradius() const
    {
        return ryradius;
    };

    /** Gets reference to the x radius of the ellipse
    */
    int& xradius()
    {
        return rxradius;
    };

    /** Gets a reference to the y radius of the ellipse
    */
    int& yradius()
    {
        return ryradius;
    };


    /** Moves the figure
    *   Parameter: _cx: new X-coordinate
    *   Parameter: _cy: new Y-coordinate
    */
    virtual bool move(int _cx, int _cy)
    {
        if (((rcx - rxradius) > 0) &&
            ((rcy - ryradius) > 0) &&
            ((rcx + rxradius) < width) &&
            ((rcy + ryradius) < height))
        {
            oxradius = rxradius;
            oyradius = ryradius;
            ocx = rcx;
            ocy = rcy;
            rcx = _cx;
            rcy = _cy;
            paint();
            return true;
        }
        else return false;
    };

    /** Moves the figure to the left one pixel
    */
    virtual bool moveLeft()
    {
        if ((rcx - rxradius) > 0)
        {
            oxradius = rxradius;
            oyradius = ryradius;
            ocx = rcx;
            ocy = rcy;
            rcx--;
            paint();
            return true;
        }
        else return false;
    };


    /** Moves the figure up one pixel
    */
    virtual bool moveUp()
    {
        if ((rcy - ryradius) > 0)
        {
            oxradius = rxradius;
            oyradius = ryradius;
            ocx = rcx;
            ocy = rcy;
            rcy--;
            paint();
            return true;
        }
        else return false;
    };


    /** Moves the figure to the right one pixel
    */
    virtual bool moveRight()
    {
        if ((rcx + rxradius) < width)
        {
            oxradius = rxradius;
            oyradius = ryradius;
            ocx = rcx;
            ocy = rcy;
            rcx++;
            paint();
            return true;
        }
        else return false;
    };

    /** Moves the figure down one pixel
    */
    virtual bool moveDown()
    {
        if ((rcy + ryradius) < height)
        {
            oxradius = rxradius;
            oyradius = ryradius;
            ocx = rcx;
            ocy = rcy;
            rcy++;
            paint();
            return true;
        }
        else return false;
    };



    /** Decreases the width one pixel
    */
    virtual bool decreaseWidth()
    {
        if (rxradius > 1)
        {
            oxradius = rxradius;
            oyradius = ryradius;
            ocx = rcx;
            ocy = rcy;
            rxradius--;
            paint();
            return true;
        }
        else return false;
    };

    /** Decreases the height one pixel
    */
    virtual bool decreaseHeight()
    {
        if (ryradius > 1)
        {
            oxradius = rxradius;
            oyradius = ryradius;
            ocx = rcx;
            ocy = rcy;
            ryradius--;
            paint();
            return true;
        }
        else return false;
    };

    /** Increases the width one pixel
    */
    virtual bool increaseWidth()
    {
        if (((rcx + rxradius) < width) && ((rcy + ryradius) < height) && ((rcx - rxradius) > 0) && ((rcy - ryradius) > 0))
        {
            oxradius = rxradius;
            oyradius = ryradius;
            ocx = rcx;
            ocy = rcy;
            rxradius++;
            paint();
            return true;
        }
        else return false;
    };

    /** Increases the height one pixel
    */
    virtual bool increaseHeight()
    {
        if (((rcx + rxradius) < width) && ((rcy + ryradius) < height) && ((rcx - rxradius) > 0) && ((rcy - ryradius) > 0))
        {
            oxradius = rxradius;
            oyradius = ryradius;
            ocx = rcx;
            ocy = rcy;
            ryradius++;
            paint();
            return true;
        }
        else return false;
    };


    /** Resizes the figure
    *   Parameter: _width: new width
    *   Parameter: _height: new height
    */
    virtual void resize(int _width, int _height)
    {
        maskFigure::resize(_width, _height);
        rxradius = (int)((float)(rxradius * _width) / (float)width);
        ryradius = (int)((float)(ryradius * _width) / (float)width);
        oxradius = rxradius;
        oyradius = ryradius;
        height = _height;
        width = _width;
        paint();
    };


    /** Paints the figure (virtual method)
    *
    */
    virtual void paint()
    {
        maskFigure::paint();
        QPoint point(ocx - oxradius, ocy - oyradius);
        QRect rect(ocx - oxradius, ocy - oyradius, oxradius*2, oyradius*2);
        drawPixmap(point, (*pixmap), rect);
        drawEllipse(rcx - rxradius, rcy - ryradius, rxradius*2, ryradius*2);
    }


    /** this method returns true if the point passed as parameter is inside the circle
    *   Parameter: _cx: X xoordinate of the point to test
    *   Parameter: _cy: X xoordinate of the point to test
    */
    virtual bool isIn(int _cx, int _cy)
    {
        float X = (float)(_cx - rcx) / (float)rxradius;
        float Y = (float)(_cy - rcy) / (float)ryradius;
        if ((X*X + Y*Y) < 1)
            return true;
        else
            return false;
    };

    /** Prints mask parameters (virtual)*/
    virtual void print(float scalex, float scaley)
    {
        std::cout << "(Origin_x,Origin_y): (" << rcx / scalex << ", "
                  << rcy / scaley << ")" << std::endl;

        std::cout << "(Radius_x,Radius_y): (" << rxradius / scalex << ", "
                  << ryradius / scaley << ")" << std::endl;
    }

protected:
    int rxradius, oxradius;  // real radius of the figure
    int ryradius, oyradius;  // real radius of the figure
};





/**@name maskCircle*/
//@{
/**
 * This class implements a circle figure to b used as mask.
 * It is used for generating geometric masks.
 */
class maskCircle: public maskEllipse
{
public:

    /** maskCircle constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    */
    maskCircle(QPaintDevice* _paintDevice, QPixmap* _pixmap, int _width, int _height): maskEllipse(_paintDevice, _pixmap, _width, _height)
    {};


    /** maskCircle constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    *   Parameter: _radius: radius of the circle
    */
    maskCircle(int _cx = 0, int _cy = 0, int _radius = 1): maskEllipse(_cx, _cy, _radius, _radius)
    {};

    /** maskCircle destructor
    */
    ~maskCircle()
    {};

    /** Gets a constant reference to the radius of the circle
    */
    const int& radius() const
    {
        return rxradius;
    };

    /** Sets the radius of the circle
    */
    void setRadius(int _radius)
    {
        oxradius = rxradius;
        oyradius = ryradius;
        rxradius = _radius;
        ryradius = _radius;
    };

    /** Decreases the width one pixel
    */
    virtual bool decreaseWidth()
    {
        maskEllipse::decreaseWidth();
        return maskEllipse::decreaseHeight();
    };

    /** Decreases the height one pixel
    */
    virtual bool decreaseHeight()
    {
        maskEllipse::decreaseWidth();
        return maskEllipse::decreaseHeight();
    };

    /** Increases the width one pixel
    */
    virtual bool increaseWidth()
    {
        maskEllipse::increaseWidth();
        return maskEllipse::increaseHeight();
    };

    /** Increases the height one pixel
    */
    virtual bool increaseHeight()
    {
        maskEllipse::increaseWidth();
        return maskEllipse::increaseHeight();
    };

};


/**@name maskRectangle*/
//@{
/**
 * This class implements circle figure to b used as mask.
 * It is used for generating geometric masks.
 */
class maskRectangle: public maskFigure
{
public:

    /** maskRectangle constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    */
    maskRectangle(QPaintDevice* _paintDevice, QPixmap* _pixmap, int _width, int _height): maskFigure(_paintDevice, _pixmap, _width, _height)
    {
        rrwidth = height / 4;
        rrheight = height / 4;
        if (rrwidth < 1) rrwidth = 20;
        if (rrheight < 1) rrheight = 20;
        orwidth = rrwidth;
        orheight = rrheight;
    };


    /** maskRectangle constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    *   Parameter: _rwidth: width of the rectangle
    *   Parameter: _rheight: height of the rectangle
    */
    maskRectangle(int _cx = 0, int _cy = 0, int _rwidth = 1, int _rheight = 1): maskFigure(_cx, _cy), rrwidth(_rwidth), rrheight(_rheight)
    {};

    /** maskRectangle destructor
    */
    ~maskRectangle()
    {};

    /** Gets a constant reference to the width of the rectangle
    */
    const int& rwidth() const
    {
        return rrwidth;
    };

    /** Gets a constant reference to the height of the rectangle
    */
    const int& rheight() const
    {
        return rrheight;
    };

    /** Gets a reference to the width of the rectangle
    */
    int& rwidth()
    {
        return rrwidth;
    };

    /** Gets a reference to the height of the rectangle
    */
    int& rheight()
    {
        return rrheight;
    };

    /** Moves the figure
    *   Parameter: _cx: new X-coordinate
    *   Parameter: _cy: new Y-coordinate
    */
    virtual bool move(int _cx, int _cy)
    {
        if (((rcx - rrwidth) > 0) &&
            ((rcy - rrheight) > 0) &&
            ((rcx + rrwidth) < width) &&
            ((rcy + rrheight) < height))
        {
            orwidth = rrwidth;
            orheight = rrheight;
            ocx = rcx;
            ocy = rcy;
            rcx = _cx;
            rcy = _cy;
            paint();
            return true;
        }
        else return false;
    };

    /** Moves the figure to the left one pixel
    */
    virtual bool moveLeft()
    {
        if ((rcx - rrwidth) > 0)
        {
            orwidth = rrwidth;
            orheight = rrheight;
            ocx = rcx;
            ocy = rcy;
            rcx--;
            paint();
            return true;
        }
        else return false;
    };


    /** Moves the figure up one pixel
    */
    virtual bool moveUp()
    {
        if ((rcy - rrheight) > 0)
        {
            orwidth = rrwidth;
            orheight = rrheight;
            ocx = rcx;
            ocy = rcy;
            rcy--;
            paint();
            return true;
        }
        else return false;
    };


    /** Moves the figure to the right one pixel
    */
    virtual bool moveRight()
    {
        if ((rcx + rrwidth) < width)
        {
            orwidth = rrwidth;
            orheight = rrheight;
            ocx = rcx;
            ocy = rcy;
            rcx++;
            paint();
            return true;
        }
        else return false;
    };

    /** Moves the figure down one pixel
    */
    virtual bool moveDown()
    {
        if ((rcy + rrheight) < height)
        {
            orwidth = rrwidth;
            orheight = rrheight;
            ocx = rcx;
            ocy = rcy;
            rcy++;
            paint();
            return true;
        }
        else return false;
    };


    /** Decreases the width one pixel
    */
    virtual bool decreaseWidth()
    {
        if (rrwidth > 1)
        {
            orwidth = rrwidth;
            orheight = rrheight;
            ocx = rcx;
            ocy = rcy;
            rrwidth--;
            paint();
            return true;
        }
        else return false;
    };

    /** Decreases the height one pixel
    */
    virtual bool decreaseHeight()
    {
        if (rrheight > 1)
        {
            orwidth = rrwidth;
            orheight = rrheight;
            ocx = rcx;
            ocy = rcy;
            rrheight--;
            paint();
            return true;
        }
        else return false;
    };

    /** Increases the width one pixel
    */
    virtual bool increaseWidth()
    {
        if (((rcx + rrwidth) < width) && ((rcy + rrheight) < height) && ((rcx - rrwidth) > 0) && ((rcy - rrheight) > 0))
        {
            orwidth = rrwidth;
            orheight = rrheight;
            ocx = rcx;
            ocy = rcy;
            rrwidth++;
            paint();
            return true;
        }
        else return false;
    };

    /** Increases the height one pixel
    */
    virtual bool increaseHeight()
    {
        if (((rcx + rrwidth) < width) && ((rcy + rrheight) < height) && ((rcx - rrwidth) > 0) && ((rcy - rrheight) > 0))
        {
            orwidth = rrwidth;
            orheight = rrheight;
            ocx = rcx;
            ocy = rcy;
            rrheight++;
            paint();
            return true;
        }
        else return false;
    };


    /** Resizes the figure
    *   Parameter: _width: new width
    *   Parameter: _height: new height
    */
    virtual void resize(int _width, int _height)
    {
        maskFigure::resize(_width, _height);
        rrwidth = (int)((float)(rrwidth * _width) / (float)width);
        rrheight = (int)((float)(rrheight * _width) / (float)width);
        orwidth = rrwidth;
        orheight = rrheight;
        height = _height;
        width = _width;
        paint();
    };


    /** Paints the figure (virtual method)
    *
    */
    virtual void paint()
    {
        maskFigure::paint();
        QPoint point(ocx - orwidth, ocy - orheight);
        QRect rect(ocx - orwidth, ocy - orheight, orwidth*2, orheight*2);
        drawPixmap(point, (*pixmap), rect);
        drawRect(rcx - rrwidth, rcy - rrheight, rrwidth*2, rrheight*2);
    }


    /** this method returns true if the point passed as parameter is inside the circle
    *   Parameter: _cx: X xoordinate of the point to test
    *   Parameter: _cy: X xoordinate of the point to test
    */
    virtual bool isIn(int _cx, int _cy)
    {
        if ((_cx > (rcx - rrwidth)) && (_cy > (rcy - rrheight)) && (_cx < (rcx + rrwidth)) && (_cy < (rcy + rrheight)))
            return true;
        else
            return false;
    };

    /** Prints mask parameters (virtual)*/
    virtual void print(float scalex, float scaley)
    {
        std::cout << "(Origin_x,Origin_y): (" << rcx / scalex << ", "
                  << rcy / scaley << ")" << std::endl;

        std::cout << "(half_width,half_height): (" << rrwidth / scalex << ", "
                  << rrheight / scaley << ")" << std::endl;
    }

protected:
    int rrwidth, orwidth;    // real width of the rectangle
    int rrheight, orheight;  // real height of the rectangle
};


/**@name maskSquare*/
//@{
/**
 * This class implements a square figure to b used as mask.
 * It is used for generating geometric masks.
 */
class maskSquare: public maskRectangle
{
public:

    /** maskSquare constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    */
    maskSquare(QPaintDevice* _paintDevice, QPixmap* _pixmap, int _width, int _height): maskRectangle(_paintDevice, _pixmap, _width, _height)
    {};


    /** maskSquare constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    *   Parameter: _width: width of the sqaure
    */
    maskSquare(int _cx = 0, int _cy = 0, int _width = 1): maskRectangle(_cx, _cy, _width, _width)
    {};

    /** maskSquare destructor
    */
    ~maskSquare()
    {};

    /** Gets a constant reference to the width of the square
    */
    const int& swidth() const
    {
        return rrwidth;
    };

    /** Gets a reference to the width of the square
    */
    int& swidth()
    {
        return rrwidth;
    };

    /** Decreases the width one pixel
    */
    virtual bool decreaseWidth()
    {
        maskRectangle::decreaseWidth();
        return maskRectangle::decreaseHeight();
    };

    /** Decreases the height one pixel
    */
    virtual bool decreaseHeight()
    {
        maskRectangle::decreaseWidth();
        return maskRectangle::decreaseHeight();
    };

    /** Increases the width one pixel
    */
    virtual bool increaseWidth()
    {
        maskRectangle::increaseWidth();
        return maskRectangle::increaseHeight();
    };

    /** Increases the height one pixel
    */
    virtual bool increaseHeight()
    {
        maskRectangle::increaseWidth();
        return maskRectangle::increaseHeight();
    };

};


/**@name maskEllipRing*/
//@{
/**
 * This class implements a ring figure to b used as mask.
 * It is used for generating geometric masks.
 */
class maskEllipRing: public maskEllipse
{
public:

    /** maskEllipRing constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    */
    maskEllipRing(QPaintDevice* _paintDevice, QPixmap* _pixmap, int _width, int _height): maskEllipse(_paintDevice, _pixmap, _width, _height)
    {
        innerXradius = rxradius - 10;
        innerYradius = rxradius - 10;
    };


    /** maskEllipse constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    *   Parameter: _xradius: radius of the ellipse
    *   Parameter: _yradius: radius of the ellipse
    */
    maskEllipRing(int _cx = 0, int _cy = 0, int _xradius = 1, int _yradius = 1): maskEllipse(_cx, _cy, _xradius, _yradius)
    {
        innerXradius = rxradius - 10;
        innerYradius = rxradius - 10;
    };

    /** maskEllipRing destructor
    */
    ~maskEllipRing()
    {};

    /** Gets a constant reference to the inner radius of the ellipse
    */
    const int& innerXRadius() const
    {
        return innerXradius;
    };

    /** Gets a constant reference to the inner y radius of the ellipse
    */
    const int& innerYRadius() const
    {
        return innerYradius;
    };

    /** Gets a reference to the inner radius of the ellipse
    */
    int& innerXRadius()
    {
        return innerXradius;
    };

    /** Gets a reference to the inner y radius of the ellipse
    */
    int& innerYRadius()
    {
        return innerYradius;
    };


    /** Decreases the area of the ring
    */
    virtual bool decreaseArea()
    {
        if ((rxradius > innerXradius) && (ryradius > innerYradius))
        {
            oxradius = rxradius;
            oyradius = ryradius;
            ocx = rcx;
            ocy = rcy;
            rxradius--;
            ryradius--;
            paint();
            return true;
        }
        else return false;
    };


    /** Increases the area of the ring
    */
    virtual bool increaseArea()
    {
        if (((rcx + rxradius) < width) && ((rcy + ryradius) < height) && ((rcx - rxradius) > 0) && ((rcy - ryradius) > 0))
        {
            oxradius = rxradius;
            oyradius = ryradius;
            ocx = rcx;
            ocy = rcy;
            rxradius++;
            ryradius++;
            paint();
            return true;
        }
        else return false;
    };


    /** Decreases the width one pixel
    */
    virtual bool decreaseWidth()
    {
        if (innerXradius > 1)
        {
            oxradius = rxradius;
            oyradius = ryradius;
            ocx = rcx;
            ocy = rcy;
            innerXradius--;
            rxradius--;
            paint();
            return true;
        }
        else return false;
    };

    /** Decreases the height one pixel
    */
    virtual bool decreaseHeight()
    {
        if (innerYradius > 1)
        {
            oxradius = rxradius;
            oyradius = ryradius;
            ocx = rcx;
            ocy = rcy;
            innerYradius--;
            ryradius--;
            paint();
            return true;
        }
        else return false;
    };

    /** Increases the width one pixel
    */
    virtual bool increaseWidth()
    {
        if (((rcx + rxradius) < width) && ((rcy + ryradius) < height) && ((rcx - rxradius) > 0) && ((rcy - ryradius) > 0))
        {
            oxradius = rxradius;
            oyradius = ryradius;
            ocx = rcx;
            ocy = rcy;
            innerXradius++;
            rxradius++;
            paint();
            return true;
        }
        else return false;
    };

    /** Increases the height one pixel
    */
    virtual bool increaseHeight()
    {
        if (((rcx + rxradius) < width) && ((rcy + ryradius) < height) && ((rcx - rxradius) > 0) && ((rcy - ryradius) > 0))
        {
            oxradius = rxradius;
            oyradius = ryradius;
            ocx = rcx;
            ocy = rcy;
            innerYradius++;
            ryradius++;
            paint();
            return true;
        }
        else return false;
    };


    /** Resizes the figure
    *   Parameter: _width: new width
    *   Parameter: _height: new height
    */
    virtual void resize(int _width, int _height)
    {
        maskFigure::resize(_width, _height);
        rxradius = (int)((float)(rxradius * _width) / (float)width);
        ryradius = (int)((float)(ryradius * _width) / (float)width);
        innerXradius = (int)((float)(innerXradius * _width) / (float)width);
        innerYradius = (int)((float)(innerYradius * _width) / (float)width);
        oxradius = rxradius;
        oyradius = ryradius;
        height = _height;
        width = _width;
        paint();
    };


    /** Paints the figure (virtual method)
    *
    */
    virtual void paint()
    {
        maskFigure::paint();
        QPoint point(ocx - oxradius, ocy - oyradius);
        QRect rect(ocx - oxradius, ocy - oyradius, oxradius*2, oyradius*2);
        drawPixmap(point, (*pixmap), rect);
        drawEllipse(rcx - rxradius, rcy - ryradius, rxradius*2, ryradius*2);
        drawEllipse(rcx - innerXradius, rcy - innerYradius, innerXradius*2, innerYradius*2);
    }


    /** this method returns true if the point passed as parameter is inside the circle
    *   Parameter: _cx: X xoordinate of the point to test
    *   Parameter: _cy: X xoordinate of the point to test
    */
    virtual bool isIn(int _cx, int _cy)
    {
        float Xout = (float)(_cx - rcx) / (float)rxradius;
        float Yout = (float)(_cy - rcy) / (float)ryradius;
        float Xin = (float)(_cx - rcx) / (float)innerXradius;
        float Yin = (float)(_cy - rcy) / (float)innerYradius;
        if ((Xout*Xout + Yout*Yout) < 1)
            if ((Xin*Xin + Yin*Yin) >= 1)
                return true;
        return false;
    };

    /** Prints mask parameters (virtual)*/
    virtual void print(float scalex, float scaley)
    {
        maskEllipse::print(scalex, scaley);
        std::cout << "(innerXradius,innerYradius): (" << innerXradius / scalex << ", "
                  << innerYradius / scaley << ")"
                  << std::endl;
    }


protected:
    int innerXradius;        // inner X radius
    int innerYradius;        // inner Y radius
};


/**@name maskCircleRing*/
//@{
/**
 * This class implements a circle ring figure to b used as mask.
 * It is used for generating geometric masks.
 */
class maskCircleRing: public maskEllipRing
{
public:

    /** maskCircleRing constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    */
    maskCircleRing(QPaintDevice* _paintDevice, QPixmap* _pixmap, int _width, int _height): maskEllipRing(_paintDevice, _pixmap, _width, _height)
    {};


    /** maskCircleRing constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    *   Parameter: _radius: outer radius of the circle
    */
    maskCircleRing(int _cx = 0, int _cy = 0, int _radius = 1): maskEllipRing(_cx, _cy, _radius, _radius)
    {};

    /** maskCircleRing destructor
    */
    ~maskCircleRing()
    {};

    /** Gets a constant reference to the inner radius of the circle
    */
    const int& innerRadius() const
    {
        return innerXradius;
    };


    /** Decreases the width one pixel
    */
    virtual bool decreaseWidth()
    {
        maskEllipRing::decreaseWidth();
        return maskEllipRing::decreaseHeight();
    };

    /** Decreases the height one pixel
    */
    virtual bool decreaseHeight()
    {
        maskEllipRing::decreaseWidth();
        return maskEllipRing::decreaseHeight();
    };

    /** Increases the width one pixel
    */
    virtual bool increaseWidth()
    {
        maskEllipRing::increaseWidth();
        return maskEllipRing::increaseHeight();
    };

    /** Increases the height one pixel
    */
    virtual bool increaseHeight()
    {
        maskEllipRing::increaseWidth();
        return maskEllipRing::increaseHeight();
    };
};




/**@name maskRectFrame*/
//@{
/**
 * This class implements rectangular frame to be used as mask.
 * It is used for generating geometric masks.
 */
class maskRectFrame: public maskRectangle
{
public:

    /** maskRectFrame constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    */
    maskRectFrame(QPaintDevice* _paintDevice, QPixmap* _pixmap, int _width, int _height): maskRectangle(_paintDevice, _pixmap, _width, _height)
    {
        innerwidth = rrwidth - 10;
        innerheight = rrheight - 10;
    };


    /** maskRectangle constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    *   Parameter: _rwidth: width of the rectangle
    *   Parameter: _rheight: height of the rectangle
    */
    maskRectFrame(int _cx = 0, int _cy = 0, int _rwidth = 1, int _rheight = 1): maskRectangle(_cx, _cy, _rwidth, _rheight)
    {
        innerwidth = rrwidth - 10;
        innerheight = rrheight - 10;
    };

    /** maskRectangle destructor
    */
    ~maskRectFrame()
    {};


    /** Decreases the area of the frame
    */
    virtual bool decreaseArea()
    {
        if ((rrwidth > innerwidth) && (rrheight > innerheight))
        {
            orwidth = rrwidth;
            orheight = rrheight;
            ocx = rcx;
            ocy = rcy;
            rrwidth--;
            rrheight--;
            paint();
            return true;
        }
        else return false;
    };


    /** Increases the area of the frame
    */
    virtual bool increaseArea()
    {
        if (((rcx + rrwidth) < width) && ((rcy + rrheight) < height) && ((rcx - rrwidth) > 0) && ((rcy - rrheight) > 0))
        {
            orwidth = rrwidth;
            orheight = rrheight;
            ocx = rcx;
            ocy = rcy;
            rrwidth++;
            rrheight++;
            paint();
            return true;
        }
        else return false;
    };


    /** Decreases the width one pixel
    */
    virtual bool decreaseWidth()
    {
        if (innerwidth > 1)
        {
            orwidth = rrwidth;
            orheight = rrheight;
            ocx = rcx;
            ocy = rcy;
            innerwidth--;
            rrwidth--;
            paint();
            return true;
        }
        else return false;
    };

    /** Decreases the height one pixel
    */
    virtual bool decreaseHeight()
    {
        if (innerheight > 1)
        {
            orwidth = rrwidth;
            orheight = rrheight;
            ocx = rcx;
            ocy = rcy;
            innerheight--;
            rrheight--;
            paint();
            return true;
        }
        else return false;
    };

    /** Increases the width one pixel
    */
    virtual bool increaseWidth()
    {
        if (((rcx + rrwidth) < width) && ((rcy + rrheight) < height) && ((rcx - rrwidth) > 0) && ((rcy - rrheight) > 0))
        {
            orwidth = rrwidth;
            orheight = rrheight;
            ocx = rcx;
            ocy = rcy;
            rrwidth++;
            innerwidth++;
            paint();
            return true;
        }
        else return false;
    };

    /** Increases the height one pixel
    */
    virtual bool increaseHeight()
    {
        if (((rcx + rrwidth) < width) && ((rcy + rrheight) < height) && ((rcx - rrwidth) > 0) && ((rcy - rrheight) > 0))
        {
            orwidth = rrwidth;
            orheight = rrheight;
            ocx = rcx;
            ocy = rcy;
            rrheight++;
            innerheight++;
            paint();
            return true;
        }
        else return false;
    };


    /** Resizes the figure
    *   Parameter: _width: new width
    *   Parameter: _height: new height
    */
    virtual void resize(int _width, int _height)
    {
        maskFigure::resize(_width, _height);
        rrwidth = (int)((float)(rrwidth * _width) / (float)width);
        rrheight = (int)((float)(rrheight * _width) / (float)width);
        innerwidth = (int)((float)(innerwidth * _width) / (float)width);
        innerheight = (int)((float)(innerheight * _width) / (float)width);
        orwidth = rrwidth;
        orheight = rrheight;
        height = _height;
        width = _width;
        paint();
    };


    /** Paints the figure (virtual method)
    *
    */
    virtual void paint()
    {
        maskFigure::paint();
        QPoint point(ocx - orwidth, ocy - orheight);
        QRect rect(ocx - orwidth, ocy - orheight, orwidth*2, orheight*2);
        drawPixmap(point, (*pixmap), rect);
        drawRect(rcx - rrwidth, rcy - rrheight, rrwidth*2, rrheight*2);
        drawRect(rcx - innerwidth, rcy - innerheight, innerwidth*2, innerheight*2);
    }


    /** this method returns true if the point passed as parameter is inside the frame
    *   Parameter: _cx: X xoordinate of the point to test
    *   Parameter: _cy: X xoordinate of the point to test
    */
    virtual bool isIn(int _cx, int _cy)
    {
        if ((_cx > (rcx - rrwidth)) && (_cy > (rcy - rrheight)) && (_cx < (rcx + rrwidth)) && (_cy < (rcy + rrheight)))
            if (!((_cx >= (rcx - innerwidth)) && (_cy >= (rcy - innerheight)) && (_cx <= (rcx + innerwidth)) && (_cy <= (rcy + innerheight))))
                return true;
        return false;
    };

    /** Prints mask parameters (virtual)*/
    virtual void print(float scalex, float scaley)
    {
        maskRectangle::print(scalex, scaley);
        std::cout << "(half_innerwidth,"
                  "half_innerheight): (" << innerwidth / scalex << ", "
                  << innerheight / scaley << ")"
                  << std::endl;
    }


protected:
    int innerwidth;   // real width of the rectangle
    int innerheight;  // real width of the rectangle
};



/**@name maskSquFrame*/
//@{
/**
 * This class implements square frame to be used as mask.
 * It is used for generating geometric masks.
 */
class maskSquFrame: public maskRectFrame
{
public:

    /** maskSquFrame constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    */
    maskSquFrame(QPaintDevice* _paintDevice, QPixmap* _pixmap, int _width, int _height): maskRectFrame(_paintDevice, _pixmap, _width, _height)
    {};


    /** maskSquFrame constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    *   Parameter: _rwidth: width of the rectangle
    *   Parameter: _rheight: height of the rectangle
    */
    maskSquFrame(int _cx = 0, int _cy = 0, int _rwidth = 1, int _rheight = 1): maskRectFrame(_cx, _cy, _rwidth, _rheight)
    {};

    /** maskSquFrame destructor
    */
    ~maskSquFrame()
    {};

    /** Decreases the width one pixel
    */
    virtual bool decreaseWidth()
    {
        maskRectFrame::decreaseWidth();
        return maskRectFrame::decreaseHeight();
    };

    /** Decreases the height one pixel
    */
    virtual bool decreaseHeight()
    {
        maskRectFrame::decreaseWidth();
        return maskRectFrame::decreaseHeight();
    };

    /** Increases the width one pixel
    */
    virtual bool increaseWidth()
    {
        maskRectFrame::increaseWidth();
        return maskRectFrame::increaseHeight();
    };

    /** Increases the height one pixel
    */
    virtual bool increaseHeight()
    {
        maskRectFrame::increaseWidth();
        return maskRectFrame::increaseHeight();
    };

};

/**@name maskEllipse*/
//@{
/**
 * This class implements circle figure to b used as mask.
 * It is used for generating geometric masks.
 */
class maskPolygon: public maskFigure
{
public:

    /** maskPolygon constructor
    *   Parameter: _cx: X-coordinate
    *   Parameter: _cy: Y-coordinate
    */
    maskPolygon(QPaintDevice* _paintDevice, QPixmap* _pixmap, int _width, int _height): maskFigure(_paintDevice, _pixmap, _width, _height)
    {
        close_polygon = false;
    };

    /** Moves the figure
    *   Parameter: _cx: new X-coordinate
    *   Parameter: _cy: new Y-coordinate
    */
    virtual bool move(int _cx, int _cy)
    {
        if (list_of_points.size() == 0) return true;
        std::cout << "Is this ever called? Contact coss@cnb.csic.es\n";
    };

    /** Moves the figure to the left one pixel
    */
    virtual bool moveLeft()
    {
        if (list_of_points.size() == 0) return true;
        if (XX(top_left_corner) > 0)
        {
            old_top_left_corner = top_left_corner;
            old_bottom_right_corner = bottom_right_corner;
            for (int i = 0; i < list_of_points.size(); i++)
                --XX(list_of_points[i]);
            --XX(top_left_corner);
            --XX(bottom_right_corner);
            paint();
            return true;
        }
        else return false;
    };

    /** Moves the figure up one pixel
    */
    virtual bool moveUp()
    {
        if (list_of_points.size() == 0) return true;
        if (YY(top_left_corner) > 0)
        {
            old_top_left_corner = top_left_corner;
            old_bottom_right_corner = bottom_right_corner;
            for (int i = 0; i < list_of_points.size(); i++)
                --YY(list_of_points[i]);
            --YY(top_left_corner);
            --YY(bottom_right_corner);
            paint();
            return true;
        }
        else return false;
    };

    /** Moves the figure to the right one pixel
    */
    virtual bool moveRight()
    {
        if (list_of_points.size() == 0) return true;
        if (XX(bottom_right_corner) < width)
        {
            old_top_left_corner = top_left_corner;
            old_bottom_right_corner = bottom_right_corner;
            for (int i = 0; i < list_of_points.size(); i++)
                ++XX(list_of_points[i]);
            ++XX(top_left_corner);
            ++XX(bottom_right_corner);
            paint();
            return true;
        }
        else return false;
    };

    /** Moves the figure down one pixel
    */
    virtual bool moveDown()
    {
        if (list_of_points.size() == 0) return true;
        if (YY(bottom_right_corner) < height)
        {
            old_top_left_corner = top_left_corner;
            old_bottom_right_corner = bottom_right_corner;
            for (int i = 0; i < list_of_points.size(); i++)
                ++YY(list_of_points[i]);
            ++YY(top_left_corner);
            ++YY(bottom_right_corner);
            paint();
            return true;
        }
        else return false;
    };

    /** Decreases the width one pixel
    */
    virtual bool decreaseWidth()
    {
        std::cout << "Resizing is not available for polygons\n";
        return true;
    };

    /** Decreases the height one pixel
    */
    virtual bool decreaseHeight()
    {
        std::cout << "Resizing is not available for polygons\n";
        return true;
    };

    /** Increases the width one pixel
    */
    virtual bool increaseWidth()
    {
        std::cout << "Resizing is not available for polygons\n";
        return true;
    };

    /** Increases the height one pixel
    */
    virtual bool increaseHeight()
    {
        std::cout << "Resizing is not available for polygons\n";
        return true;
    };

    /** Add point to the polygon */
    void add_point(int x, int y)
    {
        if (close_polygon) return;
        old_top_left_corner = top_left_corner;
        old_bottom_right_corner = bottom_right_corner;
        list_of_points.push_back(vectorR2(x, y));
        if (list_of_points.size() > 1)
        {
            XX(top_left_corner) = XMIPP_MIN(XX(top_left_corner), x);
            YY(top_left_corner) = XMIPP_MIN(YY(top_left_corner), y);
            XX(bottom_right_corner) = XMIPP_MAX(XX(bottom_right_corner), x);
            YY(bottom_right_corner) = XMIPP_MAX(YY(bottom_right_corner), y);
        }
        else
        {
            top_left_corner = vectorR2(x, y);
            bottom_right_corner = vectorR2(x, y);
            old_top_left_corner = top_left_corner;
            old_bottom_right_corner = bottom_right_corner;
        }
        paint();
    }

    /** Close polygon */
    void closePolygon()
    {
        close_polygon = true;
        list_of_points.push_back(list_of_points[0]);
        paint();
    }

    /** Move point of the polygon */
    void move_point(int x, int y)
    {
        if (close_polygon || list_of_points.size() <= 1) return;
        old_top_left_corner = top_left_corner;
        old_bottom_right_corner = bottom_right_corner;
        int imax = list_of_points.size();
        XX(list_of_points[imax-1]) = x;
        YY(list_of_points[imax-1]) = y;
        XX(top_left_corner) = XMIPP_MIN(XX(top_left_corner), x);
        YY(top_left_corner) = XMIPP_MIN(YY(top_left_corner), y);
        XX(bottom_right_corner) = XMIPP_MAX(XX(bottom_right_corner), x);
        YY(bottom_right_corner) = XMIPP_MAX(YY(bottom_right_corner), y);
        paint();
    }

    /** Clear all points */
    void clear_all_points()
    {
        int imax = list_of_points.size();
        if (imax == 0 || close_polygon) return;
        if (imax > 1)
        {
            QPoint point((int)XX(old_top_left_corner), (int)YY(old_top_left_corner));
            QRect rect((int)XX(old_top_left_corner) - 2, (int)YY(old_top_left_corner) - 2,
                       (int)(XX(old_bottom_right_corner) - XX(old_top_left_corner)) + 4,
                       (int)(YY(old_bottom_right_corner) - YY(old_top_left_corner)) + 4);
            drawPixmap(point, (*pixmap), rect);
        }
        list_of_points.clear();
    }

    /** Resizes the figure
    *   Parameter: _width: new width
    *   Parameter: _height: new height
    */
    virtual void resize(int _width, int _height)
    {
        maskFigure::resize(_width, _height);
        double factor = (double)_width / (double)width;
        for (int i = 0; i < list_of_points.size(); i++)
            list_of_points[i]*factor;
        old_top_left_corner = top_left_corner;
        old_bottom_right_corner = bottom_right_corner;
        height = _height;
        width = _width;
        paint();
    };

    /** Paints the figure (virtual method)
    *
    */
    virtual void paint()
    {
        if (list_of_points.size() == 0) return;
        maskFigure::paint();
        int imax = list_of_points.size();
        if (imax > 1)
        {
            QPoint point((int)XX(old_top_left_corner), (int)YY(old_top_left_corner));
            QRect rect((int)XX(old_top_left_corner), (int)YY(old_top_left_corner),
                       (int)(XX(old_bottom_right_corner) - XX(old_top_left_corner)),
                       (int)(YY(old_bottom_right_corner) - YY(old_top_left_corner)));
            drawPixmap(point, (*pixmap), rect);
        }
        drawPoint((int)XX(list_of_points[0]), (int)YY(list_of_points[0]));
        for (int i = 1; i < imax; i++)
            drawLine((int)XX(list_of_points[i-1]), (int)YY(list_of_points[i-1]),
                     (int)XX(list_of_points[i]), (int)YY(list_of_points[i]));
        old_top_left_corner = top_left_corner;
        old_bottom_right_corner = bottom_right_corner;
    }

    /** this method returns true if the point passed as parameter is inside the circle
    *   Parameter: _cx: X xoordinate of the point to test
    *   Parameter: _cy: X xoordinate of the point to test
    */
    virtual bool isIn(int _cx, int _cy)
    {
        if (!close_polygon) return false;
        else return point_inside_polygon(list_of_points, vectorR2(_cx, _cy));
    };

    /** Prints mask parameters (virtual)*/
    virtual void print(float scalex, float scaley)
    {
        std::cout << "List of points\n";
        for (int i = 0; i < list_of_points.size(); i++)
            std::cout << "Point " << i << " (X,Y)=("
                      << XX(list_of_points[i]) / scalex << ","
                      << YY(list_of_points[i]) / scaley << ")\n";
    }

protected:
    std::vector< Matrix1D<double> > list_of_points;
    Matrix1D<double>                top_left_corner;
    Matrix1D<double>                bottom_right_corner;
    Matrix1D<double>                old_top_left_corner;
    Matrix1D<double>                old_bottom_right_corner;
    bool                            close_polygon;
};

/**@name maskImg*/
//@{
/**
 * This class implements the GUI for selecting a mask.
 */


class QLabel;
class QMenuBar;
#ifdef QT3_SUPPORT
class Q3PopupMenu;
#else
class QPopupMenu;
#endif

class maskImg : public QWidget
{
    Q_OBJECT
    QPixmap pmScaled;  // the scaled pixmap
public:
    maskImg(QWidget *parent = 0, const char *name = 0, int wFlags = 0);
    maskImg(QWidget *parent = 0, ImageT<double> *_image = 0, maskType _typeOfMask = CIRCLE, const char *name = 0, int wFlags = 0);
    maskImg(QWidget *parent = 0, QImage *_image = 0, maskType _typeOfMask = CIRCLE, const char *name = 0, int wFlags = 0);
    ~maskImg();

    /** Apply transformations that are stored in the headers of the images? */
    bool        apply_geo;
    std::string saveasname;
    bool loadImage(const char *fileName);
    Image xmippImage;    // Xmipp Image

protected:
    void paintEvent(QPaintEvent *);
    void resizeEvent(QResizeEvent *);
    void mousePressEvent(QMouseEvent *);
    void mouseDoubleClickEvent(QMouseEvent *);
    void mouseMoveEvent(QMouseEvent *);
    void  keyPressEvent(QKeyEvent*);

private:
    void scale();
    int  alloc_context;
    bool convertEvent(QMouseEvent* e, int& x, int& y);
    const char* filename;
    QImage image;   // the loaded image
    QPixmap pm;   // the converted pixmap
    int  xmippFlag;  // stands for:
    // -1: No Xmipp image
    //  0: Xmipp image
#ifdef QT3_SUPPORT
    Q3PopupMenu *menubar;
    Q3PopupMenu *file;
    Q3PopupMenu *options;
    Q3PopupMenu *saveimage;
    Q3PopupMenu *savepixmap;
#else
    QPopupMenu *menubar;
    QPopupMenu *file;
    QPopupMenu *options;
    QPopupMenu *saveimage;
    QPopupMenu *savepixmap;
#endif
    QWidget    *helpmsg;
    QLabel     *status;
    int         circle, rect, ellip, squ, ring, ellipring, rectframe,
      squframe, si, polygon; // Menu item ids
    void Init();
    bool  xmipp2Qt(ImageT<double> &_image);
    bool  Qt2xmipp(QImage _image);
    bool  showImage();
    void updateStatus();
    bool  reconvertImage();
    int  pickx, picky;
    int  clickx, clicky;
    static  maskImg* other;
    Image       mask;
    maskType    typeOfMask;
    int         cx, cy, w, h, KEYevent;
    maskFigure* theMaskFigure;
private slots:
    void newWindow();
    void openFile();
    void Close();
    void saveImage(int);
    void giveHelp();
    void doOption(int);
    void        about();
    void        aboutXmipp();

};

#endif // MASKIMG_H
