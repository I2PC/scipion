/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
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

#include "maskimg.h"

#include <qmenubar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qlabel.h>
#include <qpainter.h>
#include <qkeycode.h>
#include <qapplication.h>

/****************************************************/

/*
  In the constructor, we just pass the standard parameters on to
  QWidget.

  The menu uses a single slot to simplify the process of adding
  more items to the options menu.
*/

void maskImg::Init()
{
    pickx = -1;
    clickx = -1;
    alloc_context = 0;

    theMaskFigure = NULL;
    menubar = new QPopupMenu();

    file = new QPopupMenu();
    menubar->insertItem("&File", file);
    file->insertItem("New window", this,  SLOT(newWindow()));
    file->insertItem("Open...", this,  SLOT(openFile()));
    file->insertItem("Save mask...", this,  SLOT(saveImage(int)));

    options =  new QPopupMenu();
    menubar->insertItem("&Mask types", options);
    circle = options->insertItem("Circle");
    ellip = options->insertItem("Ellipse");
    rect = options->insertItem("Rectangle");
    squ = options->insertItem("Square");
    ring = options->insertItem("Circular Crown");
    ellipring = options->insertItem("Elliptical Crown");
    rectframe = options->insertItem("Rectangular Frame");
    squframe = options->insertItem("Square Frame");
    polygon = options->insertItem("Polygon");
    options->setCheckable(TRUE);

    menubar->insertSeparator();

    QPopupMenu* help = new QPopupMenu();
    menubar->insertItem("&Help", help);
    help->insertItem("&About", this, SLOT(about()));
    help->insertItem("About &Xmipp", this, SLOT(aboutXmipp()));
    help->insertSeparator();
    help->insertItem("Help!", this, SLOT(giveHelp()));

    menubar->insertSeparator();
    menubar->insertItem("Quit", this,  SLOT(Close()));

    connect(options, SIGNAL(activated(int)), this, SLOT(doOption(int)));

    status = new QLabel(this);
    status->setFrameStyle(QFrame::WinPanel | QFrame::Sunken);
    status->setFixedHeight(fontMetrics().height() + 4);

    setMouseTracking(TRUE);
    options->setItemChecked(circle, true);
}


/****************************************************/
maskImg::maskImg(QWidget *parent, const char *name, int wFlags)
        : QWidget(parent, name, wFlags),
        filename(0),
        helpmsg(0)
{
    Init();
    typeOfMask = CIRCLE;
}

/****************************************************/

maskImg::maskImg(QWidget *parent, QImage *_image, maskType _typeOfMask, const char *name, int wFlags)
        : QWidget(parent, name, wFlags),
        filename(0),
        helpmsg(0)
{
    Init();
    filename = name;
    if (Qt2xmipp((QImage) *_image)) showImage();
    typeOfMask = _typeOfMask;
    if (typeOfMask == CIRCLE)
        theMaskFigure = new maskCircle(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == ELLIPSE)
        theMaskFigure = new maskEllipse(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == RECTANGLE)
        theMaskFigure = new maskRectangle(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == SQUARE)
        theMaskFigure = new maskSquare(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == RING)
        theMaskFigure = new maskCircleRing(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == ELLIPRING)
        theMaskFigure = new maskEllipRing(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == RECTFRAME)
        theMaskFigure = new maskRectFrame(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == SQUFRAME)
        theMaskFigure = new maskSquFrame(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == POLYGON)
        theMaskFigure = new maskPolygon(this, &pmScaled, width(), height() - status->height());
    theMaskFigure->paint();
}

/****************************************************/

maskImg::maskImg(QWidget *parent, ImageT<double> *_image, maskType _typeOfMask, const char *name, int wFlags)
        : QWidget(parent, name, wFlags),
        filename(0),
        helpmsg(0)
{
    Init();
    filename = name;
    if (xmipp2Qt(*_image))
        showImage();
    typeOfMask = _typeOfMask;
    if (typeOfMask == CIRCLE)
        theMaskFigure = new maskCircle(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == ELLIPSE)
        theMaskFigure = new maskEllipse(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == RECTANGLE)
        theMaskFigure = new maskRectangle(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == SQUARE)
        theMaskFigure = new maskSquare(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == RING)
        theMaskFigure = new maskCircleRing(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == ELLIPRING)
        theMaskFigure = new maskEllipRing(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == RECTFRAME)
        theMaskFigure = new maskRectFrame(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == SQUFRAME)
        theMaskFigure = new maskSquFrame(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == POLYGON)
        theMaskFigure = new maskPolygon(this, &pmScaled, width(), height() - status->height());
    theMaskFigure->paint();
}

/****************************************************/

maskImg::~maskImg()
{
    if (alloc_context)
        QColor::destroyAllocContext(alloc_context);
    if (other == this)
        other = 0;
    delete theMaskFigure;
}

/****************************************************/

/*
  This function modifies the conversion_flags when an options menu item
  is selected, then ensures all menu items are up to date, and reconverts
  the image if possibly necessary.
*/
void maskImg::doOption(int item)
{

    if (options->isItemChecked(item)) return;     // They are all radio buttons

    if (theMaskFigure != NULL)
        delete theMaskFigure;
    if (item == circle)
    {
        options->setItemChecked(circle, true);
        options->setItemChecked(rect, false);
        options->setItemChecked(ellip, false);
        options->setItemChecked(squ, false);
        options->setItemChecked(ring, false);
        options->setItemChecked(ellipring, false);
        options->setItemChecked(rectframe, false);
        options->setItemChecked(squframe, false);
        options->setItemChecked(polygon, false);
        typeOfMask = CIRCLE;
        theMaskFigure = new maskCircle(this, &pmScaled, width(), height() - status->height());
    }
    else if (item == rect)
    {
        options->setItemChecked(circle, false);
        options->setItemChecked(rect, true);
        options->setItemChecked(ellip, false);
        options->setItemChecked(squ, false);
        options->setItemChecked(ring, false);
        options->setItemChecked(ellipring, false);
        options->setItemChecked(rectframe, false);
        options->setItemChecked(squframe, false);
        options->setItemChecked(polygon, false);
        typeOfMask = RECTANGLE;
        theMaskFigure = new maskRectangle(this, &pmScaled, width(), height() - status->height());
    }
    else if (item == ellip)
    {
        options->setItemChecked(circle, false);
        options->setItemChecked(rect, false);
        options->setItemChecked(ellip, true);
        options->setItemChecked(squ, false);
        options->setItemChecked(ring, false);
        options->setItemChecked(ellipring, false);
        options->setItemChecked(rectframe, false);
        options->setItemChecked(squframe, false);
        options->setItemChecked(polygon, false);
        typeOfMask = ELLIPSE;
        theMaskFigure = new maskEllipse(this, &pmScaled, width(), height() - status->height());
    }
    else if (item == squ)
    {
        options->setItemChecked(circle, false);
        options->setItemChecked(rect, false);
        options->setItemChecked(ellip, false);
        options->setItemChecked(squ, true);
        options->setItemChecked(ring, false);
        options->setItemChecked(ellipring, false);
        options->setItemChecked(rectframe, false);
        options->setItemChecked(squframe, false);
        options->setItemChecked(polygon, false);
        typeOfMask = SQUARE;
        theMaskFigure = new maskSquare(this, &pmScaled, width(), height() - status->height());
    }
    else if (item == ring)
    {
        options->setItemChecked(circle, false);
        options->setItemChecked(rect, false);
        options->setItemChecked(ellip, false);
        options->setItemChecked(squ, false);
        options->setItemChecked(ring, true);
        options->setItemChecked(ellipring, false);
        options->setItemChecked(rectframe, false);
        options->setItemChecked(squframe, false);
        options->setItemChecked(polygon, false);
        typeOfMask = RING;
        theMaskFigure = new maskCircleRing(this, &pmScaled, width(), height() - status->height());
    }
    else if (item == ellipring)
    {
        options->setItemChecked(circle, false);
        options->setItemChecked(rect, false);
        options->setItemChecked(ellip, false);
        options->setItemChecked(squ, false);
        options->setItemChecked(ring, false);
        options->setItemChecked(ellipring, true);
        options->setItemChecked(rectframe, false);
        options->setItemChecked(squframe, false);
        options->setItemChecked(polygon, false);
        typeOfMask = ELLIPRING;
        theMaskFigure = new maskEllipRing(this, &pmScaled, width(), height() - status->height());
    }
    else if (item == rectframe)
    {
        options->setItemChecked(circle, false);
        options->setItemChecked(rect, false);
        options->setItemChecked(ellip, false);
        options->setItemChecked(squ, false);
        options->setItemChecked(ring, false);
        options->setItemChecked(ellipring, false);
        options->setItemChecked(rectframe, true);
        options->setItemChecked(squframe, false);
        options->setItemChecked(polygon, false);
        typeOfMask = RECTFRAME;
        theMaskFigure = new maskRectFrame(this, &pmScaled, width(), height() - status->height());
    }
    else if (item == squframe)
    {
        options->setItemChecked(circle, false);
        options->setItemChecked(rect, false);
        options->setItemChecked(ellip, false);
        options->setItemChecked(squ, false);
        options->setItemChecked(ring, false);
        options->setItemChecked(ellipring, false);
        options->setItemChecked(rectframe, false);
        options->setItemChecked(squframe, true);
        options->setItemChecked(polygon, false);
        typeOfMask = SQUFRAME;
        theMaskFigure = new maskSquFrame(this, &pmScaled, width(), height() - status->height());
    }
    else if (item == polygon)
    {
        options->setItemChecked(circle, false);
        options->setItemChecked(rect, false);
        options->setItemChecked(ellip, false);
        options->setItemChecked(squ, false);
        options->setItemChecked(ring, false);
        options->setItemChecked(ellipring, false);
        options->setItemChecked(rectframe, false);
        options->setItemChecked(squframe, false);
        options->setItemChecked(polygon, true);
        typeOfMask = POLYGON;
        theMaskFigure = new maskPolygon(this, &pmScaled, width(), height() - status->height());
    }
    theMaskFigure->paint();

    // And reconvert...
    repaint(image.hasAlphaBuffer()); // show image in widget
}

/****************************************************/

void maskImg::updateStatus()
{
    if (pm.size() == QSize(0, 0))
    {
        if (filename)
            status->setText("Could not load image");
        else
            status->setText("No image - select Open from File menu.");
    }
    else
    {
        QString message, moremsg;

        if (image.valid(pickx, picky))
        {
            int y_log, x_log;
            xmippImage().toLogical(picky, pickx, y_log, x_log);
            moremsg.sprintf("(%d,%d)= %.3f ",
                            x_log, y_log,
                            xmippImage(y_log, x_log));
            message += moremsg;
        }
        moremsg.sprintf("%dx%d", image.width(), image.height());
        message += moremsg;
        if (pm.size() != pmScaled.size())
        {
            moremsg.sprintf(" [%dx%d]", pmScaled.width(),
                            pmScaled.height());
            message += moremsg;
        }
        status->setText(message);
    }
}

/****************************************************/

/*
  This function exits.
*/
void maskImg::Close()
{

    if (saveasname != "")
    {
        saveImage(0);
    }
    close();

}



/****************************************************/

/*
  This function saves the image.
*/
void maskImg::saveImage(int item)
{
    //const char* fmt = saveimage->text(item);
//    QString savefilename = QFileDialog::getSaveFileName(0, 0, 0, filename);
    QString savefilename;

    if (saveasname == "")
    {
        savefilename = QFileDialog::getSaveFileName(QString::null, "*.msk", this);
    }
    else
    {
        const char* c_string = saveasname.c_str();
        savefilename = QString(c_string);
    }

    if (!savefilename.isEmpty())
    {
        QFileInfo fi(savefilename);
        if (fi.exists())
        {
            if (QMessageBox::information(this, "ImageViewer application",
                                         "The file already exist. Overwrite?",
                                         "Yes",
                                         "No") != 0)
            {
                QMessageBox::about(this, "Warning!", "Saving aborted\n");
                return;
            }
        }
        /*        if (strcmp(fmt, "Spider") != 0) {
            if ( !image.save( savefilename, fmt ) )
             QMessageBox::warning( this, "Save failed", "Error saving file" );
         } else {*/
        try
        {
            // Creates a Xmipp temporal Image

            ImageXmipp tmpImage(xmippImage().rowNumber(), xmippImage().colNumber());

            // Writes pixels.
            for (int y = 0; y < tmpImage().rowNumber(); y++)
                for (int x = 0; x < tmpImage().colNumber(); x++)
                {
                    int phy_x, phy_y;
                    phy_x = (int)((float)(x * pmScaled.width()) / (float)xmippImage().colNumber());
                    phy_y = (int)((float)(y * pmScaled.height()) / (float)xmippImage().rowNumber());
                    if (theMaskFigure->isIn(phy_x, phy_y))
                        tmpImage(y, x) = (float) 1.0;
                    else
                        tmpImage(y, x) = (float) 0.0;
                }

            // Saves Xmipp Image
            tmpImage.rename((string)((const  char *)savefilename));
            tmpImage.write();
            // print parameters
            theMaskFigure->print(pmScaled.width() / (float)xmippImage().colNumber(),
                                 pmScaled.height() / (float)xmippImage().rowNumber());
        }
        catch (Xmipp_error)
        {
            char *helptext = "Invalid image type";
            helpmsg = new QMessageBox("Error", helptext,
                                      QMessageBox::Information, QMessageBox::Ok, 0, 0, 0, 0, FALSE);
            helpmsg->show();
            helpmsg->raise();
        }
    }
}

/****************************************************/

void maskImg::newWindow()
{
    maskImg* that = new maskImg(0, "new window", WDestructiveClose);
    that->show();
}

/*
  This function is the slot for processing the Open menu item.
*/
void maskImg::openFile()
{
    QString newfilename = QFileDialog::getOpenFileName();
    if (!newfilename.isEmpty())
    {
        loadImage(newfilename) ;
        repaint();   // show image in widget
    }
}

/*****************************************/

bool maskImg::showImage()
{
    bool ok = FALSE;
    QApplication::setOverrideCursor(waitCursor);   // this might take time
    pickx = -1;
    clickx = -1;
    ok = reconvertImage();
    if (ok)
    {
        setCaption(filename);     // set window caption
        int w = pm.width();
        int h = pm.height();
        const int reasonable_width = 128;
        if (w < reasonable_width)
        {
            // Integer scale up to something reasonable
            int multiply = (reasonable_width + w - 1) / w;
            w *= multiply;
            h *= multiply;
        }
        h += status->height();
        resize(w, h);          // we resize to fit image
    }
    else
    {
        pm.resize(0, 0);       // couldn't load image
        update();
    }
    QApplication::restoreOverrideCursor();  // restore original cursor
    updateStatus();
    return ok;
}

/*****************************************/

bool maskImg::xmipp2Qt(ImageT<double> &_image)
{
    bool ok = FALSE;

    try
    {
        xmippImage = _image;
        // Creates a Qt Image to hold Xmipp Image
        QImage tmpImage(_image().colNumber(), _image().rowNumber(), 8, 256);
        _image().range_adjust(0, 255);

        // Sets Graylevel Palette.
        for (int i = 0; i < 256; i++)
        {
            QColor c;
            c.setRgb(i, i, i);
            tmpImage.setColor(i, c.rgb());
        }
        // Reads pixels.
        for (int y = 0; y < _image().rowNumber(); y++)
            for (int x = 0; x < _image().colNumber(); x++)
                tmpImage.setPixel(x, y, ((uint) _image(y, x)));

        image = tmpImage;
        xmippFlag = 0;   // Sets flag = Xmipp image.
        ok = TRUE;

    }
    catch (Xmipp_error)
    {
        ok = FALSE;
        char *helptext = "Error converting xmipp to Qt";
        helpmsg = new QMessageBox("Error", helptext,
                                  QMessageBox::Information, QMessageBox::Ok, 0, 0, 0, 0, FALSE);
        helpmsg->show();
        helpmsg->raise();
    }

    return ok;
}

/*****************************************/

bool maskImg::Qt2xmipp(QImage _image)
{
    bool ok = FALSE;
    // try to read image from standard format.

    try
    {
        image = _image;
        image.setNumColors(256);

        // Creates a Xmipp Image to hold Qt Image
        Image tmpImage(image.height(), image.width());

        // Reads pixels.
        for (int y = 0; y < image.width(); y++)
            for (int x = 0; x < image.height(); x++)
                tmpImage(x, y) = (double) image.pixelIndex(y, x);


        xmippImage = tmpImage;
        xmippFlag = 0;   // Sets flag = Xmipp image.

        ok = TRUE;

    }
    catch (Xmipp_error)
    {
        ok = FALSE;
        char *helptext = "Error converting Qt image to Xmipp image";
        helpmsg = new QMessageBox("Error", helptext,
                                  QMessageBox::Information, QMessageBox::Ok, 0, 0, 0, 0, FALSE);
        helpmsg->show();
        helpmsg->raise();
    }

    return ok;
}


/*****************************************/

/*
  This function loads an image from a file and resizes the widget to
  exactly fit the image size. If the file was not found or the image
  format was unknown it will resize the widget to fit the errorText
  message (see above) displayed in the current font.

  Returns TRUE if the image was successfully loaded.
*/

bool maskImg::loadImage(const char *fileName)
{
    filename = fileName;
    bool ok = FALSE;
    if (filename)
    {

        // try to read image from standard format.

        if (image.load(filename, 0)) ok = Qt2xmipp(image);

        if (!ok)
        {
            try
            {
                // reads Xmipp Image
                Image *tmpImage = Image::LoadImage(fileName, apply_geo);
                if (!tmpImage)
                    REPORT_ERROR(1501, (string)"Error opening file " + fileName);
                ok = TRUE;
                ok = xmipp2Qt(*tmpImage);
                delete tmpImage;
            }
            catch (Xmipp_error)
            {
                ok = FALSE;
                char *helptext = "Invalid image type";
                helpmsg = new QMessageBox("Error", helptext,
                                          QMessageBox::Information, QMessageBox::Ok, 0, 0, 0, 0, FALSE);
                helpmsg->show();
                helpmsg->raise();
            }
        }
    }

    if (ok) ok = showImage();
    if (theMaskFigure != NULL)
        delete theMaskFigure;
    if (typeOfMask == CIRCLE)
        theMaskFigure = new maskCircle(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == ELLIPSE)
        theMaskFigure = new maskEllipse(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == RECTANGLE)
        theMaskFigure = new maskRectangle(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == SQUARE)
        theMaskFigure = new maskSquare(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == RING)
        theMaskFigure = new maskCircleRing(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == ELLIPRING)
        theMaskFigure = new maskEllipRing(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == RECTFRAME)
        theMaskFigure = new maskRectFrame(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == SQUFRAME)
        theMaskFigure = new maskSquFrame(this, &pmScaled, width(), height() - status->height());
    else if (typeOfMask == POLYGON)
        theMaskFigure = new maskPolygon(this, &pmScaled, width(), height() - status->height());
    theMaskFigure->paint();
    return ok;
}

/*****************************************/

bool maskImg::reconvertImage()
{
    bool success = FALSE;

    if (image.isNull()) return FALSE;
    if (alloc_context)
    {
        QColor::destroyAllocContext(alloc_context);
        alloc_context = 0;
    }
    QApplication::setOverrideCursor(waitCursor);   // this might take time
    if (pm.convertFromImage(image))
    {
        pmScaled = QPixmap();
        scale();
        resize(width(), height());
        success = TRUE;    // load successful
    }
    else
    {
        pm.resize(0, 0);   // couldn't load image
    }
    updateStatus();
    QApplication::restoreOverrideCursor(); // restore original cursor
    return success;    // TRUE if loaded OK
}


/****************************************************/

/*
  This functions scales the pixmap in the member variable "pm" to fit the
  widget size and  puts the resulting pixmap in the member variable "pmScaled".
*/

void maskImg::scale()
{
    int h = height() - status->height();

    if (image.isNull()) return;

    QApplication::setOverrideCursor(waitCursor);   // this might take time
    if (width() == pm.width() && h == pm.height())
    {      // no need to scale if widget
        pmScaled = pm;    // size equals pixmap size
    }
    else
    {
        QWMatrix m;    // transformation matrix
        m.scale(((double)width()) / pm.width(),// define scale factors
                ((double)h) / pm.height());
        pmScaled = pm.xForm(m);    // create scaled pixmap
    }
    if (theMaskFigure != NULL)
        theMaskFigure->resize(width(), height() - status->height());
    QApplication::restoreOverrideCursor(); // restore original cursor
}


/****************************************************/

/*
  The resize event handler, if a valid pixmap was loaded it will call
  scale() to fit the pixmap to the new widget size.
*/

void maskImg::resizeEvent(QResizeEvent *)
{
    status->setGeometry(0, height() - status->height(),
                        width(), status->height());

    if (pm.size() == QSize(0, 0))      // we couldn't load the image
        return;

    int h = height() - status->height();
    if (width() != pmScaled.width() || h != pmScaled.height())
    {      // if new size,
        scale();    // scale pmScaled to window
        updateStatus();
    }
}


/****************************************************/

/*
  Handles a change in the mouse
*/

bool maskImg::convertEvent(QMouseEvent* e, int& x, int& y)
{
    if (pm.size() != QSize(0, 0))
    {
        int h = height() - status->height();
        int nx = e->x() * image.width() / width();
        int ny = (e->y()) * image.height() / h;
        if (nx != x || ny != y)
        {
            x = nx;
            y = ny;
            return true;
        }
    }
    return false;
}

/****************************************************/

/*
  Mouse press events.
*/

void maskImg::mousePressEvent(QMouseEvent *e)
{
    QPoint clickedPos = e->pos();  // extract pointer position
    if (e->button() == RightButton)
        menubar->exec(clickedPos);
    else if (e->button() == LeftButton && typeOfMask == POLYGON)
    {
        maskPolygon *polygon = dynamic_cast<maskPolygon *>(theMaskFigure);
        polygon->add_point(clickedPos.x(), clickedPos.y());
    }
}


/****************************************************/
void maskImg::mouseDoubleClickEvent(QMouseEvent *e)
{
    QPoint clickedPos = e->pos();  // extract pointer position
    if (typeOfMask == POLYGON)
    {
        maskPolygon *polygon = dynamic_cast<maskPolygon *>(theMaskFigure);
        polygon->add_point(clickedPos.x(), clickedPos.y());
        polygon->closePolygon();
    }
}

/****************************************************/

/*
  Record the pixel position of interest.
*/
void maskImg::mouseMoveEvent(QMouseEvent *e)
{
    if (convertEvent(e, pickx, picky))
        updateStatus();
    if (typeOfMask == POLYGON)
    {
        maskPolygon *polygon = dynamic_cast<maskPolygon *>(theMaskFigure);
        polygon->move_point(pickx, picky);
    }
}



/****************************************************/

/*
  Handles key press events.
*/

void maskImg::keyPressEvent(QKeyEvent* e)
{
    switch (e->key())
    {   // Look at the key code
    case Key_Escape:
        if (typeOfMask == POLYGON)
        {
            maskPolygon *polygon = dynamic_cast<maskPolygon *>(theMaskFigure);
            polygon->clear_all_points();
        }
        break;
    case Key_R:
        if (e->state() == ControlButton)
        { // If 'Ctrol R' key,
            xmippImage().moveOriginTo(-xmippImage().startingY(), -xmippImage().startingX());// sets origin at the upper left corner
        }
        break;
    case Key_O:    // Xmipp origin
        if (e->state() == ControlButton)
        { // If 'Ctrol N' key,
            xmippImage().setXmippOrigin(); // sets origin at the center of the iamge.
        }
        break;
    case Key_Q:    // Quit program
        if (e->state() == ControlButton)
        { // If 'Ctrol Q' key,
            close();                      // Exit
        }
        break;
    case Key_N:    // Natural size (original size)
        if (e->state() == ControlButton)
        { // If 'Ctrol N' key,
            resize(xmippImage().colNumber(), xmippImage().rowNumber() + status->height());
        }
        break;
    case Key_M:
    case Key_Minus:    // Half size
        if (e->state() == ControlButton)
        { // If 'Ctrol-' key
            if (width() > pm.width() / 3)
                resize(width() / 2, height() / 2 + status->height() / 2 + 1);
        }
        else if (e->state() == AltButton)
        { // If 'CtrolShift+' key,
            theMaskFigure->decreaseArea();
        }
        break;
    case Key_P:
    case Key_Plus:           // Double size
        if (e->state() == ControlButton)
        { // If 'Ctrol+' key,
            resize(width()*2, height()*2 - status->height());
        }
        else if (e->state() == AltButton)
        { // If 'CtrolShift+' key,
            theMaskFigure->increaseArea();
        }
        break;
    case Key_A:        // Aspect ratio
        if (e->state() == ControlButton)
        { // If 'Ctrol+' key,
            double ratio = (double) xmippImage().colNumber() / (double) xmippImage().rowNumber();
            resize((int)width(), (int)(width() / ratio + status->height()));
        }
        break;
    case Key_Left:    // If 'left arrow'-key,
        if (e->state() == ControlButton)  // If 'Ctrol+' key,
            theMaskFigure->decreaseWidth();
        else
            theMaskFigure->moveLeft();
        break;
    case Key_Right:    // Correspondingly...
        if (e->state() == ControlButton) // If 'Ctrol+' key,
            theMaskFigure->increaseWidth();
        else
            theMaskFigure->moveRight();
        break;
    case Key_Up:
        if (e->state() == ControlButton) // If 'Ctrol' key,
            theMaskFigure->increaseHeight();
        else
            theMaskFigure->moveUp();
        break;
    case Key_Down:
        if (e->state() == ControlButton) // If 'Ctrol' key,
            theMaskFigure->decreaseHeight();
        else
            theMaskFigure->moveDown();
        break;
    default:    // If not an interesting key,
        e->ignore();   // we don't accept the event
        return;
    }
}


/****************************************************/

/*
  Draws the portion of the scaled pixmap that needs to be updated or prints
  an error message if no legal pixmap has been loaded.
*/

void maskImg::paintEvent(QPaintEvent *e)
{
    if (pm.size() != QSize(0, 0))
    {  // is an image loaded?
        QPainter painter(this);
        painter.setClipRect(e->rect());
        painter.drawPixmap(0, 0, pmScaled);
        if (theMaskFigure != NULL)
            theMaskFigure->paint();
    }
}


/****************************************************/

/**
  Explain anything that might be confusing.
*/
void maskImg::giveHelp()
{
    if (!helpmsg)
    {
        QString helptext = "Usage: xmipp_xmask <-img> <-sel> [-sd]\n\n ";
        QStrList support = QImage::outputFormats();
        helptext += "\n\nSupported input formats:\n";
        int lastnl = helptext.length();

        support.clear();
        support.insert(0, "Spider");
        support.insert(1, "bmp");

        const char* f = support.first();
        helptext += f;
        f = support.next();
        for (; f; f = support.next())
        {
            helptext += ',';
            if (helptext.length() - lastnl > 40)
            {
                helptext += "\n  ";
                lastnl = helptext.length() - 2;
            }
            else
            {
                helptext += ' ';
            }
            helptext += f;
        }
        helptext += "\n\nCommands:\n";
        helptext += " Right-click  : Popup menu\n";
        helptext += " Ctrl Q: Quit program\n";
        helptext += " Ctrl N: Natural size of the image\n";
        helptext += " Ctrl O: Set origin to the center of the image\n";
        helptext += " Ctrl R: Restore origin to the upper left corner of the image\n";
        helptext += " Ctrl - or Ctrl M: Half the size of the image\n";
        helptext += " Ctrl + or Ctrl P: Double the size of the image\n";
        helptext += " Ctrl A: Aspect ratio\n";
        helptext += " Arrows: Moves the mask around the image\n";
        helptext += " Ctrl+Arrows: Resizes the mask\n";
        helptext += " Alt + or Ctrl P: Increases the width of the crown or frame\n";
        helptext += " Alt - or Ctrl M: Decreases the width of the crown or frame\n";
        helpmsg = new QMessageBox("Help", helptext,
                                  QMessageBox::Information, QMessageBox::Ok, 0, 0, 0, 0, FALSE);
    }
    helpmsg->show();
    helpmsg->raise();
}

void maskImg::about()
{
    QMessageBox::about(this, "XMask",
                       "Creates a binary mask out of a set of \n"
                       "Electron Microscopy Single Particle Images.\n");
}


void maskImg::aboutXmipp()
{
    QMessageBox::about(this, "Xmipp: Xmipp Image Processing Package",
                       "Biocomputing Unit.\n"
                       "National Center of Biotechnology-CSIC\n"
                       "Madrid, Spain\n"
                       "http://www.biocomp.cnb.uam.es\n");
}




/****************************************************/

maskImg* maskImg::other = 0;

/****************************************************/
