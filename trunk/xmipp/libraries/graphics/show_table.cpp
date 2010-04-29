/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.csic.es)
 *              Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Part of this module has been developed by Lorenzo Zampighi and Nelson Tang
 * Dept. Physiology of the David Geffen School of Medistd::cine
 * Univ. of California, Los Angeles.
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

/* FALTA:
 10- Cambiar el tama�o del grid por menu
        11- Agregar una barra de estado para ver los pixeles originales...
 12- Arreglar el escalado cuando la normalizacion es global...
*/

#include "show_table.h"
#include "show_2d.h"

#include <qpainter.h>
#include <qnamespace.h>
#include <qscrollbar.h>
#include <qmessagebox.h>
#include <qmenubar.h>
#include <qfontdialog.h>
#include <qcolordialog.h>

#ifdef QT3_SUPPORT
#include <q3filedialog.h>
//Added by qt3to4:
#include <QKeyEvent>
#include <QLabel>
#include <QPixmap>
#include <Q3Frame>
#include <Q3PopupMenu>
#else
#include <qfiledialog.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/* Constructor/Destructor -------------------------------------------------- */
ShowTable::ShowTable()
{
    check_file_change = false;
    init();
}

ShowTable::~ShowTable()
{
    clear();
    if (menubar != NULL) delete menubar;
    if (options != NULL) delete options;
    if (status  != NULL) delete status;
    for (int i = 0; i < tempfilenames.size(); i++)
        remove (tempfilenames[i].c_str());
    tempfilenames.clear();
}

void ShowTable::init()
{
    listSize = projXdim = projYdim = 0;
    maxWidth = maxHeight = 0;
    minPixel = maxPixel = 0;
    currScale       = 100;
    fn              = "";
    cellMarks       = NULL;
    content         = NULL;
    options         = NULL;
    menubar         = NULL;
    status          = NULL;
    QFont tmpFont("Fixed", 12);
    fontColor       = Qt::white;
    content_queue.clear();
}

void ShowTable::clear()
{
    if (content != NULL)
    {
        clearContents();
        delete[] content;
    }
    if (cellMarks != NULL) delete[] cellMarks;
    if (options   != NULL) delete options;
    if (menubar   != NULL) delete menubar;
    init();
}

void ShowTable::initContents()
{
    cellMarks       = new bool[listSize];
    content         = new QPixmap* [listSize];
    for (int i = 0; i < listSize; i++)
    {
        content[i] = NULL;
        cellMarks[i] = false;
    }
}

void ShowTable::clearContents()
{
    for (int i = 0; i < listSize; i++)
        if (content[i] != NULL)
        {
            delete content[i];
            content[i] = NULL;
        }
    content_queue.clear();
}

void ShowTable::connectTimer()
{
    timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(check_file()));
    timer->start(500);   // Check every 0.5 seconds
}

void ShowTable::annotateTime(const FileName &_fn)
{
    struct stat info;
    if (stat(_fn.c_str(), &info))
        std::cerr << "ShowTable::annotateTime: Cannot get time of file "
                  << _fn << std::endl;
    modification_time = info.st_mtime;
}

/* Initialize table shape and properties ----------------------------------- */
void ShowTable::initTable()
{
    // NumCols and NumRows are given by the user
    int maxCols, numCols = NumCols;
    int maxRows, numRows = NumRows;

    // Check size of the table
    bool vertical_scroll_needed = true;
    if (numCols*numRows >= listSize)
    {
        vertical_scroll_needed = false;
        if (numCols > listSize)
        {
            maxCols = numCols = listSize;
            maxRows = numRows = 1;
        }
        else
            maxRows = numRows = CEIL((double)listSize / numCols);
        maxCols = numCols;
    }
    else
    {
        maxCols = numCols;
        maxRows = CEIL((double)listSize / numCols);
    }

    setBackgroundMode(Qt::PaletteBase);      // set widgets background
#ifdef QT3_SUPPORT
    setSelectionMode(Q3Table::NoSelection); // do not accept Qt managed selections
#else
    setSelectionMode(QTable::NoSelection);
#endif
    setReadOnly(true);
    setNumCols(maxCols);              // set number of cols in table
    setNumRows(maxRows);       // set number of rows in table
    for (int i = 0; i < maxCols; i++)
        setColumnWidth(i, (int)(currScale / 100*projXdim));
    for (int i = 0; i < maxRows; i++)
        setRowHeight(i, (int)(currScale / 100*projYdim));

    // Set Maxwidth
    maxWidth = (int)(numCols * currScale / 100 * projXdim + 4);
    maxHeight = (int)(numRows * currScale / 100 * projYdim + 4);

    // Make room for the scroll bar
    if (vertical_scroll_needed)   maxWidth  += horizontalScrollBar()->height();

    verticalHeader()->hide();
    setLeftMargin(0);
    horizontalHeader()->hide();
    setTopMargin(0);

    // Adjust size of the status label
    adjustStatusLabel();

    // Make Table connections
    connect(this, SIGNAL(clicked(int, int, int, const QPoint &)),
            this, SLOT(contentsMousePressEvent(int, int, int,
                       const QPoint &)));

    connect(this, SIGNAL(doubleClicked(int, int, int, const QPoint &)),
            this, SLOT(contentsMouseDoubleClickEvent(int, int, int, const QPoint &)));
}
// Change colors -----------------------------------------------------------
void ShowTable::GUIchangeColor(QColor &_color, const char *_color_title)
{
    QColor tmpColor = QColorDialog::getColor(fontColor, this, _color_title);
    if (tmpColor.isValid())
    {
        _color = tmpColor;
        repaintContents();
    }
}

void ShowTable::changeFontColor()
{
    GUIchangeColor(fontColor, "Font color");
}

// Change Font -------------------------------------------------------------
void ShowTable::changeFont()
{
    QFont tmpFont;
    bool ok;
    tmpFont = QFontDialog::getFont(&ok, labelFont, this, "Font type");
    if (ok)
    {
        labelFont = tmpFont;
        repaintContents();
    }
}

/* Init Rightclick menubar ------------------------------------------------- */
void ShowTable::insertGeneralItemsInRightclickMenubar()
{
    //Colors................................................................
#ifdef QT3_SUPPORT
    Q3PopupMenu* colorMenu = new Q3PopupMenu();
#else
    QPopupMenu* colorMenu = new QPopupMenu();
#endif
    colorMenu->insertItem("&Font color", this, SLOT(changeFontColor()));
    menubar->insertItem("&Change color...", colorMenu);
    menubar->insertItem("&Change font...", this, SLOT(changeFont()));
    // Help ................................................................
#ifdef QT3_SUPPORT
    Q3PopupMenu* help = new Q3PopupMenu();
#else
    QPopupMenu* help = new QPopupMenu();
#endif
    help->insertItem("About &Xmipp", this, SLOT(aboutXmipp()));
    help->insertSeparator();
    help->insertItem("Help!", this, SLOT(giveHelp()));
    menubar->insertItem("&Help", help);

    // Quit ................................................................
    menubar->insertItem("&Quit", this,  SLOT(close()));
}


/* Change boolean option --------------------------------------------------- */
void ShowTable::changeBoolOption(int _mi, int _mic)
{
    bool show = options->isItemEnabled(_mi);
    if (show)
    {
        options->setItemEnabled(_mi,  false);
        options->setItemEnabled(_mic, true);
    }
    else
    {
        options->setItemEnabled(_mi,  true);
        options->setItemEnabled(_mic, false);
    }
    repaintContents();
}

/* Adjust label to window size --------------------------------------------- */
void ShowTable::adjustStatusLabel()
{
    if (status == NULL) status = new QLabel(this);
#ifdef QT3_STATUS
    status->setFrameStyle(Q3Frame::WinPanel | Q3Frame::Sunken);
#else
    status->setFrameStyle(QFrame::WinPanel | QFrame::Sunken);
#endif
    status->setFixedHeight(fontMetrics().height() + 4);
    status->setFixedWidth(maxWidth);
    status->move(0, maxHeight);
    status->show();
    maxHeight += status->height();
}

/* Rewrite cell painting --------------------------------------------------- */
void ShowTable::paintCell(QPainter *p, int row, int col, const QRect & cr,
                          bool selected, const QColorGroup & cg)
{
    int i = indexOf(row, col);
    if (i >= listSize) return;
    if (content[i] == NULL) producePixmapAt(i);
    insert_content_in_queue(i);
    p->drawPixmap(0, 0, *(content[i]));
    drawFrameAndLabel(p, row, col, i);
}

void ShowTable::drawFrameAndLabel(QPainter *p, int row, int col, int i,
                                  int label_position)
{
    int w = columnWidth(col);            // width of cell in pixels
    int h = rowHeight(row);           // height of cell in pixels
    int x2 = w - 1;
    int y2 = h - 1;

    // Draw frame
    if (cellMarks[i])
    {
        // if the cell is marked
        QPen pen;
        pen.setColor(Qt::red);
        pen.setWidth(3);
        p->setPen(pen);            // paint it in red
        p->drawRect(0, 0, x2, y2);        // draw rect. along cell edges
    }
    else
    {
        p->setPen(Qt::white);          // restore to normal
        p->drawLine(x2, 0, x2, y2);       // draw vertical line on right
        p->drawLine(0, y2, x2, y2);       // draw horiz. line at bottom
    }

    // Draw label
    if (cellLabel(i) != NULL)
    {
        QPen pen(fontColor);
        p->setPen(pen);
        p->setFont(labelFont);
        switch (label_position)
        {
        case 0:
            p->drawText(0, 0, x2, y2,
                        (Qt::AlignBottom | Qt::AlignRight), cellLabel(i));
            break;
        case 1:
            p->drawText(0, 0, x2, y2,
                        (Qt::AlignTop | Qt::AlignRight), cellLabel(i));
            break;
        }
    }
}

/* Scale and normalize an image -------------------------------------------- */
void ShowTable::scale_and_normalize(Matrix2D<double> &I, bool normalize,
                                    int &minGray, int &maxGray)
{
    // Care about size
    if (currScale != 100)
        I.selfScaleToSize((int)(currScale / 100*projYdim),
                          (int)(currScale / 100*projXdim));

    // Care about the normalization
    minGray = 0;
    maxGray = 255;
    if (normalize)
    {
        // There is a global normalization
        // recompute the min and max values for this image
        double min_val, max_val, slope;
        I.computeDoubleMinMax(min_val, max_val);
        if (minPixel != maxPixel)
            slope = (double)(255) / (double)(maxPixel - minPixel);
        else
            slope = 0;
        minGray = (int)(slope * (double)(min_val - minPixel));
        maxGray = (int)(slope * (double)(max_val - minPixel));
    }
}

/* Handle mouse click/keyboard events -------------------------------------- */
void ShowTable::changeMark(int row, int col)
{
    cellMarks[indexOf(row, col)] = !cellMarks[indexOf(row, col)];
    updateCell(row, col);   // show new current cell
}

void ShowTable::changeScale(double newScale)
{
    if (newScale == currScale) return;
    currScale = newScale;
    clearContents();
    initTable();
    repaintContents();
}

void ShowTable::contentsMouseDoubleClickEvent(int row, int col, int button,
        const QPoint & mousePos)
{
    changeMark(row, col);
}

void ShowTable::contentsMousePressEvent(int row, int col, int button,
                                        const QPoint & mousePos)
{
    setCurrentCell(row, col);
    if (button == Qt::RightButton)
        menubar->exec(mousePos);
}

void ShowTable::keyPressEvent(QKeyEvent* e)
{
    switch (e->key())
    {
        // Look at the key code
    case Qt::Key_F1:
        menubar->exec();
        break;
    case Qt::Key_Space:
        changeMark(currentRow(), currentColumn());
        break;
    case Qt::Key_M:
    case Qt::Key_Minus:
        if (e->state() == Qt::ControlButton) // If 'Ctrol+'-key,
            if (currScale > 10) changeScale(currScale - 10);
        break;
    case Qt::Key_P:
    case Qt::Key_Plus:
        if (e->state() == Qt::ControlButton) // If 'Ctrol+'-key,
            changeScale(currScale + 10);
        break;
    case Qt::Key_Q:
        if (e->state() == Qt::ControlButton) // If 'Ctrol Q' key,
            exit(0); // Terminate program
        break;
    default:
#ifdef QT3_SUPPORT
        Q3Table::keyPressEvent(e);
#else
        QTable::keyPressEvent(e);
#endif
        break;
    }
}

/* Open File --------------------------------------------------------------- */
void ShowTable::GUIopenFile()
{
#ifdef QT3_SUPPORT
    QString newfilename = Q3FileDialog::getOpenFileName(QString::null, "*", this, "Sel files");
#else
    QString newfilename = QFileDialog::getOpenFileName(QString::null, "*", this, "Sel files");
#endif
    if (!newfilename.isEmpty())
        openNewFile((std::string)((const char *)newfilename));
}

/* Help -------------------------------------------------------------------- */
void ShowTable::giveHelp()
{
    QString helptext = "xmipp_showsel\n";
    helptext += "\nMouse and Keyboard Commands:\n";
    helptext += " Right-click  : Popup menu\n";
    helptext += " Ctrl - or Ctrl M: Half the size of the images\n";
    helptext += " Ctrl + or Ctrl P: Double the size of the images\n";
    helptext += " Spacebar     : Mark/Unmark individual images\n";
    QMessageBox * helpmsg = new QMessageBox("Help", helptext,
                                            QMessageBox::Information, QMessageBox::Ok, 0, 0, 0, 0, FALSE);
    helpmsg->show();
    helpmsg->raise();
}

void ShowTable::aboutXmipp()
{
    QMessageBox::about(this, "Xmipp: Xmipp Image Processing Package",
                       "Biocomputing Unit.\n"
                       "National Center of Biotechnology-CSIC\n"
                       "Madrid, Spain\n"
                       "http://www.biocomp.cnb.uam.es\n");
}

/* Check file change -------------------------------------------------------- */
void ShowTable::check_file()
{
    struct stat info;
    static bool message_shown = false;
    if (stat(fn.c_str(), &info) && !message_shown)
    {
        std::cerr << "check_file: Cannot get time of file " << fn << std::endl;
        message_shown = true;
    }
    if (info.st_mtime != modification_time) reOpenFile();
}

/* Show Average and Standard deviation of a SelFile ------------------------ */
void ShowTable::showStats(SelFile &SF, bool apply_geo)
{
    try
    {
        Image _ave, _sd;
        double _minPixel, _maxPixel;
        SF.go_beginning();
        SF.get_statistics(_ave, _sd, _minPixel, _maxPixel, apply_geo);
        ImageViewer *wavg = new ImageViewer(&_ave, "Average Image");
        ImageViewer *wsd  = new ImageViewer(&_sd, "SD Image");
        wavg->show();
        wsd->show();
    }
    catch (Xmipp_error XE)
    {
        QMessageBox::about(this, "Error!",
                           "There is a problem opening files\n");
    }
}

/* Make a temporary file; returns a string which is the filename, and an
   OPEN file descriptor (passed in as a parameter) which must be closed
   by the caller AFTER writing whatever contents are needed.  The file
   will be automatically deleted when this ShowTable object is destroyed. */
std::string ShowTable::makeTempFile(int &fd)
{
    std::string tmpfname = (std::string) P_tmpdir + "/xmipp_tempfile_XXXXXX";
    char *tempfile = strdup(tmpfname.c_str());
    fd = mkstemp(tempfile);
    tmpfname = tempfile;
    free(tempfile);
    tempfilenames.push_back(tmpfname);
    return (tmpfname);
}
