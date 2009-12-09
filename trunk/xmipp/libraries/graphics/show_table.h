/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Alberto Pascual (pascual@cnb.csic.es)
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

#ifndef TABLE_H
#define TABLE_H

#include <qglobal.h>
#include <qlabel.h>
#include <qpainter.h>
#include <qtimer.h>

#ifdef QT3_SUPPORT
//Added by qt3to4:
#include <QPixmap>
#include <QMouseEvent>
#include <QKeyEvent>
#include <q3table.h>
#include <q3popupmenu.h>
#else
#include <qtable.h>
#include <qpopupmenu.h>
#endif

#include <data/funcs.h>
#include <data/matrix2d.h>
#include <data/selfile.h>

#include <iostream>
#include <list>

/**@defgroup ShowTables Show Tables
   @ingroup GraphicsLibrary */
//@{
/** Class to show tables.
    This is a general class that is inherited by real show classes
    like ShowSel or ShowVol.
*/


#ifdef QT3_SUPPORT
// MOC_SKIP_BEGIN
class ShowTable : public Q3Table
// MOC_SKIP_END
#else
class ShowTable : public QTable
#endif

{
    Q_OBJECT

protected:

    // Axis color
    QColor lableColor;
    // Font for axes
    QFont labelFont;
    // Font Color
    QColor fontColor;

    // Filename of the file being represented
    FileName    fn;

    // If the file changes and this flag is set, the representation
    // must change
    bool        check_file_change;
    // Modification time for this file as read from the directory
    time_t      modification_time;
    // Timer to check the modification of the file
    QTimer     *timer;

    // Array with the marks (red squares) of each cell
    bool       *cellMarks;
    // Array with pointer to the pixmaps of each cell
    QPixmap   **content;
    // List to know which pixmaps are older, and probably out of use
    std::list<int>   content_queue;
    // Number of cells
    int         listSize;

    // Suggested number of rows and columns
    int         NumRows,    NumCols;
    // Size of each cell
    int         projXdim,   projYdim;
    // Actual maxWidth and maxHeight of the contents
    int         maxWidth,   maxHeight;
    // Minimum and maximum values represented in the whole content
    // if both are 0 then no global normalization is applied
    double      minPixel,   maxPixel;

    // Current scale factor
    double      currScale;

#ifdef QT3_SUPPORT
    // Rightclick menubar
    Q3PopupMenu *menubar;
    // Options of the menubar
    Q3PopupMenu *options;
#else
    QPopupMenu *menubar;
    QPopupMenu *options;
#endif

    // Status label. Need not be used
    QLabel     *status;
    // Tempfiles that need to be deleted on destruction
    std::vector<std::string> tempfilenames;
public:
    /** Empty constructor */
    ShowTable();

    // Destructor
    ~ShowTable();

    /** Modify representation if the file changes.
        Call this function before initWithFile of the inherited class.*/
    void setPoll()
    {
        check_file_change = true;
    }

    /** Only way to initialize the show... widgets */
    virtual void initWithFile(int _numRows, int _numCols,
                              const FileName &_fn, double _minGray = 0, double _maxGray = 0) = 0;
    /** Open new file.
        Use this function to change the file being represented.*/
    virtual void openNewFile(const FileName &_fn) = 0;
private:
    /* Change mark of cell */
    void changeMark(int row, int col);
protected:
    /* GUI Change Color */
    void GUIchangeColor(QColor &_color, const char * _color_title);
    // Annotate modification time of the represented file
    void annotateTime(const FileName &_fn);

    /* Connect a Timer for checking the file change */
    void connectTimer();

    /* Set all information to NULL, 0, ... No check is done
       wether pointers must be deleted or not */
    virtual void init();

    /* Delete pointers and then calls to init() */
    virtual void clear();
    /* Ask for memory for contents and cellMarks. They are initialized.
       listSize must be already set */
    virtual void initContents();
    /* Remove all contents and content_queue. Memory is still allocated
       for the arrays */
    virtual void clearContents();

    /* Initialize shape and properties of the table */
    virtual void initTable();
    /* Form the right click menu bar */
    virtual void initRightclickMenubar() = 0;

    /* Insert help and Quit in the right click menu bar*/
    void insertGeneralItemsInRightclickMenubar();
    /* Change a boolean option */
    virtual void changeBoolOption(int _mi, int _mic);
    /* Adjust label to window size */
    virtual void adjustStatusLabel();

    /* Send update to a cell */
    void updateCellIdx(int i)
    {
        int row, col;
        IndextoPos(i, row, col);
        updateCell(row, col);
    }
    /* Index of cell */
    int indexOf(int row, int col) const
    {
        return (row * numCols()) + col;
    }
    /* Position of cell */
    void IndextoPos(int _index, int& _row, int& _col) const
    {
        _col = _index % NumCols;
        _row = _index / NumCols;
    }
    /* Label of cell */
    virtual const char* cellLabel(int i) const
    {
        return NULL;
    }

    /* How to repaint the cell.
       This is the main function. It calls producePixmapAt(i) if the
       corresponding pixmap pointer is NULL. */
    virtual void paintCell(QPainter *p, int row, int col, const QRect & cr,
                           bool selected, const QColorGroup & cg);
    /* Draw Red/White frame and label in the Cell.
       The position of the label can be 0=Bottom-right corner or
       1=Top-right corner */
    virtual void drawFrameAndLabel(QPainter *p, int row, int col, int i,
                                   int label_pos = 0);
    /* This is the function that all inherited classes should implement */
    virtual void producePixmapAt(int i) = 0;
    /* Scale to currScale and produce the output minGray and maxGray
       for an ImageXmipp. This function is used by producePixmapAt.*/
    virtual void scale_and_normalize(Matrix2D<double> &I, bool normalize,
                                     int &minGray, int &maxGray);
    /* How to check if the pixmap is old or not */
    virtual void insert_content_in_queue(int i)
    {};

    /* Reopen this file */
    virtual void reOpenFile()
    {};

    /* Show the average and SD of a Selfile */
    void showStats(SelFile &SF, bool apply_geo = FALSE);

    /* Make a temp file */
    std::string makeTempFile(int &fd);
private slots:
    /* Open window with help about keys */
    void giveHelp();
    /* Link to Xmipp */
    void aboutXmipp();
protected slots:
    /* GUI for opening a file */
    virtual void GUIopenFile();
    /* Change scale */
    virtual void changeScale(double newScale);
    /* Process keys */
    virtual void keyPressEvent(QKeyEvent*);
    /* Process double clicks */
    virtual void contentsMouseDoubleClickEvent(int row, int col, int button,
            const QPoint & mousePos);
    /* Process single clicks */
    virtual void contentsMousePressEvent(int row, int col, int button,
                                         const QPoint & mousePos);
    /* Process mouse move events */
    virtual void contentsMouseMoveEvent(QMouseEvent* e)
    {
#ifdef QT3_SUPPORT
        Q3Table::contentsMouseMoveEvent(e);
#else
        QTable::contentsMouseMoveEvent(e);
#endif
    };
    /* Check the change of file */
    virtual void check_file();
    /*Change background color */
    void changeFontColor();
    /* Change Font for labels */
    void changeFont();
};
//@}

#endif // TABLE_H
