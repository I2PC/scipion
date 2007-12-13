/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Alberto Pascual (pascual@cnb.uam.es)
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

#include "show_vol.h"
#include "show_tools.h"

/* Initialize with a volume file.------------------------------------------- */
void ShowVol::initWithFile(int _numRows, int _numCols,
                           const FileName &_fn, double _minGray, double _maxGray)
{
    init();
    readFile(_fn, _minGray, _maxGray);
    if (_numRows != -1 && _numCols != -1)
    {
        NumRows = _numRows;
        NumCols = _numCols;
    }
    else
    {
        NumCols = XMIPP_MIN(10, FLOOR(900.0 / XSIZE(V())));
        NumRows = XMIPP_MIN(10, FLOOR(700.0 / YSIZE(V())));
    }
    initTable();
    initRightclickMenubar();
    repaint();
    if (check_file_change) connectTimer();
}

/* Read a Volume ---------------------------------------------------------- */
void ShowVol::readFile(const FileName &_fn,
                       double _minGray, double _maxGray)
{
    FileName aux_fn = _fn;
    clear();

    fn = aux_fn;
    if (fn[fn.length()-1] == 'x')
    {
        slices = 'X';
        aux_fn = fn.substr(0, fn.length() - 1);
    }
    else if (fn[fn.length()-1] == 'y')
    {
        slices = 'Y';
        aux_fn = fn.substr(0, fn.length() - 1);
    }
    else slices = 'Z';

    wait_until_stable_size(aux_fn);
    V.read(aux_fn);
    setCaption(fn.c_str());
    annotateTime(aux_fn);

    V().setXmippOrigin();
    if (_minGray == 0 && _maxGray == 0)
        V().computeDoubleMinMax(minPixel, maxPixel);
    else
    {
        minPixel = _minGray;
        maxPixel = _maxGray;
    }
    switch (slices)
    {
    case 'X':
        listSize = XSIZE(V());
        projYdim = ZSIZE(V());
        projXdim = YSIZE(V());
        break;
    case 'Y':
        listSize = YSIZE(V());
        projYdim = ZSIZE(V());
        projXdim = XSIZE(V());
        break;
    case 'Z':
        listSize = ZSIZE(V());
        projYdim = YSIZE(V());
        projXdim = XSIZE(V());
        break;
    }
    if (listSize == 0)
        REPORT_ERROR(1, "ShowVol::readFile: Input volume is empty");
    initContents();
}

/* Init table -------------------------------------------------------------- */
void ShowVol::initTable()
{
    ShowTable::initTable();
    setFocusPolicy(NoFocus);   // no keyboard focus is accepted
    // Really set size
    setMaximumSize(maxWidth, maxHeight);
    resize(maxWidth, maxHeight);
}

/* Init Rightclick menubar ------------------------------------------------- */
void ShowVol::initRightclickMenubar()
{
    menubar = new QPopupMenu();

    // File ................................................................
    QPopupMenu * file = new QPopupMenu();
    file->insertItem("Open...", this,  SLOT(GUIopenFile()));

    // Form the menu
    menubar->insertItem("&File", file);
    insertGeneralItemsInRightclickMenubar();
}

/* Produce pixmap ---------------------------------------------------------- */
void ShowVol::producePixmapAt(int i)
{
    int log_i;
    switch (slices)
    {
    case 'X':
        log_i = i + STARTINGX(V());
        break;
    case 'Y':
        log_i = i + STARTINGY(V());
        break;
    case 'Z':
        log_i = i + STARTINGZ(V());
        break;
    }
    Image I;
    V().getSlice(log_i, I(), slices);

    int minGray, maxGray;
    scale_and_normalize(I(), true, minGray, maxGray);

    content[i] = new QPixmap;
    xmipp2Pixmap(I, content[i], minGray, maxGray);
}

/* Show voxel value -------------------------------------------------------- */
void ShowVol::contentsMouseMoveEvent(QMouseEvent* e)
{
    QPoint Pos = e->pos(); // extract pointer position
    int row = rowAt(Pos.y());
    int col = columnAt(Pos.x());
    if (row < 0 || col < 0) return;
    int k, i, j;
    switch (slices)
    {
    case 'X':
        k = Pos.y() - rowPos(row);
        i = Pos.x() - columnPos(col);
        j = indexOf(row, col);
        break;
    case 'Y':
        k = Pos.y() - rowPos(row);
        i = indexOf(row, col);
        j = Pos.x() - columnPos(col);
        break;
    case 'Z':
        k = indexOf(row, col);
        i = Pos.y() - rowPos(row);
        j = Pos.x() - columnPos(col);
        break;
    }
    updateStatus(k, i, j);
}

/* Open new file ----------------------------------------------------------- */
void ShowVol::openNewFile(const FileName &_fn)
{
    fn = _fn;
    reOpenFile();
}

/* Reopen file ------------------------------------------------------------- */
void ShowVol::reOpenFile()
{
    readFile(fn);
    initTable();
    repaint();
}

/* Update status ----------------------------------------------------------- */
void ShowVol::updateStatus(int k, int i, int j)
{
    QString message, moremsg;
    int k_log, i_log, j_log;
    V().toLogical(k, i, j, k_log, i_log, j_log);
    if (!V().outside(k_log, i_log, j_log))
    {
        moremsg.sprintf("(%d,%d,%d)=(%d,%d,%d)= %.3f",
                        k, i, j, k_log, i_log, j_log,
                        V(k_log, i_log, j_log));
        message += moremsg;
    }
    status->setText(message);
}
