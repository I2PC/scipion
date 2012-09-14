/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Alberto Pascual (pascual@cnb.csic.es)
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

#include "show_selfile.h"
#include "show_2d.h"
#include "show_tools.h"
#include <data/args.h>
#include <qmessagebox.h>

#ifdef QT3_SUPPORT
#include <q3filedialog.h>
//Added by qt3to4:
#include <QPixmap>
#include <Q3PopupMenu>
#include <QMouseEvent>
#else
#include <qfiledialog.h>
#endif

/* Empty constructor ------------------------------------------------------- */
ShowSel::ShowSel(): ShowTable()
{
    load_mode       = Normal_mode;
}

/* Init/Clear data --------------------------------------------------------- */
void ShowSel::init()
{
    labeltype       = SFLabel_LABEL;
    imgnames        = NULL;
    imgids          = NULL;
    selstatus       = NULL;
    ShowTable::init();
}

void ShowSel::clear()
{
    delete[] selstatus;
    delete[] imgnames;
    delete[] imgids;
    ShowTable::clear();
}

/* Initialize with a sel file.---------------------------------------------- */
void ShowSel::initWithFile(int _numRows, int _numCols,
                           const FileName &_fn, double _minGray, double _maxGray)
{
    init();
    selfile_fn = _fn;
    readFile(_fn, _minGray, _maxGray);
    NumRows = XMIPP_MAX(((_numRows != -1) ? _numRows : FLOOR(700.0 / projYdim)), 1);
    NumCols = XMIPP_MAX(((_numCols != -1) ? _numCols : FLOOR(900.0 / projXdim)), 1);
    initTable();
    initRightclickMenubar();
    repaint();
}

void ShowSel::initWithFile(int _numRows, int _numCols,
                           const FileName &_fn, TLoadMode _load_mode)
{
    init();
    selfile_fn = _fn;
    load_mode = _load_mode;
    readFile(_fn, 0, 0);
    NumRows = XMIPP_MAX(((_numRows != -1) ? _numRows : FLOOR(700.0 / projYdim)), 1);
    NumCols = XMIPP_MAX(((_numCols != -1) ? _numCols : FLOOR(900.0 / projXdim)), 1);
    initTable();
    initRightclickMenubar();
    repaint();
}

void ShowSel::initWithObject(int _numRows, int _numCols,
                             MetaData &_SF, const char *_title)
{
    init();
    fn = "";
    setCaption(_title);
    mdInput = _SF;
    readObject(mdInput);
    NumRows = _numRows;
    NumCols = _numCols;
    initTable();
    initRightclickMenubar();
    repaint();
}

/* Read a Selfile ---------------------------------------------------------- */
void ShowSel::readFile(const FileName &_fn, double _minGray, double _maxGray)
{
    clear();
    fn              = _fn;
    setCaption(fn.c_str());
    readSelFile(_fn, _minGray, _maxGray);
}

void ShowSel::readSelFile(const FileName &_fn, double _minGray, double _maxGray)
{
    mdInput.read(_fn);
    annotateTime(_fn);
    readObject(mdInput, _minGray, _maxGray);
}

void ShowSel::readObject(MetaData &SF, double _minGray, double _maxGray)
{
    if (showonlyactive && SF.containsLabel(MDL_ENABLED))
        SF.removeObjects(MDValueEQ(MDL_ENABLED, -1));
    listSize        = SF.size();
    if (listSize == 0)
        REPORT_ERROR(ERR_IO_SIZE, "ShowSel::readFile: Input selfile is empty");
    imgnames        = new FileName[listSize];
    imgids          = new size_t[listSize];
    selstatus       = new bool[listSize];
    initContents();
    int Zdim;
    size_t Ndim;
    getImageSize(SF, projXdim, projYdim, Zdim, Ndim);
    if (load_mode == PSD_mode && NumRows != -1 && NumCols != -1)
    {
        // Scale to make the images fit into a reasonable selfWindow
        double suggested_Xdim = XMIPP_MIN(900.0 / NumCols, projXdim);
        double suggested_Ydim = XMIPP_MIN(700.0 / NumRows, projYdim);
        double scale_X = suggested_Xdim / projXdim;
        double scale_Y = suggested_Ydim / projYdim;
        double scale = XMIPP_MIN(scale_X, scale_Y);
        projYdim = FLOOR(scale * projYdim);
        projXdim = FLOOR(scale * projXdim);
    }
    minPixel = _minGray;
    maxPixel = _maxGray;
    int i = 0;
    int enabled;
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        imgids[i] = __iter.objId;
        SF.getValue(MDL_IMAGE, imgnames[i],__iter.objId);
        if (SF.getValue(MDL_ENABLED, enabled,__iter.objId))
            selstatus[i] = enabled == 1;
        else
            selstatus[i] = true;
        i++;
    }
}

/* Compute global normalization params ------------------------------------- */
void ShowSel::compute_global_normalization_params()
{
    bool first = true;
    for (int i = 0; i < listSize; i++)
    {
        Image<double> I;
        I.read(imgnames[i]);
        if (load_mode == PSD_mode)
            xmipp2PSD(I(), I());
        double min_val, max_val;
        I().computeDoubleMinMax(min_val, max_val);
        if (first || min_val < minPixel)
            minPixel = min_val;
        if (first || max_val > maxPixel)
            maxPixel = max_val;
        first = false;
    }
}

/* Init table -------------------------------------------------------------- */
void ShowSel::initTable()
{
    ShowTable::initTable();
#ifdef QT3_SUPPORT

    setFocusPolicy(Qt::StrongFocus);   // keyboard focus is accepted
#else

    setFocusPolicy(StrongFocus);
#endif
    // Really set size
    setMaximumSize(maxWidth, maxHeight);
    resize(maxWidth, maxHeight);
}

/* Init Rightclick menubar ------------------------------------------------- */
void ShowSel::initRightclickMenubar()
{
#ifdef QT3_SUPPORT
    menubar = new Q3PopupMenu();
#else

    menubar = new QPopupMenu();
#endif

    setFileRightclickMenubar();

    // Options .............................................................
#ifdef QT3_SUPPORT

    options =  new Q3PopupMenu();
#else

    options = new QPopupMenu();
#endif

    setCommonOptionsRightclickMenubar();

    // What kind of labels
    mi_imgAsLabels = options->insertItem("Show Image Names as Labels", this,  SLOT(changeLabels()));
    mi_selAsLabels = options->insertItem("Show Sel status as Labels", this,  SLOT(changeLabels()));
    options->setItemEnabled(mi_imgAsLabels, false);
    options->setItemEnabled(mi_selAsLabels, true);
    labeltype = Filename_LABEL;
    options->insertSeparator();

    // Statistics
    options->insertItem("View average and SD Images", this,  SLOT(showStats()));
    options->insertItem("Show Sel Statistics", this,  SLOT(showSelStats()));
    options->insertSeparator();

    // Show this image
    options->insertItem("Reload all", this, SLOT(reloadAll()));
    options->insertItem("Show this image separately", this, SLOT(showThisImage()));
    options->insertSeparator();

    // Form the menu
    menubar->insertItem("&Options", options);
    menubar->insertSeparator();
    insertGeneralItemsInRightclickMenubar();
}

void ShowSel::setFileRightclickMenubar()
{
#ifdef QT3_SUPPORT
    Q3PopupMenu * file = new Q3PopupMenu();
    Q3PopupMenu * fileSave = new Q3PopupMenu();
#else

    QPopupMenu * file = new QPopupMenu();
    QPopupMenu * fileSave = new QPopupMenu();
#endif

    file->insertItem("Open...", this,  SLOT(GUIopenFile()));

    fileSave->insertItem("As discarded...",
                         this, SLOT(saveSelFileDiscarded()));
    fileSave->insertItem("As active and the rest as discarded...",
                         this, SLOT(saveSelFileActive()));
    fileSave->insertItem("In a new sel file...",
                         this, SLOT(saveSelFileNew()));
    fileSave->insertItem("Overwrite selfile with active images",
                         this, SLOT(saveSelFileNewOverwrite()));
    file->insertItem("&Save Selected Images in a Sel File", fileSave);
    menubar->insertItem("&File", file);
}

void ShowSel::setCommonOptionsRightclickMenubar()
{
    // Normalization
    mi_Individual_norm = options->insertItem("Individual Normalization",
                         this, SLOT(changeNormalize()));
    mi_Global_norm = options->insertItem("Global Normalization",
                                         this, SLOT(changeNormalize()));
    options->setItemEnabled(mi_Individual_norm, false);
    options->setItemEnabled(mi_Global_norm, true);
    options->insertSeparator();

    // Show/Hide labels
    mi_showLabel = options->insertItem("Show Labels", this,  SLOT(changeShowLabels()));
    mi_hideLabel = options->insertItem("Hide Labels", this,  SLOT(changeShowLabels()));
    if (load_mode == PSD_mode)
    {
        options->setItemEnabled(mi_showLabel, true);
        options->setItemEnabled(mi_hideLabel, false);
    }
    options->insertSeparator();

    // Select/unselect
    options->insertItem("Select All", this,  SLOT(SelectAll()));
    options->insertItem("Unselect All", this,  SLOT(unSelectAll()));
    options->insertSeparator();
}

/* Cell label -------------------------------------------------------------- */
const char * ShowSel::cellLabel(int i) const
{
    if (options->isItemEnabled(mi_showLabel))
        return NULL;
    switch (labeltype)
    {
    case SFLabel_LABEL:
        return (selstatus[i]) ? "1" : "-1";
    case Filename_LABEL:
        return imgnames[i].c_str();
    }
}

/* Produce pixmap ---------------------------------------------------------- */
void ShowSel::producePixmapAt(int i)
{
    Image<double> I;
    // Read image
    if (I.isRealImage(imgnames[i]))
    {
        // Plain Xmipp images
        if (apply_geo)
            I.readApplyGeo(imgnames[i], mdInput, imgids[i]);
        else
            I.read(imgnames[i]);
        if (load_mode == PSD_mode)
            xmipp2PSD(I(), I());
    }
    else
        // Unknown image
        I().initZeros(projYdim, projXdim);

    // Scale and normalize
    int minGray = 0, maxGray = 0;
    scale_and_normalize(I(), options->isItemEnabled(mi_Individual_norm),
                        minGray, maxGray);

    // If PSD mode, make the full selfWindow fit the current size
    if (load_mode == PSD_mode)
    {
        selfScaleToSize(LINEAR, I(), rowHeight(0), columnWidth(0));
    }

    // Convert Xmipp image to Pixmap
    content[i] = new QPixmap;
    xmipp2Pixmap(I, content[i], minGray, maxGray, 0, 0);
}

/* Grow older all contents ------------------------------------------------- */
void ShowSel::insert_content_in_queue(int i)
{
    if (listSize < NumRows*NumCols)
        return;
    // Check if the image i is already in the queue
    int jmax = content_queue.size();
    bool found = false;
    std::list<int>::iterator ptr = content_queue.begin();
    for (int j = 0; j < jmax; j++)
    {
        if ((*ptr) == i)
        {
            found = true;
            break;
        }
        ptr++;
    }
    if (found)
        content_queue.erase(ptr);
    content_queue.push_back(i);

    // If the queue is longer than 3 times the visible area then remove some
    // images
    if (jmax + 1 > 2*NumRows*NumCols)
    {
        int i_to_remove = content_queue.front();
        content_queue.pop_front();
        if (content[i_to_remove] != NULL)
        {
            delete content[i_to_remove];
            content[i_to_remove] = NULL;
        }
    }
}

/* Open new file ----------------------------------------------------------- */
void ShowSel::openNewFile(const FileName &_fn)
{
    init();
    readFile(_fn);
    if (options->isItemEnabled(mi_Individual_norm))
        compute_global_normalization_params();
    initTable();
    repaint();
}

/* Save SelFiles ----------------------------------------------------------- */
// This function saves the sel file with the selected images as discarded.
void ShowSel::saveSelFileDiscarded()
{
    MetaData SFnew;
    bool saveFile = false;
    size_t id;

    for (int i = 0; i < listSize; i++)
    {
        id = SFnew.addObject();
        SFnew.setValue(MDL_IMAGE,imgnames[i],id);
        if (cellMarks[i])
        {
            saveFile = true;
            SFnew.setValue(MDL_ENABLED,-1, id);
        }
        else
        {
            SFnew.setValue(MDL_ENABLED, selstatus[i] ? 1 : -1, id);
        }
    }
    if (saveFile)
        writeSelFile(SFnew);
    else
        QMessageBox::about(this, "Error!", "No images selected\n");
}

/* This function saves the sel file with the selected images as active and
  the rest of the sel file as discarded. */
void ShowSel::saveSelFileActive()
{
    MetaData SFnew;
    bool saveFile = false;
    size_t id;

    for (int i = 0; i < listSize; i++)
    {
        id = SFnew.addObject();
        SFnew.setValue(MDL_IMAGE,imgnames[i], id);
        if (cellMarks[i])
        {
            saveFile = true;
            SFnew.setValue(MDL_ENABLED,1, id);
        }
        else
            SFnew.setValue(MDL_ENABLED,-1, id);
    }
    if (saveFile)
        writeSelFile(SFnew);
    else
        QMessageBox::about(this, "Error!", "No images selected\n");
}

/* This function saves a new sel file with the selected images as active.*/
void ShowSel::saveSelFileNew()
{
    MetaData SFnew;
    bool saveFile = false;
    size_t id;

    for (int i = 0; i < listSize; i++)
    {
        id = SFnew.addObject();
        SFnew.setValue(MDL_IMAGE,imgnames[i], id);
        if (cellMarks[i])
        {
            saveFile = true;
            SFnew.setValue(MDL_ENABLED,1, id);
        }
    }
    if (saveFile)
        writeSelFile(SFnew);
    else
        QMessageBox::about(this, "Error!", "No images selected\n");
}

/* This function saves a  sel file with the selected images as active.
   It does not ask for filename, it overwrites present selfile*/
void ShowSel::saveSelFileNewOverwrite()
{
    MetaData SFnew;
    bool saveFile = false;
    size_t id;
    for (int i = 0; i < listSize; i++)
    {
        id = SFnew.addObject();
        SFnew.setValue(MDL_IMAGE, imgnames[i], id);
        if (cellMarks[i])
        {
            saveFile = true;
            SFnew.setValue(MDL_ENABLED, 1, id);
        }
    }
    if (saveFile)
        writeSelFile(SFnew, true);
    else
        QMessageBox::about(this, "Error!", "No images selected\n");
}

/* Save a Selfile.
   Make all possible checkings */
void ShowSel::writeSelFile(MetaData &_SF, bool overwrite)
{

    if (overwrite)
        _SF.write(selfile_fn);
    else
    {
#ifdef QT3_SUPPORT
        QString newfilename = Q3FileDialog::getSaveFileName(
                                  selfile_fn.c_str(), "*.sel", this, "Sel files");
#else

        QString newfilename = QFileDialog::getSaveFileName(
                                  selfile_fn.c_str(), "*.sel", this, "Sel files");
#endif

        if (!newfilename.isEmpty())
        {
            QFileInfo fi(newfilename)
            ;
            if (fi.extension(false) != "sel")
            {
                if (QMessageBox::information(this, "Showsel application",
                                             "The file has no ""sel"" extension. add it? ",
                                             "Yes", "No") == 0)
                    newfilename += ".sel";
            }
            fi.setFile(newfilename);
            if (fi.exists())
                if (QMessageBox::information(this, "Showsel application",
                                             "The file already exist. Overwrite?",
                                             "Yes", "No") == 0)
                    _SF.write((std::string)((const char *)newfilename));
                else
                    QMessageBox::about(this, "Warning!", "Saving aborted\n");
            else
                _SF.write((std::string)((const char *)newfilename));
        }
        else
            QMessageBox::about(this, "Warning!", "Saving aborted\n");
    }
}

// Change options ----------------------------------------------------------
void ShowSel::changeNormalize()
{
    bool indivNorm = options->isItemEnabled(mi_Individual_norm);
    if (indivNorm)
    {
        options->setItemEnabled(mi_Individual_norm, false);
        options->setItemEnabled(mi_Global_norm, true);
    }
    else
    {
        options->setItemEnabled(mi_Individual_norm, true);
        options->setItemEnabled(mi_Global_norm, false);
        if (minPixel == 0 && maxPixel == 0)
            compute_global_normalization_params();
    }
    clearContents();
    repaintContents();
}

void ShowSel::changeShowLabels()
{
    changeBoolOption(mi_showLabel, mi_hideLabel);
}
void ShowSel::changeLabels()
{
    if (options->isItemEnabled(mi_imgAsLabels))
        labeltype = Filename_LABEL;
    else if (options->isItemEnabled(mi_selAsLabels))
        labeltype = SFLabel_LABEL;
    changeBoolOption(mi_imgAsLabels, mi_selAsLabels);
}

// Show statistics ---------------------------------------------------------
void ShowSel::showStats()
{
    MetaData SFnew;
    MDRow row;
    for (int i = 0; i < listSize; i++)
    {
        if (cellMarks[i])
        {
            mdInput.getRow(row,imgids[i]);
            size_t id=SFnew.addRow(row);
            SFnew.setValue(MDL_ENABLED,1, id);
        }
    }
    if (SFnew.size())
        ShowTable::showStats(SFnew, apply_geo);
    else
        QMessageBox::about(this, "Error!", "No images selected\n");
}

// Show Sel Stats ----------------------------------------------------------
void ShowSel::showSelStats()
{
    int total = 0;
    int active = 0;
    int discarded = 0;
    int commented = 0;
    int enabled;
    MetaData SF(fn);
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        if (!SF.getValue(MDL_ENABLED, enabled,__iter.objId))
            enabled=1;
        if (enabled==1)
            active++;
        else if (enabled==-1)
            discarded++;
        total++;
    }
    QString tmpS, Str1;
    Str1 = "Sel File Name : ";
    Str1 += fn.c_str();
    Str1 += "\nTotal number of images : ";
    tmpS.setNum(total);
    Str1 += tmpS;
    Str1 += " (100.00%)\n";
    Str1 += "Active images    : ";
    tmpS.setNum(active);
    Str1 += tmpS;
    Str1 += " (";
    tmpS.setNum((float) active*100.0 / (float) total , 'f', 2);
    Str1 += tmpS;
    Str1 += "%)\n";
    Str1 += "Discarded image: ";
    tmpS.setNum(discarded);
    Str1 += tmpS;
    Str1 += " (";
    tmpS.setNum((float) discarded*100.0 / (float) total , 'f', 2);
    Str1 += tmpS;
    Str1 += "%)\n";
    QMessageBox::about((QWidget*)this, "Sel File Statistics", Str1);
}

// Reload all --------------------------------------------------------------
void ShowSel::reloadAll()
{
    clearContents();
    repaint();
}

// Show This image ---------------------------------------------------------
void ShowSel::showThisImage()
{
    int row = currentRow();
    int col = currentColumn();
    int i = indexOf(row, col);

    ImageViewer *showimg = new ImageViewer(imgnames[i].c_str(), false);
    if (load_mode == Normal_mode)
        showimg->loadImage(imgnames[i].c_str());
    else if (load_mode == PSD_mode)
        showimg->loadImage(imgnames[i].c_str(), 0, 0, ImageViewer::PSD_mode);
    showimg->show();
}

// Unselect/Select all -----------------------------------------------------
void ShowSel::SelectAll()
{
    for (int i = 0; i < listSize; i++)
        if (!cellMarks[i])
        {
            cellMarks[i] = true;
            updateCellIdx(i);
        }
}

void ShowSel::unSelectAll()
{
    for (int i = 0; i < listSize; i++)
        if (cellMarks[i])
        {
            cellMarks[i] = false;
            updateCellIdx(i);
        }
}

// Update status -----------------------------------------------------------
void ShowSel::contentsMouseMoveEvent(QMouseEvent* e)
{
    QPoint Pos = e->pos(); // extract pointer position
    int row = rowAt(Pos.y());
    int col = columnAt(Pos.x());
    if (row < 0 || col < 0)
        return;
    updateStatus(indexOf(row, col));
}

void ShowSel::updateStatus(int i)
{
    if (i > listSize)
        return;
    status->setText(imgnames[i].c_str());
}
