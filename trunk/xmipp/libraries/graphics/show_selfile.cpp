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
#include "show_ctf_estimate.h"

#include <reconstruction/ctf_estimate_from_micrograph.h>

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
    selstatus       = NULL;
    ShowTable::init();
}

void ShowSel::clear()
{
    if (selstatus != NULL) delete[] selstatus;
    if (imgnames  != NULL) delete[] imgnames;
    ShowTable::clear();
}

/* Initialize with a sel file.---------------------------------------------- */
void ShowSel::initWithFile(int _numRows, int _numCols,
                           const FileName &_fn, double _minGray, double _maxGray)
{
    init();
    selfile_fn = _fn;
    readFile(_fn, _minGray, _maxGray);
    if (_numRows != -1) NumRows = _numRows;
    else NumRows = FLOOR(700.0 / projYdim);
    if (NumRows==0) NumRows=1;
    if (_numCols != -1) NumCols = _numCols;
    else NumCols = FLOOR(900.0 / projXdim);
    if (NumCols==0) NumCols=1;
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
    NumRows = _numRows;
    NumCols = _numCols;
    readFile(_fn, 0, 0);
    if (_numRows == -1 || _numCols == -1)
    {
        NumCols = FLOOR(900.0 / projXdim);
        NumRows = FLOOR(700.0 / projYdim);
	if (NumRows==0) NumRows=1;
	if (NumCols==0) NumCols=1;
    }
    initTable();
    initRightclickMenubar();
    repaint();
}

void ShowSel::initWithObject(int _numRows, int _numCols,
                             SelFile &_SF, const char *_title)
{
    init();
    fn = "";
    setCaption(_title);
    _SF.go_first_ACTIVE();
    readObject(_SF);
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
    SelFile         SF(_fn);
    annotateTime(_fn);
    readObject(SF, _minGray, _maxGray);
}

void ShowSel::readObject(SelFile &SF, double _minGray, double _maxGray)
{
    listSize        = SF.ImgNo();
    if (!showonlyactive)   listSize += SF.ImgNo(SelLine::DISCARDED);
    if (listSize == 0)
        REPORT_ERROR(1, "ShowSel::readFile: Input selfile is empty");
    imgnames        = new FileName[listSize];
    selstatus       = new bool[listSize];
    initContents();
    SF.ImgSize(projYdim, projXdim);
    if (load_mode == PSD_mode && NumRows != -1 && NumCols != -1)
    {
        // Scale to make the images fit into a reasonable window
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
    while (!SF.eof())
    {
        if (SF.Is_ACTIVE() || !showonlyactive)
        {
            imgnames[i] = SF.get_current_file();
            selstatus[i] = SF.Is_ACTIVE();
            i++;
        }
        SF.next();
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
        if (load_mode == PSD_mode) xmipp2PSD(I(), I());
        double min_val, max_val;
        I().computeDoubleMinMax(min_val, max_val);
        if (first || min_val < minPixel) minPixel = min_val;
        if (first || max_val > maxPixel) maxPixel = max_val;
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

    // Recalculate CTF model
    if (load_mode == CTF_mode)
    {
        options->insertItem("Edit CTF model", this, SLOT(editCTFmodel()));
        options->insertItem("Recompute CTF model", this, SLOT(recomputeCTFmodel()));
        options->insertSeparator();

        options->setItemEnabled(mi_imgAsLabels, false);
        options->setItemEnabled(mi_selAsLabels, false);
    }

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
    if (load_mode == Normal_mode || load_mode == CTF_mode)
    {
        options->setItemEnabled(mi_showLabel, false);
        options->setItemEnabled(mi_hideLabel, true);
    }
    else if (load_mode == PSD_mode)
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
    if (options->isItemEnabled(mi_showLabel)) return NULL;
    if (load_mode == CTF_mode)
    {
        // Get the defocus parameters from the ctfparam file
        FileName fn_param = imgnames[i].without_extension() + ".ctfparam";
        try
        {
            XmippCTF ctf;
            ctf.read(fn_param, false);
            std::string defocus_val = integerToString(ROUND(XMIPP_MIN(ctf.DeltafU, ctf.DeltafV)), 6) + " " +
                                      integerToString(ROUND(XMIPP_MAX(ctf.DeltafU, ctf.DeltafV)), 6) + " " +
                                      integerToString(ABS(ROUND(ctf.DeltafU - ctf.DeltafV)));
            return defocus_val.c_str();
        }
        catch (Xmipp_error XE)
        {
            return ((std::string)"Cannot open " + fn_param).c_str();
        }
    }
    else
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
        I.read(imgnames[i], true, -1, apply_geo, FALSE);
        if (load_mode == PSD_mode) xmipp2PSD(I(), I());
    }
    else if (I.isComplexImage(imgnames[i]))
    {
        // FFT Xmipp images: plot log10(1+|I|^2)
        Image<std::complex<double> > If;
        If.read(imgnames[i]);
        FFT_magnitude(If(), I());
        FOR_ALL_ELEMENTS_IN_ARRAY2D(I())
            A2D_ELEM(I(), i, j) =
                log10(1 + A2D_ELEM(I(), i, j) * A2D_ELEM(I(), i, j));
    }
    else
        // Unknown image
        I().initZeros(projYdim, projXdim);

    // Scale and normalize
    int minGray = 0, maxGray = 0;
    scale_and_normalize(I(), options->isItemEnabled(mi_Individual_norm),
                        minGray, maxGray);

    // If PSD mode, make the full window fit the current size
    if (load_mode == PSD_mode || load_mode == CTF_mode)
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
    if (listSize < NumRows*NumCols) return;
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
    if (found) content_queue.erase(ptr);
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
    SelFile SFNew;
    bool saveFile = false;
    for (int i = 0; i < listSize; i++)
    {
        if (cellMarks[i])
        {
            saveFile = true;
            SFNew.insert(imgnames[i], SelLine::DISCARDED);
        }
        else
            if (selstatus[i]) SFNew.insert(imgnames[i], SelLine::ACTIVE);
            else     SFNew.insert(imgnames[i], SelLine::DISCARDED);
    }
    if (saveFile) writeSelFile(SFNew);
    else QMessageBox::about(this, "Error!", "No images selected\n");
}

/* This function saves the sel file with the selected images as active and
  the rest of the sel file as discarded. */
void ShowSel::saveSelFileActive()
{
    SelFile SFNew;
    bool saveFile = false;
    for (int i = 0; i < listSize; i++)
    {
        if (cellMarks[i])
        {
            saveFile = true;
            SFNew.insert(imgnames[i], SelLine::ACTIVE);
        }
        else
            SFNew.insert(imgnames[i], SelLine::DISCARDED);
    }
    if (saveFile) writeSelFile(SFNew);
    else QMessageBox::about(this, "Error!", "No images selected\n");
}

/* This function saves a new sel file with the selected images as active.*/
void ShowSel::saveSelFileNew()
{
    SelFile SFNew;
    bool saveFile = false;
    for (int i = 0; i < listSize; i++)
    {
        if (cellMarks[i])
        {
            saveFile = true;
            SFNew.insert(imgnames[i], SelLine::ACTIVE);
        }
    }
    if (saveFile) writeSelFile(SFNew);
    else QMessageBox::about(this, "Error!", "No images selected\n");
}

/* This function saves a  sel file with the selected images as active.
   It does not ask for filename, it overwrites present selfile*/
void ShowSel::saveSelFileNewOverwrite()
{
    SelFile SFNew;
    bool saveFile = false;
    for (int i = 0; i < listSize; i++)
    {
        if (cellMarks[i])
        {
            saveFile = true;
            SFNew.insert(imgnames[i], SelLine::ACTIVE);
        }
    }
    if (saveFile) writeSelFile(SFNew, true);
    else QMessageBox::about(this, "Error!", "No images selected\n");
}

/* Save a Selfile.
   Make all possible checkings */
void ShowSel::writeSelFile(SelFile &_SF, bool overwrite)
{

    if (overwrite)
        _SF.write((std::string)((const char *)selfile_fn.c_str()));
    else
    {
#ifdef QT3_SUPPORT
         QString newfilename = Q3FileDialog::getSaveFileName(
#else
         QString newfilename = QFileDialog::getSaveFileName(    
#endif
                                  selfile_fn.c_str(), "*.sel", this, "Sel files");
        if (!newfilename.isEmpty())
        {
            QFileInfo fi(newfilename);
            if (fi.extension(false) != "sel")
            {
                if (QMessageBox::information(this, "Showsel application",
                                             "The file has no ""sel"" extension. add it? ",
                                             "Yes", "No") == 0) newfilename += ".sel";
            }
            fi.setFile(newfilename);
            if (fi.exists())
                if (QMessageBox::information(this, "Showsel application",
                                             "The file already exist. Overwrite?",
                                             "Yes", "No") == 0) _SF.write((std::string)((const char *)newfilename));
                else QMessageBox::about(this, "Warning!", "Saving aborted\n");
            else _SF.write((std::string)((const char *)newfilename));
        }
        else  QMessageBox::about(this, "Warning!", "Saving aborted\n");
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
    if (options->isItemEnabled(mi_imgAsLabels)) labeltype = Filename_LABEL;
    else if (options->isItemEnabled(mi_selAsLabels)) labeltype = SFLabel_LABEL;
    changeBoolOption(mi_imgAsLabels, mi_selAsLabels);
}

// Show statistics ---------------------------------------------------------
void ShowSel::showStats()
{
    SelFile SFNew;
    for (int i = 0; i < listSize; i++)
        if (cellMarks[i])
            SFNew.insert(imgnames[i], SelLine::ACTIVE);
    if (SFNew.ImgNo()) ShowTable::showStats(SFNew, apply_geo);
    else QMessageBox::about(this, "Error!", "No images selected\n");
}

// Show Sel Stats ----------------------------------------------------------
void ShowSel::showSelStats()
{
    int total = 0;
    int active = 0;
    int discarded = 0;
    int commented = 0;
    SelFile SF(fn);
    while (!SF.eof())
    {
        // Get file
        if (SF.Is_ACTIVE())
            active++;
        else if (SF.Is_DISCARDED())
            discarded++;
        else if (SF.Is_COMMENT())
            commented++;
        SF.next();
        total++;
    }  // while
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
    Str1 += "Commented image: ";
    tmpS.setNum(commented);
    Str1 += tmpS;
    Str1 += " (";
    tmpS.setNum((float) commented*100.0 / (float) total , 'f', 2);
    Str1 += tmpS;
    Str1 += "%)";
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
    else if (load_mode == CTF_mode)
        showimg->loadImage(imgnames[i].c_str(), 0, 0, ImageViewer::CTF_mode);
    showimg->show();
}

// Edit CTF model ----------------------------------------------------------
void ShowSel::editCTFmodel()
{
    if (fn_assign == "" && fn_assign_sel == "")
    {
        QMessageBox::about(this, "Error!", "No Assign CTF file provided\n");
        return;
    }

    FileName fn_param;
    if (fn_assign != "")
    { // Single micrograph with or without pieces
        // Read the Assign CTF parameters
        Prog_assign_CTF_prm assign_ctf_prm;
        assign_ctf_prm.read(fn_assign);

        // Check if the CTF is computed at each particle
        FileName fn_root = assign_ctf_prm.image_fn.remove_all_extensions();

        // Get the piece name
        if (assign_ctf_prm.compute_at_particle)
        {
            int i = indexOf(currentRow(), currentColumn());
            fn_param = imgnames[i].without_extension() + ".ctfparam";
        }
        else if (assign_ctf_prm.micrograph_averaging)
        {
            // If it is the average of the micrograph
            if (assign_ctf_prm.PSD_mode == Prog_assign_CTF_prm::ARMA)
                fn_param = fn_root + "_ARMAavg.ctfparam";
            else fn_param = fn_root + "_Periodogramavg.ctfparam";
        }
        else
        {
            // If the micrograph was divided into pieces
            if (assign_ctf_prm.PSD_mode == Prog_assign_CTF_prm::ARMA)
                fn_param = fn_root + "_ARMA";
            else fn_param = fn_root + "_Periodogram";
            // Get the piece to edit
            int i = indexOf(currentRow(), currentColumn()) + 1;
            fn_param += integerToString(i, 5);
            fn_param += ".ctfparam";
        }
    }
    else
    {   // Multiple micrographs all with micrograph_averaging
        // Get the assign filename
        SelFile SF_assign;
        SF_assign.read(fn_assign_sel);
        SF_assign.jump(indexOf(currentRow(), currentColumn()));
        FileName fn_assign = SF_assign.get_current_file();

        // Read the corresponding assignment parameter file
        Prog_assign_CTF_prm assign_ctf_prm;
        assign_ctf_prm.read(fn_assign);
        FileName fn_root = assign_ctf_prm.image_fn.remove_all_extensions();

        // Get the piece name
        if (assign_ctf_prm.micrograph_averaging)
        {
            // If it is the average of the micrograph
            if (assign_ctf_prm.PSD_mode == Prog_assign_CTF_prm::ARMA)
                fn_param = fn_root + "_ARMAavg.ctfparam";
            else fn_param = fn_root + "_Periodogramavg.ctfparam";
        }
        else
        {
            REPORT_ERROR(1, "ShowSel::editCTFmodel: This function is intended"
                         " only for micrograph averages");
        }
    }

    // Edit the CTF
    system(((std::string)"xmipp_edit -i " + fn_param + " &").c_str());
}

// Recompute CTF model -----------------------------------------------------
void ShowSel::recomputeCTFmodel()
{
    if (fn_assign == "" && fn_assign_sel == "")
    {
        QMessageBox::about(this, "Error!", "No Assign CTF file provided\n");
        return;
    }

    FileName fn_psd;
    Prog_assign_CTF_prm assign_ctf_prm;
    if (fn_assign!="")
    {
	try
	{
            assign_ctf_prm.read(fn_assign);
	}
	catch (Xmipp_error XE)
	{
            std::cout << XE;
            std::cout << "It seems that " << fn_assign << " is not the parameter file"
            << " that you used to estimate the CTFs\n";
            return;
	}

	// Get the PSD name
	FileName fn_root = assign_ctf_prm.image_fn.remove_all_extensions();
	if (assign_ctf_prm.compute_at_particle)
	{
            int i = indexOf(currentRow(), currentColumn());
            fn_psd = imgnames[i].without_extension() + ".psd";
	}
	else if (assign_ctf_prm.micrograph_averaging)
	{
            // If it is the average of the micrograph
            if (assign_ctf_prm.PSD_mode == Prog_assign_CTF_prm::ARMA)
        	fn_psd = fn_root + "_ARMAavg.psd";
            else fn_psd = fn_root + "_Periodogramavg.psd";
	}
	else
	{
            // If the micrograph was divided into pieces
            if (assign_ctf_prm.PSD_mode == Prog_assign_CTF_prm::ARMA)
        	fn_psd = fn_root + "_ARMA";
            else fn_psd = fn_root + "_Periodogram";
            // Get the piece to recompute
            int i = indexOf(currentRow(), currentColumn()) + 1;
            fn_psd += integerToString(i, 5);
            fn_psd += ".psd";
	}
    }
    else
    {   // Multiple micrographs all with micrograph_averaging
        // Get the assign filename
        SelFile SF_assign;
        SF_assign.read(fn_assign_sel);
        SF_assign.jump(indexOf(currentRow(), currentColumn()));
        FileName fn_assign = SF_assign.get_current_file();

        // Read the corresponding assignment parameter file
        assign_ctf_prm.read(fn_assign);
        FileName fn_root = assign_ctf_prm.image_fn.remove_all_extensions();

        // Get the piece name
        if (assign_ctf_prm.micrograph_averaging)
        {
            // If it is the average of the micrograph
            if (assign_ctf_prm.PSD_mode == Prog_assign_CTF_prm::ARMA)
                fn_psd = fn_root + "_ARMAavg.psd";
            else fn_psd = fn_root + "_Periodogramavg.psd";
        }
        else
        {
            REPORT_ERROR(1, "ShowSel::recomputeCTFmodel: This function is intended"
                         " only for micrograph averages");
        }
    }

    // Show this image in a separate window to select the main parameters
    AssignCTFViewer *prm_selector = new AssignCTFViewer(fn_psd, assign_ctf_prm);
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
    if (row < 0 || col < 0) return;
    updateStatus(indexOf(row, col));
}

void ShowSel::updateStatus(int i)
{
    if (i > listSize) return;
    status->setText(imgnames[i].c_str());
}

// Set Assign CTF file -----------------------------------------------------
void ShowSel::setAssignCTFfile(const FileName &_fn_assign)
{
    fn_assign = _fn_assign;
}

// Set Assign CTF selfile --------------------------------------------------
void ShowSel::setAssignCTFselfile(const FileName &_fn_assign_sel)
{
    fn_assign_sel = _fn_assign_sel;
}
