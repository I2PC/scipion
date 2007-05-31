/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Alberto Pascual (pascual@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Part of this module has been developed by Lorenzo Zampighi and Nelson Tang
 * Dept. Physiology of the David Geffen School of Medicine
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include "show_som.h"
#include "show_2d.h"

#include <classification/training_vector.h>

#include <qfiledialog.h>
#include <qmessagebox.h>

/* Init/Clear data --------------------------------------------------------- */
void ShowSOM::init()
{
    SFcv        = NULL;
    hisAssigned = NULL;
    cv_errors   = NULL;
    infStr      = "";
    ShowSel::init();
}

void ShowSOM::clear()
{
    if (SFcv        != NULL) delete [] SFcv;
    if (hisAssigned != NULL) delete [] hisAssigned;
    if (cv_errors   != NULL) delete [] cv_errors;
    infStr = "";
    ShowSel::clear();
}

/* Initialize with a SOM file.---------------------------------------------- */
void ShowSOM::initWithFile(const FileName &_fn_root)
{
    init();
    readFile(_fn_root);
    initTable();
    initRightclickMenubar();
    repaint();
}

/* Read a SOM -------------------------------------------------------------- */
void ShowSOM::readFile(const FileName &_fn_root,
                       double _minGray, double _maxGray)
{
    clear();
    fn              = _fn_root;
    setCaption(fn.c_str());
    readSOMFiles(_fn_root);
    readSelFile(_fn_root + ".sel");
}

void ShowSOM::readSOMFiles(const FileName &_fn_root)
{
    FileName fn_class, fn_his, fn_err, fn_inf;
    fn_class = _fn_root + ".vs";
    fn_his  = _fn_root + ".his";
    fn_err  = _fn_root + ".err";
    fn_inf  = _fn_root + ".inf";

    // Read histogram
    if (exists(fn_his))
    {
        ifstream fh_his(fn_his.c_str());
        if (fh_his)
        {
            xmippCTVectors ts(0, true);
            fh_his >> ts;
            int imax = ts.size();
            hisAssigned = new string[imax];
            for (int i = 0; i < imax; i++)
                hisAssigned[i] = ts.theTargets[i];
        }
        fh_his.close();
    }

    // Read errors
    if (exists(fn_err))
    {
        ifstream fh_err(fn_err.c_str());
        if (fh_err)
        {
            xmippCTVectors ts(0, true);
            fh_err >> ts;
            int imax = ts.size();
            cv_errors = new string[imax];
            for (int i = 0; i < imax; i++)
                cv_errors[i] = ts.theTargets[i];
        }
        fh_err.close();
    }

    // Read inf file
    if (exists(fn_inf))
    {
        ifstream fh_inf(fn_inf.c_str());
        if (fh_inf)
        {
            string line;
            getline(fh_inf, line);
            infStr = line.c_str();
            infStr += "\n";
            while (!fh_inf.eof())
            {
                getline(fh_inf, line);
                infStr += line.c_str();
                infStr += "\n";
            }
        }
        fh_inf.close();
    }

    // Read codevectors
    if (exists(fn_class))
    {
        ifstream fh_class(fn_class.c_str());
        if (fh_class)
        {
            string line, fn;
            int dim;
            string topol, neigh;
            fh_class >> dim >> topol >> NumCols >> NumRows >> neigh;
            listSize = NumCols * NumRows;
            if (listSize == 0)
                REPORT_ERROR(1, "ShowSOM::readFile: Input file is empty");
            getline(fh_class, line);

            if (infStr == "")
            {
                infStr = "Kohonen SOM algorithm\n\n";
                infStr += "Number of variables: ";
                line = ItoA(dim);
                infStr += line.c_str();
                infStr += "\n";
                infStr += "Horizontal dimension (Xdim): ";
                line = ItoA(NumCols);
                infStr += line.c_str();
                infStr += "\n";
                infStr += "Vertical dimension (Ydim): ";
                line = ItoA(NumRows);
                infStr += line.c_str();
                infStr += "\n";
                infStr += "Topology : ";
                infStr += topol.c_str();
                infStr += "\n";
                infStr += "Neighborhood function : ";
                infStr += neigh.c_str();
                infStr += "\n";
            }

            SFcv = new SelFile[listSize];
            while (!fh_class.eof())
            {
                int row, col;
                float tmp;
                fh_class >> col >> row >> tmp;
                getline(fh_class, line);
                int i = row * NumCols + col;
                SFcv[i].insert(firstToken(line), SelLine::ACTIVE);
            }
        }
        fh_class.close();
    }
}

/* Initialize right click menubar ------------------------------------------ */
void ShowSOM::initRightclickMenubar()
{
    menubar = new QPopupMenu();
    QPopupMenu * file = new QPopupMenu();
    file->insertItem("Open...", this,  SLOT(GUIopenFile()));
    file->insertItem("Save assigned images in a sel file...",
                     this, SLOT(saveAssigned()), CTRL + Key_N);
    file->insertItem("Save assigned images in separate sel files...",
                     this, SLOT(saveAssignedSeparately()));
    menubar->insertItem("&File", file);

    // Options .............................................................
    options =  new QPopupMenu();
    setCommonOptionsRightclickMenubar();

    // What kind of labels
    mi_imgAsLabels = options->insertItem("Show Image Names as Labels", this,  SLOT(changeLabelsToImg()));
    mi_selAsLabels = options->insertItem("Show Sel status as Labels", this,  SLOT(changeLabelsToSel()));
    mi_hisAsLabels = options->insertItem("Show Histogram as Labels", this, SLOT(changeLabelsToHis()));
    mi_errAsLabels = options->insertItem("Show Quantization Errors as Labels", this, SLOT(changeLabelsToErr()));
    options->setItemEnabled(mi_imgAsLabels, true);
    options->setItemEnabled(mi_selAsLabels, true);
    options->setItemEnabled(mi_hisAsLabels, false);
    options->setItemEnabled(mi_errAsLabels, true);
    labeltype = Histogram_LABEL;
    options->insertSeparator();

    // Statistics
    options->insertItem("View average and SD of the Selected Codevectors",  this,  SLOT(showStats()));
    options->insertItem("View average and SD of the Assigned Images", this,  SLOT(showRepresentedStats()));
    options->insertItem("View average of the Assigned Images together", this, SLOT(showRepresentedAverageTogether()));
    options->insertItem("View assigned images", this,  SLOT(showRepresentedSel()), CTRL + Key_A);
    options->insertItem("View error Image", this,  SLOT(showErrorImage()), CTRL + Key_E);
    options->insertSeparator();
    options->insertItem("Show Algorithm Information", this,  SLOT(showAlgoInfo()), CTRL + Key_G);
    options->insertSeparator();

    // Insert options the menu
    menubar->insertItem("&Options", options);
    menubar->insertSeparator();

    // Inser Help and Quit
    insertGeneralItemsInRightclickMenubar();
}

/* Extract represented .---------------------------------------------------- */
void ShowSOM::extractRepresented(SelFile &SF_represented)
{
    for (int i = 0; i < listSize; i++)
    {
        if (cellMarks[i])
            SF_represented.merge(SFcv[i]);
    }
}

void ShowSOM::saveAssigned()
{
    SelFile SFNew;
    extractRepresented(SFNew);
    if (SFNew.ImgNo() != 0) writeSelFile(SFNew);
    else QMessageBox::about(this, "Error!", "No images selected\n");
}

void ShowSOM::saveAssignedSeparately()
{
    QString basename;
    for (int i = 0; i < listSize; i++)
        if (cellMarks[i] && SFcv[i].ImgNo() > 0)
        {
            if (basename.isNull())
            {
                basename =
                    QFileDialog::getSaveFileName(QString::null, QString::null,
                                                 this, "Save images",
                                                 "Enter base filename for sel files");
                if (basename.isEmpty())
                {
                    QMessageBox::about(this, "Warning!", "Saving aborted\n");
                    return;
                }
            }
            int r, c;
            IndextoPos(i, r, c);
            char buf[50];
            sprintf(buf, "_row%d_col%d", r + 1, c + 1);
            QString newfilename = basename + (QString) buf + ".sel";
            QFileInfo fi(newfilename);
            if (fi.exists() &&
                (QMessageBox::information(this, "File Exists",
                                          "File " + newfilename +
                                          " already exists. Overwrite?",
                                          "Yes", "No") != 0))
            {
                QMessageBox::about(this, "Warning!",
                                   "Skipping file " + newfilename + ".");
                continue;
            }
            SFcv[i].write(newfilename.ascii());
        }
    if (basename.isNull())
        QMessageBox::about(this, "Error!", "No images selected\n");
}



/* Change labels ----------------------------------------------------------- */
void ShowSOM::changeLabelsToImg()
{
    changeLabel(mi_imgAsLabels);
}
void ShowSOM::changeLabelsToSel()
{
    changeLabel(mi_selAsLabels);
}
void ShowSOM::changeLabelsToHis()
{
    changeLabel(mi_hisAsLabels);
}
void ShowSOM::changeLabelsToErr()
{
    changeLabel(mi_errAsLabels);
}
void ShowSOM::changeLabel(int _clicked_mi)
{
    options->setItemEnabled(mi_imgAsLabels, true);
    options->setItemEnabled(mi_selAsLabels, true);
    options->setItemEnabled(mi_hisAsLabels, true);
    options->setItemEnabled(mi_errAsLabels, true);
    options->setItemEnabled(_clicked_mi,    false);
    if (_clicked_mi == mi_imgAsLabels) labeltype = Filename_LABEL;
    else if (_clicked_mi == mi_selAsLabels) labeltype = SFLabel_LABEL;
    else if (_clicked_mi == mi_hisAsLabels) labeltype = Histogram_LABEL;
    else if (_clicked_mi == mi_errAsLabels) labeltype = Err_LABEL;
    repaintContents();
}

const char * ShowSOM::cellLabel(int i) const
{
    if (options->isItemEnabled(mi_showLabel)) return NULL;
    switch (labeltype)
    {
    case SFLabel_LABEL:
        return (selstatus[i]) ? "1" : "-1";
    case Filename_LABEL:
        return imgnames[i].c_str();
    case Histogram_LABEL:
        return hisAssigned[i].c_str();
    case Err_LABEL:
        return cv_errors[i].c_str();
    }
}

/* Show Average and SD of the represented images --------------------------- */
void ShowSOM::showRepresentedStats()
{
    SelFile SFNew;
    extractRepresented(SFNew);
    if (SFNew.ImgNo()) ShowTable::showStats(SFNew, apply_geo);
    else QMessageBox::about(this, "Error!", "No images selected\n");
}

/* Show Average of the represented images together ------------------------- */
void ShowSOM::showRepresentedAverageTogether()
{
    unsigned int numMarked = 0;

    // Count how many are marked first
    for (int i = 0; i < listSize; i++)
    {
        if (cellMarks[i])
            numMarked++;
    }
    if (numMarked == 0)
    {
        QMessageBox::about(this, "Error!", "No images selected\n");
        return;
    }

    //vector <string> cell_labels;
    QPixmap *pixmap = new QPixmap[numMarked];
    if (!pixmap)
    {
        QMessageBox::about(this, "Error!", "Out of memory creating new pixmap!\n");
        return;
    }

    // Create a blank Selfile for the averages images.
    SelFile SFAvgs;

    // Go back through cells and add images to list
    for (int i = 0; i < listSize; i++)
    {
        if (cellMarks[i])
        {
            Image _ave;
            if (SFcv[i].ImgNo())
            {
                Image _sd;
                double _minPixel, _maxPixel;
                SelFile SFNew(SFcv[i]);
                SFNew.go_beginning();
                SFNew.get_statistics(_ave, _sd, _minPixel, _maxPixel, apply_geo);
            }
            else
                _ave = Image(projXdim, projYdim);

            // Write _ave to a temp file
            // Need to convert it to ImageXmipp because Selfile can't handle plain
            // Image files when it is being read back in.
            ImageXmipp xm_ave;
            xm_ave = _ave;
            int tempfd;
            string tmpImgfile = makeTempFile(tempfd);
            // Add that image file to the SelFile and save it
            xm_ave.write(tmpImgfile);
            SFAvgs.insert(tmpImgfile);
            ::close(tempfd);
            /*
            int r, c;
            IndextoPos (i, r, c);
            char buf[20];
            sprintf (buf, "row %d col %d", r, c);
            cell_labels.push_back (buf);
            */
        }
    }

    if (SFAvgs.ImgNo())
    {
        ShowSel *showsel = new ShowSel;
        showsel->initWithObject(10, 10, SFAvgs, "Averages of assigned images");
        //for (u = 0; u < numMarked; u++)
        //  v->setCellLabel (u, hisLabels[u], cell_labels[u].c_str());
        showsel->show();
    }
    else
        QMessageBox::about(this, "Error!", "No images selected\n");
}

/* Show assigned sel ------------------------------------------------------- */
void ShowSOM::showRepresentedSel()
{
    SelFile SFNew;
    extractRepresented(SFNew);
    if (SFNew.ImgNo())
    {
        ShowSel *showsel = new ShowSel;
        showsel->initWithObject(10, 10, SFNew, "Represented images");
        showsel->apply_geo = apply_geo;
        showsel->show();
    }
    else QMessageBox::about(this, "Error!", "No images selected\n");
}

/* Show error image -------------------------------------------------------- */
void ShowSOM::showErrorImage()
{
    int row = currentRow();
    int col = currentColumn();
    if (row < 0 || col < 0) return;
    int i = indexOf(row, col);

    SelFile SFNew;
    SFNew.merge(SFcv[i]);
    if (SFNew.ImgNo())
    {
        // Compute the average of the images assigned to that cell
        Image _ave, _sd;
        double _minPixel, _maxPixel;
        SFNew.get_statistics(_ave, _sd, _minPixel, _maxPixel, apply_geo);

        // Load the cell code vector
        Image *image = Image::LoadImage(imgnames[i]), error_image;
        if (!image)
            QMessageBox::about(this, "Error", (QString) "Error loading image " +
                               imgnames[i].c_str() + "!");
        else
        {
            // Compute and show the error
            error_image() = _ave() - (*image)();
            delete image;
            ImageViewer *w = new ImageViewer(&error_image, "Error image");
            w->show();
        }
    }
    else
        QMessageBox::about(this, "Error!", "No images selected\n");
}

/* Show algorithm information ---------------------------------------------- */
void ShowSOM::showAlgoInfo()
{
    QMessageBox::about((QWidget*)this, "Algorithm Information", infStr);
}
