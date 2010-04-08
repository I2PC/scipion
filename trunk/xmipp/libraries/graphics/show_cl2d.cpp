/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include "show_cl2d.h"
#include "show_2d.h"

#include <qmessagebox.h>

#ifdef QT3_SUPPORT
#include <q3filedialog.h>
//Added by qt3to4:
#include <QPixmap>
#include <Q3PopupMenu>
#else
#include <qfiledialog.h>
#endif

/* Init/Clear data --------------------------------------------------------- */
void ShowCL2D::init()
{
    SFcv        = NULL;
    hisAssigned = NULL;
    ShowSel::init();
}

void ShowCL2D::clear()
{
    if (SFcv        != NULL) delete [] SFcv;
    if (hisAssigned != NULL) delete [] hisAssigned;
    ShowSel::clear();
}

/* Initialize with a SOM file.---------------------------------------------- */
void ShowCL2D::initWithFile(const FileName &_fn_root)
{
    init();
    readFile(_fn_root);
    initTable();
    initRightclickMenubar();
    repaint();
}

/* Read a CL2D ------------------------------------------------------------- */
void ShowCL2D::readFile(const FileName &_fn_root,
                       double _minGray, double _maxGray)
{
    clear();
    fn = _fn_root;
    setCaption(fn.c_str());
    readSelFile(_fn_root + ".sel");
    SFcv=new SelFile[listSize];
    hisAssigned=new std::string[listSize];
    for (int i=0; i<listSize; i++)
    {
        FileName fnCV=imgnames[i].without_extension();
        SFcv[i].read(fnCV+".sel");
        hisAssigned[i]=integerToString(SFcv[i].ImgNo());
    }
    NumCols=FLOOR(sqrt(listSize));
    NumRows=CEIL(sqrt(listSize));
}

/* Initialize right click menubar ------------------------------------------ */
void ShowCL2D::initRightclickMenubar()
{
#ifdef QT3_SUPPORT
    menubar = new Q3PopupMenu();
    Q3PopupMenu * file = new Q3PopupMenu();
#else
    menubar = new QPopupMenu();
    QPopupMenu * file = new QPopupMenu();
#endif
    file->insertItem("Open...", this,  SLOT(GUIopenFile()));
    file->insertItem("Save assigned images in a sel file...",
                     this, SLOT(saveAssigned()), Qt::CTRL + Qt::Key_N);
    menubar->insertItem("&File", file);

    // Options .............................................................
#ifdef QT3_SUPPORT
    options =  new Q3PopupMenu();
#else
    options = new QPopupMenu();
#endif
    setCommonOptionsRightclickMenubar();

    // What kind of labels
    mi_nameAsLabels = options->insertItem("Show Image Names as Labels", this,  SLOT(changeLabelsToName()));
    mi_hisAsLabels = options->insertItem("Show Histogram as Labels", this, SLOT(changeLabelsToHis()));
    options->setItemEnabled(mi_nameAsLabels, true);
    options->setItemEnabled(mi_hisAsLabels, false);
    labeltype = Histogram_LABEL;
    options->insertSeparator();

    // Statistics
    options->insertItem("View assigned images",  this,  SLOT(showAssigned()));
    options->insertItem("Show this image separately", this, SLOT(showThisImage()));
    options->insertSeparator();

    // Insert options the menu
    menubar->insertItem("&Options", options);
    menubar->insertSeparator();

    // Insert Help and Quit
    insertGeneralItemsInRightclickMenubar();
}

/* Extract represented .---------------------------------------------------- */
void ShowCL2D::extractRepresented(SelFile &SF_represented)
{
    for (int i = 0; i < listSize; i++)
    {
        if (cellMarks[i])
            SF_represented.merge(SFcv[i]);
    }
}

void ShowCL2D::saveAssigned()
{
    SelFile SFNew;
    extractRepresented(SFNew);
    if (SFNew.ImgNo() != 0) writeSelFile(SFNew);
    else QMessageBox::about(this, "Error!", "No images selected\n");
}

void ShowCL2D::showAssigned()
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

/* Change labels ----------------------------------------------------------- */
void ShowCL2D::changeLabelsToName()
{
    changeLabel(mi_nameAsLabels);
}

void ShowCL2D::changeLabelsToHis()
{
    changeLabel(mi_hisAsLabels);
}

void ShowCL2D::changeLabel(int _clicked_mi)
{
    options->setItemEnabled(mi_nameAsLabels, true);
    options->setItemEnabled(mi_hisAsLabels, true);
    options->setItemEnabled(_clicked_mi,    false);
    if (_clicked_mi == mi_nameAsLabels) labeltype = Name_LABEL;
    else if (_clicked_mi == mi_hisAsLabels) labeltype = Histogram_LABEL;
    repaintContents();
}

const char * ShowCL2D::cellLabel(int i) const
{
    if (options->isItemEnabled(mi_showLabel)) return NULL;
    switch (labeltype)
    {
        case Name_LABEL:
            return imgnames[i].c_str();
        case Histogram_LABEL:
            return hisAssigned[i].c_str();
    }
}

// Show This image ---------------------------------------------------------
void ShowCL2D::showThisImage()
{
    int row = currentRow();
    int col = currentColumn();
    int i = indexOf(row, col);

    ImageViewer *showimg = new ImageViewer(imgnames[i].c_str(), false);
    showimg->loadImage(imgnames[i].c_str());
    showimg->show();
}
