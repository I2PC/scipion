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

#include "filters_controller.h"
#include "filter.h"
#include "image_converter.h"

#include <data/image.h>

#include <qdialog.h>
#include <qlistbox.h>
#include <qmessagebox.h>
#include <qimage.h>

/* Constructor ------------------------------------------------------------- */
QtFiltersController::QtFiltersController(QWidget *_parent,
        const Micrograph *_M) :
        QWidget(_parent)
{
    __M = _M;

    __addFilterDialog = new QDialog(this, 0, TRUE);
    __addFilterDialog->setCaption("Add filter");
    __listFilters     = new QListBox(__addFilterDialog);
    __listFilters->setMinimumSize(300, 300);

    __listFilters->insertItem(QtInvertContrastFilter::name);
    __listFilters->insertItem(QtEnhanceContrastFilter::name);
    __listFilters->insertItem(QtSubstractBackgroundFilter::name);
    __listFilters->insertItem(QtRemoveOutlierFilter::name);
    __listFilters->insertItem(QtLowPassFilter::name);
    __listFilters->insertItem(QtHighPassFilter::name);

    connect(__listFilters, SIGNAL(selected(int)),
            this, SLOT(slotAddFilter(int)));
    connect(__listFilters, SIGNAL(selected(int)),
            __addFilterDialog, SLOT(done(int)));
}

QtFiltersController::~QtFiltersController()
{
    for (int i = 0; i < __filterList.size(); i++) delete __filterList[i];

    delete __listFilters;
    delete __addFilterDialog;
}

/* Apply the filters list -------------------------------------------------- */
void QtFiltersController::applyFilters(QImage *_img)
{
    if (__filterList.size() == 0) return;

    QtImageConverter  converter;
    Image            *xmippImg;
    QImage           *qtImg;

    xmippImg = converter.qt2xmipp(_img);

    for (int i = 0; i < __filterList.size(); i++)
        __filterList[i]->apply(xmippImg);

    qtImg    = converter.xmipp2qt(xmippImg);

    delete xmippImg;
    *_img = *qtImg;
    delete qtImg;
}


void QtFiltersController::slotAddFilter()
{
    __addFilterDialog->exec();
}

void QtFiltersController::slotAddFilter(int _f)
{
    switch (_f)
    {
    case invertContrastFilter:
        __filterList.push_back((QtFilter*)new QtInvertContrastFilter(__M));
        break;
    case enhanceContrastFilter:
        __filterList.push_back((QtFilter*)new QtEnhanceContrastFilter(__M));
        break;
    case substractBackgroundFilter:
        __filterList.push_back((QtFilter*)new QtSubstractBackgroundFilter(__M));
        break;
    case removeOutlierFilter:
        __filterList.push_back((QtFilter*)new QtRemoveOutlierFilter(__M));
        break;
    case lowpassFilter:
        __filterList.push_back((QtFilter*)new QtLowPassFilter(__M));
        break;
    case highpassFilter:
        __filterList.push_back((QtFilter*)new QtHighPassFilter(__M));
        break;
    default:
        QMessageBox::information(this, "Error",
                                 "Unable to use this filter");
    }
}

void QtFiltersController::slotCleanFilters()
{
    for (int i = 0; i < __filterList.size(); i++) delete __filterList[i];

    __filterList.clear();
}
