/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Carlos Manzanares       (cmanzana@cnb.csic.es)
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

#include "filter_menu.h"
#include "widget_micrograph.h"

/* Constructor ------------------------------------------------------------- */
QtFilterMenu::QtFilterMenu(QtWidgetMicrograph* _parent) :
        QtPopupMenuMark(_parent)
{
    insertItem("Adjust contrast", this, SLOT(slotAdjustContrast()));
    insertItem("Crop micrograph", this, SLOT(slotCrop()));
    insertItem("Add filter", this, SLOT(slotAddFilter()));
    insertItem("Clean filters", this, SLOT(slotCleanFilters()));
}

/* Add Filter -------------------------------------------------------------- */
void QtFilterMenu::slotAdjustContrast()
{
    emit signalAdjustContrast();
}

void QtFilterMenu::slotCrop()
{
    emit signalCrop();
}

void QtFilterMenu::slotAddFilter()
{
    emit signalAddFilter();
}

void QtFilterMenu::slotCleanFilters()
{
    emit signalCleanFilters();
}
