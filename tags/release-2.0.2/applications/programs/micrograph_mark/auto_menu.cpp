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

#include <cstring>

#include "auto_menu.h"
#include "widget_micrograph.h"
#include "main_widget_mark.h"

#include <data/funcs.h>

#include <qfiledialog.h>
#include <qgrid.h>
#include <qlabel.h>
#include <qpushbutton.h>
#include <qmessagebox.h>
#include <qlineedit.h>

/* Constructor ------------------------------------------------------------- */
QtAutoMenu::QtAutoMenu(QtWidgetMicrograph* _parent) :
        QtPopupMenuMark(_parent)
{
    insertItem("Configure", this, SLOT(slotConfig()));
    insertItem("Load models", this, SLOT(slotLoadModels()));
    insertItem("Learn particles", this, SLOT(slotLearnParticles()));
    insertItem("AutoSelect", this, SLOT(slotAutoSelectParticles()));
    insertItem("Save models", this, SLOT(slotSaveModels()));
}

/* Configure --------------------------------------------------------------- */
void QtAutoMenu::slotConfig()
{
    ((QtWidgetMicrograph *)parentWidget())->configure_auto();
}

/* Load models of Particles ------------------------------------------------ */
void QtAutoMenu::slotLoadModels()
{
    ((QtWidgetMicrograph *)parentWidget())->loadModels();
}

/* Learn Particles --------------------------------------------------------- */
void QtAutoMenu::slotLearnParticles()
{
    ((QtWidgetMicrograph *)parentWidget())->learnParticles();
}

/* AutoSelect Particles ---------------------------------------------------- */
void QtAutoMenu::slotAutoSelectParticles()
{
    ((QtWidgetMicrograph *)parentWidget())->automaticallySelectParticles();
}

/* Save models of Particles ------------------------------------------------ */
void QtAutoMenu::slotSaveModels()
{
    ((QtWidgetMicrograph *)parentWidget())->saveModels();
}
