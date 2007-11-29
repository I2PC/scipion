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

#include "dialog_families.h"

#include <data/micrograph.h>

#include <qlistbox.h>
#include <qpushbutton.h>
#include <qvbox.h>
#include <qhbox.h>
#include <qlabel.h>
#include <qlineedit.h>

/* Constructor ------------------------------------------------------------- */
QtDialogFamilies::QtDialogFamilies(QWidget *_parent,
                                   const char *_name,
                                   bool _modal,
                                   WFlags _f) :
        QDialog(_parent, _name, _modal, _f)
{

    __m               = NULL;
    __mTilted         = NULL;
    __vBoxLayout      = new QVBox(this);
    __vBoxLayout->setMinimumSize(300, 300);

    __addFamilyButton = new QPushButton("Add family", __vBoxLayout);
    __familyList      = new QListBox(__vBoxLayout);

    connect(__addFamilyButton, SIGNAL(clicked(void)),
            this, SLOT(slotAddFamily(void)));
    connect(__familyList, SIGNAL(selected(int)),
            this, SLOT(slotEditFamily(int)));
    connect(__familyList, SIGNAL(highlighted(int)),
            this, SLOT(slotActiveFamily(int)));
}

QtDialogFamilies::~QtDialogFamilies()
{
    delete __familyList;
    delete __addFamilyButton;
    delete __vBoxLayout;
}

/* Set the micrograph ------------------------------------------------------ */
void QtDialogFamilies::setMicrograph(Micrograph *_m)
{
    if (_m == NULL) return;

    __m = _m;
    __familyList->clear();

    if (__m->LabelNo() == 0) __m->add_label("Common");
    for (int i = 0; i < __m->LabelNo(); i++)
        __familyList->insertItem(__m->get_label(i).c_str());

    __familyList->setSelected(__familyList->count() - 1, TRUE);
    emit signalActiveFamily(__familyList->count() - 1);
}

/* Set the tilted micrograph ----------------------------------------------- */
void QtDialogFamilies::setTiltedMicrograph(Micrograph *_mTilted)
{
    if (_mTilted == NULL) return;

    __mTilted = _mTilted;

    for (int i = 0; i < __m->LabelNo(); i++)
        __mTilted->add_label(__m->get_label(i));
}

/* Find family ------------------------------------------------------------- */
int QtDialogFamilies::findFamily(const char *_familyName)
{
    for (int i = 0; i < __m->LabelNo(); i++)
        if (!strcmp(__m->get_label(i).c_str(), _familyName)) return(i);

    return(-1);
}

/* Add family -------------------------------------------------------------- */
void QtDialogFamilies::slotAddFamily()
{
    if (__m == NULL) return;

    QDialog    familyNameDialog(this, 0, TRUE);
    QHBox      qhbox(&familyNameDialog);
    qhbox.setMinimumSize(200, 20);
    QLabel     familyNameLabel("Family name: ", &qhbox);
    QLineEdit  familyNameLineEdit(&qhbox);

    familyNameLineEdit.setFocus();
    connect(&familyNameLineEdit, SIGNAL(returnPressed(void)),
            &familyNameDialog, SLOT(accept(void)));

    if (familyNameDialog.exec())
    {
        if (!familyNameLineEdit.text().isEmpty())
            slotAddFamily(familyNameLineEdit.text().ascii());
    }
}

/* Add family -------------------------------------------------------------- */
void QtDialogFamilies::slotAddFamily(const char *_familyName)
{
    if (findFamily(_familyName) != -1) return;

    __m->add_label(_familyName);
    if (__mTilted != NULL) __mTilted->add_label(_familyName);
    __familyList->insertItem(_familyName);
    __familyList->setSelected(__familyList->count() - 1, TRUE);
    emit signalActiveFamily(__familyList->count() - 1);
}

/* Edit family ------------------------------------------------------------- */
void QtDialogFamilies::slotEditFamily(int _f)
{
    if (__m == NULL) return;

    QDialog    familyNameDialog(this, 0, TRUE);
    QHBox      qhbox(&familyNameDialog);
    qhbox.setMinimumSize(200, 20);
    QLabel     familyNameLabel("Family name: ", &qhbox);
    QLineEdit  familyNameLineEdit(&qhbox);
    familyNameLineEdit.setText(__m->get_label(_f).c_str());

    familyNameLineEdit.setFocus();
    connect(&familyNameLineEdit, SIGNAL(returnPressed(void)),
            &familyNameDialog, SLOT(accept(void)));

    if (familyNameDialog.exec())
    {
        if (!familyNameLineEdit.text().isEmpty())
        {
            __m->get_label(_f) = familyNameLineEdit.text().ascii();
            if (__mTilted != NULL)
                __mTilted->get_label(_f) = familyNameLineEdit.text().ascii();
            __familyList->changeItem(familyNameLineEdit.text().ascii(), _f);
        }
    }
}

/* Active Family ----------------------------------------------------------- */
void QtDialogFamilies::slotActiveFamily(int _f)
{
    emit signalActiveFamily(_f);
}
