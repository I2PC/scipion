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

#include <cstring>

#include "file_menu.h"
#include "widget_micrograph.h"
#include "main_widget_mark.h"

#include <data/micrograph.h>
#include <data/funcs.h>
#include <data/args.h>

#include <qlabel.h>
#include <qpushbutton.h>
#include <qmessagebox.h>
#include <qlineedit.h>

#ifdef QT3_SUPPORT
#include <q3filedialog.h>
#include <q3grid.h>
//Added by qt3to4:
#include <Q3PopupMenu>
#else
#include <qfiledialog.h>
#include <qgrid.h>
#endif

/* Constructor ------------------------------------------------------------- */
QtFileMenu::QtFileMenu(QtWidgetMicrograph* _parent) :
        QtPopupMenuMark(_parent)
{
    __coordinates_are_saved = TRUE;

    insertItem("Load coords", this, SLOT(slotLoadCoords()));
    insertItem("Save coords", this, SLOT(slotSaveCoords()));
    insertItem("Save angles", this, SLOT(slotSaveAngles()));
    insertItem("Generate images", this, SLOT(slotGenerateImages()));

    insertSeparator();
#ifdef QT3_SUPPORT
    options =  new Q3PopupMenu();
#else
    options =  new QPopupMenu();
#endif
    insertItem("Change mark type", options);
    circle = options->insertItem("Circle");
    square = options->insertItem("Square");
    options->setCheckable(TRUE);
    connect(options, SIGNAL(activated(int)), this, SLOT(doOption(int)));
    setMouseTracking(TRUE);
    options->setItemChecked(circle, true);
    insertItem("Change mark radius", this, SLOT(slotChangeCircleRadius()));
    insertItem("Show families", this, SLOT(slotShowFamilies()));
    insertItem("Micrograph Info", this, SLOT(slotMicrographInfo()));
    insertSeparator();

    insertItem("Quit", this, SLOT(slotQuit()));
}


/* Change mark type -------------------------------------------------------- */
void QtFileMenu::doOption(int item)
{

    Micrograph *m = ((QtWidgetMicrograph*)parentWidget())->getMicrograph();
    if (m == NULL) return;

    if (options->isItemChecked(item)) return;     // They are all radio buttons

    if (item == circle)
    {
        options->setItemChecked(circle, true);
        options->setItemChecked(square, false);
	((QtWidgetMicrograph*)parentWidget())->changeMarkType(MARK_CIRCLE);
    }
    else if (item == square)
    {
        options->setItemChecked(circle, false);
        options->setItemChecked(square, true);
	((QtWidgetMicrograph*)parentWidget())->changeMarkType(MARK_SQUARE);
    }
}

/* Change circle radius ---------------------------------------------------- */
void QtFileMenu::slotChangeCircleRadius()
{
    Micrograph *m = ((QtWidgetMicrograph*)parentWidget())->getMicrograph();
    if (m == NULL) return;
    ((QtWidgetMicrograph*)parentWidget())->slotChangeCircleRadius();
}

/* Show families ----------------------------------------------------------- */
void QtFileMenu::slotShowFamilies()
{
    ((QtMainWidgetMark*) ((QtWidgetMicrograph*)parentWidget())->parentWidget())
        ->showFamilies();
}

/* Show Micrograph info ---------------------------------------------------- */
void QtFileMenu::slotMicrographInfo()
{
    Micrograph *m = ((QtWidgetMicrograph*)parentWidget())->getMicrograph();
    if (m == NULL) return;

    int Xdim, Ydim;
    m->size(Xdim,Ydim);
    std::string message = (std::string)"Micrograph name: "+m->micrograph_name();
    message+=(std::string)"\nSize YxX: "+integerToString(Ydim)+"x"+
        integerToString(Xdim);
    QMessageBox helpmsg((QString)"Information", (QString) message.c_str(),
        QMessageBox::Information, QMessageBox::Ok, 0, 0, 0, 0, false);
    helpmsg.exec();
}

/* Load coordinates -------------------------------------------------------- */
void QtFileMenu::slotLoadCoords()
{
    Micrograph *m = ((QtWidgetMicrograph*)parentWidget())->getMicrograph();
    if (m == NULL) return;

    try
    {
#ifdef QT3_SUPPORT
        Q3FileDialog coordFileDialog(this, 0, TRUE);
#else
        QFileDialog coordFileDialog(this, 0, TRUE);
#endif
        coordFileDialog.setCaption("Coordinates filename");
        coordFileDialog.setFilter("*.pos");
        if (coordFileDialog.exec())
        {
            // Get the family name from the filename
            FileName fn;
            fn = (char*)coordFileDialog.selectedFile().ascii();
            fn = fn.without_extension();
            FileName familyName;
            int i=fn.find(m->micrograph_name()+".");
            if (i!=-1)
                familyName=fn.substr(i+m->micrograph_name().length()+1,
                    fn.length()-(i+m->micrograph_name().length()+1));
            else familyName=fn.get_extension().c_str();
            emit signalAddFamily(familyName.c_str());

            int activeFamily = ((QtWidgetMicrograph*)parentWidget())->activeFamily();
            m->read_coordinates(activeFamily,
                                (char*)coordFileDialog.selectedFile().ascii());
            ((QtWidgetMicrograph*)parentWidget())->repaint();
            slotCoordChange();
        }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
}

/* Save coordinates -------------------------------------------------------- */
void QtFileMenu::slotSaveCoords()
{
    Micrograph *m = ((QtWidgetMicrograph*)parentWidget())->getMicrograph();
    if (m == NULL) return;

    switch (QMessageBox::information(this, "Mark",
                                     "Save active family coordinates\n"
                                     "Are you sure?",
                                     "&Yes", "&No",
                                     0,      // Enter == button 0
                                     1))
    { // Escape == button 1
    case 0: // Yes clicked, Alt-Y or Enter pressed.
        break;
    case 1: // No clicked or Alt-N pressed or ESC
        return;
        break;
    }

    int activeFamily = ((QtWidgetMicrograph*)parentWidget())->activeFamily();
    m->write_coordinates(activeFamily,  m->micrograph_name() + "." +
                         m->get_label(activeFamily) +
                         ".pos");
    __coordinates_are_saved = TRUE;
}

/* Generate images --------------------------------------------------------- */
void QtFileMenu::slotGenerateImages()
{
    Micrograph *m = ((QtWidgetMicrograph*)parentWidget())->getMicrograph();
    if (m == NULL) return;
    /*
       switch(QMessageBox::information( this, "Mark",
                                       "Generate active family images\n"
                                       "Are you sure?",
                                       "&Yes", "&No",
                                       0,      // Enter == button 0
                                       1 ) ) { // Escape == button 1
           case 0: // Yes clicked, Alt-Y or Enter pressed.
               break;
           case 1: // No clicked or Alt-N pressed or ESC
               return; break;
        }
    */
    int activeFamily = ((QtWidgetMicrograph*)parentWidget())->activeFamily();
    try
    {
        QDialog   setPropertiesDialog(this, 0, TRUE);
        setPropertiesDialog.setCaption("Generate images");
#ifdef QT3_SUPPORT
        Q3Grid     qgrid(2, &setPropertiesDialog);
#else
        QGrid     qgrid(2, &setPropertiesDialog);
#endif
        qgrid.setMinimumSize(250, 170);
        QLabel    windowSizeXLabel("X window size: ", &qgrid);
        QLineEdit windowSizeXLineEdit(&qgrid);
        windowSizeXLineEdit.setText("100");
        QLabel    windowSizeYLabel("Y window size: ", &qgrid);
        QLineEdit windowSizeYLineEdit(&qgrid);
        windowSizeYLineEdit.setText("100");
        QLabel    rootNameLabel("Root name: ", &qgrid);
        QLineEdit rootNameLineEdit(&qgrid);
        QLabel    startingIndex("Starting index: ", &qgrid);
        QLineEdit startingIndexLineEdit(&qgrid);
        startingIndexLineEdit.setText("1");
        QLabel    originalM("Original micrograph: ", &qgrid);
        QLineEdit originalMLineEdit(&qgrid);
	// Sjors & Roberto: 25jan08
	// Here we need the name of the original micrograph, as we
	// might have a 8bit version in memory because this gives
	// nicer visualization with the XVsmooth routines
	// Also see our explanation in main.cpp
	FileName fnt = (((QtWidgetMicrograph*)parentWidget())->getMicrograph())->micrograph_name();
        originalMLineEdit.setText(fnt.c_str());
        QRadioButton  computeTransmitance("Transmitance -> O.D.  ", &qgrid);
        QRadioButton  computeInverse("Invert", &qgrid);
        computeTransmitance.setChecked(FALSE);
        computeInverse.setChecked(FALSE);
        QPushButton okButton("Ok", &qgrid);
        QPushButton cancelButton("Cancel", &qgrid);

        connect(&okButton, SIGNAL(clicked(void)),
                &setPropertiesDialog, SLOT(accept(void)));
        connect(&cancelButton, SIGNAL(clicked(void)),
                &setPropertiesDialog, SLOT(reject(void)));

        if (setPropertiesDialog.exec())
        {
            m->set_window_size(windowSizeXLineEdit.text().toInt(),
                               windowSizeYLineEdit.text().toInt());
            m->set_transmitance_flag(computeTransmitance.isChecked());
            m->set_inverse_flag(computeInverse.isChecked());
            if (!rootNameLineEdit.text().isEmpty())
            {
                // Select right angle
#define MAIN_WIDGET \
    ((QtMainWidgetMark *)(parentWidget()->parentWidget()))
                double alpha = 0, psi = 0, gamma = 0;
                bool this_is_tilted;
                if (MAIN_WIDGET->there_is_tilted())
                {
                    if (((QtWidgetMicrograph *)(MAIN_WIDGET->tilted_widget())) ==
                        (QtWidgetMicrograph *) parentWidget())
                    {
                        psi   = MAIN_WIDGET->alpha_t();
                        gamma = MAIN_WIDGET->gamma_t();
                        this_is_tilted = true;
                    }
                    else
                    {
                        gamma = 0;
                        alpha = MAIN_WIDGET->alpha_u();
                        this_is_tilted = false;
                    }
                }
                m->produce_all_images(activeFamily,
                                      (char*)rootNameLineEdit.text().ascii(),
                                      startingIndexLineEdit.text().toInt(),
                                      (char*)originalMLineEdit.text().ascii(),
                                      alpha, gamma, psi/*rot,tilt,psi*/);
                MAIN_WIDGET->generated(this_is_tilted,
                                       m->get_label(activeFamily));
            }
        }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
}

/* Save angles ------------------------------------------------------------- */
void QtFileMenu::slotSaveAngles()
{
    switch (QMessageBox::information(this, "Mark",
                                     "Write micrograph angles\n"
                                     "Are you sure?",
                                     "&Yes", "&No",
                                     0,      // Enter == button 0
                                     1))
    { // Escape == button 1
    case 0: // Yes clicked, Alt-Y or Enter pressed.
        break;
    case 1: // No clicked or Alt-N pressed or ESC
        return;
        break;
    }

    ((QtMainWidgetMark *)(parentWidget()->parentWidget()))->write_angles();
}

/* Quit -------------------------------------------------------------------- */
void QtFileMenu::slotQuit()
{
    QtMainWidgetMark *M = (QtMainWidgetMark *) parentWidget()->parentWidget();
    // If there is no tilted micrograph there is nothing to do
    if (M->tilted_widget() != NULL)
    {
        M->compute_gamma();
        M->compute_alphas();
        slotSaveAngles();
    }

    if (!__coordinates_are_saved)
    {
        switch (QMessageBox::information(this, "Mark",
                                         "The document contains unsaved work\n"
                                         "Do you want to save it before exiting?",
                                         "&Save", "&Don't Save", "&Cancel",
                                         0,      // Enter == button 0
                                         2))
        { // Escape == button 2
        case 0: // Save clicked, Alt-S or Enter pressed.
            slotSaveCoords();
            exit(0);
            break;
        case 1: // Don't Save clicked or Alt-D pressed
            exit(1);
            break;
        case 2: // Cancel clicked, Alt-C or Escape pressed
            // don't do anything
            break;
        }
    }
    else exit(1);
}
