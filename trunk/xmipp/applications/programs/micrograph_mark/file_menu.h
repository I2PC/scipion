/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Carlos Manzanares       (cmanzana@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Copyright (c) 2000 , CSIC .
 *
 * Permission is granted to copy and distribute this file, for noncommercial
 * use, provided (a) this copyright notice is preserved, (b) no attempt
 * is made to restrict redistribution of this file, and (c) this file is
 * restricted by a compilation copyright.
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 *
 *****************************************************************************/

#ifndef __QT_FILE_MENU_HH__
#define __QT_FILE_MENU_HH__

#define MARK_CIRCLE 0
#define MARK_SQUARE 1

#include "popup_menu_mark.h"

#include <qradiobutton.h>

#ifdef QT3_SUPPORT
//Added by qt3to4:
#include <Q3PopupMenu>
#endif

/* Forward declarations ---------------------------------------------------- */
class QtWidgetMicrograph;

/* File menu for Mark ------------------------------------------------------ */
class QtFileMenu : public QtPopupMenuMark
{
    // For accepting signals and slots
    Q_OBJECT

    // For changing mark type
#ifdef QT3_SUPPORT
    Q3PopupMenu *options;
#else
    QPopupMenu* options;
#endif
    int circle, square;

    // Coordinates have been saved
    bool __coordinates_are_saved;
public:
    // Constructor
    QtFileMenu(QtWidgetMicrograph* _parent);

public slots:

    // Change mark type
    void doOption(int item);

    // Change circle radius
    void slotChangeCircleRadius();

    // Show families
    void slotShowFamilies();

    // Show micrograph info
    void slotMicrographInfo();

    // Load coords
    void slotLoadCoords();

    // Save coords
    void slotSaveCoords();

    // Generate images
    void slotGenerateImages();

    // Save angles
    void slotSaveAngles();

    // Quit
    void slotQuit();

    // Coordinates changed
    void slotCoordChange()
    {
        __coordinates_are_saved = FALSE;
    }
signals:
    void signalAddFamily(const char *);
};

#endif
