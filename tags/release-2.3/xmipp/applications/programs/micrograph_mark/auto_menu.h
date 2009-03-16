/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 *
 *****************************************************************************/

#ifndef __QT_AUTO_MENU_HH__
#define __QT_AUTO_MENU_HH__

#include "popup_menu_mark.h"

/* Forward declarations ---------------------------------------------------- */
class QtWidgetMicrograph;

/* File menu for Mark ------------------------------------------------------ */
class QtAutoMenu : public QtPopupMenuMark
{
    // For accepting signals and slots
    Q_OBJECT
public:
    // Constructor
    QtAutoMenu(QtWidgetMicrograph* _parent);


public slots:
    // Automatically select particles
    void slotAutoSelectParticles();

    // Save models of particles
    void slotSaveModels();

    // Learn and save
    void slotLearnSaveQuit();
signals:
    void signalAddFamily(const char *);
};

#endif
