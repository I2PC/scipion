/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Carlos Manzanares       (cmanzana@cnb.uam.es)
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

#ifndef __QT_FILE_MENU_HH__
#define __QT_FILE_MENU_HH__

/* Includes ---------------------------------------------------------------- */
#include "QtPopupMenuMark.hh"
#include <qradiobutton.h>

/* Forward declarations ---------------------------------------------------- */ 
class QtWidgetMicrograph;

/* File menu for Mark ------------------------------------------------------ */
class QtFileMenu : public QtPopupMenuMark {
   // For accepting signals and slots
   Q_OBJECT
   
   // Coordinates have been saved
   bool __coordinates_are_saved;
public:
   // Constructor
   QtFileMenu( QtWidgetMicrograph* _parent );


public slots:
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
   void slotCoordChange() {__coordinates_are_saved=FALSE;}
signals:
   void signalAddFamily( const char * );
};

#endif
