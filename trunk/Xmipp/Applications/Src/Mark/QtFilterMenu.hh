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

#ifndef __QT_FILTER_MENU_HH__
#define __QT_FILTER_MENU_HH__

/* Includes ---------------------------------------------------------------- */
#include "QtPopupMenuMark.hh"

/* Forward declarations ---------------------------------------------------- */ 
class QtWidgetMicrograph;

/* Filter menu for Mark ---------------------------------------------------- */
class QtFilterMenu : public QtPopupMenuMark {
   // For accepting signals and slots
   Q_OBJECT
   
public:
   // Constructor
   QtFilterMenu( QtWidgetMicrograph* _parent );


public slots:
   // Add filter
   void slotAdjustContrast();
   void slotAddFilter();
   void slotCleanFilters();

signals:
   void signalAdjustContrast();
   void signalAddFilter();
   void signalCleanFilters();
};

#endif
