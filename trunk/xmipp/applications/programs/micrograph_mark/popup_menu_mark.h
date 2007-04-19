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

#ifndef __QT_POPUP_MENU_MARK_HH__
#define __QT_POPUP_MENU_MARK_HH__

#include <qpopupmenu.h>

/* Forward declarations ---------------------------------------------------- */
class QtWidgetMicrograph;

/* Virtual class for Mark menus -------------------------------------------- */
class QtPopupMenuMark : public QPopupMenu {
   Q_OBJECT

public:
   // Constructor. Set reference to main_widget
   QtPopupMenuMark( QtWidgetMicrograph *_parent ) :
      QPopupMenu( (QWidget*)_parent ) {}

public slots:
};

#endif
