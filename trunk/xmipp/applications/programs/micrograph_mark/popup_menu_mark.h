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

#ifndef __QT_POPUP_MENU_MARK_HH__
#define __QT_POPUP_MENU_MARK_HH__

#ifdef QT3_SUPPORT
// MOC_SKIP_BEGIN
#include <q3popupmenu.h>
// MOC_SKIP_END
#else
#include <qpopupmenu.h>
#endif

/* Forward declarations ---------------------------------------------------- */
class QtWidgetMicrograph;

/* Virtual class for Mark menus -------------------------------------------- */
#ifdef QT3_SUPPORT
// MOC_SKIP_BEGIN
class QtPopupMenuMark : public Q3PopupMenu
// MOC_SKIP_END
#else
class QtPopupMenuMark : public QPopupMenu
#endif
{
    Q_OBJECT

public:
    // Constructor. Set reference to main_widget
    QtPopupMenuMark(QtWidgetMicrograph *_parent) :
#ifdef QT3_SUPPORT
        Q3PopupMenu((QWidget*)_parent)
#else
        QPopupMenu((QWidget*) _parent)
#endif
    {}

public slots:
};

#endif
