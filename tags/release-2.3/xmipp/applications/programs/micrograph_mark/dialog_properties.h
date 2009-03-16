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

#ifndef __QT_DIALOG_PROPERTIES_HH__
#define __QT_DIALOG_PROPERTIES_HH__

#include <qdialog.h>

/* Forward declarations ---------------------------------------------------- */
#ifdef QT3_SUPPORT
class Q3ListBox;
class Q3VBox;
#else
class QListBox;
class QVBox;
#endif
class QPushButton;
class Micrograph;
class QtWidgetMicrograph;

/* Widget for the micrograph ----------------------------------------------- */
class QtDialogProperties : public QDialog
{
    Q_OBJECT

private:
    Micrograph         *__m;
    QtWidgetMicrograph *__wm;
    int                 __coord;
    bool                __moving;
#ifdef QT3_SUPPORT
    Q3ListBox           *__familyList;
    Q3VBox              *__vBoxLayout;
#else
    QListBox           *__familyList;
    QVBox              *__vBoxLayout;
#endif
    QPushButton        *__moveButton;
    QPushButton        *__deleteButton;

public:
    // Constructor
    QtDialogProperties(Micrograph *_m, QtWidgetMicrograph *_wm,
                       int _coord, QWidget *_parent = 0,
                       const char *_name = 0, bool _modal = FALSE, Qt::WFlags _f = 0);
    ~QtDialogProperties();

public slots:
    void slotChangeFamily(int _f);
    void slotDeleteMark();
    void slotMoveMark();

signals:
    void signalDeleteMarkOther(int _coord);
    void signalChangeFamilyOther(int _coord, int _f);
};

#endif
