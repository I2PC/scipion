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

#ifndef __QT_DIALOG_FAMILIES_HH__
#define __QT_DIALOG_FAMILIES_HH__

/* Includes ---------------------------------------------------------------- */
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

/* Widget for the micrograph ----------------------------------------------- */
class QtDialogFamilies : public QDialog
{
    Q_OBJECT

private:
    Micrograph  *__m;
    Micrograph  *__mTilted;
    QPushButton *__addFamilyButton;
#ifdef QT3_SUPPORT
    Q3ListBox    *__familyList;
    Q3VBox       *__vBoxLayout;
#else
    QListBox    *__familyList;
    QVBox       *__vBoxLayout;
#endif
    

public:
    // Constructor
    QtDialogFamilies(QWidget *_parent = 0, const char *_name = 0,
                     bool _modal = FALSE, Qt::WFlags _f = 0);
    ~QtDialogFamilies();

    // Set the micrograph
    void setMicrograph(Micrograph *_m);
    void setTiltedMicrograph(Micrograph *_mTilted);

    int findFamily(const char *_familyName);

public slots:
    void slotAddFamily();
    void slotAddFamily(const char *_familyName);
    void slotEditFamily(int _f);
    void slotActiveFamily(int _f);

signals:
    void signalActiveFamily(int _f);
};

#endif
