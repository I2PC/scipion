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

/* Includes ---------------------------------------------------------------- */
#include "qdialog.h"

/* Forward declarations ---------------------------------------------------- */
class QListBox;
class QPushButton;
class QVBox;
class Micrograph;

/* Widget for the micrograph ----------------------------------------------- */
class QtDialogProperties : public QDialog {   
   Q_OBJECT

private:
   Micrograph  *__m;   
   int          __coord;
   bool         __moving;
   QListBox    *__familyList;
   QPushButton *__moveButton;
   QPushButton *__deleteButton;
   QVBox       *__vBoxLayout;

public:
   // Constructor
   QtDialogProperties( Micrograph *_m, int _coord, QWidget *_parent=0, 
                       const char *_name=0, bool _modal=FALSE, WFlags _f=0 );
   ~QtDialogProperties();
         
public slots:
   void slotChangeFamily( int _f );
   void slotDeleteMark();
   void slotMoveMark();
   
signals:
   void signalDeleteMarkOther( int _coord );
   void signalChangeFamilyOther( int _coord, int _f );
};

#endif
