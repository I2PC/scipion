/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Carlos Manzanares       (cmanzana@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/

/* Includes ---------------------------------------------------------------- */
#include "QtDialogProperties.hh"
#include "QtImageMicrograph.hh"
#include "XmippData/xmippMicrograph.hh"
#include "qlistbox.h"
#include "qpushbutton.h"
#include "qvbox.h"
#include "qhbox.h"

/* Constructor ------------------------------------------------------------- */
QtDialogProperties::QtDialogProperties( Micrograph *_m, int _coord, 
                                        QWidget *_parent, 
                                        const char *_name, 
                                        bool _modal,
                                        WFlags _f ) : 
   QDialog( _parent, _name, _modal, _f ) {

   setCaption( "Change properties" );
   __m            = _m;
   __coord        = _coord;
   __moving       = false;   
   __vBoxLayout   = new QVBox( this );
   __vBoxLayout->setMinimumSize( 300, 300 );
   __moveButton   = new QPushButton( "Move", __vBoxLayout );
   __deleteButton = new QPushButton( "Delete", __vBoxLayout );
   __familyList   = new QListBox( __vBoxLayout );

   for( int i = 0; i < __m->LabelNo(); i++ )
      __familyList->insertItem( __m->get_label(i).c_str() );   
   __familyList->setSelected( __m->coord(__coord).label, TRUE );
   
   connect( __familyList, SIGNAL(highlighted(int)),
            this, SLOT(slotChangeFamily(int)) );
   connect( __deleteButton, SIGNAL(clicked(void)),
            this, SLOT(slotDeleteMark(void)) );
   connect( __moveButton, SIGNAL(clicked(void)), 
            this, SLOT(slotMoveMark(void)) );
}

QtDialogProperties::~QtDialogProperties() {
   delete __familyList;
   delete __deleteButton;
   delete __moveButton;
   delete __vBoxLayout;
}

void QtDialogProperties::slotChangeFamily( int _f ) {
   __m->coord(__coord).label = _f;
   emit signalChangeFamilyOther( __coord, _f );
   accept();
}

void QtDialogProperties::slotDeleteMark() {
   __m->coord(__coord).valid = false;
   emit signalDeleteMarkOther( __coord );
   accept();
}

void QtDialogProperties::slotMoveMark() {
   ((QtImageMicrograph*)parentWidget())->movingMark( __coord );
   accept();
}
