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
#include "QtWidgetMicrograph.hh"
#include "QtFilterMenu.hh"
#include "QtImageMicrograph.hh"
#include "QtImageOverviewMicrograph.hh"
#include "XmippData/xmippMicrograph.hh"

/* Constructor ------------------------------------------------------------- */
QtWidgetMicrograph::QtWidgetMicrograph( QtMainWidgetMark *_mainWidget, 
                                        QtFiltersController *_f, 
                                        Micrograph *_m ) : 
   QWidget( (QWidget*) _mainWidget ) {
   __filtersController = _f;
   __m              = NULL;
   __activeFamily   = -1;
   __tilted         = FALSE;

   __gridLayout     = new QVBoxLayout( this );
   __menuBar        = new QMenuBar( this );
      __menuBar->setSeparator( QMenuBar::InWindowsStyle );
   __mImage         = new QtImageMicrograph( 0 );
   __mImageOverview = new QtImageOverviewMicrograph( this );
   __file_menu      = NULL;
   
   __mImage->setFiltersController( _f );
   __mImageOverview->setFiltersController( _f );
        
   connect( __mImageOverview, SIGNAL(signalSetCoords(int, int)),
            __mImage, SLOT(slotSetCoords(int, int)) );
   connect( __mImage, SIGNAL(signalSetCoords(int, int)),
            __mImageOverview, SLOT(slotSetCoords(int, int)) );   
   connect( __mImage, SIGNAL(signalSetWidthHeight(int, int)),
            __mImageOverview, SLOT(slotSetWidthHeight(int, int)) );
   connect( __mImage, SIGNAL(signalRepaint( void )),
            __mImageOverview, SLOT(slotRepaint( void )) );
   connect( __mImageOverview, SIGNAL(signalRepaint( void )),
            __mImage, SLOT(slotRepaint( void )) );
   connect( __mImage, SIGNAL(signalRepaint( void )),
            __mImage, SLOT(slotRepaint( void )) );
   connect( __mImageOverview, SIGNAL(signalRepaint( void )),
            __mImageOverview, SLOT(slotRepaint( void )) );
   connect( __mImage, SIGNAL(signalAddCoordOther( int, int, int )),
            this, SLOT(slotDrawEllipse( int,int,int )) );
      
   connect( this, SIGNAL(signalActiveFamily(int)),
            __mImage, SLOT(slotActiveFamily(int)) );
   connect( this, SIGNAL(signalActiveFamily(int)),
            __mImageOverview, SLOT(slotActiveFamily(int)) );
   
   setMicrograph( _m );
   
   __mImage->show();
   __gridLayout->setMenuBar( __menuBar );
   __gridLayout->addWidget( __mImageOverview );
   
   openMenus();   
}

QtWidgetMicrograph::~QtWidgetMicrograph() {
   delete __mImage;
   delete __mImageOverview;
   delete __menuBar;
   delete __gridLayout;
}

/* Set Micrograph ---------------------------------------------------------- */
void QtWidgetMicrograph::setMicrograph( Micrograph *_m ) {
   if ( _m != NULL ) {
      __m = _m;
      __mImage->setMicrograph( _m );
      __mImageOverview->setMicrograph( _m );
   }   
}

/* Open menus -------------------------------------------------------------- */
void QtWidgetMicrograph::openMenus() {
   __file_menu = new QtFileMenu( this );
   connect( __mImage, SIGNAL(signalAddCoordOther(int,int,int)),
            __file_menu, SLOT(slotCoordChange()) );   

   QtFilterMenu *filterMenu = new QtFilterMenu( this );
   
   addMenuItem( "&File", (QtPopupMenuMark *)(__file_menu) );
   addMenuItem( "F&ilters", (QtPopupMenuMark *)(filterMenu) );
   
   connect( __file_menu, SIGNAL(signalAddFamily(const char *)),
            this, SLOT(slotAddFamily(const char*)) );
   
   connect( (QObject*)filterMenu, SIGNAL(signalAddFilter()),
            (QObject*)__filtersController, SLOT(slotAddFilter()) );   
   connect( (QObject*)filterMenu, SIGNAL(signalCleanFilters()),
            (QObject*)__filtersController, SLOT(slotCleanFilters()) );
   connect( (QObject*)filterMenu, SIGNAL(signalCleanFilters()),
            this, SLOT(slotRepaint()) );

   // *** Add your own menus
}

void QtWidgetMicrograph::repaint( int t ) {
   __mImage->repaint( FALSE );
   __mImageOverview->repaint( FALSE );
}

void QtWidgetMicrograph::slotDrawEllipse(int _x, int _y, int _f) {
   __mImage->drawEllipse(_x,_y,_f);
   __mImageOverview->drawEllipse(_x,_y,_f);
}

void QtWidgetMicrograph::slotDrawLastEllipse(int _x, int _y, int _f) {
   __mImage->drawEllipse(_x,_y,_f);
   __mImageOverview->drawEllipse(_x,_y,_f);
   __mImage->drawLastEllipse(_x,_y,_f);
}

/* Active family ----------------------------------------------------------- */
void QtWidgetMicrograph::slotActiveFamily( int _f ) {
   __activeFamily = _f;
   emit signalActiveFamily( _f );
}

void QtWidgetMicrograph::slotAddFamily( const char *_familyName ) {
   emit signalAddFamily( _familyName );
}

void QtWidgetMicrograph::slotDeleteMarkOther( int _coord ) {
   __m->coord(_coord).valid = false;
   repaint();
}

void QtWidgetMicrograph::slotChangeFamilyOther( int _coord, int _f ) {
   __m->coord(_coord).label = _f;
   repaint();
}
