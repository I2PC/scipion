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

#ifndef __QT_WIDGET_MICROGRAPH_HH__
#define __QT_WIDGET_MICROGRAPH_HH__

/* Includes ---------------------------------------------------------------- */
#include "qwidget.h"
#include "qpainter.h"
#include "qlayout.h"
#include "qmenubar.h"
#include "QtImageMicrograph.hh"
#include "QtImageOverviewMicrograph.hh"
#include "QtFileMenu.hh"

/* Forward declarations ---------------------------------------------------- */ 
class QtMainWidgetMark;
class QtImageMicrograph;
class QtPopupMenuMark;
class QtFiltersController;
class Micrograph;

/* Widget for the micrograph ----------------------------------------------- */
class QtWidgetMicrograph : public QWidget {   
   Q_OBJECT
   
private:
   Micrograph                *__m;
   QtFiltersController       *__filtersController;
   int                        __activeFamily;   
   QMenuBar                  *__menuBar;      
   QtImageMicrograph         *__mImage;
   QtImageOverviewMicrograph *__mImageOverview;
   QVBoxLayout               *__gridLayout;
   QtFileMenu                *__file_menu;
   bool                       __tilted;

public:
   // Constructor
   QtWidgetMicrograph( QtMainWidgetMark *_mainWidget, 
                       QtFiltersController *_f,
                       Micrograph *_m = NULL );
   ~QtWidgetMicrograph();
   
   // Set Micrograph
   void setMicrograph( Micrograph *_m );
   
   // Get Micrograph
   Micrograph *getMicrograph() { return( __m ); }
      
   // Set this as tilted micrograph
   void setTilted() {__tilted=TRUE; __mImage->setTilted();}

   // Is tilted?
   bool isTilted() {return __tilted;}

   // Get filters controller
   QtFiltersController *getFiltersController() { return(__filtersController); }
   
   // Get active family
   int activeFamily() { return( __activeFamily); }
   
   // Get overview
   QtImageOverviewMicrograph *overview() { return( __mImageOverview ); }
   
   // Get Image
   QtImageMicrograph *image() { return( __mImage ); }   
   
   // Add menu item
   void addMenuItem( const char *_msg, const QtPopupMenuMark *_item ) {
       __menuBar->insertItem( _msg, (QPopupMenu*)_item );
   }
   
   // Draw axis
   void draw_axis(double _ang)
      {__mImageOverview->enableAxis(); __mImageOverview->draw_axis(_ang);}
   
   // Open menu.
   // Add your menus to this function
   void openMenus();
   
   void repaint( int t=FALSE );
   
public slots:
   void slotActiveFamily( int _f );
   void slotAddFamily( const char *_familyName );
   void slotDeleteMarkOther( int _coord );
   void slotChangeFamilyOther( int _coord, int _f );
   void slotRepaint() { repaint( FALSE ); }
   void slotDrawEllipse(int _x, int _y, int _f);
   void slotDrawLastEllipse(int _x, int _y, int _f);
signals:
   void signalActiveFamily( int _f );
   void signalAddFamily( const char *_familyName );
};


#endif
