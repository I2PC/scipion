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
#include <qwidget.h>
#include <qpainter.h>
#include <qlayout.h>
#include <qmenubar.h>
#include <qaccel.h>
#include <qscrollbar.h>
#include <qlabel.h>
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
   int                        __mingray;
   int                        __maxgray;
   float                      __gamma;

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
   
   // Get Filemenu
   QtFileMenu *file_menu() {return __file_menu;}
   
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
   
   // Change contrast
   void changeContrast(int _mingray, int _maxgray, float _gamma);
   
   // Repaint 
   void repaint( int t=FALSE );
   
public slots:
   void slotActiveFamily( int _f );
   void slotAddFamily( const char *_familyName );
   void slotDeleteMarkOther( int _coord );
   void slotChangeFamilyOther( int _coord, int _f );
   void slotRepaint() { repaint( FALSE ); }
   void slotDrawEllipse(int _x, int _y, int _f);
   void slotDrawLastEllipse(int _x, int _y, int _f);
   void slotQuit();
   void slotChangeContrast();
signals:
   void signalActiveFamily( int _f );
   void signalAddFamily( const char *_familyName );
};

/** Class to adjust contrast
*/
class AdjustContrastWidget : public QWidget {
   Q_OBJECT
public:
   /** Constructor */
   AdjustContrastWidget(int min, int max, float gamma, 
      QtWidgetMicrograph *_qtwidgetmicrograph,
      QWidget *parent=0, const char *name=0, int wflags=0);
private:
   QtWidgetMicrograph *__qtwidgetmicrograph;
   QScrollBar 	      *__scroll_min;
   QScrollBar 	      *__scroll_max;
   QScrollBar 	      *__scroll_gamma;
   QLabel     	      *__label_min;
   QLabel     	      *__label_max;
   QLabel     	      *__label_gamma;
private slots:
   void scrollValueChanged(int);  
};

#endif
