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

#include "image.h"
#include "main_widget_mark.h"
#include "filters_controller.h"

#include <data/micrograph.h>

#include <qpainter.h>
#include <qimage.h>

double QtImage::__zoom = 1.0;

/* Coordiantes transformations --------------------------------------------- */
void QtImage::micrographToImage( int _x, int _y, int &_rx, int &_ry ) {
   _rx = ROUND(_x / __zoom) - __x;
   _ry = ROUND(_y / __zoom) - __y;
}

void QtImage::imageToMicrograph( int _x, int _y, int &_rx, int &_ry ) {
   _rx = ROUND((_x + __x) * __zoom);
   _ry = ROUND((_y + __y) * __zoom);
}

void QtImage::exact_micrographToImage( int _x, int _y, double &_rx, double &_ry ) {
   _rx = (double)_x / __zoom - __x;
   _ry = (double)_y / __zoom - __y;
}

void QtImage::exact_imageToMicrograph( int _x, int _y, double &_rx, double &_ry ) {
   _rx = ((double)_x + __x) * __zoom;
   _ry = ((double)_y + __y) * __zoom;
}

/* Constructor ------------------------------------------------------------- */
QtImage::QtImage( QWidget *_parent, const char *_name, WFlags _f ) :
   QWidget( _parent, _name, _f ) {
   __m                 = NULL;
   __filtersController = NULL;
   __img               = new QImage( size(), 8, 256 );
   __activeFamily      = -1;
   __paint             = new QPainter( this );

   __mingray=0;
   __maxgray=255;
   __gamma=1;

   for( int i = 0; i < 256; i++ ) __img->setColor( i, qRgb( i, i, i ) );
}

QtImage::~QtImage() {
   delete __img;
   delete __paint;
}

/* Set Micrograph ---------------------------------------------------------- */
void QtImage::setMicrograph( Micrograph *_m ) {
   if ( _m != NULL ) {
      __m = _m;
      setCaption(_m->micrograph_name().c_str());
   }
}

/* Set WidgetMicrograph ---------------------------------------------------- */
void QtImage::setWidgetMicrograph( QtWidgetMicrograph *_wm ) {
   __wm=_wm;
}

/* Set Pixel --------------------------------------------------------------- */
void QtImage::setPixel(int _x, int _y, int _value) {
    __img->setPixel( _x, _y, _value);
}

/* Set the filters controller ---------------------------------------------- */
void QtImage::setFiltersController( QtFiltersController *_f ) {
   __filtersController = _f;
}

/* Get the active family --------------------------------------------------- */
void QtImage::slotActiveFamily( int _f ) {
   __activeFamily = _f;
}

/* Set Coordinates --------------------------------------------------------- */
void QtImage::slotSetCoords( int _x, int _y ) {
   __x = (int)(_x / __zoom);
   __y = (int)(_y / __zoom);
}

/* Resize event (create a new image with the new size) --------------------- */
void QtImage::resizeEvent( QResizeEvent * ) {
   __img->reset();
   __img->create( size(), 8, 256 );

   for( int i = 0; i < 256; i++ ) __img->setColor( i, qRgb( i, i, i ) );
}

/* Paint the micrograph ---------------------------------------------------- */
void QtImage::paintEvent( QPaintEvent * ) {
   loadImage();

   if ( __filtersController != NULL )
      applyFilters( __filtersController,__img );

   __paint->drawImage( 0, 0, *__img );
   loadSymbols();
   draw_axis(__axis_ang);
}

/* Change contrast --------------------------------------------------------- */
void QtImage::changeContrast(int _mingray, int _maxgray, float _gamma) {
   __mingray=_mingray;
   __maxgray=_maxgray;
   __gamma  =_gamma;

   if (__mingray < __maxgray) {
      float a =255.0/(__maxgray-__mingray);
      float ap=1.0  /(__maxgray-__mingray);
      for (int i=0; i<__mingray; i++) __img->setColor(i, qRgb(0,0,0));
      for (int i=__mingray; i<=__maxgray; i++) {
          float pregray=ap*(i-__mingray);
	  float fgray;
	  if (__gamma==1.0) fgray=pregray;
	  else              fgray=pow(pregray,__gamma);
	  int gray=(int)(255.0*fgray);
	  __img->setColor( i, qRgb(gray,gray,gray));
      }
      for (int i=__maxgray+1; i<256; i++) __img->setColor(i, qRgb(255,255,255));
   } else {
      for (int i=0; i<256; i++) __img->setColor(i, qRgb(0,0,0));
   }
   repaint(false);
}
