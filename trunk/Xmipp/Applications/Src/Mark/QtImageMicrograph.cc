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
#include "QtImageMicrograph.hh"
#include "QtColorLabel.hh"
#include "QtDialogProperties.hh"
#include "QtWidgetMicrograph.hh"
#include "XmippData/xmippMicrograph.hh"
#include "qpainter.h"
#include "../DownSample/xvsmooth.h"

/* Constructor ------------------------------------------------------------- */
QtImageMicrograph::QtImageMicrograph( QWidget *_parent = 0, 
                                      const char *_name = 0, 
                                      WFlags _f = 0 ) : 
   QtImage( _parent, _name, _f ) { 
   __x          = 0; 
   __y          = 0;
   __pressed    = false;
   __movingMark = -1;
   __tilted     = false; 
   
   emit signalSetWidthHeight( image()->width(), image()->height() );
}

/* Set the micrograph ------------------------------------------------------ */
void QtImageMicrograph::setMicrograph( Micrograph *_m ) {
   QtImage::setMicrograph( _m );
   
   emit signalSetWidthHeight( image()->width(), image()->height() );   
}


/* Load Image -------------------------------------------------------------- */
void QtImageMicrograph::loadImage() {
   int mMaxX, mMaxY, mX, mY;
   double emX, emY;
   
   if ( getMicrograph() != NULL ) getMicrograph()->size( mMaxX, mMaxY );
   else return;

   // Get the starting and finishing point in the original micrograph
   int mY0, mX0, mYF, mXF;
   imageToMicrograph( 0, 0, mX0, mY0 );
   imageToMicrograph( image()->width(), image()->height(), mXF, mYF );
   
   if (mXF-mX0<=image()->width() && mYF-mY0<=image()->height() &&
       getMicrograph()->depth()==8) {
      // Use XV for showing the image
      // Copy that piece of the micrograph in an intermidiate piece of memory
      byte *ptr,*piece=new byte[image()->width() * image()->height()];
      if (piece==NULL)
          REPORT_ERROR(1,"QtImageMicrograph::loadImage: Cannot allocate memory");
      ptr=piece;
      for (mY=mY0; mY<mYF; mY++)
         for (mX=mX0; mX<mXF; mX++)
            *ptr++=getMicrograph()->val8(mX, mY);

      // Apply xvsmooth and copy to the canvas
      byte rgb[256]; for (int i=0; i<256; i++) rgb[i]=i;
      byte *result = SmoothResize(piece,
         mXF-mX0, mYF-mY0, image()->width(), image()->height(),
         rgb, rgb, rgb, rgb, rgb, rgb, 256);
      ptr=result;
      for( int y = 0; y < image()->height(); y++ )
         for( int x = 0; x < image()->width(); x++ )
            image()->setPixel( x, y, *ptr++);
      free(result);
      delete piece;
   } else {
      // Apply bilinear interpolation
      for( int y = 0; y < image()->height(); y++ )
         for( int x = 0; x < image()->width(); x++ )
            if ( getMicrograph() != NULL ) {
               exact_imageToMicrograph( x, y, emX, emY );
               double val=0;
               if (emX>=0 && emX<mMaxX-1 && emY>=0 && emY<mMaxY-1) {
                  double wx=emX-(int)emX;
                  double wy=emY-(int)emY;
                  int    mX1=(int)emX, mX2=mX1+1;
                  int    mY1=(int)emY, mY2=mY1+1;
                  val+=(1-wy)*(1-wx)*getMicrograph()->val8(mX1,mY1)+
                       (1-wy)*   wx *getMicrograph()->val8(mX2,mY1)+
                          wy *(1-wx)*getMicrograph()->val8(mX1,mY2)+
                          wy *   wx *getMicrograph()->val8(mX2,mY2);
               }
               image()->setPixel( x, y, (unsigned int)val );
            } else image()->setPixel( x, y, 0 );
   }
}

/* Draw ellipse ------------------------------------------------------------ */
void QtImageMicrograph::drawEllipse(int _x, int _y, int _color) {
   int mX, mY;
   #define ellipse_radius 5.0
   micrographToImage( _x, _y, mX, mY );
   if ( (mX > 0 && mX < image()->width()) && 
        (mY > 0 && mY < image()->height()) ) {
      QPen pen( __col.col(_color), 2 );
      __paint->setPen( pen );
      __paint->drawEllipse(ROUND(mX-ellipse_radius/__zoom),
         ROUND(mY-ellipse_radius/__zoom),
      	 ROUND(2*ellipse_radius/__zoom), ROUND(2*ellipse_radius/__zoom) );
      __paint->flush();
   }
}

/* Load Symbols ------------------------------------------------------------ */
void QtImageMicrograph::loadSymbols() {
   if ( getMicrograph() == NULL ) return;
   for( int i = 0; i < getMicrograph()->ParticleNo(); i++ ) {
      if ( !getMicrograph()->coord(i).valid ) continue;      
      drawEllipse(getMicrograph()->coord(i).X, 
         getMicrograph()->coord(i).Y,getMicrograph()->coord(i).label);
   }
}

void QtImageMicrograph::slotDeleteMarkOther( int _coord ) {
   emit signalDeleteMarkOther( _coord );
}

void QtImageMicrograph::slotChangeFamilyOther( int _coord, int _f ) {
   emit signalChangeFamilyOther( _coord, _f );
}

void QtImageMicrograph::resizeEvent( QResizeEvent *e ) {
   QtImage::resizeEvent( e );
   
   if ( getMicrograph() == NULL ) return;
      
   emit signalSetWidthHeight( image()->width(), image()->height() );
   emit signalRepaint();   
}

void QtImageMicrograph::mousePressEvent( QMouseEvent *e ) {
   if ( e->button() == LeftButton ) __pressed = true;
}

void QtImageMicrograph::mouseReleaseEvent( QMouseEvent *e ) {
   if ( getMicrograph() == NULL ) return;

   // If moving a particle
   if ( __movingMark != -1 && __pressed ) {
      __movingMark = -1;
      __pressed    = false;
      emit signalRecalculateTiltMatrix();
      return;
   } else if ( __movingMark != -1 ) return;
   
   // If picking a new one
   int mX, mY;
   imageToMicrograph( e->pos().x(), e->pos().y(), mX, mY );
   
   if ( e->button() == RightButton ) changeProperties( mX, mY );
   else if ( __pressed == true ) {
      if ( isTilted() ) {
         cout << "Moving last particle to (X,Y)=(" << mX << "," << mY << ")\n";
	 getMicrograph()->move_last_coord_to(mX,mY);
	 __pressed = false;
         emit signalRepaint();   
         emit signalRecalculateTiltMatrix();
      } else {
	 cout << "Particle marked at (X,Y)=(" << mX << "," << mY << ")\n";
	 getMicrograph()->add_coord( mX, mY, __activeFamily );    
	 __pressed = false;
	 emit signalAddCoordOther( mX, mY, __activeFamily );
      }
   }
}

void QtImageMicrograph::mouseMoveEvent( QMouseEvent *e ) {
   if ( getMicrograph() == NULL || __movingMark == -1 || !__pressed ) return;
   
   int mX, mY;
   imageToMicrograph( e->pos().x(), e->pos().y(), mX, mY );
   getMicrograph()->coord(__movingMark).X = mX;
   getMicrograph()->coord(__movingMark).Y = mY;
   
   emit signalRepaint();   
}

void QtImageMicrograph::changeProperties( int mX, int mY ) {
   int coord = getMicrograph()->search_coord_near( mX, mY, 10 );
   if ( coord == -1 ) return;
   
   QtDialogProperties dialogProperties( getMicrograph(), coord, this, 0, TRUE );
   
   connect( &dialogProperties, SIGNAL(signalDeleteMarkOther(int)),
            this, SLOT(slotDeleteMarkOther(int)) );
   connect( &dialogProperties, SIGNAL(signalChangeFamilyOther(int,int)),
            this, SLOT(slotChangeFamilyOther(int,int)) );

   dialogProperties.exec();
}

void QtImageMicrograph::slotZoomIn() {
   if ( __zoom == 0.1 ) return;
   __zoom -= 0.1;
   
   __x = __y = 0;
   emit signalSetCoords( 0, 0 );
   emit signalRepaint();
}

void QtImageMicrograph::slotZoomOut() {
   __zoom += 0.1;
   
   __x = __y = 0;
   emit signalSetCoords( 0, 0 );
   emit signalRepaint();
}
