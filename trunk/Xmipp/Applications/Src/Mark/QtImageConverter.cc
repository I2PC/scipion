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
#include "QtImageConverter.hh"

/* Constructor ------------------------------------------------------------- */
QtImageConverter::QtImageConverter() {
}

/* qt to xmipp converter --------------------------------------------------- */
Image *QtImageConverter::qt2xmipp( QImage *qtImg ) {
   Image *xmippImg = new Image( qtImg->height(), qtImg->width() );
   
   for( int y = 0; y < qtImg->height(); y++ )
      for( int x = 0; x < qtImg->width(); x++ )
         (*xmippImg)(y, x) = qtImg->pixelIndex(x, y);

   return( xmippImg );
}

/* xmipp to qt converter --------------------------------------------------- */
QImage *QtImageConverter::xmipp2qt( Image *xmippImg ) {
   QImage *qtImg = new QImage( (*xmippImg)().ColNo(), (*xmippImg)().RowNo(), 
                               8, 256 );
   
   for( int i = 0; i < 256; i++ ) qtImg->setColor( i, qRgb( i, i, i ) );   

   for( int y = 0; y < qtImg->height(); y++ )
      for( int x = 0; x < qtImg->width(); x++ )
         qtImg->setPixel( x, y, (unsigned int)((*xmippImg)(y, x)) );
   
   return( qtImg );
}
