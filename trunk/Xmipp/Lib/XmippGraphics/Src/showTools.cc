/***************************************************************************
 *
 * Authors:      Alberto Pascual Montano (pascual@cnb.uam.es)
 *               Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#include "../showTools.hh"
#include <qlayout.h>
#include <qscrollbar.h>
#include <qpushbutton.h>
#include <qtooltip.h> 

ScrollParam::ScrollParam(float min , float max, float initial_value,
   char *caption, QWidget *parent, const char *name, int wFlags, int precision  ) : 
   QWidget( parent, name, wFlags )
{
    my_precision = (int)pow ((double)10,(double)precision);
    int tmp_min=(int) (min*my_precision);
    int tmp_max=(int) (max*my_precision);
    value  = (float)initial_value*my_precision;
    QColor col;
   // Set the window caption/title

    setCaption( caption );

    // Create a layout to position the widgets

    QBoxLayout *Layout = new QVBoxLayout( this, 10 );

    // Create a grid layout to hold most of the widgets
    QGridLayout *grid = new QGridLayout( 6, 4 );
    // This layout will get all of the stretch

    Layout->addLayout( grid, 5 );

    //title
    QLabel     *title= new QLabel( this, "title" );    
    title->setFont( QFont("times",14,QFont::Bold) );
    title->setText( "Spacing selector" );
    title->setFixedSize(title->sizeHint());
    grid->addMultiCellWidget( title, 0, 0, 0, 2);
     
    //Scroll Bar
    QScrollBar  *scroll= new QScrollBar(tmp_min, tmp_max, 1, 1, (int)value, QScrollBar::Horizontal,this,"scroll");
    scroll->setFixedWidth(200); 
    scroll->setFixedHeight(15); 
    connect( scroll, SIGNAL(valueChanged(int)), SLOT(scrollValueChanged(int)) );
    grid->addMultiCellWidget( scroll, 1, 1, 0, 2);			   

    //labels under scroll bar
    QLabel     *lab1= new QLabel( this, "lab1" );    
    lab1->setFont( QFont("times",12,QFont::Bold) );
    lab1->setNum( (float)tmp_min/my_precision);
    lab1->setFixedSize(lab1->sizeHint());
    grid->addWidget( lab1, 2, 0, AlignLeft );
    QLabel     *lab2= new QLabel( this, "lab2" );        
    lab2->setFont( QFont("times",12,QFont::Bold) );
    lab2->setNum( (float)tmp_max/my_precision);
    lab2->setFixedSize(lab2->sizeHint());
    grid->addWidget( lab2, 2, 2, AlignRight );

    //Labels for the current value
    QLabel     *lab3= new QLabel( this, "lab3" );    
    lab3->setFont( QFont("times",12,QFont::Bold) );
    lab3->setText( "Spacing: ");
    lab3->setFixedSize(lab3->sizeHint());
    grid->addWidget( lab3, 3, 0, AlignCenter );
    value_lab= new QLabel( this, "value_lab" );        
    value_lab->setFont( QFont("times",12,QFont::Bold) );
    value_lab->setNum( value/my_precision );
    grid->addWidget( value_lab, 3, 1, AlignLeft);

    //Close Button
    QPushButton *close;
    close = new QPushButton( this, "close" );	// create button 1
    close->setFont( QFont("times",12,QFont::Bold) );
    close->setText( "Close" );    
    close->setFixedSize( close->sizeHint());
    grid->addWidget( close, 4, 0, AlignVCenter ); 
    QToolTip::add( close, "Close the window" );
    connect( close, SIGNAL(clicked()), SLOT(but_close_clicked()) );
        
    QPushButton *do_it;
    do_it = new QPushButton( this, "do_it" );	// create button 3
    do_it->setFont( QFont("times",12,QFont::Bold) );
    do_it->setText( "Ok" );
    do_it->setFixedHeight( do_it->sizeHint().height());
    do_it->setFixedWidth(80);
    grid->addWidget( do_it, 4, 2, AlignVCenter ); 
    QToolTip::add( do_it, "Commit action" );
    connect( do_it, SIGNAL(clicked()), SLOT(but_ok_clicked()) );
    
    setFixedSize(350,250);
}

/****************************************************/
ScrollParam::~ScrollParam()
{
}
/****************************************************/

void ScrollParam::scrollValueChanged(int v){
   value = (float)v/my_precision;
   value_lab->setNum(value);
}

/****************************************************/

void ScrollParam::but_close_clicked()
{
  close();
}

/****************************************************/

void ScrollParam::but_ok_clicked()   
{
  emit new_value( value );
  close();  
}

ScrollParam2::ScrollParam2(float min , float max, float initial_value,
                         float initial_value2,
   char *caption, QWidget *parent, const char *name, int wFlags, int precision  ) : 
   QWidget( parent, name, wFlags )
{
    my_precision = (int)pow ((double)10,(double)precision);
    int tmp_min=(int) (min*my_precision);
    int tmp_max=(int) (max*my_precision);
    value  = (float)initial_value*my_precision;
    value2 = (float)initial_value2*my_precision;;
    QColor col;
   // Set the window caption/title

    setCaption( caption );

    // Create a layout to position the widgets

    QBoxLayout *Layout = new QVBoxLayout( this, 10 );

    // Create a grid layout to hold most of the widgets
    QGridLayout *grid = new QGridLayout( 6, 4 );
    // This layout will get all of the stretch

    Layout->addLayout( grid, 5 );

    //title
    QLabel     *title= new QLabel( this, "title" );    
    title->setFont( QFont("times",12,QFont::Bold) );
    title->setText( "Spacing and Offset X-Selector" );
    title->setFixedSize(title->sizeHint());
    grid->addMultiCellWidget( title, 0, 0, 0, 2);
     
    //Scroll Bar
    QScrollBar  *scroll= new QScrollBar(tmp_min, tmp_max, 1, 1, (int)value, QScrollBar::Horizontal,this,"scroll");
    scroll->setFixedWidth(200); 
    scroll->setFixedHeight(15); 
    connect( scroll, SIGNAL(valueChanged(int)), SLOT(scrollValueChanged(int)) );
    grid->addMultiCellWidget( scroll, 1, 1, 0, 2);			   

    //Scroll Bar2
    QScrollBar  *scroll2= new QScrollBar(tmp_min, tmp_max, 1, 1, (int)value2, QScrollBar::Horizontal,this,"scroll");
    scroll2->setFixedWidth(200); 
    scroll2->setFixedHeight(15); 
    connect( scroll2, SIGNAL(valueChanged(int)), SLOT(scrollValueChanged2(int)) );
    grid->addMultiCellWidget( scroll2, 2, 1, 0, 2);			   

    //labels under scroll bar
    QLabel     *lab1= new QLabel( this, "lab1" );    
    lab1->setFont( QFont("times",12,QFont::Bold) );
    lab1->setNum( (float)tmp_min/my_precision);
    lab1->setFixedSize(lab1->sizeHint());
    grid->addWidget( lab1, 3, 0, AlignLeft );
    QLabel     *lab2= new QLabel( this, "lab2" );	 
    lab2->setFont( QFont("times",12,QFont::Bold) );
    lab2->setNum( (float)tmp_max/my_precision);
    lab2->setFixedSize(lab2->sizeHint());
    grid->addWidget( lab2, 3, 2, AlignRight );

    //Labels for the current value
    QLabel     *lab3= new QLabel( this, "lab3" );    
    lab3->setFont( QFont("times",12,QFont::Bold) );
    lab3->setText( "Spac/tick: ");
    lab3->setFixedSize(lab3->sizeHint());
    grid->addWidget( lab3, 4, 0, AlignCenter );
    value_lab= new QLabel( this, "value_lab" );        
    value_lab->setFont( QFont("times",12,QFont::Bold) );
    value_lab->setNum( value/my_precision );
    grid->addWidget( value_lab, 4, 1, AlignLeft);
    value_lab2= new QLabel( this, "value_lab2" );        
    value_lab2->setFont( QFont("times",12,QFont::Bold) );
    value_lab2->setNum( value2/my_precision );
    grid->addWidget( value_lab2, 4, 2, AlignLeft);

    //Close Button
    QPushButton *close;
    close = new QPushButton( this, "close" );	// create button 1
    close->setFont( QFont("times",12,QFont::Bold) );
    close->setText( "Close" );    
    close->setFixedSize( close->sizeHint());
    grid->addWidget( close, 5, 0, AlignVCenter ); 
    QToolTip::add( close, "Close the window" );
    connect( close, SIGNAL(clicked()), SLOT(but_close_clicked()) );
        
    QPushButton *do_it;
    do_it = new QPushButton( this, "do_it" );	// create button 3
    do_it->setFont( QFont("times",12,QFont::Bold) );
    do_it->setText( "Ok" );
    do_it->setFixedHeight( do_it->sizeHint().height());
    do_it->setFixedWidth(80);
    grid->addWidget( do_it, 5, 2, AlignVCenter ); 
    QToolTip::add( do_it, "Commit action" );
    connect( do_it, SIGNAL(clicked()), SLOT(but_ok_clicked()) );
    
    setFixedSize(350,250);
}

/****************************************************/
ScrollParam2::~ScrollParam2()
{
}
/****************************************************/

void ScrollParam2::scrollValueChanged(int v){
   value = (float)v/my_precision;
   value_lab->setNum(value);
}

void ScrollParam2::scrollValueChanged2(int v){
   value2 = (float)v/my_precision;
   value_lab2->setNum(value2);
}

/****************************************************/

void ScrollParam2::but_close_clicked()
{
  close();
}

/****************************************************/

void ScrollParam2::but_ok_clicked()   
{
  emit new_value( value, value2 );
  close();  
}

/* Qimage -> Xmipp --------------------------------------------------------- */
void Qt2xmipp( QImage &_qimage, Image &_ximage ) {
   _ximage().resize(_qimage.height(), _qimage.width());
   for (int y = 0; y < _qimage.width(); y++)
     for (int x = 0; x < _qimage.height(); x++)
       _ximage(x,y) = (double) _qimage.pixelIndex(y, x);
}

/* Xmipp -> QImage --------------------------------------------------------- */
void xmipp2Qt(Image& _ximage, QImage &_qimage, int _min_scale,
   int _max_scale) {
   // Creates a Qt Image to hold Xmipp Image
   _qimage.create(_ximage().ColNo(), _ximage().RowNo(), 8, 256);

   // Sets Graylevel Palette.
   for (int i = 0; i < 256; i++) {
     QColor c;
     c.setRgb(i,i,i);
     _qimage.setColor(i, c.rgb());   
   }  

   double min_val, max_val;
   _ximage().compute_double_minmax(min_val,max_val);
   double a=(_max_scale-_min_scale)/(max_val-min_val);
   
   // Reads pixels.
   for (int y = 0; y < _ximage().RowNo(); y++)
     for (int x = 0; x < _ximage().ColNo(); x++)
       _qimage.setPixel(x,y,((uint) (a*(
          DIRECT_MAT_ELEM(_ximage(),y,x)-min_val)+_min_scale)));		    
}

/* Xmipp -> Pixmap --------------------------------------------------------- */
void xmipp2Pixmap(Image &xmippImage, QPixmap* pixmap,
   int _minScale, int _maxScale) {
   QImage tmpImage;
   xmipp2Qt(xmippImage,tmpImage,_minScale, _maxScale);
   pixmap->convertFromImage(tmpImage, 0); 
}
