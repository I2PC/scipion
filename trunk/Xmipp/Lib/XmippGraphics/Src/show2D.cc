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

#include "../show2D.hh"
#include "../showTools.hh"
#include <XmippData/xmippFuncs.hh>
#include <XmippInterface/xmippVTK.hh>
#include <qmenubar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpainter.h>
#include <string.h> 

/****************************************************/

/*
  In the constructor, we just pass the standard parameters on to
  QWidget.

  The menu uses a single slot to simplify the process of adding
  more items to the options menu.
*/

void ImageViewer::Init()
{
    pickx = -1;
    clickx = -1;
    alloc_context = 0;
    down = false;
    spacing = 1;

    menubar = new QPopupMenu(); 
    
    file = new QPopupMenu();
    menubar->insertItem( "&File", file );
    file->insertItem( "New window", this,  SLOT(newWindow()));
    file->insertItem( "Open...", this,  SLOT(openFile()));
    file->insertItem( "Save image...", this,  SLOT(saveImage(int)));

    // Create and setup the print button
    pi = file->insertItem( "Print", this, SLOT(printIt()));
    printer = new QPrinter;

    options =  new QPopupMenu();
    menubar->insertItem( "&Options", options );
    ss = options->insertItem( "Set Spacing" );

    menubar->insertSeparator();

    QPopupMenu* help = new QPopupMenu();
    menubar->insertItem( "&Help", help );
    help->insertItem( "&About", this, SLOT(about()));
    help->insertItem( "About &Xmipp", this, SLOT(aboutXmipp()));
    help->insertSeparator();
    help->insertItem( "Help!", this, SLOT(giveHelp()));


    menubar->insertSeparator();
    menubar->insertItem( "Quit", this,  SLOT(close()));

    connect( options, SIGNAL(activated(int)), this, SLOT(doOption(int)) );

    status = new QLabel(this);
    status->setFrameStyle( QFrame::WinPanel | QFrame::Sunken );
    status->setFixedHeight( fontMetrics().height() + 4 );

    setMouseTracking( TRUE );
    
    if (check_file_change) {
        timer = new QTimer( this );
        connect( timer, SIGNAL(timeout()), this, SLOT(check_file()) );
        timer->start( 500 ); // Check every 0.5 seconds
    } else timer=NULL;
}


/****************************************************/
ImageViewer::ImageViewer( const char *name, bool _check_file_change): 
      QWidget( NULL, name, QWidget::WDestructiveClose ),
      filename( 0 ),
      helpmsg( 0 )
{
    check_file_change=_check_file_change;
    Init();
}

/****************************************************/

ImageViewer::ImageViewer( QImage *_image, const char *name)
    : QWidget( NULL, name, QWidget::WDestructiveClose ),
      filename( 0 ),
      helpmsg( 0 )
{
    check_file_change=false;
    Init();
    filename = name;
    if (Qt2xmipp(*_image)) showImage();
}

/****************************************************/

ImageViewer::ImageViewer( Image *_image, const char *name)
    : QWidget( NULL, name, QWidget::WDestructiveClose ),
      filename( 0 ),
      helpmsg( 0 )
{
    check_file_change=false;
    Init();
    filename = name;
    if (xmipp2Qt((Image&) *_image)) showImage();
}

/****************************************************/

ImageViewer::ImageViewer( FourierImageXmipp *_FFTimage, const char *name)
    : QWidget( NULL, name, QWidget::WDestructiveClose ),
      filename( 0 ),
      helpmsg( 0 )
{
    check_file_change=false;
    Init();
    filename = name;
    // Compute the magnitude
    ImageXmipp I;
    FFT_magnitude(*_FFTimage,I());
    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(I())
       MULTIDIM_ELEM(I(),i)=
	  log10(1+MULTIDIM_ELEM(I(),i)*MULTIDIM_ELEM(I(),i));
    // Set the image and show
    if (xmipp2Qt(I)) showImage();
}

/****************************************************/

ImageViewer::~ImageViewer()
{
    if ( alloc_context )
	QColor::destroyAllocContext( alloc_context );
    if ( other == this )
	other = 0;
}

/****************************************************/

/*
	This function set the spacing constant.
*/
void ImageViewer::doOption(int item)
{
    if (item == ss) {
      float min= 0.1;
      float max= 10.0;        
      ScrollParam* param_window;	
      param_window = new ScrollParam(min,max, spacing, "Set spacing", 
         0, "new window", WDestructiveClose);
      connect( param_window, SIGNAL(new_value(float)), this, SLOT(set_spacing(float)) );
      param_window->setFixedSize(250,200);
      param_window->show();
    }
}

/****************************************************/

void ImageViewer::updateStatus()
{
    if ( pm.size() == QSize( 0, 0 ) ) {
	if ( filename )
	    status->setText("Could not load image");
	else
	    status->setText("No image - select Open from File menu.");
    } else {
	QString message, moremsg;
	
	if (image.valid(pickx,picky)) {
	    int y_log, x_log;
	    xmippImage().physical2logical(picky, pickx, y_log, x_log);
	    moremsg.sprintf("(%d,%d)=(%d,%d)= %.3f ",
	                  picky, pickx,
			  y_log, x_log,
			  xmippImage(y_log,x_log));
	    message += moremsg;
	}
	moremsg.sprintf("%dx%d", image.width(), image.height());
	message += moremsg;
	if ( pm.size() != pmScaled.size() ) {
	    moremsg.sprintf(" [%dx%d]", pmScaled.width(),
		pmScaled.height());
	    message += moremsg;
	}
	status->setText(message);
    }
}

/****************************************************/

/*
  This function saves the image.
*/
void ImageViewer::saveImage( int item )
{
//    QString savefilename = QFileDialog::getSaveFileName(0, 0, 0, filename);
    QString savefilename = QFileDialog::getSaveFileName(QString::null, "*.xmp", this);

    if ( !savefilename.isEmpty() ) {
	  try {
	    ImageXmipp tmpImage; tmpImage()=xmippImage();
 	    // Saves Xmipp Image
	    tmpImage.rename((string) ((const char *)savefilename));
            tmpImage.write(); 

	  } catch (Xmipp_error) {
	      char *helptext = "Invalid image type";
	      helpmsg = new QMessageBox( "Error", helptext,
	      	QMessageBox::Information, QMessageBox::Ok, 0, 0, 0, 0, FALSE );
    	      helpmsg->show();
    	      helpmsg->raise();
	  }
    }    
    
}


/****************************************************/

void ImageViewer::newWindow()
{
    ImageViewer* that = new ImageViewer("new window");
    that->show();
}

/*
  This function is the slot for processing the Open menu item.
*/
void ImageViewer::openFile()
{  
    QString newfilename = QFileDialog::getOpenFileName();
    if ( !newfilename.isEmpty() ) {
	loadImage( newfilename ) ;
	repaint();			// show image in widget
    }
}


/****************************************************/
//
// Called when the print button is clicked.
//

void ImageViewer::printIt()
{
    if ( printer->setup(this) ) {
	QPainter paint( printer );
 	paint.drawPixmap(0, 0, pmScaled);	
    }
}



/*****************************************/

bool ImageViewer::showImage()
{
    bool ok = FALSE;
    QApplication::setOverrideCursor( waitCursor ); // this might take time
    pickx = -1;
    clickx = -1;
    ok = reconvertImage();
    if ( ok ) {     
	setCaption( filename );			// set window caption
        int w = pm.width();
    	int h = pm.height();

    	const int reasonable_width = 128;
    	if ( w < reasonable_width ) {
    	    // Integer scale up to something reasonable
    	    int multiply = ( reasonable_width + w - 1 ) / w;
    	    w *= multiply;
    	    h *= multiply;
    	}

    	h += status->height();
    	resize( w, h ); 			    // we resize to fit image
    } else {
    	pm.resize(0,0); 			    // couldn't load image
    	update();
    }
    QApplication::restoreOverrideCursor();  // restore original cursor
    updateStatus();
//    setMenuItemFlags();
    return ok;
}


/*****************************************/

bool ImageViewer::xmipp2Qt(Image& _image )
{
    bool ok = FALSE;
		
    try { 
      xmippImage = _image;
      ::xmipp2Qt(_image,image); // Take the one in showTools
      xmippFlag = 0;	 	// Sets flag = Xmipp image.
      ok = TRUE; 
    } catch (Xmipp_error) {
      ok = FALSE;
      char *helptext = "Error converting xmipp to Qt";
      helpmsg = new QMessageBox( "Error", helptext,
    	QMessageBox::Information, QMessageBox::Ok, 0, 0, 0, 0, FALSE );
      helpmsg->show();
      helpmsg->raise();
    }
  
    return ok;
}

/*****************************************/

bool ImageViewer::Qt2xmipp( QImage &_image )
{
    bool ok = FALSE;
    // try to read image from standard format.
		
    try { 
      image = _image;      
      image.setNumColors(256);

      ::Qt2xmipp(_image,xmippImage); // Take the one in showTools
      xmippFlag = 0;	 	// Sets flag = Xmipp image.
      
      ok = TRUE;

    } catch (Xmipp_error) {
      ok = FALSE;
      char *helptext = "Error converting Qt image to Xmipp image";
      helpmsg = new QMessageBox( "Error", helptext,
    	QMessageBox::Information, QMessageBox::Ok, 0, 0, 0, 0, FALSE );
      helpmsg->show();
      helpmsg->raise();
    }
  
    return ok;
}


/*****************************************/

/*
  This function loads an image from a file and resizes the widget to
  exactly fit the image size. If the file was not found or the image
  format was unknown it will resize the widget to fit the errorText
  message (see above) displayed in the current font.

  Returns TRUE if the image was successfully loaded.
*/

bool ImageViewer::loadImage( const char *fileName ) 
{
    filename = fileName;
    bool imagic=((string)(filename)).find("imagic:")==0;
    bool ok = FALSE;
    static bool message_shown=false;
    if ( filename ) {
	
	// try to read image from standard format.	
	
	if (image.load(filename, 0)) ok = Qt2xmipp(image);	  
	
	if (!ok) {
          try { 
 	    // reads Xmipp Image
            Image tmpImage;
	    if (!imagic) wait_until_stable_size(filename);
            if (imagic) {
	       Image *p = Image::LoadImage(filename);
	       if (!p) REPORT_ERROR(1,"ImageViewer::loadImage: Unknown format");
               tmpImage() = (*p)();
               delete p;
            } else if (Is_ImageXmipp(filename)) {
               ImageXmipp p;
               p.read((FileName)filename);
               tmpImage()=p();
            } else if (Is_FourierImageXmipp(filename)) {
	       FourierImageXmipp If; If.read(filename);
	       FFT_magnitude(If,tmpImage());
	       FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(tmpImage())
                  MULTIDIM_ELEM(tmpImage(),i)=
	             log10(1+MULTIDIM_ELEM(tmpImage(),i)*
		             MULTIDIM_ELEM(tmpImage(),i));
	    } else REPORT_ERROR(1,"ImageViewer::loadImage: Unknown format");
	    tmpImage().set_Xmipp_origin();
	    ok = xmipp2Qt(tmpImage);
	  } catch (Xmipp_error) {
	    ok = FALSE;
	    if (!message_shown) {
	       char *helptext = "Invalid image type";
	       helpmsg = new QMessageBox( "Error", helptext,
		 QMessageBox::Information, QMessageBox::Ok, 0, 0, 0, 0, FALSE );
    	       helpmsg->show();
    	       helpmsg->raise();
	       message_shown=true;
	    }
	  }
	}
    }   
    
    if (ok) {
       ok = showImage();
       struct stat info;
       if (stat(filename, &info) && !imagic) {
          cerr << "loadImage: Cannot get time of file " << filename << endl;
          modification_time=0;
       } else modification_time=info.st_mtime;
    }
    
    return ok;
}

/*****************************************/
void ImageViewer::setImage( const matrix2D<double> &img ) {
    xmippImage()=img;
    showImage();
}

/*****************************************/

bool ImageViewer::reconvertImage()
{
    bool success = FALSE;

    if ( image.isNull() ) return FALSE;

    QApplication::setOverrideCursor( waitCursor ); // this might take time
    if ( pm.convertFromImage(image) )
    {
	pmScaled = QPixmap();
	scale();
	resize( width(), height() );
	success = TRUE;				// load successful
    } else {
	pm.resize(0,0);				// couldn't load image
    }
    updateStatus();
    QApplication::restoreOverrideCursor();	// restore original cursor

    return success;				// TRUE if loaded OK
}

/****************************************************/


/*
  This functions scales the pixmap in the member variable "pm" to fit the
  widget size and  puts the resulting pixmap in the member variable "pmScaled".
*/

void ImageViewer::scale()
{
    int h = height() - status->height();

    if ( image.isNull() ) return;

    QApplication::setOverrideCursor( waitCursor ); // this might take time
    if ( width() == pm.width() && h == pm.height() )
    {						// no need to scale if widget
	pmScaled = pm;				// size equals pixmap size
    } else {
       QWMatrix m;			   // transformation matrix
       m.scale(((double)width())/pm.width(),// define scale factors
	       ((double)h)/pm.height());
       pmScaled = pm.xForm( m );	   // create scaled pixmap
    }
    QApplication::restoreOverrideCursor();	// restore original cursor
}


/****************************************************/

/*
  The resize event handler, if a valid pixmap was loaded it will call
  scale() to fit the pixmap to the new widget size.
*/

void ImageViewer::resizeEvent( QResizeEvent * )
{
    status->setGeometry(0, height() - status->height(),
		        width(), status->height());

    if ( pm.size() == QSize( 0, 0 ) )		// we couldn't load the image
	return;

    int h = height() - status->height();
    if ( width() != pmScaled.width() || h != pmScaled.height())
    {						// if new size,
	scale();				// scale pmScaled to window
	updateStatus();
    }
}


/****************************************************/

/*
  Handles a change in the mouse
*/

bool ImageViewer::convertEvent( QMouseEvent* e, int& x, int& y)
{
    if ( pm.size() != QSize( 0, 0 ) ) {
	int h = height() - status->height();
	int nx = e->x() * image.width() / width();
	int ny = (e->y()) * image.height() / h;
	if (nx != x || ny != y ) {
	    x = nx;
	    y = ny;
	    return TRUE;
	}
    }
    return FALSE;
}


/****************************************************/

/*
  Mouse press events.
*/

void ImageViewer::mousePressEvent( QMouseEvent *e )
{
    QPoint clickedPos = e->pos();		// extract pointer position
    if (e->button() == RightButton) { 
      menubar->exec(clickedPos);
      down = false;
    } else {
    	if (!down) {    
  		if (convertEvent(e,xi, yi)) {
	        	repaint();
	        	xir = e->x();   
	        	yir = e->y(); 
			ox = xir;
			oy = yir;  
    			down = true;  
		}
    	}
    }      
}


/****************************************************/

/*
  Mouse release events.
*/

void ImageViewer::mouseReleaseEvent( QMouseEvent *e )
{
    if (down) {    
  	if (convertEvent(e,xf,yf)) {
	        xfr = e->x();   
	        yfr = e->y();
		if (xf > xmippImage().ColNo()) xf = xmippImage().ColNo();   
		if (xf <  0) xf = 0;
		if (yf > xmippImage().RowNo()) yf = xmippImage().RowNo();   
		if (yf <  0) yf = 0;
       		float distance = sqrt ((double) (xf -xi)*(xf -xi) + (yf -yi)*(yf -yi));	
		QString message;
	    	message.sprintf("Distance: %.3f Angstroms", distance*spacing);
		status->setText(message);
       		down = false;
	}
    }
}


/****************************************************/

/*
  Record the pixel position of interest.
*/
void ImageViewer::mouseMoveEvent( QMouseEvent *e )
{
  if (convertEvent(e,pickx,picky)) {
    updateStatus();
    if (down) {
        if ((e->x() < width()) && (e->y() < height()-(status->height()))) {
  		QPainter p(this);
    		QBrush brush( NoBrush );  // void brush
    		QPen myPen(red, 3);
    		p.setPen( myPen );	
    		p.setBrush( brush );
		p.setRasterOp(XorROP);
		p.drawLine(xir, yir, ox, oy);      
        	p.drawLine(xir, yir, e->x(), e->y());
		ox = e->x(); oy = e->y();      
		QString message;
       		float distance = sqrt ((double) (pickx -xi)*(pickx -xi) + (picky -yi)*(picky -yi));	
	    	message.sprintf("Distance: %.3f Angstroms", distance*spacing);
		status->setText(message);
	}
    }
  }
}



/****************************************************/

/*
  Handles key press events.
*/

void ImageViewer::keyPressEvent( QKeyEvent* e )
{
    switch( e->key() ) {			// Look at the key code
	case Key_R:
             if (e->state() == ControlButton) {	// If 'Ctrol R' key, 
  		  xmippImage().move_origin_to(-xmippImage().startingY(), -xmippImage().startingX());// sets origin at the upper left corner        
	     }
	     break;
	case Key_Q:
             if (e->state() == ControlButton) {	// If 'Ctrol Q' key, 
  		  exit(0); // Terminate program
	     }
	     break;
	case Key_O:				// Xmipp origin
             if (e->state() == ControlButton) {	// If 'Ctrol N' key, 
  		  xmippImage().set_Xmipp_origin(); // sets origin at the center of the iamge.        
	     }
	     break;
	case Key_N:				// Natural size (original size)
             if (e->state() == ControlButton) {	// If 'Ctrol N' key, 
    		  resize(xmippImage().ColNo(), xmippImage().RowNo() + status->height());	         
	     }
	     break;
	case Key_M:     
	case Key_Minus:				// Half size
             if (e->state() == ControlButton) {	// If 'Ctrol+' key, 
    		  resize(width()/2, height()/2 + status->height()/2 + 1);	         
	     }
	     break;
	case Key_P:
	case Key_Plus:				// Double size
             if (e->state() == ControlButton) {	// If 'Ctrol+' key, 
    		  resize(width()*2, height()*2 - status->height());
	     }
	     break;
	case Key_A:    				// Aspect ratio
             if (e->state() == ControlButton) {	// If 'Ctrol+' key, 
	          double ratio = (double) xmippImage().ColNo()/ (double) xmippImage().RowNo();		  
    		  resize(width(), (int) (width()/ratio + status->height()));	         
	     }
	     break;
	default:				// If not an interesting key,
	    e->ignore();			// we don't accept the event
	    return;	
    }
}


/****************************************************/

/*
  Draws the portion of the scaled pixmap that needs to be updated or prints
  an error message if no legal pixmap has been loaded.
*/

void ImageViewer::paintEvent( QPaintEvent *e )
{
    if ( pm.size() != QSize( 0, 0 ) ) {		// is an image loaded?
	QPainter painter(this);
	painter.setClipRect(e->rect());
 	painter.drawPixmap(0, 0, pmScaled);
   }
}


/****************************************************/

/*
  Explain anything that might be confusing.
*/
void ImageViewer::giveHelp()
{
    if (!helpmsg) {
	QString helptext = "Usage: xmipp_iv [filename]\n\n ";
    	QStrList support = QImage::outputFormats();
	helptext += "\n\nSupported input formats:\n";
	int lastnl = helptext.length();

    	support.clear();
    	support.insert(0, "Spider");
    	support.insert(1, "bmp");
	
	const char* f = support.first();
	helptext += f;
	f = support.next();
	for (; f; f = support.next()) {
	    helptext += ',';
	    if ( helptext.length() - lastnl > 40 ) {
		helptext += "\n  ";
		lastnl = helptext.length() - 2;
	    } else {
		helptext += ' ';
	    }
	    helptext += f;
	}
	helptext += "\n\nCommands:\n";
	helptext += " Right-click : Popup menu\n";
	helptext += " Ctrl Q : Quit\n";
	helptext += " Ctrl N : Natural size of the image\n";
	helptext += " Ctrl O : Set origin to the center of the image\n";
	helptext += " Ctrl R : Restore origin to the upper left corner of the image\n";
	helptext += " Ctrl - or Ctrl M: Half the size of the image\n";
	helptext += " Ctrl + or Ctrl P: Double the size of the image\n";
	helptext += " Ctrl A : Aspect ratio\n";
	helpmsg = new QMessageBox( "Help", helptext,
	    QMessageBox::Information, QMessageBox::Ok, 0, 0, 0, 0, FALSE );
    }
    helpmsg->show();
    helpmsg->raise();
}


void ImageViewer::about()
{
    QMessageBox::about( this, "IV (Image Viewer)",
			"Visualizes an Image in Spider format. \n");
}


void ImageViewer::aboutXmipp()
{
    QMessageBox::about( this, "Xmipp: Xmipp Image Processing Package",
    				"Biocomputing Unit.\n"
				"National Center of Biotechnology-CSIC\n"
				"Madrid, Spain\n"
				"http://www.biocomp.cnb.uam.es\n");
}

/****************************************************/

void ImageViewer::set_spacing(float _spacing) {
  spacing = _spacing;
}


/****************************************************/

ImageViewer* ImageViewer::other = 0;

/****************************************************/

void ImageViewer::check_file() {
   struct stat info;
   static bool message_shown=false;
   if (stat(filename, &info) && !message_shown) {
      cerr << "check_file: Cannot get time of file "<< filename << endl;
      message_shown=true;
   }
   if (info.st_mtime!=modification_time) {
      loadImage(filename);
      repaint();
      updateStatus();
   }
}
