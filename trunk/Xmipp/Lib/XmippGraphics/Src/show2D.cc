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
#include <XmippData/xmippHistograms.hh>
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
    minGray=maxGray=0;
    fft_show_mode=0;
   
    pickx = -1;
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
    ravg = options->insertItem( "Radial average" );
    line_setup = options->insertItem ( "Line setup" ); 
    profile = options->insertItem( "Profile" );
    sfft = options->insertItem( "Set FFT show mode" );
    options->setItemEnabled(sfft,false);

    // Add CTF actions
    editctfmodel=options->insertItem( "Edit CTF model");
    options->setItemEnabled(editctfmodel,false);

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

    xi=yi=xf=yf=0;
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
    apply_geo=true;
    isFourierImage=false;
    load_mode=ImageViewer::Normal_mode;
    Init();
}

/****************************************************/

ImageViewer::ImageViewer( QImage *_image, const char *name)
    : QWidget( NULL, name, QWidget::WDestructiveClose ),
      filename( 0 ),
      helpmsg( 0 )
{
    check_file_change=false;
    apply_geo=true;
    isFourierImage=false;
    load_mode=ImageViewer::Normal_mode;
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
    apply_geo=true;
    isFourierImage=false;
    load_mode=ImageViewer::Normal_mode;
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
    apply_geo=true;
    isFourierImage=true;
    load_mode=ImageViewer::Normal_mode;
    Init();
    filename = name;
    xmippImageFourier=(*_FFTimage)();
    CenterFFT(xmippImageFourier,true);
    ImageXmipp I;
    generateFFTImage(I());
    if (xmipp2Qt(I)) showImage();
}

void ImageViewer::generateFFTImage(matrix2D<double> &out) {
    matrix2D<int> Isubs;
    Isubs.init_zeros(YSIZE(xmippImageFourier),XSIZE(xmippImageFourier));
    out.init_zeros(YSIZE(xmippImageFourier),XSIZE(xmippImageFourier));
    double min_positive;
    bool first=true;
    FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(xmippImageFourier) {
       double ampl, phase, val, eps;
       eps=0;
       val=min_positive;
       switch (fft_show_mode) {
          case 0:
	     ampl=abs(MULTIDIM_ELEM(xmippImageFourier,i));
	     if (ampl!=0)
		val=MULTIDIM_ELEM(out,i)=10*log10(ampl*ampl);
	     else MULTIDIM_ELEM(Isubs,i)=1;
	     break;
	  case 1:
	     val=MULTIDIM_ELEM(out,i)=real(MULTIDIM_ELEM(xmippImageFourier,i));
	     break;
	  case 2:
	     val=MULTIDIM_ELEM(out,i)=imag(MULTIDIM_ELEM(xmippImageFourier,i));
	     break;
	  case 3:
	     val=MULTIDIM_ELEM(out,i)=abs(MULTIDIM_ELEM(xmippImageFourier,i));
	     break;
	  case 4:
	     ampl=abs(MULTIDIM_ELEM(xmippImageFourier,i));
	     val=MULTIDIM_ELEM(out,i)=ampl*ampl;
	     break;
	  case 5:
	     val=MULTIDIM_ELEM(out,i)=arg(MULTIDIM_ELEM(xmippImageFourier,i));
	     break;
       }
       if (val<min_positive || first) {min_positive=val; first=false;}
    }
    
    // Substitute 0s by something a little bit smaller
    if (fft_show_mode==0) {
       FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(Isubs)
	  if (MULTIDIM_ELEM(Isubs,i)==1) MULTIDIM_ELEM(out,i)=min_positive-1;
    }
    
    out.set_Xmipp_origin();
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
       ScrollParam* param_window;	
       param_window = new ScrollParam(0.1, 10, spacing, "Set spacing", "Spacing",
          0, "new window", WDestructiveClose);
       connect( param_window, SIGNAL(new_value(float)), this, SLOT(set_spacing(float)) );
       param_window->setFixedSize(250,150);
       param_window->show();
    } else if (item == sfft) {
       ExclusiveParam* param_window;	
       vector<string> list_values;
       list_values.push_back("10*log10(abs(z)^2)");
       list_values.push_back("real(z)");
       list_values.push_back("imag(z)");
       list_values.push_back("abs(z)");
       list_values.push_back("abs(z)^2");
       list_values.push_back("phase(z)");

       param_window = new ExclusiveParam(list_values, fft_show_mode, "Set FFT show mode", 
          0, "new window", WDestructiveClose);
       connect( param_window, SIGNAL(new_value(int)), this, SLOT(set_fft_show_mode(int)) );
       param_window->setFixedSize(250,200);
       param_window->show();
    } else if (item == ravg) {
       matrix1D<double> radial_profile;
       matrix1D<int> center_of_rot(2), radial_count;
       radial_average(xmippImage(),center_of_rot,radial_profile,radial_count);
       radial_profile.show_with_gnuplot("Radius","Radial average");
       radial_profile.edit();
    } else if (item == line_setup) {
       refineProfileLine();
    } else if (item == profile) {
       matrix1D<double> profile;
       xmippImage().profile(
          xi+STARTINGX(xmippImage()),yi+STARTINGY(xmippImage()),
          xf+STARTINGX(xmippImage()),yf+STARTINGY(xmippImage()),100,profile);
       profile.show_with_gnuplot("Length (%)","Profile");
       profile.edit();
    } else if (item == editctfmodel) {
       FileName fn_param=((FileName)filename).without_extension()+".ctfparam";
       string command=(string)"xmipp_edit -i "+fn_param+" &";
       cout << command << endl;
       system(command.c_str());
    }
}

/* Refine profile line ----------------------------------------------------- */
void ImageViewer::refineProfileLine() {
   vector<float> min, max;
   min.push_back(0); max.push_back(XSIZE(xmippImage())-1);
   min.push_back(0); max.push_back(YSIZE(xmippImage())-1);
   min.push_back(0); max.push_back(XSIZE(xmippImage())-1);
   min.push_back(0); max.push_back(YSIZE(xmippImage())-1);

   vector<float> initial_value;
      initial_value.push_back(xi); initial_value.push_back(yi);
      initial_value.push_back(xf); initial_value.push_back(yf);
   vector<char *> prm_name;
      prm_name.push_back("X initial"); prm_name.push_back("Y initial");
      prm_name.push_back("X final");   prm_name.push_back("Y final");

   ScrollParam* param_window = new ScrollParam(min,max,initial_value,prm_name,
      "Setup profile line", 0, "new window", WDestructiveClose,0);

   // Connect its output to my input (set_spacing)
   connect( param_window, SIGNAL(new_value(vector<float>)),
           this,          SLOT(set_profile_line(vector<float>)) );

   // Show
   param_window->setFixedSize(200,300);
   param_window->show();
   
   // Repaint and draw the old line
   repaint();
   drawLine(xir,yir,xfr,yfr);
}

void ImageViewer::set_profile_line(vector<float> prm) {
   xi=(int)prm[0];
   yi=(int)prm[1];
   xf=(int)prm[2];
   yf=(int)prm[3];
   
   // Remove current drawn line
   drawLine(xir,yir,xfr,yfr);
   
   // Convert coordinates to real coordinates
   int h = height() - status->height();
   xir=ROUND((double) xi*width()/image.width());
   yir=ROUND((double) yi*h      /image.height());
   xfr=ROUND((double) xf*width()/image.width());
   yfr=ROUND((double) yf*h      /image.height());
   
   // Draw the new line
   drawLine(xir,yir,xfr,yfr);
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
    	} else if (w>900 || h>650) {
           double suggested_X=MIN(w,900.0);
           double suggested_Y=MIN(h,650.0);
           double scale=MIN(suggested_X/w,suggested_Y/h);
           w=ROUND(scale*w);
           h=ROUND(scale*h);
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

bool ImageViewer::xmipp2Qt(Image& _image, bool treat_differently_left_right )
{
    bool ok = FALSE;
		
    try { 
      xmippImage = _image;
      // Take the one in showTools
      if (minGray==0 && maxGray==0)
	 ::xmipp2Qt(_image,image,0,255,0,0,treat_differently_left_right); 
      else
	 ::xmipp2Qt(_image,image,0,255,minGray,maxGray,
            treat_differently_left_right);
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

bool ImageViewer::loadImage( const char *fileName,
   double _minGray, double _maxGray, TLoadMode _load_mode ) 
{
    filename = fileName;
    load_mode=_load_mode;
    bool imagic=((string)(filename)).find("imagic:")==0;
    bool ok = FALSE;
    static bool message_shown=false;
    bool treat_separately_left_right=false;
    if ( filename ) {
	
	// try to read image from standard format.	
	if (image.load(filename, 0)) ok = Qt2xmipp(image);	  
	
	if (!ok) {
          try { 
 	    // reads Xmipp Image
            Image tmpImage;
	    if (!imagic) wait_until_stable_size(filename);
            if (imagic) {
               isFourierImage=false;
	       Image *p = Image::LoadImage(filename);
	       if (!p) REPORT_ERROR(1,"ImageViewer::loadImage: Unknown format");
               tmpImage() = (*p)();
               delete p;
            } else if (Is_ImageXmipp(filename)) {
               isFourierImage=false;
               ImageXmipp p;
               p.read((FileName)filename,FALSE,FALSE,apply_geo,FALSE);
               if (load_mode==ImageViewer::PSD_mode) {
                  // It is only the ARMA model
                  xmipp2PSD(p(),p());
               } else if (load_mode==ImageViewer::CTF_mode) {
                  // It is ARMA and CTF together
                  xmipp2CTF(p(),p());
                  treat_separately_left_right=true;
                  options->setItemEnabled(editctfmodel,true);
               }
               tmpImage()=p();
            } else if (Is_FourierImageXmipp(filename)) {
               isFourierImage=true;
	       FourierImageXmipp If; If.read(filename);
               xmippImageFourier=If();
               CenterFFT(xmippImageFourier,true);
	       generateFFTImage(tmpImage());
	       options->setItemEnabled(sfft,true);
	    } else REPORT_ERROR(1,"ImageViewer::loadImage: Unknown format");
            
	    tmpImage().set_Xmipp_origin();
	    minGray=_minGray;
	    maxGray=_maxGray;
	    ok = xmipp2Qt(tmpImage,treat_separately_left_right);
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
	int nx = ROUND((double)e->x() * image.width()  / width());
	int ny = ROUND((double)e->y() * image.height() / h);
        x = CLIP(nx,0,image.width());
        y = CLIP(ny,0,image.height());
        return true;
    }
    return false;
}

/****************************************************/

/*
  Mouse press events.
*/

void ImageViewer::mousePressEvent( QMouseEvent *e ) {
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
	      old_xfr = xir;
	      old_yfr = yir;  
    	      down = true;  
	   }
    	}
    }      
}

/****************************************************/
/*  Mouse release events. */
void ImageViewer::mouseReleaseEvent( QMouseEvent *e ) {
    if (down) {
  	if (convertEvent(e,xf,yf)) {
	   xfr = e->x();   
	   yfr = e->y();
           drawLine(xir,yir, old_xfr, old_yfr);
           drawLine(xir,yir, xfr, yfr);
           old_xfr=xfr;
           old_yfr=yfr;
       	   float distance = sqrt ((double) (xf-xi)*(xf-xi) + (yf-yi)*(yf-yi)); 
	   QString message;
	   message.sprintf("Distance: %.3f Angstroms", distance*spacing);
	   status->setText(message);
       	   down = false;
	}
    }
}

/****************************************************/
/* Record the pixel position of interest. */
void ImageViewer::mouseMoveEvent( QMouseEvent *e ) {
   if (convertEvent(e,pickx,picky)) {
      updateStatus();
      if (down) {
         if ((e->x() < width()) && (e->y() < height()-(status->height()))) {
            xf=pickx;
            yf=picky;
            xfr=e->x();
            yfr=e->y();
            drawLine(xir,yir, old_xfr, old_yfr);
            drawLine(xir,yir, xfr, yfr);
            old_xfr=xfr;
            old_yfr=yfr;

 	    QString message;
            float distance = sqrt ((double) (xf-xi)*(xf-xi) + (yf-yi)*(yf-yi));     
            message.sprintf("Distance: %.3f Angstroms", distance*spacing);
	    status->setText(message);
	 }
      }
   }
}

/****************************************************/
/* Handles key press events. */
void ImageViewer::keyPressEvent( QKeyEvent* e ) {
    switch( e->key() ) {			// Look at the key code
	case Key_F1:
	   menubar->exec();
	   break;
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
             if (e->state() == ControlButton) {	// If 'Ctrol-' key,
                  // Aspect ratio of the original image
                  double ratio_original=(double)xmippImage().ColNo()/xmippImage().RowNo();
                  int new_width=width()/2;
                  int new_height=ROUND(new_width/ratio_original);
    		  resize(new_width, new_height+status->height());
	     }
	     break;
	case Key_P:
	case Key_Plus:				// Double size
             if (e->state() == ControlButton) {	// If 'Ctrol+' key, 
                  double ratio_original=(double)xmippImage().ColNo()/xmippImage().RowNo();
                  int new_width=width()*2;
                  int new_height=ROUND(new_width/ratio_original);
    		  resize(new_width, new_height+status->height());
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

void ImageViewer::paintEvent( QPaintEvent *e ) {
    if ( pm.size() != QSize( 0, 0 ) ) {		// is an image loaded?
       QPainter painter(this);
       painter.setClipRect(e->rect());
       painter.drawPixmap(0, 0, pmScaled);
   }
}

void ImageViewer::drawLine(int x1, int y1, int x2, int y2) {
   QPainter painter(this);
   QBrush brush( NoBrush );
   QPen myPen(red, 3);
   painter.setPen( myPen );      
   painter.setBrush( brush );
   painter.setRasterOp(XorROP);
   painter.drawLine(x1, y1, x2, y2);
}

/****************************************************/
/*
  Explain anything that might be confusing.
*/
void ImageViewer::giveHelp() {
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

void ImageViewer::about() {
    QMessageBox::about( this, "IV (Image Viewer)",
			"Visualizes an Image in Spider format. \n");
}

void ImageViewer::aboutXmipp() {
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

void ImageViewer::set_fft_show_mode(int _fft_show_mode) {
  fft_show_mode = _fft_show_mode;
  if (isFourierImage) {
     ImageXmipp I;
     generateFFTImage(I());
     if (xmipp2Qt(I)) showImage();
  }
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
      loadImage(filename,0,0,load_mode);
      repaint();
      updateStatus();
   }
}
