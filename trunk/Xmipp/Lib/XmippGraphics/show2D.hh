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

#ifndef SHOWIMG_H
#define SHOWIMG_H

#include <qwidget.h>
#include <qmultilinedit.h>
#include <qimage.h>
#include <qlabel.h>
#include <qprinter.h>
#include <qkeycode.h>
#include <qapplication.h>
#include <qtimer.h>
#include <XmippData/xmippImages.hh>
#include <sys/stat.h>
#include <XmippData/xmippFFT.hh>

/**@name Show 2D*/
//@{

/** ImageViewer.
    This class is in the one in charge of opening a window with the 2D image.
    
    Example of use:
    \begin{verbatim}
       ImageViewer *showimg = new ImageViewer(argv[i], poll);
       showimg->loadImage( argv[i] );
       showimg->show();
    \end{verbatim}
*/
class ImageViewer : public QWidget
{
    Q_OBJECT
public:
    /** Empty constructor. */
    ImageViewer( const char *name=0, bool _check_file_change=false);
    
    /** Constructor with a pointer to an Image */
    ImageViewer( Image *_image=0, const char *name=0);

    /** Constructor with a pointer to a Fourier Xmipp Image */
    ImageViewer( FourierImageXmipp *_FFTimage=0, const char *name=0);

    /** Constructor with a pointer to QImage */
    ImageViewer( QImage *_image=0, const char *name=0);

    ~ImageViewer();

    /** Load image from file.
        The image may be in any standard format or Xmipp format.
	The gray limits are used for common normalization. Set them to 0
	if you don't want to use this option.*/
    bool loadImage( const char *fileName,
       double _minGray=0, double _maxGray=0 );
    
    /** Set image from matrix */
    void setImage( const matrix2D<double> &img);

    /** Flag whether or not to apply header transformation **/ 
    bool        apply_geo;
 
protected:
    void	paintEvent( QPaintEvent * );
    void	resizeEvent( QResizeEvent * );
    void	mousePressEvent( QMouseEvent * );
    void	mouseReleaseEvent( QMouseEvent * );
    void	mouseMoveEvent( QMouseEvent * );
    void 	keyPressEvent( QKeyEvent* );
private:
    void	scale();
    int		alloc_context;
    bool	convertEvent( QMouseEvent* e, int& x, int& y );
    const char* filename;
    Image	xmippImage;		// Xmipp Image
    double      minGray;                // Minimum value of the image
    double      maxGray;                // Maximum value of the image
    QImage	image;			// the loaded image
    QPixmap	pm;			// the converted pixmap
    QPixmap	pmScaled;		// the scaled pixmap
    int 	xmippFlag;		// stands for:
    					// -1: No Xmipp image
					//  0: Xmipp image
    bool        check_file_change;
    time_t      modification_time;      // of the file
    QPopupMenu *menubar;
    QPopupMenu *file;
    QPopupMenu *options;
    QPopupMenu *saveimage;
    QWidget    *helpmsg;
    QPrinter   *printer;    
    QLabel     *status;
    QTimer     *timer;
    int 	ss, si, pi, ravg;
    void	Init();
    bool 	xmipp2Qt(Image& _image);
    bool 	Qt2xmipp(QImage &_image);
    bool 	showImage();
    void	updateStatus();
    bool 	reconvertImage();
    int		pickx, picky;
    int		clickx, clicky;
    bool	may_be_other;
    static 	ImageViewer* other;
    bool 	down;
    int 	xi, yi, xf, yf;
    int 	xir, yir, xfr, yfr, ox, oy; 
    float       spacing;    
    QMultiLineEdit *ed1;

private slots:
    void	newWindow();
    void	openFile();
    void	saveImage(int);
    void	giveHelp();
    void	doOption(int);
    void 	printIt();
    void 	set_spacing(float _spacing);
    void        about();
    void        aboutXmipp();
    void        check_file();
};

//@}
#endif // SHOWIMG_H
