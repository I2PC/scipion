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
    typedef enum {
       Normal_mode,
       PSD_mode,
       CTF_mode
    } TLoadMode;

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
       double _minGray=0, double _maxGray=0,
       TLoadMode load_mode=ImageViewer::Normal_mode );
    
    /** Set this image as a Fourier image */
    void set_Fourier_flag() {isFourierImage=true;}

    /** Flag whether or not to apply header transformation **/ 
    bool apply_geo;
 
protected:
    virtual void paintEvent( QPaintEvent * );
    void	resizeEvent( QResizeEvent * );
    void	mousePressEvent( QMouseEvent * );
    void	mouseReleaseEvent( QMouseEvent * );
    void	mouseMoveEvent( QMouseEvent * );
    void 	keyPressEvent( QKeyEvent* );
protected:
    void	scale();
    int		alloc_context;
    bool	convertEvent( QMouseEvent* e, int& x, int& y );
    const char* filename;
    Image	xmippImage;		// Xmipp Image
    matrix2D< complex<double> > xmippImageFourier; // Fourier image
    double      minGray;                // Minimum value of the image
    double      maxGray;                // Maximum value of the image
    QImage	image;			// the loaded image
    QPixmap	pm;			// the converted pixmap
    QPixmap	pmScaled;		// the scaled pixmap
    int 	xmippFlag;		// stands for:
    					// -1: No Xmipp image
					//  0: Xmipp image
    bool        check_file_change;
    int         fft_show_mode;          // 0 10*log10(abs(z)^2)
                                        // 1 real(z)
					// 2 imag(z)
					// 3 abs(z)
					// 4 abs(z)^2
					// 5 phase(z)
    bool        isFourierImage;         // True if the image is in Fourier
    TLoadMode   load_mode;
    time_t      modification_time;      // of the file
    QPopupMenu *menubar;
    QPopupMenu *file;
    QPopupMenu *options;
    QPopupMenu *saveimage;
    QWidget    *helpmsg;
    QPrinter   *printer;    
    QLabel     *status;
    QTimer     *timer;
    int 	ss, si, pi, ravg, profile, sfft, line_setup, editctfmodel;
    void	Init();
    bool 	xmipp2Qt(Image& _image, bool treat_differently_left_right=false);
    bool 	Qt2xmipp(QImage &_image);
    bool 	showImage();
    void 	generateFFTImage(matrix2D<double> &out);
    void	updateStatus();
    bool 	reconvertImage();
    void        refineProfileLine();
    void        drawLine(int x1, int y1, int x2, int y2);
    int		pickx, picky;
    bool	may_be_other;
    static 	ImageViewer* other;
    bool 	down;
    int 	xi, yi, xf, yf;
    int 	xir, yir, xfr, yfr, old_xfr, old_yfr; 
    float       spacing;    
    QMultiLineEdit *ed1;

protected slots:
    void	newWindow();
    void	openFile();
    void	saveImage(int);
    void	giveHelp();
    void	doOption(int);
    void 	printIt();
    void 	set_spacing(float _spacing);
    void 	set_profile_line(vector<float> prm);
    void 	set_fft_show_mode(int _fft_show_mode);
    void        about();
    void        aboutXmipp();
    void        check_file();
};

//@}
#endif // SHOWIMG_H
