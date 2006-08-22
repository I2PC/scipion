/***************************************************************************
 *
 * Author:     Javier Rodríguez Fernández (javrodri@gmail.com)
 *             Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Universidad San Pablo CEU (Monteprincipe, Madrid)
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

#include "../showPlotter.hh"
#include <qpainter.h>
#include <qfiledialog.h>
#include <qstyle.h>
#include <qstatusbar.h>
#include <qmime.h>
#include <qapplication.h>
#include <math.h>

/* Empty constructor ------------------------------------------------------- */
PlotSettings::PlotSettings() {
    minX = 0.0; 
    maxX = 0.2; 
    numXTicks = 5;
    
    minY = 0.0;  
    maxY = 0.2;  
    numYTicks = 5;
}

/* Scrolling --------------------------------------------------------------- */
// Increments min, max, by the interval between two ticks times a given number
// to implement scrolling in Plotter::keyPressEvent()
void PlotSettings::scroll(int dx,int dy) {
   double stepX = spanX() / numXTicks;
   minX += dx*stepX;
   maxX += dx*stepX;
   
   double stepY = spanY() / numYTicks;
   minY += dy*stepY;
   maxY += dy*stepY;
}

/* Adjust ------------------------------------------------------------------ */
// called from mouseReleaseEvent, to round the min & max, values to 'nice'
// values and to determinate the number of ticks apropiate for each axis
void PlotSettings::adjust() {
    adjustAxis(minX,maxX);
    adjustAxis(minY,maxY);
}

void PlotSettings::adjustAxis(double &min,double &max) {
    const int MinTicks = 4;
    // a kind of maximum for the step value,
    double grossStep = (max-min) / MinTicks;
    // this to take the decimal logarithm of the gross step, to be nice for
    double step = pow(10,floor(log10(grossStep)));
    // the user, for ex. gross Step = 236, log 236 = 2.37291...,
    // we round it down to 2, and obtain 10²=100, as the candidate.
    // and calculate the other two candidates ex. 5*100=500 and 2*100=200
    if      (5*step < grossStep)  step *=5;
    else if (2*step < grossStep)  step *=2;
    
    // Obtained by rounding the original min down to the nearest multiple
    // of the step
    min = floor(min/step) * step;
    // rounding up to the nearest multiple of the step. 
    max = ceil(max/step) * step; 
}

/*****************************************************************************/
/* Empty constructor ------------------------------------------------------- */
Plotter::Plotter(QWidget *parent, const char *name) :QMainWindow(parent,name) {
   // To use white component to provide easy printing
   setBackgroundMode(PaletteBase); 

   // It makes the layout manager responsible for the widget willing to
   // grow or shrink
   setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding);
   
   // Makes the widget accept focus by clicking or by pressing Tab
   setFocusPolicy(StrongFocus);

   // Initially there is no rubber band
   rubberBandIsShown= false;
    
   // Add XmippGraphics/images to the MIME Path
   QMimeSourceFactory::defaultFactory()->addFilePath(
      xmippBaseDir()+"/Lib/XmippGraphics/images");

   // Initialize save button
   saveButton = new QToolButton (this);
   saveButton->setIconSet(QPixmap::fromMimeSource("filesave20.png"));
   saveButton->adjustSize();
   saveButton->setFixedSize(15,15);
   connect( saveButton, SIGNAL( clicked() ),this, SLOT( saveToFile() ) );
      
   //Buttons to zoomin&out 
   zoomInButton = new QToolButton(this);
   zoomInButton->setIconSet(QPixmap::fromMimeSource("zoomin20.png"));
   zoomInButton->adjustSize();
   zoomInButton->setFixedSize(20,20);
   zoomInButton->hide(); 
   connect(zoomInButton, SIGNAL(clicked()),this, SLOT(zoomIn()));
   
   zoomOutButton = new QToolButton(this);
   zoomOutButton ->setIconSet(QPixmap::fromMimeSource("zoomout20.png"));
   zoomOutButton->adjustSize();
   zoomOutButton->setFixedSize(20,20);
   zoomOutButton->hide(); 
   connect(zoomOutButton, SIGNAL(clicked()),this, SLOT(zoomOut()));

   //Number of posible zooms
   zoomStack.resize(1);
   curZoom=0; 
   zoomStack[curZoom] = PlotSettings(); 

   // Create status bar
   createStatusBar();

   // To catch mouse move
   setMouseTracking(true);

   // Set title
   setCaption(name);
   setGeometry(290,0,250,250);

   show();
}

/* Destructor -------------------------------------------------------------- */
Plotter::~Plotter() {
   delete locationLabelX;
   delete locationLabelY;
   delete formulaLabel;
   delete modLabel;
   delete zoomInButton;
   delete zoomOutButton;
   delete saveButton;
   curveMap.clear();
   curveActive.clear();
}

/* Zoom -------------------------------------------------------------------- */
/* The zoom stack is represented by two member variables
   zoomStack holds the different zoom settings as vector<PlotSettings>
   curZoom holds the current PlotSettings index in the zoomstack
   the zoomIN & OUT buttons are hidden until we call show on them in
   the zoomIn or zoomOut slots
*/
//Enables the zoomOut button depending if the graph can be
// zoomed out any more or not.
void Plotter::zoomOut() {
    if (curZoom>0) {
       --curZoom;
       zoomOutButton->setEnabled (curZoom >0); 
       zoomInButton->setEnabled(true);
       zoomInButton->show();
       refreshCurves();
    }
}

//Enables the zoomIn button depending if the graph can be zoomed in
// any more or not.
void Plotter::zoomIn() {
    if (curZoom<1) {
       ++curZoom;
       zoomInButton->setEnabled(curZoom == 0); 
       zoomOutButton->setEnabled(true);
       zoomOutButton->show();
       refreshCurves();
    }
}

/* Set curve --------------------------------------------------------------- */
void Plotter::setCurveData(int id, const matrix1D<double> &Y) {
   matrix2D<double> data(XSIZE(Y),2);
   FOR_ALL_ELEMENTS_IN_MATRIX1D(Y) {
      data(i,0)=i;
      data(i,1)=Y(i);
   }
   setCurveData(id,data);
}

void Plotter::setCurveData(int id, const matrix1D<double> &X,
   const matrix1D<double> &Y) {
   if (XSIZE(X)!=XSIZE(Y))
      REPORT_ERROR(1,"Plotter::setCurveData: X and Y have different sizes");
   matrix2D<double> data(XSIZE(Y),2);
   FOR_ALL_ELEMENTS_IN_MATRIX1D(Y) {
      data(i,0)=X(i);
      data(i,1)=Y(i);
   }
   setCurveData(id,data);
}

void Plotter::setCurveData(int id, const matrix2D<double> &data) {
    curveMap[id]=data;   
    curveActive[id]=true;
    
    // Determine the bounding box of this data
    double data_minX, data_maxX, data_minY, data_maxY;
    data_minX=data_maxX=data(0,0);
    data_minY=data_maxY=data(0,1);
    for (int i=1; i<YSIZE(data); i++) {
       data_minX=MIN(data_minX,data(i,0));
       data_maxX=MAX(data_minX,data(i,0));
       data_minY=MIN(data_minY,data(i,1));
       data_maxY=MAX(data_minY,data(i,1));
    }

    if (id==0) {
       zoomStack[curZoom].minX=data_minX;
       zoomStack[curZoom].maxX=data_maxX;
       zoomStack[curZoom].minY=data_minY;
       zoomStack[curZoom].maxY=data_maxY;
    } else {
       zoomStack[curZoom].minX=MIN(zoomStack[curZoom].minX,data_minX);
       zoomStack[curZoom].maxX=MAX(zoomStack[curZoom].maxX,data_maxX);
       zoomStack[curZoom].minY=MIN(zoomStack[curZoom].minY,data_minY);
       zoomStack[curZoom].maxY=MAX(zoomStack[curZoom].maxY,data_maxY);
    }
    zoomStack[curZoom].adjust();
    
    refreshCurves();
}

/* Delete Curve ------------------------------------------------------------ */
void Plotter::deleteCurve(int id) {
    curveMap.erase(id);
    curveActive.erase(id);
    refreshCurves();
}

/* Size hints -------------------------------------------------------------- */
QSize Plotter::minimumSizeHint() const {
    return QSize(4 * Margin, 4 * Margin);
}

QSize Plotter::sizeHint() const {
    // We return an ideal size proportion 4:3 aspect ratio
    return QSize(8 * Margin, 6 * Margin);
}

/* Paint event ------------------------------------------------------------- */
void Plotter::paintEvent(QPaintEvent *event) {
    // rect() returns an array of QRects that define the region to repaint
    QMemArray <QRect> rects = event ->region().rects();
    for (int i=0; i<(int)rects.size(); i++)
       // Copy each rectangular area from the pixmap to the widget
       bitBlt(this, rects[i].topLeft(),&pixmap,rects[i]);
       // (dest, destPos, source widget, sourceRect->
       // rectangle in the source that should be copied)
       
    QPainter painter(this);
    
    if (rubberBandIsShown) {
       // Ensure good contrast with the white background.
       painter.setPen(colorGroup().dark());
       painter.drawRect(rubberBandRect.normalize());
    }
    
    if (hasFocus())
       style().drawPrimitive (QStyle::PE_FocusRect, &painter, 
          rect(), colorGroup(), QStyle::Style_FocusAtBorder,
          colorGroup().dark());
}

/* Resize event ------------------------------------------------------------ */
void Plotter::resizeEvent(QResizeEvent *)  {
    int x = width()- (zoomInButton->width() + zoomOutButton->width() +10);
    zoomInButton->move(x,5);
    zoomOutButton->move(x + zoomInButton->width() +5, 5);
    //the buttons are separated by a 5 pixel gap and with 5 pixel offset
    // from the rop and right edges
    refreshCurves();
}

/* Mouse press event ------------------------------------------------------- */
/* Event to control mouse buttons:
   - Left Button  : Rubber band to know the zoom area is showed
   - Right Button : Show the position of the mouse at StatusBar */
void Plotter::mousePressEvent(QMouseEvent *event) {
    if (event->button() == LeftButton) {
       rubberBandIsShown = true;
       rubberBandRect.setTopLeft(event->pos());
       rubberBandRect.setBottomRight(event->pos());
       updateRubberBandRegion();
       setCursor(crossCursor); // and we change the cursor to a crosshair
    }
    
    if (event->button() == RightButton) {
       setCursor(crossCursor); // and we change the cursor to a crosshair
       rubberBandRect.setTopLeft(event->pos());// in the cursor position
       rubberBandRect.setBottomRight(event->pos());
       updateRubberBandRegion();       
    }
}

/* Mouse move event -------------------------------------------------------- */
/* Event to perform zoom if left button is pressed, or update StatusBar if
   right button is pressed. */
void Plotter::mouseMoveEvent(QMouseEvent *event) {
   updateRubberBandRegion(); //we repaint the rubber band in the new position
   rubberBandRect.setBottomRight(event->pos());
   updateRubberBandRegion();
   
   QRect rect = rubberBandRect.normalize();
   rect.moveBy(-Margin, -Margin);

   PlotSettings prevPosSettings = zoomStack [curZoom];
   double dx = prevPosSettings.spanX() / (width()  -2 * Margin);
   double dy = prevPosSettings.spanY() / (height() -2 * Margin);
   PlotSettings Possettings;
   Possettings.minX  = prevPosSettings.minX + dx * rect.left();
   Possettings.maxX  = prevPosSettings.minX + dx * rect.right();
   Possettings.minY  = prevPosSettings.maxY - dy * rect.bottom();
   Possettings.maxY  = prevPosSettings.maxY - dy * rect.top();
   Possettings.adjust();
   Possettings.posX = (prevPosSettings.minX) + rect.left()   * dx;
   Possettings.posY = (prevPosSettings.maxY) - rect.bottom() * dy;
       
   updateCellIndicators(Possettings.posX,Possettings.posY);
   
   //When the user moves the cursor while holding the left button
   if (event->state() & LeftButton) {
      // Repaint the rubber band in the new position
      updateRubberBandRegion(); 
      rubberBandRect.setBottomRight(event->pos());
      updateRubberBandRegion();
   }
   if (event->state() & RightButton) {
      rubberBandRect.setBottomRight(event->pos());
      updateRubberBandRegion();
   }
}

/* Mouse move status ------------------------------------------------------- */
//To refresh the status bar
void Plotter::mouseMoveStatus(QMouseEvent *event) {
  if (event->state())
      updateCellIndicators(event->x(),event->y());
}

/*  Update Status ---------------------------------------------------------- */
void Plotter::updateCellIndicators(double x, double y) {
   locationLabelX->setText(QString::number(x,'f',3));
   locationLabelY->setText(QString::number(y,'f',3));
   statusBar()->message("");
}

// To show messages in StatusBar
void Plotter::updateModLabel(const QString & message) {
   modLabel->setText(message);
}

/* Mouse release event ----------------------------------------------------- */
void Plotter::mouseReleaseEvent(QMouseEvent *event) {
   if (event->button() == LeftButton) {  
      //When the user releases the button we erase the rubberband
      rubberBandIsShown = false;
      updateRubberBandRegion();
      //Restore the standard arrow cursor
      unsetCursor();
       
      QRect rect = rubberBandRect.normalize();
      /* If the rubberband is at least 4x4, we perform the zoom, 
          If it is smaller, it's likely that the user clicked the 
 .       widget by mistake or to give focus, so we do nothing */
      if (rect.width()<4 || rect.height()<4) return;   

      rect.moveBy(-Margin, -Margin);
      PlotSettings prevSettings = zoomStack [curZoom];
      PlotSettings settings;
      double dx = prevSettings.spanX() / (width() -2 * Margin);
      double dy = prevSettings.spanY() / (height() -2 * Margin);
      settings.minX  = prevSettings.minX + dx * rect.left();
      settings.maxX  = prevSettings.minX + dx * rect.right();
      settings.minY  = prevSettings.maxY - dy * rect.bottom();
      settings.maxY  = prevSettings.maxY - dy * rect.top();
      settings.adjust();
       
      zoomStack.resize(curZoom + 1);
      zoomStack.push_back(settings);
      zoomIn();
   }

   if (event->button() == RightButton){
      QRect rect = rubberBandRect.normalize();
      rect.moveBy(-Margin, -Margin);
       
      PlotSettings prevPosSettings = zoomStack [curZoom];
      PlotSettings Possettings;
      double dx = prevPosSettings.spanX() / (width() -2 * Margin);
      double dy = prevPosSettings.spanY() / (height() -2 * Margin);
      Possettings.minX  = prevPosSettings.minX + dx * rect.left();
      Possettings.maxX  = prevPosSettings.minX + dx * rect.right();
      Possettings.minY  = prevPosSettings.maxY - dy * rect.bottom();
      Possettings.maxY  = prevPosSettings.maxY - dy * rect.top();
      Possettings.adjust();
      Possettings.posX  = prevPosSettings.minX + rect.left()   * dx;
      Possettings.posY  = prevPosSettings.maxY - rect.bottom() * dy;
       
      updateCellIndicators(Possettings.posX,Possettings.posY);
      unsetCursor();
   }
}

/* Key press event --------------------------------------------------------- */
/* Event to control Key Pressing:
   - Plus  (+)    : To zoomIn
   - Minus (-)    : To zoomOut
   - Left  Arrow  : Move left
   - Right Arrow  : Move right
   - Up    Arrow  : Move up
   - Down  Arrow  : Move down */
void Plotter::keyPressEvent(QKeyEvent *event) {
   PlotSettings Settings = zoomStack [curZoom];       
   switch ( event->key()) {
      case Key_Q:
         if (event->state() == ControlButton) // If 'Ctrol Q' key, 
  	    exit(0); // Terminate program
         break;
      case Key_Plus:
         zoomIn();
         break;
      case Key_Minus:
         zoomOut();
         break;
      case Key_Left:
         if (Settings.minX>0.01) {
            zoomStack[curZoom].scroll(-1,0);
            refreshCurves();
         }
         break;
      case Key_Right:
         zoomStack[curZoom].scroll(+1,0);
         refreshCurves();
         break;
      case Key_Down:
         if (Settings.minY>0.01) {
            zoomStack[curZoom].scroll(0,-1);
            refreshCurves();
         }
         break;
      case Key_Up:
         zoomStack[curZoom].scroll(0,+1);
         refreshCurves();
         break;
      default:
         QWidget::keyPressEvent(event);
   }
}

/* Wheel event ------------------------------------------------------------- */
/* Event to control mouse wheel to perform zoom */
void Plotter::wheelEvent(QWheelEvent * event) {
   //delta() returns the distance the wheel was rotated in eights of degree
   int numDegrees = event->delta() / 8; 
   int numTicks = numDegrees / 15; //mice tipically work in steps of 15 degrees
    
   if (event->orientation() == Horizontal)
      zoomStack[curZoom].scroll(numTicks,0);
   else
      zoomStack[curZoom].scroll(0,numTicks);
   refreshCurves();
}

/* Update Rubber band (zoom) region ---------------------------------------- */
void Plotter::updateRubberBandRegion() {
    QRect rect = rubberBandRect.normalize();
    update(rect.left(),rect.top(),rect.width(),1);
    update(rect.left(),rect.top(),1, rect.height());
    update(rect.left(),rect.bottom(),rect.width(),1);
    update(rect.right(),rect.top(),1,rect.height());
}

/* Refresh pixmap ---------------------------------------------------------- */
void Plotter::refreshPixmap() {
    pixmap.resize(size()); // Resize the pixmap to the same size as the widget 
    pixmap.fill(this,0,0); // Fill it with the widget erase color
    QPainter painter(&pixmap,this); //To draw on the pixmap
    drawGrid(&painter); // to perform the drawing
    drawCurves(&painter); // to perform the drawing
    update(); // To schedule a paint event for the whole widget
}

/* Draw grid --------------------------------------------------------------- */
void Plotter::drawGrid(QPainter *painter) {
    QRect rect(Margin, Margin , width()-2 *Margin, height() -2 *Margin);
    PlotSettings settings = zoomStack[curZoom];
    QPen quiteDark;
    QPen light;
    light.setColor(white);
    quiteDark.setColor(black);
    
    for(int i =0;i<=settings.numXTicks;i++) {
        int x = rect.left() + (i * (rect.width()-1) / settings.numXTicks);
        
        painter->setPen(quiteDark);
        painter->drawLine(x,rect.top(),x,rect.bottom());
        painter->drawLine(x,rect.bottom(),x,rect.bottom() + 5);

        //it is to draw the numbers corresponding to the tick mark
        double label = settings.minX + (i*settings.spanX() / settings.numXTicks);
        painter->drawText(x - 50,rect.bottom() + 5, 100, 15,
                          AlignHCenter | AlignTop,
                          QString::number(label,'f',2));

        /* COSS: This is showing the inverse of the X axis. Valid only for
                 Fourier plots*/
        /*
        double uplabel; 
        if (label > 0.001) uplabel = 1/label;
        else uplabel = 0;

        //To paint the values uptop
        painter->drawText(x - 50,rect.top() - 15, 100, 15,
                          AlignHCenter | AlignTop,
                        QString::number(uplabel,'f',1));
        */
    }

    for (int j=0; j<=settings.numYTicks;++j) {
        int y = rect.bottom() - (j * (rect.height()-1) / settings.numYTicks);
        double label = settings.minY + (j *settings.spanY() / settings.numYTicks);
        painter->setPen(quiteDark);
        painter->drawLine(rect.left(),y,rect.right(),y);
        painter->drawLine(rect.left()-5,y,rect.left(),y);
        painter->drawText(rect.left()-Margin,y-10,Margin-5,20,AlignRight |AlignVCenter, QString::number(label,'f',2));
     }
     painter->drawRect(rect);
}

/* Draw curves ------------------------------------------------------------- */
void Plotter::drawCurves(QPainter *painter) {
   QPen pen1 (red,   1, SolidLine);
   QPen pen2 (blue , 1, DotLine);
   QPen pen3 (green, 1, DashLine); 
   QPen pen4 (black, 1, DashDotLine); 

   QPen penForIds[4] = {pen1,pen2,pen3,pen4};
   PlotSettings settings = zoomStack[curZoom];

   // Set the clip region as a rectangle that 
   // contains the curves, so it will ignore drawing operations outside
   // the area
   QRect rect(Margin ,Margin , width()-2 *Margin, height() -2*Margin);
   painter->setClipRect(rect.x() +1,rect.y()+1,rect.width() -2, rect.height()-2); 
   
   map <int, matrix2D<double> > ::const_iterator it = curveMap.begin();
   while(it != curveMap.end()) {
      // The first member of the it value gives us the ID
      int id = (*it).first;
      if (curveActive[id]) {
         const matrix2D<double> &data = (*it).second;
         int numPoints = 0;
         int maxPoints = YSIZE(data);
         QPointArray points(maxPoints);

         for (int i=0; i<maxPoints; ++i) {
            // Convert a coordinate pair from plotter coordinate to widget
            // coordinate and stores them in points variable
            // If the users zooms in a lot, we could easily end up with
            // number that cannot be represented as 16-bit signal integers
            double dx = data(i,0) - settings.minX;
            double dy = data(i,1) - settings.minY;
            double x = rect.left()+(dx*(rect.width()-1)/settings.spanX());
            double y = rect.bottom()-(dy*(rect.height()-1)/settings.spanY());
            if (fabs(x)<32768 && fabs(y)<32768)
               points[numPoints++] = QPoint((int)x,(int)y);
         }
         
         points.truncate(numPoints);
         painter->setPen(penForIds[(uint) id % 4]);
         painter->drawPolyline(points);
   }

   // Proceed with the next curve
   ++it;
   }
}

void Plotter::refreshCurves() {
   // Resize the pixmap to have the same size as the widget 
   pixmap.resize(size());
   pixmap.fill(this,0,0); //and fill it with the widget erase color
   QPainter painter(&pixmap,this); //To draw on the pixmap
   drawGrid(&painter); // to perform the drawing
   drawCurves(&painter); //to perform the drawing
   update(); // to schedule a paint event for the whole widget
}

void Plotter::createStatusBar() {
   locationLabelX = new QLabel("X Value",this);
   locationLabelX->setAlignment(AlignHCenter);
   locationLabelX->setMinimumSize(locationLabelX->sizeHint());
   
   locationLabelY = new QLabel("Y Value",this);
   locationLabelY->setAlignment(AlignHCenter);
   locationLabelY->setMinimumSize(locationLabelY->sizeHint());
   
   statusBar()->addWidget(locationLabelX);
   statusBar()->addWidget(locationLabelY);
   statusBar()->addWidget(saveButton);
}

void Plotter::saveToFile() {       
   int lastSlashPos;
   int length;
   FileName file;
   
   QString fileName = QFileDialog::getSaveFileName("","*.png",this);
   if (!fileName.isEmpty()) {
       length = fileName.length();
       //We get the positon of the last Slash
       lastSlashPos = fileName.findRev('/');
       fileName = fileName.right(length-lastSlashPos-1);
       pixmap.save(fileName,"png");
    }
}

// Set plot settings -------------------------------------------------------
void Plotter::setPlotSettings(const PlotSettings &new_settings) {
   zoomStack.clear();
   curZoom=0;
   zoomStack.push_back(new_settings);
   update();
}

// Copy --------------------------------------------------------------------
void Plotter::copy(const Plotter &plotter) {
   curveMap=plotter.curveMap;
   curveActive=plotter.curveActive;
   setPlotSettings(plotter.zoomStack[0]);
}
