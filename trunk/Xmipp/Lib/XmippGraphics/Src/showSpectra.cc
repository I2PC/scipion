/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Alberto Pascual (pascual@cnb.uam.es)
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

#include "../showSpectra.hh"
#include <qcolordialog.h>
#include <qfontdialog.h>
#include <qmessagebox.h>
#include "../showTools.hh"

/* Init/Clear data --------------------------------------------------------- */
void ShowSpectra::init() {
    V = NULL;
    offX = 10; offY = 10;
    spacing=3;
    x_tick_off=1;
    backColor = black;
    axisColor = white;
    curveColor = green;
    QFont tmpFont("Clean", 6);
    axisFont = tmpFont; 
    ShowSel::init();
}

void ShowSpectra::clear() {
    if (V != NULL) delete V;
    ShowSel::clear();
}

/* Init with vectors ------------------------------------------------------- */
void ShowSpectra::initWithVectors(int _numRows, int _numCols,
   xmippCTVectors *_V, const char *_title) {
   init();
   V=_V;
   fn="";
   initFromVectors();
   setCaption(_title);
   NumRows = _numRows;
   NumCols = _numCols;
   initTable();
   initRightclickMenubar();
   repaint();
}

/* Read a Spectra ---------------------------------------------------------- */
void ShowSpectra::readFile(const FileName &_fn) _THROW {
    clear();
    fn = _fn;
    setCaption(fn.c_str());
    readDatFile(_fn);
}

void ShowSpectra::readDatFile(const FileName &_fn) _THROW {
    ifstream fh_in(_fn.c_str());
    if (!fh_in)
       REPORT_ERROR(1,(string)"ShowSpectra::readFile: Cannot open"+_fn);
    V=new xmippCTVectors(fh_in);
    fh_in.close();

    annotateTime(_fn);
    initFromVectors();
}

/* Read vectors ------------------------------------------------------------ */
void ShowSpectra::initFromVectors() {
    listSize = V->size();
    if (listSize==0)
       REPORT_ERROR(1,"ShowSpectra::readFile: Input file is empty");
    imgnames        = new FileName[listSize];
    selstatus       = new bool[listSize];
    initContents();

    // Determine min, max and average
    minPixel=MAXFLOAT;
    maxPixel=MINFLOAT;
    for (long i = 0; i < V->size(); i++)
       for (int j = 0; j < V->theItems[0].size(); j++) {
          double val=V->theItems[i][j];
          if (minPixel>val) minPixel=val;
	  if (maxPixel<val) maxPixel=val;
       }
    projXdim=projYdim=100;
    for (long i=0; i<V->size(); i++) {
      	imgnames[i] = V->theTargets[i];
	selstatus[i]= SelLine::ACTIVE;
    }
}

/* Init Rightclick menubar ------------------------------------------------- */
void ShowSpectra::initRightclickMenubar() {    
   menubar = new QPopupMenu(); 
   setFileRightclickMenubar(); 

   // Options .............................................................
   options =  new QPopupMenu();
   setCommonOptionsRightclickMenubar();
   labeltype=Filename_LABEL;
   // spectra common options (pased to spectraSOM too)
   setCommonSpectraOptionsRightclickMenubar();

   // Statistics
   options->insertItem( "View average and SD Images",  this,  SLOT(showStats()));
   options->insertItem( "Show average and SD Spectra", this,  SLOT(showSpectraStats()));
   options->insertSeparator();
    
   // Insert options the menu
   menubar->insertItem( "&Options", options );    
   menubar->insertSeparator();


   // Inser Help and Quit
   insertGeneralItemsInRightclickMenubar();
}

/* Produce pixmap ---------------------------------------------------------- */
void ShowSpectra::paintCell(QPainter *p, int row, int col,const QRect & cr,
   bool selected, const QColorGroup & cg) {
   int scprojXdim=(int)((currScale*projXdim)/100.0);
   int scprojYdim=(int)((currScale*projYdim)/100.0);
   int i=indexOf(row,col);
   if (i >= listSize) return;
   int N=V->theItems[i].size();

   if (indexOf(row,col) >= listSize) return; 
   QPixmap background( columnWidth(col), rowHeight(row) );
   background.fill(backColor);
   p->drawPixmap(0, 0, background);
   p->setFont(axisFont);

   // Get minimum and Maximum of this spectrum
   double myMinValue, myMaxValue;
   if (!options->isItemEnabled(mi_Individual_norm)) {
      // Individual normalization
      myMaxValue = MINFLOAT;
      myMinValue = MAXFLOAT;
      for (int l=0; l<N ; l++) {
	  double val=V->theItems[i][l];
	  if (myMinValue>val) myMinValue=val;
	  if (myMaxValue<val) myMaxValue=val;
      }
   } else {
      // Global normalization
      myMaxValue = maxPixel;
      myMinValue = minPixel;
   }


   // Calculate slope and intersection of linear transformation
   double slope=0;
   if (myMinValue !=myMaxValue) 
      slope=(scprojXdim-2*offX)/(myMaxValue-myMinValue);
  
   // Paint X axis .........................................................
   QPen pen(axisColor);
   p->setPen(pen);
   p->drawLine(offX, scprojYdim-offY, scprojXdim-offX, scprojYdim-offY);
   double scaleX = (scprojXdim-2*offX)/N;
   int    istep = spacing;
   for (int l=x_tick_off-1; l<=N; l+=istep) {
	 int x = offX + (int)(l*scaleX);
	 // Draw legend
         if (!options->isItemEnabled(mi_showXlegend)) {
           QString tmpS;
           tmpS.setNum(l+1);
           p->drawText(x, (int) (scprojYdim-offY+ p->font().pixelSize()), tmpS);
	 }
	 // Draw X grid lines
	 if (!options->isItemEnabled(mi_showgrid))
   	    p->drawLine(x, offY, x, scprojYdim-offY);
	 else
   	    p->drawLine(x, scprojYdim-offY, x,scprojYdim-offY-3);
   }

   // Paint Y axis .........................................................
   p->drawLine(offX, offY, offX, scprojYdim-offY);
   for (int l=0; l<=3; l++) {
      // Draw legend
      if (!options->isItemEnabled(mi_showYlegend)) {
         QString tmpS;
         if (l == 0) tmpS.setNum(myMinValue, 'f', 2);
         else	     tmpS.setNum((l*(myMaxValue/3.0)), 'f', 2);
         p->drawText(2, (int) (scprojYdim-offY-(l*((scprojYdim-2*offY)/3.0))), tmpS);
      }
      
      // Draw Y grid lines
      if (!options->isItemEnabled(mi_showgrid))
         p->drawLine(offX, (int) (scprojYdim-offY-(l*((scprojYdim-2*offY)/3.0))),
	   scprojXdim-offX, (int) (scprojYdim-offY-(l*((scprojYdim-2*offY)/3.0))));
      else
         p->drawLine(offX, (int) (scprojYdim-offY-(l*((scprojYdim-2*offY)/3.0))),
	          2+offX, (int) (scprojYdim-offY-(l*((scprojYdim-2*offY)/3.0))));
   }

   // Paint curves .........................................................
   pen.setColor(curveColor);
   pen.setWidth(3);
   p->setPen( pen );

   int      x= offX;
   double myY= offY + slope*(V->theItems[i][0] - myMinValue);
   int      y= scprojYdim - (int) myY;
   p->moveTo(x, y);
   p->drawPoint(x, y);
   for (int l=1; l<N; l++) {
       x = offX + (int)(l*scaleX);
       myY = offY + slope*(V->theItems[i][l] - myMinValue);
       y = scprojYdim - (int) myY;
       p->lineTo(x, y);
   }
   // Draw Frame and label
   drawFrameAndLabel(p,row,col,i,1);
}

/* Open new file ----------------------------------------------------------- */
void ShowSpectra::openNewFile(const FileName &_fn) {
   init();
   readFile(_fn);
   initTable();
   repaint();
}

/* Change options ----------------------------------------------------------- */
void ShowSpectra::changeGrid()    {changeBoolOption(mi_showgrid, mi_hidegrid);}
void ShowSpectra::changeXlegend() {changeBoolOption(mi_showXlegend, mi_hideXlegend);}
void ShowSpectra::changeYlegend() {changeBoolOption(mi_showYlegend, mi_hideYlegend);}

// Show Spectra Stats ------------------------------------------------------
void ShowSpectra::showSpectraStats() {
   long counter=0;
   // get number of marked spectra
   for (long i = 0; i<listSize; i++) if (cellMarks[i]) counter++; 
   if (counter < 2)
       QMessageBox::about( this, "Error!", "No enough spectra selected\n");
   else {
       xmippCTVectors mySpectra(0, true);
       mySpectra.theItems.resize(counter);
       mySpectra.theTargets.resize(counter);
       long myIndex = 0;
       for (long i=0; i<listSize; i++)
   	 if (cellMarks[i]) {
            mySpectra.theItems[myIndex] = V->theItems[i]; 
            mySpectra.theTargets[myIndex] = V->theTargets[i]; 
            myIndex++;
         } 

       xmippCTVectors *myVect = new xmippCTVectors(0, true);
       *myVect = mySpectra.getStatVector();  
       ShowSpectra *myST = new ShowSpectra;
       myST->initWithVectors(1, 2, myVect, "Average and SD");
       myST->show();
   }
}

// Change colors -----------------------------------------------------------
void ShowSpectra::GUIchangeColor(QColor &_color, const char *_color_title) {
   QColor tmpColor = QColorDialog::getColor(axisColor, this, _color_title);
   if (tmpColor.isValid()) {
      _color = tmpColor;
      repaintContents();
   }
}

void ShowSpectra::changeBackColor()  {GUIchangeColor(backColor, "Back color");}
void ShowSpectra::changeAxisColor()  {GUIchangeColor(axisColor, "Axis color");}
void ShowSpectra::changeCurveColor() {GUIchangeColor(curveColor,"Curve color");}

// Change Font -------------------------------------------------------------
void ShowSpectra::changeFont() {
   QFont tmpFont; bool ok;
   tmpFont = QFontDialog::getFont(&ok, axisFont, this, "Font type");
   if (ok) {
   	axisFont = tmpFont;
        repaintContents();
   }
}
// Change Grid spacing in X axis -------------------------------------
void ShowSpectra::changeXstep() {
      int N=V->theItems[0].size();
      ScrollParam2* param_window;
      param_window = new ScrollParam2(1 //min
                                    ,N //max
				    , spacing //init value
				    , x_tick_off // init value
				    , "Set spacing",
         0, "new window", WDestructiveClose,0);
      connect( param_window, SIGNAL(new_value(float,float)), this,
      SLOT(set_spacing(float,float)) );
      param_window->setFixedSize(250,200);
      param_window->show();

}
/****************************************************/

void ShowSpectra::set_spacing(float _spacing, float _x_tick_off) {
  spacing = (int) _spacing;
  x_tick_off = (int) _x_tick_off;
    clearContents();
    repaintContents();
  
}

void ShowSpectra::setCommonSpectraOptionsRightclickMenubar()
{
   // Show grid
   mi_showgrid = options->insertItem( "Show Grid", this,  SLOT(changeGrid()));
   mi_hidegrid = options->insertItem( "Hide Grid", this,  SLOT(changeGrid()));
   options->setItemEnabled(mi_showgrid, true);
   options->setItemEnabled(mi_hidegrid, false);
   options->insertSeparator();

   // Show X legend
   mi_showXlegend = options->insertItem( "Show X legend", this,  SLOT(changeXlegend()));
   mi_hideXlegend = options->insertItem( "Hide X legend", this,  SLOT(changeXlegend()));
   options->insertItem( "&Change X-axis step...", this, SLOT(changeXstep()));
   options->setItemEnabled(mi_showXlegend, false);
   options->setItemEnabled(mi_hideYlegend, true);
   options->insertSeparator();

   // Show Y legend
   mi_showYlegend = options->insertItem( "Show Y legend", this,  SLOT(changeYlegend()));
   mi_hideYlegend = options->insertItem( "Hide Y legend", this,  SLOT(changeYlegend()));
   options->setItemEnabled(mi_showYlegend, false);
   options->setItemEnabled(mi_hideYlegend, true);
   options->insertSeparator();

   // Colors and change font
   QPopupMenu* colorMenu = new QPopupMenu();
      colorMenu->insertItem( "&Background", this, SLOT(changeBackColor()));
      colorMenu->insertItem( "&Curve", this, SLOT(changeCurveColor()));
      colorMenu->insertItem( "&Axis", this, SLOT(changeAxisColor()));
   menubar->insertItem( "&Colors", colorMenu);
   menubar->insertSeparator();
}
