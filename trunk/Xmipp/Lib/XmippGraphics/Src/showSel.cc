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

#include "../showSel.hh"
#include "../showTools.hh"
#include <XmippInterface/xmippVTK.hh>
#include <qmessagebox.h>
#include <qfiledialog.h>

/* Init/Clear data --------------------------------------------------------- */
void ShowSel::init() {
    labeltype       = SFLabel_LABEL;
    imgnames        = NULL;
    selstatus       = NULL;
    ShowTable::init();
}

void ShowSel::clear() {
    if (selstatus != NULL) delete[] selstatus;
    if (imgnames  != NULL) delete[] imgnames;
    ShowTable::clear();
}

/* Initialize with a sel file.---------------------------------------------- */
void ShowSel::initWithFile( int _numRows, int _numCols,
   const FileName &_fn) {
   init();
   readFile(_fn);
   NumRows = _numRows;
   NumCols = _numCols;
   initTable();
   initRightclickMenubar();
   repaint();
}

void ShowSel::initWithObject(int _numRows, int _numCols,
   SelFile &_SF, const char *_title) {
   init();
   fn="";
   setCaption(_title);
   _SF.go_first_ACTIVE();
   readObject(_SF);
   NumRows = _numRows;
   NumCols = _numCols;
   initTable();
   initRightclickMenubar();
   repaint();
}

/* Read a Selfile ---------------------------------------------------------- */
void ShowSel::readFile(const FileName &_fn) _THROW {
    clear();
    fn              = _fn;
    setCaption(fn.c_str());
    readSelFile(_fn);
}

void ShowSel::readSelFile(const FileName &_fn) _THROW {
    SelFile         SF(_fn);
    annotateTime(_fn);
    readObject(SF);
}

void ShowSel::readObject(SelFile &SF) {
    listSize        = SF.ImgNo()+SF.ImgNo(SelLine::DISCARDED);
    if (listSize==0)
       REPORT_ERROR(1,"ShowSel::readFile: Input selfile is empty");
    imgnames        = new FileName[listSize];
    selstatus       = new bool[listSize];
    initContents();
    SF.ImgSize(projYdim,projXdim);
    int i=0;
    while (!SF.eof()) {
      	imgnames[i] = SF.get_current_file();
	selstatus[i]= SF.Is_ACTIVE();
	SF.next();
	i++;
    }
}

/* Compute global normalization params ------------------------------------- */
void ShowSel::compute_global_normalization_params() {
   bool first=true;
   for (int i=0; i<listSize; i++) {
      ImageXmipp I(imgnames[i]);
      double min_val, max_val;
      I().compute_double_minmax(min_val,max_val);
      if (first || min_val<minPixel) minPixel=min_val;
      if (first || max_val>maxPixel) maxPixel=max_val;
      first=false;
   }
}

/* Init table -------------------------------------------------------------- */
void ShowSel::initTable() {
   ShowTable::initTable();
   setFocusPolicy( StrongFocus ); // keyboard focus is accepted
   // Really set size
   setMaximumSize(maxWidth,maxHeight);    
   resize(maxWidth,maxHeight);
}

/* Init Rightclick menubar ------------------------------------------------- */
void ShowSel::initRightclickMenubar() {    
   menubar = new QPopupMenu(); 
   setFileRightclickMenubar(); 

   // Options .............................................................
   options =  new QPopupMenu();
   setCommonOptionsRightclickMenubar();

   // What kind of labels
   mi_imgAsLabels = options->insertItem( "Show Image Names as Labels", this,  SLOT(changeLabels()));
   mi_selAsLabels = options->insertItem( "Show Sel status as Labels", this,  SLOT(changeLabels()));
   options->setItemEnabled(mi_imgAsLabels, false);
   options->setItemEnabled(mi_selAsLabels, true);
   labeltype=Filename_LABEL;
   options->insertSeparator();

   // Statistics
   options->insertItem( "View average and SD Images", this,  SLOT(showStats()));
   options->insertItem( "Show Sel Statistics", this,  SLOT(showSelStats()));
   options->insertSeparator();
    
   // Form the menu
   menubar->insertItem( "&Options", options );    
   menubar->insertSeparator();
   insertGeneralItemsInRightclickMenubar();
}

void ShowSel::setFileRightclickMenubar() {
  QPopupMenu * file = new QPopupMenu();
      file->insertItem( "Open...", this,  SLOT(GUIopenFile()));
      QPopupMenu * fileSave = new QPopupMenu();
      fileSave->insertItem( "As discarded...",
         this, SLOT(saveSelFileDiscarded()));
      fileSave->insertItem( "As active and the rest as discarded...",
         this, SLOT(saveSelFileActive()));
      fileSave->insertItem( "In a new sel file...",
         this, SLOT(saveSelFileNew()));
      file->insertItem( "&Save Selected Images in a Sel File", fileSave);
   menubar->insertItem( "&File", file );
}

void ShowSel::setCommonOptionsRightclickMenubar() {
   // Normalization
   mi_Individual_norm = options->insertItem( "Individual Normalization",
      this, SLOT(changeNormalize()));
   mi_Global_norm = options->insertItem( "Global Normalization",
      this, SLOT(changeNormalize()));
   options->setItemEnabled(mi_Individual_norm, false);
   options->setItemEnabled(mi_Global_norm, true);
   options->insertSeparator();

   // Show/Hide labels
   mi_showLabel = options->insertItem( "Show Labels", this,  SLOT(changeShowLabels()));
   mi_hideLabel = options->insertItem( "Hide Labels", this,  SLOT(changeShowLabels()));
   options->setItemEnabled(mi_showLabel, false);
   options->setItemEnabled(mi_hideLabel, true);
   options->insertSeparator();

   // Select/unselect
   options->insertItem( "Select All", this,  SLOT(SelectAll()));
   options->insertItem( "Unselect All", this,  SLOT(unSelectAll()));
   options->insertSeparator();
}

/* Cell label -------------------------------------------------------------- */
const char * ShowSel::cellLabel(int i) const {
   if (options->isItemEnabled(mi_showLabel)) return NULL;
   switch (labeltype) {
      case SFLabel_LABEL:
         return (selstatus[i])? "1":"-1";
      case Filename_LABEL:
         return imgnames[i].c_str();
   }
}

/* Produce pixmap ---------------------------------------------------------- */
void ShowSel::producePixmapAt(int i) {
    ImageXmipp I;
    // Read image
    if (Is_ImageXmipp(imgnames[i]))
       // Plain Xmipp images
       I.read(imgnames[i]);
    else if (Is_FourierImageXmipp(imgnames[i])) {
       // FFT Xmipp images: plot log10(1+|I|^2)
       FourierImageXmipp If;
       If.read(imgnames[i]);
       FFT_magnitude(If,I());
       FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(I())
          MULTIDIM_ELEM(I(),i)=
	     log10(1+MULTIDIM_ELEM(I(),i)*MULTIDIM_ELEM(I(),i));
    } else
       // Unknown image
       I().init_zeros(projYdim,projXdim);

    // Scale and normalize
    int minGray, maxGray;
    scale_and_normalize(I(), options->isItemEnabled(mi_Individual_norm),
       minGray, maxGray);

    // Convert Xmipp image to Pixmap
    content[i]=new QPixmap;
    xmipp2Pixmap(I,content[i],minGray,maxGray);
}

/* Grow older all contents ------------------------------------------------- */
void ShowSel::insert_content_in_queue(int i) {
   if (listSize<NumRows*NumCols) return;
   // Check if the image i is already in the queue
   int jmax=content_queue.size();
   bool found=false;
   list<int>::iterator ptr=content_queue.begin();
   for (int j=0; j<jmax; j++) {
      if ((*ptr)==i) {found=true; break;}
      ptr++;
   }
   if (found) content_queue.erase(ptr);
   content_queue.push_back(i);
   
   // If the queue is longer than 3 times the visible area then remove some
   // images
   if (jmax+1>2*NumRows*NumCols) {
      int i_to_remove=content_queue.front();
      content_queue.pop_front();
      if (content[i_to_remove]!=NULL) {
         delete content[i_to_remove];
         content[i_to_remove]=NULL;
      }
   }
}

/* Open new file ----------------------------------------------------------- */
void ShowSel::openNewFile(const FileName &_fn) {
   init();
   readFile(_fn);
   if (options->isItemEnabled(mi_Individual_norm))
      compute_global_normalization_params();
   initTable();
   repaint();
}

/* Save SelFiles ----------------------------------------------------------- */
// This function saves the sel file with the selected images as discarded.
void ShowSel::saveSelFileDiscarded() {
   SelFile SFNew;
   bool saveFile=false;
   for (int i=0; i<listSize; i++) {
      if (cellMarks[i]) {
   	 saveFile=true;
   	 SFNew.insert(imgnames[i],SelLine::DISCARDED);
      } else
   	 if (selstatus[i]) SFNew.insert(imgnames[i],SelLine::ACTIVE);
   	 else		   SFNew.insert(imgnames[i],SelLine::DISCARDED);
   }
   if (saveFile) writeSelFile(SFNew);
   else QMessageBox::about( this, "Error!", "No images selected\n");
}

/* This function saves the sel file with the selected images as active and 
  the rest of the sel file as discarded. */
void ShowSel::saveSelFileActive() {
   SelFile SFNew;
   bool saveFile=false;
   for (int i=0; i<listSize; i++) {
      if (cellMarks[i]) {
   	 saveFile=true;
   	 SFNew.insert(imgnames[i],SelLine::ACTIVE);
      } else
   	 SFNew.insert(imgnames[i],SelLine::DISCARDED);
   }
   if (saveFile) writeSelFile(SFNew);
   else QMessageBox::about( this, "Error!", "No images selected\n");
}

/* This function saves a new sel file with the selected images as active.*/
void ShowSel::saveSelFileNew() {
   SelFile SFNew;
   bool saveFile=false;
   for (int i=0; i<listSize; i++) {
      if (cellMarks[i]) {
   	 saveFile=true;
   	 SFNew.insert(imgnames[i],SelLine::ACTIVE);
      }
   }
   if (saveFile) writeSelFile(SFNew);
   else QMessageBox::about( this, "Error!", "No images selected\n");
}

/* Save a Selfile.
   Make all possible checkings */
void ShowSel::writeSelFile(SelFile &_SF) {
   QString newfilename = QFileDialog::getSaveFileName(
      QString::null, "*.sel", this, "Sel files");
   if (!newfilename.isEmpty() ) {
      QFileInfo fi(newfilename);
      if (fi.extension(false) != "sel") {
   	 if ( QMessageBox::information( this, "Showsel application",
   		"The file has no ""sel"" extension. add it? ", 
   		"Yes", "No") == 0) newfilename += ".sel";
      }
      fi.setFile(newfilename);
      if (fi.exists())
   	 if ( QMessageBox::information( this, "Showsel application",
   		"The file already exist. Overwrite?", 
   		"Yes", "No") == 0) _SF.write((string) newfilename);
   	 else QMessageBox::about( this, "Warning!", "Saving aborted\n");			   
      else _SF.write((string) newfilename);
   } else  QMessageBox::about( this, "Warning!", "Saving aborted\n");
}

// Change options ----------------------------------------------------------
void ShowSel::changeNormalize() {  
    bool indivNorm = options->isItemEnabled(mi_Individual_norm); 
    if (indivNorm) {
       options->setItemEnabled(mi_Individual_norm, false);
       options->setItemEnabled(mi_Global_norm, true);
    } else {
       options->setItemEnabled(mi_Individual_norm, true);
       options->setItemEnabled(mi_Global_norm, false);
       if (minPixel==0 && maxPixel==0)
          compute_global_normalization_params();
    }
    clearContents();
    repaintContents();
}

void ShowSel::changeShowLabels() {changeBoolOption(mi_showLabel, mi_hideLabel);}
void ShowSel::changeLabels()     {changeBoolOption(mi_imgAsLabels, mi_selAsLabels);}

// Show statistics ---------------------------------------------------------
void ShowSel::showStats() {  
    SelFile SFNew;
    for (int i=0; i<listSize; i++)
        if (cellMarks[i])
	  SFNew.insert(imgnames[i], SelLine::ACTIVE);
    if (SFNew.ImgNo()) ShowTable::showStats(SFNew);
    else QMessageBox::about( this, "Error!", "No images selected\n");
}

// Show Sel Stats ----------------------------------------------------------
void ShowSel::showSelStats() {
   int total = 0; 
   int active = 0;
   int discarded = 0;
   int commented = 0;
   SelFile SF(fn);
   while (!SF.eof()) {
      // Get file  
         if (SF.Is_ACTIVE())
	    active++;
	 else if (SF.Is_DISCARDED())
	    discarded++;
	 else if (SF.Is_COMMENT())
	    commented++;   
      SF.next();
      total++;
   }  // while 
   QString tmpS, Str1;
   Str1 = "Sel File Name : ";
   Str1 += fn.c_str();
   Str1 += "\nTotal number of images : ";
   tmpS.setNum(total);
   Str1 += tmpS;
   Str1 += " (100.00%)\n";
   Str1 += "Active images 	  : ";
   tmpS.setNum(active);
   Str1 += tmpS;
   Str1 += " (";
   tmpS.setNum((float) active*100.0/(float) total , 'f', 2);
   Str1 += tmpS;
   Str1 += "%)\n";
   Str1 += "Discarded image: ";
   tmpS.setNum(discarded);
   Str1 += tmpS;
   Str1 += " (";
   tmpS.setNum((float) discarded*100.0/(float) total , 'f', 2);
   Str1 += tmpS;
   Str1 += "%)\n";
   Str1 += "Commented image: ";
   tmpS.setNum(commented);
   Str1 += tmpS;
   Str1 += " (";
   tmpS.setNum((float) commented*100.0/(float) total , 'f', 2);
   Str1 += tmpS;
   Str1 += "%)";
   QMessageBox::about( (QWidget*)this, "Sel File Statistics", Str1);
}

// Unselect/Select all -----------------------------------------------------
void ShowSel::SelectAll() {  
  for (int i = 0; i < listSize; i++) 
     if (!cellMarks[i]) {
       cellMarks[i] = true;
       updateCellIdx(i);
     }
}

void ShowSel::unSelectAll() {  
  for (int i = 0; i < listSize; i++) 
     if (cellMarks[i]) {
       cellMarks[i] = false;
       updateCellIdx(i);
     }
}

// Update status -----------------------------------------------------------
void ShowSel::contentsMouseMoveEvent( QMouseEvent* e) {
    QPoint Pos = e->pos(); // extract pointer position
    int row=rowAt(Pos.y());
    int col=columnAt(Pos.x());
    if (row<0 || col<0) return;
    updateStatus(indexOf(row,col));
}

void ShowSel::updateStatus(int i) {
   if (i>listSize) return;
   status->setText(imgnames[i].c_str());
}
