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
#include "../show2D.hh"
#include "../showTools.hh"
#include "../showAssignCTF.hh"
#include <Reconstruction/Programs/Prog_assign_CTF.hh>
#include <XmippData/xmippArgs.hh>
#include <qmessagebox.h>
#include <qfiledialog.h>

/* Empty constructor ------------------------------------------------------- */
ShowSel::ShowSel(): ShowTable() {
    load_mode       = Normal_mode;
}

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
   const FileName &_fn, double _minGray, double _maxGray) {
   init();
   readFile(_fn, _minGray, _maxGray);
   if (_numRows!=-1 && _numCols!=-1) {
      NumRows = _numRows;
      NumCols = _numCols;
   } else {
      NumCols=FLOOR(900.0/projXdim);
      NumRows=FLOOR(700.0/projYdim);
   }
   initTable();
   initRightclickMenubar();
   repaint();
}

void ShowSel::initWithFile( int _numRows, int _numCols,
   const FileName &_fn, TLoadMode _load_mode) {
   init();
   load_mode=_load_mode;
   NumRows = _numRows;
   NumCols = _numCols;
   readFile(_fn, 0, 0);
   if (_numRows==-1 || _numCols==-1) {
      NumCols=FLOOR(900.0/projXdim);
      NumRows=FLOOR(700.0/projYdim);
   }
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
void ShowSel::readFile(const FileName &_fn, double _minGray, double _maxGray)
    {
    clear();
    fn              = _fn;
    setCaption(fn.c_str());
    readSelFile(_fn, _minGray, _maxGray);
}

void ShowSel::readSelFile(const FileName &_fn, double _minGray, double _maxGray) {
    SelFile         SF(_fn);
    annotateTime(_fn);
    readObject(SF, _minGray, _maxGray);
}

void ShowSel::readObject(SelFile &SF, double _minGray, double _maxGray) {
    listSize        = SF.ImgNo();
    if(!showonlyactive)   listSize += SF.ImgNo(SelLine::DISCARDED);
    if (listSize==0)
       REPORT_ERROR(1,"ShowSel::readFile: Input selfile is empty");
    imgnames        = new FileName[listSize];
    selstatus       = new bool[listSize];
    initContents();
    SF.ImgSize(projYdim,projXdim);
    if (load_mode==PSD_mode && NumRows!=-1 && NumCols!=-1) {
       // Scale to make the images fit into a reasonable window
       double suggested_Xdim=MIN(900.0/NumCols,projXdim);
       double suggested_Ydim=MIN(700.0/NumRows,projYdim);
       double scale_X=suggested_Xdim/projXdim;
       double scale_Y=suggested_Ydim/projYdim;
       double scale=MIN(scale_X, scale_Y);
       projYdim=FLOOR(scale*projYdim);
       projXdim=FLOOR(scale*projXdim);
    }
    minPixel=_minGray; maxPixel=_maxGray;
    int i=0;
    while (!SF.eof()) {
        if(SF.Is_ACTIVE() || !showonlyactive) {
      	   imgnames[i] = SF.get_current_file();
	   selstatus[i]= SF.Is_ACTIVE();
	   i++;
        }
	SF.next();
    }
}

/* Compute global normalization params ------------------------------------- */
void ShowSel::compute_global_normalization_params() {
   bool first=true;
   for (int i=0; i<listSize; i++) {
      ImageXmipp I(imgnames[i]);
      if      (load_mode==PSD_mode) xmipp2PSD(I(),I());
      else if (load_mode==CTF_mode) xmipp2CTF(I(),I());
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
   
   // Show this image
   options->insertItem( "Reload all", this, SLOT(reloadAll()));
   options->insertItem( "Show this image separately", this, SLOT(showThisImage()));
   options->insertSeparator();
   
   // Recalculate CTF model
   if (load_mode==CTF_mode) {
      options->insertItem( "Edit CTF model", this, SLOT(editCTFmodel()));
      options->insertItem( "Recompute CTF model", this, SLOT(recomputeCTFmodel()));
      options->insertSeparator();
      
      options->setItemEnabled(mi_imgAsLabels,false);
      options->setItemEnabled(mi_selAsLabels,false);
   }
   
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
   if (load_mode==Normal_mode || load_mode==CTF_mode) {
      options->setItemEnabled(mi_showLabel, false);
      options->setItemEnabled(mi_hideLabel, true);
   } else if (load_mode==PSD_mode) {
      options->setItemEnabled(mi_showLabel, true);
      options->setItemEnabled(mi_hideLabel, false);
   }
   options->insertSeparator();

   // Select/unselect
   options->insertItem( "Select All", this,  SLOT(SelectAll()));
   options->insertItem( "Unselect All", this,  SLOT(unSelectAll()));
   options->insertSeparator();
}

/* Cell label -------------------------------------------------------------- */
const char * ShowSel::cellLabel(int i) const {
   if (options->isItemEnabled(mi_showLabel)) return NULL;
   if (load_mode==CTF_mode) {
      // Get the defocus parameters from the ctfparam file
      FileName fn_param=imgnames[i].without_extension()+".ctfparam";
      XmippCTF ctf; ctf.read(fn_param,false);
      string defocus_val=ItoA(ROUND(MIN(ctf.DeltafU,ctf.DeltafV)),6)+" "+
                         ItoA(ROUND(MAX(ctf.DeltafU,ctf.DeltafV)),6)+" "+
                         ItoA(ABS(ROUND(ctf.DeltafU-ctf.DeltafV)));
      return defocus_val.c_str();
   }
   else
      switch (labeltype) {
         case SFLabel_LABEL:
            return (selstatus[i])? "1":"-1";
         case Filename_LABEL:
            return imgnames[i].c_str();
      }
}

/* Produce pixmap ---------------------------------------------------------- */
void ShowSel::producePixmapAt(int i) {
    bool treat_separately_left_right=false;
    ImageXmipp I;
    // Read image
    if (imgnames[i].find("imagic:")!=-1) {
       Image *img = Image::LoadImage(imgnames[i]);
       I()=(*img)();
       delete img;
    } else if (Is_ImageXmipp(imgnames[i])) {
       // Plain Xmipp images
       I.read(imgnames[i],FALSE,FALSE,apply_geo,FALSE);
       if      (load_mode==PSD_mode) xmipp2PSD(I(),I());
       else if (load_mode==CTF_mode) {
          treat_separately_left_right=true;
          xmipp2CTF(I(),I());
       }
    } else if (Is_FourierImageXmipp(imgnames[i])) {
       // FFT Xmipp images: plot log10(1+|I|^2)
       FourierImageXmipp If;
       If.read(imgnames[i]);
       FFT_magnitude(If(),I());
       FOR_ALL_ELEMENTS_IN_MULTIDIM_ARRAY(I())
          MULTIDIM_ELEM(I(),i)=
	     log10(1+MULTIDIM_ELEM(I(),i)*MULTIDIM_ELEM(I(),i));
    } else
       // Unknown image
       I().init_zeros(projYdim,projXdim);

    // Scale and normalize
    int minGray=0, maxGray=0;
    scale_and_normalize(I(), options->isItemEnabled(mi_Individual_norm),
       minGray, maxGray);

    // If PSD mode, make the full window fit the current size
    if (load_mode==PSD_mode || load_mode==CTF_mode)
       I().self_scale_to_size(rowHeight(0),columnWidth(0));

    // Convert Xmipp image to Pixmap
    content[i]=new QPixmap;
    xmipp2Pixmap(I,content[i],minGray,maxGray,0,0,treat_separately_left_right);
}

/* Resize event ------------------------------------------------------------ */
void ShowSel::resizeEvent(QResizeEvent *event) {
   if (load_mode==Normal_mode)
      QTable::resizeEvent(event);
   else {
      int Xdim=(event->size().width()-4)/NumCols;
      int Ydim=(event->size().height()-4-status->height())/NumRows;
      for (int i=0; i<NumCols; i++) setColumnWidth(i,Xdim);
      for (int i=0; i<NumRows; i++) setRowHeight(i,Ydim);
      
      // Reset all flags so that images are reloaded
      clearContents();
      
      // Repaint
      maxWidth=NumCols*Xdim+4;
      maxHeight=NumRows*Ydim+4;
      adjustStatusLabel();
      repaintContents();
   }
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
   		"Yes", "No") == 0) _SF.write((string)  ((const char *)newfilename));
   	 else QMessageBox::about( this, "Warning!", "Saving aborted\n");			   
      else _SF.write((string)  ((const char *)newfilename));
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
void ShowSel::changeLabels()     {
   if      (options->isItemEnabled(mi_imgAsLabels)) labeltype=Filename_LABEL;
   else if (options->isItemEnabled(mi_selAsLabels)) labeltype=SFLabel_LABEL;
   changeBoolOption(mi_imgAsLabels, mi_selAsLabels);
}

// Show statistics ---------------------------------------------------------
void ShowSel::showStats() {  
    SelFile SFNew;
    for (int i=0; i<listSize; i++)
        if (cellMarks[i])
	  SFNew.insert(imgnames[i], SelLine::ACTIVE);
    if (SFNew.ImgNo()) ShowTable::showStats(SFNew,apply_geo);
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

// Reload all --------------------------------------------------------------
void ShowSel::reloadAll() {
   clearContents();
   repaint();
}

// Show This image ---------------------------------------------------------
void ShowSel::showThisImage() {
   int row=currentRow();
   int col=currentColumn();
   int i=indexOf(row,col);

   ImageViewer *showimg = new ImageViewer(imgnames[i].c_str(), false);
   if (load_mode==Normal_mode)
      showimg->loadImage(imgnames[i].c_str());
   else if (load_mode==PSD_mode)
      showimg->loadImage(imgnames[i].c_str(), 0, 0, ImageViewer::PSD_mode);
   else if (load_mode==CTF_mode)
      showimg->loadImage(imgnames[i].c_str(), 0, 0, ImageViewer::CTF_mode);
   showimg->show();
}

// Edit CTF model ----------------------------------------------------------
void ShowSel::editCTFmodel() {
   if (fn_assign=="") {
      QMessageBox::about( this, "Error!", "No Assign CTF file provided\n");
      return;
   }
   // Read the Assign CTF parameters
   Prog_assign_CTF_prm assign_ctf_prm;
   assign_ctf_prm.read(fn_assign);
   
   // Check if the CTF is computed at each particle
   FileName fn_root=assign_ctf_prm.image_fn.remove_all_extensions();
   FileName fn_param;

   // Get the piece name
   if (assign_ctf_prm.compute_at_particle) {
      int i=indexOf(currentRow(),currentColumn());
      fn_param=imgnames[i].without_extension()+".ctfparam";
   } else if (assign_ctf_prm.micrograph_averaging) {
      // If it is the average of the micrograph
      if (assign_ctf_prm.PSD_mode==Prog_assign_CTF_prm::ARMA)
           fn_param=fn_root+"_ARMAavg.ctfparam";
      else fn_param=fn_root+"_Periodogramavg.ctfparam";
   } else {
      // If the micrograph was divided into pieces
      if (assign_ctf_prm.PSD_mode==Prog_assign_CTF_prm::ARMA)
           fn_param=fn_root+"_ARMA";
      else fn_param=fn_root+"_Periodogram";
      // Get the piece to edit
      int i=indexOf(currentRow(),currentColumn())+1;
      fn_param+=ItoA(i,5);
      fn_param+=".ctfparam";
   }

   // Edit the CTF
   system(((string)"xmipp_edit -i "+fn_param+" &").c_str());
}

// Recompute CTF model -----------------------------------------------------
void ShowSel::recomputeCTFmodel() {
   if (fn_assign=="") {
      QMessageBox::about( this, "Error!", "No Assign CTF file provided\n");
      return;
   }

   // Read the Assign CTF parameters
   Prog_assign_CTF_prm assign_ctf_prm;
   assign_ctf_prm.read(fn_assign);
   
   // Get the PSD name
   FileName fn_root=assign_ctf_prm.image_fn.remove_all_extensions();
   FileName fn_psd;
   if (assign_ctf_prm.compute_at_particle) {
      int i=indexOf(currentRow(),currentColumn());
      fn_psd=imgnames[i].without_extension()+".psd";
   } else if (assign_ctf_prm.micrograph_averaging) {
      // If it is the average of the micrograph
      if (assign_ctf_prm.PSD_mode==Prog_assign_CTF_prm::ARMA)
           fn_psd=fn_root+"_ARMAavg.psd";
      else fn_psd=fn_root+"_Periodogramavg.psd";
   } else {
      // If the micrograph was divided into pieces
      if (assign_ctf_prm.PSD_mode==Prog_assign_CTF_prm::ARMA)
           fn_psd=fn_root+"_ARMA";
      else fn_psd=fn_root+"_Periodogram";
      // Get the piece to recompute
      int i=indexOf(currentRow(),currentColumn())+1;
      fn_psd+=ItoA(i,5);
      fn_psd+=".psd";
   }

   // Show this image in a separate window to select the main parameters
   AssignCTFViewer *prm_selector=new AssignCTFViewer(fn_psd,assign_ctf_prm);
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
