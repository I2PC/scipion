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

#include "../showSpectraSOM.hh"
#include <qmessagebox.h>

// Init/Clear -------------------------------------------------------------- */
void ShowSpectraSOM::init() {
    Vdat        = NULL;
    SFcv        = NULL;
    SFcvs       = NULL;
    hisAssigned = NULL;
    cv_errors   = NULL;
    infStr      = "";
    ShowSpectra::init();
}

void ShowSpectraSOM::clear() {
    if (Vdat        != NULL) delete    Vdat;
    if (SFcv        != NULL) delete [] SFcv;
    if (SFcvs       != NULL) delete [] SFcvs;
    if (hisAssigned != NULL) delete [] hisAssigned;
    if (cv_errors   != NULL) delete [] cv_errors;
    infStr = "";
    ShowSpectra::clear();
}

/* Initialize with a SOM file.---------------------------------------------- */
void ShowSpectraSOM::initWithFile(const FileName &_fn_root,
   const FileName &_fn_dat) {
   init();
   readFile(_fn_root);
   
   ifstream fh_in(_fn_dat.c_str());
   if (!fh_in)
      REPORT_ERROR(1,(string)"ShowSpectra::readFile: Cannot open"+_fn_dat);
   Vdat=new xmippCTVectors(fh_in);
   fh_in.close();
   
   initTable();
   initRightclickMenubar();
   repaint();
}

/* Read a SOM -------------------------------------------------------------- */
void ShowSpectraSOM::readFile(const FileName &_fn_root) _THROW {
    clear();
    fn              = _fn_root;
    setCaption(fn.c_str());
    readSOMFiles(_fn_root);
    readDatFile(_fn_root+".cod");
}

void ShowSpectraSOM::readSOMFiles(const FileName &_fn_root) _THROW {
    FileName fn_class, fn_his, fn_err, fn_inf;
    fn_class=_fn_root+".vs";
    fn_his  =_fn_root+".his";
    fn_err  =_fn_root+".err";
    fn_inf  =_fn_root+".inf";
    
    // Read histogram
    if (exists(fn_his)) {
     	ifstream fh_his(fn_his.c_str());
     	if (fh_his) {
     	   xmippCTVectors ts(0, true);
     	   fh_his >> ts;
	   int imax=ts.size();
     	   hisAssigned = new string[imax];
     	   for (int i=0; i<imax; i++)
               hisAssigned[i] = ts.theTargets[i];		   
   	}
	fh_his.close();
    }

    // Read errors
    if (exists(fn_err)) {
     	ifstream fh_err(fn_err.c_str());
     	if (fh_err) {
     	   xmippCTVectors ts(0, true);
     	   fh_err >> ts;
	   int imax=ts.size();
     	   cv_errors = new string[imax];
     	   for (int i=0; i<imax; i++)
               cv_errors[i] = ts.theTargets[i];		   
   	}
	fh_err.close();
    }

    // Read inf file
    if (exists(fn_inf)) {
       ifstream fh_inf(fn_inf.c_str());
       if (fh_inf) {
          string line;
          getline(fh_inf, line);
          infStr = line.c_str(); infStr += "\n";
          while (!fh_inf.eof()) {
              getline(fh_inf, line);
              infStr += line.c_str(); infStr += "\n";
          }	
        }
        fh_inf.close();
    }
    
    // Read codevectors
    if (exists(fn_class)) {
       ifstream fh_class(fn_class.c_str());
       if (fh_class) {
      	  string line, fn;
          int dim;
	  string topol, neigh;
          fh_class >> dim >> topol >> NumCols >> NumRows >> neigh; 
	  listSize = NumCols*NumRows;
	  if (listSize==0)
	     REPORT_ERROR(1,"ShowSpectraSOM::readFile: Input file is empty");
          getline(fh_class, line);

	  if (infStr=="") {
	     infStr = "Kohonen SOM algorithm\n\n";
	     infStr += "Number of variables: ";
	     line = ItoA(dim);
	     infStr += line.c_str();
	     infStr += "\n";
	     infStr += "Horizontal dimension (Xdim): ";
	     line = ItoA(NumCols);
	     infStr += line.c_str();
	     infStr += "\n";
	     infStr += "Vertical dimension (Ydim): ";
	     line = ItoA(NumRows);
	     infStr += line.c_str();
	     infStr += "\n";
	     infStr += "Topology : ";
	     infStr += topol.c_str();
	     infStr += "\n";
	     infStr += "Neighborhood function : ";
	     infStr += neigh.c_str();
	     infStr += "\n";
	  }

          SFcv=new vector<string>[listSize];
          SFcvs=new vector<int>[listSize];
          int j=0;
	  while (!fh_class.eof()) {
	       int row, col;
	       float tmp;
	       fh_class >> col >> row >> tmp;
	       getline(fh_class, line);
	       int i=row*NumCols+col;
	       SFcv[i].push_back(first_token(line)); 
	       SFcvs[i].push_back(j);
               j++; 
	  }
       }
       fh_class.close();
    }
}

/* Initialize right click menubar ------------------------------------------ */
void ShowSpectraSOM::initRightclickMenubar() {
   menubar = new QPopupMenu(); 
   QPopupMenu * file = new QPopupMenu();
      file->insertItem( "Open...", this,  SLOT(GUIopenFile()));
      file->insertItem( "Save assigned images in a sel file...",
         this, SLOT(saveAssigned()), CTRL+Key_N);
   menubar->insertItem( "&File", file );

   // Options .............................................................
   options =  new QPopupMenu();
   setCommonOptionsRightclickMenubar();

   setCommonSpectraOptionsRightclickMenubar();

   // What kind of labels
   mi_imgAsLabels = options->insertItem( "Show Image Names as Labels", this,  SLOT(changeLabelsToImg()));
   mi_selAsLabels = options->insertItem( "Show Sel status as Labels", this,  SLOT(changeLabelsToSel()));
   mi_hisAsLabels = options->insertItem( "Show Histogram as Labels", this, SLOT(changeLabelsToHis()));
   mi_errAsLabels = options->insertItem( "Show Errors as Labels", this, SLOT(changeLabelsToErr()));
   options->setItemEnabled(mi_imgAsLabels, true);
   options->setItemEnabled(mi_selAsLabels, true);
   options->setItemEnabled(mi_hisAsLabels, false);
   options->setItemEnabled(mi_errAsLabels, true);
   labeltype=Histogram_LABEL;
   options->insertSeparator();

   // Statistics
   options->insertItem( "View average and SD Represented Spectra",  this,  SLOT(showRepresentedSpectraStats()));
   options->insertItem( "View average and SD Represented Images", this,  SLOT(showRepresentedImagesStats()));
   options->insertItem( "View Represented Spectra", this,  SLOT(showRepresentedSpectra()));
   options->insertItem( "View Represented Images", this,  SLOT(showRepresentedSel()));
   options->insertItem( "View error spectrum", this,  SLOT(showErrorSpectrum()), CTRL+Key_E );
   options->insertSeparator();
   options->insertItem( "Show Algorithm Information", this,  SLOT(showAlgoInfo()), CTRL+Key_G );
   options->insertSeparator();

   // Insert options the menu
   menubar->insertItem( "&Options", options );    
   menubar->insertSeparator();

   // Inser Help and Quit
   insertGeneralItemsInRightclickMenubar();
}

/* Extract represented.----------------------------------------------------- */
void ShowSpectraSOM::extractRepresented(SelFile &SF_represented) {
   for (int i=0; i<listSize; i++)
      if (cellMarks[i]) {
         int jmax=SFcv[i].size();
	 for (int j=0; j<jmax; j++)
   	    SF_represented.insert((SFcv[i])[j],SelLine::ACTIVE);
      }
}

void ShowSpectraSOM::extractRepresented(xmippCTVectors &_v_represented) {
   // Count the number of represented vectors
   int counter=0;
   for (int i=0; i<listSize; i++)
      if (cellMarks[i]) counter+=SFcvs[i].size();

   // Resize the represented vectors
   _v_represented.theItems.resize(counter);
   _v_represented.theTargets.resize(counter);

   // Now copy the items indicated by SFcvs
   long myIndex = 0;
   for (long i=0; i<listSize; i++)
     if (cellMarks[i]) {
         int jmax=SFcvs[i].size();
	 for (int j=0; j<jmax; j++) {
            _v_represented.theItems[myIndex] = Vdat->theItems[(SFcvs[i])[j]]; 
            _v_represented.theTargets[myIndex] = Vdat->theTargets[(SFcvs[i])[j]]; 
            myIndex++;
         }
     } 
}

/* Save assigned images ---------------------------------------------------- */
void ShowSpectraSOM::saveAssigned() {
   SelFile SFNew;
   extractRepresented(SFNew);
   if (SFNew.ImgNo()) writeSelFile(SFNew);
   else QMessageBox::about( this, "Error!", "No images selected\n");
}

/* Change labels ----------------------------------------------------------- */
void ShowSpectraSOM::changeLabelsToImg() {changeLabel(mi_imgAsLabels);}
void ShowSpectraSOM::changeLabelsToSel() {changeLabel(mi_selAsLabels);}
void ShowSpectraSOM::changeLabelsToHis() {changeLabel(mi_hisAsLabels);}
void ShowSpectraSOM::changeLabelsToErr() {changeLabel(mi_errAsLabels);}
void ShowSpectraSOM::changeLabel(int _clicked_mi) {
   options->setItemEnabled(mi_imgAsLabels, true);
   options->setItemEnabled(mi_selAsLabels, true);
   options->setItemEnabled(mi_hisAsLabels, true);
   options->setItemEnabled(mi_errAsLabels, true);
   options->setItemEnabled(_clicked_mi,    false);
   if      (_clicked_mi==mi_imgAsLabels) labeltype=Filename_LABEL;
   else if (_clicked_mi==mi_selAsLabels) labeltype=SFLabel_LABEL;
   else if (_clicked_mi==mi_hisAsLabels) labeltype=Histogram_LABEL;
   else if (_clicked_mi==mi_errAsLabels) labeltype=Err_LABEL;
   repaintContents();
}

const char * ShowSpectraSOM::cellLabel(int i) const {
   if (options->isItemEnabled(mi_showLabel)) return NULL;
   switch (labeltype) {
      case SFLabel_LABEL:
         return (selstatus[i])? "1":"-1";
      case Filename_LABEL:
         return imgnames[i].c_str();
      case Histogram_LABEL:
         return hisAssigned[i].c_str();
      case Err_LABEL:
         return cv_errors[i].c_str();
   }
}

/* Show Average and SD of the represented images --------------------------- */
void ShowSpectraSOM::showRepresentedImagesStats() {
    SelFile SFNew;
    extractRepresented(SFNew);
    if (SFNew.ImgNo()) ShowTable::showStats(SFNew);
    else QMessageBox::about( this, "Error!", "No images selected\n");
}

void ShowSpectraSOM::showRepresentedSpectraStats() {
   // Extract represented vectors
   xmippCTVectors V_represented(0,true);
   extractRepresented(V_represented);
   // Get Avg and SD
   xmippCTVectors *myVect = new xmippCTVectors(0, true);
   *myVect = V_represented.getStatVector();
   // Represent
   ShowSpectra *myST = new ShowSpectra;
   myST->initWithVectors(1, 2, myVect, "Average and SD of represented spectra");
   myST->show();
}

/* Show assigned sel ------------------------------------------------------- */
void ShowSpectraSOM::showRepresentedSel() {
    SelFile SFNew;
    extractRepresented(SFNew);
    if (SFNew.ImgNo()) {
      ShowSel *showsel=new ShowSel;
      showsel->initWithObject(10,10,SFNew,"Represented images");
      showsel->show();
   }
    else QMessageBox::about( this, "Error!", "No images selected\n");
}

/* Show assigned spectra --------------------------------------------------- */
void ShowSpectraSOM::showRepresentedSpectra() {
   xmippCTVectors *V_represented=new xmippCTVectors(0,true);
   extractRepresented(*V_represented);
   ShowSpectra *myST = new ShowSpectra;
   myST->initWithVectors(6, 6, V_represented, "Represented spectra");
   myST->show();
}

/* Show error spectrum ----------------------------------------------------- */
void ShowSpectraSOM::showErrorSpectrum() {  
   // Extract represented vectors
   xmippCTVectors V_represented(0,true);
   int row=currentRow();
   int col=currentColumn();
   if (row<0 || col<0) return;
   int i=indexOf(row,col);
   int represented_images=SFcvs[i].size();
   if (represented_images==0) {
      QMessageBox::about( this, 
         "Error!", "This vector does not represent any spectrum\n");
      return;
   }

   V_represented.theItems.resize(represented_images);
   V_represented.theTargets.resize(represented_images);
   for (int j=0; j<represented_images; j++) {
      V_represented.theItems[j] = Vdat->theItems[(SFcvs[i])[j]]; 
      V_represented.theTargets[j] = Vdat->theTargets[(SFcvs[i])[j]]; 
   }

   // Get Avg
   xmippCTVectors *myVect = new xmippCTVectors(0, true);
   *myVect = V_represented.getStatVector();
   myVect->deleteRow(1);
   myVect->theTargets[0]="Error spectrum";

   // Compute the difference with the code vector
   int jmax=myVect->theItems[0].size();
   for (int j=0; j<jmax; j++)
      myVect->theItems[0][j]-=V->theItems[i][j];
   
   // Represent
   ShowSpectra *myST = new ShowSpectra;
   myST->initWithVectors(1, 1, myVect, "Error spectrum of the represented spectra");
   myST->show();
}

/* Show algorithm information ---------------------------------------------- */
void ShowSpectraSOM::showAlgoInfo() {
   QMessageBox::about( (QWidget*)this, "Algorithm Information", infStr);
}
