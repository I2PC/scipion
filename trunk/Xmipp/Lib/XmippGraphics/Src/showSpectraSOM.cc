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
    hisAssigned = NULL;
    cv_errors   = NULL;
    infStr      = "";
    ShowSpectra::init();
}

void ShowSpectraSOM::clear() {
    if (Vdat        != NULL) delete    Vdat;
    if (SFcv        != NULL) delete [] SFcv;
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
	  while (!fh_class.eof()) {
	       int row, col;
	       float tmp;
	       fh_class >> col >> row >> tmp;
	       getline(fh_class, line);
	       int i=row*NumCols+col;
	       SFcv[i].push_back(first_token(line)); 
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
   options->insertItem( "View average and SD Images",  this,  SLOT(showStats()));
   options->insertItem( "Show average and SD Represented Images", this,  SLOT(showRepresentedStats()));
   options->insertItem( "View assigned images", this,  SLOT(showRepresentedSel()), CTRL+Key_A );
   options->insertItem( "View error Image", this,  SLOT(showErrorImage()), CTRL+Key_E );
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
void ShowSpectraSOM::extractRepresented(vector<string> &SF_represented) {
   for (int i=0; i<listSize; i++)
      if (cellMarks[i]) {
         int jmax=SFcv[i].size();
	 for (int j=0; j<jmax; j++)
   	    SF_represented.push_back((SFcv[i])[j]);
      }
}

void ShowSpectraSOM::saveAssigned() {
   QMessageBox::about( this, "Error!", "Not implemented\n");
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
void ShowSpectraSOM::showRepresentedStats() {
   QMessageBox::about( this, "Error!", "Not implemented\n");
}

/* Show assigned sel ------------------------------------------------------- */
void ShowSpectraSOM::showRepresentedSel() {
   QMessageBox::about( this, "Error!", "Not implemented\n");
}

/* Show error image -------------------------------------------------------- */
void ShowSpectraSOM::showErrorImage() {  
   QMessageBox::about( this, "Error!", "Not implemented\n");
}

/* Show algorithm information ---------------------------------------------- */
void ShowSpectraSOM::showAlgoInfo() {
   QMessageBox::about( (QWidget*)this, "Algorithm Information", infStr);
}
