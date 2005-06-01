/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Part of this module has been developed by Lorenzo Zampighi and Nelson Tang
 * Dept. Physiology of the David Geffen School of Medicine
 * Univ. of California, Los Angeles.
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
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include "../xmippFuncs.hh"
#include "../xmippSelFiles.hh"
#include "../xmippImages.hh"
#include "../xmippImagic.hh"

/*****************************************************************************/   
/* SEL FILE LINE        	             	      	             	     */
/*****************************************************************************/
/* Copy Constructor */
SelLine::SelLine(const SelLine &l) {
   line_type = l.line_type;
   text      = l.text;
   label     = l.label;
}

SelLine& SelLine::operator = (const SelLine &SL) {
   if (this!=&SL) {
      line_type = SL.line_type;
      text      = SL.text;
      label     = SL.label;
   }
   return *this;
}

// Returns FALSE if the comparison is nonsense, ie, if line_types are not the
// same or the lines are to be deleted or are not assigned.
// Returns TRUE if l1<l2 otherwise FALSE
int operator < (const SelLine &l1, const SelLine &l2) {
   if      (l1.line_type<l2.line_type) return 1;
   else if (l1.line_type>l2.line_type) return 0;
   else                                return l1.text<l2.text;
}

ostream& operator << (ostream& o, SelLine &SFL) {
   switch (SFL.line_type) {
      case (SelLine::DATALINE):
         o << SFL.text << " " << SFL.label << endl; break;
      case (SelLine::COMMENT):
         o << SFL.text << endl; break;
   }
   return o;
}

istream& operator >> (istream& o, SelLine &SFL) _THROW {
   string   line;
   char     img_name[MAX_FILENAME_LENGTH];
   int      no_elements_read;
   int      label;

   // Get line
   getline(o,line);

   // Initialise target
   SFL.line_type=SelLine::NOT_ASSIGNED;
   SFL.text="";
   SFL.label=SelLine::DISCARDED;
   if (line.length()==0) return o;

   // Check if comment or empty line
   if (line[0]=='#' || line[0]=='\0' || line[0]==';') {
      line[line.length()-1]='\0';
      SFL.line_type=SelLine::COMMENT;
      SFL.text=line;

   // Check if a true "filename label" line
   } else {
      no_elements_read=sscanf(line.c_str(),"%s %d",img_name,&label);
      // *** THE SSCANF CAN BE REPLACED BY A STRING I/O OPERATION
      if (no_elements_read==2) {
         SFL.line_type=SelLine::DATALINE;
         SFL.text=img_name;
         SFL.label=(label>=0) ? SelLine::ACTIVE:SelLine::DISCARDED;
      } else
         REPORT_ERROR(1552,"Format error when reading Selection line");
   }
   return o;
}

/*****************************************************************************/   
/* SEL FILE     	      	             	      	             	     */
/*****************************************************************************/
/* Constructor ------------------------------------------------------------- */
SelFile::SelFile() {
   fn_sel       = "Unnamed";
   no_imgs      = 0;
   current_line = text_line.begin();
}

/* Copy Constructor -------------------------------------------------------- */
SelFile::SelFile(const SelFile &SF) {
   fn_sel       = SF.fn_sel;
   text_line    = SF.text_line;
   no_imgs      = SF.no_imgs;
   current_line = SF.current_line;
}

/* Clear ------------------------------------------------------------------- */
void SelFile::clear() {
   fn_sel="Unnamed";
   text_line.erase(text_line.begin(),text_line.end());
   no_imgs=0;
   current_line = text_line.begin();
}

/* Assignment -------------------------------------------------------------- */
SelFile& SelFile::operator = (const SelFile &SF) {
   if (this!=&SF) {
      fn_sel       = SF.fn_sel;
      text_line    = SF.text_line;
      no_imgs      = SF.no_imgs;
      current_line = SF.current_line;
   }
   return *this;
}

/* Show Sel file ----------------------------------------------------------- */
ostream& operator << (ostream& o, SelFile &SF) {
   vector<SelLine>::iterator current = SF.text_line.begin();
   vector<SelLine>::iterator last    = SF.text_line.end();
   while (current != last) {
      o << *current;
      current++;
   }
   return o;
}

/* Clean ------------------------------------------------------------------- */
void SelFile::clean() {
   vector<SelLine>::iterator current = text_line.begin();
   vector<SelLine>::iterator last    = text_line.end();
   vector<SelLine>::iterator temp;
   while (current != last) {
      if ((*current).line_type==SelLine::DATALINE &&
         (*current).label==SelLine::DISCARDED) {
         temp=current; temp++;
         text_line.erase(current);
         current=temp;
      } else current++;
   }
   current_line=text_line.begin();
}

/* Clean comments ---------------------------------------------------------- */
void SelFile::clean_comments() {
   vector<SelLine>::iterator current = text_line.begin();
   vector<SelLine>::iterator last    = text_line.end();
   vector<SelLine>::iterator temp;
   while (current != last) {
      if ((*current).line_type==SelLine::COMMENT) {
         temp=current; temp++;
         text_line.erase(current);
         current=temp;
      } else current++;
   }
   current_line=text_line.begin();
}

/* Read -------------------------------------------------------------------- */
void SelFile::read(const FileName &sel_name, int overriding) _THROW {
   SelLine   temp;
   ifstream  fh_sel;
   int       line_no=1;
   
   // Empties current SelFile
   if (overriding) clear();
   
   // Open file
   if (sel_name.find (IMAGIC_TAG) == 0) {
      // Read Imagic selfile
      const FileName hed_fname = sel_name.substr(IMAGIC_TAG_LEN);
      ImageImagicInfo info = ImagicGetImgInfo (hed_fname);
      no_imgs = info.num_img;
      temp.line_type = SelLine::DATALINE;
      temp.label = SelLine::ACTIVE;
      for (unsigned i = 0; i < no_imgs; i++) {
          temp.text = ImagicMakeName (hed_fname.c_str(), i);
          text_line.push_back (temp);
      }
   } else {
      // Read normal selfile
      fh_sel.open(sel_name.c_str(), ios::in);
      if (!fh_sel)
         REPORT_ERROR(1551,(string)"SelFile::read: File "+sel_name+" not found");

      // Read each line and keep it in the list of the SelFile object
      fh_sel.peek();
      while (!fh_sel.eof()) {
         try {
            fh_sel >> temp;
         } catch (Xmipp_error) {
            cout << "Sel file: Line " << line_no << " is skipped due to an error\n";
         }
         switch (temp.line_type) {
            case (SelLine::NOT_ASSIGNED): break; // Line with an error
            case (SelLine::DATALINE):
	       if (temp.label!=SelLine::DISCARDED) no_imgs++;
	       text_line.push_back(temp);
	       break;
            case (SelLine::COMMENT):
	       text_line.push_back(temp);
	       break;
          }
          line_no++;
          fh_sel.peek();
      }
     
      // Close file
      fh_sel.close();
   }
   
   // Set "pointer" to the beginning of the file
   if (overriding) fn_sel=sel_name;
   go_first_ACTIVE();
}
/* Merge ------------------------------------------------------------------- */
void SelFile::merge(const FileName &sel_name) {
   SelFile SF(sel_name);
   *this=*this+SF;
   go_first_ACTIVE();
}

/* Write ------------------------------------------------------------------- */
void SelFile::write(const FileName &sel_name) _THROW {
   ofstream    fh_sel;
   vector<SelLine>::iterator current = text_line.begin();
   vector<SelLine>::iterator last    = text_line.end();
   
   if (strcmp(sel_name.c_str(),"")!=0) fn_sel=sel_name;
      // Don't use sel_name=="" because it wastes memory

   if (sel_name.find (IMAGIC_TAG) == 0) {
      // Write Imagic selfile
      const FileName hed_fname = sel_name.substr(IMAGIC_TAG_LEN);
      vector<Image *> imgs;
      for ( ; current != last; current++) {
          Image *img;
          if (current->Is_data() && (current->get_label() == SelLine::ACTIVE) &&
	     (img = Image::LoadImage (current->get_text())))
	     imgs.push_back (img);
      }
      if (!ImagicWriteImagicFile (hed_fname, imgs))
         REPORT_ERROR (1553, "Error writing selfile to Imagic file "+sel_name);
      for (vector<Image *>::iterator i = imgs.begin(); i != imgs.end(); i++)
          delete (*i);
   } else {
      // Write Xmipp selfile
      // Open file
      fh_sel.open(fn_sel.c_str(), ios::out);
      if (!fh_sel)
         REPORT_ERROR(1553,"SelFile::write: File "+fn_sel+" cannot be written");

      // Read each line and keep it in the list of the SelFile object
      while (current != last) fh_sel << *(current++);
     
      // Close file
      fh_sel.close();
   }
}

/* Merging with another selfile -------------------------------------------- */
void SelFile::merge(SelFile &SF) {
   vector<SelLine>::iterator current = SF.text_line.begin();
   vector<SelLine>::iterator last    = SF.text_line.end();
   vector<SelLine>::iterator found;

   SelLine discrepancy;
   discrepancy.line_type=SelLine::COMMENT;
   discrepancy.text="# There were discrepancy in the tags for next line, the "
      "ACTIVE state is kept";

   while (current != last) {
      if ((*current).line_type!=SelLine::DATALINE) {current++; continue;}
      if ((found=find((*current).text))==text_line.end()) {
         // New image not found in the whole Sel File. 
         // Add it if it is not discarded
         if ((*current).label!=SelLine::DISCARDED) {
            text_line.push_back(*current);
            no_imgs++;
         }
      } else
         // New image is found, check that its line is not going
         // to be removed, if it is add it again; else, check if
         // there is a discrepancy between them
         if ((*found).label!=(*current).label) {
            if ((*found).label<(*current).label) {
               (*found).label=SelLine::ACTIVE;
               no_imgs++;
            }
            text_line.insert(found,1,discrepancy);
         }
      current++;
   }
}

/* Merging, operator + ----------------------------------------------------- */
// If the same file is in both Sel Files the label in the first is kept
SelFile SelFile::operator + (SelFile &SF) {
   SelFile result;
   result=*this;
   result.merge(SF);
   return result;
}

/* Split randomly in two equally large selfiles ---------------------------- */
void SelFile::split_in_two(SelFile &SF1,SelFile &SF2) {
  SelFile  SFtmp;
  SF1=*this;
  SFtmp=SF1.randomize();
  SF1.clear();
  int N=SFtmp.ImgNo();
  SF1.reserve(N);
  SF2.reserve(N);
  int half=N/2;
  SFtmp.go_beginning();
  for (int i=0;i<N; i++) {
    if (i<half) SF1.insert(SFtmp.current());
    else SF2.insert(SFtmp.current());
    if (i<N-1) SFtmp.NextImg();
  }  
  SFtmp=SF1.sort_by_filenames();
  SF1=SFtmp;
  SFtmp=SF2.sort_by_filenames();
  SF2=SFtmp;
}

/* Select only part of the selfile for parallel MPI-runs ------------------ */
void SelFile::mpi_select_part(int rank, int size, int &num_img_tot) {

  (*this).clean_comments();
  (*this).clean();
  num_img_tot = (*this).ImgNo();
  int remaining = num_img_tot % size;
  int Npart = (int)(num_img_tot-remaining)/size;
  int myFirst, myLast;
  if ( rank < remaining ) {
    myFirst = rank * (Npart + 1);
    myLast = myFirst + Npart;
  } else {
    myFirst = rank * Npart + remaining;
    myLast = myFirst + Npart - 1;
  }
  // Now discard all images in Selfile that are outside myFirst-myLast
  (*this).go_beginning();
  SelFile  SFpart=*this;
  SFpart.clear();
  for (int nr=myFirst; nr<=myLast; nr++) {
    (*this).go_beginning();
    (*this).jump_lines(nr);
    SFpart.insert((*this).current());
  }
  *this=SFpart;

}

/* Adjust to label --------------------------------------------------------- */
void SelFile::adjust_to_label(SelLine::Label label) {
   if (current_line==text_line.end()) return;
   while ((*current_line).line_type!=SelLine::DATALINE ||
      (*current_line).label!=label) {
      current_line++;
      if (current_line==text_line.end()) return;
   }
}

/* Next Image with a certain label ----------------------------------------- */
string SelFile::NextImg(SelLine::Label label) {
   adjust_to_label(label);
   if (current_line!=text_line.end()) return (*current_line++).text;
   else                               return "";
}

/* Jump over a certain number of data lines (disregarding any label) ------- */
bool SelFile::jump_lines(int how_many) {
  for (int i=0; i<how_many; i++) {
    if (current_line!=text_line.end()) current_line++;
    else return false;
  }
  return true;
}

/* Jump over images with a certain label ----------------------------------- */
void SelFile::jump(int how_many, SelLine::Label label) {
   adjust_to_label(label);
   for (int i=0; i<how_many; i++)
      if (current_line!=text_line.end()) {current_line++; adjust_to_label(label);}
}

/* Find an image (inside the list) ----------------------------------------- */
// It returns a pointer to past-last element if the image is not inside
vector<SelLine>::iterator find(vector<SelLine> &text, string &img_name) {
   vector<SelLine>::iterator current = text.begin();
   vector<SelLine>::iterator last    = text.end();

   while (current != last) {
      if ((*current).line_type==SelLine::DATALINE &&
         (*current).text==img_name) return current;
      current++;
   }
   return current;
}

/* Find an image (inside the Sel File) ------------------------------------- */
// It returns a pointer to past-last element if the image is not inside
// *** THIS SHOULD USE THE PREVIOUS FUNCTION BUT I CANNOT MAKE IT TO COMPILE
vector<SelLine>::iterator SelFile::find(string img_name) 
{
   vector<SelLine>::iterator current = text_line.begin();
   vector<SelLine>::iterator last    = text_line.end();

   while (current != last) {
      if ((*current).line_type==SelLine::DATALINE &&
         (*current).text==img_name) return current;
      current++;
   }
   return current;
}

/* Number of images with a certain label ----------------------------------- */
// If the label is 0 it means any valid image
int SelFile::ImgNo(SelLine::Label label) const {
   int N=0;
   vector<SelLine>::const_iterator current = text_line.begin();
   vector<SelLine>::const_iterator last    = text_line.end();
   while (current != last) {
      if ((*current).line_type==SelLine::DATALINE &&
         (*current).label==label) N++;
      current++;
   }
   return N;
}

/* Number of lines within file --------------------------------------------- */
int SelFile::LineNo() {
   int N=0;
   vector<SelLine>::iterator current = text_line.begin();
   vector<SelLine>::iterator last    = text_line.end();
   while (current != last) {
      N++; current++;
   }
   return N;
}

/* Image size -------------------------------------------------------------- */
void SelFile::ImgSize(int &Ydim, int &Xdim) _THROW {
   vector<SelLine>::iterator aux=current_line;
   go_first_ACTIVE();
   FileName fn_img=(*current_line).text;
   if (fn_img.find("imagic:")!=-1) {
      Image *img = Image::LoadImage(fn_img);
      Ydim=(*img)().ydim;
      Xdim=(*img)().xdim;
      delete img;
   } else if (Is_ImageXmipp(fn_img)) {
      ImageXmipp img;
      img.read(fn_img);
      Ydim=img().ydim;
      Xdim=img().xdim;
   } else if (Is_FourierImageXmipp(fn_img)) {
      FourierImageXmipp img;
      img.read(fn_img);
      Ydim=img().ydim;
      Xdim=img().xdim;
   } else
      REPORT_ERROR(1,"SelFile::ImgSize: First Active file is not an image");
      
  current_line=aux;
}

/* File Extension ---------------------------------------------------------- */
FileName SelFile::FileExtension() {
   vector<SelLine>::iterator aux=current_line;
   go_first_ACTIVE();
   FileName ext=(*current_line).text;
   ext=ext.get_extension();
   current_line=aux;
   return ext;
}

/* Statistics -------------------------------------------------------------- */
void SelFile::get_statistics(Image& _ave, Image& _sd, double& _min,
   double& _max, bool apply_geo) {
   _min = MAXFLOAT; _max = 0; bool first = true; int n = 0;   
   // Calculate Mean
   go_beginning();
   while ((!eof())) {
        string image_name = NextImg();
	if (image_name == "") continue;
        Image *image = Image::LoadImage(image_name,apply_geo); // reads image 
        double min, max, avg, stddev;      
        (*image)().compute_stats(avg, stddev, min, max);
        if (_min > min) _min = min;  
        if (_max < max) _max = max; 
        if (first) {    
           _ave = *image;
	   first = false;	 
        } else {
           _ave() += (*image)();
	}
	delete image;
        n++;
   }
   if (n > 0)
      _ave() /= n;
   _sd = _ave;  
   _sd().init_zeros();
   // Calculate SD
   go_beginning(); 
   while ((!eof())) {
        string image_name = NextImg();
	if (image_name == "") continue;
        Image *image = Image::LoadImage (image_name,apply_geo); // reads image 
        Image tmpImg;
        tmpImg() = (((*image)() - _ave()));   
        tmpImg() *= tmpImg();       
        _sd() += tmpImg();
	delete image;
   }
   _sd() /= (n-1);
   _sd().SQRTnD(); 
}


/* Maximum filename length ------------------------------------------------- */
int SelFile::MaxFileNameLength() {
   vector<SelLine>::iterator aux=current_line;
   int max_length=0;
   go_first_ACTIVE();
   while (!eof()) {
     FileName fn=NextImg();
     max_length=MAX(max_length,fn.length());
   }
   current_line=aux;
   return max_length;
}

/* Get current filename ---------------------------------------------------- */
string SelFile::get_current_file() {
   if (current_line==text_line.end()) return "";
   if ((*current_line).line_type!=SelLine::DATALINE)  return "";
   return (*current_line).text;
}

/* Get filename number i --------------------------------------------------- */
string SelFile::get_file_number(int i) {
   if (i<0) return "";
   vector<SelLine>::iterator current = text_line.begin();
   vector<SelLine>::iterator last    = text_line.end();

   int currenti=0;
   while (current != last) {
      if ((*current).line_type==SelLine::DATALINE &&
         (*current).label==SelLine::ACTIVE) currenti++;
      if (currenti>i) return (*current).text;
      current++;
   }
   return "";
}

/* Remove a certain file --------------------------------------------------- */
void SelFile::remove(string img_name) {
   vector<SelLine>::iterator aux=find(img_name);
   vector<SelLine>::iterator temp;
   if (aux!=text_line.end()) {
      if (aux==current_line) {temp=current_line; temp++;}
      else              temp=current_line;
      if ((*aux).line_type==SelLine::DATALINE) no_imgs--;
      text_line.erase(aux);
      current_line=temp;
   }
}

/* Remove current line ----------------------------------------------------- */
void SelFile::remove_current() {
   if (current_line!=text_line.end()) {
      vector<SelLine>::iterator temp;
      temp=current_line;
      temp++;
      if ((*current_line).line_type==SelLine::DATALINE) no_imgs--;
      text_line.erase(current_line);
      current_line=temp;
   }
}

/* Append a file or change label ------------------------------------------- */
void SelFile::set(string img_name, SelLine::Label label) {
   SelLine temp;
   vector<SelLine>::iterator aux=find(img_name);
   if (aux==text_line.end()) {
      temp.line_type=SelLine::DATALINE;
      temp.text=img_name;
      temp.label=label;
      text_line.push_back(temp);
      if (label!=SelLine::DISCARDED) no_imgs++;
   } else {
      if ((*aux).label!=label) {
         (*aux).label=label;
         if (label!=SelLine::DISCARDED) no_imgs++;
      }
   }
}

/* Append a file or change label ------------------------------------------- */
void SelFile::set_current(SelLine::Label label) {
   if ((*current_line).label!=label) {
      (*current_line).label=label;
      if (label!=SelLine::DISCARDED) no_imgs++;
   }
}

/* Change current filename ------------------------------------------- */
void SelFile::set_current_filename(const FileName &fn_new) {
   if ((*current_line).line_type==SelLine::DATALINE) {
      (*current_line).text=fn_new;
   }
}

/* Insert image before current line ---------------------------------------- */
void SelFile::insert(string img_name, SelLine::Label label) {
   SelLine temp;
   temp.line_type=SelLine::DATALINE;
   temp.text=img_name;
   temp.label=label;
   if (label!=SelLine::DISCARDED) no_imgs++;
   
   // Insert and updates current_line
   current_line=text_line.insert(current_line,temp);
   current_line++;
}

/* Insert line before current line ----------------------------------------- */
void SelFile::insert(const SelLine &_selline) _THROW {
   if (_selline.line_type!=SelLine::DATALINE &&
       _selline.line_type!=SelLine::COMMENT)
      REPORT_ERROR(1552,"SelFile::insert(SelLine): SelLine type not valid");
   if (_selline.line_type==SelLine::DATALINE)
     if (_selline.label!=SelLine::DISCARDED &&
	 _selline.label!=SelLine::ACTIVE)
      REPORT_ERROR(1552,"SelFile::insert(SelLine): SelLine label not valid");

   // Insert and updates current_line
   current_line=text_line.insert(current_line,_selline);
   current_line++;
}

/* Insert a comment before current line ------------------------------------ */
void SelFile::insert_comment(string comment) {
   SelLine temp;
   temp.line_type=SelLine::COMMENT;
   temp.text="# " + comment;
   temp.label=SelLine::DISCARDED;
   
   // Insert and updates current_line
   current_line=text_line.insert(current_line,temp);
   current_line++;
}

/* Sort -------------------------------------------------------------------- */
SelFile SelFile::sort_by_filenames() {
   SelFile result(*this);
   sort(result.text_line.begin(),result.text_line.end());
   result.current_line=result.text_line.begin();
   return result;
}

/* Randomize --------------------------------------------------------------- */
SelFile SelFile::randomize() {
   SelFile  result, aux;
   int      i,j;
   int      rnd_indx;

   randomize_random_generator();
   if (no_imgs==0) return aux;
   aux=*this;
   for (i=no_imgs; i>0; i--) {
      // Jump a random number from the beginning
      rnd_indx=(int) rnd_unif(0,i);
      aux.go_first_ACTIVE(); aux.jump(rnd_indx);
      result.text_line.push_back(*(aux.current_line));
      (*aux.current_line).line_type=SelLine::NOT_CONSIDERED;
   }

   // Adjust remaining fields
   result.no_imgs=no_imgs;
   result.current_line=result.text_line.begin();
   return result;
}

/* Discard randomly a set of images ---------------------------------------- */
SelFile SelFile::random_discard(int N) {
   SelFile  result;
   int      i, j, rnd_indx;
   
   SelLine::Label label=SelLine::ACTIVE;
   result=*this;
   N=min(N,no_imgs);
   for (i=0; i<N; i++) {
      // Jump a random number from the beginning
      rnd_indx=(int) rnd_unif(0,result.no_imgs);
      result.go_first_ACTIVE(); result.jump(rnd_indx,label);

      // Discard that image
      (*(result.current_line)).label=SelLine::DISCARDED;

      // Decrease the number of images such that next time 
      result.no_imgs--;
   }
   
   result.go_beginning();
   return result;
}

/* Compare ----------------------------------------------------------------- */
// Only img_files with the active label are considered
SelFile compare(SelFile &SF1, SelFile &SF2) {
   vector<SelLine>     only_in_SF1;
   vector<SelLine>     only_in_SF2;
   vector<SelLine>     in_both;
   SelFile           result;
   SelLine           temp;
   int               SF1_discarded=0, SF2_discarded=0;
   char              str[10];
   
   // Search in File 1
   vector<SelLine>::iterator current = SF1.text_line.begin();
   vector<SelLine>::iterator last    = SF1.text_line.end();
   vector<SelLine>::iterator last_SF = SF2.text_line.end();
   vector<SelLine>::iterator found;
   
   while (current != last) {
      // Skip if not active
      if ((*current).line_type!=SelLine::DATALINE) {current++; continue;}
      if ((*current).label==SelLine::DISCARDED)
         {SF1_discarded++; current++; continue;}

      // Try to find this archive into Sel File 2
      found=SF2.find((*current).text);
      if (found==last_SF) only_in_SF1.push_back(*current);
      else
         if      ((*found).label==SelLine::DISCARDED)
                     only_in_SF1.push_back(*current);
         else        in_both.push_back(*current);
      current++;      
   }
   
   // Search in File 2
   current = SF2.text_line.begin();
   last    = SF2.text_line.end();
   
   while (current != last) {
      // Skip if not active
      if ((*current).line_type!=SelLine::DATALINE) {current++; continue;}
      if ((*current).label==SelLine::DISCARDED)
         {SF2_discarded++; current++; continue;}

      // Try to find this archive into Sel File 2
      found=find(in_both,(*current).text);
      if (found!=in_both.end()) {current++; continue;}
      only_in_SF2.push_back(*current);
      current++;
   }
   
   // Write Statistics
   temp.line_type=SelLine::COMMENT;
   temp.label=SelLine::DISCARDED;
   temp.text="# Statistics of comparison";
   result.text_line.push_back(temp);
   temp.text="# -------------------------------------------------------------";
   result.text_line.push_back(temp);
   sprintf(str,"%6d",SF1.no_imgs);
   temp.text="# File 1: " + SF1.fn_sel + "(VALID: " + str;
   sprintf(str,"%6d",SF1_discarded);
   temp.text += (string) " DISCARDED: " + str + ")";
   result.text_line.push_back(temp);
   sprintf(str,"%6d",SF2.no_imgs);
   temp.text="# File 2: " + SF2.fn_sel + "(VALID: " + str;
   sprintf(str,"%6d",SF2_discarded);
   temp.text += (string) " DISCARDED: " + str + ")";
   result.text_line.push_back(temp);
   temp.text="";
   result.text_line.push_back(temp);
   sprintf(str,"%6d",in_both.size());
   temp.text=(string)"# Matching Files: " + str;
   result.text_line.push_back(temp);
   sprintf(str,"%6d",only_in_SF1.size());
   temp.text=(string)"# Only in file 1: " + str;
   result.text_line.push_back(temp);
   sprintf(str,"%6d",only_in_SF2.size());
   temp.text=(string)"# Only in file 2: " + str;
   result.text_line.push_back(temp);
   temp.text="# -------------------------------------------------------------";
   result.text_line.push_back(temp);

   // Write files in both
   temp.text="";
   result.text_line.push_back(temp);
   temp.text="# Files in both .sel files";
   result.text_line.push_back(temp);
   current = in_both.begin();
   last    = in_both.end();
   while (current != last)
      result.text_line.push_back(*current++);
  
   // Write files only in Sel File 1
   temp.text="";
   result.text_line.push_back(temp);
   temp.text="# Files only in the first file";
   result.text_line.push_back(temp);
   current = only_in_SF1.begin();
   last    = only_in_SF1.end();
   while (current != last)
      result.text_line.push_back(*current++);
  
   // Write files only in Sel File 2
   temp.text="";
   result.text_line.push_back(temp);
   temp.text="# Files only in the second file";
   result.text_line.push_back(temp);
   current = only_in_SF2.begin();
   last    = only_in_SF2.end();
   while (current != last)
      result.text_line.push_back(*current++);

   // Adjust the remaining fields
   result.no_imgs=in_both.size() +
      only_in_SF1.size() + only_in_SF2.size();
   result.current_line=result.text_line.begin();
  
   return result;
}

/* For all ----------------------------------------------------------------- */
void SelFile::for_all(void (*f)(FileName, FileName), string _ext,
   SelLine::Label _label) {
   vector<SelLine>::iterator current = text_line.begin();
   vector<SelLine>::iterator last    = text_line.end();
   while (current != last) {
      if ((*current).line_type==SelLine::DATALINE && (*current).label==_label) {
         if (_ext!="") (*f)((*current).text,(*current).text+"."+_ext);
         else          (*f)((*current).text,(*current).text);
      }
      current++;
   }
}
