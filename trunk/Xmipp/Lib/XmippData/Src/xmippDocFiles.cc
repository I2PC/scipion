/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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

#include <fstream.h>
#include <sstream>
#include <stdio.h>
#include "../xmippDocFiles.hh"
#include "../xmippArgs.hh"

/*****************************************************************************/   
/* DOC FILE LINE        	             	      	             	     */
/*****************************************************************************/
/* Copy Constructor */
DocLine::DocLine(const DocLine &DL) {
   line_type = DL.line_type;
   text      = DL.text;
   key       = DL.key;
   data      = DL.data;
}

DocLine& DocLine::operator = (const DocLine &DL) {
   if (this!=&DL) {
      line_type = DL.line_type;
      text      = DL.text;
      key       = DL.key;
      data      = DL.data;
   }
   return *this;
}

double& DocLine::operator [] (int i) _THROW {
   if (i+1>data.size())
      REPORT_ERROR(1604,"Trying to access to non-existing element " +ItoA(i)+
         " of a document line");
   return data[i];
}

double DocLine::operator [] (int i) const _THROW {
   if (i+1>data.size())
      REPORT_ERROR(1604,"Trying to access to non-existing element " +ItoA(i)+
         " of a document line");
   return data[i];
}

void DocLine::set(int i, double val) {
   // Make sure there is enough memory
   if (i+1>data.size()) data.reserve(i+1);
   
   // Pad with zeros for the non-existing indexes in between
   int space_needed=i+1-data.size();
   for (int k=0; k<space_needed; k++) data.push_back(0);
   
   // Set required data
   data[i]=val;
}

void DocLine::set(const matrix1D<double> &v) {
   data.clear();
   if (line_type!=DATALINE) {
      line_type=DATALINE;
      key=0;
   }

   data.reserve(XSIZE(v));
   for (int i=STARTINGX(v); i<=FINISHINGX(v); i++)
       data.push_back(VEC_ELEM(v,i));
}

void DocLine::clear() {
   line_type = NOT_ASSIGNED;
   text      = "";
   key       = 0;
   data.clear();
}

ostream& operator << (ostream& o, DocLine &DL) {
   char    aux[30];
   switch (DL.line_type) {
      case (DocLine::DATALINE):
         // Print a data line
         sprintf(aux,"%5d ",DL.key);          o << aux;
         sprintf(aux,"%-2d",DL.data.size());  o << aux;
         int imax;
         imax=DL.data.size();
         for (int i=0; i<imax; i++) {
            sprintf(aux," % 10.5f",DL.data[i]); o << aux;
         }
         o << endl;
         break;

      case (DocLine::COMMENT):
         // Print a comment
         o << DL.text << endl; break;
   }
   return o;
}

void DocLine::read(istream& in) _THROW {
   string   line;
   int      key, param_no;

   // Get line
   getline(in,line);

   // Initialise target
   line_type=DocLine::NOT_ASSIGNED;
   text="";
   key=0;
   data.clear();

   // Check if comment or empty line
   int charpos1=line.find_first_not_of(" ");
   if (line[0]=='\0' || line[charpos1]=='#' || line[charpos1]==';') {
      line_type=DocLine::COMMENT;
      text=line;
      data.clear();
      key=0;

   // Read a true document file line
   } else {
      line_type=DocLine::DATALINE;
      text="";
      key=AtoI(first_token(line),1602,"Error reading key");
      param_no=AtoI(next_token(),1602,"Error reading number parameters");
      string auxline=line;
      try {
         // Try unfixed mode first
         read_float_list(NULL,param_no,data,1602,"Error reading doc file line");
      } catch (Xmipp_error) {
         // If doesn't work the try fix
         data.clear();
         data.reserve(param_no);
         for (int i=0; i<param_no; i++) {
            data.push_back(AtoF(line.substr(8+i*12,12)));
	 }
      }
   }
}

/*****************************************************************************/   
/* DOC FILE     	      	             	      	             	     */
/*****************************************************************************/
/* Copy Constructor -------------------------------------------------------- */
DocFile::DocFile(DocFile &DF) {
   fn_doc       = DF.fn_doc;
   m            = DF.m;
   no_lines     = DF.no_lines;
   first_key    = DF.first_key;
   current_line = DF.current_line;
}

/* Clear ------------------------------------------------------------------- */
void DocFile::clear() {
   fn_doc="";
   m.clear();
   no_lines=0;
   first_key=1;
   current_line = m.begin();
}

/* Assignment -------------------------------------------------------------- */
DocFile& DocFile::operator = (const DocFile &DF) {
   if (this!=&DF) {
      fn_doc       = DF.fn_doc;
      m            = DF.m;
      no_lines     = DF.no_lines;
      first_key    = DF.first_key;
      current_line = DF.current_line;
   }
   return *this;
}

/* Assignment from matrix -------------------------------------------------- */
DocFile& DocFile::operator = (const matrix2D<double> &A) {
   clear();
   DocLine temp;
   for (int i=STARTINGY(A); i<=FINISHINGY(A); i++) {
      temp.clear();
      temp.line_type=DocLine::DATALINE;
      temp.data.resize(XSIZE(A));
      for (int j=STARTINGX(A); j<=FINISHINGX(A); j++)
         temp.data[j-STARTINGX(A)]=MAT_ELEM(A,i,j);
      m.push_back(temp);
   }
   fn_doc       = "";
   no_lines     = A.RowNo();
   renum();
   go_beginning();
   return *this;
}

/* Show Doc file ----------------------------------------------------------- */
ostream& operator << (ostream& o, DocFile &DF){
   vector<DocLine>::iterator current = DF.m.begin();
   vector<DocLine>::iterator last    = DF.m.end();
   while (current != last) {
      o << *current;
      current++;
   }
   return o;
}

/* Show a given line ------------------------------------------------------- */
void DocFile::show_line(ostream &o, int key) {
   if (key==-1) {
      if (current_line==m.end()) o << "Current line is at the end of file\n";
      else                       o << *current_line;
   } else {
      vector<DocLine>::iterator line = find(key);
      if (line==m.end()) o << "Key " << key << " not found\n";
      else o << *line;
   }
}

/* Debug Doc file ----------------------------------------------------------- */
void DocFile::debug (){
   vector<DocLine>::iterator current = m.begin();
   vector<DocLine>::iterator last    = m.end();
   while (current != last) {
      if ((*current).line_type==DocLine::DATALINE ||
          (*current).line_type==DocLine::COMMENT)
         cout << *current;
      else {
         char aux[30];
         string str="";
         cout << "Special line\n";
         cout << "  Type: " << (*current).line_type << endl;
         cout << "  Key:  " << (*current).key << endl;
         cout << "  Text: " << (*current).text << endl;
         cout << "  Data: ";
         for (int i=0; i<(*current).data.size(); i++) {
            sprintf(aux," % 11.5f",(*current).data[i]); str += aux;
         }
         cout << str << endl;
      }
      current++;
   }
}

/* Find a key -------------------------------------------------------------- */
// It returns a pointer to past-last element if the key is not inside
vector<DocLine>::iterator DocFile::find(int _key) {
   vector<DocLine>::iterator current = m.begin();
   vector<DocLine>::iterator last    = m.end();

   while (current != last) {
      if ((*current).line_type==DocLine::DATALINE && (*current).key==_key)
         return current;
      current++;
   }
   return current;
}

/* Adjust to data line ----------------------------------------------------- */
void DocFile::adjust_to_data_line() {
   if (current_line==m.end()) return;
   while ((*current_line).line_type!=DocLine::DATALINE) {
      current_line++;
      if (current_line==m.end()) return;
   }
}

/* Renumerate keys --------------------------------------------------------- */
void DocFile::renum() {
   vector<DocLine>::iterator current = m.begin();
   vector<DocLine>::iterator last    = m.end();
   int act_key=first_key;

   while (current != last) {
      if ((*current).line_type==DocLine::DATALINE) (*current).key=act_key++;
      current++;
   }
}

/* Read -------------------------------------------------------------------- */
void DocFile::read(FileName _name, int overriding) _THROW {
   DocLine   temp;
   ifstream  fh_doc;
   int       line_no=1;
   
   // Empties current DocFile
   if (overriding) clear();
   
   // Open file
   fh_doc.open(_name.c_str(), ios::in);
   if (!fh_doc)
      REPORT_ERROR(1601,"DocFile::read: File "+_name+" not found");

   // Read each line and keep it in the list of the DocFile object
   fh_doc.peek();
   while (!fh_doc.eof()) {
      #ifndef _NO_EXCEPTION
      try {
      #endif
         temp.read(fh_doc);
      #ifndef _NO_EXCEPTION
      }
      catch (Xmipp_error) {
         cout << "Doc File: Line " << line_no << " is skipped due to an error\n";
      }
      #endif
      switch (temp.line_type) {
         case (DocLine::NOT_ASSIGNED): break; // Line with an error
         case (DocLine::DATALINE):
            no_lines++;
            m.push_back(temp);
            break;
         case (DocLine::COMMENT):
            m.push_back(temp);
            break;
      }
      line_no++;
      fh_doc.peek();
   }
   
   // Close file
   fh_doc.close();
   
   // Set "pointer" to the beginning of the file
   if (overriding) fn_doc=_name;
   go_first_data_line();
   renum();
}

/* Write ------------------------------------------------------------------- */
void DocFile::write(FileName _name) _THROW {
   ofstream    fh_doc;
   vector<DocLine>::iterator current = m.begin();
   vector<DocLine>::iterator last    = m.end();

   if (_name!="") fn_doc=_name;
   renum();

   // Open file
   fh_doc.open(fn_doc.c_str(), ios::out);
   if (!fh_doc)
      REPORT_ERROR(1603,"DocFile::write: File "+fn_doc+" cannot be written");

   // Read each line and keep it in the list of the SelFile object
   while (current != last) fh_doc << *(current++);

   // Close file
   fh_doc.close();
}

/* Jump over data lines ---------------------------------------------------- */
void DocFile::jump(int how_many) {
   adjust_to_data_line();
   for (int i=0; i<how_many; i++)
      if (current_line!=m.end()) {current_line++; adjust_to_data_line();}
}

/* Locate ------------------------------------------------------------------ */
// It returns a pointer to the next element if the key is not inside
void DocFile::locate(int _key) {
   vector<DocLine>::iterator last = m.end();

   current_line = m.begin();
   while (current_line != last) {
      if ((*current_line).line_type==DocLine::DATALINE &&
          (*current_line).key>=_key)
         return;
      current_line++;
   }
}

/* First data line column number ------------------------------------------- */
int DocFile::FirstLine_ColNo() {
   vector<DocLine>::iterator aux = current_line;
   go_first_data_line();
   int retval=current_line->data.size();
   current_line=aux;
   return retval;
}

/* Last key of the file ---------------------------------------------------- */
int DocFile::get_last_key() {
   vector<DocLine>::iterator last  = m.end();
   vector<DocLine>::iterator first = m.begin();
   do {
      if (last==first) return first_key-1;
      last--;
      if (last->line_type==DocLine::DATALINE) return last->key;
   } while (TRUE);
}

/* Element access to any line ---------------------------------------------- */
double DocFile::operator ()(int _key, int i) _THROW {
   vector<DocLine>::iterator aux = find(_key);
   if (aux==m.end())
      REPORT_ERROR(1604,"DocFile::operator(): The given key (" + ItoA(_key) + ") is not in the file");
   return (*aux)[i];
}

/* Getting angles ---------------------------------------------------------- */
void DocFile::get_angles(int _key, double &rot, double &tilt, double &psi,
   const string &ang1, const string &ang2, const string &ang3) _THROW {
   vector<DocLine>::iterator aux = find(_key);
   if (aux==m.end())
      REPORT_ERROR(1604,"DocFile::operator(): The given key (" + ItoA(_key) + ") is not in the file");
   switch (ang1[0]) {
      case 'r': rot  = (*aux)[0]; break;
      case 't': tilt = (*aux)[0]; break;
      case 'p': psi  = (*aux)[0]; break;
   }
   switch (ang2[0]) {
      case 'r': rot  = (*aux)[1]; break;
      case 't': tilt = (*aux)[1]; break;
      case 'p': psi  = (*aux)[1]; break;
   }
   switch (ang3[0]) {
      case 'r': rot  = (*aux)[2]; break;
      case 't': tilt = (*aux)[2]; break;
      case 'p': psi  = (*aux)[2]; break;
   }
}
/* Getting angles (second triad) ------------------------------------------- */
void DocFile::get_angles1(int _key, double &rot, double &tilt, double &psi,
   const string &ang1, const string &ang2, const string &ang3) _THROW {
   vector<DocLine>::iterator aux = find(_key);
   if (aux==m.end())
      REPORT_ERROR(1604,"DocFile::operator(): The given key (" + ItoA(_key) + ") is not in the file");
   switch (ang1[0]) {
      case 'r': rot  = (*aux)[3]; break;
      case 't': tilt = (*aux)[3]; break;
      case 'p': psi  = (*aux)[3]; break;
   }
   switch (ang2[0]) {
      case 'r': rot  = (*aux)[4]; break;
      case 't': tilt = (*aux)[4]; break;
      case 'p': psi  = (*aux)[4]; break;
   }
   switch (ang3[0]) {
      case 'r': rot  = (*aux)[5]; break;
      case 't': tilt = (*aux)[5]; break;
      case 'p': psi  = (*aux)[5]; break;
   }
}
/* Getting angles (third triad) ------------------------------------------- */
void DocFile::get_angles2(int _key, double &rot, double &tilt, double &psi,
   const string &ang1, const string &ang2, const string &ang3) _THROW {
   vector<DocLine>::iterator aux = find(_key);
   if (aux==m.end())
      REPORT_ERROR(1604,"DocFile::operator(): The given key (" + ItoA(_key) + ") is not in the file");
   switch (ang1[0]) {
      case 'r': rot  = (*aux)[6]; break;
      case 't': tilt = (*aux)[6]; break;
      case 'p': psi  = (*aux)[6]; break;
   }
   switch (ang2[0]) {
      case 'r': rot  = (*aux)[7]; break;
      case 't': tilt = (*aux)[7]; break;
      case 'p': psi  = (*aux)[7]; break;
   }
   switch (ang3[0]) {
      case 'r': rot  = (*aux)[8]; break;
      case 't': tilt = (*aux)[8]; break;
      case 'p': psi  = (*aux)[8]; break;
   }
}

/* Setting angles ---------------------------------------------------------- */
void DocFile::set_angles(int _key, double rot, double tilt, double psi,
   const string &ang1, const string &ang2, const string &ang3) _THROW {
   vector<DocLine>::iterator aux = find(_key);
   if (aux==m.end())
      REPORT_ERROR(1604,"DocFile::operator(): The given key (" + ItoA(_key) + ") is not in the file");
   switch (ang1[0]) {
      case 'r': (*aux)[0]=rot ; break;
      case 't': (*aux)[0]=tilt; break;
      case 'p': (*aux)[0]=psi ; break;
   }
   switch (ang2[0]) {
      case 'r': (*aux)[1]=rot ; break;
      case 't': (*aux)[1]=tilt; break;
      case 'p': (*aux)[1]=psi ; break;
   }
   switch (ang3[0]) {
      case 'r': (*aux)[2]=rot ; break;
      case 't': (*aux)[2]=tilt; break;
      case 'p': (*aux)[2]=psi ; break;
   }
}

/* Element setting in current line ----------------------------------------- */
void DocFile::set(int i, double val) {
   if (current_line!=m.end()) (*current_line).set(i,val);
   else {
      DocLine temp;
      temp.line_type=DocLine::DATALINE;
      temp.set(i,val);
      m.push_back(temp);
      current_line=m.end();
   }
}
   
/* Element setting in any line --------------------------------------------- */
void DocFile::set(int _key, int i, double val) {
   vector<DocLine>::iterator aux = find(_key);
   if (aux==m.end())
      REPORT_ERROR(1604,"DocFile::set: The given key (" + ItoA(_key) + ") is not in the file");

   (*aux).set(i,val);
}

/* Remove a certain key or range ------------------------------------------- */
void DocFile::remove(int _key0, int _keyF) {
   vector<DocLine>::iterator last = m.end();

   if (_keyF==-1) _keyF=_key0;
   int old_key=(*current_line).key;

   current_line = m.begin();
   while (current_line != last) {
      if ((*current_line).line_type!=DocLine::DATALINE) continue;
      if ((*current_line).key>=_key0 && (*current_line).key<=_keyF)
         remove_current();
      else current_line++;
   }
   
   locate(old_key);
}

/* Remove current line ----------------------------------------------------- */
void DocFile::remove_current() {
   if (current_line!=m.end()) {
      vector<DocLine>::iterator temp;
      temp=current_line;
      temp++;
      if ((*current_line).line_type==DocLine::DATALINE) no_lines--;
      m.erase(current_line);
      current_line=temp;
   }
}

/* Insert empty data line before current line ------------------------------- */
int DocFile::insert_data_line(int no_lines_to_insert) {
   DocLine temp;
   temp.line_type=DocLine::DATALINE;

   vector<DocLine>::iterator first_line_inserted;
   for (int i=0; i<no_lines_to_insert; i++) {
       current_line=m.insert(current_line,temp);    // Pointer to the just 
                                                    // inserted line
       if (i==0) first_line_inserted=current_line;  // Annotate it
       current_line++;                              // Go on pointing to the
                                                    // previous current line
   }
   renum();
   no_lines += no_lines_to_insert;
   return (*first_line_inserted).key;
}

/* Insert data line before current line ------------------------------------ */
int DocFile::insert_data_line(const matrix1D<double> &v) {
   DocLine temp;
   temp.set(v); 

   current_line=m.insert(current_line,temp);
   renum();
   int retval=(*current_line).key;
   current_line++;
   no_lines++;
   return retval;
}

/* Insert a comment before current line ------------------------------------ */
void DocFile::insert_comment(string comment) {
   DocLine temp;
   temp.line_type=DocLine::COMMENT;
   temp.text=" ; " + comment;
   
   // Insert and updates current_line
   current_line=m.insert(current_line,temp);
   current_line++;
}

/* Insert document line before current line -------------------------------- */
int DocFile::insert_line(const DocLine &DL) {
   int retval=0;
   current_line=m.insert(current_line,DL);
   if (DL.line_type==DocLine::DATALINE) {
      renum();
      retval=(*current_line).key;
      no_lines++;
   }
   current_line++;
   return retval;
}

/* Append empty data lines ------------------------------------------------- */
int DocFile::append_data_line(int no_lines_to_insert) {
   DocLine temp;
   temp.line_type=DocLine::DATALINE;
   int retval=get_last_key()+1;
   int act_key=retval;

   for (int i=0; i<no_lines_to_insert; i++, act_key++) {
       temp.key=act_key;
       m.push_back(temp);
   }
   no_lines += no_lines_to_insert;
   return retval;
}

/* Append a data line ------------------------------------------------------ */
int DocFile::append_data_line(const matrix1D<double> &v) {
   DocLine temp;
   temp.set(v);
   temp.key=get_last_key()+1;
   m.push_back(temp);
   no_lines++;
   return temp.key;
}

/* Append a comment line --------------------------------------------------- */
void DocFile::append_comment(const string &comment) {
   DocLine temp;
   temp.line_type=DocLine::COMMENT;
   temp.text=" ; " + comment;
   m.push_back(temp);
}

/* Append angles ----------------------------------------------------------- */
int DocFile::append_angles(double rot,  double tilt,  double psi,
   const string &ang1, const string &ang2, const string &ang3){
   matrix1D<double> aux(3);
   if      (ang1[0]=='r')  VEC_ELEM(aux,0)=rot;
   else if (ang1[0]=='t')  VEC_ELEM(aux,0)=tilt;
   else if (ang1[0]=='p')  VEC_ELEM(aux,0)=psi;

   if      (ang2[0]=='r')  VEC_ELEM(aux,1)=rot;
   else if (ang2[0]=='t')  VEC_ELEM(aux,1)=tilt;
   else if (ang2[0]=='p')  VEC_ELEM(aux,1)=psi;

   if      (ang3[0]=='r')  VEC_ELEM(aux,2)=rot;
   else if (ang3[0]=='t')  VEC_ELEM(aux,2)=tilt;
   else if (ang3[0]=='p')  VEC_ELEM(aux,2)=psi;
   
   return append_data_line(aux);
}

/* Append angles ----------------------------------------------------------- */
int DocFile::append_angles(double rot,  double tilt,  double psi,
                           double rot1, double tilt1, double psi1,
   const string &ang1, const string &ang2, const string &ang3) {
   matrix1D<double> aux(6);
   if      (ang1[0]=='r')  VEC_ELEM(aux,0)=rot;
   else if (ang1[0]=='t')  VEC_ELEM(aux,0)=tilt;
   else if (ang1[0]=='p')  VEC_ELEM(aux,0)=psi;

   if      (ang2[0]=='r')  VEC_ELEM(aux,1)=rot;
   else if (ang2[0]=='t')  VEC_ELEM(aux,1)=tilt;
   else if (ang2[0]=='p')  VEC_ELEM(aux,1)=psi;

   if      (ang3[0]=='r')  VEC_ELEM(aux,2)=rot;
   else if (ang3[0]=='t')  VEC_ELEM(aux,2)=tilt;
   else if (ang3[0]=='p')  VEC_ELEM(aux,2)=psi;
   
   if      (ang1[0]=='r')  VEC_ELEM(aux,3)=rot1;
   else if (ang1[0]=='t')  VEC_ELEM(aux,3)=tilt1;
   else if (ang1[0]=='p')  VEC_ELEM(aux,3)=psi1;

   if      (ang2[0]=='r')  VEC_ELEM(aux,4)=rot1;
   else if (ang2[0]=='t')  VEC_ELEM(aux,4)=tilt1;
   else if (ang2[0]=='p')  VEC_ELEM(aux,4)=psi1;

   if      (ang3[0]=='r')  VEC_ELEM(aux,5)=rot1;
   else if (ang3[0]=='t')  VEC_ELEM(aux,5)=tilt1;
   else if (ang3[0]=='p')  VEC_ELEM(aux,5)=psi1;
   
   return append_data_line(aux);
}
/* Append angles ----------------------------------------------------------- */
int DocFile::append_angles(double rot,  double tilt,  double psi,
                           double rot1, double tilt1, double psi1,
                           double rot2, double tilt2, double psi2,
   const string &ang1, const string &ang2, const string &ang3) {
   matrix1D<double> aux(9);
   if      (ang1[0]=='r')  VEC_ELEM(aux,0)=rot;
   else if (ang1[0]=='t')  VEC_ELEM(aux,0)=tilt;
   else if (ang1[0]=='p')  VEC_ELEM(aux,0)=psi;

   if      (ang2[0]=='r')  VEC_ELEM(aux,1)=rot;
   else if (ang2[0]=='t')  VEC_ELEM(aux,1)=tilt;
   else if (ang2[0]=='p')  VEC_ELEM(aux,1)=psi;

   if      (ang3[0]=='r')  VEC_ELEM(aux,2)=rot;
   else if (ang3[0]=='t')  VEC_ELEM(aux,2)=tilt;
   else if (ang3[0]=='p')  VEC_ELEM(aux,2)=psi;
   
   if      (ang1[0]=='r')  VEC_ELEM(aux,3)=rot1;
   else if (ang1[0]=='t')  VEC_ELEM(aux,3)=tilt1;
   else if (ang1[0]=='p')  VEC_ELEM(aux,3)=psi1;

   if      (ang2[0]=='r')  VEC_ELEM(aux,4)=rot1;
   else if (ang2[0]=='t')  VEC_ELEM(aux,4)=tilt1;
   else if (ang2[0]=='p')  VEC_ELEM(aux,4)=psi1;

   if      (ang3[0]=='r')  VEC_ELEM(aux,5)=rot1;
   else if (ang3[0]=='t')  VEC_ELEM(aux,5)=tilt1;
   else if (ang3[0]=='p')  VEC_ELEM(aux,5)=psi1;
   
   if      (ang1[0]=='r')  VEC_ELEM(aux,6)=rot2;
   else if (ang1[0]=='t')  VEC_ELEM(aux,6)=tilt2;
   else if (ang1[0]=='p')  VEC_ELEM(aux,6)=psi2;

   if      (ang2[0]=='r')  VEC_ELEM(aux,7)=rot2;
   else if (ang2[0]=='t')  VEC_ELEM(aux,7)=tilt2;
   else if (ang2[0]=='p')  VEC_ELEM(aux,7)=psi2;

   if      (ang3[0]=='r')  VEC_ELEM(aux,8)=rot2;
   else if (ang3[0]=='t')  VEC_ELEM(aux,8)=tilt2;
   else if (ang3[0]=='p')  VEC_ELEM(aux,8)=psi2;

   return append_data_line(aux);
}

/* Append a document line -------------------------------------------------- */
int DocFile::append_line(DocLine &DL) {
   int retval=0;
   if (DL.line_type==DocLine::DATALINE) {
      no_lines++;
      retval=DL.key=get_last_key()+1;
   }
   m.push_back(DL);
   return retval;
}

/* Clean comments ---------------------------------------------------------- */
void DocFile::clean_comments() {
   vector<DocLine>::iterator last    = m.end();

   current_line=m.begin();
   while (current_line != last) {
      if ((*current_line).line_type==DocLine::COMMENT) remove_current();
      else current_line++;
   }
   current_line=m.begin();
}

/* Randomize --------------------------------------------------------------- */
DocFile DocFile::randomize() {
   DocFile  result, aux;
   int      i,j;
   int      rnd_indx;

   randomize_random_generator();
   if (no_lines==0) return aux;
   aux=*this;
   for (i=no_lines; i>0; i--) {
      // Jump a random number from the beginning
      rnd_indx=(int) rnd_unif(0,i);
      aux.go_first_data_line(); aux.jump(rnd_indx);
      result.m.push_back(*(aux.current_line));
      (*aux.current_line).line_type=DocLine::NOT_CONSIDERED;
   }

   // Adjust remaining fields
   result.no_lines=no_lines;
   result.current_line=result.m.begin();
   result.renum();
   return result;
}

/* Discard randomly a set of images ---------------------------------------- */
DocFile DocFile::random_discard(int N) {
   DocFile  result;
   int      i, j, rnd_indx;
   
   result=*this;
   randomize_random_generator();
   N=min(N,no_lines);
   for (i=0; i<N; i++) {
      // Jump a random number from the beginning
      rnd_indx=(int) rnd_unif(0,result.no_lines);
      result.go_first_data_line(); result.jump(rnd_indx);
      result.remove_current();
   }
   
   result.go_beginning();
   return result;
}

/* Column --> Vector ------------------------------------------------------- */
matrix1D<double> DocFile::col(int _Col) {
   matrix1D<double> result(no_lines);
   vector<DocLine>::iterator current = m.begin();
   vector<DocLine>::iterator last    = m.end();
   int i=0;
   while (current != last) {
      if ((*current).line_type==DocLine::DATALINE)
         VEC_ELEM(result,i++)=(*current)[_Col];
      current++;
   }
   return result;
}

/* Row --> Vector ---------------------------------------------------------- */
matrix1D<double> DocFile::row(int _key) {
   matrix1D<double> result;
   vector<DocLine>::iterator aux = find(_key);
   if (aux==m.end()) return result;
   
   result.resize((*aux).data.size());
   result.setRow();
   for (int i=0; i<result.xdim; i++)
       VEC_ELEM(result,i)=(*aux).data[i];

   return result;
}

/* Vector --> Column ------------------------------------------------------- */
void DocFile::setCol(int _col, matrix1D<double> &v) _THROW {
   go_first_data_line();
   
   for (int i=STARTINGX(v); i<=FINISHINGX(v); i++) {
      set(_col,VEC_ELEM(v,i));
      next_data_line();
   }
   renum();
   
   if (v.get_dim()<m.size())
      REPORT_ERROR(1605,"DocFile::setCol: Column assignment not complete");
}

/* For all ----------------------------------------------------------------- */
void DocFile::for_all_lines(
   void (*f)(const matrix1D<double> &, matrix1D<double> &),
   int key0, int keyF) {
   int current_key;

   // Look for starting point
   if (key0!=-1) locate(key0);
   else          go_beginning();
   
   // While not at the end do
   // Notice that the key range termination condition is inside the loop
   while (!eof()) {
      // Out of range?
      current_key=get_current_key();
      if (keyF!=-1 && current_key>keyF) break;
      
      // If current line is a data line
      if (current_key!=0) {
         // Get current data line, transform it and insert it before
         // the current position
         matrix1D<double> v_in=row(current_key);
         matrix1D<double> v_out;
         f(v_in,v_out);
         insert_data_line(v_out);
         
         // Besides removing the current line, the following removing moves
         // the current line pointer to the next line (either it is
         // a comment or a data line). That is why an adjust to next
         // data line is needed
         remove_current();
         adjust_to_data_line();
         
         // This renumeration is needed because after the insertion the
         // number of lines has been increased by 1, and later a line
         // has been removed leaving a "hole" in the file.
         renum();
      } else
         // Move to next data line
         next_data_line();
   }
   renum();
   go_beginning();
}

/* For a column ------------------------------------------------------------ */
void DocFile::for_column(double (*f)(double), int _col, int _key0, int _keyF) {
   current_line=m.begin();
   vector<DocLine>::iterator last=m.end();
   
   while (current_line!=last) {
      if ((*current_line).line_type==DocLine::DATALINE)
         if (_col==-1)
            for (int i=0; i<(*current_line).data.size(); i++)
                (*current_line).data[i]=f((*current_line).data[i]);
         else
            if (_col<(*current_line).data.size() &&
                (_key0==-1 || (*current_line).key>=_key0) && 
                (_keyF==-1 || (*current_line).key<=_keyF))
                   (*current_line).data[_col]=f((*current_line).data[_col]);
      current_line++;
   }
   renum();
   go_beginning();
}

/* Read document file with Euler angles ------------------------------------ */
int read_Euler_document_file(FileName fn, string ang1, string ang2,
    string ang3, DocFile &DF) {
      DocFile DFaux(fn);
      DocLine DL1, DL2;

      // Set DL2 as a data line and go to the beginning of file
      DL2.set_type(DocLine::DATALINE);
      DFaux.go_beginning();

      // Macro to assign the angle from DL1 in the correcto place of DL2
      // The angle order in DL2 is (rot, tilt, psi)
      #define assign_in_correct_place_of_DL2(angle_descr,angle_index) \
         switch (angle_descr[0]) { \
            case ('r'): DL2.set(0,DL1[angle_index]); break; \
            case ('t'): DL2.set(1,DL1[angle_index]); break; \
            case ('p'): DL2.set(2,DL1[angle_index]); break; \
         }

      // Read the whole file
      while (!DFaux.eof()) {
         DL1=DFaux.get_current_line();
         
         // If DL1 type is neither DATALINE nor COMMENT the line is skipped!!
         if (DL1.Is_data()) {
            // Reorder the angles and insert
            assign_in_correct_place_of_DL2(ang1,0);
            assign_in_correct_place_of_DL2(ang2,1);
            assign_in_correct_place_of_DL2(ang3,2);
            DF.append_line(DL2);
         } else if (DL1.Is_comment())
            // Insert the comment
            DF.append_line(DL1);
         
         // Next line
         DFaux.next();
      }

      return DF.dataLineNo();
}

/* Select images ----------------------------------------------------------- */
void select_images(DocFile &DF, SelFile &SF, int col,
   bool en_limit0, double limit0, bool en_limitF, double limitF) {
   SF.go_beginning();
   DF.go_first_data_line();
   while (!SF.eof()) {
      if (SF.Is_ACTIVE()) {
         if (en_limit0 && DF(col)<limit0) SF.set_current(SelLine::DISCARDED);
         if (en_limitF && DF(col)>limitF) SF.set_current(SelLine::DISCARDED);
      }
      SF.next();
      DF.next_data_line();
   }
}
