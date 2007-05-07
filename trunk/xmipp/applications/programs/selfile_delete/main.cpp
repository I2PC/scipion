/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
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

#include <data/selfile.h>

#include <cstdlib>

/* Prototypes -============================================================= */

void Usage (char **argv);


int main (int argc, char *argv[]) {

/* Input Parameters ======================================================== */
FileName       sel_file;   // selection file

/* Parameters ============================================================== */
   try {
       if (argc != 2) {Usage(argv); exit(0);} else {
       		sel_file = argv[1];
       }
   }
   catch (Xmipp_error XE) {cout << XE; Usage(argv);}

try {

/* Perform copy or move =================================================== */

  // Finds last slash
  string org_path;
  int break_point = -1;
  for(int i = sel_file.size()- 1; i >= 0; i--) {
    if (sel_file[i] == '/') {
    	break_point = i;
	break;
    }
  }

  // Copy only the path  	
  if (break_point >=0) {
    org_path.resize(break_point+1);
    for(int j = 0; j <= break_point; j++)
      org_path[j] = sel_file[j];
  }

   SelFile SF(sel_file);
   string comStr;
   while (!SF.eof()) {
      // Get file
      SelLine line= SF.current();
      if (line.Is_data()) { 		//The SelLine is not a comment
       FileName in_name = line.get_text();
       comStr = "rm -f " + org_path + in_name;       	

       if (!system(comStr.c_str()))
	   cout << " file " << org_path << in_name << " removed " << endl;
      }
      SF.next();
   }  // while

   	
   // now remove sel file
   comStr = "rm -f " + sel_file;			
   system(comStr.c_str());

} catch (Xmipp_error XE) {cout << XE;}
   exit(0);
} //main

/* ------------------------------------------------------------------------- */
/* Help Message for this Program                                             */
/* ------------------------------------------------------------------------- */
void Usage (char **argv) {
  printf (
     "Usage: %s [Purpose and Parameters]"
     "\nPurpose: Remove the images in a sel file (and the sel file) "
     "\nParameter Values: (note space before value)"
     "\nI/O parameters"
     "\n    input_file    input sel file"
     "\n  "
     "\nExample: "
     "\n    rmsel c3u.sel "
     "\n    (will remove all images in c3u.sel file ) "

     "\n"
     ,argv[0]);
}
