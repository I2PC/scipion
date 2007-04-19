/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Roberto Marabini (roberto@mipg.upenn.edu)
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

#include <fstream>

#include <data/docfile.h>
#include <interface/spider.h>

#include "recons_spider.h"

// SIRT --------------------------------------------------------------------
void SIRT_Spider(SelFile &SF, double lambda, double no_it, int radius,
   const FileName &fn_root, const FileName &fn_ext,
   const FileName &fn_recons_root, const FileName &fn_batch) {
      // Generate Spider selection file
      DocFile DF;
      generate_Spider_count(SF.ImgNo(),DF);
      FileName fn_spider_sel=fn_root+"sel."+fn_ext;
      DF.write(fn_spider_sel);

      // Save angles
      FileName fn_spider_ang=fn_root+"ang."+fn_ext;
      extract_angles(SF,DF,"psi","tilt","rot");
      DF.write(fn_spider_ang);

      // Generate appropiate Spider batch file
      ofstream fh_batch;
      FileName full_fn_batch=fn_batch+"."+fn_ext;
      fh_batch.open(full_fn_batch.c_str(),ios::out);
      if (!fh_batch)
         REPORT_ERROR(3005,(string)"single_recons_test:: Could not open "
	    "the file "+full_fn_batch+" for output");

      SF.go_first_ACTIVE();
      fh_batch << "bp rp\n";
      fh_batch << SF.get_current_file() << endl;
      fh_batch << fn_spider_sel << endl;
      fh_batch << radius << endl;
      fh_batch << fn_spider_ang << endl;
      fh_batch << "n\n";
      fh_batch << fn_recons_root << endl;
      fh_batch << "(" << lambda << ",0)\n";
      fh_batch << "(" << no_it << ",1)\n";
      fh_batch << "(0.,2.)\n";
      fh_batch << "(" << 0.95*1/(1+6*lambda) << ")\n";
      fh_batch << "en\n";
      fh_batch.close();

      // Call Spider
      string command_line=(string)"spider "+fn_ext+" "+fn_batch;
      cerr << "Reconstructing with SIRT Spider ...\n";
      system(command_line.c_str());

      // Rename Spider output
      command_line=(string)"mv "+fn_recons_root+"."+fn_ext+" "+
         fn_recons_root+".vol";
      cerr << "Renaming Spider output ...\n";
      system(command_line.c_str());
}
