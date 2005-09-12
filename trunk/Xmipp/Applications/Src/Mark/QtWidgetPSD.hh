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

#ifndef __QT_WIDGET_PSD_HH__
#define __QT_WIDGET_PSD_HH__

#include <qwidget.h>
#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippMicrograph.hh>
#include <Reconstruction/Programs/Prog_assign_CTF.hh>

/* Widget for the PSDs --------------------------------------------------- */
class QtWidgetPSD : public QWidget {   
   Q_OBJECT
public:
   /// Show the CTF instead the PSD
   bool ctf_mode;

   /// Filename of the assign CTF parameters file
   FileName fn_assign_CTF;

   /// Empty constructor
   QtWidgetPSD() {ctf_mode=false;}

   /// Destructor
   ~QtWidgetPSD();

   /// Set CTF mode
   void set_CTF_mode() {ctf_mode=true;}

   /// Set the assign CTF parameters file
   void set_assign_CTF_file(Micrograph &m,
      const FileName &_fn_assign_CTF); 
public:
   // Assign CTF parameters
   Prog_assign_CTF_prm assign_ctf_prm;

   // Filenames to remove when this object is destroyed
   vector<FileName> files_to_remove;
};

#endif
