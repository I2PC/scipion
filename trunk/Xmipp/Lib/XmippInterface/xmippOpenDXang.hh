/* author: Roberto Marabini
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
//Draw spheres in openDX
//lines will be follow soon
#ifndef _XMIPPopenDXang_HH
#   define _XMIPPopenDXang_HH

#include <stdlib.h>
#include <fstream>  
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <XmippData/xmippFuncs.hh>
#include <XmippData/xmippMatrices1D.hh>
/**@name openDXang */
//@{
/**  openDXang: 
      This is a class to show grid point using openDXang(tm) (It is usefull to
      view Euler angles. Use the networks for these purpose
      \URL[distance.net]{../Extra_Docs/distances.net}. (To measuure distances select: Options->viewcontrol and then mode->pick, click in teh image select picks->pick_2 and pick again in the image.)
    */
      
class openDXang {
  private:
    ofstream fh_out_header;
    ofstream fh_out_data;
    
    int number_of_elements; 
    
  public:
    /**Open file and write header. Default filenames are myopenDXang.general and
    myopentest.data"
      */      
   void openDXangFile(FileName openDXangfilename="myopenDXang.dx" );
   
   /** Add one grid point 
     */ 
   void Add_Item(const matrix1D<double> RotTiltPsi=vector_R3(0.,0.,0.));
   /** write a footer and Close file, also write the total number of elements
       */
   ~openDXang();   
   };
//@}
#endif
