/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Roberto Marabini (roberto@mipg.upenn.edu)
 *
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

#include "../xmippOpenDXang.hh"
void openDXang::openDXangFile(FileName openDXangname)
{  
 
  FileName header_name, data_name;
  
  header_name=openDXangname.add_extension("general");
  data_name=openDXangname.add_extension("data");
   
 //header
   openDXang::fh_out_header.open(header_name.c_str(), ios::out); //trunc file
   if ( fh_out_header.fail() ) {cerr << "Cant open file: " 
                              << header_name.c_str()<<endl;
                         exit( 0 );
                        }

//do not change the size of this header without updating the destructor
   openDXang::number_of_elements=0;
   
   openDXang::fh_out_header << "file = " <<openDXangname.c_str()
   			     << "\n"
   			     << "points = XXXXX" 
			     << "\n"
			     << "format = ascii\n"  
			     << "interleaving = field\n"
			     << "field = locations\n"
			     << "structure = 3-vector\n"
			     << "type = float\n\n"
                      	     << "end" << endl;
   
   
		      
//data
   openDXang::fh_out_data.open(data_name.c_str(), ios::out); //trunc file
   if ( fh_out_data.fail() ) {cerr << "Cant open file: " 
                              << data_name.c_str()<<endl;
                         exit( 0 );
                        }

		  
}/* openDXang */   


void openDXang::Add_Item(const matrix1D<double> RotTiltPsi)
{
double X,Y,Z;
X = XX(RotTiltPsi);Y = YY(RotTiltPsi);Z = ZZ(RotTiltPsi);
openDXang::number_of_elements++;
openDXang::fh_out_data << " "<< X << " " << Y << " " << Z << endl;


}

openDXang::~openDXang()
   {
     
   
//write right number of point 

  fh_out_header.seekp( 20, ios::beg);
  openDXang::fh_out_header.width(5);
  openDXang::fh_out_header << openDXang::number_of_elements;
   
  
}   


