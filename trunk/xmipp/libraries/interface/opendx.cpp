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

#include "opendx.h"

void openDX::openDXFile(FileName openDXname)
{
   openDX::fh_out.open(openDXname.c_str(), ios::out); //trunc file
   if ( fh_out.fail() ) {cerr << "Cant open file: "
                              << openDXname.c_str()<<endl;
                         exit( 0 );
                        }

//do not change the size of this header without updating the destructor
   openDX::number_of_elements=0;
   openDX::fh_out << "# The following example describes an irregular grid."
                      << endl;
   openDX::fh_out << "object 1 class array type float rank 1 shape 3 "
                      << "items XXXXX data follows" << endl;
}/* openDX */


void openDX::Add_Item(const matrix1D<double> XYZ)
{
double X,Y,Z;
X = XX(XYZ);Y = YY(XYZ);Z = ZZ(XYZ);
openDX::number_of_elements++;
openDX::fh_out << " "<< X << " " << Y << " " << Z << endl;
}//add sphere

openDX::~openDX()
   {
   openDX::fh_out<< "# The data, which is in a one-to-one"
                     << " correspondence with the positions\n"
                     << "object 2 class array type float rank"
                     << " 0 items " <<  openDX::number_of_elements
                     << " data follows";
   for (int ii=0; ii<openDX::number_of_elements; ii++)
       {
       if(ii%10==0)
          openDX::fh_out << endl;
       openDX::fh_out <<"1 ";
       }

   openDX::fh_out << "\nattribute \"dep\" string \"positions\""
                      << "\n# the field, with three components: "
                      << "\n#\"positions\" and \"data\"\n";
   openDX::fh_out << "object \"irregular positions\" class field\n"
                      << "component \"positions\" value 1\n"
                      << "component \"data\" value 2\n"
                      << "end\n";
//write right number of point

  fh_out.seekp(106, ios::beg);
  openDX::fh_out.width(5);
  openDX::fh_out << openDX::number_of_elements;

}


