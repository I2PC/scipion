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

#include "opendxang.h"

void openDXang::openDXangFile(FileName openDXangname_aux)
{

    FileName data_name;
    openDXangname = openDXangname_aux;
    data_name = openDXangname.add_extension("data");

    openDXang::number_of_elements = 0;

//data
    openDXang::fh_out_data.open(data_name.c_str(), std::ios::out); //trunc file
    if (fh_out_data.fail())
    {
        std::cerr << "Cant open file: "
        << data_name.c_str() << std::endl;
        exit(0);
    }


}/* openDXang */


void openDXang::Add_Item(const Matrix1D<double> RotTiltPsi)
{
    double Xp, Yp, Zp;
    double X, Y, Z;

    X = Xp = XX(RotTiltPsi);
    Y = Yp = YY(RotTiltPsi);
    Z = Zp = ZZ(RotTiltPsi);
    if (Y < -90. || Y >= 90)
        Euler_up_down(Xp, Yp, Zp, X, Y, Z);


    openDXang::number_of_elements++;
    openDXang::fh_out_data << " " << X << " " << Y << " " << Z << std::endl;


}

openDXang::~openDXang()
{

    FileName header_name;
    header_name = openDXangname.add_extension("general");

//header
    openDXang::fh_out_header.open(header_name.c_str(), std::ios::out); //trunc file
    if (fh_out_header.fail())
    {
        std::cerr << "Cant open file: "
        << header_name.c_str() << std::endl;
        exit(0);
    }

//write right number of point
    openDXang::fh_out_header << "file = " << openDXang::openDXangname.c_str() << ".data"
    << "\n"
    << "points = " << openDXang::number_of_elements
    << "\n"
    << "format = ascii\n"
    << "interleaving = field\n"
    << "field = locations\n"
    << "structure = 3-vector\n"
    << "type = float\n\n"
    << "end" << std::endl;



}


