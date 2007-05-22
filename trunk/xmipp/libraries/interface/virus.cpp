/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@mipg.upenn.edu)
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

#include "virus.h"

#include <data/args.h>

void VirusEulerMatrices::read(const FileName &fn)
{
    ifstream  fh_Euler;
    int       line_no = 0;
    string    line;
    int       i;

    //alloc Memory
    for (i = 0;i < Vir_Eq_Views;i++)
        E_Matrices[i].resize(3, 3);

    // Open file
    fh_Euler.open(fn.c_str(), ios::in);
    if (!fh_Euler)
        REPORT_ERROR(1601, "VirusEulerMatrices::read: File " + fn + " not found");

    fh_Euler.peek();
    while (!fh_Euler.eof())
    {
        try
        {
            getline(fh_Euler, line);
            if (line[0] == 0)    continue;
            if (line[0] == '#')  continue;
            if (line[0] == '\n') continue;

            E_Matrices[(int)(line_no/3)](line_no % 3, 0) = AtoF(first_token(line));
            E_Matrices[(int)(line_no/3)](line_no % 3, 1) = AtoF(next_token());
            E_Matrices[(int)(line_no/3)](line_no % 3, 2) = AtoF(next_token());
//cout<<line<<endl;
        }
        catch (Xmipp_error)
        {
            cout << "Euler File: Line " << line_no << " is skipped due to an error\n";
        }
        line_no++;
        fh_Euler.peek();
//if(line_no%3==0)
//{
//cout<<E_Matrices[(int)((line_no-1)/3)];
//exit(1);
//}
    }/* while */

    // Close file
    fh_Euler.close();
//   fh_Euler=fn;

}/*  VirusEulerMatrices::read */
