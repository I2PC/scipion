/***************************************************************************
 *
 * Authors:     Carlos Oscar Sanchez Sorzano (coss@cnb.uam.es)
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

#include "pdb.h"
#include <fstream>
#include <data/args.h>
#include <data/matrix2d.h>

/* Atom charge ------------------------------------------------------------- */
int atomCharge(const string &atom)
{
    switch (atom[0])
    {
    case 'H':
        return 1;
        break;
    case 'C':
        return 6;
        break;
    case 'N':
        return 7;
        break;
    case 'O':
        return 8;
        break;
    case 'P':
        return 15;
        break;
    case 'S':
        return 16;
        break;
    case 'F': // Iron
        return 26;
        break;
    default:
        return 0;
    }
}

/* Atom radius ------------------------------------------------------------- */
double atomRadius(const string &atom)
{
    switch (atom[0])
    {
    case 'H':
        return 0.25;
        break;
    case 'C':
        return 0.70;
        break;
    case 'N':
        return 0.65;
        break;
    case 'O':
        return 0.60;
        break;
    case 'P':
        return 1.00;
        break;
    case 'S':
        return 1.00;
        break;
    case 'F': // Iron
        return 1.40;
        break;
    default:
        return 0;
    }
}

/* Compute geometry -------------------------------------------------------- */
void computePDBgeometry(const std::string &fnPDB,
   Matrix1D<double> &centerOfMass,
   Matrix1D<double> &limit0, Matrix1D<double> &limitF)
{
    // Initialization
    centerOfMass.initZeros(3);
    limit0.initZeros(3);
    limitF.initZeros(3);
    limit0.init_constant(1e30);
    limitF.init_constant(-1e30);
    double total_mass = 0;

    // Open the file
    std::ifstream fh_pdb;
    fh_pdb.open(fnPDB.c_str());
    if (!fh_pdb)
        REPORT_ERROR(1, (string)"computePDBgeometry:"
                     "Cannot open " + fnPDB + " for reading");

    // Process all lines of the file�
    while (!fh_pdb.eof())
    {
        // Read a ATOM line
        std::string line;
        getline(fh_pdb, line);
        if (line == "") continue;
        std::string kind = firstToken(line);
        if (kind != "ATOM") continue;

        // Extract atom type and position
        // Typical line:
        // ATOM    909  CA  ALA A 161      58.775  31.984 111.803  1.00 34.78
        std::string dummy = nextToken();
        std::string atom_type = nextToken();
        dummy = nextToken();
        dummy = nextToken();
        dummy = nextToken();
        double x = textToFloat(nextToken());
        double y = textToFloat(nextToken());
        double z = textToFloat(nextToken());

        // Update center of mass and limits
        if (x < XX(limit0)) XX(limit0) = x;
        else if (x > XX(limitF)) XX(limitF) = x;
        if (y < YY(limit0)) YY(limit0) = y;
        else if (y > YY(limitF)) YY(limitF) = y;
        if (z < ZZ(limit0)) ZZ(limit0) = z;
        else if (z > ZZ(limitF)) ZZ(limitF) = z;
        int weight=atomCharge(atom_type);
        total_mass += weight;
        XX(centerOfMass) += weight * x;
        YY(centerOfMass) += weight * y;
        ZZ(centerOfMass) += weight * z;
    }

    // Finish calculations
    centerOfMass /= total_mass;
    limit0 -= centerOfMass;
    limitF -= centerOfMass;

    // Close file
    fh_pdb.close();
}

/* Apply geometry ---------------------------------------------------------- */
void applyGeometry(const string &fn_in, const string &fn_out,
    const Matrix2D<double> &A, bool centerPDB)
{
    Matrix1D<double> centerOfMass, limit0, limitF;
    if (centerPDB)
    	computePDBgeometry(fn_in, centerOfMass,limit0, limitF);

    // Open files
    std::ifstream fh_in;
    fh_in.open(fn_in.c_str());
    if (!fh_in)
        REPORT_ERROR(1, (std::string)"applyGeometry:"
                     "Cannot open " + fn_in + " for reading");
    std::ofstream fh_out;
    fh_out.open(fn_out.c_str());
    if (!fh_out)
        REPORT_ERROR(1, (std::string)"applyGeometry:"
                     "Cannot open " + fn_out + " for writing");

    // Process all lines of the file
    while (!fh_in.eof())
    {
        // Read a ATOM line
        std::string line, auxLine;
        getline(fh_in, line);
        if (line == "") 
	{
	    fh_out << "\n";
	    continue;
	}
	auxLine=line.c_str();
        std::string kind = firstToken(auxLine);
        if (kind != "ATOM") 
	{
	    fh_out << line << std::endl;
	    continue;
	}

        // Extract atom type and position
        // Typical line:
        // ATOM    909  CA  ALA A 161      58.775  31.984 111.803  1.00 34.78
        std::string t1 = nextToken();
        std::string t2 = nextToken();
        std::string t3 = nextToken();
        std::string t4 = nextToken();
        std::string t5 = nextToken();
        double x = textToFloat(nextToken());
        double y = textToFloat(nextToken());
        double z = textToFloat(nextToken());
        std::string t6 = nextToken();
        std::string t7 = nextToken();

    	Matrix1D<double> v(4);
        if (centerPDB)
	{
	    VECTOR_R3(v,x-XX(centerOfMass),
		y-YY(centerOfMass),z-ZZ(centerOfMass));
	}
	else
	{
    	    VECTOR_R3(v,x,y,z);
	}
	v(3)=1;
	v=A*v;
	fh_out << "ATOM " << t1 << " " << t2 << " " << t3 << " " << t4 << " " 
	       << t5 << " " << XX(v) << " " << YY(v) << " " 
	       << ZZ(v) << " " << t6 << " " << t7 << std::endl;
    }

    // Close files
    fh_in.close();
    fh_out.close();
}
