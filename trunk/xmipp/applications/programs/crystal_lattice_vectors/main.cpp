/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
 *              Roberto Marabini (roberto@csam.temple.edu)
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

#include <data/geometry.h>
#include <data/args.h>

void Usage();

int main(int argc, char **argv)
{
    Matrix1D<double> a, b;
    double crystal_size, /*taya,*/ rot, tilt, psi;
    string op;
    /*   bool   MRC;*/

    // Get parameters .......................................................
    try
    {
        a = get_vector_param(argc, argv, "-lattice_a");
        b = get_vector_param(argc, argv, "-lattice_b");
        crystal_size = AtoF(get_param(argc, argv, "-crystal_size", "1024"));
//      taya = AtoF(get_param(argc,argv,"-taya","0"));
        int i = position_param(argc, argv, "-ang");
        if (i == -1)
        {
            rot = tilt = psi = 0;
        }
        else
        {
            if (i + 3 >= argc) REPORT_ERROR(1, "Not enough parameters after -ang");
            rot  = AtoF(argv[i+1]);
            tilt = AtoF(argv[i+2]);
            psi  = AtoF(argv[i+3]);
        }
        op = get_param(argc, argv, "-op", "RF");
//      MRC = check_param(argc,argv,"-MRC");
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        Usage();
        exit(1);
    }

    // Process ..............................................................
    try
    {
        double gamma;
        matrix2D<double> E2D, vp(2, 2), Vp(2, 2), v(2, 2);
        double x, y;
        /* compute projection crystal vectors in projection coordinate system
           Solve system of equations:

            (a*)^t  a = 1, (a*)^t  b = 0
            (b*)^t  a = 0, (b*)^t  b = 1

            ( (a*)^t ) ( a b ) = I
            ( (b*)^t )

            Renaming

            V^t v = I

            Taking transpose

            v^t V = I

            Thus, V=(v^t)^(-1)

            NOTE 1: asume column vectors
            NOTE 2: a and b are the crystal vectors in real space
        */
#define SHOW(str,v,i) x=v(0,i); y=v(1,i); cout << str << "("<< x \
    << "," << y << ") (" << sqrt(x*x+y*y) << ")\n";
        if (op == "RF")
        {
            Euler_angles2matrix(rot, tilt, psi, E2D);
            E2D.resize(2, 2);
            v.setCol(0, a);
            v.setCol(1, b);
            vp = E2D * v;
            Vp = crystal_size * (vp.transpose()).inv();
            SHOW("a  =", v, 0);
            SHOW("b  =", v, 1);
            SHOW("ap =", vp, 0);
            SHOW("bp =", vp, 1);
            SHOW("ap*=", Vp, 0);
            SHOW("bp*=", Vp, 1);
            a = Vp.Col(0);
            a.resize(3);
            b = Vp.Col(1);
            b.resize(3);
            double gammap = RAD2DEG(atan2(vector_product(a, b).module(),
                                          dot_product(a, b)));
            cout << "Angle from ap* to bp*="
            << gammap << endl;
        }
        else if (op == "FR")
        {
            Euler_angles2matrix(rot, tilt, psi, E2D);
            E2D.resize(2, 2);
            Vp.setCol(0, a);
            Vp.setCol(1, b);
            vp = crystal_size * (Vp.transpose()).inv();
            v = E2D.inv() * vp;
            a = v.Col(0);
            a.resize(3);
            b = v.Col(1);
            b.resize(3);
            gamma = RAD2DEG(atan2(vector_product(a, b).module(),
                                  dot_product(a, b)));
            /*         if (MRC) {
                        rot=taya+(180-gamma)-rot;
                        Euler_angles2matrix(rot,tilt,psi,E2D); E2D.resize(2,2);
                 v=E2D.inv()*vp;
                     cout << "rot: " << rot <<endl;
                     cout << "tilt: " << tilt <<endl;
                     cout << "psi: " << psi <<endl;
                     cout << "gamma: " << rot <<endl;

                     }
            */
            SHOW("ap*=", Vp, 0);
            SHOW("bp*=", Vp, 1);
            SHOW("ap =", vp, 0);
            SHOW("bp =", vp, 1);
            SHOW("a  =", v, 0);
            SHOW("b  =", v, 1);
            cout << "Angle from a to b: " << gamma << endl;
        }
        else if (op == "PR")
        {
            Euler_angles2matrix(rot, tilt, psi, E2D);
            E2D.resize(2, 2);
            vp.setCol(0, a);
            vp.setCol(1, b);
            v = E2D.inv() * vp;
            SHOW("ap =", vp, 0);
            SHOW("bp =", vp, 1);
            SHOW("a  =", v, 0);
            SHOW("b  =", v, 1);
        }
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
    }
    exit(0);
}

void Usage()
{
    cout << "Purpose: This program allows you to compute the lattice vectors\n"
    << "         in the reciprocal space (Real-->Fourier & Fourier-->Real)\n"
    << "Usage: lattice_vectors [options]\n"
    << "  -lattice_a [xa,ya]             : First lattice vector\n"
    << "  -lattice_b [xb,yb]             : Second lattice vector\n"
    << " [-crystal_size <size=1024>      : Crystal size\n"
    << " [-ang <rot=0> <tilt=0> <psi=0>] : Projection angle\n"
    << " [-op <op=RF>]                   : RF -> Real space to Fourier space\n"
    << "                                   PR -> Projected space to Real space\n"
    << "                                   FR -> Fourier space to Real space\n"
//        << " [-MRC]                          : data coming from Spectra\n"
    ;
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Lattice_vectors {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Lattice_vectors/Help/lattice_vectors.html";
      help="Computes the lattice vectors of a crystal after projection
            in a certain direction";
      OPEN MENU menu_lattice_vectors;
      COMMAND LINES {
         +usual: xmipp_lattice_vectors
                    -lattice_a "["$XA","$YA"]"
                    -lattice_b "["$XB","$YB"]"
                   [-crystal_size $CRYSTAL_SIZE]
                   [-ang $ROT $TILT $PSI]
                   [-op $OP]
                   [-MRC]
      }
      PARAMETER DEFINITIONS {
         $XA   {label="X component of vector A"; type=float;}
         $YA   {label="Y component of vector A"; type=float;}
         $XB   {label="X component of vector B"; type=float;}
         $YB   {label="Y component of vector B"; type=float;}
         $CRYSTAL_SIZE {
            help="Size in pixels of the crystal image";
            label="Crystal size";
            type=natural;
            by default=1024;
         }
         OPT(-ang) {label="Projection direction";}
            $ROT  {label="Rotational angle"; type=float; by default=0;}
            $TILT {label="Tilting angle"; type=float; by default=0;}
            $PSI  {label="In-plane rotational angle"; type=float; by default=0;}
         OPT(-MRC) {label="Give results for MRC";}
         $OP {
            label="Operation";
            type=Exclusion {
               "Real -> Fourier" {RF}
               "Projection -> Real" {PR}
               "Fourier -> Real" {FR}
            };
         }
      }
   }

   MENU menu_lattice_vectors {
      "Lattice vectors"
      {$XA $YA}
      {$XB $YB}
      "Other parameters"
      OPT($CRYSTAL_SIZE)
      OPT($OP)
      OPT(-MRC)
      "Projection angles"
      OPT(-ang)
   }
*/
