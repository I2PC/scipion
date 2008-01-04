/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.uam.es (1999)
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

#include <data/args.h>
#include <data/selfile.h>
#include <data/geometry.h>
#include <data/histogram.h>
#include <interface/spider.h>
#include <interface/opendxang.h>

#include <fstream>

void Usage();

int main(int argc, char *argv[])
{
    std::string     ang1 = "rot", ang2 = "tilt", ang3 = "psi";
    DocFile         angles;
    FileName        fn_ang, fn_sel, fn_hist, fn_ps, fn_DX;
    int             steps;
    int             tell;
    float           R, r, rmax, wmax = -99.e99;
    float           rot_view;
    float           tilt_view;
    int             up_down_correction, colw;
    bool            solid_sphere;


// Check the command line ==================================================
    try
    {
        fn_sel = getParameter(argc, argv, "-sel", "");
        fn_ang = getParameter(argc, argv, "-ang", "");
        fn_hist = getParameter(argc, argv, "-hist", "");
        fn_ps = getParameter(argc, argv, "-ps", "");
        fn_DX = getParameter(argc, argv, "-DX", "");
        steps = textToInteger(getParameter(argc, argv, "-steps", "100"));
        tell = checkParameter(argc, argv, "-show_process");
        R = textToFloat(getParameter(argc, argv, "-R", "60"));
        rmax = textToFloat(getParameter(argc, argv, "-r", "1.5"));
        rot_view = textToFloat(getParameter(argc, argv, "-rot_view",  "0"));
        tilt_view = textToFloat(getParameter(argc, argv, "-tilt_view", "30"));
        up_down_correction = checkParameter(argc, argv, "-up_down_correction");
        solid_sphere = checkParameter(argc, argv, "-solid_sphere");
        colw = textToInteger(getParameter(argc, argv, "-wcol", "-1"));

        // Angle order
        int i;
        if ((i = paremeterPosition(argc, argv, "-order")) != -1)
        {
            if (i + 3 >= argc)
            {
                std::cout << "Angular distribution: Not enough parameters behind -ang\n";
                Usage();
                exit(1);
            }
            ang1 = argv[i+1];
            ang2 = argv[i+2];
            ang3 = argv[i+3];
        }

        // Check they are "rot", "tilt", and "psi"
        checkAngle(ang1);
        checkAngle(ang2);
        checkAngle(ang3);
        if (ang1[1] == ang2[1] || ang1[1] == ang3[1] || ang2[1] == ang3[1])
            REPORT_ERROR(1, "Angular distribution: There is an angle twice in the angle order");

        // Check there is some input
        if (fn_ang == "" && fn_sel == "")
            REPORT_ERROR(1, "Angular distribution: There is no input information");
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage();
        exit(1);
    }

    try
    {
// Get angles ==============================================================
        if (fn_ang != "")
            angles.read(fn_ang);
        else
        {
            SelFile selfile(fn_sel);
            extract_angles(selfile, angles, ang1, ang2, ang3);
        }
        int AngleNo = angles.dataLineNo();
        if (AngleNo == 0)
            EXIT_ERROR(1, "Angular distribution: Input files doesn't contain angular information");

        if (colw >= 0)
        {
            // Find maximum weight
            for (int i = 0; i < AngleNo; i++) if (angles(i + 1, colw) > wmax) wmax = angles(i + 1, colw);
        }

// Build vector tables ======================================================
#define GET_ANGLES(i) \
    angles.get_angles(i,rot,tilt,psi,ang1,ang2,ang3); \
    if (up_down_correction && ABS(tilt)>90) \
        Euler_up_down(rot,tilt,psi,rot,tilt,psi);

        double rot, tilt, psi;
        std::vector< Matrix1D<double> > v, v_ang;
        v.reserve(AngleNo);
        v_ang.reserve(AngleNo);
        for (int i = 0; i < AngleNo; i++)
        {
            Matrix1D<double> aux(3);
            Matrix1D<double> aux_ang(6);


            GET_ANGLES(i + 1);
            Euler_direction(rot, tilt, psi, aux);
            v.push_back(aux);


            aux_ang = vectorR3(rot, tilt, psi);
            v_ang.push_back(aux_ang);

        }

//Show distribution with OpenDx ==============================================
        openDXang DX;
        DX.openDXangFile(fn_DX);

        for (int i = 0; i < AngleNo; i++)
        {


            DX.Add_Item(v_ang[i]);
        }


// Compute histogram of distances =============================================
        if (fn_hist != "")
        {
            Matrix1D<double> dist;

#define di VEC_ELEM(dist,i)
#define dj VEC_ELEM(dist,j)

#define SHOW {\
        GET_ANGLES(i+1); \
        std::cout << i << " " << rot << " " << tilt << " v[i]=" \
        << v[i].transpose() << std::endl; \
        GET_ANGLES(j+1); \
        std::cout << j << " " << rot << " " << tilt << " v[j]=" \
        << v[j].transpose() << std::endl; \
        std::cout << " d= " << d << std::endl << std::endl; \
    }

            // Compute minimum distance table
            dist.initZeros(AngleNo);
            for (int i = 0; i < AngleNo; i++)
                for (int j = i + 1; j < AngleNo; j++)
                {
                    double d = spherical_distance(v[i], v[j]);
                    if (di == 0 || d < di)
                    {
                        di = d;
                        if (tell) SHOW;
                    }
                    if (dj == 0 || d < dj)
                    {
                        dj = d;
                        if (tell) SHOW;
                    }
                }

            histogram1D dist_hist;
            double min, max;
            dist.computeDoubleMinMax(min, max);
            dist_hist.init(min, max, steps);
            for (int i = 0; i < AngleNo; i++) dist_hist.insert_value(di);
            dist_hist.write(fn_hist);
        }



// Show distribution as triangles ==========================================
        if (fn_ps != "")
        {
            std::ofstream fh_ps;
            fh_ps.open(fn_ps.c_str(), std::ios::out);
            if (!fh_ps)
                EXIT_ERROR(1, (std::string)"Ang_distribution: Cannot open " + fn_ps + " for output");

            fh_ps << "%%!PS-Adobe-2.0\n";
            fh_ps << "%% Creator: Angular Distribution\n";
            fh_ps << "%% Title: Angular distribution of " << fn_sel << "\n";
            fh_ps << "%% Pages: 1\n";

#define TO_PS(x,y) \
    tmp=y; \
    y=400.0f-x*250.0f/60; \
    x=300.0f+tmp*250.0f/60;

            Matrix1D<double> p0(4), p1(4), p2(4), p3(4), origin(3);
            Matrix2D<double> A, euler_view;
            Euler_angles2matrix(rot_view, tilt_view, 0.0f, euler_view);
            origin.initZeros();
            double tmp;
            for (int i = 0; i < AngleNo; i++)
            {


                // Triangle size depedent on w
                if (colw >= 0)
                {
                    r = angles(i + 1, colw);
                    r *= rmax / wmax;
                }
                else r = rmax;


                // Initially the triangle is on the floor of the projection plane
                VECTOR_R3(p0,    0   ,      0        , 0);
                VECTOR_R3(p1,    0   , r*2 / 3*SIND(60), 0);
                VECTOR_R3(p2, r / 2*0.6, -r*1 / 3*SIND(60), 0);
                VECTOR_R3(p3, -r / 2*0.6, -r*1 / 3*SIND(60), 0);

                // Convert to homogeneous coordinates
                p0(3) = 1;
                p1(3) = 1;
                p2(3) = 1;
                p3(3) = 1;

                // Compute Transformation matrix
                GET_ANGLES(i + 1);
                Euler_angles2matrix(rot, tilt, psi, A);

                // We go from the projeciton plane to the universal coordinates
                A = A.transpose();

                // Convert to homogeneous coordinates and apply a translation
                // to the sphere of radius R
                A.resize(4, 4);
                A(0, 3) = R * XX(v[i]);
                A(1, 3) = R * YY(v[i]);
                A(2, 3) = R * ZZ(v[i]);
                A(3, 3) = 1;

                // Convert triangle coordinates to universal ones
                p0 = A * p0;
                p1 = A * p1;
                p2 = A * p2;
                p3 = A * p3;

                // Check if this triangle must be drawn
                if (solid_sphere)
                {
                    Matrix1D<double> view_direction, p0p;
                    euler_view.getRow(2, view_direction);
                    p0p = p0;
                    p0p.resize(3);
                    if (point_plane_distance_3D(p0, origin, view_direction) < 0)
                        continue;
                }

                // Project this triangle onto the view plane and write in PS
                Matrix1D<double> pp(3);
                Uproject_to_plane(p1, euler_view, pp);
                TO_PS(XX(pp), YY(pp));
                fh_ps << "newpath\n";
                fh_ps << XX(pp) << " " << YY(pp) << " moveto\n";

                Uproject_to_plane(p2, euler_view, pp);
                TO_PS(XX(pp), YY(pp));
                fh_ps << XX(pp) << " " << YY(pp) << " lineto\n";

                Uproject_to_plane(p3, euler_view, pp);
                TO_PS(XX(pp), YY(pp));
                fh_ps << XX(pp) << " " << YY(pp) << " lineto\n";

                Uproject_to_plane(p1, euler_view, pp);
                TO_PS(XX(pp), YY(pp));
                fh_ps << XX(pp) << " " << YY(pp) << " lineto\n";

                fh_ps << "closepath\nstroke\n";
            }
            fh_ps << "showpage\n";
            fh_ps.close();
        }

    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }

}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    std::cout << "Usage:\n";
    std::cout << "   ang_distribution <options>\n";
    std::cout << "   Where <options> are:\n";
    std::cout << "      (-sel <sel_file> |            : selection file with the set of images\n"
    << "       -ang <ang_file>              : Spider document file with the angles\n"
    << "      [-order <ang1> <ang2> <ang3>]): where ang1, ang2 and ang3 are\n"
    << "                                      either psi, psi, tilt, tilt,\n"
    << "                                      rot or rot. The default\n"
    << "                                      order is rot, tilt and psi.\n"
    << "      [-DX <-DX file out>           : DX file\n"
    << "      [-hist <doc_file>]            : histogram of distances\n"
    << "      [-steps <stepno=100>]         : number of divisions in the histogram\n"
    << "      [-show_process]               : show distances.\n"
    << "      [-ps <PS file out>]           : PS file with the topological sphere\n"
    << "      [-R <big_sphere_radius=60>]   : sphere radius for the PS file\n"
    << "      [-r <triangle side=1.5>]      : triangle size for the PS file\n"
    << "      [-rot_view <rot angle=0>]     : rotational angle for the view\n"
    << "      [-tilt_view <tilt angle=30>]  : tilting angle for the view\n"
    << "      [-wcol <column number=-1>]    : generate triangles with size depending on \n"
    << "                                      number in corresponding column of the docfile\n"
    << "      [-up_down_correction]         : correct angles so that a semisphere\n"
    << "                                      is shown\n"
    << "      [-solid_sphere]               : projections in the back plane are\n"
    << "                                      not shown\n"
    ;
}

/* ------------------------------------------------------------------------- */
/* Menus                                                                     */
/* ------------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Ang_distribution {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Ang_distribution/Help/ang_distribution.html";
      help="Show angular distribution";
      OPEN MENU Ang_distribution;
      COMMAND LINES {
         + Angle_Set:
             ang_distribution -ang $ANGFILE
                [-ang $ANG1 $ANG2 $ANG3] [$ANGL1 $ANGL2 $ANGL3]
                [-hist $HISTFILE [-steps $STEPNO] [-show_process]]
                [-ps $PSFILE [-R $RADIUS] [-r $SIDE] [-rot_view $ROT]
                   [-tilt_view $TILT] [-up_down_correction]]
         + SelFile:
             ang_distribution -sel $SELFILE OPT(-hist) OPT(-ps)
      }
      PARAMETER DEFINITIONS {
         $ANGFILE {
            label="DocFile with angles";
            type=FILE EXISTING;
         }
         $SELFILE {
            label="Selection file";
            type=FILE EXISTING;
         }
         #include "angles.mnu"
         OPT(-hist) {label="Generate Distance Histogram";}
            $HISTFILE {
               label="Output histogram file";
               type=FILE;
            }
            $STEPNO {
               label="Histogram steps";
               help="The histogram is divided in this number of steps";
               type=NATURAL;
               by default=100;
            }
            OPT(-show_process) {label="Show process information";}
         OPT(-ps) {label="Generate Topological sphere";}
            $PSFILE {
               label="OutputPostscript file";
               type=FILE;
            }
            $RADIUS {
               label="Topological sphere radius";
               type=FLOAT [0...];
               by default=60;
            }
            $SIDE {
               label="Triangle side length";
               type=FLOAT [0...];
               by default=1.5;
            }
            $ROT {
               label="Rotational angle for Point of view";
               help="degrees";
               type=FLOAT [0...360];
               by default=0;
            }
            $TILT {
               label="Tilting angle for Point of view";
               help="degrees";
               type=FLOAT [0...360];
               by default=30;
            }
            OPT(-up_down_correction) {label= "Up-Down correction";}
      }
   }

   MENU Ang_distribution {
      "I/O Parameters"
      $ANGFILE
      OPT($ANGL1)
      $SELFILE
      "Distance Histogram Parameters"
      OPT(-hist)
      "Topological Sphere Parameters"
      OPT(-ps)
   }
*/
