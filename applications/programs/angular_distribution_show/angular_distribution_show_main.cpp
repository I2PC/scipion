/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (1999)
 *             Roberto Marabini        added bild option (2008)
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
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <data/args.h>
#include <data/geometry.h>
#include <data/histogram.h>
#include <data/metadata.h>
#include <data/xmipp_program.h>
#include <interface/spider.h>
#include <fstream>

class ProgAngularDistributionShow: public XmippProgram
{
public:
    FileName fnIn, fnOut;
    String type;
    double R, rmax;
    int shift_center, solid_sphere, steps;
    double rot_view, tilt_view;
    bool up_down_correction;

    void readParams()
    {
        fnIn = getParam("-i");
        fnOut = getParam("-o");
        type = getParam("-o",1);
        if (type=="chimera")
        {
            R=getDoubleParam("-o",2);
            rmax=getDoubleParam("-o",3);
            shift_center=getIntParam("-o",4);
        }
        else if (type=="ps")
        {
            R=getDoubleParam("-o",2);
            rmax=getDoubleParam("-o",3);
            rot_view=getDoubleParam("-o",4);
            tilt_view=getDoubleParam("-o",5);
            solid_sphere=getIntParam("-o",6);
        }
        else if (type=="histogram")
            steps = getIntParam("-o",2);
        up_down_correction = checkParam("--up_down_correction");
    }

    void defineParams()
    {
        addUsageLine("Shows which is the angular distribution of a set of projections.");
        addUsageLine("+There are three kinds of outputs: a bild file (chimera understands this format), a distance histogram and a postscript file. ");
        addUsageLine("+In postcript, each projection is represented by a small isosceles triangles (the longest part is pointing to Y in the projection plane). ");
        addUsageLine("+In chimera, each projection direction has a sphere whose radius is proportional to the number of images assigned to it");
        addUsageLine("+The histogram output is the histogram of the minimum distance (in degrees) among projections. The distance is measured on the sphere surface. This gives an idea of how compact is the angular distribution.");
        addParamsLine(" -i <metadata>                 : Metadata file with the angles");
        addParamsLine(" -o <file> <type>              : Output file");
        addParamsLine("    where <type>");
        addParamsLine("          chimera <R=60> <rmax=1.5> <shift_center=0>               : R=sphere radius, rmax=maximum point radius, shift_center=shift in pixels applied to all coordinates");
        addParamsLine("                                                                   : shift_center should be half the volume size and R should be 2/3 the volume size");
        addParamsLine("          ps <R=60> <rmax=1.5> <rot=0> <tilt=30> <solid_sphere=0>  : R=sphere radius, rmax=maximum point radius, rot and tilt defines the point of view, solid_sphere=0 (=no) or 1 (=yes)");
        addParamsLine("                                                                : If the sphere is solid, projections in the back plane are not shown");
        addParamsLine("          histogram <stepno=100> : Number of divisions in the histogram");
        addParamsLine("[--up_down_correction]         : correct angles so that a semisphere is shown");
    }

    void run()
    {
        // Get angles ==============================================================
        MetaData angles;
        angles.read(fnIn);
        size_t AngleNo = angles.size();
        if (AngleNo == 0 || !angles.containsLabel(MDL_ANGLE_ROT))
            REPORT_ERROR(ERR_MD_BADLABEL, "Input file doesn't contain angular information");

        double maxWeight = -99.e99;
        MultidimArray<double> weight;
        weight.initZeros(AngleNo);
        if (angles.containsLabel(MDL_WEIGHT))
        {
            // Find maximum weight
            int i=0;
            FOR_ALL_OBJECTS_IN_METADATA(angles)
            {
                double w;
                angles.getValue(MDL_WEIGHT,w,__iter.objId);
                DIRECT_A1D_ELEM(weight,i++)=w;
                maxWeight=XMIPP_MAX(w,maxWeight);
            }
        }
        double imaxWeight=1.0/maxWeight;

        // Build vector tables ======================================================
        std::vector< Matrix1D<double> > v, v_ang;
        v.reserve(AngleNo);
        v_ang.reserve(AngleNo);
        Matrix1D<double> aux(3);
        Matrix1D<double> aux_ang(3);
        FOR_ALL_OBJECTS_IN_METADATA(angles)
        {
            double rot, tilt, psi;
            angles.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
            angles.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
            angles.getValue(MDL_ANGLE_PSI,psi,__iter.objId);
            if (up_down_correction && fabs(tilt)>90)
                \
                Euler_up_down(rot,tilt,psi,rot,tilt,psi);
            Euler_direction(rot, tilt, psi, aux);
            v.push_back(aux);
            VECTOR_R3(aux_ang, rot, tilt, psi);
            v_ang.push_back(aux_ang);
        }

        if (type=="histogram")
        {
            // Compute minimum distance table
            MultidimArray<double> dist;
            dist.initZeros(AngleNo);
            for (size_t i = 0; i < AngleNo; i++)
            {
                const Matrix1D<double> &vi=v[i];
                for (size_t j = i + 1; j < AngleNo; j++)
                {
                    const Matrix1D<double> &vj=v[j];
                    // Since the two vectors are in the unit sphere, the spherical distance
                    // is simply the arc cosine
                    double d = acos(XX(vi)*XX(vj) + YY(vi)*YY(vj) + ZZ(vi)*ZZ(vj));
                    d=RAD2DEG(d);
                    if (DIRECT_A1D_ELEM(dist,i) == 0 || d < DIRECT_A1D_ELEM(dist,i))
                    	DIRECT_A1D_ELEM(dist,i) = d;
                    if (DIRECT_A1D_ELEM(dist,j) == 0 || d < DIRECT_A1D_ELEM(dist,j))
                    	DIRECT_A1D_ELEM(dist,j) = d;
                }
            }

            Histogram1D dist_hist;
            compute_hist(dist, dist_hist, steps);
            dist_hist.write(fnOut);
        }
        else if (type=="chimera")
        {
            std::ofstream fh_bild;
            fh_bild.open(fnOut.c_str(), std::ios::out);
            if (!fh_bild)
                REPORT_ERROR(ERR_IO_NOWRITE, fnOut);
            fh_bild << ".color 1 0 0" << std::endl;

            int imax=v.size();
            for (int i=0; i<imax; i++)
            {
                double r=rmax;
                if (maxWeight>0)
                    r *= DIRECT_A1D_ELEM(weight,i)*imaxWeight;
                fh_bild
                << ".sphere "
                << R*XX(v[i])  + shift_center << " "
                << R*YY(v[i])  + shift_center << " "
                << R*ZZ(v[i])  + shift_center << " "
                << r
                <<"\n";
            }
            fh_bild.close();
        }
        else if (type=="ps")
        {
            std::ofstream fh_ps;
            fh_ps.open(fnOut.c_str(), std::ios::out);
            if (!fh_ps)
                REPORT_ERROR(ERR_IO_NOWRITE, fnOut);

            fh_ps << "%%!PS-Adobe-2.0\n";
            fh_ps << "%% Creator: Angular Distribution\n";
            fh_ps << "%% Title: Angular distribution of " << fnIn << "\n";
            fh_ps << "%% Pages: 1\n";

#define TO_PS(x,y) \
        tmp=y; \
        y=400.0f-x*250.0f/60; \
        x=300.0f+tmp*250.0f/60;

            Matrix1D<double> p0(4), p1(4), p2(4), p3(4), view_direction, pp(3);
            Matrix2D<double> A, euler_view;
            Euler_angles2matrix(rot_view, tilt_view, 0.0f, euler_view);
            euler_view.getRow(2, view_direction);
            double tmp;
            int imax=v.size();
            for (int i=0; i<imax; i++)
            {
                double r=rmax;
                if (maxWeight>0)
                    r *= DIRECT_A1D_ELEM(weight,i)*imaxWeight;

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
                const Matrix1D<double> &v_angi=v_ang[i];
                const Matrix1D<double> &vi=v[i];
                Euler_angles2matrix(VEC_ELEM(v_angi,1), VEC_ELEM(v_angi,2),
                                    VEC_ELEM(v_angi,3), A, true);
                A = A.transpose(); // We go from the projection plane to the universal coordinates

                // Apply a translation to the sphere of radius R
                MAT_ELEM(A, 0, 3) = R * XX(vi);
                MAT_ELEM(A, 1, 3) = R * YY(vi);
                MAT_ELEM(A, 2, 3) = R * ZZ(vi);

                // Convert triangle coordinates to universal ones
                p0 = A * p0;
                p1 = A * p1;
                p2 = A * p2;
                p3 = A * p3;

                // Check if this triangle must be drawn
                if (solid_sphere)
                {
                	// Point-plane distance
                	double d=XX(p0)*XX(view_direction)+
                			 YY(p0)*YY(view_direction)+
                			 ZZ(p0)*ZZ(view_direction);
                    if (d < 0)
                        continue;
                }

                // Project this triangle onto the view plane and write in PS
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
};

int main(int argc, char *argv[])
{
    ProgAngularDistributionShow prm;
    prm.read(argc,argv);
    return prm.tryRun();
}
