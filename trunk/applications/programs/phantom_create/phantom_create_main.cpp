/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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

#include <data/phantom.h>
#include <data/xmipp_program.h>

class ProgPhantomCreate: public XmippProgram
{
protected:
    FileName          fn_phantom;
    FileName          fn_vol;
    Phantom           phantom;
    Image<double>     vol;

    void defineParams()
    {
        addUsageLine("Create phantom volume from a feature description file.");
        addUsageLine("+You may define a mathematical phantom from its geometrical features");
        addUsageLine("+(cubes, cylinders, ...) and translate it into a voxel volume with this program");
        addUsageLine("+The structure of the description file is:");
        addUsageLine("+# General Volume Parameters:",true);
        addUsageLine("+#      Xdim      Ydim      Zdim   Background_Density   Scale",true);
        addUsageLine("+      <Xdim>    <Ydim>    <Zdim>     <back_den=0>    <scale=1>",true);
        addUsageLine("The *scale* is by default 1. This scale factor is a number which will multiply all features in order to reduce them (scale<1) or enlarge them (scale>1)");
        addUsageLine("+# Feature Parameters:",true);
        addUsageLine("+Sphere, with center at (x0,y0,z0) and radius in pixels [radius]");
        addUsageLine("+sph +/- den x0 y0 z0 radius",true);
        addUsageLine("+Blob with center at (x0,y0,z0) and radius in pixels [radius], alpha is the tampering factor and order is the order of the Bessel functions (usually 2)");
        addUsageLine("+blo +/- den x0 y0 z0 radius alpha order",true);
        addUsageLine("+Gaussian with center at (x0,y0,z0) and sigma in pixels [sigma]");
        addUsageLine("+gau +/- den x0 y0 z0 sigma",true);
        addUsageLine("+Cylinder, initially with base at plane XY and height in Z (-h/2 to h/2) and then it is moved to (x0,y0,z0) (its center) and rotated after (tilt,rot,psi)");
        addUsageLine("+cyl +/- den x0 y0 z0 xradius yradius height rot tilt psi",true);
        addUsageLine("+double cylinder (two cylinders), initially with base at plane XY and height in Z (-h/2 to h/2), then it is moved to (x0,y0,z0) and rotated after (tilt,rot,psi)");
        addUsageLine("+dcy +/- den x0 y0 z0 radius height separation rot tilt psi",true);
        addUsageLine("+Cube of size (xdim,ydim,zdim) whose center is moved to (x0,y0,z0) and rotated");
        addUsageLine("+cub +/- den x0 y0 z0 xdim ydim zdim rot tilt psi",true);
        addUsageLine("+Ellipsoid with 3 radius and then moved and rotated");
        addUsageLine("+ell +/- den x0 y0 z0 xradius yradius zradius rot tilt psi",true);
        addUsageLine("+Cone, initially with base at plane XY and height in Z (-h/2 to h/2), and then the (0,0,0) point is moved to (x0,y0,z0) and rotated");
        addUsageLine("+con +/- den x0 y0 z0 radius height rot tilt psi",true);
        addUsageLine("+");
        addUsageLine("+*Overlapping* : features can overlap and add their densities at overlapping positions, use + to allow the sum of overlapping densities or = if you want that the density of the latest feature prevail at a given position");
        addUsageLine("+");
        addUsageLine("+*Index symmetry* : due to the Xmipp definition of center, volumes are symmetrical (with respect to the size, ie, are defined between =-xi= and =+xi=) for odd dimensions (for instance, 65x65x65, in this case the volume is defined between =[-32,+32]= for the three directions)");
        addUsageLine("+");
        addUsageLine("+*Euler angles* : Euler angles cause figures to be rotated, here you have some useful Euler combinations to make a cone (this figure has been selected as it seems to point in some direction) point in the given direction");
        addUsageLine("+|  *Direction*  |  *rot*  |  *tilt*  |  *psi*  |");
        addUsageLine("+|  -X  |  180  |  90  |  0  |");
        addUsageLine("+|  +X  |  0  |  90  |  0  |");
        addUsageLine("+|  -Y  |  270  |  90  |  0  |");
        addUsageLine("+|  +Y  |  90  |  90  |  0  |");
        addUsageLine("+|  -Z  |  0  |  180  |  0  |");
        addUsageLine("+|  +Z  |  0  |  0  |  0  |");
        addUsageLine("+");
        addUsageLine("+*Graphical design* : There is a graphical phantom [[http://sites.google.com/site/phan3d][designer]] which was published in this [[http://www.ncbi.nlm.nih.gov/pubmed/15217810][article]]. The program needs OpenGL libraries installed.");
        addUsageLine("+");
        addUsageLine("+Here you have an example file with 12 balls in two rings:");
        addUsageLine("<pre>");
        addUsageLine("# Phantom description file, (generated with phantom help)");
        addUsageLine("# General Volume Parameters: ");
        addUsageLine("#      Xdim      Ydim      Zdim   Background");
        addUsageLine("        65        65        65        0");
        addUsageLine("# Feature Parameters:");
        addUsageLine("#Type +/= Density X_Center Y_Center Z_Center Radius");
        addUsageLine("# Bottom ring -----------------------------------------------------");
        addUsageLine(" sph  +      1      15         0      -8.5      7.5");
        addUsageLine(" sph  +      1       7.5      13      -8.5      7.5");
        addUsageLine(" sph  +      1       7.5     -13      -8.5      7.5");
        addUsageLine(" sph  +      1      -7.5      13      -8.5      7.5");
        addUsageLine(" sph  +      1      -7.5     -13      -8.5      7.5");
        addUsageLine(" sph  +      1     -15         0      -8.5      7.5");
        addUsageLine("# Top ring --------------------------------------------------------");
        addUsageLine(" sph  +      1      15         0       8.5      7.5");
        addUsageLine(" sph  +      1       7.5      13       8.5      7.5");
        addUsageLine(" sph  +      1       7.5     -13       8.5      7.5");
        addUsageLine(" sph  +      1      -7.5      13       8.5      7.5");
        addUsageLine(" sph  +      1      -7.5     -13       8.5      7.5");
        addUsageLine(" sph  +      1     -15         0       8.5      7.5");
        addUsageLine("</pre>");
        addParamsLine("-i <description_file> : Input file with the mathematical features");
        addParamsLine("-o <output_file>      : Output volume in voxels");
        addExampleLine("xmipp_phantom_create -i phantom.descr -o phantom.vol");
    }

    void readParams()
    {
        fn_phantom = getParam("-i");
        fn_vol = getParam("-o");
    }

public:
    void run()
    {
        phantom.read(fn_phantom);
        phantom.draw_in(vol());
        vol.write(fn_vol);
    }
}
;//end of class ProgPhantomCreate

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgPhantomCreate program;
    program.read(argc, argv);
    return program.tryRun();
}
