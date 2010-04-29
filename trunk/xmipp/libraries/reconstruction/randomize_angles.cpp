/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

#include "randomize_angles.h"

/* Read parameters from command line. -------------------------------------- */
void RandomizeAngles::read(int argc, char **argv)
{
    fn_sel  = getParameter(argc, argv,  "-sel");
    fn_doc  = getParameter(argc, argv,  "-doc","");
    fh_out  = getParameter(argc, argv,  "-o");
    fn_sym  = getParameter(argc, argv,  "-sym","i3");
    verbose = checkParameter(argc, argv,"-verbose");
}

/* Show -------------------------------------------------------------------- */
void RandomizeAngles::show()
{
    std::cout   << "input selfile:             " << fn_sel  << std::endl
                << "output file:               " << fh_out << std::endl
                << "docfile :                  " << fn_doc  << std::endl
                << "outputfile dimensions:     " << dim << std::endl
                << "symmetry:                  " << fn_sym << std::endl
                << "verbose:                   " << verbose << std::endl
                ;
}

/* Usage ------------------------------------------------------------------- */
void RandomizeAngles::usage()
{
    std::cerr   << "   -sel filename             : Selfile with image names\n"
                << "   -o filename               : Name for output docfile\n"
                << "   -doc  filename  (optional): Auxiliary file with projection angles\n"
                << "   -sym symmetry flag (i3)   : symmetry description flag\n"
                << "   -verbose                  : set verbose mode on\n"
                ;
}

void RandomizeAngles::processAngles()
{
    double rot, tilt, psi, xoff, yoff, flip, weight,rotp,tiltp,psip;
    rot = tilt = psi = xoff = yoff = flip = weight = 0.;
    SF.go_beginning();
    Projection proj;
    int Ydim, Xdim;
    SF.ImgSize(Ydim, Xdim);
    DocFile DFout;
    DFout.append_comment("Headerinfo columns: rot (1) , tilt (2),\
 psi (3), Xoff (4), Yoff (5), Weight (6), Flip (7)");
    Matrix1D<double> docline;
    docline.initZeros(7);
    int repaint = ceil((double)SF.ImgNo()/60);

    init_progress_bar(SF.ImgNo());
    int imgno=0;//progress bar
    while (!SF.eof())
    {
        FileName fn_img = SF.NextImg();
        if (imgno++%repaint==0) progress_bar(imgno);
        if (fn_img=="") break;
        if (fn_doc == "")
        {
            proj.read(fn_img, false); //true means apply shifts
            rot  = proj.rot();
            tilt = proj.tilt();
            psi  = proj.psi();
            xoff = proj.Xoff();
            yoff = proj.Yoff();
            flip = proj.flip();
            weight = proj.weight();
        } else
        {
            get_angles_for_image(fn_img, rot, tilt, psi, xoff, yoff, flip, weight);
        }
        Matrix2D<double> euler(3, 3), temp;
        Euler_angles2matrix(rot, tilt, psi, euler);
        int irandom;
        irandom=rnd_unif(0, SL.SymsNo()+1);
        //#define DEBUG
        #ifdef DEBUG
            for (int h=0; h < 100000; h++)
            {
                irandom=rnd_unif(0, SL.SymsNo()+1);
                fprintf(stderr,"%02d\n",irandom);
            }
            exit(1);
        #endif
        #undef DEBUG
        temp = euler *
               R_repository[irandom].inv();
        //temp = euler;
        Euler_matrix2angles(temp, rotp, tiltp, psip);
        docline(0) = rotp;
        docline(1) = tiltp;
        docline(2) = psip;
        docline(3) = xoff;
        docline(4) = yoff;
        docline(5) = weight;
        docline(6) = flip;
        DFout.append_comment(fn_img);
        DFout.append_data_line(docline);
    }
    progress_bar(SF.ImgNo());
    DFout.write(fh_out);
}

void RandomizeAngles::get_angles_for_image(const FileName &fn, double &rot,
    double &tilt, double &psi, double &xoff, double &yoff, double &flip,
    double &weight)
{
    if (DFangles.search_comment(fn))
    {
        rot    = DFangles(col_rot);
        tilt   = DFangles(col_tilt);
        psi    = DFangles(col_psi);
        xoff   = DFangles(col_xoff);
        yoff   = DFangles(col_yoff);
        if (col_flip < 0)
            flip   = 0.;
        else
            flip   = DFangles(col_flip);
        if (col_weight < 0)
            weight = 0.;
        else
            weight = DFangles(col_weight);
    }
    else
    {
        REPORT_ERROR(1, (std::string)"RandomizeAnles: Cannot find " + fn + " in docfile " + fn_doc);
    }
    
}
void RandomizeAngles::symmetryInicialization()
{
    /** vector with symmetry matrices */
    Matrix2D<double>  R(4, 4),L(4,4);
    Matrix2D<double>  Identity(3,3);
    Identity.initIdentity();
    //symmetryMatrixVertex.resize(12, 5);
    R_repository.push_back(Identity);
    for (int isym = 0; isym < SL.SymsNo(); isym++)
    {
        SL.get_matrices(isym, L, R);
        R.resize(3, 3);
        R_repository.push_back(R);
    }

}
/* Main program ------------------------------------------------------------ */
void RandomizeAngles::run()
{
    double accuracy=1e-6;
    randomize_random_generator();
    show();
    //read doc and sel files file
    if (fn_doc != "")
        DFangles.read(fn_doc);
        col_rot    = DFangles.getColNumberFromHeader("rot")    - 1;
        col_tilt   = DFangles.getColNumberFromHeader("tilt")   - 1;
        col_psi    = DFangles.getColNumberFromHeader("psi")    - 1;
        col_xoff   = DFangles.getColNumberFromHeader("Xoff")   - 1;
        col_yoff   = DFangles.getColNumberFromHeader("Yoff")   - 1;
        col_flip   = DFangles.getColNumberFromHeader("Flip")   - 1;
        col_weight = DFangles.getColNumberFromHeader("Weight") - 1;

    SF.read(fn_sel);
    
    //load icosahedron vertex
    symmetry=SL.read_sym_file(fn_sym, accuracy);
    symmetryInicialization();
    processAngles();
}


