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

#include "virus_vertex.h"
#include <data/metadata_extension.h>

/* Read parameters from command line. -------------------------------------- */
void VirusVertex::read(int argc, char **argv)
{
    fn_doc  = getParameter(argc, argv,  "-i");
    fn_root = getParameter(argc, argv,  "-root");
    virusRadius = textToFloat(getParameter(argc, argv, "-radius"));
    minVirusRadius = textToFloat(getParameter(argc, argv, "-min_radius"));
    dim = textToInteger(getParameter(argc, argv, "-dim"));
    fn_sym  = getParameter(argc, argv, "-sym","i3");
    verbose = checkParameter(argc, argv,"-verbose");
    removeCloseVertex = checkParameter(argc, argv,"-remove");
}

/* Show -------------------------------------------------------------------- */
void VirusVertex::show()
{
    std::cout
    << "input docfile:             " << fn_doc  << std::endl
    << "output files root:         " << fn_root << std::endl
    << "virus radius:              " << virusRadius << std::endl
    << "minimun virus radius:      " << minVirusRadius << std::endl
    << "outputfile dimensions:     " << dim << std::endl
    << "symmetry:                  " << fn_sym << std::endl
    << "verbose:                   " << verbose << std::endl
    << "remove:                    " << removeCloseVertex<< std::endl
    ;
}

/* Usage ------------------------------------------------------------------- */
void VirusVertex::usage()
{
    std::cerr
    << "   -i  filename              : List of images, and optionally projection angles\n"
    << "   -root filename            : Root name for output images\n"
    << "   -radius float_number      : virus radius (distance from center\n"
    << "                               to vertex in pixels)\n"
    << "   -min_radius  float_number : extract vertex located at a radius greater\n"
    << "                               than min_radius\n"
    << "   -dim integer              : size of the extracted images (pixels)\n"
    << "   -sym symmetry flag (i3)   : symmetry description flag\n"
    << "   -verbose                  : set verbose mode on\n"
    << "   -remove                   : remove vertexs if closer than 'dim'\n"
    ;
}
void VirusVertex::loadIcosahedronVertex()
{
    vertices_vectors.push_back(vectorR3(0., 0., 1.));
    vertices_vectors.push_back(vectorR3(0.723606900230461, -0.525731185781806, 0.447213343087301));
    vertices_vectors.push_back(vectorR3(0.723606900230461, 0.525731185781806, 0.447213343087301));
    vertices_vectors.push_back(vectorR3(-0.276393239417711, 0.850650928976665, 0.447213343087301));
    vertices_vectors.push_back(vectorR3(-0.8944273172062, 0., 0.447213343087301));
    vertices_vectors.push_back(vectorR3(-0.276393239417711, -0.850650928976665, 0.447213343087301));
    vertices_vectors.push_back(vectorR3(0.8944273172062, 0., -0.447213343087301));
    vertices_vectors.push_back(vectorR3(0.276393242471372, 0.850650927984471, -0.447213343087301));
    vertices_vectors.push_back(vectorR3(-0.723606898343194, 0.525731188379405, -0.447213343087301));
    vertices_vectors.push_back(vectorR3(-0.723606898343194, -0.525731188379405, -0.447213343087301));
    vertices_vectors.push_back(vectorR3(0.276393242471372, -0.850650927984471, -0.447213343087301));
    vertices_vectors.push_back(vectorR3(0., 0., -1.));

    Matrix2D<double>  A(3, 3);
    if (symmetry  == pg_I || symmetry  == pg_I2)
    {
        Euler_angles2matrix(0,-31.7174745559,0, A);
    }
    else if (symmetry  == pg_I1)
    {//OK
        Euler_angles2matrix(0,-31.7174745559+90.,0, A);
    }
    else if (symmetry  == pg_I3)
    {//OK
        A.initIdentity();
    }
    else if (symmetry  == pg_I4)
    {//OK
        Euler_angles2matrix(0,-31.7174745559 *2.0,0, A);
    }
    else if (symmetry  == pg_I5)
    {//OK
        std::cerr << "ERROR: Symmetry pg_I5 not implemented" << std::endl;
        exit(0);
    }
    else
    {//OK
        std::cerr << "Unknown symmetry" << std::endl;
        exit(0);
    }
    for (int i = 0; i < 12; i++)
    {
        vertices_vectors[i]=A*vertices_vectors[i];
    }// for i
#define CHIMERA
    #ifdef CHIMERA
    std::ofstream filestr;
    filestr.open ("sym.bild");
    for (int i = 0; i < 12; i++)
    {
        filestr    << ".color red"
        << std::endl
        << ".sphere "
        << XX(vertices_vectors[i]) << " "
        << YY(vertices_vectors[i]) << " "
        << ZZ(vertices_vectors[i]) << " "
        << " .05"
        << std::endl
        ;
    }
    filestr.close();
#endif
    #undef CHIMERA
}

void VirusVertex::processAngles()
{
    double rot, tilt, psi, xoff, yoff, flip, weight,rotp,tiltp,psip;
    rot = tilt = psi = xoff = yoff = flip = weight = 0.;
    Projection proj,proj_aux;
    int Zdim, Ydim, Xdim;
    size_t Ndim;
    getImageSize(DF, Xdim, Ydim, Zdim, Ndim);
    MetaData DFout;
    int repaint = ceil((double)DF.size()/60);

    init_progress_bar(DF.size());
    int imgno=0;
    FileName fn_img;
    FOR_ALL_OBJECTS_IN_METADATA(DF)
    {
        DF.getValue(MDL_IMAGE,fn_img,__iter.objId);
        if (imgno++%repaint==0)
            progress_bar(imgno);
        proj.read(fn_img, false); //true means apply shifts
        if (!DF.containsLabel(MDL_ANGLE_ROT))
        {
            rot  = proj.rot();
            tilt = proj.tilt();
            psi  = proj.psi();
            xoff = proj.Xoff();
            yoff = proj.Yoff();
            flip = proj.flip();
            weight = proj.weight();
        }
        else
        {
                DF.getValue(MDL_ANGLE_ROT,rot,__iter.objId);
                DF.getValue(MDL_ANGLE_TILT,tilt,__iter.objId);
                DF.getValue(MDL_ANGLE_PSI,psi,__iter.objId);
                DF.getValue(MDL_SHIFT_X,xoff,__iter.objId);
                DF.getValue(MDL_SHIFT_Y,yoff,__iter.objId);
                DF.getValue(MDL_FLIP,flip,__iter.objId);
                DF.getValue(MDL_WEIGHT,weight,__iter.objId);
                proj.setEulerAngles(rot,tilt,psi);
                proj.setShifts(xoff,yoff);
                proj.setFlip(flip);
                proj.setWeight(weight);
        }
        Matrix2D<double> A;
        Matrix1D<double> projected_point(3);
        Matrix1D<double> projected_point_2D(2);
        std::vector <Matrix1D<double> > proj_vectors;
        for (int i = 0; i < 12; i++)
        {
            Uproject_to_plane(vertices_vectors[i],rot,tilt,psi,projected_point);
            XX(projected_point_2D)=XX(projected_point);
            YY(projected_point_2D)=YY(projected_point);
            proj_vectors.push_back(projected_point_2D);
        }
        bool scissor;//if true extract this vertex
        for (int i = 0; i < 12; i++)
        {
            scissor=1;
            //Remove vertex if they are close to another vertex
            if (removeCloseVertex)
            {
                for (int j = 0; j < 12; j++)
                {
                    if(i==j)
                        continue;
                    if((proj_vectors[i]-proj_vectors[j]).module()*virusRadius < dim)
                    {
                        scissor=0;
                        break;
                    }
                }
            }
            if (!scissor)
                continue;
            int radius = proj_vectors[i].module()*virusRadius;
            if (radius > minVirusRadius)
            {
                proj_aux=proj;
                FileName fn_tmp;
                proj_aux.setShifts(xoff + XX(proj_vectors[i])*virusRadius,
                                   yoff + YY(proj_vectors[i])*virusRadius);
                proj_aux.getTransformationMatrix(A,true);
                if (!A.isIdentity())
                	selfApplyGeometry(BSPLINE3,proj_aux(),A, IS_INV, DONT_WRAP);
                fn_tmp = fn_img.withoutExtension() + "_";
                fn_tmp.compose(fn_tmp, i, "xmp");
                proj_aux.setShifts(0,0);//shift already applied, same for flip? I do not think so
                proj_aux().selfWindow((Ydim-dim)/2,(Xdim-dim)/2,(Ydim+dim)/2,(Xdim+dim)/2,0);
                Matrix2D<double> Identity(3,3);
                Identity.initIdentity();
                int irandom;
                irandom=rnd_unif(0, 4);
                Matrix2D<double> euler(3, 3), temp;
                Euler_angles2matrix(rot, tilt, psi, euler);
                temp = euler *
                       R_repository[symmetryMatrixVertex(i,irandom)].inv();
                //temp = euler;
                Euler_matrix2angles(temp, rotp, tiltp, psip);
                proj_aux.setEulerAngles(rotp,tiltp,psip);
                proj_aux.write(fn_tmp);

                size_t objId = DFout.addObject();
                DFout.setValue(MDL_IMAGE,fn_tmp,objId);
                DFout.setValue(MDL_ANGLE_ROT,rotp,objId);
                DFout.setValue(MDL_ANGLE_TILT,tiltp,objId);
                DFout.setValue(MDL_ANGLE_PSI,psip,objId);
                DFout.setValue(MDL_SHIFT_X,0,objId);
                DFout.setValue(MDL_SHIFT_Y,0,objId);
                DFout.setValue(MDL_FLIP,flip,objId);
                DFout.setValue(MDL_WEIGHT,weight,objId);
            }
        }
    }
    progress_bar(DF.size());
    FileName fn_tmp1;
    fn_tmp1= fn_doc.withoutExtension()+"_out.doc";
    DFout.write(fn_tmp1);
}

void VirusVertex::assignSymmetryMatricesToVertex()
{
    /** vector with symmetry matrices */
    Matrix2D<double>  R(4, 4),L(4,4);
    Matrix2D<double>  Identity(3,3);
    Identity.initIdentity();
    symmetryMatrixVertex.resize(12, 5);
    R_repository.push_back(Identity);
    for (int isym = 0; isym < SL.SymsNo(); isym++)
    {
        SL.get_matrices(isym, L, R);
        R.resize(3, 3);
        R_repository.push_back(R);
//std::cerr<< R << std::endl;
    }

    //#define CREATEICOSAHEDRALPHANTOM
#ifdef CREATEICOSAHEDRALPHANTOM
    std::ofstream filestr;
    double alpha, beta, gamma;
    filestr.open ("ico.feat");
    filestr
    << "# Phantom description file, (generated with phantom help)\n"
    << "# General Volume Parameters:\n"
    << "#      Xdim      Ydim      Zdim   Background_Density Scale\n"
    << "      2    2    2    0  128\n"
    << "# Feature Parameters: \n";

    for (int i = 0; i < 12; i++)
    {
        Euler_direction2angles(vertices_vectors[i], alpha, beta, gamma);
        filestr
        << "cyl + 1 "
        << XX(vertices_vectors[i])*0.8 << " "
        << YY(vertices_vectors[i])*0.8 << " "
        << ZZ(vertices_vectors[i])*0.8 << " "
        << " .1 .1 .2 "
        << alpha << " "
        << beta  << " "
        << gamma << " "
        << std::endl
        ;
        // cyl  <+/=> <den>    <x0>     <y0>     <z0>    <xradius> <yradius> <height>               <rot> <tilt> <psi>
        // cyl + 1 15   0  0  5  5 15    0 90  0    ; Cylinder  in X
    }
    Matrix1D<double>  myvector(3);
    for (int j = 0; j < R_repository.size(); j++)
    {
        myvector = vectorR3(0., 0.15, 0.9).transpose() * R_repository[j];
        filestr
        << "sph + 1 "
        << XX(myvector)*0.9 << " "
        << YY(myvector)*0.9 << " "
        << ZZ(myvector)*0.9 << " "
        << " .05 "
        << std::endl;
    }


    filestr.close();

#endif

    Matrix1D<double> r(3);
    int k,sixty;
    sixty=0;
    for (int i = 0; i < 12; i++)
    {
        k=0;
        for (int j = 0; j < R_repository.size(); j++)
        {
            if( (vertices_vectors[i]-vectorR3(0., 0., 1.).transpose() * R_repository[j]).module() < 0.0001)
            {
                symmetryMatrixVertex(i,k)=j;
                k++;
                sixty++;
            }
        }
    }
    if (sixty>60)
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, (std::string)"assignSymmetryMatricesToVertex: more than 60 symmetries " );
}
/* Main program ------------------------------------------------------------ */
void VirusVertex::run()
{
    double accuracy=1e-6;
    randomize_random_generator();
    show();
    //read doc and sel files file
    DF.read(fn_doc);

    //load icosahedron vertex
    symmetry=SL.read_sym_file(fn_sym, accuracy);
    loadIcosahedronVertex();
    assignSymmetryMatricesToVertex();

    //Assign symmetry matrices to vertex
    //process one set of angles
    processAngles();
}


