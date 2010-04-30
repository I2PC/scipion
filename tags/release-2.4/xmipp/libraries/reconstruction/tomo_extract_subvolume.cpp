/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2010)
 *
 * Unidad de Bioinformatica del Centro Nacional de Biotecnologia , CSIC
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
#include "tomo_extract_subvolume.h"

//#define DEBUG

// Read arguments ==========================================================
void Prog_tomo_extract_subvolume_prm::read(int argc, char **argv)
{
    // Read command line
    if (checkParameter(argc, argv, "-more_options"))
    {
	usage();
	extendedUsage();
    }
    fn_doc = getParameter(argc, argv, "-doc", "");
    fn_sym = getParameter(argc, argv, "-sym", "c1");
    fn_root = getParameter(argc, argv, "-o", "subvolume");
    size = textToInteger(getParameter(argc, argv, "-size"));
    if (checkParameter(argc, argv, "-center"))
    {
        int i = paremeterPosition(argc, argv, "-center");
        if (i + 3 >= argc)
            REPORT_ERROR(1, "Not enough parameters after -center");
        center_ref.resize(3);
        XX(center_ref) = textToFloat(argv[i+1]);
        YY(center_ref) = textToFloat(argv[i+2]);
        ZZ(center_ref) = textToFloat(argv[i+3]);
    }
    mindist = textToFloat(getParameter(argc, argv, "-mindist","1"));
    verb = textToInteger(getParameter(argc, argv, "-verb", "1"));


}

// Show ====================================================================
void Prog_tomo_extract_subvolume_prm::show()
{

    if (verb > 0)
    {
        // To screen
        std::cerr << " -----------------------------------------------------------------" << std::endl;
        std::cerr << " | Read more about this program in the following publication:    |" << std::endl;
        std::cerr << " |  Scheres ea. in preparation                                   |" << std::endl;
        std::cerr << " |                                                               |" << std::endl;
        std::cerr << " |   *** Please cite it if this program is of use to you! ***    |" << std::endl;
        std::cerr << " -----------------------------------------------------------------" << std::endl;
        std::cerr << "--> Maximum-likelihood multi-reference refinement " << std::endl;
        std::cerr << "  Input docfile            : " << fn_doc << " (" << nr_exp_images << ")" << std::endl;
        std::cerr << "  Symmetry group           : " << fn_sym << std::endl;
        std::cerr << "  Reference center         : " << XX(center_ref)<<" "<<YY(center_ref)<<" "<<ZZ(center_ref)<<std::endl;
        std::cerr << "  Size of subvolumes       : " << size << " pixels"<<std::endl;
        std::cerr << "  Minimum distance         : " << mindist << " pixels" <<std::endl;
        std::cerr << " -----------------------------------------------------------------" << std::endl;
        if (SL.SymsNo()>0)
        {
            std::cerr << "  Number of symmetry-related subvolumes= "<<centers_subvolumes.size()<<std::endl;
            std::cerr << "    With the following centers= "<<std::endl;
            for (int i = 0; i < centers_subvolumes.size(); i++)
            {
                std::cerr << "    * subvolume nr "<<i+1<<" (x,y,z)= "<<XX(centers_subvolumes[i])<<" "<<YY(centers_subvolumes[i])<<" "<<ZZ(centers_subvolumes[i])<<std::endl;
            }
        }
    }

}

// Usage ===================================================================
void Prog_tomo_extract_subvolume_prm::usage()
{
    std::cerr << "Usage:  tomo_extract_subvolume [options] " << std::endl;
    std::cerr << "   -doc <docfile>              : Docfile with input images \n";
    std::cerr << "   -center <x> <y> <z>         : X,Y,Z-coords of subvolume center in reference \n";
    std::cerr << "   -size <int>                 : Size of output subvolumes \n";
    std::cerr << " [ -sym <symgroup=c1> ]        : Symmetry group \n";
    std::cerr << " [ -mindist <1.> ]             : Minimum distance between subvolumes \n";
    //std::cerr << " [ -more_options ]             : Show all possible input parameters \n";
}

// Extended usage ===================================================================
void Prog_tomo_extract_subvolume_prm::extendedUsage()
{
    std::cerr << "Additional options: " << std::endl;
    std::cerr << std::endl;
    exit(1);
}

// Set up a lot of general stuff
// This side info is general, i.e. in parallel mode it is the same for
// all processors! (in contrast to produce_Side_info2)
void Prog_tomo_extract_subvolume_prm::produceSideInfo()
{


#ifdef  DEBUG
    std::cerr<<"Start produceSideInfo"<<std::endl;
#endif

    // Read in Docfile
    DF.read(fn_doc);
    nr_exp_images = DF.dataLineNo();

    // Setup symmetry
    if (!SL.isSymmetryGroup(fn_sym, symmetry, sym_order))
        REPORT_ERROR(3005, (std::string)"tomo_extract_subvolume: Invalid symmetry" +  fn_sym);
    SL.read_sym_file(fn_sym);

    // Check which symmetry operators give unique output points
    centers_subvolumes.clear();
    rotations_subvolumes.clear();
    Matrix2D<double>  L(4, 4), R(4, 4);
    Matrix2D<double>  I(3,3);
    Matrix1D<double>  newcenter(3), distcenter(3);
    double dist;
    bool is_uniq;
    I.initIdentity();

    centers_subvolumes.push_back(center_ref);
    rotations_subvolumes.push_back(I);

    for (int isym = 0; isym < SL.SymsNo(); isym++)
    {
        SL.get_matrices(isym, L, R);
        L.resize(3,3);
        R.resize(3,3);
        newcenter = L * (center_ref.transpose() * R).transpose();
        is_uniq=true;
        for (int i = 0; i < centers_subvolumes.size(); i++)
        {
            distcenter = centers_subvolumes[i] - newcenter;
            dist = sqrt(distcenter.sum2());
            if (dist < mindist)
            {
                is_uniq=false;
                break;
            }
        }
        if (is_uniq)
        {
            centers_subvolumes.push_back(newcenter);               
            rotations_subvolumes.push_back(R);
        }
    }

}

void Prog_tomo_extract_subvolume_prm::processImages(int imgno_start, int imgno_end)
{

    FileName fn_img, fn_out;
    VolumeXmipp vol, volout;
    DocLine DL, DLout;
    Matrix1D<double> center(3), intcenter(3), doccenter(3);
    Matrix2D<double> A(3,3), R(3,3), I(3,3);
    double rot, tilt, psi, rotp, tiltp, psip, x0, xF;
    I.initIdentity();
    DFout.clear();
    SPEED_UP_temps;

    for (int imgno=imgno_start; imgno <imgno_end; imgno++)
    {
        DF.locate(imgno+1);
        DL = DF.get_current_line();
        DF.previous();
        if (DF.get_current_line().Is_comment())
        {
            fn_img = ((DF.get_current_line()).get_text()).erase(0, 3);
        }
        else
        {
            REPORT_ERROR(1,"BUG: no comment in DF where expected....");
        }

        rot = DL[0];
        tilt = DL[1];
        psi = DL[2];
        XX(doccenter) = DL[3];
        YY(doccenter) = DL[4];
        ZZ(doccenter) = DL[5];

        // Read volume
        vol.read(fn_img);
        vol().setXmippOrigin();

        x0 = FIRST_XMIPP_INDEX(size);
        xF = LAST_XMIPP_INDEX(size);
        
        // Tomo_Extract each of the unique subvolumes
        for (int i = 0; i < centers_subvolumes.size(); i++)
        {
            center=centers_subvolumes[i];

            // 1. rotate center
            Euler_angles2matrix(-psi,-tilt,-rot,A);
            M3x3_BY_V3x1(center, A, center);
            // 2. translate center
            center -= doccenter;
            // 3. Apply possible non-integer center to volume
            vol().translate(-center, volout(), DONT_WRAP);
            //4. Window operation and write subvolume to disc
            volout().window(x0,x0,x0,xF,xF,xF);
            fn_out=fn_img.without_extension();
            fn_out+="_sub";
            fn_out.compose(fn_out,i+1,"vol");
            volout.write(fn_out);

            // 5. Calculate output angles: apply symmetry rotation to rot,tilt and psi
            DLout=DL;
            Euler_apply_transf(rotations_subvolumes[i], I, rot, tilt, psi, rotp, tiltp, psip);
            DLout[0] = rotp;
            DLout[1] = tiltp;
            DLout[2] = psip;

            // 6. Output translations will be zero because subvolumes are centered by definition
            DLout[3] = 0.;
            DLout[4] = 0.;
            DLout[5] = 0.;
            DFout.append_comment(fn_out);
            DFout.append_line(DLout);
        }
        
    }

}
