/***************************************************************************
 *
 * Authors:    Sjors Scheres           scheres@cnb.csic.es (2004)
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

#include <data/metadata_extension.h>
#include "break_symmetry.h"
#include "reconstruct_wbp.h"

// Read ===================================================================
void Prog_Break_Sym_prm::read(int argc, char **argv)
{
    //Read Break_Sym parameters
    fn_sel = getParameter(argc, argv, "-i");
    fn_root = getParameter(argc, argv, "-o", "breaksym");
    fn_vol = getParameter(argc, argv, "-vol");
    fn_mask = getParameter(argc, argv, "-mask", "");
    // Fill volume selfile
    if (!fn_vol.isMetaData())
    {
        SFvol.addObject();
        SFvol.setValue(MDL_IMAGE,fn_vol);
    }
    else
        SFvol.read(fn_vol);
    Nvols = SFvol.size();

    fn_sym = getParameter(argc, argv, "-sym");
    mask_radius = textToFloat(getParameter(argc, argv, "-mask_radius", "-1"));
    eps = textToFloat(getParameter(argc, argv, "-eps", "5e-5"));
    verb = textToInteger(getParameter(argc, argv, "-verb", "1"));
    Niter = textToInteger(getParameter(argc, argv, "-iter", "100"));
    // Hidden
    istart = textToInteger(getParameter(argc, argv, "-istart", "1"));

    // Read stuff into memory
    SF.read(fn_sel);
    ImgSize(SF, dim, dim);
    SL.read_sym_file(fn_sym);
    Image<double> vol;
    Nvols = 0;
    FOR_ALL_OBJECTS_IN_METADATA(SFvol)
    {
        FileName fn_img;
        SFvol.getValue(MDL_IMAGE,fn_img);
        vol.read(fn_img);
        vol().setXmippOrigin();
        vols.push_back(vol());
        Nvols++;
    }
    if (fn_mask != "")
    {
        vol.read(fn_mask);
        mask().resize(vol());
        symmetrize(SL, vol(), mask());
    }
}

// Usage ===================================================================
void Prog_Break_Sym_prm::usage()
{
    std::cerr << "Usage:  Break_Sym [options] " << std::endl;
    std::cerr << "   -i <selfile>                : Selfile with input images \n"
    << "   -vol <volume/selfile>       : Initial reference volume \n"
    << "                               :  OR selfile with multiple reference volumes\n"
    << " [ -sym <symfile> ]            : Symmetry to be broken \n"
    << " [ -mask <volume> ]            : Limit correlation to masked region \n"
    << " [ -o <rootname=\"breaksym\"> ]  : Output rootname \n"
    << " [ -iter <int=100> ]           : Maximum number of iterations \n"
    << " [ -mask_radius <rad.> ]       : Radius for masking volume (default is no masking)\n";

}
// Show ======================================================================
void Prog_Break_Sym_prm::show()
{
    if (verb > 0)
    {
        // To screen
        std::cerr << " =================================================================" << std::endl;
        std::cerr << " Break symmetry " << std::endl;
        std::cerr << " =================================================================" << std::endl;
        if (Nvols == 1)
            std::cerr << " Initial reference volume : " << fn_vol << std::endl;
        else
        {
            std::cerr << " Selfile with references  : " << fn_vol << std::endl;
            std::cerr << "   with # of volumes      : " << Nvols << std::endl;
        }
        std::cerr << " Experimental images:     : " << SF.getFilename() << " (" << SF.size() << ")" << std::endl;
        std::cerr << " Symmetry file:           : " << fn_sym << std::endl;
        std::cerr << " Output rootname          : " << fn_root << std::endl;
        std::cerr << " Convergence criterion    : " << eps << std::endl;
        if (mask_radius > 0)
            std::cerr << " Mask radius              : " << mask_radius << std::endl;
        std::cerr << " =================================================================" << std::endl;

        // Also open and fill history file
        fh_hist.open((fn_root + ".hist").c_str(), std::ios::out);
        if (!fh_hist)
            REPORT_ERROR(3008, (std::string)"Prog_Break_Sym: Cannot open file " + fn_root + ".hist");

        fh_hist << " =================================================================" << std::endl;
        fh_hist << " Break symmetry " << std::endl;
        fh_hist << " =================================================================" << std::endl;
        if (Nvols == 1)
            fh_hist << " Initial reference volume : " << fn_vol << std::endl;
        else
        {
            fh_hist << " Selfile with references  : " << fn_vol << std::endl;
            fh_hist << "   with # of volumes      : " << Nvols << std::endl;
        }
        fh_hist << " Experimental images:     : " << SF.getFilename() << " (" << SF.size() << ")" << std::endl;
        fh_hist << " Symmetry file:           : " << fn_sym << std::endl;
        fh_hist << " Output rootname          : " << fn_root << std::endl;
        fh_hist << " Convergence criterion    : " << eps << std::endl;
        if (mask_radius > 0)
            fh_hist << " Mask radius              : " << mask_radius << std::endl;
        fh_hist << " =================================================================" << std::endl;
    }
}

// Projection of the reference (blob) volume =================================
void Prog_Break_Sym_prm::process_one_image(Image<double> &img, int &opt_vol,
        int &opt_sym, double &maxcorr)
{
    double rot, tilt, psi, newrot, newtilt, newpsi;
    double Xsq, Asq, corr;
    double opt_rot, opt_tilt, opt_psi;
    Projection proj, projmask;
    Matrix2D<double> L(4, 4), R(4, 4);

    rot = img.rot();
    tilt = img.tilt();
    psi = img.psi();
    maxcorr = -999.;

    Image<double> Itmp;
    FileName fnt;

    if (fn_mask != "")
    {
        projectVolume(mask(), projmask, dim, dim, rot, tilt, psi);
        projmask().setXmippOrigin();
    }
    else
    {
        projmask().resize(dim, dim);
        projmask().setXmippOrigin();
        projmask().initConstant(1.);
    }

    // Pre-calculate Xsq:
    Xsq = 0.;
    FOR_ALL_ELEMENTS_IN_ARRAY2D(img())
    {
        if (i*i + j*j > dim*dim / 4)
            A2D_ELEM(projmask(), i, j) = 0.;
        if (A2D_ELEM(projmask(), i, j) > 0.)
            Xsq += A2D_ELEM(img(), i, j) * A2D_ELEM(img(), i, j);
    }
    if (verb > 1)
    {
        Itmp() = projmask();
        Itmp.write("mask.xmp");
    }

    for (int volno = 0; volno < Nvols; volno++)
    {
        projectVolume(vols[volno], proj, dim, dim, rot, tilt, psi);
        proj().setXmippOrigin();
        corr = Asq = 0.;
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(proj())
        {
            if (dAij(projmask(), i, j) > 0.)
            {
                corr += dAij(proj(), i, j) * dAij(img(), i, j);
                Asq += dAij(proj(), i, j) * dAij(proj(), i, j);
            }
        }
        corr /= sqrt(Asq) * sqrt(Xsq);
        if (verb > 1)
            std::cout << " " << corr;
        if (corr > maxcorr)
        {
            maxcorr = corr;
            opt_rot = rot;
            opt_tilt = tilt;
            opt_psi = psi;
            opt_sym = 0;
            opt_vol = volno;
        }
        if (verb > 1)
        {
            Itmp() = proj();
            Itmp.write("proj.xmp");
            Itmp() = img();
            Itmp.write("img.xmp");
        }
        for (int isym = 0; isym < SL.SymsNo(); isym++)
        {
            SL.get_matrices(isym, L, R);
            L.resize(3, 3);
            R.resize(3, 3);
            Euler_apply_transf(L, R, rot, tilt, psi, newrot, newtilt, newpsi);
            projectVolume(vols[volno], proj, dim, dim, newrot, newtilt, newpsi);
            proj().setXmippOrigin();
            if (verb > 1)
            {
                Itmp() = proj();
                fnt.compose("sym", isym, "xmp");
                Itmp.write(fnt);
            }

            corr = Asq = 0.;
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(proj())
            {
                if (dAij(projmask(), i, j) > 0.)
                {
                    corr += dAij(proj(), i, j) * dAij(img(), i, j);
                    Asq += dAij(proj(), i, j) * dAij(proj(), i, j);
                }
            }
            corr /= sqrt(Asq) * sqrt(Xsq);
            if (verb > 1)
                std::cout << " " << corr;
            if (corr > maxcorr)
            {
                maxcorr = corr;
                opt_rot = newrot;
                opt_tilt = newtilt;
                opt_psi = newpsi;
                opt_sym = isym + 1;
                opt_vol = volno;
            }
        }
    }

    // Read image again because of applied shifts
    img.read(img.name());
    img.setEulerAngles(opt_rot, opt_tilt, opt_psi);
    img.write(img.name());
    if (verb > 1)
    {
        std::cout << std::endl;
        std::cout << opt_rot << " " << opt_tilt << " " << opt_psi << std::endl;
    }
}

// Projection of the reference (blob) volume =================================
void Prog_Break_Sym_prm::process_selfile(MetaData &SF, std::vector<MetaData> &SFout,
    double &avecorr)
{
    Image<double> img;
    MetaData SFtmp;
    FileName fn_img;
    double maxcorr;
    int c, nn, imgno, opt_vol, opt_sym;

    // Initialize
    nn = SF.size();
    if (verb > 0)
        init_progress_bar(nn);
    c = XMIPP_MAX(1, nn / 60);

    SFout.clear();
    for (int volno = 0; volno < Nvols; volno++)
        SFout.push_back(SFtmp);

    avecorr = 0.;
    imgno = 0;
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        SF.getValue(MDL_IMAGE,fn_img);

        // read image, only applying the shifts
        img.read(fn_img, false, false, false, true);
        img().setXmippOrigin();

        process_one_image(img, opt_vol, opt_sym, maxcorr);
        SFout[opt_vol].addObject();
        SFout[opt_vol].setValue(MDL_IMAGE,fn_img);
        SFout[opt_vol].setValue(MDL_ENABLED,true);

        avecorr += maxcorr;
        if (verb > 0)
            if (imgno % c == 0)
                progress_bar(imgno);
        imgno++;
    }
    if (verb > 0)
        progress_bar(nn);
    avecorr /= (double)imgno;
}


// Reconstruction using the ML-weights ==========================================
void Prog_Break_Sym_prm::reconstruction(int argc, char **argv, MetaData &SF,
                                        int iter, int volno)
{
    Image<double>          new_vol;
    FileName               fn_tmp, fn_insel;

    // Setup selfile for reconstruction
    fn_tmp = fn_root + "_it";
    fn_tmp.compose(fn_tmp, iter, "");
    if (Nvols > 1)
    {
        fn_tmp += "_vol";
        fn_tmp.compose(fn_tmp, volno + 1, "");
    }
    fn_insel = fn_tmp + ".sel";
    SF.write(fn_insel);

    Prog_WBP_prm wbp_prm;
    if (verb > 0)
        std::cerr << "--> WBP reconstruction " << std::endl;

    // read command line (fn_sym, angular etc.)
    wbp_prm.read(argc, argv);
    if (mask_radius > 0)
        wbp_prm.diameter = (int)mask_radius * 2;

    wbp_prm.fn_sel = fn_insel;
    wbp_prm.fn_sym = "";
    wbp_prm.show();
    wbp_prm.fn_out = fn_tmp + ".vol";
    wbp_prm.produce_Side_info();
    wbp_prm.apply_2Dfilter_arbitrary_geometry(wbp_prm.SF, new_vol());
    new_vol.write(wbp_prm.fn_out);
    vols[volno] = new_vol();

    if (verb > 0)
        std::cerr << " =================================================================" << std::endl;
}

