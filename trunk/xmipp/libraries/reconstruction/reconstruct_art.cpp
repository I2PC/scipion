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

#include "reconstruct_art.h"
#include "art_crystal.h"
#include "art_xray.h"
#include "denoise.h"
#include "fourier_filter.h"
#include <data/wavelet.h>
#include <sys/time.h>


ProgReconsART::ProgReconsART()
{
    isMpi = false;
    artRecons = NULL;
}
ProgReconsART::~ProgReconsART()
{
    delete artRecons;
}

void ProgReconsART::setIO(const FileName &fn_in, const FileName &fn_out)
{}

void ProgReconsART::defineParams()
{
    addUsageLine("Generate 3D reconstructions from projections using the ART algorithm (even with SIRT).");
    addUsageLine("+A history file with the information of the reconstruction process is created. You can ");
    addUsageLine("+also give symmetry elements to specify different points of view reusing the projections ");
    addUsageLine("+you have. The reconstruction is performed using some basis (blobs or voxels) in ");
    addUsageLine("+different grids (BCC (by default), FCC or CC).");

    ARTReconsBase::defineParams(this, isMpi);

#ifndef RELEASE_MODE

    addParamsLine(" == Special Parameters for X-rays == ");
    XrayARTRecons::defineParams(this);
#endif

    //    addParamsLine(" == Special Parameters for crystals == ");
    //    CrystalARTRecons::defineParams(this);

    addExampleLine("+++* Basis volume definition: The basis volume is defined internally such that it should cover ",false);
    addExampleLine("+++the space cube corresponding to =(-Xdim/2,-Xdim/2,-Xdim/2)= to =(X dim /2,Xdim/2,Xdim/2)= where Xdim  ",false);
    addExampleLine("+++is the X dimension of the projections.", false);
    addExampleLine("+++* Voxel volume definition: The final reconstructed volume is defined in such a way that it ", false);
    addExampleLine("+++covers the whole basis volume. For these two definitions and the fact that the basis has got ", false);
    addExampleLine("+++some radius (normally 2) the final reconstructions are a little wider than the original ",false);
    addExampleLine("+++projections (usually =2*radius+1=)",false);
    addExampleLine("+++* Crystal reconstructions: Crystal projections must be the deformed projections of the ",false);
    addExampleLine("+++deformed volume. I mean, the lattice vectors of the volume to reconstruct need not to lie ",false);
    addExampleLine("+++on a grid point inside the BCC, FCC or CC grid. If we ant to force that they lie on a grid ",false);
    addExampleLine("+++point we have to deform the volume. When he have this deformed volume and we project it, the ",false);
    addExampleLine("+++two 3D lattice vectors project to two different 2D lattice vectors, that in principle there ",false);
    addExampleLine("+++is no relationship between them. We must force this two 2D lattice vectors to align with ",false);
    addExampleLine("+++=a=(Xdim,0)= and =b=(0,Ydim)= by deforming the projection. This is the second deformation. ",false);
    addExampleLine("+++This complex projections are either generated from phantoms using the program Project or ",false);
    addExampleLine("+++from lists of Fourier spots coming out from MRC as APH files by applying the program ",false);
    addExampleLine("+++Spots2RealSpace2D.",false);
    addExampleLine("+++* Surface constraints: These constraints force the volume to be 0 at known places. ",false);
    addExampleLine("+++This has got the beneficial consequence that mass tends to be compacted where it should be.",false);
    addExampleLine("+++* Grids: This program defines three grids for working with data. The typical simple cubic grid, ",false);
    addExampleLine("+++the face-centered grid (fcc) and the body-centered cubic grid (bcc). We define these grids as follows: ",false);
    addExampleLine("+++|  Simple Cubic Grid  |  %ATTACHURL%/image004.gif  |",false);
    addExampleLine("+++|  Face-centered cubic grid  |  %ATTACHURL%/image008.gif  |",false);
    addExampleLine("+++|  Body-centered cubic grid  |  %ATTACHURL%/image006.gif  |",false);
    addExampleLine("+++where %ATTACHURL%/image010.gif is the set of integers and %ATTACHURL%/image002.gif is a positive ",false);
    addExampleLine("+++real number called sampling distance. See =Gabor H., Geometry of Digital Spaces, Birkhauser, 1998, Massachussets=",false);
    addExampleLine("+++Spots2RealSpace2D. %BR%",false);
    addExampleLine("+++In the case of the simple cubic grid the voxels are cubes of volume %ATTACHURL%/image013.gif ³",false);
    addExampleLine("+++. The voxels associated to the bcc grid are truncated octahedra, polyhedra of 8 hexagonal ",false);
    addExampleLine("+++faces and 6 square faces, of volume 4 %ATTACHURL%/image013.gif ³ .Finally, the voxels associated ",false);
    addExampleLine("+++to the fcc grid are rhombic dodacehedra, polyhedra of 12 identical rhombic faces, of volume 2 ",false);
    addExampleLine("+++%ATTACHURL%/image013.gif ³. However, in practice the above definitions are implemented using a ",false);
    addExampleLine("+++combination of simple cubic grids. Clearly, a simple cubic grid defines the complex simple cubic ",false);
    addExampleLine("+++grid and thus the relative size used in the implementation is exactly the %ATTACHURL%/image013.gif ³. ",false);
    addExampleLine("+++For defining the bcc grid two simple cubic grids are used, it can be seen from the definition for ",false);
    addExampleLine("+++the bcc above that the valid positions for even values of =z=, =x= and =y= have to be also even, ",false);
    addExampleLine("+++in the case of odd values of =z= , =x= and =y= have to be odd. Consequently, the relative size used ",false);
    addExampleLine("+++in the BCC grid is equivalent to 2 %ATTACHURL%/image013.gif above. For defining the fcc grid four ",false);
    addExampleLine("+++simple cubic grids are used, it can be seen from the definition for the fcc above that the valid ",false);
    addExampleLine("+++positions for even values of =z=, the sum of =x= and =y= has to be even; in the case of odd values ",false);
    addExampleLine("+++of =z= the sum of =x= and =y= has to be odd. Consequently, the relative size used in the FCC grid ",false);
    addExampleLine("+++is equivalent to 2 %ATTACHURL%/image013.gif above. %BR%",false);
    addExampleLine(" ",false); // For style twiki reasons
    addExampleLine("Basic reconstructing commands:",false);
    addExampleLine("reconstruct_art -i projections.sel -o artrec");
    addExampleLine("Create the equivalent noisy reconstruction:",false);
    addExampleLine("reconstruct_art -i projections.sel -o artrec --noisy_reconstruction");
    //    addExampleLine("Reconstruct using SIRT parallelization algorithm for five iterations:",false);
    //    addExampleLine("reconstruct_art -i projections.sel -o artrec -n 5 --parallel_mode SIRT");
    addExampleLine("Save the basis information at each iteration:",false);
    addExampleLine("reconstruct_art -i projections.sel -o artrec -n 3 --save_basis");
}

void ProgReconsART::readParams()
{
    //    if (checkParam("--crystal"))
    //        artRecons = new CrystalARTRecons;
    //    else

#ifndef RELEASE_MODE
    if (checkParam("--xray"))
        artRecons = new XrayARTRecons;
    else
#endif

        artRecons = new SinPartARTRecons;

    artRecons->readParams(this);

    if (artRecons->artPrm.threads > 1 && isMpi)
        REPORT_ERROR(ERR_ARG_BADCMDLINE, "Threads not compatible in mpi version.");
    if (!isMpi && artRecons->artPrm.parallel_mode != BasicARTParameters::ART)
        REPORT_ERROR(ERR_ARG_BADCMDLINE, "If --parallel_mode is passed, then mpi version must be used.");

}

void ProgReconsART::show()
{
    if (verbose > 0)
    {
        std::cout << " =====================================================================" << std::endl;
        std::cout << " ART reconstruction method " << std::endl;
        std::cout << " =====================================================================" << std::endl;
        std::cout << " Projections file             : "  << artRecons->artPrm.fn_sel << std::endl;
        std::cout << " Output rootname              : "  << artRecons->artPrm.fn_root << std::endl;
        std::cout << " Iterations                   : "  << artRecons->artPrm.no_it << std::endl;
        std::cout << " Lambda                       : "  << artRecons->artPrm.lambda_list(0) << std::endl;
        std::cout << std::endl;
        std::cout << " More info in the history file: "  << artRecons->artPrm.fn_root + ".hist" << std::endl;
        std::cout << " ---------------------------------------------------------------------\n" << std::endl;
    }
}

void ProgReconsART::run()
{
    BasicARTParameters &artPrm = artRecons->artPrm;

    Image<double> vol_voxels;
    GridVolume vol_basis;
    // Configure time clock
    time_config();

    struct timeval start_time, end_time;
    long int init_usecs, process_usecs, finish_usecs;

    gettimeofday(&start_time, NULL);

    show();
    // Produce side information and initial volume
    artRecons->preIterations(vol_basis);
    if (verbose > 0)
    {
        std::cout << " ---------------------------------------------------------------------" << std::endl;
        std::cout << " Projections                  : "  << artRecons->artPrm.numIMG << std::endl;
        std::cout << " ---------------------------------------------------------------------" << std::endl;
    }
    // Show parameters and initiate history
    artRecons->initHistory(vol_basis);

    gettimeofday(&end_time, NULL);

    init_usecs = (end_time.tv_sec-start_time.tv_sec)*1000000+(end_time.tv_usec-start_time.tv_usec);

    gettimeofday(&start_time,NULL);

    // Iterations
    artRecons->iterations(vol_basis);

    gettimeofday(&end_time,NULL);

    process_usecs = (end_time.tv_sec-start_time.tv_sec)*1000000+(end_time.tv_usec-start_time.tv_usec);

    gettimeofday(&start_time,NULL);

    // Finish iterations
    artRecons->postIterations(vol_basis);

    // Write final volume
    int Xoutput_volume_size=(artPrm.Xoutput_volume_size==0) ?
                            artPrm.projXdim:artPrm.Xoutput_volume_size;
    int Youtput_volume_size=(artPrm.Youtput_volume_size==0) ?
                            artPrm.projYdim:artPrm.Youtput_volume_size;
    int Zoutput_volume_size=(artPrm.Zoutput_volume_size==0) ?
                            artPrm.projXdim:artPrm.Zoutput_volume_size;

    //   int min_distance = ceil( ( 2 * artPrm.grid_relative_size ) / artPrm.basis.blob.radius) + 1;

    artPrm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
                                Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size,artPrm.threads);

    vol_voxels.write(artPrm.fn_root+".vol");

    if (artPrm.tell&TELL_SAVE_BASIS)
        vol_basis.write(artPrm.fn_root+".basis");
    artPrm.fh_hist->close();

    gettimeofday(&end_time,NULL);

    finish_usecs = (end_time.tv_sec-start_time.tv_sec)*1000000+(end_time.tv_usec-start_time.tv_usec);

    std::cout << "INIT_TIME: " << (double)init_usecs/(double)1000000 << std::endl;
    std::cout << "PROCESS_TIME: " << (double)process_usecs/(double)1000000 << std::endl;
    std::cout << "FINISH_TIME: " << (double)finish_usecs/(double)1000000 << std::endl;

}

