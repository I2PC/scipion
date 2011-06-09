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
#include "denoise.h"
#include "fourier_filter.h"
#include <data/wavelet.h>
#include <sys/time.h>


ProgReconsART::ProgReconsART()
{}
ProgReconsART::~ProgReconsART()
{}

void ProgReconsART::setIO(const FileName &fn_in, const FileName &fn_out)
{}

void ProgReconsART::defineParams()
{
    ARTReconsBase::defineParams(this);
    addParamsLine(" == Special Parameters for crystals == ");
    CrystalARTRecons::defineParams(this);
}

void ProgReconsART::readParams()
{
    if (checkParam("--crystal"))
        artRecons = new CrystalARTRecons;
    else
        artRecons = new ARTReconsBase;

    artRecons->readParams(this);
    //    artPrm.readParams(this);
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

    // Produce side information and initial volume
    artRecons->produceSideInfo(vol_basis);

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
    artRecons->finishIterations(vol_basis);

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

