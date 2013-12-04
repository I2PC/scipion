/***************************************************************************
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Joaquin Oton (joton@cnb.csic.es)
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

#include <algorithm>
#include "base_art_recons.h"
#include "recons_misc.h"
#include "fourier_filter.h"


void ARTReconsBase::readParams(XmippProgram * program)
{
    artPrm.readParams(program);
}

void ARTReconsBase::preProcess(GridVolume &vol_basis0, int level, int rank)
{
    artPrm.produceSideInfo(vol_basis0, level, rank);
}

void ARTReconsBase::iterations(GridVolume &vol_basis, int rank)
{
    // Some variables .......................................................
    int        ART_numIMG;              // Normalizing factor
    GridVolume *ptr_vol_out;            // Pointer to output volume
    GridVolume vol_basis_out;           // Output volume (only useful in SIRT)
    Image<double> vol_voxels;             // This one is useful only in the
    // case of saving intermediate volumes
    int Xoutput_volume_size, Youtput_volume_size, Zoutput_volume_size;
    double aux_tilt;

    // Projection related ...................................................
    Projection      read_proj;          // Projection read from file
    Projection      theo_proj;          // Projection from the
    // reconstruction
    Projection      alig_proj;          // Projection with the correlation
    // maps needed for translation alignment
    Projection      corr_proj;          // Image with the correction
    // factors for unitary basis
    Projection      diff_proj;          // Difference between the
    // theoretical and real image

    //    preIterations(vol_basis);

    // Reconstruction results ...............................................
    double          mean_error;
    double          global_mean_error,global_mean_error_1stblock;

    // Initialize residual image vector for wlsART ..........................
    if (artPrm.WLS)
    {
        artPrm.residual_imgs.clear();
        for (int iact_proj = 0; iact_proj < artPrm.numIMG ; iact_proj++)
        {
            read_proj.read(artPrm.IMG_Inf[iact_proj].fn_proj);
            read_proj().setXmippOrigin();
            read_proj().initZeros();
            artPrm.residual_imgs.push_back(read_proj);
        }
    }

    // Some initialization ..................................................
    Xoutput_volume_size = (artPrm.Xoutput_volume_size==0) ?
                          artPrm.projXdim : artPrm.Xoutput_volume_size;
    Youtput_volume_size = (artPrm.Youtput_volume_size==0) ?
                          artPrm.projYdim : artPrm.Youtput_volume_size;
    Zoutput_volume_size = (artPrm.Zoutput_volume_size==0) ?
                          artPrm.projXdim : artPrm.Zoutput_volume_size;

    // POCS constraints .....................................................
    POCSClass POCS(&artPrm, Zoutput_volume_size, Youtput_volume_size,
                   Xoutput_volume_size);

    // Variability analysis .................................................
    VariabilityClass VC(&artPrm, Zoutput_volume_size, Youtput_volume_size,
                        Xoutput_volume_size);

    // Noisy reconstruction .................................................
    GridVolume vol_basis_noisy; // Output volume (only for ART)
    Projection noisy_projection;
    MetaData SF_noise, SF_signal;
    if (artPrm.noisy_reconstruction)
    {
        vol_basis_noisy = vol_basis;
        vol_basis_noisy.initZeros();
    }

    // If SIRT set normalizing factor and create output volume ..............
    if (artPrm.parallel_mode==BasicARTParameters::SIRT )
    {
        ART_numIMG = artPrm.numIMG;
        vol_basis_out = vol_basis;         // Copy the structure of vol_basis
        ptr_vol_out = &vol_basis_out;      // Pointer to output volume
        if (artPrm.variability_analysis)
            ptr_vol_out->initZeros();
    }
    else if ( artPrm.parallel_mode==BasicARTParameters::pfSIRT ||
              artPrm.parallel_mode==BasicARTParameters::pSIRT ||
              artPrm.parallel_mode==BasicARTParameters::pSART ||
              artPrm.parallel_mode==BasicARTParameters::pCAV )
    {
        // In those cases, normalization must be done at the top level program once we
        // have the reconstruction values. Have a look ar /Applications/Src/MPIArt/mpi_art.cc for
        // an example
        ART_numIMG = 1;
        ptr_vol_out = &vol_basis_out;      // Pointer to output volume
        vol_basis_out = vol_basis;         // Copy the structure of vol_basis
        // and possible initial values
    }
    else if (artPrm.eq_mode == CAV)
    {
        ART_numIMG = 1;                    // Normalizing factor = Total no. images
        ptr_vol_out = &vol_basis_out;      // Pointer to output volume
        vol_basis_out = vol_basis;         // Copy the structure of vol_basis
        // and possible initial values
    }
    else
    {
        ART_numIMG = 1;                    // No normalizing factor
        ptr_vol_out = &vol_basis;          // Output volume is the same as
        // input one
    }
    // Now iterate ..........................................................
    ProcessorTimeStamp time0;                    // For measuring the elapsed time
    annotate_processor_time(&time0);
    int images=0;
    double mean_error_2ndblock,pow_residual_imgs;
    bool iv_launched=false;
    for (int it = 0; it < artPrm.no_it; it++)
    {
        // Initialization of some variables
        global_mean_error = 0;
        global_mean_error_1stblock = 0;
        POCS.newIteration();
        VC.newIteration();
        if (rank==-1)
        {
            std::cout << "Running iteration " << it << " with lambda= " << artPrm.lambda(it)<< "\n" << std::endl;
            if (!(artPrm.tell&TELL_SHOW_ERROR))
                init_progress_bar(artPrm.numIMG);
        }

        // For each projection -----------------------------------------------
        for (int act_proj = 0; act_proj < artPrm.numIMG ; act_proj++)
        {
            POCS.newProjection();

            // Select next projection ........................................
            int iact_proj; // Index inside the sorting information for act_proj
            if (artPrm.tell&TELL_MANUAL_ORDER)
            {
                int proj_number;
                std::cout << "Introduce next projection to study: ";
                std::cin  >> proj_number;
                int sym_number=-1;
                if (artPrm.SL.symsNo()!=0)
                {
                    std::cout << "Introduce symmetry to study: ";
                    std::cin  >> sym_number;
                }
                iact_proj=0;
                while (iact_proj<artPrm.numIMG)
                {
                    if (artPrm.IMG_Inf[iact_proj].fn_proj.getNumber()==proj_number &&
                        artPrm.IMG_Inf[iact_proj].sym==sym_number)
                        break;
                    iact_proj++;
                }
            }
            else
            {
                iact_proj = artPrm.ordered_list(act_proj);
            }

            ReconsInfo &imgInfo = artPrm.IMG_Inf[iact_proj];

            read_proj.read(imgInfo.fn_proj, artPrm.apply_shifts, DATA, &imgInfo.row);
            read_proj().setXmippOrigin();
            read_proj.setEulerAngles(imgInfo.rot, imgInfo.tilt, imgInfo.psi);

            // If noisy reconstruction
            if (artPrm.noisy_reconstruction)
            {
                init_random_generator(imgInfo.seed);
                noisy_projection().resize(read_proj());
                noisy_projection().initRandom(0, 1, RND_GAUSSIAN);
                noisy_projection().setXmippOrigin();
                noisy_projection.setEulerAngles(imgInfo.rot,
                                                imgInfo.tilt,
                                                imgInfo.psi);
                if ( it == 0 && imgInfo.sym==-1 )
                {
                    FileName fn_noise;
                    MDRow row;
                    fn_noise.compose(read_proj.name().getPrefixNumber(),artPrm.fn_root+"_noise_proj.stk");

                    noisy_projection.write(fn_noise);
                    row.setValue(MDL_IMAGE, fn_noise);
                    row.setValue(MDL_ENABLED, 1);
                    row.setValue(MDL_ANGLE_PSI, read_proj.psi());
                    row.setValue(MDL_ANGLE_ROT, read_proj.rot());
                    row.setValue(MDL_ANGLE_TILT, read_proj.tilt());
                    SF_noise.addRow(row);

                    row.setValue(MDL_IMAGE, read_proj.name());
                    SF_signal.addRow(row);
                }
            }

            //skipping if  tilt greater than max_tilt
            //tilt is in between 0 and 360
            aux_tilt=read_proj.tilt();
            if((aux_tilt > artPrm.max_tilt && aux_tilt < 180.-artPrm.max_tilt) ||
               (aux_tilt > artPrm.max_tilt + 180 && aux_tilt < 360.-artPrm.max_tilt))
            {
                std::cout << "Skipping Proj no: " << iact_proj
                << "tilt=" << read_proj.tilt()  << std::endl;
                continue;
            }

            // Projection extension? .........................................
            if (artPrm.proj_ext!=0)
            {
                read_proj().selfWindow(
                    STARTINGY (read_proj())-artPrm.proj_ext,
                    STARTINGX (read_proj())-artPrm.proj_ext,
                    FINISHINGY(read_proj())+artPrm.proj_ext,
                    FINISHINGX(read_proj())+artPrm.proj_ext);
                noisy_projection().resize(read_proj());
            }

            //Skip if desired
            if (artPrm.tell&TELL_ONLY_SYM)
            {
                if( imgInfo.sym != -1)
                {
                    std::cout << "Skipping Proj no: " << iact_proj
                    << " with symmetry no: " << imgInfo.sym
                    <<  std::endl;
                    continue;
                }
                else
                    std::cout << "NO Skipping Proj no: " << iact_proj
                    << " with symmetry no: " << imgInfo.sym
                    <<  std::endl;
            }

            // For wlsART: use alig_proj for residual image!!
            if (artPrm.WLS)
                alig_proj=artPrm.residual_imgs[iact_proj];

            // Is there a mask ...............................................
            MultidimArray<int> mask;
            const MultidimArray<int> *maskPtr=NULL;
            if (artPrm.goldmask<1e6 || artPrm.shiftedTomograms)
            {
                mask.resize(read_proj());
                FOR_ALL_ELEMENTS_IN_ARRAY2D(read_proj())
                {
                    mask(i,j)=1;
                    if ((read_proj(i,j)<artPrm.goldmask && artPrm.goldmask<1e6) ||
                        (ABS(read_proj(i,j))<1e-5 && artPrm.shiftedTomograms))
                        mask(i,j)=0;
                }
                maskPtr=&mask;
            }

            // Apply the reconstruction algorithm ............................
            // Notice that the following function is specific for each art type
            singleStep(vol_basis, ptr_vol_out,
                       theo_proj, read_proj, imgInfo.sym , diff_proj,
                       corr_proj, alig_proj,
                       mean_error, ART_numIMG, artPrm.lambda(it),
                       images, imgInfo.fn_ctf, maskPtr,
                       artPrm.refine);

            if (artPrm.WLS)
            {
                artPrm.residual_imgs[iact_proj]=alig_proj;
                global_mean_error_1stblock += diff_proj().sum2()/(XSIZE(diff_proj())*YSIZE(diff_proj()));
            }

            global_mean_error += mean_error;
            if (artPrm.noisy_reconstruction)
            {
                double noise_mean_error;
                singleStep(vol_basis_noisy, &vol_basis_noisy,
                           theo_proj, noisy_projection, imgInfo.sym,
                           diff_proj,  corr_proj, alig_proj,
                           noise_mean_error, ART_numIMG, artPrm.lambda(it),
                           images, imgInfo.fn_ctf,maskPtr,
                           false);
            }

            // Force symmetry in the volume ..................................
            // so far only crystallographic
            if (artPrm.sym_each &&
                act_proj%artPrm.sym_each==0 &&
                (act_proj!=0 || artPrm.sym_each==1))
            {
                //                apply_symmetry(vol_basis,ptr_vol_out, artPrm.grid_type);
                if (artPrm.noisy_reconstruction)
                    //                    apply_symmetry(vol_basis_noisy,&vol_basis_noisy, artPrm.grid_type);
                    ;
            }

            // Apply POCS ....................................................
            POCS.apply(vol_basis,it,images);
            if (artPrm.noisy_reconstruction)
                POCS.apply(vol_basis_noisy,it,images);

            // Variability analysis ..........................................
            if (artPrm.variability_analysis)
                VC.newUpdateVolume(ptr_vol_out,read_proj);

            // Apply anisotropic diffusion ...................................
            if (it>=1 && artPrm.diffusionWeight>-1)
            {
                artPrm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
                                            Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
                Matrix1D<double> alpha;
                alpha=vectorR3(1.0,1.0,artPrm.diffusionWeight);
                double regError=tomographicDiffusion(vol_voxels(),alpha,
                                                     artPrm.lambda(it)/100);
                if (artPrm.tell&TELL_SHOW_ERROR)
                    std::cout << "Regularization error = " << regError << std::endl;
                *artPrm.fh_hist << "Regularization error = " << regError << std::endl;
                artPrm.basis.changeFromVoxels(vol_voxels(),vol_basis,artPrm.grid_type,
                                              artPrm.grid_relative_size, NULL, NULL, artPrm.R, artPrm.threads);
            }

            // Apply sparsity constraint .....................................
            if (it>=1 && artPrm.sparseEps>0 &&
                artPrm.parallel_mode!=BasicARTParameters::SIRT &&
                artPrm.parallel_mode!=BasicARTParameters::pSIRT &&
                artPrm.parallel_mode!=BasicARTParameters::pfSIRT)
            {
                artPrm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
                                            Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
                forceDWTSparsity(vol_voxels(),artPrm.sparseEps);
                artPrm.basis.changeFromVoxels(vol_voxels(),vol_basis,artPrm.grid_type,
                                              artPrm.grid_relative_size, NULL, NULL, artPrm.R, artPrm.threads);
            }

            // Show results ..................................................
            *artPrm.fh_hist << imgInfo.fn_proj << ", sym="
            << imgInfo.sym << "\t\t" << mean_error;
            if (POCS.apply_POCS)
                *artPrm.fh_hist << "\tPOCS:" << POCS.POCS_mean_error;
            *artPrm.fh_hist << std::endl;
            if (artPrm.tell&TELL_SHOW_ERROR)
            {
                std::cout << imgInfo.fn_proj << ", sym="
                << imgInfo.sym << "\t\t" << mean_error;
                if (POCS.apply_POCS)
                    std::cout        << "\tPOCS:" << POCS.POCS_mean_error;
                std::cout << std::endl;
            }
            else if (act_proj%XMIPP_MAX(1,artPrm.numIMG/60)==0)
                progress_bar(act_proj);

            if (artPrm.tell&TELL_STATS)
            {
                std::cout << "   read      ";
                read_proj().printStats();
                std::cout << "\n   theo      ";
                theo_proj().printStats();
                std::cout << "\n   corr      ";
                corr_proj().printStats();
                std::cout << "\n   alig      ";
                alig_proj().printStats();
                std::cout << "\n   diff      ";
                diff_proj().printStats();
                std::cout << "\n   subvol(0) ";
                (*ptr_vol_out)(0)().printStats();
                std::cout << "\n";
            }

            if (artPrm.tell&TELL_SAVE_AT_EACH_STEP)
            {
                std::cout << "Stats PPPdiff.xmp: ";
                diff_proj().printStats();
                std::cout << std::endl;
                std::cout << "Stats PPPtheo.xmp: ";
                theo_proj().printStats();
                std::cout << std::endl;
                std::cout << "Stats PPPread.xmp: ";
                read_proj().printStats();
                std::cout << std::endl;
                std::cout << "Stats PPPcorr.xmp: ";
                corr_proj().printStats();
                std::cout << std::endl;
                diff_proj.write("PPPdiff.xmp");
                theo_proj.write("PPPtheo.xmp");
                read_proj.write("PPPread.xmp");
                corr_proj.write("PPPcorr.xmp");
                if(act_proj!=0)
                    alig_proj.write("PPPalign.xmp");
                vol_basis.write("PPPbasis.basis");
                artPrm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
                                            Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
                std::cout << "Stats PPPvol.vol : ";
                vol_voxels().printStats();
                std::cout << std::endl;
                vol_voxels.write("PPPvol.vol");

                if (!iv_launched)
                {
                    if (!system("xmipp_show -img PPPdiff.xmp PPPtheo.xmp PPPread.xmp PPPcorr.xmp -dont_apply_geo -poll &"))
                    	REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
                    if (!system("xmipp_show -vol PPPvol.vol -poll &"))
                    	REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
                    iv_launched=true;
                }
                std::cout << "\nHit any key and enter\n";
                char c;
                std::cin >> c;
            }

            // Save intermediate
            if (((artPrm.tell&TELL_SAVE_INTERMIDIATE) | (artPrm.tell&TELL_IV)) &&
                artPrm.save_intermidiate_every!=0 &&
                act_proj%artPrm.save_intermidiate_every==0)
            {
                if (artPrm.tell&TELL_SAVE_INTERMIDIATE)
                    std::cout << "\nSaving intermediate ...\n"
                    << "Converting basis volume to voxels ...\n";
                // Save reconstructed volume
                artPrm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
                                            Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
                if (artPrm.tell&TELL_SAVE_INTERMIDIATE)
                    vol_voxels.write(artPrm.fn_root+"it"+integerToString(it)+"proj"+
                                     integerToString(act_proj,5)+".vol");
                else
                    vol_voxels.write("PPPvol.vol");

                // Launch viewer
                if (!iv_launched)
                {
                    if (!system("xmipp_show -i PPPvol.vol --poll &"))
                    	REPORT_ERROR(ERR_UNCLASSIFIED,"Cannot open shell");
                    iv_launched=true;
                }
            }

            // Check if algorithm must stop via stop_at
            if (++images==artPrm.stop_at)
                break;
        }

        if (!(artPrm.tell&TELL_SHOW_ERROR))
            progress_bar(artPrm.numIMG);

        // Update residual images for WLS
        if (artPrm.WLS)
        {
            double kappa=artPrm.kappa(it);
            updateResidualVector( artPrm, vol_basis, kappa,
                                  mean_error_2ndblock, pow_residual_imgs);
        }

        // Prepare for next global iteration ---------------------------------
        // Calculate norm and print
        if (rank==-1)
        {
            *artPrm.fh_hist << "Finished - Iteration " << it << std::endl;
            if (artPrm.WLS)
                std::cout        << "   Weighted error: " << global_mean_error
                << " 1st block: "        << global_mean_error_1stblock
                << " 2nd block: "        << mean_error_2ndblock
                << " residual: "         << pow_residual_imgs << std::endl;
            else
                std::cout << "   Global mean squared error: "
                <<  global_mean_error/artPrm.numIMG << std::endl;

            if (artPrm.WLS)
                *artPrm.fh_hist << "   Weighted error: " << global_mean_error
                << " 1st block: "        << global_mean_error_1stblock
                << " 2nd block: "        << mean_error_2ndblock
                << " residual: "         << pow_residual_imgs << std::endl;
            else
                *artPrm.fh_hist << "   Global mean squared error: "
                << global_mean_error/artPrm.numIMG << std::endl;

            if (POCS.apply_POCS)
            {
                std::cout        << "   POCS Global mean squared error: "
                << POCS.POCS_global_mean_error/POCS.POCS_N << std::endl;
                *artPrm.fh_hist << "   POCS Global mean squared error: "
                << POCS.POCS_global_mean_error/POCS.POCS_N << std::endl;
            }
        }

        // Convert volume and write if not last iteration
        // If in SIRT mode move the volume to the reference one
        if ((artPrm.parallel_mode==BasicARTParameters::SIRT
             && !artPrm.variability_analysis)||
            artPrm.parallel_mode==BasicARTParameters::pSIRT ||
            artPrm.parallel_mode==BasicARTParameters::pfSIRT ||
            artPrm.parallel_mode==BasicARTParameters::pSART ||
            artPrm.parallel_mode==BasicARTParameters::pCAV ||
            artPrm.eq_mode==CAV )
        {
            vol_basis=vol_basis_out;
            if (artPrm.sparseEps>0)
            {
                artPrm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
                                            Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
                selfScaleToSize(BSPLINE3,vol_voxels(),
                                (size_t)NEXT_POWER_OF_2(XSIZE(vol_voxels())),
                                (size_t)NEXT_POWER_OF_2(YSIZE(vol_voxels())),
                                (size_t)NEXT_POWER_OF_2(ZSIZE(vol_voxels())));
                Image<double> vol_wavelets, vol_wavelets_abs;
                set_DWT_type(DAUB12);
                DWT(vol_voxels(),vol_wavelets());
                vol_wavelets_abs()=vol_wavelets();
                vol_wavelets_abs().selfABS();
                double *begin=MULTIDIM_ARRAY(vol_wavelets_abs());
                double *end=MULTIDIM_ARRAY(vol_wavelets_abs())+
                            MULTIDIM_SIZE(vol_wavelets_abs());
                std::sort(begin,end);
                double threshold1=DIRECT_MULTIDIM_ELEM(vol_wavelets_abs(),
                                                       (long int)((1-artPrm.sparseEps)*MULTIDIM_SIZE(vol_wavelets_abs())));
                std::cout << "Threshold=" << threshold1 << std::endl;
                vol_wavelets().threshold("abs_below", threshold1, 0.0);
                IDWT(vol_wavelets(),vol_voxels());
                selfScaleToSize(BSPLINE3,vol_voxels(),
                                Xoutput_volume_size,
                                Youtput_volume_size,
                                Zoutput_volume_size);
                artPrm.basis.changeFromVoxels(vol_voxels(),vol_basis,artPrm.grid_type,
                                              artPrm.grid_relative_size, NULL, NULL, artPrm.R, artPrm.threads);
            }
        }

        if ( (artPrm.tell & TELL_SAVE_INTERMIDIATE) && it != artPrm.no_it-1)
        {
            if (rank==-1)
                std::cout << "Converting basis volume to voxels ...\n";
            artPrm.basis.changeToVoxels(vol_basis, &(vol_voxels()),
                                        Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
            FileName fn_tmp = artPrm.fn_out.insertBeforeExtension(formatString("_it%06d", it));
            vol_voxels.write(fn_tmp);

            if (artPrm.tell & TELL_SAVE_BASIS)
                vol_basis.write(fn_tmp.removeLastExtension().addExtension("basis"));
        }

        // Check if algorithm must stop via stop_at
        if (images==artPrm.stop_at)
            break;
    }

    // Times on screen
    if (rank==-1)
    {
        std::cout << "\nTime of " << artPrm.no_it << " iterations: \n";
        print_elapsed_time(time0);
    }

    // Finish variability analysis
    VC.finishAnalysis();

    // Save the noisy reconstruction
    if (artPrm.noisy_reconstruction)
    {
        Image<double> vol_voxels_noisy;
        artPrm.basis.changeToVoxels(vol_basis_noisy, &(vol_voxels_noisy()),
                                    Zoutput_volume_size, Youtput_volume_size, Xoutput_volume_size);
        vol_voxels_noisy.write(artPrm.fn_out.insertBeforeExtension("_noise"));
        SF_noise.write(artPrm.fn_root+"_noise_proj.xmd");
        SF_signal.write(artPrm.fn_root+"_signal_proj.xmd");
    }

}

void ARTReconsBase::singleStep(GridVolume & vol_in, GridVolume *vol_out, Projection & theo_proj, Projection & read_proj, int sym_no, Projection & diff_proj, Projection & corr_proj, Projection & alig_proj, double & mean_error, int numIMG, double lambda, int act_proj, const FileName & fn_ctf, const MultidimArray<int> *maskPtr, bool refine)
{}

void ARTReconsBase::postProcess(GridVolume &vol_basis)
{}

void ARTReconsBase::initHistory(const GridVolume &vol_basis0)
{
    // Show general information ................................................
    *artPrm.fh_hist << " ============================================================\n";
    *artPrm.fh_hist << " ART RECONSTRUCTION HISTORY: \n";
    *artPrm.fh_hist << " ============================================================\n\n";
    *artPrm.fh_hist << " Parameters -------------------------------------------------\n";
    *artPrm.fh_hist << " Projections file     : " << artPrm.fn_sel << std::endl;
    *artPrm.fh_hist << " CTF file             : " << artPrm.fn_ctf << std::endl;
    *artPrm.fh_hist << " Goldmask             : " << artPrm.goldmask << std::endl;
    *artPrm.fh_hist << " Unmatched projectors : " << artPrm.unmatched << std::endl;
    *artPrm.fh_hist << " Sampling             : " << artPrm.sampling << std::endl;
    *artPrm.fh_hist << artPrm.basis << std::endl;
    switch (artPrm.grid_type)
    {
    case FCC:
        *artPrm.fh_hist << " Grid Type: FCC ";
        break;
    case BCC:
        *artPrm.fh_hist << " Grid Type: BCC ";
        break;
    case CC:
        *artPrm.fh_hist << " Grid Type: CC ";
        break;
    }
    *artPrm.fh_hist << " Shifted tomograms:" << artPrm.shiftedTomograms << std::endl;
    *artPrm.fh_hist << " Ray length: " << artPrm.ray_length << std::endl;
    *artPrm.fh_hist << "\n Radius of the interest sphere= " << artPrm.R
    << " pixels" << std::endl;
    *artPrm.fh_hist << " Grid unit=" << artPrm.grid_relative_size
    << " pixels" << std::endl;
    if (artPrm.proj_ext==0)
        *artPrm.fh_hist << " No Projection extension\n";
    else
        *artPrm.fh_hist << " Projection extension frame: " << artPrm.proj_ext
        << " pixels\n";
    if (artPrm.Zoutput_volume_size==0)
        *artPrm.fh_hist << " Output volume size as input images\n";
    else
        *artPrm.fh_hist << " Output volume size (ZxYxX): "
        << artPrm.Zoutput_volume_size << "x" << artPrm.Youtput_volume_size
        << "x" << artPrm.Xoutput_volume_size << std::endl;
    *artPrm.fh_hist << " Iterations:    lambda=" << artPrm.lambda_list << std::endl
    << "                No. Iter=" << artPrm.no_it << std::endl;
    if (artPrm.WLS)
    {
        *artPrm.fh_hist << " Perform weighted least-squares ART "<< std::endl;
        *artPrm.fh_hist << " Iterations:    kappa=" << artPrm.kappa_list << std::endl
        << "                No. Iter=" << artPrm.no_it << std::endl;
    }
    *artPrm.fh_hist << " Parallel mode: ";
    switch (artPrm.parallel_mode)
    {
    case BasicARTParameters::ART:
        *artPrm.fh_hist << "ART\n";
        break;
    case BasicARTParameters::SIRT:
        *artPrm.fh_hist << "SIRT\n";
        break;
    case BasicARTParameters::pSIRT:
        *artPrm.fh_hist << "Parallel SIRT\n";
        break;
    case BasicARTParameters::pfSIRT:
        *artPrm.fh_hist << "Parallel False SIRT\n";
        break;
    case BasicARTParameters::pSART:
        *artPrm.fh_hist << "SART, block size=" << artPrm.block_size << std::endl;
        break;
    case BasicARTParameters::pAVSP:
        *artPrm.fh_hist << "AVSP\n";
        break;
    case BasicARTParameters::pBiCAV:
        *artPrm.fh_hist << "BiCAV, block size=" << artPrm.block_size << std::endl;
        break;
    case BasicARTParameters::pCAV:
        *artPrm.fh_hist << "CAV (global algorithm)\n";
        break;
    }
    *artPrm.fh_hist << " Equation mode: ";
    if (artPrm.eq_mode==CAV)
        *artPrm.fh_hist << "CAV (Projection update)\n";
    else if (artPrm.eq_mode==CAVK)
        *artPrm.fh_hist << "CAVK\n";
    else if (artPrm.eq_mode==CAVARTK)
        *artPrm.fh_hist << "CAVARTK\n";
    else
        *artPrm.fh_hist << "ARTK\n";
    *artPrm.fh_hist << " Projections: Ydim x Xdim = " << artPrm.projYdim
    << " x " << artPrm.projXdim << std::endl;
    *artPrm.fh_hist << " Surface mask: " << artPrm.fn_surface_mask << std::endl;
    *artPrm.fh_hist << " POCS freq: " << artPrm.POCS_freq << std::endl;
    *artPrm.fh_hist << " Known volume: " << artPrm.known_volume << std::endl;
    *artPrm.fh_hist << " Positivity: " <<artPrm.positivity << std::endl;
    *artPrm.fh_hist << " Apply shifts in images headers: " << artPrm.apply_shifts << std::endl;
    *artPrm.fh_hist << " Symmetry file: " << artPrm.fn_sym << std::endl;
    *artPrm.fh_hist << " Force Symmetry each: " << artPrm.sym_each << " projections"
    <<std::endl;
    *artPrm.fh_hist << " Maximun absolute tilt angle: " << artPrm.max_tilt << " degrees"
    <<std::endl;
    *artPrm.fh_hist << " Generating symmetry group: " << !artPrm.do_not_generate_subgroup << std::endl;
    *artPrm.fh_hist << " Do not use symmetrized projections: "
    << !artPrm.do_not_use_symproj << std::endl;
    *artPrm.fh_hist << " Number of total projections (including symmetrized): "
    << artPrm.numIMG << std::endl;
    *artPrm.fh_hist << " Number of different projections: " << artPrm.trueIMG << std::endl;
    *artPrm.fh_hist << " Forcing symmetry: " << artPrm.force_sym << std::endl;
    *artPrm.fh_hist << " Stop at: " << artPrm.stop_at << std::endl;
    if (artPrm.random_sort)
        *artPrm.fh_hist << " Random sort" << std::endl;
    else
        *artPrm.fh_hist << " Sort with last " << artPrm.sort_last_N << " images\n";
    *artPrm.fh_hist << " Variability analysis: " << artPrm.variability_analysis << std::endl;
    *artPrm.fh_hist << " Refine projections: " << artPrm.refine << std::endl;
    *artPrm.fh_hist << " Sparsity epsilon: " << artPrm.sparseEps << std::endl;
    *artPrm.fh_hist << " Diffusion weight: " << artPrm.diffusionWeight << std::endl;
    if (artPrm.SL.symsNo()!=0)
    {
        Matrix2D<double> L(4,4),R(4,4); // A matrix from the list
        *artPrm.fh_hist << " Symmetry matrices -------\n";
        for (int j=0; j<artPrm.SL.symsNo(); j++)
        {
            artPrm.SL.getMatrices(j,L,R);
            *artPrm.fh_hist << " Left  Symmetry matrix " << j << std::endl << L;
            *artPrm.fh_hist << " Right Symmetry matrix " << j << std::endl << R << std::endl;
        }
    }
    *artPrm.fh_hist << " Saving intermidiate at every "
    << artPrm.save_intermidiate_every << " projections\n";
    if(artPrm.ref_trans_step > 0.)
    {
        *artPrm.fh_hist << " Refine translational alignement after "
        << artPrm.ref_trans_after << " projection presentations\n"
        << " Maximun allowed shift is: " << artPrm.ref_trans_step << "\n\n";
    }

    // Show angles .............................................................
    // Prepare info structure for showing
    MetaData MD;
    double dfrot, dftilt, dfpsi;
    size_t id;
    for (int i=0; i<artPrm.numIMG; i++)
    {
        id=MD.addObject();
        MD.setValue(MDL_ANGLE_ROT, artPrm.IMG_Inf[i].rot,id);
        MD.setValue(MDL_ANGLE_TILT, artPrm.IMG_Inf[i].tilt,id);
        MD.setValue(MDL_ANGLE_PSI, artPrm.IMG_Inf[i].psi,id);
        MD.setValue(MDL_SYMNO, artPrm.IMG_Inf[i].sym,id);
    }

    // Now show
    *artPrm.fh_hist << " Projection angles -----------------------------------------\n";
    FOR_ALL_OBJECTS_IN_METADATA(MD)
    {
        MD.getValue(MDL_ANGLE_ROT, dfrot,__iter.objId);
        MD.getValue(MDL_ANGLE_TILT, dftilt,__iter.objId);
        MD.getValue(MDL_ANGLE_PSI, dfpsi,__iter.objId);

        *artPrm.fh_hist << "rot= "<<dfrot<<" tilt= "<<dftilt<<" psi= "<<dfpsi<<std::endl;
    }
    *artPrm.fh_hist << " -----------------------------------------------------------\n";

    // Show Initial volume and volume structure ................................
    if (artPrm.fn_start != "")
        *artPrm.fh_hist << " Starting from file: " << artPrm.fn_start << std::endl;
    else
        *artPrm.fh_hist << " Starting from a zero volume\n";

    *artPrm.fh_hist << " Grid structure ------\n" << vol_basis0.grid();

    // Show extra information ..................................................
    if (typeid(*this) != typeid(ARTReconsBase))
        *artPrm.fh_hist << " Extra: ";
    print(*artPrm.fh_hist);
}

void ARTReconsBase::print(std::ostream &o) const
{}
std::ostream & operator<< (std::ostream &o, const ARTReconsBase& artRecons)
{
    artRecons.print(o);
    return o;
}

void ARTReconsBase::applySymmetry(GridVolume &vol_in, GridVolume *vol_out,int grid_type)
{
    REPORT_ERROR(ERR_NOT_IMPLEMENTED, "ARTReconsBase::applySymmetry: Function not implemented for single particles");
}


void SinPartARTRecons::preProcess(GridVolume & vol_basis0, int level, int rank)
{
    ARTReconsBase::preProcess(vol_basis0, level, rank);

    // As this is a threaded implementation, create structures for threads, and
    // create threads
    if( artPrm.threads > 1 )
    {
        th_ids = (pthread_t *)malloc( artPrm.threads * sizeof( pthread_t));

        // Initialize the structures which will contain the parameters passed to different
        // threads
        project_threads = (project_thread_params *) malloc ( artPrm.threads * sizeof( project_thread_params ) );

        // Initialize barrier to wait for working threads and the master thread.
        barrier_init( &project_barrier, (artPrm.threads+1) );

        // Threads are created in a waiting state. They can only run when master thread unlocks them by calling
        // barrier_wait()

        for( int c = 0 ; c < artPrm.threads ; c++ )
        {
            project_threads[c].thread_id = c;
            project_threads[c].threads_count = artPrm.threads;
            project_threads[c].destroy = false;

            pthread_create( (th_ids+c), NULL, project_SimpleGridThread<double>, (void *)(project_threads+c) );
        }
    }
}


void SinPartARTRecons::singleStep(GridVolume &vol_in, GridVolume *vol_out,
                                  Projection &theo_proj, Projection &read_proj,
                                  int sym_no,
                                  Projection &diff_proj, Projection &corr_proj, Projection &alig_proj,
                                  double &mean_error, int numIMG, double lambda, int act_proj,
                                  const FileName &fn_ctf, const MultidimArray<int> *maskPtr,
                                  bool refine)
{
    // Prepare to work with CTF ................................................
    FourierFilter ctf;
    ImageOver *footprint = (ImageOver *) & artPrm.basis.blobprint;
    ImageOver *footprint2 = (ImageOver *) & artPrm.basis.blobprint2;
    bool remove_footprints = false;
    double weight, sqrtweight;

    if (fn_ctf != "" && !artPrm.unmatched)
    {
        if (artPrm.basis.type != Basis::blobs)
            REPORT_ERROR(ERR_VALUE_INCORRECT, "ART_single_step: This way of correcting for the CTF "
                         "only works with blobs");
        // It is a description of the CTF
        ctf.FilterShape = ctf.FilterBand = CTF;
        ctf.ctf.enable_CTFnoise = false;
        ctf.ctf.read(fn_ctf);
        ctf.ctf.Tm /= BLOB_SUBSAMPLING;
        ctf.ctf.produceSideInfo();

        // Create new footprints
        footprint = new ImageOver;
        footprint2 = new ImageOver;
        remove_footprints = true;

        // Enlarge footprint, bigger than necessary to avoid
        // aliasing
        *footprint = artPrm.basis.blobprint;
        (*footprint)().setXmippOrigin();
        double blob_radius = artPrm.basis.blob.radius;
        int finalsize = 2 * CEIL(30 + blob_radius) + 1;
        footprint->window(
            FIRST_XMIPP_INDEX(finalsize), FIRST_XMIPP_INDEX(finalsize),
            LAST_XMIPP_INDEX(finalsize), LAST_XMIPP_INDEX(finalsize));

        // Generate mask to the size of the footprint, correct phase
        // and apply CTF
        ctf.generateMask((*footprint)());
        ctf.correctPhase();
        ctf.applyMaskSpace((*footprint)());

        // Remove unnecessary regions
        finalsize = 2 * CEIL(15 + blob_radius) + 1;
        footprint->window(
            FIRST_XMIPP_INDEX(finalsize), FIRST_XMIPP_INDEX(finalsize),
            LAST_XMIPP_INDEX(finalsize), LAST_XMIPP_INDEX(finalsize));
#ifdef DEBUG

        Image<double> save;
        save() = (*footprint)();
        save.write("PPPfootprint.xmp");
#endif

        // Create footprint2
        *footprint2 = *footprint;
        (*footprint2)() *= (*footprint2)();
    }

    // Project structure .......................................................
    // The correction image is reused in this call to store the normalising
    // projection, ie, the projection of an all-1 volume
    Matrix2D<double> *A = NULL;
    if (artPrm.print_system_matrix)
        A = new Matrix2D<double>;
    corr_proj().initZeros();

    project_GridVolume(vol_in, artPrm.basis, theo_proj,
                       corr_proj, YSIZE(read_proj()), XSIZE(read_proj()),
                       read_proj.rot(), read_proj.tilt(), read_proj.psi(), FORWARD, artPrm.eq_mode,
                       artPrm.GVNeq, A, maskPtr, artPrm.ray_length, artPrm.threads);

    if (fn_ctf != "" && artPrm.unmatched)
    {
        ctf.generateMask(theo_proj());
        ctf.applyMaskSpace(theo_proj());
    }

    // Print system matrix
    if (artPrm.print_system_matrix)
    {
        std::cout << "Equation system (Ax=b) ----------------------\n";
        std::cout << "Size: "<< (*A).mdimx <<"x"<<(*A).mdimy<< std::endl;
        for (size_t i = 0; i < MAT_YSIZE(*A); i++)
        {
            bool null_row = true;
            for (size_t j = 0; j < MAT_XSIZE(*A); j++)
                if (MAT_ELEM(*A, i, j) != 0)
                {
                    null_row = false;
                    break;
                }
            if (!null_row)
            {
                std::cout << "pixel=" << integerToString(i, 3) << " --> "
                << DIRECT_MULTIDIM_ELEM(read_proj(), i) << " = ";
                for (size_t j = 0; j < MAT_XSIZE(*A); j++)
                    std::cout << MAT_ELEM(*A, i, j) << " ";
                std::cout << std::endl;
            }
        }
        std::cout << "---------------------------------------------\n";
        delete A;
    }

    // Refine ..............................................................
    if (refine)
    {
        Matrix2D<double> M;
        /*
        Image<double> save;
        save()=theo_proj(); save.write("PPPtheo.xmp");
        save()=read_proj(); save.write("PPPread.xmp");
        */
        alignImages(theo_proj(),read_proj(),M);
        //save()=read_proj(); save.write("PPPread_aligned.xmp");
        std::cout << M << std::endl;
        read_proj().rangeAdjust(theo_proj());
        //save()=read_proj(); save.write("PPPread_aligned_grey.xmp");
        //std::cout << "Press any key\n";
        //char c; std::cin >> c;
    }

    // Now compute differences .................................................
    double applied_lambda = lambda / numIMG; // In ART mode, numIMG=1

    mean_error = 0;
    diff_proj().resize(read_proj());

    // Weighted least-squares ART for Maximum-Likelihood refinement
    if (artPrm.WLS)
    {
        weight = read_proj.weight() / artPrm.sum_weight;
        sqrtweight = sqrt(weight);

        FOR_ALL_ELEMENTS_IN_ARRAY2D(IMGMATRIX(read_proj))
        {
            // Compute difference image and error
            IMGPIXEL(diff_proj, i, j) = IMGPIXEL(read_proj, i, j) - IMGPIXEL(theo_proj, i, j);
            mean_error += IMGPIXEL(diff_proj, i, j) * IMGPIXEL(diff_proj, i, j);

            // Subtract the residual image (stored in alig_proj!)
            IMGPIXEL(diff_proj, i, j) = sqrtweight * IMGPIXEL(diff_proj, i, j) - IMGPIXEL(alig_proj, i, j);

            // Calculate the correction and the updated residual images
            IMGPIXEL(corr_proj, i, j) =
                applied_lambda * IMGPIXEL(diff_proj, i, j) / (weight * IMGPIXEL(corr_proj, i, j) + 1.);
            IMGPIXEL(alig_proj, i, j) += IMGPIXEL(corr_proj, i, j);
            IMGPIXEL(corr_proj, i, j) *= sqrtweight;

        }
        mean_error /= XSIZE(diff_proj()) * YSIZE(diff_proj());
        mean_error *= weight;

    }
    else
    {
        long int Nmean=0;
        FOR_ALL_ELEMENTS_IN_ARRAY2D(IMGMATRIX(read_proj))
        {
            if (maskPtr!=NULL)
                if ((*maskPtr)(i,j)<0.5)
                    continue;
            // Compute difference image and error
            IMGPIXEL(diff_proj, i, j) = IMGPIXEL(read_proj, i, j) - IMGPIXEL(theo_proj, i, j);
            mean_error += IMGPIXEL(diff_proj, i, j) * IMGPIXEL(diff_proj, i, j);
            Nmean++;

            // Compute the correction image
            IMGPIXEL(corr_proj, i, j) = XMIPP_MAX(IMGPIXEL(corr_proj, i, j), 1);
            IMGPIXEL(corr_proj, i, j) =
                applied_lambda * IMGPIXEL(diff_proj, i, j) / IMGPIXEL(corr_proj, i, j);
        }
        mean_error /= Nmean;
    }

    // Backprojection of correction plane ......................................
    project_GridVolume(*vol_out, artPrm.basis, theo_proj,
                       corr_proj, YSIZE(read_proj()), XSIZE(read_proj()),
                       read_proj.rot(), read_proj.tilt(), read_proj.psi(), BACKWARD, artPrm.eq_mode,
                       artPrm.GVNeq, NULL, maskPtr, artPrm.ray_length, artPrm.threads);

    // Remove footprints if necessary
    if (remove_footprints)
    {
        delete footprint;
        delete footprint2;
    }
}


void SinPartARTRecons::postProcess(GridVolume & vol_basis)
{
    // Destroy created threads. This is done in a tricky way. At this point, threads
    // are "slept" waiting for a barrier to be reached by the master thread to continue
    // projecting/backprojecting a new projection. Here we set the flag destroy=true so
    // the threads won't process a projections but will return.
    if( artPrm.threads > 1 )
    {
        for( int c = 0 ; c < artPrm.threads ; c++ )
        {
            project_threads[c].destroy = true;
        }

        // Trigger threads "self-destruction"
        barrier_wait( &project_barrier );

        // Wait for effective threads death
        for( int c = 0 ; c < artPrm.threads ; c++ )
        {
            pthread_join(*(th_ids+c),NULL);
        }

        // Destroy barrier and mutex, as they are no longer needed.
        // Sjors, 27march2009
        // NO, dont do this, because running a second threaded-ART
        // in a single program will not have mutexes anymore...
        //pthread_mutex_destroy( &project_mutex );
        //barrier_destroy( &project_barrier );
    }
}



