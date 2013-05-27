/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2002)
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
#include <data/morphology.h>
#include <data/filters.h>
#include "volume_segment.h"

// Read arguments ==========================================================
void ProgVolumeSegment::readParams()
{
    voxel_mass=dalton_mass=aa_mass=-1;
    sampling_rate=-1;
    threshold=-1;
    wang_radius=-1;
    en_threshold=false;
    do_prob=false;
    otsu=false;

    fn_vol = getParam("-i");
    fn_mask = getParam("-o");
    method=getParam("--method");
    if (method=="voxel_mass")
        voxel_mass = getDoubleParam("--method", 1);
    else if (method=="dalton_mass")
    {
        dalton_mass = getDoubleParam("--method", 1);
        sampling_rate = getDoubleParam("--method", 2);
    }
    else if (method=="aa_mass")
    {
        aa_mass = getDoubleParam("--method", 1);
        sampling_rate = getDoubleParam("--method", 2);
    }
    else if (method=="threshold")
    {
        en_threshold = true;
        threshold = getDoubleParam("--method", 1);
    }
    else if (method=="otsu")
        otsu=true;
    else if (method=="prob")
    {
        do_prob = true;
        wang_radius = getIntParam("--method", 1);
    }
}

// Show ====================================================================
void ProgVolumeSegment::show() const
{
    std::cout
    << "Input file   : " << fn_vol        << std::endl
    << "Voxel mass   : " << voxel_mass    << std::endl
    << "Dalton mass  : " << dalton_mass   << std::endl
    << "AA mass      : " << aa_mass       << std::endl
    << "Sampling rate: " << sampling_rate << std::endl
    << "Output mask  : " << fn_mask       << std::endl
    << "Enable thres.: " << en_threshold  << std::endl
    << "Threshold    : " << threshold     << std::endl
    << "Otsu         : " << otsu          << std::endl
    << "Wang radius  : " << wang_radius   << std::endl
    << "Probabilistic: " << do_prob       << std::endl
    ;
}

// usage ===================================================================
void ProgVolumeSegment::defineParams()
{
    addParamsLine("   -i <volume>              : Volume to segment");
    addParamsLine("  [-o <mask=\"\">]          : Output mask");
    addParamsLine("   --method <method>        : Segmentation method");
    addParamsLine("      where <method>");
    addParamsLine("            voxel_mass  <mass>        : Mass in voxels");
    addParamsLine("            dalton_mass <mass> <Tm=1> : Mass in daltons");
    addParamsLine("                                      : Tm is the sampling rate (A/pixel)");
    addParamsLine("            aa_mass     <mass> <Tm=1> : Mass in aminoacids");
    addParamsLine("                                      : Tm is the sampling rate (A/pixel)");
    addParamsLine("            threshold   <th>          : Thresholding");
    addParamsLine("            otsu                      : Otsu's method segmentation");
    addParamsLine("            prob        <radius=-1>   : Probabilistic solvent mask (typical value 3)");
    addParamsLine("                                      : Radius [pix] is used for B.C. Wang cone smoothing");
}

// Produce side information ================================================
void ProgVolumeSegment::produce_side_info()
{
    V.read(fn_vol);
    if (method=="dalton_mass" || method=="aa_mass")
    {
        if ((dalton_mass == -1 && aa_mass == -1) || sampling_rate == -1)
            REPORT_ERROR(ERR_VALUE_INCORRECT, "No way to compute voxel mass");
        double sampling_rate3 = sampling_rate * sampling_rate * sampling_rate;
        if (dalton_mass != -1)
            voxel_mass = dalton_mass * 1.207 / sampling_rate3;
        else
            voxel_mass = aa_mass * 110 * 1.207 / sampling_rate3;
        std::cout << std::endl << "Derived voxel_mass=" << voxel_mass << std::endl;
    }
}

// Count voxels ============================================================
// Segment with a given threshold and compute the number of voxels in the
// biggest piece
//#define DEBUG
double segment_threshold(const Image<double> *V_in, Image<double> *V_out,
                         double threshold, bool do_prob)
{
    Image<double> aux;

    // Binarize input volume
    (*V_out)() = (*V_in)();
    (*V_out)().threshold("below", threshold, threshold);
    (*V_out)().binarize(threshold);

#ifdef DEBUG

    std::cout << threshold << std::endl;
    Image<double> save;
    save() = (*V_in)();
    save.write("PPP0.vol");
    save() = (*V_out)();
    save.write("PPP1.vol");
#endif

    if (!do_prob)
    {
        // Apply morphological opening to input volume
        aux().initZeros((*V_out)());
        opening3D((*V_out)(), aux(), 18, 0, 1);
        closing3D(aux(), (*V_out)(), 18, 0, 1);
#ifdef DEBUG

        save() = (*V_out)();
        save.write("PPP2.vol");
#endif
    }

    // Count the number of different objects
    int no_comp = labelImage3D((*V_out)(), aux());
    Matrix1D<double> count(no_comp + 1);
    const MultidimArray<double> &maux=aux();
    FOR_ALL_ELEMENTS_IN_ARRAY3D(maux)
    VEC_ELEM(count,(int)A3D_ELEM(maux,k, i, j))++;

#ifdef DEBUG

    std::cout << count << std::endl << std::endl;
    std::cout << "Press any key\n";
    char c;
    std::cin >> c;
#endif

    // Pick the maximum
    count(0) = 0; // We don't want to pick the background
    int imax;
    count.maxIndex(imax);

    // Select the mask with only that piece
    const MultidimArray<double> &mVout=(*V_out)();
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mVout)
    A3D_ELEM(mVout,k, i, j) = A3D_ELEM(maux,k, i, j) == imax;

    return count(imax);
}

void wang_smoothing(const Image<double> *V_in, Image<double> *V_out, int radius)
{

    int r2;
    double radius2 = radius * radius;
    double sumw, weight;

    (*V_out)().initZeros((*V_in)());

    const MultidimArray<double> &mVin=(*V_in)();
    const MultidimArray<double> &mVout=(*V_out)();
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mVin)
    {
        sumw = 0.;
        for (int kp = k - radius; kp < k + radius; kp++)
        {
            if (kp > STARTINGZ(mVin) && kp < FINISHINGZ(mVin))
            {
                for (int ip = i - radius; ip < i + radius; ip++)
                {
                    if (ip > STARTINGY(mVin) && ip < FINISHINGY(mVin))
                    {
                        for (int jp = j - radius; jp < j + radius; jp++)
                        {
                            if (jp > STARTINGX(mVin) && jp < FINISHINGX(mVin))
                            {
                                r2 = (kp - k) * (kp - k) + (ip - i) * (ip - i) + (jp - j) * (jp - j);
                                if ((r2 < radius2) && (A3D_ELEM(mVin, kp, ip, jp) > 0.))
                                {
                                    weight = 1. - sqrt(r2 / radius2);
                                    A3D_ELEM(mVout, k, i, j) += weight * XMIPP_MAX(0., A3D_ELEM(mVin, kp, ip, jp));
                                    sumw += weight;
                                }
                            }
                        }
                    }
                }
            }
        }
        if (sumw > 0.)
            A3D_ELEM(mVout, k, i, j) /= sumw;
    }
}


void probabilistic_solvent(Image<double> *V_in, Image<double> *V_out)
{
    // Calculate mean and sigma for protein and solvent regions
    // according to the traditional segmentation
    double Np, sump, sum2p, Ns, sums, sum2s;
    double avgp, sigp, avgs, sigs, aux, solv_frac, prot_frac;
    double p_prot, p_solv;

    MultidimArray<double> &mVin=(*V_in)();
    MultidimArray<double> &mVout=(*V_out)();
    mVin.setXmippOrigin();
    mVout.setXmippOrigin();

    Np = sump = sum2p = Ns = sums = sum2s = 0.;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mVin)
    {
        aux = A3D_ELEM(mVin, k, i, j);
        if (A3D_ELEM(mVout, k, i, j) < 0.5)
        {
            sums += aux;
            sum2s += aux * aux;
            Ns += 1.;
        }
        else
        {
            sump += aux;
            sum2p += aux * aux;
            Np += 1.;
        }
    }
    if (Np > 0. && Ns > 0.)
    {
        avgs = sums / Ns;
        sigs = sum2s / Ns - avgs * avgs;
        avgp = sump / Np;
        sigp = sum2p / Np - avgp * avgp;
        prot_frac = Np / (Np + Ns);
        solv_frac = 1. - prot_frac;
    }
    else
    {
        REPORT_ERROR(ERR_NUMERICAL, "empty solvent or protein region");
    }

    // Terwilliger-like calculation of P(x|solv) & P(x|prot)
    // Bayes: P(prot|x)= P(x|prot)/{P(x|prot)+P(x|solv)}
    double isigs=-1/(2.0*sigs);
    double isigp=-1/(2.0*sigp);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(mVin)
    {
    	double voxelVal=A3D_ELEM(mVin, k, i, j);
        aux = voxelVal - avgs;
        p_solv = solv_frac * exp(aux * aux * isigs);
        aux = voxelVal - avgp;
        p_prot = prot_frac * exp(aux * aux * isigp);
        A3D_ELEM(mVout, k, i, j) = p_prot / (p_prot + p_solv);
    }
}

// Really segment ==========================================================
void ProgVolumeSegment::segment(Image<double> &mask)
{
    double th_min, th_max, val_min, val_max;
    V().computeDoubleMinMax(val_min, val_max);
    th_min = val_min;
    th_max = val_max;

    bool ok = false;
    if (!otsu)
    {
        if (!en_threshold)
        {
            // Perform a bracketing search until the mass is
            // within a 0.1% of the desired mass
            do
            {
                double th_med = (th_min + th_max) * 0.5;
                double mass_med = segment_threshold(&V, &mask, th_med, do_prob);
                std::cout << "Threshold= " << th_med
                << " mass of the main piece= " << mass_med << std::endl;
                if (ABS(mass_med - voxel_mass) / voxel_mass < 0.001)
                {
                    ok = true;
                    break;
                }
                if ((th_max - th_min) / (val_max - val_min) < 0.0001)
                    break;
                if (mass_med < voxel_mass)
                {
                    th_max = th_med;
                }
                else
                {
                    th_min = th_med;
                }
            }
            while (true);
        }
        else
        {
            // Perform a single thresholding
            double mass_med = segment_threshold(&V, &mask, threshold, do_prob);
            std::cout << "Threshold= " << threshold
            << " mass of the main piece= " << mass_med << std::endl;
            ok = true;
        }
    }

    if (otsu)
    {
        mask()=V();
        double th=EntropyOtsuSegmentation(mask());
        std::cout << "Threshold " << th << std::endl;
        ok=true;
    }

    if (do_prob)
    {
        // Wang-Leslie like modification of the input volume
        if (wang_radius >= 3)
        {
            Image<double> Vwang;
            wang_smoothing(&V, &Vwang, wang_radius);
            V = Vwang;
        }
        // Terwilliger-like calculation of P(solv|x) through P(x|solv) & P(x|prot)
        probabilistic_solvent(&V, &mask);
    }

    // Save mask if necessary
    if (fn_mask != "" && (ok || do_prob))
        mask.write(fn_mask);
    if (!ok && !do_prob)
        std::cout << "Segment: Cannot find an appropriate threshold\n";
}

void ProgVolumeSegment::run()
{
    Image<double> mask;
    show();
    produce_side_info();
    segment(mask);
}
