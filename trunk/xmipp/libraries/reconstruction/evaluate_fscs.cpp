/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include <data/args.h>
#include <data/selfile.h>

#include "evaluate_fscs.h"
#include "volume_foms.h"

#include <fstream>

/* Read Evaluate parameters from command line ============================== */
void Prog_Evaluate_FSCs_Parameters::read(int argc, char **argv)
{
    fn_phantom       = getParameter(argc, argv, "-p");
    fn_recons        = getParameter(argc, argv, "-r");
    fn_out           = getParameter(argc, argv, "-o", "");
    sampling_rate    = AtoF(getParameter(argc, argv, "-sampling_rate", "1"));
    action           = ESTIMATE_SINGLE_FSC;
    if (checkParameter(argc, argv, "-estimate_average_resolution"))
        action = ESTIMATE_AVERAGE_RESOLUTION;
    else if (checkParameter(argc, argv, "-estimate_average_FSC"))
        action = ESTIMATE_AVERAGE_FSC;
    else if (checkParameter(argc, argv, "-compare_two_sets"))
    {
        action = COMPARE_TWO_SETS;
        fn_recons2 = getParameter(argc, argv, "-compare_two_sets");
    }
}

/* Evaluate usage ========================================================== */
void Prog_Evaluate_FSCs_Parameters::usage()
{
    cerr << "Usage: evaluate_FSCs\n"
    << "   -p <phantom filename>          : Xmipp volume\n"
    << "   -r <reconstruction filename>   : Xmipp volume or selfile with volumes\n"
    << "  [-o <output filename>]          : for the FSC\n"
    << "  [-sampling_rate <Tm=1>]         : Angstroms/pixel\n"
    << "  [-estimate_average_resolution]  : provide a selfile for -r\n"
    << "  [-estimate_average_FSC]         : provide a selfile for -r\n"
    << "  [-compare_two_sets <selfile>]   : second set of volumes\n"
    ;
}

/* Show parameters ========================================================= */
ostream & operator << (ostream &out, const Prog_Evaluate_FSCs_Parameters &prm)
{
    out << "Phantom         : " << prm.fn_phantom        << endl
    << "Reconstruction  : " << prm.fn_recons         << endl
    << "Output file     : " << prm.fn_out            << endl
    << "Sampling rate   : " << prm.sampling_rate     << endl
    ;
    out << "Action          : ";
    switch (prm.action)
    {
    case ESTIMATE_SINGLE_FSC:
        out << "Estimate single FSC\n";
        break;
    case ESTIMATE_AVERAGE_RESOLUTION:
        out << "Estimate average resolution\n";
        break;
    case ESTIMATE_AVERAGE_FSC:
        out << "Estimate average FSC\n";
        break;
    case COMPARE_TWO_SETS:
        out << "Compare two sets: "
        << prm.fn_recons2 << endl;
        break;
    }
    return out;
}

/* Produce side information ================================================ */
void Prog_Evaluate_FSCs_Parameters::produce_side_info()
{
    phantom.read(fn_phantom);
    switch (action)
    {
    case ESTIMATE_SINGLE_FSC:
        reconstruction.read(fn_recons);
        break;
    case ESTIMATE_AVERAGE_RESOLUTION:
    case ESTIMATE_AVERAGE_FSC:
        SF_recons.read(fn_recons);
        break;
    case COMPARE_TWO_SETS:
        SF_recons.read(fn_recons);
        SF_recons2.read(fn_recons2);
        break;
    }
}

/* Compute average resolution ============================================== */
void Prog_Evaluate_FSCs_Parameters::compute_average_resolution(
    double &avg_resol, double &stddev_resol, Matrix1D<double> &resol)
{
    SF_recons.go_first_ACTIVE();
    resol.initZeros(SF_recons.ImgNo());

    Matrix1D<double> frequency, FSC;
    int i = 0;
    cerr << "Estimating average resolution ...\n";
    init_progress_bar(XSIZE(resol));
    while (!SF_recons.eof())
    {
        reconstruction.read(SF_recons.NextImg());
        resol(i) =
            compute_FSC(phantom, reconstruction, sampling_rate, frequency, FSC);
        i++;
        progress_bar(i);
    }
    progress_bar(XSIZE(resol));
    double min_resol, max_resol;
    resol.computeStats(avg_resol, stddev_resol, min_resol, max_resol);
}

/* Compute average resolution ============================================== */
void Prog_Evaluate_FSCs_Parameters::compute_average_FSC(
    Matrix1D<double> &frequency,
    Matrix1D<double> &avg_FSC, Matrix1D<double> &min_FSC,
    Matrix1D<double> &max_FSC)
{
    SF_recons.go_first_ACTIVE();
    int N = SF_recons.ImgNo();

    Matrix1D<double> FSC;
    int n = 0;
    cerr << "Estimating average FSC ...\n";
    init_progress_bar(N);
    while (!SF_recons.eof())
    {
        reconstruction.read(SF_recons.NextImg());
        compute_FSC(phantom, reconstruction, sampling_rate, frequency, FSC);
        if (n == 0)
        {
            avg_FSC.initZeros(FSC);
            max_FSC = min_FSC = FSC;
        }
        avg_FSC   += FSC;
        FOR_ALL_ELEMENTS_IN_MATRIX1D(FSC)
        if (FSC(i) < min_FSC(i)) min_FSC(i) = FSC(i);
        else if (FSC(i) > max_FSC(i)) max_FSC(i) = FSC(i);
        progress_bar(++n);
    }
    progress_bar(N);

    avg_FSC /= N;
}

/* Compare two sets ======================================================== */
void Prog_Evaluate_FSCs_Parameters::compare_two_sets(
    Matrix1D<double> &frequency,
    Matrix1D<double> &avg_diff_FSC, Matrix1D<double> &stddev_diff_FSC)
{
    SF_recons.go_first_ACTIVE();
    SF_recons2.go_first_ACTIVE();
    int N = SF_recons.ImgNo();

    VolumeXmipp reconstruction2;
    Matrix1D<double> FSC1, FSC2;
    int n = 0;
    cerr << "Estimating average FSC ...\n";
    init_progress_bar(N);
    while (!SF_recons.eof())
    {
        reconstruction.read(SF_recons.NextImg());
        reconstruction2.read(SF_recons2.NextImg());
        compute_FSC(phantom, reconstruction,  sampling_rate, frequency, FSC1);
        compute_FSC(phantom, reconstruction2, sampling_rate, frequency, FSC2);
        if (n == 0)
        {
            avg_diff_FSC.initZeros(FSC1);
            stddev_diff_FSC.initZeros(FSC1);
        }
        FOR_ALL_ELEMENTS_IN_MATRIX1D(FSC1)
        {
            double diff = FSC2(i) - FSC1(i);
            avg_diff_FSC(i) += diff;
            stddev_diff_FSC(i) += diff * diff;
        }
        progress_bar(++n);
    }
    progress_bar(N);

    if (N > 0) avg_diff_FSC /= N;
    if (N > 1)
        FOR_ALL_ELEMENTS_IN_MATRIX1D(FSC1)
        stddev_diff_FSC(i) = sqrt(ABS(stddev_diff_FSC(i) / (N - 1) -
                                      N / (N - 1) * avg_diff_FSC(i) * avg_diff_FSC(i)));
}

/* Main routine ============================================================ */
void ROUT_Evaluate_FSCs(Prog_Evaluate_FSCs_Parameters &prm)
{
    Matrix1D<double> frequency, FSC, min_FSC, max_FSC, stddev_FSC;
    double avg_resol, stddev_resol, resolution, avg_FSC;
    int i;
    ofstream fh_out;

    switch (prm.action)
    {
        // Estimate single FSC
    case ESTIMATE_SINGLE_FSC:
        resolution = compute_FSC(prm.phantom, prm.reconstruction,
                                 prm.sampling_rate, frequency, FSC);
        cout << "The resolution (FSC<0.5) of this volume is "
        << resolution << endl;

        // Compute the average FSC until the resolution
        avg_FSC = 0;
        i = 0;
        while (frequency(i) < resolution) avg_FSC += FSC(i++);
        if (i != 0) avg_FSC /= i;
        cout << "The average FSC until the resolution is "
        << avg_FSC << endl;

        // Write the FSC
        fh_out.open(prm.fn_out.c_str());
        if (!fh_out)
            REPORT_ERROR(1, (string)"Evaluate FSCs: Cannot open file" +
                         prm.fn_out + " for output");
        fh_out << "# Freq(1/A) FSC\n";
        FOR_ALL_ELEMENTS_IN_MATRIX1D(frequency)
        fh_out << frequency(i) << " " << FSC(i) << endl;
        fh_out.close();
        break;

        // Estimate average resolution
    case ESTIMATE_AVERAGE_RESOLUTION:
        prm.compute_average_resolution(avg_resol, stddev_resol, FSC);
        cout << "Average resolution=" << avg_resol << "+-" << stddev_resol
        << endl << endl;
        cout << "Vector of all resolutions involved=\n" << FSC.transpose()
        << endl;
        break;

        // Estimate average resolution
    case ESTIMATE_AVERAGE_FSC:
        prm.compute_average_FSC(frequency, FSC, min_FSC, max_FSC);

        // Write the FSC
        fh_out.open(prm.fn_out.c_str());
        if (!fh_out)
            REPORT_ERROR(1, (string)"Evaluate FSCs: Cannot open file" +
                         prm.fn_out + " for output");
        fh_out << "# Freq(1/A) FSC min_FSC max_FSC\n";
        FOR_ALL_ELEMENTS_IN_MATRIX1D(frequency)
        fh_out << frequency(i) << " " << FSC(i) << " "
        << min_FSC(i) << " " << max_FSC(i) << endl;
        fh_out.close();
        break;

        // Compare two sets
    case COMPARE_TWO_SETS:
        prm.compare_two_sets(frequency, FSC, stddev_FSC);

        // Write the FSC
        fh_out.open(prm.fn_out.c_str());
        if (!fh_out)
            REPORT_ERROR(1, (string)"Evaluate FSCs: Cannot open file" +
                         prm.fn_out + " for output");
        fh_out << "# Freq(1/A) diff_FSC stddev_diff_FSC avg-stddev avg+stddev\n";
        FOR_ALL_ELEMENTS_IN_MATRIX1D(frequency)
        fh_out << frequency(i) << " " << FSC(i) << " "
        << stddev_FSC(i) << " " << FSC(i) - stddev_FSC(i) << " "
        << FSC(i) + stddev_FSC(i) << endl;
        fh_out.close();
        break;
    }
}
