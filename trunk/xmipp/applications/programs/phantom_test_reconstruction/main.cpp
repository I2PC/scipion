/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (1999)
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

#include <fstream>

#include <reconstruction/phantom_test_reconstruction.h>

void Usage();
void test_training(Recons_test_Parameters &prm, int nvol,
                   const FileName &fn_hist, const std::string &training_FOM);
void test_all_FOMs(Recons_test_Parameters &prm, int nvol,
                   const FileName &fn_hist);

int main(int argc, char *argv[])
{
    FileName                fn_test_params, fn_root;
    bool                    only_training;
    std::string             training_FOM;
    bool                    dont_rewrite;
    Recons_test_Parameters  recons_prm;
    int                     nvol;

// Check the command line ..................................................
    try
    {
        fn_test_params = getParameter(argc, argv, "-i");
        only_training  = checkParameter(argc, argv, "-training");
        training_FOM   = getParameter(argc, argv, "-training", "scL21");
        dont_rewrite   = checkParameter(argc, argv, "-dont_rewrite");

        recons_prm.read(fn_test_params);

        if (dont_rewrite) nvol = 1;
        else              nvol = -1;
    }
    catch (Xmipp_error &XE)
    {
        std::cout << XE;
        Usage();
        exit(1);
    }
    std::cout << recons_prm << std::endl;

// Make tests ..............................................................
    try
    {
        fn_root = fn_test_params.without_extension();
        randomize_random_generator();
        if (only_training) test_training(recons_prm, nvol, fn_root + ".hist", training_FOM);
        else               test_all_FOMs(recons_prm, nvol, fn_root + ".hist");
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
    exit(0);
}

/* Show results ------------------------------------------------------------ */
void show_measured_point(std::ostream &out, Recons_test_Parameters prm, int i)
{
    out << "Test " << i;
    if (prm.recons_method == use_WBP)
        out << ": threshold= " << prm.WBP_threshold[i] << " ";
    else
    {
        if (!prm.succesive_params)
        {
            out << ": lambda= " << prm.lambda0[i]
                << "  no_it= " << prm.no_it0[i];
            if (prm.lambda0[i] != prm.lambdaF[i] || prm.no_it0[i] != prm.no_itF[i])
                out << " lambdaF= " << prm.lambdaF[i] << " no_itF= "
                    << prm.no_itF[i] << " ";
        }
        else
        {
            for (int j = 0; j < prm.lambda0.size(); j++)
            {
                out << "\n    lambda= " << prm.lambda0[j]
                    << "  no_it= " << prm.no_it0[j] << " ";
                if (prm.lambda0[j] != prm.lambdaF[j]
                    || prm.no_it0[j] != prm.no_itF[j])
                    out << " lambdaF= " << prm.lambdaF[j] << " no_itF= "
                        << prm.no_itF[j] << " ";
            }
            out << "\n    ";
        }
    }
}

/* Tests on SCL2 ----------------------------------------------------------- */
void test_training(Recons_test_Parameters &prm, int nvol,
                   const FileName &fn_hist, const std::string &training_FOM)
{
    Matrix1D<double> training_avgs;
    Matrix1D<double> training_devs;
    Matrix1D<double> training_N;
    EVALUATE_results results;

    // Open History file
    std::ofstream fh_hist;
    fh_hist.open(fn_hist.c_str(), std::ios::out);
    if (!fh_hist)
        REPORT_ERROR(1, "recons_test: Cannot open file " + fn_hist + " for output");

    // Some initialisation
    int TestNo;
    if (prm.recons_method == use_WBP) TestNo = prm.WBP_threshold.size();
    else if (prm.recons_method == use_FOURIER) TestNo = prm.MeasNo;
    else if (!prm.succesive_params) TestNo = prm.lambda0.size();
    else TestNo = 1;
    training_avgs.resize(TestNo);
    training_devs.resize(TestNo);
    training_N.resize(TestNo);

    // Perform measures
    fh_hist << "Training on " << training_FOM << std::endl;
    for (int i = 0; i < TestNo; i++)
    {
        single_measure_on_FOM(prm, i, nvol,
                              training_avgs(i), training_devs(i), training_N(i),
                              results, training_FOM);
        show_measured_point(fh_hist, prm, i);
        if (!prm.evaluate)
            fh_hist << " ---> " << training_avgs(i) << "+-" << training_devs(i)
                    << " (" << training_N(i) << ")";
        fh_hist << "std::endl";
        fh_hist.flush();
    }
    fh_hist.close();

    // Show results
    if (prm.evaluate)
    {
        std::cout << "****************************************************************\n";
        std::cout << "Reconstruction test results\n";
        std::cout << "****************************************************************\n";
        for (int i = 0; i < TestNo; i++)
        {
            show_measured_point(std::cout, prm, i);
            std::cout << " ---> " << training_avgs(i) << "+-" << training_devs(i)
                      << " (" << training_N(i) << ")";
            std::cout << std::endl;
        }
    }
}

/* Tests all FOMs ---------------------------------------------------------- */
void test_all_FOMs(Recons_test_Parameters &prm, int nvol,
                   const FileName &fn_hist)
{
    EVALUATE_results results;

    // Open History file
    std::ofstream fh_hist;
    fh_hist.open(fn_hist.c_str(), std::ios::out);
    if (!fh_hist)
        REPORT_ERROR(1, "recons_test: Cannot open file " + fn_hist + " for output");

    // Some initialisation
    int TestNo;
    if (prm.recons_method == use_WBP) TestNo = prm.WBP_threshold.size();
    else if (prm.recons_method == use_FOURIER) TestNo = prm.MeasNo;
    else if (!prm.succesive_params) TestNo = prm.lambda0.size();
    else TestNo = 1;
    FOMs foms_mean(TestNo), foms_stddev(TestNo);

    // Perform measures
    for (int i = 0; i < TestNo; i++)
    {
        single_measure_on_all_FOMs(prm, i, nvol, foms_mean, foms_stddev,
                                   results);
        show_measured_point(fh_hist, prm, i);
        std::cout << std::endl;
        if (prm.evaluate)
            show_stats(fh_hist, i, foms_mean, foms_stddev);
    }
    fh_hist.close();

    // Show results
    if (prm.evaluate)
    {
        std::cout << "****************************************************************\n";
        std::cout << "Reconstruction test results\n";
        std::cout << "****************************************************************\n";
        for (int i = 0; i < TestNo; i++)
        {
            show_measured_point(std::cout, prm, i);
            std::cout << std::endl;
            show_stats(std::cout, i, foms_mean, foms_stddev);
        }
    }
}

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    printf("Error in the arguments\n");
    printf("Usage: \n");
    printf("       recons_test <options>\n");
    printf("-i <test params>         : a file with the test parameters\n");
    printf("                           see the Manual help for more information\n");
    printf("[-training]              : tests are performed only on scL2(1)\n");
    printf("[-dont_rewrite]          : tests are performed using different filenames\n");
}
