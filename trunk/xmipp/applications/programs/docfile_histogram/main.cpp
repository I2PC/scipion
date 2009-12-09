/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (1999)
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
#include <data/histogram.h>
#include <data/docfile.h>

void Usage();

int main(int argc, char *argv[])
{
    DocFile          DF;
    FileName         fn_out;
    FileName         fn_img;
    int              generate_img;
    int              col, col2;         // Columns for histogram
    double           m, M, m2, M2;      // range for histogram
    int              automatic_range, automatic_range2;
    int              StepsNo, StepsNo2;
    histogram1D      hist;
    histogram2D      hist2;
    Matrix1D<double> C, C2;
    double           percentil;

// Check the command line ==================================================
    try
    {
        int i;
        DF.read(getParameter(argc, argv, "-i"));
        fn_out = getParameter(argc, argv, "-o", "");
        percentil = textToFloat(getParameter(argc, argv, "-percentil", "50"));
        col = textToInteger(getParameter(argc, argv, "-col", "0"));
        col2 = textToInteger(getParameter(argc, argv, "-col2", "-1"));

        StepsNo = textToInteger(getParameter(argc, argv, "-steps", "100"));
        if ((i = paremeterPosition(argc, argv, "-range")) != -1)
        {
            if (i + 2 >= argc)
                EXIT_ERROR(1, "DocFile Histogram: Not enough parameters behind -range\n");
            m = textToFloat(argv[i+1]);
            M = textToFloat(argv[i+2]);
            automatic_range = false;
        }
        else automatic_range = true;

        StepsNo2 = textToInteger(getParameter(argc, argv, "-steps2", "100"));
        if ((i = paremeterPosition(argc, argv, "-range2")) != -1)
        {
            if (i + 2 >= argc)
                EXIT_ERROR(1, "DocFile Histogram: Not enough parameters behind -range2\n");
            m2 = textToFloat(argv[i+1]);
            M2 = textToFloat(argv[i+2]);
            automatic_range2 = false;
        }
        else automatic_range2 = true;

        // Check columns are possible
        if (col < 0 || col >= DF.FirstLine_colNumber())
            EXIT_ERROR(1, "DocFile Histogram: Column for histogram not valid");
        if ((col2 < 0 || col2 >= DF.FirstLine_colNumber()) && col2 != -1)
            EXIT_ERROR(1, "DocFile Histogram: Column 2 for histogram not valid");

        // Check if the 2D histogram must be an image
        if (col2 != -1)
        {
            fn_img = getParameter(argc, argv, "-img", "");
            generate_img = (fn_img != "");
        }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
        Usage();
        exit(1);
    }

// Compute Histogram =======================================================
    try
    {
        double avg, stddev, dummy;
        // 1D histograms -----------------------------------------------------
        if (col2 == -1)
        {
            C = DF.col(col);
            if (automatic_range) C.computeDoubleMinMax(m, M);
            compute_hist(C, hist, m, M, StepsNo);
            std::cerr << "Min: " << m << " max: " << M
            << " Steps: " << StepsNo << std::endl;
            C.computeStats(avg, stddev, dummy, dummy);
            std::cerr << "Mean: " << avg << " Stddev: " << stddev << std::endl;
            std::cerr << "Percentil (" << percentil << "): "
            << hist.percentil(percentil) << std::endl;
            if (fn_out != "") hist.write(fn_out);
            else            std::cout << hist;

            // 2D histograms -----------------------------------------------------
        }
        else
        {
            C = DF.col(col);
            C2 = DF.col(col2);
            if (automatic_range)
            {
                C.computeDoubleMinMax(m, M);
                C2.computeDoubleMinMax(m2, M2);
            }
            compute_hist(C, C2, hist2, m, M, m2, M2, StepsNo, StepsNo2);
            std::cerr << "Min1: "   << m        << " max1: " << M
            << " Steps1: " << StepsNo  << std::endl;
            C.computeStats(avg, stddev, dummy, dummy);
            std::cerr << "Mean: " << avg << " Stddev: " << stddev << std::endl;
            std::cerr << "Min2: "   << m2       << " max2: " << M2
            << " Steps2: " << StepsNo2 << std::endl;
            C2.computeStats(avg, stddev, dummy, dummy);
            std::cerr << "Mean: " << avg << " Stddev: " << stddev << std::endl;
            if (fn_out != "") hist2.write(fn_out);
            else            std::cout << hist2;
            if (generate_img)
            {
                ImageXmipp I;
                I() = hist2;
                I.write(fn_img);
            }
        }
    }
    catch (Xmipp_error XE)
    {
        std::cout << XE;
    }
}
/* Usage ------------------------------------------------------------------- */
void Usage()
{
    printf("Purpose:\n");
    printf("   This program performs the histogram of a column or a \n"
           "   couple of columns in a document file\n");
    printf("Usage:\n");
    printf("   docfile_histogram <options>\n");
    printf("   Where <options> are:\n");
    printf("      -i <docfile>                 : Spider document file\n"
           "      [-o <text file>]             : Output histogram\n"
           "      [-img <img_file>]            : only for 2D histograms\n"
           "      [-col <n=0>]                 : by default, histogram1D for \n"
           "                                     the first column\n"
           "      [-range <m> <M>]             : range for the first column\n"
           "                                     by default, it is automatic\n"
           "      [-steps <N=100>]             : number of subdivisions for column 1\n"
           "      [-percentil <p=50>]          : Only for 1D histograms\n"
           "      [-col2 <n>]                  : if specified, then a 2D histogram\n"
           "                                     is calculated\n"
           "      [-range2 <m> <M>]            : range for the second column\n"
           "                                     by default, it is automatic\n"
           "      [-steps2 <N=100>]            : number of subdivisions for column 2\n"
          );
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM DocFile_Histogram {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/DocFile_Histogram/Help/docfile_histogram.html";
      help="Computes 1D and 2D histograms of Document Files";
      OPEN MENU DocFile_Histogram;
      COMMAND LINES {
         + usual: docfile_histogram -i $DOCFILE [-o $TEXTFILE] [-img $IMGFILE]
            [-col $N1]  [-range $MIN1 $MAX1]  [-steps $STEPS1]
            [-col2 $N2] [-range2 $MIN2 $MAX2] [-steps2 $STEPS2]
      }
      PARAMETER DEFINITIONS {
         $DOCFILE {label="Input Document File"; type=FILE EXISTING;}
         $TEXTFILE {label="Output histogram"; type=FILE;}
         $IMGFILE {label="Output image";
                   help="Only for 2D histograms"; type=FILE;}
         $N1 {label="Column 1"; type=NATURAL;}
         OPT(-range) {label="Range for the first column";}
            $MIN1 {label="Minimum"; type=FLOAT;}
            $MAX1 {label="Maximum"; type=FLOAT;}
         $STEPS1 {label="Steps for column 1"; type=NATURAL;
            by default=100;}
         $N2 {label="Column 2"; type=NATURAL;}
         OPT(-range2) {label="Range for the second column";}
            $MIN2 {label="Minimum"; type=FLOAT;}
            $MAX2 {label="Maximum"; type=FLOAT;}
         $STEPS2 {label="Steps for column 2"; type=NATURAL;
            by default=100;}
      }
   }
   MENU DocFile_Histogram {
      "I/O Parameters"
      $DOCFILE
      OPT(-o)
      "First Column Parameters"
      OPT(-col)
      OPT(-range)
      OPT(-steps)
      "2D histograms"
      OPT(-img)
      OPT(-col2)
      OPT(-range2)
      OPT(-steps2)
   }
*/
