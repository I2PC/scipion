/***************************************************************************
 *
 * Authors:    Sjors Scheres (scheres@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include <vector>
#include <iomanip>
#include <sstream>

#include <data/args.h>
#include <data/image.h>
#include <data/selfile.h>
#include <data/docfile.h>

void Usage();

int main(int argc, char **argv)
{
    FileName     fn_logsel, fn_sel, fn_oroot, fn_img, fn_img2, fn_tmp;
    ImageXmipp   img, sumimg;
    int          n, nlog, nimg, num, l, fac, istart, m, weight, sumweight, ln;
    SelFile      SFlog, SFo, SFa;
    DocFile      DFi, DFo, DFa, DFtmp;
    DocLine      DLi;
    Matrix1D<double> indataline(3), outdataline(3);
    vector<DocFile> allDFs;
    double       noise, offset, LL, R;
    int          argcp;
    char         **argvp, *copyp;
    string       line;
    bool         out_assign = false;

    // Read arguments --------------------------------------------------------
    try
    {
        fn_logsel = getParameter(argc, argv, "-i", NULL, 1, "MLalign2D_combine: Input logselfile not found");
        SFlog.read(fn_logsel);
        fn_oroot   = getParameter(argc, argv, "-o", "");
        if (fn_oroot == "") fn_oroot = fn_logsel.without_extension() + ".log";
    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        Usage();
        exit(1);
    }

    try
    {
        allDFs.clear();
        SFa.clear();
        DFa.clear();
        SFlog.go_beginning();
        nlog = 0;
        // Read in all logfiles in the logselfile
        while (!SFlog.eof())
        {
            fn_tmp = SFlog.NextImg();
            DFi.read(fn_tmp);
            if (nlog > 0)
                if (DFi.dataLineNo() != nimg) REPORT_ERROR(1, "MLalign2D-log files are not of equal length!");
            DFi.go_beginning();
            allDFs.push_back(DFi);
            nimg = DFi.dataLineNo();
            nlog++;
            // Check whether document files with asigned angles exist and append
            fn_tmp = fn_tmp.remove_extension("log");
            fn_tmp = fn_tmp.add_extension("doc");
            ln = DFa.dataLineNo();
            if (exists(fn_tmp))
            {
                out_assign = true;
                DFa.append(fn_tmp);
                DFa.search(ln);
                DFa.next();
                // Remove header from next file
                DFa.remove_current();
                DFa.search(ln);
                DFa.next();
                DFa.remove_current();
            }
        }
        if (out_assign)
        {
            fn_tmp = fn_oroot + ".doc";
            DFa.write(fn_tmp);
            // Check whether selfiles with reference assignment exist and append
            for (int ref = 1;ref <= nimg;ref++)
            {
                SFa.clear();
                SFlog.go_beginning();
                while (!SFlog.eof())
                {
                    fn_tmp = SFlog.NextImg();
                    fn_tmp = fn_tmp.remove_extension("log");
                    fn_tmp += "_ref";
                    fn_tmp.compose(fn_tmp, ref, "sel");
                    if (exists(fn_tmp)) SFa.append(fn_tmp);
                }
                SFa.sort_by_filenames();
                fn_tmp = fn_oroot + "_ref";
                fn_tmp.compose(fn_tmp, ref, "sel");
                SFa.write(fn_tmp);
            }
        }

        // Calculate weighted average of all images in the reference selfiles
        SFo.reserve(nimg);
        l = 0;
        DFo = allDFs[0];
        while (l < nimg)
        {
            sumweight = 0;
            DFo.go_beginning();
            m = 4 + l * 2;// skip headers + go to relevant line
            for (int mm = 0;mm < m;mm++) DFo.next();
            if ((DFo.get_current_line()).Is_comment())
            {
                fn_img = ((DFo.get_current_line()).get_text()).erase(0, 3);
                sumimg.read(fn_img);
                sumimg().initZeros();
                sumimg.weight() = 0.;
            }
            else REPORT_ERROR(1, "MLalign2D-log files does not have expected format");
            fn_img = fn_oroot + "_ref";
            fn_img.compose(fn_img, l + 1, "xmp");
            DFo.remove_current();
            DFo.previous();
            DFo.insert_comment(fn_img);
            outdataline.initZeros();
            for (n = 0; n < nlog; n++)
            {
                allDFs[n].go_beginning();
                if (allDFs[n].get_current_line().Is_comment()) line = ((allDFs[n].get_current_line()).get_text()).erase(0, 3);
                else REPORT_ERROR(1, "MLalign2D-log files does not have expected format");
                copyp = NULL;
                argvp = NULL;
                generateCommandLine(line, argcp, argvp, copyp);
                weight = AtoI(getParameter(argcp, argvp, "images=", ""));
                sumweight += weight;
                allDFs[n].go_beginning();
                for (int mm = 0;mm < m;mm++) allDFs[n].next();
                if (allDFs[n].get_current_line().Is_comment()) fn_img2 = ((allDFs[n].get_current_line()).get_text()).erase(0, 3);
                else REPORT_ERROR(1, "MLalign2D-log files does not have expected format");
                img.read(fn_img2);
                if (img.weight() < 0.) REPORT_ERROR(1, "MLwp<0: Weights are not set correctly!");
                sumimg() += img() * img.weight();
                sumimg.weight() += img.weight();
                allDFs[n].adjust_to_data_line();
                DLi = allDFs[n].get_current_line();
                outdataline(0) += DLi[0] * weight;
                outdataline(1) += DLi[1] * weight;
                if (DLi[2] > outdataline(2)) outdataline(2) = DLi[2];
            }
            outdataline(0) /= sumweight;
            outdataline(1) /= sumweight;
            DFo.adjust_to_data_line();
            DFo.remove_current();
            DFo.previous();
            DFo.insert_data_line(outdataline);
            if (sumimg.weight() > 0) sumimg() /= sumimg.weight();
            else sumimg().initZeros();
            SFo.insert(fn_img);
            sumimg.write(fn_img);
            l++;
        }
        fn_img = fn_oroot + ".sel";
        SFo.write(fn_img);


        // Calculate weighted average of the parameters in the MLalign2D-log file header
        noise = 0.;
        offset = 0.;
        LL = 0.;
        R = 0.;
        for (n = 0; n < nlog; n++)
        {

            allDFs[n].go_beginning();
            if (allDFs[n].get_current_line().Is_comment()) line = ((allDFs[n].get_current_line()).get_text()).erase(0, 3);
            else REPORT_ERROR(1, "MLalign2D-log files does not have expected format");
            copyp = NULL;
            argvp = NULL;
            generateCommandLine(line, argcp, argvp, copyp);
            LL += AtoF(getParameter(argcp, argvp, "LL=", ""));
            weight = AtoI(getParameter(argcp, argvp, "images=", ""));
            R += weight * AtoF(getParameter(argcp, argvp, "R=", "0"));
            allDFs[n].next();
            if (allDFs[n].get_current_line().Is_comment()) line = ((allDFs[n].get_current_line()).get_text()).erase(0, 3);
            else REPORT_ERROR(1, "MLalign2D-log files does not have expected format");
            generateCommandLine(line, argcp, argvp, copyp);
            noise += weight * AtoF(getParameter(argcp, argvp, "-noise")) * AtoF(getParameter(argcp, argvp, "-noise"));
            offset += weight * AtoF(getParameter(argcp, argvp, "-offset")) * AtoF(getParameter(argcp, argvp, "-offset"));
            istart = AtoI(getParameter(argcp, argvp, "-istart"));
        }
        noise /= sumweight;
        noise = sqrt(noise);
        offset /= sumweight;
        offset = sqrt(offset);
        R /= sumweight;
        DFo.go_beginning();
        DFo.remove_current();
        DFo.previous();
        line = "MLalign2D-logfile: Number of images= " + floatToString(sumweight);
        if (LL != 0.) line += " LL= " + floatToString(LL, 12, 5);
        if (R != 0.) line += " R= " + floatToString(R, 12, 5);
        DFo.insert_comment(line);
        DFo.remove_current();
        DFo.previous();
        line = "-noise " + floatToString(noise, 10, 7) + " -offset " + floatToString(offset, 10, 7) + " -istart " + ItoA(istart);
        DFo.insert_comment(line);
        DFo.write(fn_oroot + ".log");

    }
    catch (Xmipp_error XE)
    {
        cout << XE;
    }
    exit(0);

} //main

/* Usage ------------------------------------------------------------------- */
void Usage()
{
    cerr << "Purpose:\n";
    cerr << "    Combine weighted average images from parallel runs of xmipp_MLalign2D \n";
    cerr << "Usage: MLalign2D_combine <parameters>\n"
    << "   -i <logselfile>  : selection file containing MLalign2D-log files of parallel runs\n"
    << "   -o <rootname>    : rootname for output files (default: logsel-rootname)\n"

    << endl;
}


