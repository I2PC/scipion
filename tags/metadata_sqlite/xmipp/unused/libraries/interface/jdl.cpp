/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Roberto Marabini (roberto@mipg.upenn.edu)
 *
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
 *                                      <
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "jdl.h"

#include <data/selfile.h>

#include <fstream>
#include <stdlib.h>

/* Clear ------------------------------------------------------------------- */
void JDLFile::clear()
{
    jdl_type = "Job";
    jdl_jobtype = "normal";
    jdl_virtualorganization = "biomed";
    job_root = "";
    job_sh.clear();
    input_files.clear();
    input_selfiles.clear();
    input_wildfiles.clear();
    output_files.clear();
    output_selfiles.clear();
    output_wildfiles.clear();
    running_programs.clear();
}

/* Write files ------------------------------------------------------------- */
void JDLFile::write()
{
    // Check integrity
    if (publisher_data_out == "" && publisher_data_in != "")
        publisher_data_out = publisher_data_in;
    if (local_output_dir == "")
        REPORT_ERROR(1, "JDLFile::write: doesn't have a local output directory");

    // Get a random filename for the data tar
    FileName fn_tmpdir;
    fn_tmpdir.init_random(4);
    FileName fn_data_in;
    fn_data_in.init_random(4);
    fn_data_in = (std::string)"xmipp" + job_root + "in" + fn_data_in;
    FileName fn_data_out;
    fn_data_out.init_random(4);
    fn_data_out = (std::string)"xmipp" + job_root + "out" + fn_data_out;

    // Write the local script ...............................................
    std::ofstream fh_sh;
    fh_sh.open((job_root + "_local_in.sh").c_str());
    if (!fh_sh)
        REPORT_ERROR(1, (std::string)"JDLFile::write: Cannot open file " +
                     job_root + "_local_in.sh for output");

    // Create temporary directory with all the data
    fh_sh << "#!/bin/sh\n";
    fh_sh << "# Create temporary directory\n";
    fh_sh << "mkdir " << fn_tmpdir << std::endl;
    if (local_output_dir != "." && local_output_dir != "..")
        fh_sh << "mkdir " << fn_tmpdir << "/" << local_output_dir << std::endl;
    fh_sh << "cd " << fn_tmpdir << std::endl;
    fh_sh << std::endl;

    // Search for all programs
    fh_sh << "# Create links to all programs\n";
    int imax = running_programs.size();
    for (int i = 0; i < imax; i++)
        fh_sh << "cp " << " `which " << running_programs[i] << "` .\n";

    // Add all input files
    fh_sh << "# Create links to all input files\n";
    imax = input_files.size();
    for (int i = 0; i < imax; i++)
        fh_sh << "cp " << input_files[i] << " .\n";
    fh_sh << std::endl;

    // Add all files in selfiles
    fh_sh << "# Create links to all input files in selfiles\n";
    imax = input_selfiles.size();
    for (int i = 0; i < imax; i++)
    {
        SelFile SF;
        SF.read(input_selfiles[i]);
        while (!SF.eof())
        {
            FileName fn_img = SF.NextImg();
            if (fn_img=="") break;
            fh_sh << "cp " << fn_img << " .\n";
        }
    }
    fh_sh << std::endl;

    // Build the input tar
    fh_sh << "# Create the tar file and remove the temporary directory\n";
    fh_sh << "tar zchf " << fn_data_in << ".tgz .\n";
    fh_sh << "mv " << fn_data_in << ".tgz ..\n";
    fh_sh << "cd ..\n";
    fh_sh << "rm -rf " << fn_tmpdir << std::endl;
    imax = input_wildfiles.size();
    for (int i = 0; i < imax; i++)
        fh_sh << "tar zrf " << fn_data_in << ".tgz " << input_wildfiles[i] << std::endl;
    fh_sh << std::endl;

    // Publish the data
    fh_sh << "# Publish the data\n";
    if (publisher_data_in != "")
    {
        fh_sh << "lcg-cr --vo " << jdl_virtualorganization
        << " -l lfn:" << fn_data_in << ".tgz"
        << " -d " << publisher_data_in
        << " file:`pwd`/" << fn_data_in << ".tgz\n"
        ;
    }
    else
    {
        // Rank the Computing Elements according to free CPUs ................
        // Create the jdl to get the rank of free CPUs
        std::string rankCE_jdl = (std::string)"xmipp" + job_root + "rankCE.jdl";
        std::string rankCE_list = (std::string)"xmipp" + job_root + "rankCE_list_match.txt";
        std::string rankSE_list = (std::string)"xmipp" + job_root + "rankSE_list_match.txt";
        std::ofstream fh_rank_jdl;
        fh_rank_jdl.open(rankCE_jdl.c_str());
        if (!fh_rank_jdl)
            REPORT_ERROR(1, "JDLFile::write: Cannot open file " +
                         rankCE_jdl + " for output");
        fh_rank_jdl
        << "Type = \"" << jdl_type << "\";\n"
        << "JobType = \"" << jdl_jobtype << "\";\n"
        << "VirtualOrganisation = \"" << jdl_virtualorganization << "\";\n"
        << "Executable = \"" << job_root << "remote.sh\";\n"
        << "StdOutput = \"" << job_root << "stdout.txt\";\n"
        << "StdError = \"" << job_root << "stderr.txt\";\n"
        << "InputSandbox = {\"" << job_root << "remote.sh\"};\n"
        << "OutputSandbox = {\"" << job_root << "stdout.txt\","
        << "\"" << job_root << "stderr.txt\"};\n"
        << "Rank = other.GlueCEStateFreeCPUs;\n"
        << "Requirements =(other.GlueCEPolicyMaxCPUTime >= 1440);\n"
        ;
        fh_rank_jdl.close();

        // Get the rank and select storage element
        fh_sh
        << "# Get the best storage element and publish the data\n"
        << "edg-job-list-match --rank " << rankCE_jdl << " > "
        << rankCE_list << std::endl
        << "lcg-infosites --vo " << jdl_virtualorganization << " closeSE > "
        << rankSE_list << std::endl
        << "i=11\n"
        << "until ! [ \"$se\" == \"\" ];\n"
        << "do\n"
        << "   ce=`cat " << rankCE_list << " |head -$i |tail -1 |awk -F\"/\" '{print $1}'`\n"
        << "   line=`grep -n $ce " << rankSE_list << " |head -1| awk -F\":\" '{print $1}'`\n"
        << "   if ! [ \"$line\" = \"\"  ];\n"
        << "   then\n"
        << "      let line=$line+1\n"
        << "      se=`head -$line " << rankSE_list << "|tail -1 | awk '{print $6}'`\n"
        << "      lcg-cr --vo " << jdl_virtualorganization
        << " -l lfn:" << fn_data_in << ".tgz -d $se file:`pwd`/" << fn_data_in << ".tgz\n"
        << "   if [ `lcg-lr --vo " << jdl_virtualorganization << " lfn:" << fn_data_in << ".tgz` ];\n"
        << "   then\n"
        << "      break\n"
        << "   else\n"
        << "      se=\"\" \n"
        << "   fi\n"
        << "   fi\n"
        << "   let i=$i+2\n"
        << "done\n"
        << "rm " << rankSE_list << std::endl
        << "rm " << rankCE_list << std::endl
        << "rm " << rankCE_jdl  << std::endl
        << std::endl
        << "   if [ $se == \"\" ];\n"
        << "   then\n"
        << "     echo \" Fail: lcg_lr: \"No such file or directory\"\"  > error_publish.txt \n"
        << "   fi\n"
        ;
    }
    fh_sh.close();
    system(((std::string)"chmod u+x " + job_root + "_local_in.sh").c_str());

    // Write the remote script ..............................................
    fh_sh.open((job_root + "remote.sh").c_str());
    if (!fh_sh)
        REPORT_ERROR(1, (std::string)"JDLFile::write: Cannot open file " +
                     job_root + "remote.sh for output");

    // Get the input data and programs
    fh_sh << "#!/bin/sh\n\n";
    fh_sh << "export PATH=`pwd`:$PATH\n\n";
    fh_sh << "# Get the input data and programs\n";
    fh_sh << "lcg-cp --vo " << jdl_virtualorganization
    << " lfn:" << fn_data_in << ".tgz"
    << " file:`pwd`/" << fn_data_in << ".tgz\n";
    fh_sh << "tar -zxf " << fn_data_in << ".tgz\n";
    fh_sh << "rm " << fn_data_in << ".tgz\n";
    imax = running_programs.size();
    for (int i = 0; i < imax; i++)
        fh_sh << "chmod 700 " << running_programs[i] << "\n";
    fh_sh << std::endl;

    // Write user commands
    fh_sh << "# User commands\n";
    imax = job_sh.size();
    for (int i = 0; i < imax; i++)
        fh_sh << job_sh[i] << std::endl;
    fh_sh << std::endl;

    // Pack the output
    fh_sh << "# Prepare an output directory\n";
    fh_sh << "mkdir " << fn_tmpdir << std::endl << std::endl;
    fh_sh << "cd " << fn_tmpdir << std::endl;

    // Add all output files
    fh_sh << "# Create links to all output files\n";
    imax = output_files.size();
    for (int i = 0; i < imax; i++)
        fh_sh << "cp " << output_files[i] << " .\n";
    fh_sh << std::endl;

    // Add all files in selfiles
    fh_sh << "# Create links to all output files in selfiles\n";
    imax = output_selfiles.size();
    for (int i = 0; i < imax; i++)
    {
        SelFile SF;
        SF.read(output_selfiles[i]);
        while (!SF.eof())
        {
            FileName fn_img = SF.NextImg();
            if (fn_img=="") break;
            fh_sh << "cp " << fn_img << " .\n";
        }
    }
    fh_sh << std::endl;

    // Build the output tar
    fh_sh << "# Create the tar file and remove the temporary directory\n";
    fh_sh << "tar cf " << fn_data_out << ".tar .\n";
    fh_sh << "mv " << fn_data_out << ".tar ..\n";
    fh_sh << "cd ..\n";
    fh_sh << "rm -rf " << fn_tmpdir << std::endl;
    imax = output_wildfiles.size();
    for (int i = 0; i < imax; i++)
        fh_sh << "tar rf " << fn_data_out << ".tar " << output_wildfiles[i] << std::endl;
    fh_sh << "gzip " << fn_data_out << ".tar\n";
    fh_sh << std::endl;
    fh_sh << std::endl;

    // Publish the data
    fh_sh << "# Publish the data\n"
    ;
    if (publisher_data_out != "")
    {
        fh_sh << "lcg-cr --vo " << jdl_virtualorganization
        << " -l lfn:" << fn_data_out << ".tar.gz"
        << " -d " << publisher_data_out
        << " file:`pwd`/" << fn_data_out << ".tar.gz\n"
        ;
    }
    else
    {
        fh_sh << "i=1\n"
        << "until ! [ \"$se\" == \"\" ];\n"
        << "do\n"
        << "se=`lcg-infosites --vo " << jdl_virtualorganization
        << " se | grep --regexp=\"[i,f][t,r]\" | head -$i | tail -1 | awk '{print $4}'`\n"
        << " lcg-cr --vo " << jdl_virtualorganization
        << " -l lfn:" << fn_data_out << ".tar.gz -d $se file:`pwd`/" << fn_data_out << ".tar.gz\n"
        << " if [ `lcg-lr --vo " << jdl_virtualorganization << " lfn:" << fn_data_in << ".tar.gz` ];\n"
        << "   then\n"
        << "      break\n"
        << "   else\n"
        << "      se=\"\" \n"
        << "   fi\n"
        << "   let i=$i+1\n"
        << "done\n"
        ;
    }
    fh_sh.close();
    system(((std::string)"chmod u+x " + job_root + "remote.sh").c_str());

    // Write JDL file .......................................................
    std::ofstream fh_jdl;
    fh_jdl.open((job_root + ".jdl").c_str());
    if (!fh_jdl)
        REPORT_ERROR(1, (std::string)"JDLFile::write: Cannot open file " +
                     job_root + ".jdl for output");
    fh_jdl << "Type = \"" << jdl_type << "\";\n"
    << "JobType = \"" << jdl_jobtype << "\";\n"
    << "VirtualOrganisation = \"" << jdl_virtualorganization << "\";\n"
    << "Executable = \"" << job_root << "remote.sh\";\n"
    << "StdOutput = \"" << job_root << "stdout.txt\";\n"
    << "StdError = \"" << job_root << "stderr.txt\";\n"
    << "InputSandbox = {\"" << job_root << "remote.sh\"};\n"
    << "OutputSandbox = {\"" << job_root << "stdout.txt\","
    << "\"" << job_root << "stderr.txt\","
    << "\"xmipp" << job_root << "_error_publish_out.txt\"};\n"
    << "InputData = {\"lfn:" << fn_data_in << ".tgz\"};\n"
    << "DataAccessProtocol= {\"gsiftp\"};\n"
    << "Requirements =(other.GlueCEPolicyMaxCPUTime >= 1440);\n"
    ;
    fh_jdl.close();

    // Write the local out script ...........................................
    fh_sh.open((job_root + "_local_out.sh").c_str());
    if (!fh_sh)
        REPORT_ERROR(1, (std::string)"JDLFile::write: Cannot open file " +
                     job_root + "_local_out.sh for output");

    // Get the output data
    fh_sh << "#!/bin/csh\n";
    fh_sh << "# Get the output data and programs\n";
    fh_sh << "lcg-cp --vo " << jdl_virtualorganization
    << " lfn:" << fn_data_out << ".tar.gz"
    << " file:`pwd`/" << fn_data_out << ".tgz\n";
    fh_sh << "lcg-del -a --vo " << jdl_virtualorganization
    << " lfn:" << fn_data_in << ".tgz\n";
    fh_sh << "lcg-del -a --vo " << jdl_virtualorganization
    << " lfn:" << fn_data_out << ".tar.gz\n";
    if (local_output_dir == ".")
    {
        fh_sh << "tar -zxf " << fn_data_out << ".tar.gz\n";
    }
    else if (local_output_dir == "..")
    {
        fh_sh << "export CURRENT_DIR=`pwd`\n";
        fh_sh << "cd ..\n";
        fh_sh << "tar -zxf $CURRENT_DIR/" << fn_data_out << ".tgz\n";
        fh_sh << "cd $CURRENT_DIR\n";
    }
    else
    {
        fh_sh << "if ( !(-e " << local_output_dir << ") ) then\n"
        << "   mkdir " << local_output_dir << "\n"
        << "endif\n";
        fh_sh << "cd " << local_output_dir << std::endl;
        fh_sh << "tar -zxf ../" << fn_data_out << ".tgz\n";
        fh_sh << "cd ..\n";
    }
    fh_sh << "rm " << fn_data_out << ".tgz\n";
    fh_sh << "rm " << fn_data_in << ".tgz\n";

    fh_sh.close();
    system(((std::string)"chmod u+x " + job_root + "_local_out.sh").c_str());
}
