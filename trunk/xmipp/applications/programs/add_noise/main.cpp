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

#include <data/progs.h>
#include <data/args.h>

class Add_noise_parameters: public Prog_parameters
{
public:
    double noise_min, noise_max;
    double noise_avg, noise_stddev;
    double df, limit0, limitF;
    bool   gaussian,  uniform, student, do_limit0, do_limitF;

    void read(int argc, char **argv)
    {
        Prog_parameters::read(argc, argv);
        gaussian = uniform = student = false;
        do_limit0 = checkParameter(argc, argv, "-limit0");
        if (do_limit0)
            limit0 =  textToFloat(getParameter(argc, argv, "-limit0"));
        do_limitF = checkParameter(argc, argv, "-limitF");
        if (do_limitF)
            limitF =  textToFloat(getParameter(argc, argv, "-limitF"));

        if (checkParameter(argc, argv, "-gaussian"))
        {
            gaussian = true;
            int i = paremeterPosition(argc, argv, "-gaussian");
            if (i + 1 >= argc) REPORT_ERROR(1, "Not enough parameters after -gaussian");
            noise_stddev = textToFloat(argv[i+1]);
            if (i + 2 < argc)
            {
                noise_avg = textToFloat(argv[i+2]);
            }
            else noise_avg = 0;
        }
        else if (checkParameter(argc, argv, "-student"))
        {
            student = true;
            int i = paremeterPosition(argc, argv, "-student");
            if (i + 2 >= argc) REPORT_ERROR(1, "Not enough parameters after -student");
            df = textToFloat(argv[i+1]);
            noise_stddev = textToFloat(argv[i+2]);
            if (i + 3 < argc)
            {
                noise_avg = textToFloat(argv[i+3]);
            }
            else noise_avg = 0;
        }
        else if (checkParameter(argc, argv, "-uniform"))
        {
            uniform = true;
            int i = paremeterPosition(argc, argv, "-uniform");
            if (i + 2 >= argc) REPORT_ERROR(1, "Not enough parameters after -uniform");
            noise_min = textToFloat(argv[i+1]);
            noise_max = textToFloat(argv[i+2]);
        }
        else
            REPORT_ERROR(1, "Unknown noise type");
    }

    void show()
    {
        Prog_parameters::show();
        if (gaussian)
            std::cout << "Noise avg=" << noise_avg << std::endl
            << "Noise stddev=" << noise_stddev << std::endl;
        else if (student)
            std::cout << "Degrees of freedom= "<<df<< std::endl
                      << "Noise avg=" << noise_avg << std::endl
                      << "Noise stddev=" << noise_stddev << std::endl;
        else if (uniform)
            std::cout << "Noise min=" << noise_min << std::endl
                      << "Noise max=" << noise_max << std::endl;
        if (do_limit0)
            std::cout << "Crop noise histogram below=" << limit0 << std::endl;
        if (do_limitF)
            std::cout << "Crop noise histogram above=" << limitF << std::endl;
    }

    void usage()
    {
        Prog_parameters::usage();
        std::cerr 
            << "  [-gaussian <stddev> [<avg>=0]] : Gaussian noise parameters\n"
            << "  [-student <df> <stddev> [<avg>=0]] : t-student noise parameters\n"
            << "  [-uniform  <min> <max>]   : Uniform noise parameters\n"
            << "  [-limit0 <float> ]        : Crop noise histogram below this value \n"
            << "  [-limitF <float> ]        : Crop noise histogram above this value \n";
    }
};

bool process_img(ImageXmipp &img, const Prog_parameters *prm)
{
    Add_noise_parameters *eprm = (Add_noise_parameters *) prm;
    if (eprm->gaussian)
        img().addNoise(eprm->noise_avg, eprm->noise_stddev, "gaussian");
    else if (eprm->student)
        img().addNoise(eprm->noise_avg, eprm->noise_stddev, "student", eprm->df);
    else if (eprm->uniform)
        img().addNoise(eprm->noise_min, eprm->noise_max, "uniform");
    if (eprm->do_limit0)
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(img())
        {
            dMij(img(),i,j) = XMIPP_MAX(dMij(img(),i,j),eprm->limit0);
        }
    if (eprm->do_limitF)
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX2D(img())
        {
            dMij(img(),i,j) = XMIPP_MIN(dMij(img(),i,j),eprm->limitF);
        }

    return true;
}

bool process_vol(VolumeXmipp &vol, const Prog_parameters *prm)
{
    Add_noise_parameters *eprm = (Add_noise_parameters *) prm;
    if (eprm->gaussian)
        vol().addNoise(eprm->noise_avg, eprm->noise_stddev, "gaussian");
    else if (eprm->student)
        vol().addNoise(eprm->noise_avg, eprm->noise_stddev, "student", eprm->df);
    else if (eprm->uniform)
        vol().addNoise(eprm->noise_min, eprm->noise_max, "uniform");
    if (eprm->do_limit0)
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(vol())
        {
            dVkij(vol(),k,i,j) = XMIPP_MAX(dVkij(vol(),k,i,j),eprm->limit0);
        }
    if (eprm->do_limitF)
        FOR_ALL_DIRECT_ELEMENTS_IN_MATRIX3D(vol())
        {
            dVkij(vol(),k,i,j) = XMIPP_MIN(dVkij(vol(),k,i,j),eprm->limitF);
        }
    return true;
}

int main(int argc, char **argv)
{
    Add_noise_parameters prm;
    randomize_random_generator();
    SF_main(argc, argv, &prm, (void*)&process_img, (void*)&process_vol);
}

/* Menus ------------------------------------------------------------------- */
/*Colimate:
   PROGRAM Add_noise {
      url="http://www.cnb.uam.es/~bioinfo/NewXmipp/Applications/Src/Add_noise/Help/add_noise.html";
      help="Add noise to volumes and images";
      OPEN MENU menu_add_noise;
      COMMAND LINES {
 + usual: xmipp_add_noise
               #include "prog_line.mnu"
               [-gaussian $STDDEV [$AVG]]
               [-uniform  $MIN $MAX]
      }
      PARAMETER DEFINITIONS {
        #include "prog_vars.mnu"
        OPT(-gaussian) {label="Add gaussian noise";}
           $STDDEV {label="Standard deviation"; type=float;}
           $AVG    {label="Average";            type=float;}
        OPT(-uniform) {label="Add uniform noise";}
           $MIN    {label="Minimum value";      type=float;}
           $MAX    {label="Maximum value";      type=float;}
      }
   }

   MENU menu_add_noise {
      #include "prog_menu.mnu"
      "Noise parameters"
      OPT(-gaussian)
      OPT(-uniform)
   }
*/
