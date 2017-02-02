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

#include <data/xmipp_program.h>
#include <data/xmipp_funcs.h>


class ProgTestWork: public XmippProgram
{
protected:
    int workTime, randMin, randMax;

    double param1, param2;
    double df, limit0, limitF;
    bool   do_limit0, do_limitF;
    std::string noise_type;

    void defineParams()
    {
        addUsageLine("Testing paralellization with some sleeping program.");
        addParamsLine("  [--time <sec=5> ]              : Number of seconds to sleep ");
        addParamsLine("  [--rand <min=0> <max=2> ]       : Randon variation (uniform distributed between min and max");
        addParamsLine("  [--tag <string> ]               : Tag to monitor. ");

    }

    void readParams()
    {
    }

    void run()
    {
        size_t ms = getIntParam("--time") * 1000; // time in mili-seconds

        Timer t;
        t.tic();

        while (true)
        {
            // Just waste some cpu cycles
            for (int i = 0; i < 100000000; ++i);

            if (t.elapsed() >= ms)
                break;
        }
    }

}
;//end of class TestWork

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgTestWork program;
    program.read(argc, argv);
    //init random seed
    randomize_random_generator();
    program.tryRun();
    return 0;
}

