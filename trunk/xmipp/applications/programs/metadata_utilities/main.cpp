/*
 * main.cpp
 *
 *  Created on: Sep 28, 2010
 *      Author: roberto
 */
#include <data/argsparser.h>
#include <data/program.h>
class ProgMetadataUtilities: public XmippProgram
{
private:

protected:
    void defineParams()
    {
        addUsageLine ("Perform several operation over the metadata files");
        addParamsLine("   --union  <md1> <md2>             : union of md1 and md2");
        addParamsLine("     alias -u;");
        addParamsLine("or -intersection <md1> <md2>        : Intersection of md1 and md2");
        addParamsLine("     alias -i;");
        addParamsLine("or -subtraction <md1> <md2>         : Subtraction of md1 and md2");
        addParamsLine("   -o  <w1>                   : Cutoff freq (<1/2 or A)");
        addParamsLine("     alias -s;");
    }

    void readParams()
    {
    	std::cerr << "readParams" <<std::endl;
        //int round_shifts = checkParam("-round_shifts");
    }
public:
    void run()
    {
        std::cerr << "running" <<std::endl;
    }
};

int main(int argc, char **argv)
{
    std::cerr << "hello" <<std::endl;
    ProgMetadataUtilities program;
    program.read(argc, argv);
    program.run();
}
