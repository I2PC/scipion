/***************************************************************************
 *
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

#include <data/progs.h>
#include <data/args.h>
#include <data/strings.h>
#include <iomanip>
#include <iostream>
#include <queue>

///Constant values for weighting the keywords match
#define NAME_MATCH 5
#define KEYS_MATCH 3
#define DESC_MATCH 1

class ProgApropos: public XmippProgram
{
protected:
    StringVector keywords;
    bool update;
    bool list;
    int maxlen;

    void defineParams()
    {
        //usage
        addUsageLine("Search for Xmipp programs that are related to some keywords.");
        addUsageLine("Useful for when does not remeber a program name. gaussian noise header");
        //examples
        addExampleLine("Search for program containing the keyword 'header'", false);
        addExampleLine("   xmipp_apropos -k header");
        addExampleLine("Search for program containing the keywords 'noise' and 'gaussian'", false);
        addExampleLine("   xmipp_apropos -k \"noise gaussian\"");
        addExampleLine("List all xmipp programs", false);
        addExampleLine("   xmipp_apropos --list");
        //params
        addParamsLine(" -k <key>   :keyword to search programs matching");
        addParamsLine("                   :if you want to search for more than one keyword,");
        addParamsLine("                   :use quotes. See example above");
        addParamsLine("   alias --keyword;");
        addParamsLine("or -u        : Update the database with programs info");
        addParamsLine("   alias --update;");
        addParamsLine("or -l          : List all Xmipp programs");
        addParamsLine("   alias --list;");
    }

    void readParams()
    {
        //Get keywords from cmd line
        if (checkParam("--keyword"))
        {
            String keys = getParam("--keyword");
            toLower(keys);
            //StringVector vector;
            tokenize(keys, keywords);
        }
        update = checkParam("--update");
        list = checkParam("--list");
    }

    void show()
    {}

    int getRank(const String &line, const String &key, int weight) const
    {
        String lineLower = line;
        toLower(lineLower);
        if (lineLower.find(key) != String::npos)
            return weight;
        return 0;
    }

    void printProgram(DbProgram * prog)
    {
        size_t endline;
        endline = prog->description.find_first_of('\n');
        std::cout
        << std::left << std::setw(maxlen) << prog->name
        << "- " << prog->description.substr(0, endline) << std::endl;
    }

public:
    void run()
    {
        std::vector<DbProgram*> progs;
        //std::vector<DbProgram*> progsRank;

        std::vector<DbProgram*>::iterator it, it2;
        StringVector::iterator keyIt;

        FileName dbName = xmippBaseDir().append("/programs.db");

        if (update)
        {
            createDB();
            exit(0);
        }

        if (!exists(dbName))
            createDB();
        XmippDB db(dbName);
        db.beginTrans();
        db.selectPrograms(progs);
        db.commitTrans();

        String line;
        DbProgram *prog;
        int len;
        maxlen = 20; //take max program name

        for (int i  = 0; i < progs.size(); ++i)
        {
            prog = progs[i];
            maxlen = XMIPP_MAX(prog->name.length(), maxlen);
            if (list)
                prog->rank = 1;
            else
            {
                prog->rank = 0;
                for (keyIt = keywords.begin(); keyIt < keywords.end(); ++keyIt)
                {
                    prog->rank += getRank(prog->name, *keyIt, NAME_MATCH);
                    prog->rank += getRank(prog->keywords, *keyIt, KEYS_MATCH);
                    prog->rank += getRank(prog->description, *keyIt, DESC_MATCH);
                }
                //Order by insertion sort
                for (int j = i - 1; j >= 0 && prog->rank > progs[j]->rank; --j)
                {
                    //Swap values
                    progs[j + 1] = progs[j];
                    progs[j] = prog;
                }
            }
        }

        //Print out results
        maxlen += 3;
        for (int i  = 0; i < progs.size(); ++i)
        {
            prog = progs[i];
            if (prog->rank)
                printProgram(prog);
        }
    }

    void createDB()
    {
        ///Create db
        FileName dirPath = xmippBaseDir().append("/bin/");
        FileName dbName = xmippBaseDir().append("/programs.db");
        XmippDB db(dbName);
        db.beginTrans();
        db.createCategoryTable();
        db.createProgramTable();
        db.commitTrans();
        std::vector<FileName> files;

        String progCmd, mpiPrefix;
        FILE * input;
        getdir(dirPath, files);
        char readbuffer[256];

        int pfd[2];
        if (pipe(pfd) == -1)
        {
            perror("pipe");
            exit(EXIT_FAILURE);
        }
        //Redirect std::cerr writing end of the pipe
        //dup2(pfd[1], 2);
        std::vector<FileName>::const_iterator iter;
        int times = 0;
        for (iter = files.begin(); iter != files.end(); ++iter)
        {
            if (iter->find("xmipp_") == 0 &&
                iter->rfind("j") != iter->length()-1 &&
                iter->find("xmipp_test") == std::string::npos)
            {
                mpiPrefix = iter->find("_mpi_") == std::string::npos  ? "" : "mpirun -np 2 ";
                int fdOut;
                progCmd = mpiPrefix + *iter + " --xmipp_write_definition";
                std::cerr << progCmd << std::endl;
                                dup2(1, fdOut);//save std::cout
                                dup2(pfd[1], 1);
                                input = popen(progCmd.c_str(), "r");
                                int nbytes;
                                bool ok = false;
                                pclose(input);
                                dup2(fdOut, 1);//restore std::cout
            }
        }
    }
}
;//end of class ProgApropos

/* MAIN -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    ProgApropos program;
    program.read(argc, argv);
    program.tryRun();
}

