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
#include <data/image_collection.h>
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

    void defineParams()
    {
        addUsageLine("Search for Xmipp programs that are related to some keywords.");
        addUsageLine("Useful for when does not remeber a program name. gaussian noise header");
        addUsageLine("Examples:");
        addUsageLine("   xmipp_apropos -k header");
        addUsageLine("   xmipp_apropos -k \"noise gaussian\"");

        addParamsLine(" --keyword <key>   :keyword to search programs matching");
        addParamsLine("                   :if you want to search for more than one keyword,");
        addParamsLine("                   :use quotes. See example above");
        addParamsLine("   alias -k;");
        addParamsLine("or --update        : Update the database with programs info");
        addParamsLine("   alias -up;");
        addParamsLine("or --list          : List all Xmipp programs");
        addParamsLine("   alias -l;");

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

public:
    void run()
    {
        std::vector<DbProgram*> progs;
        std::priority_queue<DbProgram*> progsRank;

        std::vector<DbProgram*>::iterator it;
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
        int len, maxlen = 0; //take max program name

        for (it = progs.begin(); it < progs.end(); ++it)
        {
            prog = *it;
            if (list)
                prog->rank = 3;
            else
            {
                prog->rank = 0;
                for (keyIt = keywords.begin(); keyIt < keywords.end(); ++keyIt)
                {
                    prog->rank += getRank(prog->name, *keyIt, NAME_MATCH);
                    prog->rank += getRank(prog->keywords, *keyIt, KEYS_MATCH);
                    prog->rank += getRank(prog->description, *keyIt, DESC_MATCH);
                }
            }
            if (prog->rank > 2)
            {
                maxlen = XMIPP_MAX(prog->name.length(), maxlen);
                progsRank.push(prog);
            }
        }

        //Print out results
        size_t endline;
        maxlen += 3;
        while (!progsRank.empty())
        {
            prog = progsRank.top();
            progsRank.pop();
            endline = prog->description.find_first_of('\n');
            std::cout
            << std::left << std::setw(maxlen) << prog->name
            << "- " << prog->description.substr(0, endline) << std::endl;
            delete prog;
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

        std::string progCmd;
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
                iter->find("_mpi_") == std::string::npos &&
                iter->rfind("j") != iter->length()-1 &&
                iter->find("xmipp_test") == std::string::npos)
            {
                int fdOut;
                progCmd = *iter + " --xmipp_write_definition";
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
    try
    {
        ProgApropos program;
        program.read(argc, argv);
        program.run();
    }
    catch (XmippError xe)
    {
        std::cerr << xe;
    }
}

