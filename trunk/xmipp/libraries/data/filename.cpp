/***************************************************************************
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "filename.h"
#include "funcs.h"

// Constructor with root, number and extension .............................
void FileName::compose(const std::string &str, int no, const std::string &ext)
{
    *this = (FileName) str;
    if (no != -1)
    {

        char aux_str[FILENAMENUMBERLENGTH+1];
        std::string tmp_fileformat;
        tmp_fileformat = (std::string) "%0" +
                         integerToString(FILENAMENUMBERLENGTH)+
                         (std::string)"d";
        sprintf(aux_str, tmp_fileformat.c_str(), no);
        *this += aux_str;
    }

    if (ext != "")
        *this += (std::string)"." + ext;
}

// Constructor: prefix number and filename, mainly for selfiles..
void FileName::compose(int no , const std::string &str)
{
    *this = (FileName) str;
    if (no != -1)
    {

        char aux_str[FILENAMENUMBERLENGTH+1];
        std::string tmp_fileformat;
        tmp_fileformat = (std::string) "%0" +
                         integerToString(FILENAMENUMBERLENGTH)+
                         (std::string)"d@";
        sprintf(aux_str, tmp_fileformat.c_str(), no);
        *this = aux_str + str;
    }
    else
        *this = str;


}

// Is in stack ............................................................
bool FileName::isInStack() const
{
    return find("@") != std::string::npos;
}

// Decompose ..............................................................
void FileName::decompose(int &no, std::string &str) const
{
    size_t idx = find('@');
    no = textToInteger(substr(0,idx));
    str = substr(idx+1,length()-idx);
}

// Get the root name of a filename .........................................
FileName FileName::getRoot() const
{
    int skip_directories = find_last_of("/") + 1;
    int point = find_first_of(".", skip_directories);
    if (point == -1)
        point = length();
    int root_end = find_last_not_of("0123456789", point - 1);
    if (root_end + 1 != point)
        if (point - root_end > FILENAMENUMBERLENGTH)
            root_end = point - FILENAMENUMBERLENGTH - 1;
    return (FileName) substr(0, root_end + 1);
}

// Convert to lower case characters .........................................
FileName FileName::toLowercase() const
{
    FileName result = *this;
    for(unsigned int i=0;i<result.length();i++)
        result[i] = tolower(result[i]);
    return result;
}

// Convert to upper case characters .........................................
FileName FileName::toUppercase() const
{
    FileName result = *this;
    for(unsigned int i=0;i<result.length();i++)
        result[i] = toupper(result[i]);
    return result;
}

// Is substring present?
bool FileName::contains(const std::string& str) const
{
    int point = rfind(str);
    if (point > -1)
        return true;
    else
        return false;
}

// Get substring before first instance of str
FileName FileName::beforeFirstOf(const std::string& str) const
{
    int point = find_first_of(str);
    if (point > -1)
        return substr(0, point);
    else
        return *this;
}

// Get substring before last instance of str
FileName FileName::beforeLastOf(const std::string& str) const
{
    int point = find_last_of(str);
    if (point > -1)
        return substr(0, point);
    else
        return *this;
}

// Get substring after first instance of str
FileName FileName::afterFirstOf(const std::string& str) const
{
    int point = find_first_of(str);
    if (point > -1)
        return substr(point + 1);
    else
        return *this;
}

// Get substring after last instance of str
FileName FileName::afterLastOf(const std::string& str) const
{
    int point = find_last_of(str);
    if (point > -1)
        return substr(point + 1);
    else
        return *this;
}

// Get the base name of a filename .........................................
std::string FileName::getBaseName() const
{
    std::string basename = "";
    std::string myname = *this;
    int myindex = 0;
    for (int p = myname.size() - 1; p >= 0; p--)
    {
        if (myname[p] == '/')
        {
            myindex = p + 1;
            break;
        }
    }
    for (int p = myindex; p < myname.size(); p++)
    {
        if (myname[p] != '.')
            basename += myname[p];
        else
            break;
    }
    return basename;
}

// Get number from file ....................................................
int FileName::getNumber() const
{
    int skip_directories = find_last_of("/") + 1;
    int point = find_first_of(".", skip_directories);
    if (point == -1)
        point = length();
    int root_end = find_last_not_of("0123456789", point - 1);
    if (root_end + 1 != point)
    {
        if (point - root_end > FILENAMENUMBERLENGTH)
            root_end = point - FILENAMENUMBERLENGTH - 1;
        std::string aux = substr(root_end + 1, point - root_end + 1);
        return atoi(aux.c_str());
    }
    else
        return -1;
}

// Get the extension of a filename .........................................
std::string FileName::getExtension() const
{
    int skip_directories = find_last_of("/") + 1;
    int first_point = find_first_of(".", skip_directories);
    if (first_point == -1)
        return "";
    else
        return substr(first_point + 1);
}

// Init random .............................................................
void FileName::initRandom(int length)
{
    randomize_random_generator();
    *this = "";
    for (int i = 0; i < length; i++)
        *this += 'a' + FLOOR(rnd_unif(0, 26));
}

// Add at beginning ........................................................
FileName FileName::addPrefix(const std::string &prefix) const
{
    FileName retval = *this;
    int skip_directories = find_last_of("/") + 1;
    return retval.insert(skip_directories, prefix);
}

// Add at the end ..........................................................
FileName FileName::addExtension(const std::string &ext) const
{
    if (ext == "")
        return *this;
    else
    {
        FileName retval = *this;
        retval = retval.append((std::string)"." + ext);
        return retval;
    }
}

// Remove last extension ...................................................
FileName FileName::withoutExtension() const
{
    FileName retval = *this;
    return retval.substr(0, rfind("."));
}

// Remove root .............................................................
FileName FileName::withoutRoot() const
{
    return without(getRoot());
}

// Insert before extension .................................................
FileName FileName::insertBeforeExtension(const std::string &str) const
{
    int point = -1;
    bool done = false;
    do
    {
        point = find(".", point + 1);
        if (point == -1)
        {
            point = length();
            done = true;
        }
        else if (point == length() - 1)
            done = true;
        else if ((*this)[point+1] == '.' || (*this)[point+1] == '/')
            done = false;
        else
            done = true;
    }
    while (!done);
    FileName retval = *this;
    return retval.insert(point, str);
}

// Remove an extension wherever it is ......................................
FileName FileName::removeExtension(const std::string &ext) const
{
    int first = find((std::string)"." + ext);
    if (first == -1)
        return *this;
    else
    {
        FileName retval = *this;
        return retval.erase(first, 1 + ext.length());
    }
}

// Remove all extensions....................................................
FileName FileName::removeAllExtensions() const
{
    int first = rfind("/");
    first = find(".", first + 1);
    if (first == -1)
        return *this;
    else
        return substr(0, first);
}

FileName FileName::getFileFormat() const
{
    int first;
    FileName result;
    if (find("#") != std::string::npos)
        return "raw";
    else if ( (first = rfind(":"))!=std::string::npos)
        result = substr(first + 1) ;
    else if ( (first = rfind("."))!=std::string::npos)
        result = substr(first + 1);
    else
        result="spi";
    return result.toLowercase();

}

FileName FileName::removeFileFormat() const
{
    if ( find("#", 0) > -1 )
        REPORT_ERROR(ERR_IO,"Not implemented for raw data");
    size_t found=rfind(":");
    if (found!=std::string::npos)
        return substr(0, found);
    return *this;
}

bool FileName::isMetaData(bool failIfNotExists) const
{
    //file names containing @, : or % are not metadatas
    size_t found=this->find('@');
    if (found!=std::string::npos)
        return false;
    found=this->find(':');
    if (found!=std::string::npos)
        return false;
    found=this->find('#');
    if (found!=std::string::npos)
        return false;
    FileName ext = getFileFormat();
    //
    if (ext=="sel" || ext=="xmd" || ext=="doc")
    {
        return true;
    }
    else
    {
        std::ifstream infile(data(), std::ios_base::in);
        std::string line;

        if (infile.fail())
        {
            if (failIfNotExists)
                REPORT_ERROR( ERR_IO_NOTEXIST, (std::string) "File " + *this + " does not exist." );
            else
                return false;
        }

        // Search for xmipp_3,
        getline(infile, line, '\n');
        int pos = line.find("XMIPP_3 * ");

        if (pos != std::string::npos) // xmipp_3 token found
            return true;
        else
            return false;
    }
}

// Substitute one extension by other .......................................
FileName FileName::substituteExtension(const std::string &ext1,
                                        const std::string &ext2) const
{
    int first = find((std::string)"." + ext1);
    if (first == -1)
        return *this;
    else
    {
        FileName retval = *this;
        return retval.replace(first, 1 + ext1.length(), (std::string)"." + ext2);
    }
}

// Remove a substring ......................................................
FileName FileName::without(const std::string &str) const
{
    if (str.length() == 0)
        return *this;
    int pos = find(str);
    if (pos == -1)
        return *this;
    else
    {
        FileName retval = *this;
        return retval.erase(pos, str.length());
    }
}

// Remove until prefix .....................................................
FileName FileName::removeUntilPrefix(const std::string &str) const
{
    if (str.length() == 0)
        return *this;
    int pos = find(str);
    if (pos == -1)
        return *this;
    else
    {
        FileName retval = *this;
        return retval.erase(0, pos + str.length());
    }
}

// Remove directories ......................................................
FileName FileName::removeDirectories(int keep) const
{
    int last_slash = rfind("/");
    int tokeep = keep;
    while (tokeep > 0)
    {
        last_slash = rfind("/", last_slash - 1);
        tokeep--;
    }
    if (last_slash == -1)
        return *this;
    else
        return substr(last_slash + 1, length() - last_slash);
}
