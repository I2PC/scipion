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
void FileName::compose(const String &str, size_t no, const String &ext)
{
    if (no == ALL_IMAGES || no == (size_t) -1)
        REPORT_ERROR(ERR_DEBUG_TEST, "Don't compose with 0 or -1 index, now images index start at 1");

    *this = (FileName) str;

    if (no != ALL_IMAGES)
        this->append(formatString("%06lu", no));

    if (ext != "")
        *this += (String)"." + ext;
}

// Constructor: prefix number and filename, mainly for selfiles..
void FileName::compose(size_t no , const String &str)
{
    if (no == ALL_IMAGES || no == (size_t) -1)
        REPORT_ERROR(ERR_DEBUG_TEST, "Don't compose with 0 or -1 index, now images index start at 1");

    if (no != ALL_IMAGES)
        formatStringFast(*this,"%06lu@%s", no, str.c_str());
    else
        *this = str;
}

// Constructor: prefix number, filename root and extension, mainly for selfiles..
void FileName::compose(size_t no , const String &str , const String &ext)
{
    if (no == ALL_IMAGES || no == (size_t) -1)
        REPORT_ERROR(ERR_DEBUG_TEST, "Don't compose with 0 or -1 index, now images index start at 1");

    if (no != ALL_IMAGES)
        formatStringFast(*this,"%06lu@%s.%s", no, str.c_str(), ext.c_str());
    else
        *this = str;
}

// Constructor: string  and filename, mainly for metadata blocks..
void FileName::compose(const String &blockName , const String &str)
{
    *this = (FileName)( blockName + (String)"@" + str);
}

// Is in stack ............................................................
bool FileName::isInStack() const
{
    return find("@") != String::npos;
}

// Decompose ..............................................................
void FileName::decompose(size_t &no, String &str) const
{
    char buffer[1024];
    int ok = sscanf(c_str(), "%lu@%s", &no,&buffer);
    if (ok!=2)
    {
        no = ALL_IMAGES;
        str = *this;
        return;
    }
    str = buffer;
}

// Get decomposed filename .......................................
String FileName::getDecomposedFileName() const
{
    String str;
    size_t no;
    decompose(no, str);
    return str;
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
bool FileName::contains(const String& str) const
{
    int point = rfind(str);
    if (point > -1)
        return true;
    else
        return false;
}

// Get substring before first instance of str
FileName FileName::beforeFirstOf(const String& str) const
{
    int point = find_first_of(str);
    if (point > -1)
        return substr(0, point);
    else
        return *this;
}

// Get substring before last instance of str
FileName FileName::beforeLastOf(const String& str) const
{
    int point = find_last_of(str);
    if (point > -1)
        return substr(0, point);
    else
        return *this;
}

// Get substring after first instance of str
FileName FileName::afterFirstOf(const String& str) const
{
    int point = find_first_of(str);
    if (point > -1)
        return substr(point + 1);
    else
        return *this;
}

// Get substring after last instance of str
FileName FileName::afterLastOf(const String& str) const
{
    int point = find_last_of(str);
    if (point > -1)
        return substr(point + 1);
    else
        return *this;
}

// Get the base name of a filename .........................................
String FileName::getBaseName() const
{
    String basename = "";
    String myname = *this;
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
        String aux = substr(root_end + 1, point - root_end + 1);
        return atoi(aux.c_str());
    }
    else
        return -1;
}

// Get the extension of a filename .........................................
String FileName::getExtension() const
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

// Init Unique .............................................................
void FileName::initUniqueName(const char *templateStr)
{
    int fd;
    char filename[L_tmpnam];
    strcpy(filename, templateStr);
    filename[L_tmpnam-1]=0;

    if ((fd = mkstemp(filename)) == -1)
    {
        perror("FileName::Error generating tmp lock file");
        exit(1);
    }
    close(fd);
    *this = filename;
}

// Add at beginning ........................................................
FileName FileName::addPrefix(const String &prefix) const
{
    FileName retval = *this;
    int skip_directories = find_last_of("/") + 1;
    return retval.insert(skip_directories, prefix);
}

// Add at the end ..........................................................
FileName FileName::addExtension(const String &ext) const
{
    if (ext == "")
        return *this;
    else
    {
        FileName retval = *this;
        retval = retval.append((String)"." + ext);
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
FileName FileName::insertBeforeExtension(const String &str) const
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
FileName FileName::removeExtension(const String &ext) const
{
    int first = find((String)"." + ext);
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
    if (find("#") != String::npos)
        return "raw";
    else if ( (first = rfind(":"))!=String::npos)
        result = substr(first + 1) ;
    else if ( (first = rfind("."))!=String::npos)
        result = substr(first + 1);
    return result.toLowercase();
}

FileName FileName::removeFileFormat() const
{
    size_t found=rfind("#");
    if (found!=String::npos)
        return substr(0, found);
    found=rfind(":");
    if (found!=String::npos)
        return substr(0, found);
    return *this;
}

String FileName::getBlockName() const
{
    size_t first = rfind("@");
    String result = "";
    if ( first != npos && isalpha(this->at(0)))
        result = substr(0, first);
    return result;

}

FileName FileName::removeBlockName() const
{
    size_t first = rfind("@");
    if ( first != npos && isalpha(this->at(0)))
        return substr(first + 1);
    return *this;
}

bool FileName::isMetaData(bool failIfNotExists) const
{
    //file names containing @, : or % are not metadatas
    size_t found = this->find('@');
    if (find_first_of("@:#") != npos)
        return false;
    //check if file exists
    if (!exists(*this))
        REPORT_ERROR(ERR_IO_NOTFILE,(String)"file: " + *this + (String)" does not exist");
    //This is dangerous and should be removed
    //in next version. only star1 files should be OK
    //ROB
    //FIXME
    FileName ext = getFileFormat();
    if (ext=="sel" || ext=="xmd" || ext=="doc")
    {
        return true;
    }
    else
        return isStar1(failIfNotExists);
}

bool FileName::isStar1(bool failIfNotExists) const
{
    std::ifstream infile(data(), std::ios_base::in);
    String line;

    if (infile.fail())
    {
        if (failIfNotExists)
            REPORT_ERROR( ERR_IO_NOTEXIST, formatString("File '%s' does not exist.", c_str()));
        else
            return false;
    }

    // Search for xmipp_3,
    char cline[128];
    infile.getline(cline, 128);
    line = cline;
    int pos = line.find("XMIPP_STAR_1 *");

    return (pos != npos); // xmipp_star_1 token found
}

// Substitute one extension by other .......................................
FileName FileName::substituteExtension(const String &ext1,
                                       const String &ext2) const
{
    int first = find((String)"." + ext1);
    if (first == -1)
        return *this;
    else
    {
        FileName retval = *this;
        return retval.replace(first, 1 + ext1.length(), (String)"." + ext2);
    }
}

// Remove a substring ......................................................
FileName FileName::without(const String &str) const
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
FileName FileName::removeUntilPrefix(const String &str) const
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

void FileName::copyFile(const FileName & target) const
{
    std::ifstream f1 (this->c_str(), std::fstream::binary);
    std::ofstream f2 (target.c_str(),std::fstream::trunc|std::fstream::binary);
    f2<<f1.rdbuf();
}



typedef struct stat Stat;
//return 0 on success
int do_mkdir(const char *path, mode_t mode)
{
    Stat            st;
    int             status = 0;

    if (stat(path, &st) != 0)
    {
        /* Directory does not exist */
        if (mkdir(path, mode) != 0)
            status = -1;
    }
    else if (!S_ISDIR(st.st_mode))
    {
        errno = ENOTDIR;
        status = -1;
    }

    return(status);
}

/**
** mkpath - ensure all directories in path exist
** Algorithm takes the pessimistic view and works top-down to ensure
** each directory in path exists, rather than optimistically creating
** the last element and working backwards. Return null if fails
*/
int mkpath(const FileName &path, mode_t mode)
{
    char           *pp;
    char           *sp;
    int             status;
    char           *copypath = strdup(path.c_str());
    if(copypath==NULL)
        REPORT_ERROR(ERR_MEM_BADREQUEST,"Mkpath: Canot alloc memory");

    status = 0;
    pp = copypath;
    while (status == 0 && (sp = strchr(pp, '/')) != 0)
    {
        if (sp != pp)
        {
            /* Neither root nor double slash in path */
            *sp = '\0';
            status = do_mkdir(copypath, mode);
            *sp = '/';
        }
        pp = sp + 1;
    }
    if (status == 0)
        status = do_mkdir(path.c_str(), mode);
    free(copypath);
    return (status);
}
