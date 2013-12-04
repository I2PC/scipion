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

#include <dirent.h>
#include <unistd.h>
#include <algorithm>
#include "xmipp_filename.h"
#include "xmipp_funcs.h"
#include "xmipp_image_macros.h"
#include "xmipp_image_generic.h"


String FileNameVersion=METADATA_XMIPP_STAR;

void setMetadataVersion(String version)
{
    FileNameVersion=version;
}

String getMetadataVersion(void)
{
    return FileNameVersion;
}

// Constructor with root, number and extension .............................
void FileName::compose(const String &str, size_t no, const String &ext)
{

    if (no == ALL_IMAGES || no == (size_t) -1)
        REPORT_ERROR(ERR_DEBUG_TEST, "Don't compose with 0 or -1 index, now images index start at 1");

    *this = (FileName) str;

    if (no != ALL_IMAGES)
        this->append(formatString("%06lu", no));

    if (!ext.empty())
        *this += (String) "." + ext;
}

// Constructor: prefix number and filename, mainly for selfiles..
void FileName::compose(size_t no, const String &str)
{
    if (no == ALL_IMAGES || no == (size_t) -1)
        REPORT_ERROR(ERR_DEBUG_TEST, "Don't compose with 0 or -1 index, now images index start at 1");

    if (no != ALL_IMAGES)
    {
        size_t first = str.rfind(AT);
        if (first != npos)
        {
            std::vector<String> prefixes;
            int nPref = splitString(str.substr(0, first),",",prefixes, false);

            if (isalpha(prefixes[nPref-1].at(0)))
                formatStringFast(*this, "%06lu,%s", no, str.c_str());
        }
        else
            formatStringFast(*this, "%06lu@%s", no, str.c_str());
    }
    else
        *this = str;
}

// Constructor: prefix number, filename root and extension, mainly for selfiles..
void FileName::compose(size_t no, const String &str, const String &ext)
{
    if (no == ALL_IMAGES || no == (size_t) -1)
        REPORT_ERROR(ERR_DEBUG_TEST, "Don't compose with 0 or -1 index, now images index start at 1");

    if (no != ALL_IMAGES)
        formatStringFast(*this, "%06lu@%s.%s", no,
                         str.c_str(), ext.c_str());
    else
        *this = str;
}

// Constructor: string  and filename, mainly for metadata blocks..
void FileName::compose(const String &blockName, const String &str)
{
    if (blockName.empty())
        *this = str;
    else
        formatStringFast(*this, "%s@%s", blockName.c_str(), str.c_str());
}

// Constructor: string, number and filename, mainly for numered metadata blocks..
void FileName::composeBlock(const String &blockName, size_t no, const String &root, const String &ext)
{
    formatStringFast(*this, "%s%06lu@%s", blockName.c_str(), no, root.c_str());
    if (ext != "")
        *this += (String) "." + ext;
}

// Is in stack ............................................................
bool FileName::isInStack() const
{
    return find(AT) != String::npos;
}

// Decompose ..............................................................
void FileName::decompose(size_t &no, String &str) const
{
    char buffer[1024];
    unsigned long int auxNo;
    int ok = sscanf(c_str(), "%lu@%s", &auxNo, buffer);
    no=auxNo;
    if (ok != 2)
    {
        no = ALL_IMAGES;
        str = *this;
        return;
    }
    else if (no == 0)
        REPORT_ERROR(ERR_INDEX_OUTOFBOUNDS, formatString("FileName::decompose: Incorrect index number at filename %s; It must start at %lu",c_str(),FIRST_IMAGE));

    str = buffer;
}

// Get decomposed filename .......................................
FileName FileName::getDecomposedFileName() const
{
    String str;
    size_t no;
    decompose(no, str);
    return str;
}

// Get the root name of a filename .........................................
// TODO: Check if it is really needed
FileName FileName::getRoot() const
{
    size_t skip_directories = find_last_of("/") + 1;
    size_t point = find_first_of(".", skip_directories);
    if (point == npos)
        point = length();
    size_t root_end = find_last_not_of("0123456789", point - 1);
    if (root_end + 1 != point)
        if (point - root_end > FILENAMENUMBERLENGTH)
            root_end = point - FILENAMENUMBERLENGTH - 1;
    return substr(0, root_end + 1);
}

// Convert to lower case characters .........................................
FileName FileName::toLowercase() const
{
    FileName result = *this;
    for (size_t i = 0; i < result.length(); i++)
        result[i] = tolower(result[i]);
    return result;
}

// Convert to upper case characters .........................................
FileName FileName::toUppercase() const
{
    FileName result = *this;
    for (size_t i = 0; i < result.length(); i++)
        result[i] = toupper(result[i]);
    return result;
}

// Is substring present?
bool FileName::contains(const String& str) const
{
    return find(str) != npos;
}

// Get substring before first instance of str
//TODO: Check behaviour
FileName FileName::beforeFirstOf(const String& str) const
{
    size_t point = find_first_of(str);
    return (point != npos ? (FileName)substr(0, point) : *this);
}

// Get substring before last instance of str
//TODO: Check behaviour
FileName FileName::beforeLastOf(const String& str) const
{
    size_t point = find_last_of(str);
    return point != npos ? (FileName)substr(0, point) : *this;
}

// Get substring after first instance of str
//TODO: Check behaviour
FileName FileName::afterFirstOf(const String& str) const
{
    size_t point = find_first_of(str);
    return point != npos ? (FileName)substr(point + 1) : *this;
}

// Get substring after last instance of str
//TODO: Check behaviour
FileName FileName::afterLastOf(const String& str) const
{
    size_t point = find_last_of(str);
    return point != npos ? (FileName)substr(point + 1) : *this;
}

// Get the base name of a filename .........................................
FileName FileName::getBaseName() const
{
    FileName baseName = removeLastExtension();
    return baseName.afterLastOf("/");
}


// Get the dir of a filename .........................................
FileName FileName::getDir() const
{
    size_t pos = find_last_of("/");
    return (FileName)( pos != npos ? substr(0, pos+1) : "");
}

// Get the extension of a filename .........................................
String FileName::getExtension() const
{
    size_t posA = find_last_of("/");
    size_t posB = find_last_of(".");
    if (posB==npos)
        return "";
    if (posA==npos)
        return substr(posB+1);
    if (posB>posA)
        return substr(posB+1);
    return "";
}

// Has image extension .....................................................
bool FileName::hasImageExtension() const
{
    String ext = getFileFormat();
    return (ext=="img" || ext=="hed" || ext=="inf" || ext=="raw" || ext=="mrc" ||
            ext=="map" || ext=="spi" || ext=="xmp" || ext=="tif" || ext=="dm3" ||
            ext=="spe" || ext=="em"  || ext=="pif" || ext=="ser" || ext=="stk" ||
            ext=="mrcs"|| ext=="jpg");
}

// Has image extension .....................................................
bool FileName::hasStackExtension() const
{
    String ext = getFileFormat();
    return (ext=="stk" || ext=="spi" || ext=="xmp" || ext=="mrcs" || ext=="mrc" ||
            ext=="img" || ext=="hed" || ext=="pif" || ext=="tif"  || ext=="dm3" ||
            ext=="ser" || ext=="st");
}

// Has image extension .....................................................
bool FileName::hasVolumeExtension() const
{
    String ext = getFileFormat();
    return (ext=="vol" || ext=="spi" || ext=="xmp" || ext=="mrc" || ext=="map" ||
            ext=="em"  || ext=="pif" || ext=="inf" || ext=="raw");
}

// Has image extension .....................................................
bool FileName::hasMetadataExtension() const
{
    String ext = getFileFormat();
    return (ext == "sel"    || ext == "xmd" || ext == "doc" ||
            ext == "ctfdat" || ext == "ctfparam" || ext == "pos" ||
            ext == "sqlite" || ext == "xml");
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
void FileName::initUniqueName(const char *templateStr, const String &fnDir)
{
#ifndef __MINGW32__
    int fd;
    const int len=512;
    char filename[len];
    if (fnDir!="")
        strcpy(filename,(fnDir+"/").c_str());
    else
        filename[0]=0;
    strcat(filename, templateStr);
    filename[len - 1] = 0;
    if ((fd = mkstemp(filename)) == -1)
    {
        perror("FileName::Error generating tmp lock file");
        exit(1);
    }
    close(fd);
    *this = filename;
#endif
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
        retval = retval.append((String) "." + ext);
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
    return removeSubstring(getRoot());
}

// Insert before extension .................................................
FileName FileName::insertBeforeExtension(const String &str) const
{
    FileName retval = *this;
    size_t pos = find_last_of('.');
    return  pos != npos ? retval.insert(pos, str) : retval.append(str);
}

// Remove an extension wherever it is ......................................
FileName FileName::removeExtension(const String &ext) const
{
    FileName retval = *this;
    size_t first = find((String) "." + ext);
    return (first == npos) ? retval: retval.erase(first, 1 + ext.length());
}

// Remove the last extension ................................................
FileName FileName::removeLastExtension() const
{
    FileName retval = *this;
    size_t first = find_last_of('.');
    return (first == npos) ? retval : retval.substr(0, first);
}

// Remove all extensions....................................................
FileName FileName::removeAllExtensions() const
{
    FileName retval = *this;
    size_t first = find_last_of('/');
    first = find_first_of('.', first + 1);
    return (first == npos) ? retval: retval.substr(0, first);
}

FileName FileName::removeFilename() const
{
    size_t first = find_last_of('/');
    return (first == npos) ? "" : substr(0, first);
}

String FileName::getFileFormat() const
{
    size_t first;
    FileName result;
    if (find(NUM) != npos)
        return "raw";
    else if ((first = rfind(COLON)) != npos)
        result = substr(first + 1);
    else if ((first = rfind(".")) != npos)
        result = substr(first + 1);
    return result.toLowercase();
}

size_t FileName::getFileSize() const
{
    Stat info;
    if (stat(c_str(), &info))
    {
        char cCurrentPath[FILENAME_MAX];
        char *success=getcwd(cCurrentPath, sizeof(cCurrentPath));
        if (success==NULL)
        	cCurrentPath[0]='\0';
        REPORT_ERROR(ERR_UNCLASSIFIED,formatString("FileName::getFileSize: Cannot get size of file %s/%s",cCurrentPath,this->c_str()));
    }
    return info.st_size;
}

FileName FileName::removeFileFormat() const
{
    size_t found = rfind(NUM);
    if (found != String::npos)
        return substr(0, found);
    found = rfind(COLON);
    if (found != String::npos)
        return substr(0, found);
    return *this;
}

// Get number from file base name ....................................................
int FileName::getNumber() const
{
    size_t skip_directories = find_last_of("/") + 1;
    size_t point = find_first_of(".", skip_directories);
    if (point == npos)
        point = length();
    size_t root_end = find_last_not_of("0123456789", point - 1);
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

// Get number from file ....................................................
size_t FileName::getPrefixNumber(size_t pos) const
{
    size_t first = rfind(AT);
    size_t result = ALL_IMAGES;
    if (first != npos)
    {
        std::vector<String> prefixes;
        size_t nPref = splitString(substr(0, first),",",prefixes, false);

        if (pos > nPref-1)
            REPORT_ERROR(ERR_ARG_INCORRECT, formatString("getPrefixNumber: Selected %lu position greater than %lu positions \n"
                         " detected in %s filename.",pos+1, nPref, this->c_str()));

        if (isdigit(prefixes[pos].at(0)))
            result = textToSizeT(prefixes[pos].c_str());
    }
    return result;
}

String FileName::getBlockName() const
{
    size_t first = rfind(AT);
    String result = "";
    if (first != npos)
    {
        result = substr(0, first);
        if ((first = result.find(COMMA)) != npos) // Assign and compare at the same time
          result = result.substr(first+1);

        /* using isdigit instead of isalpha allows to
         * detect as blockname rootnames starting by "/"
         */
        if (result.empty() || isdigit(result[0]))
            result = "";
    }
    return result;

}

FileName FileName::removeBlockName() const
{
    size_t first = rfind(AT);

    if (first != npos)
    {
        String block = substr(0, first);
        size_t second = block.find(COMMA);

        if (second == npos)
        {
          if (!isdigit(block[0])){
            return substr(first + 1);
          }
        }
        else
        {
          String prefix = block.substr(0, second);
          block = block.substr(second + 1);
          if (!isdigit(block[0]))
            return prefix + substr(first);
        }
    }
    return *this;
}

FileName FileName::removePrefixNumber() const
{
    size_t first = rfind(AT);

    if (first != npos)
    {
        std::vector<String> prefixes;
        int nPref = splitString(substr(0, first),",",prefixes, false);

        if (isdigit(prefixes[nPref-1].at(0)))
            return substr(first + 1);
        else if (nPref > 1) // isalpha and we remove the ","
            return substr(first - prefixes[nPref-1].size());
    }
    return *this;
}

FileName FileName::removeAllPrefixes() const
{
    size_t first = rfind(AT);
    if (first != npos)
        return substr(first + 1);
    return *this;
}

bool FileName::isMetaData(bool failIfNotExists) const
{
    //check empty string
    if (empty())
        REPORT_ERROR(ERR_ARG_INCORRECT, "FileName::isMetaData: Empty string is not a MetaData");
    //file names containing : or % are not metadatas
    //size_t found = this->find('@');
    if (find_first_of(":#") != npos)
        return false;

    //check if file exists
    if (failIfNotExists && !existsTrim())
        REPORT_ERROR(ERR_IO_NOTFILE, formatString("FileName::isMetaData: File: '%s' does not exist", c_str()));
    //This is dangerous and should be removed
    //in next version. only star1 files should be OK
    //ROB
    //FIXME
    return (hasMetadataExtension() || isStar1(failIfNotExists));
}

bool FileName::isStar1(bool failIfNotExists) const
{
    std::ifstream infile( this->removeAllPrefixes().c_str(), std::ios_base::in);
    String line;

    if (infile.fail())
    {
        if (failIfNotExists)
            REPORT_ERROR( ERR_IO_NOTEXIST, formatString("File '%s' does not exist.", this->removeAllPrefixes().c_str()));
        else
            return false;
    }

    // Search for xmipp_3,
    char cline[128];
    infile.getline(cline, 128);
    infile.close();
    line = cline;
    size_t pos = line.find(METADATA_XMIPP_STAR);
    return (pos != npos); // xmipp_star_1 token found
}

// Replace one substring by other .......................................
FileName FileName::replaceSubstring(const String &subOld, const String &subNew) const
{
    size_t pos = find(subOld);
    if (pos == npos)
        return *this;

    FileName result = *this;
    result.replace(pos, subOld.length(), subNew);
    return result;
}

// Substitute one extension by other .......................................
FileName FileName::replaceExtension(const String &newExt) const
{
    return removeLastExtension() + "." + newExt;
}

// Remove a substring ......................................................
FileName FileName::removeSubstring(const String &sub) const
{
    return replaceSubstring(sub, "");
}

// Remove until prefix .....................................................
FileName FileName::removeUntilPrefix(const String &prefix) const
{
    size_t pos = find(prefix);
    if (pos == npos)
        return *this;
    FileName result = *this;
    return result.erase(0, pos + prefix.length());
}

// Remove directories ......................................................
FileName FileName::removeDirectories(int keep) const
{
    size_t last_slash = rfind("/");
    int tokeep = keep;
    while (tokeep > 0)
    {
        last_slash = rfind("/", last_slash - 1);
        tokeep--;
    }
    if (last_slash == npos)
        return *this;
    else
        return substr(last_slash + 1, length() - last_slash);
}

void FileName::copyFile(const FileName & target) const
{
    std::ifstream f1(this->c_str(), std::fstream::binary);
    std::ofstream
    f2(target.c_str(), std::fstream::trunc | std::fstream::binary);
    f2 << f1.rdbuf();
}

/* Check if a file exists -------------------------------------------------- */
bool FileName::exists() const
{
    return fileExists(getDecomposedFileName().removeFileFormat());
}
/* Delete  file exists -------------------------------------------------- */
void FileName::deleteFile() const
{
    FileName temp = this->removeFileFormat().removeAllPrefixes();
    if (temp.exists())
        unlink(temp.c_str());
}
/* Check if a file exists remove leading @ and tailing : */
bool FileName::existsTrim() const
{
    FileName auxF(*this);
    size_t found = find_first_of(AT);

    if (found != String::npos)
        auxF =  substr(found+1);

    found = auxF.find_first_of(NUM);

    if ( found != String::npos)
        auxF = auxF.substr(0, found);
    found = auxF.find_first_of(COLON);

    if (found != String::npos)
        auxF = auxF.substr(0, found);
    return fileExists(auxF.c_str());
}

/* List of files within a directory ---------------------------------------- */
void FileName::getFiles(std::vector<FileName> &files) const
{
    files.clear();

    DIR *dp;
    struct dirent *dirp;
    if ((dp  = opendir(c_str())) == NULL)
        REPORT_ERROR(ERR_IO_NOTEXIST,*this);

    while ((dirp = readdir(dp)) != NULL)
        if (strcmp(dirp->d_name,".")!=0 && strcmp(dirp->d_name,"..")!=0)
            files.push_back(FileName(dirp->d_name));
    closedir(dp);
    std::sort(files.begin(),files.end());
}

/* Is directory ------------------------------------------------------------ */
bool FileName::isDir() const
{
    Stat st_buf;
    if (stat (c_str(), &st_buf) != 0)
        REPORT_ERROR(ERR_UNCLASSIFIED,(String)"Cannot determine status of filename "+ *this);
    return (S_ISDIR (st_buf.st_mode));
}

/* Wait until file has a stable size --------------------------------------- */
void FileName::waitUntilStableSize(size_t time_step)
{
    size_t idx;
    FileName basicName;
    decompose(idx, basicName);

    if (!exists())
        return;
    Stat info1, info2;
    if (stat(basicName.c_str(), &info1))
        REPORT_ERROR(ERR_UNCLASSIFIED,
                     (String)"FileName::waitUntilStableSize: Cannot get size of file " + *this);
    off_t size1 = info1.st_size;
    do
    {
        usleep(time_step);
        if (stat(basicName.c_str(), &info2))
            REPORT_ERROR(ERR_UNCLASSIFIED,
                         (String)"FileName::waitUntilStableSize: Cannot get size of file " + *this);
        off_t size2 = info2.st_size;
        if (size1 == size2)
            break;
        size1 = size2;
    }
    while (true);
    return;
}

/* Create empty file ------------------------------------------------------- */
void FileName::createEmptyFile(size_t size, size_t block_size)
{
    unsigned char * buffer = (unsigned char*) calloc(sizeof(unsigned char),
                             block_size);
    if (buffer == NULL)
        REPORT_ERROR(ERR_MEM_NOTENOUGH, "create_empty_file: No memory left");
    FILE * fd = fopen(c_str(), "w");
    if (fd == NULL)
        REPORT_ERROR(ERR_IO_NOTOPEN, (String)"FileName::createEmptyFile: Cannot open file" + *this);
    for (size_t i = 0; i < size / block_size; i++)
        fwrite(buffer, sizeof(unsigned char), block_size, fd);
    fwrite(buffer, sizeof(unsigned char), size % block_size, fd);
    fclose(fd);
}

void FileName::createEmptyFileWithGivenLength(size_t length) const
{
    FILE* fMap = fopen(c_str(),"wb");
    if (!fMap)
        REPORT_ERROR(ERR_IO_NOWRITE, *this);
    if (length>0)
    {
        char c=0;
        if ((fseek(fMap, length-1, SEEK_SET) == -1) || (fwrite(&c,1,1,fMap) != 1))
            REPORT_ERROR(ERR_IO_NOWRITE,"FileName::createEmptyFileWithGivenLength: Cannot create empty file");
    }
    fclose(fMap);
}

/** Auxiliary function used to create a tree of directories
 * return 0 on success
 */
int do_mkdir(const char *path, mode_t mode)
{
    Stat st;
    int status = 0;

    if (stat(path, &st) != 0)
    {
        /* Directory does not exist */
#ifndef __MINGW32__
        if (mkdir(path, mode) != 0)
#else

            if (mkdir(path) != 0)
#endif

                status = -1;
    }
    else if (!S_ISDIR(st.st_mode))
    {
        errno = ENOTDIR;
        status = -1;
    }

    return (status);
}

int FileName::makePath(mode_t mode) const
{
    char *pp;
    char *sp;
    int status;
    char *copypath = strdup(c_str());
    if (copypath == NULL)
        REPORT_ERROR(ERR_MEM_BADREQUEST,"FileName::makePath: Canot alloc memory");

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
        status = do_mkdir(c_str(), mode);
    free(copypath);
    return status;
}

/* Exit program if filename is not empry and file does not exist ----------- */
void FileName::assertExists()
{
    if (!empty() && !exists())
    {
        std::cerr << "FileName::assertExists: control file" << *this
        << " doesn't exist, exiting..." << std::endl;
        exit(ERR_IO_NOTEXIST);
    }
    //TODO: Maybe change to report error???
    //REPORT_ERROR(ERR_IO_NOTEXIST, (String)"FileName::assertExists: control file" + *this " doesn't exist, exiting...");
}

/* Get the Xmipp Base directory -------------------------------------------- */
char * getXmippPath()
{
    char* path = getenv("XMIPP_HOME");
    if (path == NULL)
        REPORT_ERROR(ERR_VALUE_EMPTY, "getXmippPath::Variable XMIPP_HOME is not defined");
    return path;

}

void copyImage(const FileName & source, const FileName & target)
{
    ImageGeneric img(source);
    img.write(target);
}

void deleteFile(const FileName &fn)
{
	fn.deleteFile();
}

void FileLock::lock(int _fileno)
{
#ifndef __MINGW32__
    if (islocked)
        unlock();

    if (_fileno != 0)
        filenum = _fileno;

    fl.l_type = F_WRLCK;
    fcntl(filenum, F_SETLKW, &fl);
    islocked = true;
#endif

}

void FileLock::lock(FILE * hdlFile)
{
    if (islocked)
        unlock();

    if (hdlFile != NULL)
        this->filenum = fileno(hdlFile);

#ifdef __MINGW32__

    HANDLE hFile = (HANDLE)_get_osfhandle(filenum);
    DWORD dwLastPos = SetFilePointer(hFile, 0, NULL, FILE_END);
    if (LockFile(hFile, 0, 0, dwLastPos, 0) != NULL)
        REPORT_ERROR(ERR_IO_LOCKED,"File cannot be locked.");
#else

    fl.l_type = F_WRLCK;
    fcntl(filenum, F_SETLKW, &fl);
#endif

    islocked = true;
}


void FileLock::unlock()
{
    if (islocked)
    {
#ifdef __MINGW32__
        HANDLE hFile = (HANDLE)_get_osfhandle(filenum);
        DWORD dwLastPos = SetFilePointer(hFile, 0, NULL, FILE_END);
        if (UnlockFile(hFile, 0, 0, dwLastPos, 0) != NULL)
            REPORT_ERROR(ERR_IO_LOCKED,"File cannot be unlocked.");
#else

        fl.l_type = F_UNLCK;
        fcntl(filenum, F_SETLK, &fl);
#endif

        islocked = false;
    }
}

