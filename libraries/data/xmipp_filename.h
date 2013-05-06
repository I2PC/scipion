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

#ifndef XMIPP_FILENAME_H_
#define XMIPP_FILENAME_H_

#include <vector>
#include <sys/stat.h>

//#ifdef HAVE_UNISTD_H
#include <unistd.h>
//#endif /* HAVE_UNISTD_H */

#include "xmipp_strings.h"

#ifndef __MINGW32__
#include <fcntl.h>
#endif

#define FILENAMENUMBERLENGTH 6

typedef struct stat Stat;

//extern variable with magic Word  used to define metadata version
// xmipp XMIPP_STAR_1
// sqlite SQLite format 3

#define METADATA_XMIPP_STAR "# XMIPP_STAR_1"
#define METADATA_XMIPP_SQLITE "SQLite format 3"

//@{
/** Filenames.
 *
 * This class allows you a lot of usual and common manipulations with filenames.
 * See filename conventions for a detailed explanation of the Filenames dealed
 * here, although most of the functions work with the more general model
 * "name.extension"
 */
class FileName: public String
{
public:
    /// @name Filename constructors
    /// @{

    /** Empty constructor
     *
     * The empty constructor is inherited from the string class, so an empty
     * FileName is equal to "".
     *
     * @code
     * FileName fn_blobs;
     * @endcode
     */
    inline FileName(): String("")
    {}

    /** Constructor from string
     *
     * The constructor from a string allows building complex expressions based
     * on the string class. Notice that in the following example the type
     * casting to string is very important, if not, the operation is just a
     * pointer movement instead of a string concatenation.
     *
     * @code
     * FileName fn_blobs((String) "art00001" + ".blobs");
     * @endcode
     */
    FileName(const String& str): String(str)
    {}

    /** Constructor from char*
     */
    FileName(const char* str) : String(str)
    {}

    /** Assignment constructor
     */
    FileName& operator=(const FileName& op)
    {
        return (FileName&) std::string::operator=(op);
    }

    /** Copy constructor
     */
    FileName(const FileName& fn) : String(fn)
    {}

    /** Constructor from root, number and extension
     *
     * The number and extension are optional.
     *
     * @code
     * FileName fn_proj("g1ta000001.xmp"); // fn_proj = "g1ta000001.xmp"
     * FileName fn_proj("g1ta",1,"xmp"); // fn_proj = "g1ta000001.xmp"
     * FileName fn_proj("g1ta",1); // fn_proj = "g1ta000001"
     * @endcode
     */
    FileName(const char* str, int no, const String& ext = "")
    {
        compose(str, no, ext);
    }

    /** Constructor from root and extension
     *
     * None of the parameters is optional
     *
     * @code
     * FileName fn_proj("g1ta00001", "xmp"); // fn_proj = "g1ta00001.xmp"
     * @endcode
     */
    FileName(const char* str, const String& ext): String(str + ext)
    {}
    //@}

    /// @name Composing/Decomposing the filename
    /// @{

    /**
     * Convert to String
     */
    inline String getString() const
    {
        return (String)(*this);
    }

    /** Compose from root, number and extension
     *
     * @code
     * fn_proj.compose("g1ta", 1, "xmp");  // fn_proj = "g1ta000001.xmp"
     * @endcode
     */
    void compose(const String& str, const size_t no, const String& ext = "");

    /** Prefix with number @. Mainly for selfiles
     *
     * @code
     * fn_proj.compose(1,"g1ta.xmp");  // fn_proj = "000001@g1ta000001.xmp"
     * @endcode
     */
    void compose(const size_t no, const String& str);


    // Constructor: prefix number, filename root and extension, mainly for selfiles..
    /** Prefix with number and extension @.
     *
     * @code
     * fn_proj.compose(1, "g1ta", "xmp");  // fn_proj = "000001@g1ta000001.xmp"
     * @endcode
     */
    void compose(size_t no , const String &str , const String &ext);

    // Constructor: string  and filename, mainly for metadata blocks..
    /** Prefix with string (blockname) .
     *
     * @code
     * fn_proj.compose(b00001, "g1ta.xmp");  // fn_proj = "b000001@g1ta.xmp"
     * @endcode
     */
    void compose(const String &blockName , const String &str);

    /** Constructor: string, number, rootfilename and extension:
     *  mainly for numered metadata blocks..
     *
     * @code
     * fn_proj.composeBlock("bb",5, "g1ta","xmp");  // fn_proj = "bb000005@g1ta.xmp"
     * fn_proj.composeBlock("bb",5, "g1ta.xmp");    // fn_proj = "bb000005@g1ta.xmp"
     * @endcode
    */
    void composeBlock(const String &blockName, size_t no, const String &str, const String &ext="");


    /** True if this filename belongs to a stack
     */
    bool isInStack() const;

    /** Decompose filenames with @. Mainly from selfiles
     *
     * @code
     * fn_proj.decompose(no,filename);  // fn_proj = "000001@g1ta000001.xmp"
     *                                  // no=1
     *                                  // filename = "g1ta000001.xmp"
     * @endcode
     */
    void decompose(size_t &no, String& str) const;

    /** Get decomposed filename from filenames with @. Mainly from selfiles
         *
         * @code
         * filename = fn_proj.decomposedFileName();  // fn_proj = "000001@g1ta000001.xmp"
         *                                          // no=1
         *                                          // filename = "g1ta000001.xmp"
         * @endcode
         */
    FileName getDecomposedFileName() const;

    /** Get the root from a filename
     *
     * @code
     * FileName fn_root, fn_proj("g1ta00001");
     * fn_root = fn_proj.get_root();
     * @endcode
     */
    FileName getRoot() const;

    /** Get the base name from a filename
     */
    FileName getBaseName() const;

    /** Get the number from a filename
     *
     * If there is no number a -1 is returned.
     *
     * @code
     * FileName proj("g1ta00001");
     * int num = proj.getNumber();
     * @endcode
     */
    int getNumber() const;

    /** Get the number from a stack filename
        *
        * If there is no number a 0 is returned.
        *
        * @code
        * FileName proj("24@images.stk");
        * size_t num = proj.getStackNumber();
        * @endcode
        */
    size_t getStackNumber() const;

    FileName getDir() const;

    /** Get the last extension from filename
     *
     * The extension is returned without the dot. If there is no extension "" is
     * returned.
     *
     * @code
     * String ext = fn_proj.get_extension();
     * @endcode
     */
    String getExtension() const;

    /** Get image format identifier
     *
     * @code
     * fn_proj = "g1ta00001.xmp";
     * fn_proj = fn_proj.get_file_format(); // fn_proj == "xmp"
     * fn_proj = "g1ta00001.nor:spi";
     * fn_proj = fn_proj.get_file_format(); // fn_proj == "spi"
     * fn_proj = "input.file#120,120,55,1024,float";
     * fn_proj = fn_proj.get_file_format(); // fn_proj == "raw"
     * @endcode
     */
    String getFileFormat() const;

    /** True if the extension of this filename is of an image type
     */
    bool hasImageExtension() const;

    /** True if the extension of this filename is of a stack type
     */
    bool hasStackExtension() const;

    /** True if the extension of this filename is of a stack type
     */
    bool hasVolumeExtension() const;

    /** True if the extension of this filename is of a metadata file
     */
    bool hasMetadataExtension() const;

    /** Get the size of the file.
     */
    size_t getFileSize() const;


    /** Get blockName from filename
     * @code
     * fn_meta="block1@md1.doc"
     * String blockName;
     * blockName=fn_meta.getblockName();//blockName="block1"
     * @endcode
     */

    String getBlockName() const;

    /** Remove leading number from filename
     * @code
     * fn_meta="1@md1.doc"
     * String blockName;
     * filename=fn_meta.removeSliceNumber();//filename="md1.doc"
     * @endcode
     */

    FileName removeSliceNumber() const;

    /** Remove blockName from filename
     * @code
     * fn_meta="block1@md1.doc"
     * String blockName;
     * filename=fn_meta.getblockName();//filename="md1.doc"
     * @endcode
     */

    FileName removeBlockName() const;

    /** Remove blockName or slice number from filename
    * @code
    * fn_meta="block1@md1.doc"
    * String blockName;
    * filename=removeBlockNameOrSliceNumber();//filename="md1.doc"
    * @endcode
    */
    FileName removeBlockNameOrSliceNumber() const;

    /** Random name
     *
     * Generate a random name of the desired length.
     */
    void initRandom(int length);

    /** Unique name
     *
     * Generate a unique name replacing each 'X' with a character from
     * the portable filename character set. The characters are chosen such
     * that the resulting name does not duplicate the name of an existing file.
     */
    void initUniqueName(const char * templateStr = "xmippTemp_XXXXXX", const String &fnDir="");
    //@}

    ///@name Filename utilities
    //@{
    /** Change all characters for lowercases
     *
     * @code
     * FileName fn_proj("g1tA00001");
     * fn_proj = fn_proj.to_lowercase(); // fn_proj = "g1ta00001"
     * @endcode
     */
    FileName toLowercase() const;

    /** Change all characters for uppercases
     *
     * @code
     * FileName fn_proj("g1tA00001");
     * fn_proj = fn_proj.to_uppercase(); // fn_proj = "G1Ta00001"
     * @endcode
     */
    FileName toUppercase() const;

    /** Check whether the filename contains the argument substring
     *
     * @code
     * FileName fn_proj("g1ta00001.raw#d=f");
     * if (fn_proj.contains("raw) )  // true
     * @endcode
     */
    bool contains(const String& str) const;

    /** Return substring before first instance of argument (as in Bsoft)
     *
      * @code
     * FileName fn_proj("g1ta00001.raw#d=f");
     * fn_proj = fn_proj.before_first_of("#"); // fn_proj = "g1ta00001.raw"
     * @endcode
     */
    FileName beforeFirstOf(const String& str) const;

    /** Return substring before last instance of argument (as in Bsoft)
     *
      * @code
     * FileName fn_proj("g1ta00001.raw#d=f");
     * fn_proj = fn_proj.before_last_of("#"); // fn_proj = "g1ta00001.raw"
     * @endcode
     */
    FileName beforeLastOf(const String& str) const;

    /** Return substring after first instance of argument (as in Bsoft)
     *
      * @code
     * FileName fn_proj("g1ta00001.raw#d=f");
     * fn_proj = fn_proj.after_first_of("#"); // fn_proj = "d=f"
     * @endcode
     */
    FileName afterFirstOf(const String& str) const;

    /** Return substring after last instance of argument (as in Bsoft)
     *
      * @code
     * FileName fn_proj("g1ta00001.raw#d=f");
     * fn_proj = fn_proj.after_last_of("#"); // fn_proj = "d=f"
     * @endcode
     */
    FileName afterLastOf(const String& str) const;

    /** Add string at the beginning
     *
     * If there is a path then the prefix is added after the path.
     *
     * @code
     * fn_proj = "imgs/g1ta00001";
     * fn_proj.add_prefix("h"); // fn_proj == "imgs/hg1ta00001"
     *
     * fn_proj = "g1ta00001";
     * fn_proj.add_prefix("h"); // fn_proj == "hg1ta00001"
     * @endcode
     */
    FileName addPrefix(const String& prefix) const;

    /** Add extension at the end.
     *
     * The "." is added. If teh input extension is "" then the same name is
     * returned, with nothing added.
     *
     * @code
     * fn_proj = "g1ta00001";
     * fn_proj.add_extension("xmp"); // fn_proj == "g1ta00001.xmp"
     * @endcode
     */
    FileName addExtension(const String& ext) const;

    /** Remove last extension, if any
     *
     * @code
     * fn_proj = "g1ta00001.xmp";
     * fn_proj = fn_proj.without_extension(); // fn_proj == "g1ta00001"
     *
     * fn_proj = "g1ta00001";
     * fn_proj = fn_proj.without_extension(); // fn_proj == "g1ta00001"
     * @endcode
     */
    FileName withoutExtension() const;

    /** Remove the root
     *
     * @code
     * fn_proj = "g1ta00001.xmp";
     * fn_proj = fn_proj.without_root(); // fn_proj == "00001.xmp"
     *
     * fn_proj = "g1ta00001";
     * fn_proj = fn_proj.without_root(); // fn_proj == "00001"
     * @endcode
     */
    FileName withoutRoot() const;

    /** Insert before first extension
     *
     * If there is no extension, the insertion is performed at the end.
     *
     * @code
     * fn_proj = "g1ta00001.xmp";
     * fn_proj = fn_proj.insert_before_extension("pp");
     * // fn_proj == "g1ta00001pp.xmp"
     *
     * fn_proj = "g1ta00001";
     * fn_proj = fn_proj.insert_before_extension("pp");
     * // fn_proj=="g1ta00001pp"
     * @endcode
     */
    FileName insertBeforeExtension(const String& str) const;

    /** Remove a certain extension
     *
     * It doesn't matter if there are several extensions and the one to be
     * removed is in the middle. If the given extension is not present in the
     * filename nothing is done.
     *
     * @code
     * fn_proj = "g1ta00001.xmp.bak";
     * fn_proj = fn_proj.remove_extension("xmp");
     * // fn_proj == "g1ta00001.bak"
     * @endcode
     */
    FileName removeExtension(const String& ext) const;

    FileName removeLastExtension() const;


    /**Extract the directory portion from a filename string
     * (e.g. /foo/bar.txt -> /foo).
     *
     */
    FileName removeFilename() const;

    /** Remove all extensions
     */
    FileName removeAllExtensions() const;

    /** Remove file format
     * @code
     * fn_proj = "g1ta00001.xmp";
     * fn_proj = fn_proj.get_file_format(); // fn_proj == "xmp"
     * fn_proj = "g1ta00001.nor:spi";
     * fn_proj = fn_proj.get_file_format(); // fn_proj == "spi"
     * fn_proj = "input.file#d=f#x=120,120,55#h=1024";
     * fn_proj = fn_proj.get_file_format(); // fn_proj == "raw"
     * @endcode
     */
    FileName removeFileFormat() const;

    /** Is this file a MetaData file?
     * Returns true if the get_file_format extension == "sel", "doc" or "xmd"
     * Otherwise, the file is opened and checked for the occurence of "XMIPP_3 *" in its first line
     */
    bool isMetaData(bool failIfNotExists=true) const;

    /** True if it is a Star 1 file. */
    bool isStar1(bool failIfNotExists) const;

    /** Replace a the first ocurrence of substring 'subOld' by 'subNew'
     * If 'subOld' not present, not changes are done.
     * The search is done from left to right
     *
     * @code
     * fn_proj = "g1ta00001.xmp.bak";
     * fn_proj = fn_proj.replace(".xmp", ".bor");
     * // fn_proj == "g1ta00001.bor.bak"
     *
     * fn_proj = "g1ta00001.xmp.bak";
     * fn_proj = fn_proj.replace(".bor", "");
     * // fn_proj == "g1ta00001.bak"
     * @endcode
     * */
    FileName replaceSubstring(const String& subOld, const String& subNew) const;

    /** Substitute last extension by a new one
     *
     * @code
     * fn_proj = "g1ta00001.xmp.bak";
     * fn_proj = fn_proj.replaceExtension("bor");
     * // fn_proj == "g1ta00001.xmp.bor"
     * @endcode
     */
    FileName replaceExtension(const String& newExt) const;

    /** Remove a substring from the filename
     *
     * If the substring is not present the same FileName is returned, if it is
     * there the substring is removed.
     * It is equivalent to replaceSubstring(sub, "");
     */
    FileName removeSubstring(const String& sub) const;

    /** Remove until prefix
     *
     * Remove the starting string until the given prefix, inclusively. For
     * instance /usr/local/data/ctf-image00001.fft with ctf- yields
     * image00001.fft. If the prefix is not found nothing is done.
     */
    FileName removeUntilPrefix(const String& prefix) const;

    /** Remove all directories
     *
     * Or if keep>0, then keep the lowest keep directories
     */
    FileName removeDirectories(int keep = 0) const;
    /**copy one file
     *
     * s
     */
    void copyFile(const FileName & target) const;

    // This funtion is in funcs, cannot be here because need metadata_generic and metadata_generic
    // need filename
    //    /* Copy one image
    //     *
    //     */
    //    void copyImage(const FileName & target) const;

    /** True if the filename exists, check from current directory if is relative
     *
     * @code
     * FileName fn("g1ta0001");
     * if (fn.exists())
     *     std::cout << "The file exists" << std::endl;
     * @endcode
     */
    bool exists() const;

    /** Delete the file if exists */
    void deleteFile() const;

    /** True if the file exists in the current directory
     *  Remove leading xx@ and tailing :xx
     *
     * @code
     * if (exists("g1ta00001"))
     *     std::cout << "The file exists" << std::endl;
     * @endcode
     */
    bool existsTrim() const;

    /** True if the path is a directory */
    bool isDir() const;

    /** Return the list of files if this filename is a directory*/
    void getFiles(std::vector<FileName> &files) const;

    /** Waits until the given filename has a stable size
     *
     * The stable size is defined as having the same size within two samples
     * separated by time_step (microsecs).
     *
     * An exception is throw if the file exists but its size cannot be stated.
     */
    void waitUntilStableSize(size_t time_step = 250000);

    /** Write a zero filled file with the desired size.
     *
     * The file is written by blocks to speed up, you can modify the block size.
     * An exception is thrown if any error happens
     */
    void createEmptyFile(size_t size, size_t block_size = 102400);

    /** Create empty file with a given length. */
    void createEmptyFileWithGivenLength(size_t length = 0) const;

    /**
     * makePath - ensure all directories in path exist
     * Algorithm takes the pessimistic view and works top-down to ensure
     * each directory in path exists, rather than optimistically creating
     * the last element and working backwards. Return null if fails
     */
    int makePath(mode_t mode = 0755) const;

    /** This function raised an ERROR if the filename if not empty and if
     * the corresponding file does not exist.
     * This may be useful to have a better (killing) control on (mpi-regulated) jobs
     *
     * @code
     *   exit_if_not_exists("control_file.txt");
     * @endcode
     *
     * This function is not ported to Python.
     */
    void assertExists();
    //@}
};

/** Returns the base directory of the Xmipp installation
 */
char* getXmippPath();

/** Check if the file exists using the stat function
 */
inline bool fileExists( const char *filename )
{
    Stat buffer ;
    // do not use stat since updates the atime flag. This may be dangerous
    // for mpi processes in which the same file is accessed many times
    // from different computers even if not writting is made 
    //return ( stat(filename, &buffer) == 0 ? true : false );
    return ( access (filename, F_OK)  == 0 ? true : false);

}

/** Check if the file exists using the stat function
 */
inline bool fileExists( const FileName &filename )
{
    return fileExists(filename.c_str());
}

/** Copy one image
 */
void copyImage(const FileName & source,const FileName & target);

//Functions to make easy to change the magicWord used to identify the
//file format of metadata
extern String FileNameVersion;
void setMetadataVersion(String version);
String getMetadataVersion(void);

/** This class is used for comparing filenames.
 *
 * Example: "g0ta00001.xmp" is less than "g0ta00002.xmp"
 *
 * This class is needed to define a std::map as
 * map<FileName,FileName,FileNameComparison> myMap;
 *
 * This function is not ported to Python.
 */
class FileNameComparison
{
public:
    inline bool operator ()(const FileName &fn1, const FileName &fn2)
    {
        return fn1<fn2;
    }
};


/** Lock/unlock files
 *
 */
class FileLock
{
protected:


#ifndef __MINGW32__
    struct flock fl;
#endif

    bool islocked;
    int  filenum;

public:

    /// constructor
    FileLock()
    {
#ifndef __MINGW32__
        fl.l_type   = F_WRLCK;  /* F_RDLCK, F_WRLCK, F_UNLCK    */
        fl.l_whence = SEEK_SET; /* SEEK_SET, SEEK_CUR, SEEK_END */
        fl.l_start  = 0;        /* Offset from l_whence         */
        fl.l_len    = 0;        /* length, 0 = to EOF           */
        fl.l_pid    = getpid(); /* our PID                      */

        islocked = false;
        filenum = NULL;
#endif
    }

    FileLock(int fileno)
    {
        FileLock();
        this->filenum = fileno;
    }

    /// Lock file
    void lock(int fileno = NULL);
    void lock(FILE* file);

    ///Unlock
    void unlock();

};

//@}

#endif /* FILENAME_H_ */
