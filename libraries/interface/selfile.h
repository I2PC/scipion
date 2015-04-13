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
#ifndef SELFILE_H
#define SELFILE_H

#include <string>
#include <vector>

#include <data/xmipp_funcs.h>
#include <data/xmipp_image.h>

// Forward declaration
class SelFile;

/// @defgroup SelFiles Selection files (Selfiles)
/// @ingroup DataLibrary
//@{
/** Line of a selection file.
 *
 * The selection file is a collection of (STL list) of Selection Lines. This
 * class needn't be accessed in common programs since the SelFile class offers
 * most of the used functions. However, the class is shown in case you may need
 * to access specifically to any of the Selection Line functions alone. The
 * SelLines can be either comments or data line, they might be of other
 * internal used types (not assigned or to be discarded), but you may assume
 * that if a line is not a comment nor a data line then it can be skipped.
 */
class SelLine
{
public:
    typedef enum
    {
        DISCARDED = -1,
        ACTIVE = 1
    } Label;

    typedef enum
    {
        NOT_CONSIDERED = -1,
        NOT_ASSIGNED = 0,
        DATALINE = 1,
        COMMENT = 2
    } Line_Type;

private:
    Line_Type line_type;
    std::string text;
    Label label;
    int number;

    friend class SelFile;

public:
    /** Empty Constructor.
     *
     * The selection line is created with no type (neither comment or data).
     * You must use the function @p set_type to assign a type.
     */
    SelLine() : line_type(NOT_ASSIGNED)
    {}

    /** Copy constructor.
     */
    SelLine(const SelLine& line);

    /** Assignment.
     */
    SelLine& operator=(const SelLine& line);

    /** Another function for assigment.
     */
    void assign(const SelLine& line);

    /** Get text of this line.
     */
    const std::string& get_text()
    {
        return text;
    }

    /** Get label of this line.
     */
    short get_label()
    {
        return label;
    }

    /** Get number of this line.
     */
    int get_number()
    {
        return number;
    }

    /** True if current line is a comment.
     */
    int Is_comment()
    {
        return line_type == COMMENT;
    }

    /** True if current line is data.
     */
    int Is_data()
    {
        return line_type == DATALINE;
    }

    /** Set type.
     * Only the comment flag is set, the possible data is not lost. The comment
     * text is not touched. The valid types are DATALINE, COMMENT, NOT_ASSIGNED,
     * and NOT_CONSIDERED.
     */
    void set_type(Line_Type type)
    {
        line_type = type;
    }

    /** Set number of this line.
     */
    void set_number(int n)
    {
        number = n;
    }

    /** Lesser than.
     * l1 is lesser than l2 if line_type(l1)<line_type(l2) or if they are equal
     * if the text of l1 is lesser than the text of l2. The order of the line
     * types are NOT_CONSIDERED, NOT_ASSIGNED, DATALINE, COMMENT
     */
    friend bool operator<(const SelLine& l1, const SelLine& l2);

    friend SelFile compare(SelFile& SF1, SelFile& SF2, const int mode);
    friend std::vector< SelLine >::iterator find(std::vector< SelLine >& text,
            const std::string& img_name);

    /** Show a SelLine.
     */
    friend std::ostream& operator<<(std::ostream& o, const SelLine& line);

    /** Read a Selection Line.
     *
     * An exception is thrown if the line doesn't meet the SelFile
     * specifications.
     */
    friend std::istream& operator>>(std::istream& i, SelLine& lin);
};

/** Selection File
 *
 * The SelFile is an object which keeps in memory all the information
 * associated to a @p .sel file.
 */
class SelFile
{
    // Filename of the .sel
    FileName fn_sel;

    // The whole file is inside
    std::vector< SelLine > text_line;

    // Number of valid images
    int no_imgs;

    // Pointer to current line
    std::vector< SelLine >::iterator current_line;

    // Find a specific image name in a list of sellines
    std::vector< SelLine >::iterator find(const std::string& img_name);

    // Move the current_line to the next line with this label
    void adjust_to_label(SelLine::Label label);
public:
    /// @name Selfile constructors
    //@{

    /** Empty constructor.
     *
     * There is no file associated yet.
     *
     * @code
     * SelFile sel;
     * @endcode
     */
    SelFile();

    /** Constructor with filename, read from disk.
     *
     * The given name is loaded (method @p read) as a selection file.
     *
     * @code
     * SelFile sel("g1t.sel");
     * @endcode
     */
    SelFile(const FileName &name)
    {
        read(name);
    }

    /** Copy constructor.
     *
     * @code
     * SelFile sel2(sel1);
     * @endcode
     */
    SelFile(const SelFile& sel);

    /** Reserve memory for N entries.
     *
     * It doesn't matter if entries are comments or data. The current line is
     * set to the beginning of the SelFile
     */
    void reserve(int n)
    {
        text_line.reserve(n);
        current_line = text_line.begin();
    }

    /** Empties the object.
     *
     * @code
     * sel.clear();
     * @endcode
     */
    void clear();
    //@}

    /// @name Some operators
    //@{

    /** Assignment.
     *
     * @code
     * SelFile sel2 = sel1;
     * @endcode
     */
    SelFile& operator=(const SelFile& sel);

    /** Another function for assigment.
     */
    void assign(const SelFile& sel);

    /** Show a selection file.
     *
     * Shows all the lines either they are comments, active images or discarded
     * images. A new line is printed at the end.
     *
     * @code
     * std::cout << sel;
     * @endcode
     */
    friend std::ostream& operator<<(std::ostream& o, const SelFile& sel);
    //@}

    /// @name Managing files in disk
    //@{

    /** Read a file from disk.
     *
     * The old information on the variable is overwritten. An exception is
     * thrown if the file doesn't exist. Lines which do not fit the comment
     * structure or the "image-label" structure are ignored. The image name is
     * limited to MAX_FILENAME_LENGTH characters. After reading the selfile
     * pointer is moved to the first ACTIVE image.
     *
     * @code
     * sel.read("g2t.sel");
     * @endcode
     */
    void read(const FileName& name, int overrinding = 1);

    /** Append a file from disk to an already read one.
     *
     * The old information on the variable is not lost. All lines in the
     * selection file to be read are appened at the end of the already read one
     * without any kind of check.
     *
     * @code
     * sel.read("g1t.sel");
     * sel.append("g2t.sel");
     * @endcode
     */
    void append(const FileName& name)
    {
        read(name, 0);
    }

    /** Merge a file from disk with an already read one.
     *
     * All lines (except the comments) in the selection file to be read are
     * either added at the end of the already read one if they are not present
     * at them, either ignored if they are already present, or marked with a
     * comment if the corresponding image name is present in both files but
     * with different labels (active, discarded), in this case the image
     * remains active but a comment in the preceeding line informs you of the
     * situation.
     *
     * @code
     * sel.read("g1t.sel");
     * sel.merge("g2t.sel");
     * @endcode
     */
    void merge(const FileName& name);

    /** Merge this file with another selfile.
     */
    void merge(SelFile& sel);

    /** Merge two already read files.
     *
     * @code
     * sel1.read("g1t.sel");
     * sel2.merge("g2t.sel");
     * sel1 = sel1 + sel2;
     * @endcode
     */
    SelFile operator+(SelFile& sel);

    /** Split this file into two random halves.
     */
    void split_in_two(SelFile& sel1, SelFile& sel2);

    /** Split this file into N random parts of approximately the same size.
     */
    void split_in_N(int N, std::vector< SelFile >& parts);

    /** For MPI-parallelized runs: select relevant part of selfile for that
     * rank.
     */
    void mpi_select_part(int rank, int size, int& num_img_tot);

    /** For MPI-parallelized runs: select relevant part of selfile for that
     * rank. 
     * I just do not undertand the previous function
     */
    void mpi_select_part2(int rank, int size, int& num_img_tot,int mpi_job_size);

    /** Choose a subset of images.
     */
    void chooseSubset(int firstImage, int lastImage, SelFile &SFsubset);

    /** Write a selection file to disk.
     *
     * If you give a name then it becomes like a "Save as ..." and from this
     * point on the name of the selection file has changed.
     *
     * @code
     * sel.write(); // Save
     * sel.write("g3t.sel"); // Save as ...
     * @endcode
     */
    void write(const FileName& sel_name = "");
    //@}

    /// @name Moving the current line pointer
    //@{

    /** Go to the beginning of the file.
     *
     * Moves the pointer to the first line of the file either it is a comment,
     * an active image or a discarded one.
     *
     * @code
     * sel.go_beginning();
     * @endcode
     */
    void go_beginning()
    {
        current_line = text_line.begin();
    }

    /** Go to the first ACTIVE image.
     *
     * Moves the pointer to the first active image in the file.
     *
     * @code
     * sel.go_first_ACTIVE();
     * @endcode
     */
    void go_first_ACTIVE()
    {
        go_beginning();
        adjust_to_label(SelLine::ACTIVE);
    }

    /** Get the line i in the file
     *
     * Get the line i in the file for reading. If it is a comment
     * the comment text is returned. If it is an image, the image filename.
     * The first line in the file is number 0.
     *
     * @code
     * std::cout << sel(0) << std::endl;
     * @endcode
     */
    const std::string& operator()(int i)
    {
        return text_line[i].text;
    }

    /** Returns the name of the next image with a certain label.
     *
     * The default label is ACTIVE, ie, by default this function returns the
     * name of the next ACTIVE image. But you can give as label DISCARDED and
     * the function will return the name of the next DISCARDED image starting
     * at the current position of the current_line "pointer". If the file is at
     * the end of the selection file, "" is returned. After this function the
     * "pointer" is actually pointing to the next line following the returned
     * image name.
     *
     * @code
     * name = sel.NextImg();  // Next active image
     * name = sel.NextImg(SelLine::DISCARDED); // Next discarded image
     * @endcode
     */
    const std::string& NextImg(SelLine::Label label = SelLine::ACTIVE);

    /** Move the current pointer to the next image, disregarding its label.
     *
     * It doesn't matter if next image is ACTIVE or DISCARDED, this function
     * moves the current pointer to it.
     *
     * @code
     * sel.go_beginning();
     * while (!sel.eof())
     * {
     *     std::cout << sel.current();
     *     sel.next();
     * }
     * @endcode
     */
    void next()
    {
        current_line++;
    }

    /** Jump over a number of lines with a given label.
     *
     * Starting from the current_line "pointer" this function skips a given
     * number of entries with a certain label. For instance, jump over 1 active
     * image is to jump to the next active image. Jump over 2 active images is
     * to jump to the next of the next active image, and so on. You can give as
     * label DISCARDED, too. The number of images to jump must always be
     * positive, the jump cannot be done backwards.
     *
     * @code
     * sel.jump(2); // Jump over 2 active images
     * sel.jump(2, SelLine::ACTIVE); // The same
     * sel.jump(2, SelLine::DISCARDED) // Jump over 2 discarded images
     * @endcode
     */
    void jump(int count, SelLine::Label label = SelLine::ACTIVE);

    /** Jump over a number of lines disregarding the label.
     *
     * Returns false if the end of the line is reached before the number of
     * requested line jumps
     */
    bool jump_lines(int count);

    /** Move "pointer" to a certain image filename.
     *
     * This function searches for an image name within the file, and locate
     * the current line "pointer" pointing to that line. If the image name is
     * not present (it is not the same "not present" and "discarded") in the
     * selection file, then the pointer is moved to the end of the selection
     * file. You can check this situation using eof(). It doesn't matter if the
     * current line "pointer" before the function call is after the line where
     * the image name is, this function makes a search all over the file,
     * regardless the previous situation of the current line "pointer".
     *
     * @code
     * sel.search("g1ta0001");
     * @endcode
     */
    void search(const std::string& img_name)
    {
        current_line = find(img_name);
    }

    /** True if current line "pointer" is at the end of file.
     *
     * @code
     * if (sel.eof())
     *     std::cout << "The selection file is over\n";
     * @endcode
     */
    int eof()
    {
        return current_line == text_line.end();
    }
    //@}

    /// @name Getting information
    //@{

    /** Returns the name of the file.
     */
    FileName name() const
    {
        return fn_sel;
    }

    /** True if current line is a data line and it is active.
     */
    int Is_ACTIVE() const
    {
        return current_line->Is_data() &&
               current_line->get_label() == SelLine::ACTIVE;
    }

    /** True if current line is a data line and it is active.
     */
    int Is_DISCARDED() const
    {
        return current_line->Is_data() &&
               current_line->get_label() == SelLine::DISCARDED;
    }

    /** True if current line is a comment.
     */
    int Is_COMMENT() const
    {
        return current_line->Is_comment();
    }

    /** Returns current line as a Sel Line.
     */
    const SelLine& current()
    {
        return *current_line;
    }

    /** Another function to get a SelLine.
     */
    void get_current(SelLine& _SL)
    {
        _SL = current();
    }

    /** True if the image name is inside the selection file.
     *
     * The current line "pointer" is not modified. If an image is discarded in
     * the selection file, this function still will say that it exists,
     * although it is discarded.
     *
     * @code
     * if (sel.exists("g1ta0001"))
     *     std::cout << "g1ta0001 exists in the selection file\n";
     * @endcode
     */
    int exists(const std::string& img_name)
    {
        return find(img_name) != text_line.end();
    }

    /** Number of images inside a selection file with a certain label.
     *
     * This function returns the number of images inside the selection file
     * with a given label. By default this label is ACTIVE.
     *
     * @code
     * std::cout << "There are " << sel.ImgNo() << " active images\n";
     * std::cout << "There are " << sel.ImgNo(SelLine::ACTIVE)
     *     << " active images\n";
     * std::cout << "There are " << sel.ImgNo(SelLine::DISCARDED)
     *     << " discarded images\n";
     * @endcode
     */
    int ImgNo(SelLine::Label label = SelLine::ACTIVE) const;

    /** Returns the number of lines within a file.
     *
     * This function gives the total number of lines (including comments)
     * within a file.
     *
     * @code
     * std::cout << "There are " << sel.LineNo() << " lines in this file\n";
     * @endcode
     */
    int LineNo();

    /** Returns the size of the images inside.
     *
     * The filenames within a selection file are supposed to be for SPIDER
     * images, this function opens one of the images (an active one) and
     * returns the size of that image, supposed to be the same for the rest of
     * the images in the selection file.
     *
     * An exception is thrown if the first valid image in the selfile, doesn't
     * exist in the disk or it is not a XMIPP image.
     *
     * @code
     * sel.ImgSize(y, x);
     * @endcode
     */
    void ImgSize(int& Ydim, int& Xdim);

    /** Returns the extension of the files inside.
     *
     * This function returns the extension of the first active file.
     */
    FileName FileExtension();

    /** Returns the maximum length of an active filename inside the selfile.
     *
     * The current pointer is not moved.
     */
    int MaxFileNameLength();

    /** Get the filename of the current line.
     *
     * If the current line "pointer" is at the end of the file or is pointing
     * to a comment then an empty string is returned.
     *
     * @code
     * fn = sel.get_current_file();
     * @endcode
     */
    const std::string get_current_file();

    /** Get current line.
     */
    SelLine get_current_line()
    {
        return *current_line;
    }

    /** Get the filename at the ACTIVE line number i.
     *
     * The first file is number 0. If i is greater than the total number
     * of ACTIVE files, then "" is returned.
     *
     * @code
     * fn = sel.get_file_number(i);
     * @endcode
     */
    const std::string get_file_number(int i);
    //@}

    /// @name Modifying the selection file
    //@{

    /** Removes an image from the selection file.
     *
     * This function searches for an image in the selection file, if it is
     * found then the corresponding line is deleted. If the image is actually
     * being pointed by the current line, then the current line is now the
     * following line.
     *
     * @code
     * sel.remove("g1ta0001");
     * @endcode
     */
    void remove(const std::string& img_name);

    /** Removes actual line.
     *
     * This function removes the current line, either it is a comment or an
     * image. The current line "pointer" is moved to the following line in the
     * file.
     *
     * @code
     * sel.remove_current();
     * @endcode
     */
    void remove_current();

    /** Set label of an image.
     *
     * This function searches for an image inside the selection file and sets
     * its label to the given label. If the image is not found in the file,
     * then it is added at the end with the given label. The current line
     * pointer is not modified.
     *
     * @code
     * sel.set("g1ta0001", SelLine::ACTIVE);
     * @endcode
     */
    void set(const std::string& img_name, SelLine::Label label);

    /** Set the label of the current file.
     *
     * The same as the previous function but the label is set to the file
     * currently pointed.
     */
    void set_current(SelLine::Label label);

    /** Change current filename.
     *
     * This function changes the current filename to a new one if it is not a
     * comment. If it is a comment line, nothing is done.
     */
    void set_current_filename(const FileName& fn_new);

    /** Insert image before current line.
     *
     * There is no checking for the previous existence of the img. The current
     * line is still pointing to the same line as it was before entering the
     * function.
     *
     * @code
     * sel.insert("g1ta0000");
     * sel.insert("g1ta0000", SelLine::DISCARDED);
     * @endcode
     */
    void insert(const std::string& img_name, SelLine::Label label = SelLine::ACTIVE);

    /** Insert line before current line.
     *
     * It is checked that the line is either a comment or data line, in this
     * case that the label is right, too. The current line is still pointing
     * to the same line as it was before entering the function.
     *
     * @code
     * sel.insert("g1ta0000", SelLine::ACTIVE);
     * @endcode
     */
    void insert(const SelLine& _selline);

    /** Insert a comment before the current line.
     *
     * Comments must not start with any special character since a "#" is
     * automatically added at the beginning of the line. The current line is
     * still pointing to the same line as it was before entering the function.
     *
     * @code
     * sel.insert_comment("This is a comment");
     * @endcode
     */
    void insert_comment(const std::string& comment);

    /** Deletes all DISCARDED images from the selection file.
     *
     * The current line "pointer" is moved to the beginning of the file.
     *
     * @code
     * sel.clean();
     * @endcode
     */
    void clean();

    /** Deletes all comments from the selection file.
     *
     * The current line "pointer" is moved to the beginning of the file.
     *
     * @code
     * sel.clean_comments();
     * @endcode
     */
    void clean_comments();
    //@}

    /// @name Helpful procedures
    //@{

    /** Sort images in ascending order.
     *
     * All images are sorted in ascending order either they are active or
     * discarded. All comments are gathered at the end of the resulting
     * selection file. The current line of the resulting selection file is
     * placed at the beginning of the file.
     *
     * @code
     * sel2 = sel1.sort_by_filenames();
     * @endcode
     */
    SelFile sort_by_filenames();

    /** Alter order in sel file.
     *
     * A new selection file with all the images of the actual object (either
     * they are active or discarded) is created, but this time all images are
     * in a random order. The comments of the original selection file are lost
     * in the new copy. The current line of the resulting selection file is
     * placed at the beginning of the file.
     *
     * @code
     * sel2 = sel1.randomize();
     * @endcode
     */
    SelFile randomize();

    /** Random subset.
     *
     * A new selection file is created with the number of images specified in
     * subsetN. The new subset can be created with or without replacement.
     */
    SelFile randomSubset(int subsetN, bool withReplacement=true);

    /** Discard randomly N images.
     *
     * A set of N images are discarded from the actual selection file. If N is
     * equal or greater than the actual number of images within the file, all
     * images are discarded. Comments are kept at their original positions. The
     * current line of the resulting selection file is placed at the beginning
     * of the file.
     *
     * @code
     * sel2 = sel1.random_discard(3);
     * @endcode
     */
    SelFile random_discard(int N);

    /** Compare two selection files.
     *
     * The result is another selection file. At the beginning of it there is
     * information about the number of active and discarded images on both
     * input selection files, about the number of matching files (a file is
     * said to match if it is active in both selection files), the number of
     * active files which are only in the first selection file, and the number
     * of active files which are only in the second. Then goes the list of
     * matching files, the list of files only in SF1 and the list of files
     * only in SF2. There are comments enough to know where things start and
     * finish, and what the numbers are at the beginning. If a file is active
     * in a file and discarded in the other, then it is said to match and it is
     * kept as active, a preceeding comment warns of this situation.
     *
     * if mode<0, output file will contain all information
     * if mode=0, output will be selfile with overlapping images
     * if mode=1, output will be selfile with images only in file 1
     * if mode=2, output will be selfile with images only in file 2
     *
     * @code
     * sel3 = compare(sel1, sel2);
     * @endcode
     */
    friend SelFile compare(SelFile& SF1, SelFile& SF2, const int mode);
    //@}
};
//@}
#endif
