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
/*****************************************************************************/
/* SEL FILES                                                                 */
/*****************************************************************************/
#ifndef _XMIPPSELFILES_HH
   #define _XMIPPSELFILES_HH

#include <string>
#include <vector>
#include "xmippFuncs.hh"
#include "xmippImages.hh"

// Possible labels for images inside a .Sel file
// Forward declarations ====================================================
class SelFile;

/**@name Selection Files */
//@{

/*****************************************************************************/
/* LINE OF A SEL FILE                                                        */
/*****************************************************************************/
/** Line of a selection file.
    The selection file is a collection of (STL list) of Selection
    Lines. This class needn't be accessed in common programs since the
    SelFile class offers most of the used functions. However, the
    class is shown in case you may need to access specifically to any
    of the Selection Line functions alone. The SelLines can be either
    comments or data line, they might be of other internal used types
    (not assigned or to be discarded), but you may assume that if a
    line is not a comment nor a data line then it can be skipped.
*/
class SelLine {
public:
   typedef enum {DISCARDED=-1, ACTIVE=1} Label;
   typedef enum {NOT_CONSIDERED=-1, NOT_ASSIGNED=0, DATALINE=1, COMMENT=2}
      Line_Type;

private:
   Line_Type line_type;
   string text;                           // if a comment, the whole line
                                          // if a file, the filename only
   Label label;                           // if a file, the label associated

   friend class SelFile;
public:
   /**@name Constructors */
   //@{
   /** Empty Constructor.
       The selection line is created with no type (neither comment or data).
       You must use the function \Ref{set_type} to assign a type. */
   SelLine(): line_type(NOT_ASSIGNED) {};
   
   /** Copy constructor. */
   SelLine(const SelLine &l);
   
   /** Assignment. */
   SelLine& operator = (const SelLine &SL);
   //@}
   
   /**@name Structure information */
   //@{
   /// Get text of this line
   string get_text() {return text;}
   
   /// Get label of this line
   short get_label() {return label;}
   
   /** True if current line is a comment. */
   int Is_comment() {return line_type==COMMENT;}

   /** True if current line is a comment. */
   int Is_data() {return line_type==DATALINE;}

   /** Set type.
       Only the comment flag is set, the possible data is not lost.
       The comment text is not touched. The valid types are DATALINE,
       COMMENT, NOT_ASSIGNED, and NOT_CONSIDERED. */
   void set_type(Line_Type _line_type) {line_type=_line_type;}
   //@}

   /**@name Useful operations */
   //@{
   /** l1<l2.
       l1 is lesser than l2 if line_type(l1)<line_type(l2) or
       if they are equal if the text of l1 is lesser than the text
       of l2. The order of the line types are NOT_CONSIDERED,
       NOT_ASSIGNED, DATALINE, COMMENT */
   friend int operator < (const SelLine &l1, const SelLine &l2);

   // Compare two .sel files
   friend SelFile compare(SelFile &SF1, SelFile &SF2);
   
   // Find an image inside a list
   friend vector<SelLine>::iterator find(vector<SelLine> &text,
      string &img_name);
   //@}
   
   /**@name I/O */
   //@{
   /** Show a Selection Line */
   friend ostream& operator << (ostream& o, SelLine &SFL);
   
   /** Read a Selection Line.
       An exception is thrown if the line doesn't meet the Selection
       File specifications. */
   friend istream& operator >> (istream& o, SelLine &SFL) _THROW;
   //@}
};

/*****************************************************************************/
/* SEL FILE                                                                  */
/*****************************************************************************/
/** Selection Files.
   The sel file is an object which keeps in memory all the information
   associated to a .sel file (see \Ref{Selection Files} for detailed information
   about the file structure). With respect to the object structure you can
   think of the selection file as a file loaded in memory with an internal
   "pointer" (in fact is an STL iterator) pointing to
   the current line considered at this moment. You can move forward
   this pointer or set it to the beginning of the file, jump over any
   number of images with a certain label. There are only two valid labels:
   ACTIVE and DISCARDED.
*/   
class SelFile {
   // Structure ------------------------------------------------------------
   FileName                  fn_sel;        // Filename of the .sel
   vector<SelLine>           text_line;     // The whole file is inside
   int                       no_imgs;       // Number of valid images inside
   vector<SelLine>::iterator current_line;  // "pointer" to current line
                                            // within file

   // Private procedures ---------------------------------------------------
   vector<SelLine>::iterator find(        // "pointer" to that image inside
      string img_name);                   // .sel file. If the image is not
                                          // there then a "pointer" to past-
                                          // last element is returned
   void adjust_to_label                   // given any current_line position
      (SelLine::Label label);             // this function moves until
                                          // the next entry with the label
                                          // selected. If label==0 then to
                                          // the next not DISCARDED

public:
   // Procedures -----------------------------------------------------------
   /**@name Constructors*/
   //@{
   /** Empty constructor.
       There is no file associated yet.
       \\Ex: SelFile SF; */
   SelFile();
   
   /** Constructor with filename, read from disk.
       The given name is loaded (method read) as a selection file.
       \\Ex: SelFile SF("g1t.sel");
       @see read*/
   SelFile(FileName sel_name) {read(sel_name);}
   
   /** Copy constructor.
       \\Ex: SelFile SF2(SF1); */
   SelFile(const SelFile &SF);

   /** Reserve memory for N entries.
       It doesn't matter if entries are comments or data. It
       is very important to reserve memory if you don't want a
       memory explosion. The current line is set to the beginning of the
       selfile*/
   void reserve(int N) {text_line.reserve(N); current_line=text_line.begin();}

   /** Empties the object.
       \\Ex: SF.clear(); */
   void clear();
   //@}

   // Some basic operators ..................................................
   /**@name Some operators*/
   //@{
   /** Assignment.
       \\ Ex: SF2=SF1; */
   SelFile& operator =(const SelFile &SF);
   
   /** Show a selection file.
       Shows all the lines either they are comments, active images or
       discarded images. A new line is printed at the end.
       \\ Ex: cout << SF; */
   friend ostream& operator << (ostream& o, SelFile &SF);
   //@}

   // Managing files from disk .............................................
   /**@name Managing files in disk*/
   //@{
   /** Read a file from disk.
       The old information on the variable is overwritten. An exception
       is thrown if the file doesn't exist. Lines which do not fit the
       comment structure or the "image-label" structure are ignored.
       The image name is limited to MAX_FILENAME_LENGTH characters.
       After reading the selfile pointer is moved to the first ACTIVE
       image. 
       \\Ex: SF.read("g2t.sel"); */
   void read(const FileName &sel_name, int overrinding=1) _THROW;
   
   /** Append a file from disk to an already read one.
       The old information on the variable is not lost. All lines in
       the selection file to be read are appened at the end of the
       already read one without any kind of check.
       \\ Ex: SF.read("g1t.sel"); SF.append("g2t.sel"); */
   void append(const FileName &sel_name) {read(sel_name,0);}
   
   /** Merge a file from disk with an already read one.
       All lines (except the comments) in the selection file to be read
       are either added at the end of the already read one if they
       are not present at them, either ignored if they are already
       present, or marked with a comment if the corresponding
       image name is present in both files but with different labels
       (active, discarded), in this case the image remain active
       but a comment in the preceeding line informs you of the situation.
       \\Ex: SF.read("g1t.sel"); SF.merge("g2t.sel"); */
   void merge(const FileName &sel_name);
   
   /** Merge this file with another selfile. */
   void merge(SelFile &SF);

   /** Merge two already read files.
       \\ Ex: SF1.read("g1t.sel"); SF2.merge("g2t.sel"); SF1=SF1+SF2; */
   SelFile  operator +(SelFile &SF);
   
   /** Write a selection file to disk.
       If you give a name then it becomes like a "Save as ..." and from
       this point on the name of the selection file has changed.
       \\Ex: SF.write(); ---> Save
       \\Ex: SF.write("g3t.sel") ---> Save as ... */
   void write(const FileName &sel_name="") _THROW;
   //@}

   // Moving the current_line "pointer" ....................................
   /**@name Moving the current line "pointer"
      See the get information section to know how to extract the information
      from the current selected line */
   //@{
   /** Go to the beginning of the file.
       Moves the pointer to the first line of the file either it is
       a comment, an active image or a discarded one.
       \\Ex: SF.go_beginning(); */
   void go_beginning() {current_line=text_line.begin();}
   
   /** Go to the first ACTIVE image.
       Moves the pointer to the first active image in the file.
       \\Ex: SF.go_first_ACTIVE(); */
   void go_first_ACTIVE() {go_beginning(); adjust_to_label(SelLine::ACTIVE);}
   
   /** Returns the name of the next image with a certain label.
       The default label is ACTIVE, ie, by default this function returns
       the name of the next ACTIVE image. But you can give as label
       DISCARDED and the fucntion will return the name of the next DISCARDED
       image starting at the current position of the current_line "pointer".
       If the file is at the end of the selection file, "" is returned.
       After this function the "pointer" is actually pointing to the next
       line following the returned image name.
       \\ Ex: fn=SF.NextImg(); ---> Next active image
       \\ Ex: fn=SF.NextImg(SelLine::DISCARDED); ---> Next discarded image */
   string NextImg(SelLine::Label label=SelLine::ACTIVE);
   
   /** Move the current pointer to the next image, disregarding its label.
       It doesn't matter if next image is ACTIVE or DISCARDED, this function
       moves the current pointer to it.
       \\ Ex: 
       \begin{verbatim}
       SF.go_beginning();
       while (!SF.eof()) {
          cout << SF.current();
          SF.next();
       }
       \end{verbatim} */
   void next() {current_line++;}

   /** Jump over a number of lines with a given label.
       Starting from the current_line "pointer" this function skips 
       a given number of entries with a certain label. For instance,
       jump over 1 active image is to jump to the next active image.
       Jump over 2 active images is to jump to the next of the next
       active image, and so on. You can give as label DISCARDED, too.
       The number of images to jump must always be positive, the jumpings
       cannot be done backwards.
       \\Ex: SF.jump(2); ---> Jump over 2 active images
       \\Ex: SF.jump(2,SelLine::ACTIVE); ---> The same
       \\Ex: SF.jump(2,SelLine::DISCARDED) ---> Jump over 2 discarded images
       @see get_current_file */
   void jump(int how_many, SelLine::Label label=SelLine::ACTIVE);
   
   /** Move "pointer" to a certain image filename.
       This function searches for an image name within the file, and
       locate the current line "pointer" pointing to that line. If
       the image name is not present (it is not the same "not present" and
       "discarded") in
       the selection file, then the pointer is moved to the end
       of the selection file. You can check this situation using eof().
       It doesn't matter if the current line "pointer" before the function
       call is after the line where the image name is, this function
       makes a search all over the file, regardless the previous
       situation of the current line "pointer".
       \\Ex: SF.search("g1ta0001");
       @see exists, get_current_file*/
   void search(string img_name) {current_line=find(img_name);}

   /** True if current line "pointer" is at the end of file.
       \\ Ex: if (SF.eof()) cout << "The selection file is over\n"; */
   int eof() {return current_line==text_line.end();}
   //@}

   // Getting information ..................................................
   /**@name Getting information*/
   //@{
   /** Returns the name of the file */
   FileName name() const {return fn_sel;}

   /** True if current line is a data line and it is active */
   int Is_ACTIVE() const
       {return current_line->Is_data() &&
               current_line->get_label()==SelLine::ACTIVE;}

   /** True if current line is a data line and it is active */
   int Is_DISCARDED() const
       {return current_line->Is_data() &&
               current_line->get_label()==SelLine::DISCARDED;}

   /** True if current line is a comment */
   int Is_COMMENT() const
       {return current_line->Is_comment();}

   /** Returns current line as a Sel Line. */
   const SelLine & current() {return *current_line;}

   /** True if the image name is inside the selection file.
       The current line "pointer" is not modified. If an image is
       discarded in the selection file, this function still will say
       that it exists, although it is discarded.
       \\Ex: if (SF.exists("g1ta0001")) cout << "g1ta0001 exists in the
             selection file\n"; */
   int exists(string img_name)
       {return find(img_name)!=text_line.end();}
   
   /** Number of images inside a selection file with a certain label.
       This function returns the number of images inside the
       selection file with a given label. By default this label
       is ACTIVE.
       \\Ex: cout << "There are " << SF.ImgNo() << " active images\n";
       \\Ex: cout << "There are " << SF.ImgNo(SelLine::ACTIVE) << " active images\n";
       \\Ex: cout << "There are " << SF.ImgNo(SelLine::DISCARDED) << " discarded images\n";*/
   int ImgNo(SelLine::Label label=SelLine::ACTIVE) const;

   /** Returns the number of lines within a file.
       This function gives the total number of lines (including comments)
       within a file.
       \\ Ex: cout << "There are " << SF.LineNo() << " lines in this file\n"; */
   int LineNo();
   
   /** Returns the size of the images inside.
       The filenames within a selection file are supposed to be for
       SPIDER images, this function opens one of the images (an active one)
       and returns the size of that image, supposed to be the same for
       the rest of the images in the selection file.
       
       An exception is thrown if the first valid image in the sel file,
       doesn't exist in the disk or it is not a XMIPP image.
       \\ Ex: SF.ImgSize(Ydim,Xdim); */
   void ImgSize(int &Ydim, int &Xdim) _THROW;
   
   /** Returns the extension of the files inside.
       This function returns the extension of the first active file.*/
   FileName FileExtension();
   
   /** Returns the maximum length of an active filename inside the selfile.
       The current pointer is not moved. */
   int MaxFileNameLength();

   /** Get the filename of the current line.
       If the current line "pointer" is at the end of the file or 
       is pointing to a comment then "" is returned.
       \\Ex: fn=SF.get_current_file(); */
   string get_current_file();
   
   // Get current line
   SelLine get_current_line() {return *current_line;}

   /** Gets statistics of the active images in the sel file
       it returns the average image, the minimum and maximum.
       \\Ex: SF.get_statistics(aveImg, min, max); */
   void get_statistics(Image& _ave, Image& _sd, double& _min, double& _max, bool apply_geo=FALSE);

   //@}

   // Modifying lines ......................................................
   /**@name Modifying the selection file*/
   //@{
   /** Removes an image from the selection file.
       This function searches for an image in the selection file, if it
       is found then the corresponding line is deleted. If the image
       is actually being pointed by the current line, then the current
       line is now the following line.
       \\ Ex: SF.remove("g1ta0001"); */
   void remove(string img_name);
   
   /** Removes actual line.
       This function removes the current line, either it is a comment or
       an image. The current line "pointer" is moved to the following
       line in the file.
       \\ Ex: SF.remove_current(); */
   void remove_current();
   
   /** Set label of an image.
       This function searches for an image inside the selection file and
       sets its label to the given label. If the image is not found in
       the file, then it is added at the end with the given label.
       The current line pointer is not modified.
       \\ Ex: SF.set("g1ta0001",SelLine::ACTIVE); */
   void set(string img_name, SelLine::Label label);

   /** Set the label of the current file.
       The same as the previous function but the label is set to
       the file currently pointed. */
   void set_current(SelLine::Label label);

   /** Change current filename.
       This function changes the current filename to a new one if it
       is not a comment. If it is a comment line, nothing is done. */
   void set_current_filename(const FileName &fn_new);

   /** Insert image before current line.
       There is no checking for the previous existence of the img.
       The current line is still pointing to the same line as it was
       before entering the function.
       \\ Ex: SF.insert("g1ta0000");
       \\ Ex: SF.insert("g1ta0000",SelLine::DISCARDED); */
   void insert(string img_name, SelLine::Label label=SelLine::ACTIVE);

   /** Insert line before current line.
       It is checked that the line is either a comment or data line,
       in this case that the label is right, too.
       The current line is still pointing to the same line as it was
       before entering the function.
       \\ Ex: SF.insert("g1ta0000",SelLine::ACTIVE); */
   void insert(const SelLine &_selline) _THROW;

   /** Insert a comment before the current line.
       Comments must not start with any special character since a "#" is
       automatically added at the beginning of the line.
       The current line is still pointing to the same line as it was
       before entering the function.
       \\ Ex: SF.insert_comment("This is a comment"); */
   void insert_comment(string comment);
   
   /** Deletes all DISCARDED images from the selection file.
       The current line "pointer" is moved to the beginning of the file.
       \\Ex: SF.clean(); */
   void clean();
   
   /** Deletes all comments from the selection file.
       The current line "pointer" is moved to the beginning of the file.
       \\Ex: SF.clean_comments(); */
   void clean_comments();
   //@}

   // Helpful procedures ...................................................
   /**@name Helpful procedures*/
   //@{
   /** Sort images in ascending order.
       All images are sorted in ascending order either they are active or
       discarded. All comments are gathered at the end of the resulting
       selection file. The current line of the resulting selection file
       is placed at the beginning of the file.
       \\ Ex: SF2=SF1.sort_by_filenames(); */
   SelFile sort_by_filenames();
   
   /** Alter order in sel file.
       A new selection file with all the images of the actual object
       (either they are active or discarded) is created, but this
       time all images are in a random order. The comments of the
       original selection file are lost in the new copy.
       The current line of the resulting selection file
       is placed at the beginning of the file.
       \\ SF2=SF1.randomize(); */
   SelFile randomize();
   
   /** Discard randomly N images.
       A set of N images are discarded from the actual selection file.
       If N is equal or greater than the actual number of images
       within the file, all images are discarded. Comments are kept
       at their original positions.
       The current line of the resulting selection file
       is placed at the beginning of the file.
       \\ Ex: SF2=SF1.random_discard(3); */
   SelFile random_discard(int N);

   /** Compare two selection files.
       The result is another selection file. At the beginning of it there
       is information about the number of active and discarded images on
       both input selection files, about the number of matching files
       (a file is said to match if it is active in both selection files),
       the number of active files which are only in the first selection
       file, and the number of active files which are only in the second.
       Then goes the list of matching files, the list of files only
       in SF1 and the list of files only in SF2. There are comments enough
       to know where things start and finish, and what the numbers are at
       the beginning. If a file is active in a file and discarded in the
       other, then it is said to match and it is kept as active, a
       preceeding comment warns of this situation.
       \\ Ex: SF3=compare(SF1,SF2);*/
   friend SelFile compare(SelFile &SF1, SelFile &SF2);

   /** Apply a function to all images with a certain label.
       The function must take an input image name and an output image
       name. Then transform the input image into the output one according
       to its functionality, and finally save the result as an Xmipp image.
       The output names are formed by adding an extension to the input ones,
       this extension can be empty and then the input and output images
       are the same. By default, the function is only applied to ACTIVE images
       and the output image is the same as the input one.
       \\Ex: \begin{verbatim}
             void function1(FileName _fn_in, FileName _fn_out) {
                ImageXmipp A(_fn_in);
                A()=A().reverseY();
                A.save(_fn_out);
             }
             
             ...
             SF.for_all(&function1);               --> Same I/O image, apply
                                                   --> only to ACTIVE images
             SF.for_all(&function1,"out");         --> Add ".out"
             SF.for_all(&function1,"",SelLine::DISCARDED);
                                                   --> Same I/O, apply only to
                                                   --> DISCARDED images
            \end{verbatim}*/
       void for_all(void (*f)(FileName, FileName), string _ext="",
          SelLine::Label _label=SelLine::ACTIVE);
   //@}
};

//@}
#endif
