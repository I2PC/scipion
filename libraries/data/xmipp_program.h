/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

#ifndef PROGRAM_H_
#define PROGRAM_H_

#include "argsprinter.h"
#include "xmipp_error.h"
#include "xmipp_strings.h"
#include "metadata.h"
#include "xmipp_image.h"
#include "xmipp_program_sql.h"


/** @defgroup Programs2 Basic structure for Xmipp programs
 *  @ingroup DataLibrary
 *
 * General program routines.
 *
 * Here you will find some routines which might help you to write programs
 * with a common structure where only a small part changes. These routines
 * make the heavy work for you and you only have to supply the changing
 * part from one to another.
 *
 * @{
 */

/** This class represent an Xmipp Program.
 * It have some of the basic functionalities of
 * the programs like argument parsing, checking, usage printing.
 */
class XmippProgram
{
private:
    /** Initialization function */
    void init();

    /** Function to check built-ins actions like --more, --help...etc */
    bool checkBuiltIns();

    /** Write Program info to DB */
    void writeToDB();
    /** Write protocol header information */
    void writeToProtocol();
    /** Write bash autocomplete options */
    void writeToAutocomplete();

    /** Create program GUI */
    /** By default, a simple Tk GUI is create based on parameters definition.
     * If an specific program implements a more specialized GUI, should redefine this function.
     */
    virtual void createGUI();
    /** Create Wiki for help */
    void createWiki();

    /** Variables related to progress notification */
    size_t progressTotal;
    size_t progressStep;
    size_t progressLast;

protected:
    /** Define Commons */
    void defineCommons();
    /** Value to store possible error codes */
    int errorCode;
    /// Program definition and arguments parser
    ProgramDef * progDef;

    /** Default comments */
    std::map<String, CommentList> defaultComments;

    ///Original command line arguments
    int argc;
    const char ** argv;

public:
    /** Flag to check whether to run or not*/
    bool doRun;
    /** Flag to mark if the program should run without arguments
     * or just print the usage(default behavior)
     */
    bool runWithoutArgs;
    /** @name Functions to be implemented by subclasses.
     * @{
     */
    /** Add the comments for a given default parameter */
    void processDefaultComment(const char *param, const char *left);
    /** Set default comment */
    void setDefaultComment(const char *param, const char *comment);
    /** Initialize comments for -v, ...*/
    virtual void initComments();
    /** Function in which the param of each Program are defined. */
    virtual void defineParams();
    /** Function in which each program will read parameters that it need.
     * If some error occurs the usage will be printed out.
     * */
    virtual void readParams();
    /** @name Program definitions
     * Following functions will be used in defineParams()
     * for define the arguments of
     * a program. Very useful for checking the command line parameters
     * and for standard usage message.
     * @{
     */
    /** Set the program name */
    void setProgramName(const char * name);
    /** Add usage line */
    void addUsageLine(const char * line, bool verbatim=false);
    /** Clear usage */
    void clearUsage();
    /** Add examples */
    void addExampleLine(const char * example, bool verbatim=true);
    /** Add other programs.
     * Separated by commas and without xmipp_
     */
    void addSeeAlsoLine(const char * seeAlso);
    /** Add keywords to the program definition */
    void addKeywords(const char * keywords);

    /** @} */

    /** Get the argument of this param, first start at 0 position */
    const char * getParam(const char * param, int arg = 0);
    const char * getParam(const char * param, const char * subparam, int arg = 0);
    int getIntParam(const char * param, int arg = 0);
    int getIntParam(const char * param, const char * subparam, int arg = 0);
    double getDoubleParam(const char * param, int arg = 0);
    double getDoubleParam(const char * param, const char * subparam, int arg = 0);
    /** Get arguments supplied to param as a list */
    void getListParam(const char * param, StringVector &list);
    /** Get the number of arguments supplied to the param */
    int getCountParam(const char * param);
    /** Check if the param was supplied to command line */
    bool checkParam(const char * param);
    /** Return true if the program is defined */
    bool existsParam(const char * param);

    /** Add a params definition line*/
    void addParamsLine(const String &line);
    void addParamsLine(const char * line);

    /** Get Parameter definition */
    ParamDef * getParamDef(const char * param) const;

    /// Verbosity level
    int verbose;
    /** debug flag and seed for randomization */
    int debug, seed;

    /** @name Public common functions
     * The functions in this section are available for all programs
     * and should not be reimplemented in derived class, since they
     * are thought to have the same behaivor and it depends on the
     * definition of each program.
     * @{
     */
    /** Returns the name of the program
     * the name of the program is defined
     * by each subclass of this base class.
     */
    const char * name() const;
    /** Print the usage of the program, reduced version */
    virtual void usage(int verb = 0) const;
    /** Print help about specific parameter */
    virtual void usage(const String & param, int verb = 2);

    /** Return the version of the program */
    int version() const;

    /** Show parameters */
    virtual void show() const;

    /** Read the command line arguments
     * If an error occurs while reading arguments,
     * the error message will be printed and the
     * usage of the program showed. So you don't need
     * to do that in readParams();
     * */
    virtual void read(int argc, const char ** argv, bool reportErrors = true);

    /** Read the command line arguments
     * A convenience wrapper
     * */
    virtual void read(int argc, char ** argv, bool reportErrors = true);

    /** Read arguments from an string.
     * This function should do the same as reading arguments
     * but first convert the string to arguments.
     */
    void read(const String &argumentsLine);

    /** @} */
    /** This function will be start running the program.
     * it also should be implemented by derived classes.
     */
    virtual void run();
    /** function to exit the program
     * can be usefull redefined for mpi programs
     */
    virtual void quit(int exit_code = 0) const;
    /** Call the run function inside a try/catch block
     * The function will return the error code when
     * 0 means success
     * and a value greater than 0 represents the error type
     * */
    virtual int tryRun();

    /** Functions related to progress notification */

    /** Set the total amount of work and initialize
     * the progress to 0
     * The update step is calculated as step = XMIPP_MAX(1, total/stepBin)
     * so, the updates only are done at step multiples
     */
    void initProgress(size_t total, size_t stepBin = 60);

    /** Notify progress on work */
    void setProgress(size_t value = 0);

    /** Notify end of work */
    void endProgress();



    /** @name Constructors
     * @{
     */
    /** Constructor */
    XmippProgram();
    /** Constructor for read params */
    XmippProgram(int argc, const char ** argv);
    /** Destructor */
    virtual ~XmippProgram();
    /** @} */

}
;//end of class XmippProgram

/** Special class of XmippProgram that performs some operation related with processing images.
 * It can receive a file with images(MetaData) or a single image.
 * The function processImage is virtual here and needs to be implemented by derived classes.
 * Optionally can be implemented preProcess and postProcess to perform some customs actions.
 */
class XmippMetadataProgram: public virtual XmippProgram
{
private:
    /// Input and output metadatas
    MetaData * mdIn;
    MetaData mdOut; //TODO: can be treated by reference as mdIn for
    // uses from another programs...
public:
    /// The input metadata should not be used
    /// if there is a very very special case
    /// you can use this function
    MetaData * getInputMd()
    {
        return mdIn;
    }
    MetaData * getOutputMd()
    {
        return &mdOut;
    }

public:
    //Image<double>   img;
    /// Filenames of input and output Metadata
    FileName fn_in, fn_out, baseName, pathBaseName, oextBaseName;
    /// Apply geo
    bool apply_geo;
    /// Output dimensions
    size_t ndimOut, zdimOut, ydimOut, xdimOut;
    DataType datatypeOut;
    /// Number of input elements
    size_t mdInSize;

protected:
    /// Metadata writing mode: OVERWRITE, APPEND
    WriteModeMetaData mode;

    /// Iterator over input metadata
    MDIterator * iter;
    /// Filenames of input and output Images
    //FileName        fnImg, fnImgOut;
    /// Output extension and root
    FileName oext, oroot;
    /// MDLabel to be used to read/write images, usually will be MDL_IMAGE
    MDLabel image_label;


    // BEHAVIOR CONTROL FLAGS //

    /// Indicate that a unique final output is produced
    bool produces_an_output; // Default false (only -o param is used)
    /// Indicate that the unique final output file is a Metadata
    bool produces_a_metadata; // Default false (if true, then produces_an_output is set true)
    /// Indicate that an output is produced for each image in the input
    bool each_image_produces_an_output; // Default false (both -o --oroot params are used)
    /// Provide the program with the param --dont_apply_geo to allow the user deciding whether
    /// or not applying the transformation info as stored in the input metadata
    bool allow_apply_geo; // Default false
    /// Input Metadata will treat a stack file as a set of images instead of a unique file
    bool decompose_stacks; // Default true
    /// Delete previous output stack file prior to process images
    bool delete_output_stack; // Default true
    /// Get the input image file  dimensions to further operations
    bool get_image_info; // Default true
    /// Save the associated output metadata when output file is a stack
    bool save_metadata_stack; // Default false
    /// Include the original input image filename in the output stack
    bool track_origin; // Default false
    /// Keep input metadata columns
    bool keep_input_columns; // Default false
    /// Remove disabled images from the input selfile
    bool remove_disabled; // Default true
    /// Show process time bar
    bool allow_time_bar; // Default true

    // DEDUCED FLAGS
    /// Input is a metadata
    bool input_is_metadata;
    /// Input is a single image
    bool single_image;
    /// Input is a stack
    bool input_is_stack;
    /// Output is a stack
    bool output_is_stack;
    // Create empty output stack file prior to process images
    bool create_empty_stackfile; //
    //check whether to delete or not the input metadata
    bool delete_mdIn;

    /// Some time bar related counters
    size_t time_bar_step, time_bar_size, time_bar_done;

    virtual void initComments();
    virtual void defineParams();
    virtual void readParams();
    virtual void preProcess();
    virtual void postProcess();
    virtual void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut) = 0;
    virtual void show();
    /** Do some stuff before starting processing
     * in a parallel environment usually this only be executed
     * by master.
     */
    virtual void startProcessing();
    virtual void finishProcessing();
    virtual void showProgress();
    /** This method will be used to distribute the images to process
     * it will set the objectId and objectIndexto read from input metadata
     * or -1 if there are no more images to process.
     * This method will be useful for parallel task distribution
     */
    virtual bool getImageToProcess(size_t &objId, size_t &objIndex);

    /** Define the label param */
    virtual void defineLabelParam();

public:
    XmippMetadataProgram();

    /** Call the read function inside a try/catch block
     * The function will return the error code when
     * 0 means success
     * and a value greater than 0 represents the error type
     * */
    virtual int tryRead(int argc, const char ** argv, bool reportErrors = true);

    /** Initialization of variables should be done here
     */
    virtual void init();

    /** Setup of common XmippMetadataProgram arguments
     *  to be called from another program.
     */
    virtual void setup(MetaData *md, const FileName &o="", const FileName &oroot="",
                       bool applyGeo=false, MDLabel label=MDL_IMAGE);


    /** Destructor
     */
    virtual ~XmippMetadataProgram()
    {
        if (delete_mdIn)
            delete mdIn;
    }

    void setMode(WriteModeMetaData _mode)
    {
        mode = _mode;
    }

    /// Prepare rowout
    void setupRowOut(const FileName &fnImgIn, const MDRow &rowIn, const FileName &fnImgOut, MDRow &rowOut) const;

    void run1();
    void run2();
    virtual void run();
}
;// end of class XmippMetadataProgram

/** This class will serve as an interface for python scripts
 * useful for command line parsing and help message printing
 */
class XmippProgramGeneric: public XmippProgram
{
public:
    bool definitionComplete;
    ///Constructor
    XmippProgramGeneric();
    void endDefinition();
    virtual void read(int argc, const char ** argv, bool reportErrors = true);

protected:
    void defineParams();
    void readParams();
    void show();
    void run();
}
;// end of class XmippProgramGeneric

/** This macro will be useful for define the main and run
 * an XmippProgram
 */
#define RUN_XMIPP_PROGRAM(progName) \
  int main(int argc, char** argv) { \
      progName program; program.read(argc, argv);\
      return program.tryRun();}


///Declare all programs
/** @} */

#endif /* PROGRAM_H_ */
