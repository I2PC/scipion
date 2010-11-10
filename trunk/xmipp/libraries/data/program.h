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

#include "argsparser.h"
#include "error.h"
#include "strings.h"
#include "metadata.h"
#include "image.h"
#include "program_sql.h"

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
 * the programs like argument parsing, cheking, usage printing.
 */
class XmippProgram
{
private:
    /** Initialization function */
    void init();

    /** Function to check built-ins actions like --more, --help...etc */
    void checkBuiltIns();

    /** Write Program info to DB */
    void writeToDB(const FileName &dbName);

    /** Create Tk GUI */
    void createGUI();

protected:
    /// Program definition and arguments parser
    ProgramDef * progDef;

    ///Original command line arguments
    int argc;
    char ** argv;

    /** @name Functions to be implemented by subclasses.
     * @{
     */
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
    void addUsageLine(const char * line);
    /** Clear usage */
    void clearUsage();

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
    /** Get params supplied as a vector */
    void getListParam(const char * param, StringVector &list);
    /** Check if the param was supplied to command line */
    bool checkParam(const char * param);
    /** Return true if the program is defined */
    bool existsParam(const char * param);

public:
    /** Add a params definition line*/
    void addParamsLine(const char * line);
    /// Verbosity level
    int verbose;

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
    virtual void usage(const std::string & param, int verb = 2);

    /** Return the version of the program */
    int version() const;

    /** Show parameters */
    virtual void show() const;

    /** Read the command line arguments
     * If an error occurs while reading arguments,
     * the error message will be printed and the
     * usage of the program showed. So the is no need
     * to do that in readParams();
     * */
    void read(int argc, char ** argv, bool reportErrors = true);

    /** @} */
    /** This function will be start running the program.
     * it also should be implemented by derived classes.
     */
    virtual void run() = 0;
    /** Call the run function inside a try/catch block */
    void tryRun();

    /** @name Constructors
     * @{
     */
    /** Constructor */
    XmippProgram();
    /** Constructor for read params */
    XmippProgram(int argc, char ** argv);
    /** Destructor */
    virtual ~XmippProgram();
    /** @} */

};

/** Special class of XmippProgram that performs some operation related with processing images.
 * It can receive a file with images(MetaData) or a single image.
 * The function processImage is virtual here and needs to be implemented by derived classes.
 * Optionally can be implemented preProcess and postProcess to perform some customs actions.
 */
class XmippMetadataProgram: public XmippProgram
{
protected:
    Image<double>   img;
    /// Filenames of input and output Metadata
    FileName        fn_in, fn_out;
    /// Metadata Blockname (several metadata may go in the same file)
    std::string blockName;
    /// Input and output metadatas
    MetaData        mdIn, mdOut;
    /// Filenames of input and output Images
    FileName        fnImg, fnImgOut;
    /// Output extension and root
    FileName oext, oroot;
    /// Set this to true for allow the option
    /// of applicating the transformation as stored in the header
    bool apply_geo;
    /// Set this to true when the program produce a result metadata
    bool produces_an_output;
    /// Use this flag for not writing at every image
    bool each_image_produces_an_output;
    /// Use this flag for not producing a time bar
    bool allow_time_bar;

    virtual void defineParams();
    virtual void readParams();
    virtual void preProcess();
    virtual void postProcess();
    virtual void processImage() = 0;
    virtual void show();

public:
    XmippMetadataProgram();
    virtual void run();
}
;// end of class XmippMetadataProgram

/** @} */

#endif /* PROGRAM_H_ */
