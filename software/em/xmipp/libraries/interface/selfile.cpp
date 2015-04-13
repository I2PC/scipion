/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Part of this module has been developed by Lorenzo Zampighi and Nelson Tang
 * Dept. Physiology of the David Geffen School of Medistd::cine
 * Univ. of California, Los Angeles.
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
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <algorithm>

#include <data/xmipp_funcs.h>
#include "selfile.h"
#include <data/xmipp_image.h>

/*****************************************************************************/
/* SEL FILE LINE                                                 */
/*****************************************************************************/
/* Copy Constructor */
SelLine::SelLine(const SelLine &l)
{
    line_type = l.line_type;
    text      = l.text;
    label     = l.label;
    number    = l.number;
}

SelLine& SelLine::operator = (const SelLine &SL)
{
    if (this != &SL)
    {
        line_type = SL.line_type;
        text      = SL.text;
        label     = SL.label;
        number    = SL.number;
    }
    return *this;
}

void SelLine::assign(const SelLine &SL)
{
    *this = SL;
}

// Returns false if the comparison is nonsense, ie, if line_types are not the
// same or the lines are to be deleted or are not assigned.
// Returns TRUE if l1<l2 otherwise false
bool operator < (const SelLine &l1, const SelLine &l2)
{
    if (l1.line_type < l2.line_type)
        return 1;
    else if (l1.line_type > l2.line_type)
        return 0;
    else
        return l1.text < l2.text;
}

std::ostream& operator << (std::ostream& o, const SelLine &SFL)
{
    switch (SFL.line_type)
    {
    case (SelLine::DATALINE):
                    o << SFL.text << " " << SFL.label << std::endl;
        break;
    case (SelLine::COMMENT):
                    o << SFL.text << std::endl;
        break;
    default:
    	break;
    }
    return o;
}

std::istream& operator >> (std::istream& o, SelLine &SFL)
{
    std::string line;
    char     img_name[1024];
    int      no_elements_read;
    int      label;

    // Get line
    getline(o, line);

    // Initialise target
    SFL.line_type = SelLine::NOT_ASSIGNED;
    SFL.text = "";
    SFL.label = SelLine::DISCARDED;
    if (line.length() == 0)
        return o;

    // Check if comment or empty line
    if (line[0] == '#' || line[0] == '\0' || line[0] == ';')
    {
        line[line.length()-1] = '\0';
        SFL.line_type = SelLine::COMMENT;
        SFL.text = line;

        // Check if a true "filename label" line
    }
    else
    {
        no_elements_read = sscanf(line.c_str(), "%s %d", img_name, &label);
        // *** THE SSCANF CAN BE REPLACED BY A STRING I/O OPERATION
        if (no_elements_read == 2)
        {
            SFL.line_type = SelLine::DATALINE;
            SFL.text = img_name;
            SFL.label = (label >= 0) ? SelLine::ACTIVE : SelLine::DISCARDED;
            SFL.number = label;
        }
        else
            REPORT_ERROR(ERR_SELFILE, "Format error when reading Selection line");
    }
    return o;
}

/*****************************************************************************/
/* SEL FILE                                                     */
/*****************************************************************************/
/* Constructor ------------------------------------------------------------- */
SelFile::SelFile()
{
    fn_sel       = "Unnamed";
    no_imgs      = 0;
    current_line = text_line.begin();
}

/* Copy Constructor -------------------------------------------------------- */
SelFile::SelFile(const SelFile &SF)
{
    fn_sel       = SF.fn_sel;
    text_line    = SF.text_line;
    no_imgs      = SF.no_imgs;
    current_line = SF.current_line;
}

/* Clear ------------------------------------------------------------------- */
void SelFile::clear()
{
    fn_sel = "Unnamed";
    text_line.erase(text_line.begin(), text_line.end());
    no_imgs = 0;
    current_line = text_line.begin();
}

/* Assignment -------------------------------------------------------------- */
SelFile& SelFile::operator = (const SelFile &SF)
{
    if (this != &SF)
    {
        fn_sel       = SF.fn_sel;
        text_line    = SF.text_line;
        no_imgs      = SF.no_imgs;
        current_line = SF.current_line;
    }
    return *this;
}

/* Another function for assignment ----------------------------------------- */
void SelFile::assign(const SelFile &SF)
{
    *this = SF;
}

/* Show Sel file ----------------------------------------------------------- */
std::ostream& operator << (std::ostream& o, const SelFile &SF)
{
    std::vector<SelLine>::const_iterator current = SF.text_line.begin();
    std::vector<SelLine>::const_iterator last    = SF.text_line.end();
    while (current != last)
    {
        o << *current;
        current++;
    }
    return o;
}

/* Clean ------------------------------------------------------------------- */
void SelFile::clean()
{
    std::vector<SelLine>::iterator current = text_line.begin();
    while (current != text_line.end())
    {
        if ((*current).line_type == SelLine::DATALINE &&
            (*current).label == SelLine::DISCARDED)
        {
            text_line.erase(current);
        }
        else
            current++;
    }
    current_line = text_line.begin();
}

/* Clean comments ---------------------------------------------------------- */
void SelFile::clean_comments()
{
    std::vector<SelLine>::iterator current = text_line.begin();
    std::vector<SelLine>::iterator last    = text_line.end();
    std::vector<SelLine>::iterator temp;
    while (current != last)
    {
        if ((*current).line_type == SelLine::COMMENT)
        {
            temp = current;
            temp++;
            text_line.erase(current);
            current = temp;
        }
        else
            current++;
    }
    current_line = text_line.begin();
}

/* Read -------------------------------------------------------------------- */
void SelFile::read(const FileName &sel_name, int overriding)
{
    SelLine   temp;
    std::ifstream  fh_sel;
    int       line_no = 1;

    // Empties current SelFile
    if (overriding)
        clear();

    // Open file
    else
    {
        // Read normal selfile
        fh_sel.open(sel_name.c_str(), std::ios::in);
        if (!fh_sel)
            REPORT_ERROR(ERR_IO_NOTEXIST, sel_name);

        // Read each line and keep it in the list of the SelFile object
        fh_sel.peek();
        while (!fh_sel.eof())
        {
            try
            {
                fh_sel >> temp;
            }
            catch (XmippError)
            {
                std::cout << "Sel file: Line " << line_no << " is skipped due to an error\n";
            }
            switch (temp.line_type)
            {
            case (SelLine::NOT_ASSIGNED): break; // Line with an error
            case (SelLine::DATALINE):
                            if (temp.label != SelLine::DISCARDED)
                                no_imgs++;
                text_line.push_back(temp);
                break;
            case (SelLine::COMMENT):
                            text_line.push_back(temp);
                break;
            default:
            	break;
            }
            line_no++;
            fh_sel.peek();
        }

        // Close file
        fh_sel.close();
    }

    // Set "pointer" to the beginning of the file
    if (overriding)
        fn_sel = sel_name;
    go_first_ACTIVE();
}
/* Merge ------------------------------------------------------------------- */
void SelFile::merge(const FileName &sel_name)
{
    SelFile SF(sel_name);
    *this = *this + SF;
    go_first_ACTIVE();
}

/* Write ------------------------------------------------------------------- */
void SelFile::write(const FileName &sel_name)
{
    std::ofstream    fh_sel;
    std::vector<SelLine>::iterator current = text_line.begin();
    std::vector<SelLine>::iterator last    = text_line.end();

    if (strcmp(sel_name.c_str(), "") != 0)
        fn_sel = sel_name;
    // Don't use sel_name=="" because it wastes memory
#ifdef NEVEREVERROB
    if (sel_name.find(IMAGIC_TAG) == 0)
    {
        // Write Imagic selfile
        const FileName hed_fname = sel_name.substr(IMAGIC_TAG_LEN);
        std::vector<Image *> imgs;
        for (; current != last; current++)
        {
            Image *img;
            if (current->Is_data() && (current->get_label() == SelLine::ACTIVE) &&
                (img = Image::LoadImage(current->get_text())))
                imgs.push_back(img);
        }
        if (!ImagicWriteImagicFile(hed_fname, imgs))
            REPORT_ERROR(1553, "Error writing selfile to Imagic file " + sel_name);
        for (std::vector<Image *>::iterator i = imgs.begin(); i != imgs.end(); i++)
            delete(*i);
    }
    else
#endif
    {
        // Write Xmipp selfile
        // Open file
        fh_sel.open(fn_sel.c_str(), std::ios::out);
        if (!fh_sel)
            REPORT_ERROR(ERR_IO_NOWRITE, fn_sel);

        // Read each line and keep it in the list of the SelFile object
        while (current != last)
            fh_sel << *(current++);

        // Close file
        fh_sel.close();
    }
}

/* Merging with another selfile -------------------------------------------- */
void SelFile::merge(SelFile &SF)
{
    std::vector<SelLine>::iterator current = SF.text_line.begin();
    std::vector<SelLine>::iterator last    = SF.text_line.end();
    std::vector<SelLine>::iterator found;

    SelLine discrepancy;
    discrepancy.line_type = SelLine::COMMENT;
    discrepancy.text = "# There were discrepancy in the tags for next line, the "
                       "ACTIVE state is kept";

    while (current != last)
    {
        if ((*current).line_type != SelLine::DATALINE)
        {
            current++;
            continue;
        }
        if ((found = find((*current).text)) == text_line.end())
        {
            // New image not found in the whole Sel File.
            // Add it if it is not discarded
            if ((*current).label != SelLine::DISCARDED)
            {
                text_line.push_back(*current);
                no_imgs++;
            }
        }
        else
            // New image is found, check that its line is not going
            // to be removed, if it is add it again; else, check if
            // there is a discrepancy between them
            if ((*found).label != (*current).label)
            {
                if ((*found).label < (*current).label)
                {
                    (*found).label = SelLine::ACTIVE;
                    no_imgs++;
                }
                text_line.insert(found, 1, discrepancy);
            }
        current++;
    }
}

/* Merging, operator + ----------------------------------------------------- */
// If the same file is in both Sel Files the label in the first is kept
SelFile SelFile::operator + (SelFile &SF)
{
    SelFile result;
    result = *this;
    result.merge(SF);
    return result;
}

/* Split randomly in two equally large selfiles ---------------------------- */
void SelFile::split_in_two(SelFile &SF1, SelFile &SF2)
{
    SelFile  SFtmp;
    SF1 = *this;
    SFtmp = SF1.randomize();
    SF1.clear();
    int N = SFtmp.ImgNo();
    SF1.reserve(N);
    SF2.reserve(N);
    int half = N / 2;
    SFtmp.go_beginning();
    for (int i = 0;i < N; i++)
    {
        if (i < half)
            SF1.insert(SFtmp.current());
        else
            SF2.insert(SFtmp.current());
        if (i < N - 1)
            SFtmp.NextImg();
    }
    SFtmp = SF1.sort_by_filenames();
    SF1 = SFtmp;
    SFtmp = SF2.sort_by_filenames();
    SF2 = SFtmp;
}

/* Split randomly in N equally large selfiles ------------------------------ */
void SelFile::split_in_N(int N, std::vector<SelFile> &SF)
{
    // Randomize input data
    SelFile  SFtmp, SFrnd;
    SFrnd = *this;
    SFtmp = SFrnd.randomize();
    SFtmp.go_beginning();
    int Nimg = SFtmp.ImgNo();
    SF.clear();

    // Create space for all SelFiles
    for (int n = 0; n < N; n++)
    {
        SelFile *ptr_SF = new SelFile;
        ptr_SF->reserve(CEIL(Nimg / N));
        SF.push_back(*ptr_SF);
    }

    // Distribute images
    int n = 0;
    for (int i = 0;i < Nimg; i++)
    {
        SF[n].insert(SFtmp.current());
        n = (n + 1) % N;
        if (i < Nimg - 1)
            SFtmp.NextImg();
    }

    // Sort the Selfiles
    for (int n = 0; n < N; n++)
        SF[n] = SF[n].sort_by_filenames();
}

/* Select only part of the selfile for parallel MPI-runs ------------------ */
void SelFile::mpi_select_part(int rank, int size, int &num_img_tot)
{

    (*this).clean_comments();
    (*this).clean();
    num_img_tot = (*this).ImgNo();
    int remaining = num_img_tot % size;
    int Npart = (int)(num_img_tot - remaining) / size;
    int myFirst, myLast;
    if (rank < remaining)
    {
        myFirst = rank * (Npart + 1);
        myLast = myFirst + Npart;
    }
    else
    {
        myFirst = rank * Npart + remaining;
        myLast = myFirst + Npart - 1;
    }
    // Now discard all images in Selfile that are outside myFirst-myLast
    (*this).go_beginning();
    SelFile  SFpart = *this;
    SFpart.clear();
    for (int nr = myFirst; nr <= myLast; nr++)
    {
        (*this).go_beginning();
        (*this).jump_lines(nr);
        SFpart.insert((*this).current());
    }
    *this = SFpart;

}
/* Select only part of the selfile for parallel MPI-runs ------------------ */
void SelFile::mpi_select_part2(int jobNumber, 
                               int numberJobs, 
                               int &totalNumImg,
                               int mpi_job_size)
{   // jobNumber process number
    // total number of processes
    // total number of images

    (*this).clean_comments();
    (*this).clean();
    totalNumImg = (*this).ImgNo();
    int myFirst = jobNumber * mpi_job_size;
    int myLast  = myFirst + mpi_job_size-1;
    while ((myLast+1) > totalNumImg)
    {
        myLast = totalNumImg-1;
    }
    // Now discard all images in Selfile that are outside myFirst-myLast
    (*this).go_beginning();
    SelFile  SFpart = *this;
    SFpart.clear();
    for (int nr = myFirst; nr <= myLast; nr++)
    {
        (*this).go_beginning();
        (*this).jump_lines(nr);
        SFpart.insert((*this).current());
    }
    *this = SFpart;
}

/* Choose subset ----------------------------------------------------------- */
void SelFile::chooseSubset(int firstImage, int lastImage, SelFile &SFsubset)
{
    SFsubset.clear();
    go_beginning();
    jump(firstImage);
    for (int i=firstImage; i<=lastImage; i++)
    {
        if (!eof()) SFsubset.insert(current());
        next();
    }
}

/* Adjust to label --------------------------------------------------------- */
void SelFile::adjust_to_label(SelLine::Label label)
{
    if (current_line == text_line.end())
        return;
    while ((*current_line).line_type != SelLine::DATALINE ||
           (*current_line).label != label)
    {
        current_line++;
        if (current_line == text_line.end())
            return;
    }
}

/* Next Image with a certain label ----------------------------------------- */
const std::string& SelFile::NextImg(SelLine::Label label)
{
    adjust_to_label(label);
    static const std::string emptyString;
    if (current_line != text_line.end())
        return (*current_line++).text;
    else
        return emptyString;
}

/* Jump over a certain number of data lines (disregarding any label) ------- */
bool SelFile::jump_lines(int how_many)
{
    for (int i = 0; i < how_many; i++)
    {
        if (current_line != text_line.end())
            current_line++;
        else
            return false;
    }
    return true;
}

/* Jump over images with a certain label ----------------------------------- */
void SelFile::jump(int how_many, SelLine::Label label)
{
    adjust_to_label(label);
    for (int i = 0; i < how_many; i++)
        if (current_line != text_line.end())
        {
            current_line++;
            adjust_to_label(label);
        }
}

/* Find an image (inside the list) ----------------------------------------- */
// It returns a pointer to past-last element if the image is not inside
std::vector<SelLine>::iterator find(std::vector<SelLine> &text,
    const std::string &img_name)
{
    std::vector<SelLine>::iterator current = text.begin();
    std::vector<SelLine>::iterator last    = text.end();

    while (current != last)
    {
        if ((*current).line_type == SelLine::DATALINE &&
            (*current).text == img_name)
            return current;
        current++;
    }
    return current;
}

/* Find an image (inside the Sel File) ------------------------------------- */
// It returns a pointer to past-last element if the image is not inside
// *** THIS SHOULD USE THE PREVIOUS FUNCTION BUT I CANNOT MAKE IT TO COMPILE
std::vector<SelLine>::iterator SelFile::find(const std::string &img_name)
{
    std::vector<SelLine>::iterator current = text_line.begin();
    std::vector<SelLine>::iterator last    = text_line.end();

    while (current != last)
    {
        if ((*current).line_type == SelLine::DATALINE &&
            (*current).text == img_name)
            return current;
        current++;
    }
    return current;
}

/* Number of images with a certain label ----------------------------------- */
// If the label is 0 it means any valid image
int SelFile::ImgNo(SelLine::Label label) const
{
    int N = 0;
    std::vector<SelLine>::const_iterator current = text_line.begin();
    std::vector<SelLine>::const_iterator last    = text_line.end();
    while (current != last)
    {
        if ((*current).line_type == SelLine::DATALINE &&
            (*current).label == label)
            N++;
        current++;
    }
    return N;
}

/* Number of lines within file --------------------------------------------- */
int SelFile::LineNo()
{
    int N = 0;
    std::vector<SelLine>::iterator current = text_line.begin();
    std::vector<SelLine>::iterator last    = text_line.end();
    while (current != last)
    {
        N++;
        current++;
    }
    return N;
}


/* File Extension ---------------------------------------------------------- */
FileName SelFile::FileExtension()
{
    std::vector<SelLine>::iterator aux = current_line;
    go_first_ACTIVE();
    FileName ext = (*current_line).text;
    ext = ext.getExtension();
    current_line = aux;
    return ext;
}


/* Maximum filename length ------------------------------------------------- */
int SelFile::MaxFileNameLength()
{
    std::vector<SelLine>::iterator aux = current_line;
    size_t max_length = 0;
    go_first_ACTIVE();
    while (!eof())
    {
        FileName fn = NextImg();
        max_length = XMIPP_MAX(max_length, fn.length());
    }
    current_line = aux;
    return max_length;
}

/* Get current filename ---------------------------------------------------- */
const std::string SelFile::get_current_file()
{
    if (current_line == text_line.end())
        return "";
    if ((*current_line).line_type != SelLine::DATALINE)
        return "";
    return (*current_line).text;
}

/* Get filename number i --------------------------------------------------- */
const std::string SelFile::get_file_number(int i)
{
    if (i < 0)
        return "";
    std::vector<SelLine>::iterator current = text_line.begin();
    std::vector<SelLine>::iterator last    = text_line.end();

    int currenti = 0;
    while (current != last)
    {
        if ((*current).line_type == SelLine::DATALINE &&
            (*current).label == SelLine::ACTIVE)
            currenti++;
        if (currenti > i)
            return (*current).text;
        current++;
    }
    return "";
}

/* Remove a certain file --------------------------------------------------- */
void SelFile::remove(const std::string &img_name)
{
    std::vector<SelLine>::iterator aux = find(img_name);
    std::vector<SelLine>::iterator temp;
    if (aux != text_line.end())
    {
        if (aux == current_line)
        {
            temp = current_line;
            temp++;
        }
        else
            temp = current_line;
        if ((*aux).line_type == SelLine::DATALINE)
            no_imgs--;
        text_line.erase(aux);
        current_line = temp;
    }
}

/* Remove current line ----------------------------------------------------- */
void SelFile::remove_current()
{
    if (current_line != text_line.end())
    {
        std::vector<SelLine>::iterator temp;
        temp = current_line;
        temp++;
        if ((*current_line).line_type == SelLine::DATALINE)
            no_imgs--;
        text_line.erase(current_line);
        current_line = temp;
    }
}

/* Append a file or change label ------------------------------------------- */
void SelFile::set(const std::string& img_name, SelLine::Label label)
{
    SelLine temp;
    std::vector<SelLine>::iterator aux = find(img_name);
    if (aux == text_line.end())
    {
        temp.line_type = SelLine::DATALINE;
        temp.text = img_name;
        temp.label = label;
        text_line.push_back(temp);
        if (label != SelLine::DISCARDED)
            no_imgs++;
    }
    else
    {
        if ((*aux).label != label)
        {
            (*aux).label = label;
            if (label != SelLine::DISCARDED)
                no_imgs++;
        }
    }
}

/* Append a file or change label ------------------------------------------- */
void SelFile::set_current(SelLine::Label label)
{
    if ((*current_line).label != label)
    {
        (*current_line).label = label;
        if (label != SelLine::DISCARDED)
            no_imgs++;
    }
}

/* Change current filename ------------------------------------------- */
void SelFile::set_current_filename(const FileName &fn_new)
{
    if ((*current_line).line_type == SelLine::DATALINE)
    {
        (*current_line).text = fn_new;
    }
}

/* Insert image before current line ---------------------------------------- */
void SelFile::insert(const std::string& img_name, SelLine::Label label)
{
    SelLine temp;
    temp.line_type = SelLine::DATALINE;
    temp.text = img_name;
    temp.label = label;
    if (label != SelLine::DISCARDED)
        no_imgs++;

    // Insert and updates current_line
    current_line = text_line.insert(current_line, temp);
    current_line++;
}

/* Insert line before current line ----------------------------------------- */
void SelFile::insert(const SelLine &_selline)
{
    if (_selline.line_type != SelLine::DATALINE &&
        _selline.line_type != SelLine::COMMENT)
        REPORT_ERROR(ERR_SELFILE, "SelFile::insert(SelLine): SelLine type not valid");
    if (_selline.line_type == SelLine::DATALINE)
        if (_selline.label != SelLine::DISCARDED &&
            _selline.label != SelLine::ACTIVE)
            REPORT_ERROR(ERR_SELFILE, "SelFile::insert(SelLine): SelLine label not valid");

    // Sjors 18sep06: added next line
    if (_selline.label != SelLine::DISCARDED)
        no_imgs++;

    // Insert and updates current_line
    current_line = text_line.insert(current_line, _selline);
    current_line++;
}

/* Insert a comment before current line ------------------------------------ */
void SelFile::insert_comment(const std::string& comment)
{
    SelLine temp;
    temp.line_type = SelLine::COMMENT;
    temp.text = "# " + comment;
    temp.label = SelLine::DISCARDED;

    // Insert and updates current_line
    current_line = text_line.insert(current_line, temp);
    current_line++;
}

/* Sort -------------------------------------------------------------------- */
SelFile SelFile::sort_by_filenames()
{
    SelFile result(*this);
    sort(result.text_line.begin(), result.text_line.end());
    result.current_line = result.text_line.begin();
    return result;
}

/* Randomize --------------------------------------------------------------- */
SelFile SelFile::randomize()
{
    SelFile  result, aux;
    int      i;
    int      rnd_indx;

    randomize_random_generator();
    if (no_imgs == 0)
        return aux;
    aux = *this;
    for (i = no_imgs; i > 0; i--)
    {
        // Jump a random number from the beginning
        rnd_indx = (int) rnd_unif(0, i);
        aux.go_first_ACTIVE();
        aux.jump(rnd_indx);
        result.text_line.push_back(*(aux.current_line));
        (*aux.current_line).line_type = SelLine::NOT_CONSIDERED;
    }

    // Adjust remaining fields
    result.no_imgs = no_imgs;
    result.current_line = result.text_line.begin();
    return result;
}


/* Discard randomly a set of images ---------------------------------------- */
SelFile SelFile::random_discard(int N)
{
    SelFile  result;
    int      i, rnd_indx;

    SelLine::Label label = SelLine::ACTIVE;
    result = *this;
    N = std::min(N, no_imgs);
    for (i = 0; i < N; i++)
    {
        // Jump a random number from the beginning
        rnd_indx = (int) rnd_unif(0, result.no_imgs);
        result.go_first_ACTIVE();
        result.jump(rnd_indx, label);

        // Discard that image
        (*(result.current_line)).label = SelLine::DISCARDED;

        // Decrease the number of images such that next time
        result.no_imgs--;
    }

    result.go_beginning();
    return result;
}

/* Compare ----------------------------------------------------------------- */
// Only img_files with the active label are considered
SelFile compare(SelFile &SF1, SelFile &SF2, const int mode)
{
    std::vector<SelLine>     only_in_SF1;
    std::vector<SelLine>     only_in_SF2;
    std::vector<SelLine>     in_both;
    SelFile           result;
    SelLine           temp;
    int               SF1_discarded = 0, SF2_discarded = 0;
    char              str[10];

    // Search in File 1
    std::vector<SelLine>::iterator current = SF1.text_line.begin();
    std::vector<SelLine>::iterator last    = SF1.text_line.end();
    std::vector<SelLine>::iterator last_SF = SF2.text_line.end();
    std::vector<SelLine>::iterator found;

    while (current != last)
    {
        // Skip if not active
        if ((*current).line_type != SelLine::DATALINE)
        {
            current++;
            continue;
        }
        if ((*current).label == SelLine::DISCARDED)
        {
            SF1_discarded++;
            current++;
            continue;
        }

        // Try to find this archive into Sel File 2
        found = SF2.find((*current).text);
        if (found == last_SF)
            only_in_SF1.push_back(*current);
        else
            if ((*found).label == SelLine::DISCARDED)
                only_in_SF1.push_back(*current);
            else
                in_both.push_back(*current);
        current++;
    }

    // Search in File 2
    current = SF2.text_line.begin();
    last    = SF2.text_line.end();

    while (current != last)
    {
        // Skip if not active
        if ((*current).line_type != SelLine::DATALINE)
        {
            current++;
            continue;
        }
        if ((*current).label == SelLine::DISCARDED)
        {
            SF2_discarded++;
            current++;
            continue;
        }

        // Try to find this archive into Sel File 2
        found = find(in_both, (*current).text);
        if (found != in_both.end())
        {
            current++;
            continue;
        }
        only_in_SF2.push_back(*current);
        current++;
    }

    // Write Statistics
    if (mode < 0)
    {
	temp.line_type = SelLine::COMMENT;
	temp.label = SelLine::DISCARDED;
	temp.text = "# Statistics of comparison";
	result.text_line.push_back(temp);
	temp.text = "# -------------------------------------------------------------";
	result.text_line.push_back(temp);
	sprintf(str, "%6d", SF1.no_imgs);
	temp.text = "# File 1: " + SF1.fn_sel + "(VALID: " + str;
	sprintf(str, "%6d", SF1_discarded);
	temp.text += (std::string) " DISCARDED: " + str + ")";
	result.text_line.push_back(temp);
	sprintf(str, "%6d", SF2.no_imgs);
	temp.text = "# File 2: " + SF2.fn_sel + "(VALID: " + str;
	sprintf(str, "%6d", SF2_discarded);
	temp.text += (std::string) " DISCARDED: " + str + ")";
	result.text_line.push_back(temp);
	temp.text = "";
	result.text_line.push_back(temp);
	sprintf(str, "%6lu", (unsigned long int)in_both.size());
	temp.text = (std::string)"# Matching Files: " + str;
	result.text_line.push_back(temp);
	sprintf(str, "%6lu", (unsigned long int)only_in_SF1.size());
	temp.text = (std::string)"# Only in file 1: " + str;
	result.text_line.push_back(temp);
	sprintf(str, "%6lu", (unsigned long int)only_in_SF2.size());
	temp.text = (std::string)"# Only in file 2: " + str;
	result.text_line.push_back(temp);
	temp.text = "# -------------------------------------------------------------";
	result.text_line.push_back(temp);

	// Write files in both
	temp.text = "";
	result.text_line.push_back(temp);
	temp.text = "# Files in both .sel files";
	result.text_line.push_back(temp);
    }
    if (mode<0 || mode==0) 
    {
	current = in_both.begin();
	last    = in_both.end();
	while (current != last)
	    result.text_line.push_back(*current++);
    }

    if (mode<0)
    {
	// Write files only in Sel File 1
	temp.text = "";
	result.text_line.push_back(temp);
	temp.text = "# Files only in the first file";
	result.text_line.push_back(temp);
    }
    if (mode<0 || mode==1) 
    {
	current = only_in_SF1.begin();
	last    = only_in_SF1.end();
	while (current != last)
	    result.text_line.push_back(*current++);
    }

    if (mode<0)
    {
	// Write files only in Sel File 2
	temp.text = "";
	result.text_line.push_back(temp);
	temp.text = "# Files only in the second file";
	result.text_line.push_back(temp);
    }
    if (mode<0 || mode==2) 
    {
	current = only_in_SF2.begin();
	last    = only_in_SF2.end();
	while (current != last)
	    result.text_line.push_back(*current++);
    }
    // Adjust the remaining fields
    if (mode<0) 
	result.no_imgs = in_both.size() + only_in_SF1.size() + only_in_SF2.size();
    else if (mode==0) 
	result.no_imgs = in_both.size();
    else if (mode==1)
	result.no_imgs = only_in_SF1.size();
    else if (mode==2)
	result.no_imgs = only_in_SF2.size();
    result.current_line = result.text_line.begin();

    return result;
}
