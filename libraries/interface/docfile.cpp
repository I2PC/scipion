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

#include <fstream>
#include <sstream>
#include <cstdio>

#include "docfile.h"
#include <data/args.h>

DocLine::DocLine(const DocLine& line)
{
    line_type = line.line_type;
    text = line.text;
    key = line.key;
    data = line.data;
}

DocLine& DocLine::operator=(const DocLine& line)
{
    if (this != &line)
    {
        line_type = line.line_type;
        text = line.text;
        key = line.key;
        data = line.data;
    }

    return *this;
}

double& DocLine::operator[](size_t i)
{
    if (i+1 > data.size())
    {
        std::cout << "Docline=" << *this << std::endl;
        REPORT_ERROR(ERR_DOCFILE, "Trying to access to non-existing element " +
                     integerToString(i) + " of a document line");
    }

    return data[i];
}

double DocLine::operator[](size_t i) const
{
    if (i+1 > data.size())
    {
        std::cout << "Docline=" << *this << std::endl;
        REPORT_ERROR(ERR_DOCFILE, "Trying to access to non-existing element " +
                     integerToString(i) + " of a document line");
    }

    return data[i];
}

void DocLine::set(size_t i, double val)
{
    // Make sure there is enough memory
    if (i + 1 > data.size())
        data.reserve(i + 1);

    // Pad with zeros for the non-existing indexes in between
    int space_needed = i + 1 - data.size();
    for (int k = 0; k < space_needed; k++)
        data.push_back(0);

    // Set required data
    data[i] = val;
}

void DocLine::set(const Matrix1D< double >& v)
{
    data.clear();
    if (line_type != DATALINE)
    {
        line_type = DATALINE;
        key = 0;
    }

    data.reserve(VEC_XSIZE(v));
    FOR_ALL_ELEMENTS_IN_MATRIX1D(v)
    data.push_back(VEC_ELEM(v, i));
}

void DocLine::clear()
{
    line_type = NOT_ASSIGNED;
    text = "";
    key = 0;
    data.clear();
}

std::ostream& operator<<(std::ostream& o, const DocLine& line)
{
    char aux[30];
    switch (line.line_type)
    {
    case (DocLine::DATALINE):
                    // Print a data line
                    sprintf(aux, "%5d ", line.key);
        o << aux;
        sprintf(aux, "%-2lu", (unsigned long int)(line.data.size()));
        o << aux;

        int imax;
        imax = line.data.size();
        for (int i = 0; i < imax; i++)
{
            sprintf(aux, " % 10.5f", line.data[i]);
            o << aux;
        }

        o << std::endl;
        break;

    case (DocLine::COMMENT):
                    // Print a comment
                    o << line.text << std::endl;
        break;
    default:
    	break;
    }

    return o;
}

void DocLine::read(std::istream& in)
{
    std::string line;
    int param_no;

    // Get line
    getline(in, line);

    // Initialize target
    line_type = DocLine::NOT_ASSIGNED;
    text = "";
    key = 0;
    data.clear();

    // Check if comment or empty line
    int charpos1 = line.find_first_not_of(" \t");
    if (line[0] == '\0' || line[charpos1] == '#' || line[charpos1] == ';')
    {
        line_type = DocLine::COMMENT;
        text = line;
        data.clear();
        key = 0;
    }
    // Read a true document file line
    else
    {
        line_type = DocLine::DATALINE;
        text = "";
        size_t i = 0;

        key = textToInteger(nextToken(line, i));
        param_no = textToInteger(nextToken(line, i));
        std::string auxline = line;

        try
        {
            // Try unfixed mode first
            readFloatList(line, i, param_no, data);
        }
        catch (XmippError e)
        {
            // Try fixed mode then
            data.clear();
            data.reserve(param_no);
            for (int i = 0; i < param_no; i++)
            {
                data.push_back(textToFloat(line.substr(8 + i*12, 12)));
            }
        }
    }
}

DocFile::DocFile(const DocFile& doc)
{
    fn_doc = doc.fn_doc;
    m = doc.m;
    no_lines = doc.no_lines;
    first_key = doc.first_key;
    current_line = doc.current_line;
}

void DocFile::clear()
{
    fn_doc = "";
    m.clear();
    no_lines = 0;
    first_key = 1;
    current_line = m.begin();
}

DocFile& DocFile::operator=(const DocFile& doc)
{
    if (this != &doc)
    {
        fn_doc = doc.fn_doc;
        m = doc.m;
        no_lines = doc.no_lines;
        first_key = doc.first_key;
        current_line = doc.current_line;
    }

    return *this;
}

DocFile& DocFile::operator=(const Matrix2D< double >& A)
{
    clear();
    DocLine temp;

    for (size_t i = 0; i <MAT_YSIZE(A); i++)
    {
        temp.clear();
        temp.line_type = DocLine::DATALINE;
        temp.data.resize(MAT_XSIZE(A));

        for (size_t j = 0; j < MAT_XSIZE(A); j++)
            temp.data[j] = MAT_ELEM(A, i, j);

        m.push_back(temp);
    }

    fn_doc = "";
    no_lines = MAT_YSIZE(A);
    renum();
    go_beginning();

    return *this;
}

//NOTE: MPI PROGRAMS USE << TP PASS DOCFILES
//FROM MASTER TO SLAVE
//PLEASE DO NOT ALTER THE OUTPUR FORMAT
std::ostream& operator<<(std::ostream& o, const DocFile& doc)
{
    std::vector< DocLine >::const_iterator current = doc.m.begin();
    std::vector< DocLine >::const_iterator last = doc.m.end();

    while (current != last)
    {
        o << *current;
        current++;
    }

    return o;
}

void DocFile::show_line(std::ostream& o, int key)
{
    if (key == -1)
    {
        if (current_line == m.end())
            o << "Current line is at the end of file\n";
        else
            o << *current_line;
    }
    else
    {
        std::vector< DocLine >::iterator line = find(key);
        if (line == m.end())
            o << "Key " << key << " not found\n";
        else
            o << *line;
    }
}

void DocFile::debug()
{
    std::vector< DocLine >::iterator current = m.begin();
    std::vector< DocLine >::iterator last = m.end();

    while (current != last)
    {
        if ((*current).line_type == DocLine::DATALINE ||
            (*current).line_type == DocLine::COMMENT)
            std::cout << *current;
        else
        {
            char aux[30];
            std::string str = "";

            std::cout << "Special line\n";
            std::cout << "  Type: " << (*current).line_type << std::endl;
            std::cout << "  Key:  " << (*current).key << std::endl;
            std::cout << "  Text: " << (*current).text << std::endl;
            std::cout << "  Data: ";
            for (size_t i = 0; i < (*current).data.size(); i++)
            {
                sprintf(aux, " % 11.5f", (*current).data[i]);
                str += aux;
            }
            std::cout << str << std::endl;
        }

        current++;
    }
}

std::vector< DocLine >::iterator DocFile::find(int k)
{
    std::vector< DocLine >::iterator current = m.begin();
    std::vector< DocLine >::iterator last = m.end();

    while (current != last)
    {
        if ((*current).line_type == DocLine::DATALINE && (*current).key == k)
            return current;

        current++;
    }

    return current;
}

void DocFile::adjust_to_data_line()
{
    if (current_line == m.end())
        return;

    while ((*current_line).line_type != DocLine::DATALINE)
    {
        current_line++;

        if (current_line == m.end())
            return;
    }
}

void DocFile::renum()
{
    std::vector< DocLine >::iterator current = m.begin();
    std::vector< DocLine >::iterator last = m.end();
    int act_key = first_key;

    while (current != last)
    {
        if ((*current).line_type == DocLine::DATALINE)
            (*current).key = act_key++;

        current++;
    }
}

void DocFile::read(const FileName& name, int overriding)
{
    DocLine temp;
    std::ifstream in;
    int line_no = 1;

    // Empties current DocFile
    if (overriding)
        clear();

    // Open file
    in.open(name.c_str(), std::ios::in);
    if (!in)
        REPORT_ERROR(ERR_DOCFILE, "DocFile::read: File " + name + " not found");

    // Read each line and keep it in the list of the DocFile object
    in.peek();
    while (!in.eof())
    {
        try
        {
            temp.read(in);
        }
        catch (XmippError e)
        {
            std::cout << "Doc File: Line " << line_no <<
            " is skipped due to an error\n";
        }

        switch (temp.line_type)
        {
        case (DocLine::NOT_ASSIGNED):
                        break; // Line with an error

        case (DocLine::DATALINE):
                        no_lines++;
            m.push_back(temp);
            break;

        case (DocLine::COMMENT):
                        m.push_back(temp);
            break;
        default:
        	break;
        }

        line_no++;
        in.peek();
    }

    // Close file
    in.close();

    // Set "pointer" to the beginning of the file
    if (overriding)
        fn_doc = name;

    go_first_data_line();
    renum();
}

void DocFile::write(const FileName& name)
{
    std::ofstream out;
    std::vector< DocLine >::iterator current = m.begin();
    std::vector< DocLine >::iterator last = m.end();

    if (name != "")
        fn_doc = name;

    renum();

    // Open file
    out.open(fn_doc.c_str(), std::ios::out);
    if (!out)
        REPORT_ERROR(ERR_IO_NOWRITE, "DocFile::write: File " + fn_doc +
                     " cannot be written");

    // Read each line and keep it in the list of the SelFile object
    while (current != last)
        out << *(current++);

    // Close file
    out.close();
}

void DocFile::jump(int count)
{
    adjust_to_data_line();

    for (int i = 0; i < count; i++)
        if (current_line != m.end())
        {
            current_line++;
            adjust_to_data_line();
        }
}

int DocFile::search_comment(std::string comment)
{
    go_beginning();
    comment = " ; " + comment;

    while (!eof())
    {
        if ((*current_line).Is_comment())
        {
            if (strcmp(comment.c_str(), ((*current_line).get_text()).c_str())
                == 0)
            {
                adjust_to_data_line();
                return 1;
            }
        }

        next();
    }

    return 0;
}

int DocFile::remove_multiple_strings(std::string pattern)
{
    bool found_once = false;
    go_beginning();
    while (!eof())
    {
        if ((*current_line).Is_comment())
        {
            if (((*current_line).get_text()).find(pattern) <
                ((*current_line).get_text()).size() )
            {
                if (found_once)
                {
                    remove_current();
                }
                found_once = true;
            }
        }
        next();
    }

    return 0;

}


void DocFile::get_selfile(MetaData& sel)
{
    go_beginning();

    if ((*current_line).Is_comment())
        if (strstr(((*current_line).get_text()).c_str(), "Headerinfo") == NULL)
            REPORT_ERROR(ERR_DOCFILE,
                         "DocFile::get_selfile: Docfile is of non-NewXmipp type!");

    sel.clear();
    next();
    FileName img;
    while (!eof())
    {
        if (strstr(((*current_line).get_text()).c_str(), " ; ") != NULL)
        {
            img = (*current_line).get_text();
            sel.setValue(MDL_IMAGE,img.removeSubstring(" ; "), sel.addObject());
        }

        next();
    }
}

void DocFile::locate(int k)
{
    std::vector< DocLine >::iterator last = m.end();

    current_line = m.begin();
    while (current_line != last)
    {
        if ((*current_line).line_type == DocLine::DATALINE &&
            (*current_line).key >= k)
            return;

        current_line++;
    }
}

int DocFile::getColNumberFromHeader(const char * pattern)
{
    std::string header;
    go_beginning();
    if ((*current_line).Is_comment())
    {
        header = (*current_line).get_text();
        if (strstr(header.c_str(), "Headerinfo") == NULL)
            REPORT_ERROR(ERR_DOCFILE,"DocFile:: docfile is of non-NewXmipp type!");
        else
        {
            std::vector<std::string> tokens;
            tokenize(header,tokens," \t()");
            for (size_t i = 0; i < tokens.size(); i++)
            {
                if (strstr(tokens[i].c_str(), pattern) != NULL)
                {
                    return textToInteger(tokens[i+1]);
                }
            }
        }
    }
    return -1;

}


int DocFile::FirstLine_colNumber()
{
    std::vector< DocLine >::iterator aux = current_line;
    go_first_data_line();

    int ret = current_line->data.size();
    current_line = aux;

    return ret;
}

int DocFile::get_last_key()
{
    std::vector< DocLine >::iterator last = m.end();
    std::vector< DocLine >::iterator first = m.begin();

    do
    {
        if (last == first)
            return first_key - 1;

        last--;
        if (last->line_type == DocLine::DATALINE)
            return last->key;

    }
    while (true);
}

double DocFile::operator()(int k, int i)
{
    std::vector< DocLine >::iterator aux = find(k);

    if (aux == m.end())
        REPORT_ERROR(ERR_DOCFILE, "DocFile::operator(): The given key (" + integerToString(k)
                     + ") is not in the file");

    return (*aux)[i];
}

void DocFile::get_angles(int k, double& rot, double& tilt, double& psi,
                         const std::string& ang1, const std::string& ang2, const std::string& ang3)
{
    std::vector< DocLine >::iterator it = find(k);

    if (it == m.end())
        REPORT_ERROR(ERR_VALUE_INCORRECT,
                     "DocFile::get_angles(): The given key (" + integerToString(k)
                     + ") is not in the file");

    switch (ang1[0])
    {
    case 'r':
        rot = (*it)[0];
        break;

    case 't':
        tilt = (*it)[0];
        break;

    case 'p':
        psi = (*it)[0];
        break;
    }

    switch (ang2[0])
    {
    case 'r':
        rot = (*it)[1];
        break;

    case 't':
        tilt = (*it)[1];
        break;

    case 'p':
        psi = (*it)[1];
        break;
    }

    switch (ang3[0])
    {
    case 'r':
        rot = (*it)[2];
        break;

    case 't':
        tilt = (*it)[2];
        break;

    case 'p':
        psi = (*it)[2];
        break;
    }
}

void DocFile::get_angles1(int k, double& rot, double& tilt, double& psi,
                          const std::string& ang1, const std::string& ang2, const std::string& ang3)
{
    std::vector< DocLine >::iterator it = find(k);

    if (it == m.end())
        REPORT_ERROR(ERR_DOCFILE, "DocFile::get_angles1(): The given key (" +
                     integerToString(k) + ") is not in the file");

    switch (ang1[0])
    {
    case 'r':
        rot = (*it)[3];
        break;

    case 't':
        tilt = (*it)[3];
        break;

    case 'p':
        psi = (*it)[3];
        break;
    }

    switch (ang2[0])
    {
    case 'r':
        rot = (*it)[4];
        break;

    case 't':
        tilt = (*it)[4];
        break;

    case 'p':
        psi = (*it)[4];
        break;
    }

    switch (ang3[0])
    {
    case 'r':
        rot = (*it)[5];
        break;

    case 't':
        tilt = (*it)[5];
        break;

    case 'p':
        psi = (*it)[5];
        break;
    }
}

void DocFile::get_angles2(int k, double& rot, double& tilt, double &psi,
                          const std::string& ang1, const std::string& ang2, const std::string& ang3)
{
    std::vector< DocLine >::iterator it = find(k);
    if (it == m.end())
        REPORT_ERROR(ERR_DOCFILE, "DocFile::get_angles2(): The given key (" + integerToString(k)
                     + ") is not in the file");

    switch (ang1[0])
    {
    case 'r':
        rot = (*it)[6];
        break;

    case 't':
        tilt = (*it)[6];
        break;

    case 'p':
        psi = (*it)[6];
        break;
    }

    switch (ang2[0])
    {
    case 'r':
        rot = (*it)[7];
        break;

    case 't':
        tilt = (*it)[7];
        break;

    case 'p':
        psi = (*it)[7];
        break;
    }

    switch (ang3[0])
    {
    case 'r':
        rot = (*it)[8];
        break;

    case 't':
        tilt = (*it)[8];
        break;

    case 'p':
        psi = (*it)[8];
        break;
    }
}

void DocFile::set_angles(int k, double rot, double tilt, double psi,
                         const std::string& ang1, const std::string& ang2, const std::string& ang3)
{
    std::vector< DocLine >::iterator it = find(k);
    if (it == m.end())
        REPORT_ERROR(ERR_DOCFILE, "DocFile::set_angles(): The given key (" +
                     integerToString(k) + ") is not in the file");

    switch (ang1[0])
    {
    case 'r':
        (*it)[0] = rot;
        break;

    case 't':
        (*it)[0] = tilt;
        break;

    case 'p':
        (*it)[0] = psi;
        break;
    }

    switch (ang2[0])
    {
    case 'r':
        (*it)[1] = rot;
        break;

    case 't':
        (*it)[1] = tilt;
        break;

    case 'p':
        (*it)[1] = psi;
        break;
    }

    switch (ang3[0])
    {
    case 'r':
        (*it)[2] = rot;
        break;

    case 't':
        (*it)[2] = tilt;
        break;

    case 'p':
        (*it)[2] = psi;
        break;
    }
}

void DocFile::get_image(int key, Image<double> &I, bool apply_geo)
{
    current_line = m.begin();
    current_line += 2*key;
    DocLine DL = get_current_line();
    previous();
    FileName fn_img;
    if (get_current_line().Is_comment())
        fn_img = ((get_current_line()).get_text()).erase(0, 3);
    else
        REPORT_ERROR(ERR_DOCFILE,"The docfile provided is not of type Alignment");

    // Read actual image
    I.read(fn_img);
    I().setXmippOrigin();

    // Store translation in header and apply it to the actual image
    I.setRot(DL[0]);
    I.setTilt(DL[1]);
    I.setPsi(DL[2]);
    I.setXoff(DL[3]);
    I.setYoff(DL[4]);
    I.setFlip(0.);

    if (apply_geo)
    {
        Matrix2D<double> A;
        I.getTransformationMatrix(A, true);
        if (!A.isIdentity())
            selfApplyGeometry(BSPLINE3, I(), A, IS_INV, WRAP);
    }
}

FileName DocFile::get_imagename(int key)
{
    current_line = m.begin();
    current_line += 2*key;
    DocLine DL = get_current_line();
    previous();
    FileName fn_img;
    if (get_current_line().Is_comment())
        fn_img = ((get_current_line()).get_text()).erase(0, 3);
    else
        REPORT_ERROR(ERR_DOCFILE,"The docfile provided is not of type Alignment");
    return fn_img;
}

void DocFile::set(int i, double val)
{
    if (current_line != m.end())
        (*current_line).set(i, val);
    else
    {
        // Add a new line
        DocLine tmp;
        tmp.line_type = DocLine::DATALINE;
        tmp.set(i, val);
        m.push_back(tmp);
        current_line = m.end();
    }
}

void DocFile::set(int k, int i, double val)
{
    std::vector< DocLine >::iterator it = find(k);
    if (it == m.end())
        REPORT_ERROR(ERR_DOCFILE, "DocFile::set(): The given key (" + integerToString(k)
                     + ") is not in the file");

    it->set(i, val);
}

void DocFile::remove(int key0, int keyF)
{
    std::vector< DocLine >::iterator last = m.end();

    if (keyF == -1)
        keyF = key0;

    int old = (*current_line).key;

    current_line = m.begin();
    while (current_line != last)
    {
        if (current_line->line_type != DocLine::DATALINE)
            continue;

        if (current_line->key >= key0 && current_line->key <= keyF)
            remove_current();
        else
            current_line++;
    }

    locate(old);
}

void DocFile::remove_current()
{
    if (current_line != m.end())
    {
        std::vector< DocLine >::iterator it = current_line;
        it++;

        if (current_line->line_type == DocLine::DATALINE)
            no_lines--;

        m.erase(current_line);
        current_line = it;
    }
}

int DocFile::insert_data_line(int count)
{
    DocLine tmp;
    tmp.line_type = DocLine::DATALINE;

    std::vector< DocLine >::iterator it;
    for (int i = 0; i < count; i++)
    {
        current_line = m.insert(current_line, tmp);

        if (i == 0)
            it = current_line;

        current_line++;
    }

    renum();
    no_lines += count;

    return it->key;
}

int DocFile::insert_data_line(const Matrix1D< double >& v)
{
    DocLine tmp;
    tmp.set(v);

    current_line = m.insert(current_line, tmp);
    renum();

    int ret = current_line->key;
    current_line++;
    no_lines++;

    return ret;
}

void DocFile::insert_comment(std::string comment)
{
    DocLine tmp;
    tmp.line_type = DocLine::COMMENT;
    tmp.text = " ; " + comment;

    // Insert and update current_line
    current_line = m.insert(current_line, tmp);
    current_line++;
}

int DocFile::insert_line(const DocLine& line)
{
    int ret = 0;
    current_line = m.insert(current_line, line);

    if (line.line_type == DocLine::DATALINE)
    {
        renum();
        ret = current_line->key;
        no_lines++;
    }

    current_line++;
    return ret;
}

int DocFile::append_data_line(int count)
{
    DocLine tmp;
    tmp.line_type = DocLine::DATALINE;

    int ret = get_last_key() + 1;
    int act_key = ret;

    for (int i = 0; i < count; i++, act_key++)
    {
        tmp.key = act_key;
        m.push_back(tmp);
    }

    no_lines += count;

    return ret;
}

int DocFile::append_data_line(const Matrix1D< double >& v)
{
    DocLine tmp;
    tmp.set(v);
    tmp.key = get_last_key() + 1;
    m.push_back(tmp);
    no_lines++;

    return tmp.key;
}

void DocFile::append_comment(const std::string& comment)
{
    DocLine tmp;
    tmp.line_type = DocLine::COMMENT;

    tmp.text = " ; " + comment;
    m.push_back(tmp);
}

int DocFile::append_angles(double rot, double tilt, double psi,
                           const std::string& ang1, const std::string& ang2, const std::string& ang3)
{
    Matrix1D< double > aux(3);

    if (ang1[0] == 'r')
        VEC_ELEM(aux, 0) = rot;
    else if (ang1[0] == 't')
        VEC_ELEM(aux, 0) = tilt;
    else if (ang1[0] == 'p')
        VEC_ELEM(aux, 0) = psi;

    if (ang2[0] == 'r')
        VEC_ELEM(aux, 1) = rot;
    else if (ang2[0] == 't')
        VEC_ELEM(aux, 1) = tilt;
    else if (ang2[0] == 'p')
        VEC_ELEM(aux, 1) = psi;

    if (ang3[0] == 'r')
        VEC_ELEM(aux, 2) = rot;
    else if (ang3[0] == 't')
        VEC_ELEM(aux, 2) = tilt;
    else if (ang3[0] == 'p')
        VEC_ELEM(aux, 2) = psi;

    return append_data_line(aux);
}

int DocFile::append_angles(double rot, double tilt, double psi,
                           double rot1, double tilt1, double psi1,
                           const std::string& ang1, const std::string& ang2, const std::string& ang3)
{
    Matrix1D< double > aux(6);

    if (ang1[0] == 'r')
        VEC_ELEM(aux, 0) = rot;
    else if (ang1[0] == 't')
        VEC_ELEM(aux, 0) = tilt;
    else if (ang1[0] == 'p')
        VEC_ELEM(aux, 0) = psi;

    if (ang2[0] == 'r')
        VEC_ELEM(aux, 1) = rot;
    else if (ang2[0] == 't')
        VEC_ELEM(aux, 1) = tilt;
    else if (ang2[0] == 'p')
        VEC_ELEM(aux, 1) = psi;

    if (ang3[0] == 'r')
        VEC_ELEM(aux, 2) = rot;
    else if (ang3[0] == 't')
        VEC_ELEM(aux, 2) = tilt;
    else if (ang3[0] == 'p')
        VEC_ELEM(aux, 2) = psi;

    if (ang1[0] == 'r')
        VEC_ELEM(aux, 3) = rot1;
    else if (ang1[0] == 't')
        VEC_ELEM(aux, 3) = tilt1;
    else if (ang1[0] == 'p')
        VEC_ELEM(aux, 3) = psi1;

    if (ang2[0] == 'r')
        VEC_ELEM(aux, 4) = rot1;
    else if (ang2[0] == 't')
        VEC_ELEM(aux, 4) = tilt1;
    else if (ang2[0] == 'p')
        VEC_ELEM(aux, 4) = psi1;

    if (ang3[0] == 'r')
        VEC_ELEM(aux, 5) = rot1;
    else if (ang3[0] == 't')
        VEC_ELEM(aux, 5) = tilt1;
    else if (ang3[0] == 'p')
        VEC_ELEM(aux, 5) = psi1;

    return append_data_line(aux);
}

int DocFile::append_angles(double rot, double tilt, double psi,
                           double rot1, double tilt1, double psi1,
                           double rot2, double tilt2, double psi2,
                           const std::string& ang1, const std::string& ang2, const std::string& ang3)
{
    Matrix1D< double > aux(9);

    if (ang1[0] == 'r')
        VEC_ELEM(aux, 0) = rot;
    else if (ang1[0] == 't')
        VEC_ELEM(aux, 0) = tilt;
    else if (ang1[0] == 'p')
        VEC_ELEM(aux, 0) = psi;

    if (ang2[0] == 'r')
        VEC_ELEM(aux, 1) = rot;
    else if (ang2[0] == 't')
        VEC_ELEM(aux, 1) = tilt;
    else if (ang2[0] == 'p')
        VEC_ELEM(aux, 1) = psi;

    if (ang3[0] == 'r')
        VEC_ELEM(aux, 2) = rot;
    else if (ang3[0] == 't')
        VEC_ELEM(aux, 2) = tilt;
    else if (ang3[0] == 'p')
        VEC_ELEM(aux, 2) = psi;

    if (ang1[0] == 'r')
        VEC_ELEM(aux, 3) = rot1;
    else if (ang1[0] == 't')
        VEC_ELEM(aux, 3) = tilt1;
    else if (ang1[0] == 'p')
        VEC_ELEM(aux, 3) = psi1;

    if (ang2[0] == 'r')
        VEC_ELEM(aux, 4) = rot1;
    else if (ang2[0] == 't')
        VEC_ELEM(aux, 4) = tilt1;
    else if (ang2[0] == 'p')
        VEC_ELEM(aux, 4) = psi1;

    if (ang3[0] == 'r')
        VEC_ELEM(aux, 5) = rot1;
    else if (ang3[0] == 't')
        VEC_ELEM(aux, 5) = tilt1;
    else if (ang3[0] == 'p')
        VEC_ELEM(aux, 5) = psi1;

    if (ang1[0] == 'r')
        VEC_ELEM(aux, 6) = rot2;
    else if (ang1[0] == 't')
        VEC_ELEM(aux, 6) = tilt2;
    else if (ang1[0] == 'p')
        VEC_ELEM(aux, 6) = psi2;

    if (ang2[0] == 'r')
        VEC_ELEM(aux, 7) = rot2;
    else if (ang2[0] == 't')
        VEC_ELEM(aux, 7) = tilt2;
    else if (ang2[0] == 'p')
        VEC_ELEM(aux, 7) = psi2;

    if (ang3[0] == 'r')
        VEC_ELEM(aux, 8) = rot2;
    else if (ang3[0] == 't')
        VEC_ELEM(aux, 8) = tilt2;
    else if (ang3[0] == 'p')
        VEC_ELEM(aux, 8) = psi2;

    return append_data_line(aux);
}

int DocFile::append_line(DocLine& line)
{
    int ret = 0;
    if (line.line_type == DocLine::DATALINE)
    {
        no_lines++;
        ret = line.key = get_last_key() + 1;
    }

    m.push_back(line);

    return ret;
}

void DocFile::clean_comments()
{
    std::vector< DocLine >::iterator last = m.end();
    current_line = m.begin();

    while (current_line != last)
    {
        if (current_line->line_type == DocLine::COMMENT)
            remove_current();
        else
            current_line++;
    }

    current_line = m.begin();
}

DocFile DocFile::randomize()
{
    DocFile result, aux;
    int i;
    int rnd_indx;

    randomize_random_generator();
    if (no_lines == 0)
        return aux;

    aux = *this;

    for (i = no_lines; i > 0; i--)
    {
        // Jump a random number from the beginning
        rnd_indx = (int) rnd_unif(0, i);
        aux.go_first_data_line();
        aux.jump(rnd_indx);

        result.m.push_back(*aux.current_line);
        aux.current_line->line_type = DocLine::NOT_CONSIDERED;
    }

    // Adjust remaining fields
    result.no_lines = no_lines;
    result.current_line = result.m.begin();
    result.renum();

    return result;
}

void DocFile::perturb_column(int col, double sigma)
{

    randomize_random_generator();

    std::vector< DocLine >::iterator current = m.begin();
    std::vector< DocLine >::iterator last = m.end();

    while (current != last)
    {
        if (current->line_type == DocLine::DATALINE)
            (*current)[col] += rnd_gaus(0., sigma);
        current++;
    }

}

//merge
void DocFile::merge(const FileName& name, int mode, int sumcol)
{
    DocFile DFaux(name);
    merge(DFaux, mode, sumcol);
}

//merge
void DocFile::merge(DocFile& DF, int mode, int sumcol)
{

    DocLine  DL, DLold;
    MetaData  SF;
    FileName fn_img;
    double   w;

    DF.get_selfile(SF);
    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        SF.getValue(MDL_IMAGE,fn_img, __iter.objId);
        DF.search_comment(fn_img);
        DL=DF.get_current_line();
        if (search_comment(fn_img))
        {
            switch (mode)
            {
            case(DOCMERGE_KEEP_OLD):
                {
                    // just keep what's there and do nothing
                    break;
                }
            case(DOCMERGE_KEEP_NEW):
                {
                    //Replace current data line with the new one
                    (*current_line) = DL;
                    break;
                }
            case(DOCMERGE_SUM_COLUMN):
                {
                    // Just sum column
                    w = (*current_line).data[sumcol] + DF(sumcol);
                    (*current_line).set(sumcol, w);
                    break;
                }
            case(DOCMERGE_ERROR):
                            std::cerr<<"image name = "<<fn_img;
                REPORT_ERROR(ERR_DOCFILE,"Image occured in two docfiles to be merged");
            }
        }
        else
{
            // Just add new line at the end of the file
            current_line = m.end();
            append_comment(fn_img);
            append_line(DL);
        }
    }
}

DocFile DocFile::random_discard(int n)
{
    DocFile result;
    int i, rnd_indx;

    result = *this;
    randomize_random_generator();
    n = std::min(n, no_lines);

    for (i = 0; i < n; i++)
    {
        // Jump a random number from the beginning
        rnd_indx = (int) rnd_unif(0, result.no_lines);
        result.go_first_data_line();
        result.jump(rnd_indx);
        result.remove_current();
    }
    result.go_beginning();

    return result;
}

Matrix1D< double > DocFile::col(int c)
{
    Matrix1D< double > result(no_lines);

    std::vector< DocLine >::iterator current = m.begin();
    std::vector< DocLine >::iterator last = m.end();
    int i = 0;

    while (current != last)
    {
        if (current->line_type == DocLine::DATALINE)
            VEC_ELEM(result, i++) = (*current)[c];

        current++;
    }

    return result;
}

Matrix1D< double > DocFile::row(int k)
{
    Matrix1D< double > result;
    std::vector< DocLine >::iterator it = find(k);
    if (it == m.end())
        return result;

    result.resize(it->data.size());
    result.setRow();

    for (size_t i = 0; i < VEC_XSIZE(result); i++)
        VEC_ELEM(result, i) = it->data[i];

    return result;
}

void DocFile::setCol(int c, Matrix1D< double >& v)
{
    go_first_data_line();

    for (size_t i = 0; i <VEC_XSIZE(v); i++)
    {
        set(c, VEC_ELEM(v, i));
        next_data_line();
    }
    renum();

    if (VEC_XSIZE(v) < m.size())
        REPORT_ERROR(ERR_DOCFILE, "DocFile::setCol(): Column assignment not complete");
}

int read_Euler_document_file(FileName name, std::string ang1, std::string ang2,
                             std::string ang3, DocFile& doc)
{
    DocFile aux(name);
    DocLine line1, line2;

    // Set line2 as a data line and go to the beginning of file
    line2.set_type(DocLine::DATALINE);
    aux.go_beginning();

    // Macro to assign the angle from line1 in the correct place of line2
    // The angle order in line2 is (rot, tilt, psi)
#define assign_in_correct_place_of_line2(angle_descr,angle_index) \
    switch (angle_descr[0]) \
    { \
    case ('r'): line2.set(0, line1[angle_index]); break; \
    case ('t'): line2.set(1, line1[angle_index]); break; \
    case ('p'): line2.set(2, line1[angle_index]); break; \
    }

    // Read the whole file
    while (!aux.eof())
    {
        line1 = aux.get_current_line();

        // If line1 type is neither DATALINE nor COMMENT the line is skipped!!
        if (line1.Is_data())
        {
            // Reorder the angles and insert
            assign_in_correct_place_of_line2(ang1, 0);
            assign_in_correct_place_of_line2(ang2, 1);
            assign_in_correct_place_of_line2(ang3, 2);

            doc.append_line(line2);

        }
        else if (line1.Is_comment())
            // Insert the comment
            doc.append_line(line1);

        // Next line
        aux.next();
    }

    return doc.dataLineNo();
}

void select_images(DocFile& doc, MetaData& sel, int col, bool en_limit0,
                   double limit0, bool en_limitF, double limitF)
{
    doc.go_first_data_line();
    int enabled;
    FOR_ALL_OBJECTS_IN_METADATA(sel)
    {
        sel.getValue(MDL_ENABLED,enabled, __iter.objId);
        if (enabled == 1)
        {
            if ((en_limit0 && doc(col) < limit0) || (en_limitF && doc(col) > limitF))
                sel.setValue(MDL_ENABLED, -1, __iter.objId);
        }

        doc.next_data_line();
    }
}

void get_subset_docfile(DocFile& DFin, MetaData& SF, DocFile& DFout)
{
    DocLine DL;
    FileName fn_tmp;

    DFout.clear();
    DFin.go_beginning();
    DL = DFin.get_current_line();
    if (DL.Is_comment())
        fn_tmp = DL.get_text();
    if (strstr(fn_tmp.c_str(), "Headerinfo") == NULL)
        REPORT_ERROR(ERR_DOCFILE,"Input docfile is not of NewXmipp-style");
    else
        // append the same header to DFout
        DFout.append_comment(fn_tmp.removeSubstring(" ; "));

    FOR_ALL_OBJECTS_IN_METADATA(SF)
    {
        SF.getValue(MDL_IMAGE,fn_tmp, __iter.objId);
        if (DFin.search_comment(fn_tmp))
        {
            DFout.append_comment(fn_tmp);
            DL=DFin.get_current_line();
            DFout.append_line(DL);
        }
        else
            REPORT_ERROR(ERR_DOCFILE, (std::string)"Docfile: Cannot find " + fn_tmp + " in docfile ");
    }
}

