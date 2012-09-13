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
#include "funcs.h"
#include "args.h"

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <complex>
#include <fstream>
#include <typeinfo>

using namespace std;

/* Numerical functions ----------------------------------------------------- */
// Kaiser-Bessel constructor
KaiserBessel::KaiserBessel(double alpha_, int K_, double r_, double v_,
			   int N_, double vtable_, int ntable_) 
    : alpha(alpha_), v(v_), r(r_), N(N_), K(K_), vtable(vtable_), 
      ntable(ntable_) 
{
    // Default values are alpha=1.25, K=6, r=0.5, v = K/2
    if (0.f == v) v = double(K)/2;
    if (0.f == vtable) vtable = v;
    alphar = alpha*r;
    fac = static_cast<double>(2.*PI)*alphar*v;
    vadjust = 1.0f*v;
    facadj = static_cast<double>(2.*PI)*alphar*vadjust;
    build_I0table();
}

// Kaiser-Bessel I0 window function
double KaiserBessel::i0win(double x) const 
{
    double val0 = double(bessi0(facadj));
    double absx = fabs(x);
    if (absx > vadjust) return 0.f;
    double rt = sqrt(1.f - pow(absx/vadjust, 2));
    double res = bessi0(facadj*rt)/val0;
    return res;
}

// Tabulate I0 window for speed
void KaiserBessel::build_I0table() 
{
    i0table.resize(ntable+1); // i0table[0:ntable]
    int ltab = int(ROUND(double(ntable)/1.25f));
    fltb = double(ltab)/(K/2);
    //double val0 = gsl_sf_bessel_I0(facadj);
    double val0 = bessi0(facadj);
    for (int i=ltab+1; i <= ntable; i++) 
	i0table[i] = 0.f;
    for (int i=0; i <= ltab; i++) 
    {
	double s = double(i)/fltb/N;
	if (s < vadjust) {
	    double rt = sqrt(1.f - pow(s/vadjust, 2));
	    //i0table[i] = gsl_sf_bessel_I0(facadj*rt)/val0;
	    i0table[i] = bessi0(facadj*rt)/val0;
	} else {
	    i0table[i] = 0.f;
	}
    }
}

// Compute the maximum error in the table 
double KaiserBessel::I0table_maxerror() 
{
    double maxdiff = 0.f;
    for (int i = 1; i <= ntable; i++) 
    {
	double diff = fabs(i0table[i] - i0table[i-1]);
	if (diff > maxdiff) 
	    maxdiff = diff;
    }
    return maxdiff;
}

// Kaiser-Bessel Sinh window function
double KaiserBessel::sinhwin(double x) const 
{
    double val0 = sinh(fac)/fac;
    double absx = fabs(x);
    if (0.0 == x) 
    {
	double res = 1.0f;
	return res;
    } 
    else if (absx == alphar) 
    {
	return 1.0f/val0;
    } 
    else if (absx < alphar) 
    {
	double rt = sqrt(1.0f - pow((x/alphar), 2));
	double facrt = fac*rt;
	double res = (sinh(facrt)/facrt)/val0;
	return res;
    } 
    else 
    {
	double rt = sqrt(pow((x/alphar),2) - 1.f);
	double facrt = fac*rt;
	double res = (sin(facrt)/facrt)/val0;
	return res;
    }
}


// Solve second degree equation. ax^2+bx+c=0 -------------------------------
int solve_2nd_degree_eq(float a, float b, float c, float &x1, float &x2,
                        float prec)
{
    // Degenerate case?
    if (ABS(a) < prec)
        if (ABS(b) < prec)
            return -1;
        else
        {
            x1 = -c / b;
            return 1;
        }

    // Normal case
    float d = b * b - 4 * a * c;
    if (d < 0)
        return 0;
    else
    {
        x1 = (-b + sqrt(d)) / (2 * a);
        x2 = (-b - sqrt(d)) / (2 * a);
        return 2;
    }
}

/* Gaussian value ---------------------------------------------------------- */
double gaussian1D(double x, double sigma, double mu)
{
    x -= mu;
    return 1 / sqrt(2*PI*sigma*sigma)*exp(-0.5*((x / sigma)*(x / sigma)));
}

double gaussian2D(double x, double y, double sigmaX, double sigmaY,
                  double ang, double muX, double muY)
{
    // Express x,y in the gaussian internal coordinates
    x -= muX;
    y -= muY;
    double xp = cos(ang) * x + sin(ang) * y;
    double yp = -sin(ang) * x + cos(ang) * y;

    // Now evaluate
    return 1 / sqrt(2*PI*sigmaX*sigmaY)*exp(-0.5*((xp / sigmaX)*(xp / sigmaX) +
                                            (yp / sigmaY)*(yp / sigmaY)));
}

/* Print a boolean value --------------------------------------------------- */
void print(ostream &o, const bool b)
{
    if (b)
        o << "TRUE";
    else
        o << "FALSE";
}

/* Random functions -------------------------------------------------------- */
int idum;

// Uniform distribution ....................................................
void  init_random_generator(int seed)
{
    idum = -1;
    ran1(&idum);
    if (seed != -1)
        for (int i = 0; i < seed; i++)
            ran1(&idum);
}

void  randomize_random_generator()
{
    static  unsigned int seed;
    int rand_return;

    srand(seed);
    rand_return = rand();

    time_t t;
    time(&t);
    rand_return = abs(rand_return);
    idum = (-(int)(t % 10000)
            - (int)(rand_return % 10000));
    ran1(&idum);
    seed = (unsigned int)rand_return;
}

float rnd_unif()
{
    return ran1(&idum);
}
float rnd_unif(float a, float b)
{
    if (a == b)
        return a;
    else
        return a + (b - a)*ran1(&idum);
}

// Gaussian distribution ...................................................
float rnd_gaus()
{
    return gasdev(&idum);
}
float rnd_gaus(float a, float b)
{
    if (b == 0)
        return a;
    else
        return b*gasdev(&idum) + a;
}
float gaus_within_x0(float x0, float mean, float stddev)
{
    float z0 = (x0 - mean) / stddev;
    return erf(ABS(z0) / sqrt(2.0));
}

float gaus_outside_x0(float x0, float mean, float stddev)
{
    float z0 = (x0 - mean) / stddev;
    return erfc(ABS(z0) / sqrt(2.0));
}

float gaus_up_to_x0(float x0, float mean, float stddev)
{
    if (x0 > mean)
        return 1.0 -gaus_outside_x0(x0, mean, stddev) / 2;
    else if (x0 == mean)
        return 0.5;
    else
        return gaus_outside_x0(x0, mean, stddev) / 2;
}

float gaus_from_x0(float x0, float mean, float stddev)
{
    if (x0 > mean)
        return gaus_outside_x0(x0, mean, stddev) / 2;
    else if (x0 == mean)
        return 0.5;
    else
        return 1.0 -gaus_outside_x0(x0, mean, stddev) / 2;
}

float gaus_outside_probb(float p, float mean, float stddev)
{
    // Make a Bolzano search for the right value
    float p1, p2, pm, x1, x2, xm;
    x1 = mean;
    x2 = mean + 5 * stddev;
    do
    {
        xm = (x1 + x2) / 2;
        p1 = gaus_outside_x0(x1, mean, stddev);
        p2 = gaus_outside_x0(x2, mean, stddev);
        pm = gaus_outside_x0(xm, mean, stddev);
        if (pm > p)
            x1 = xm;
        else
            x2 = xm;
    }
    while (ABS(pm - p) / p > 0.005);
    return xm;
}

// See Numerical Recipes, Chap. 6.3
float student_within_t0(float t0, float degrees_of_freedom)
{
    return 1 -betai(degrees_of_freedom / 2, 0.5,
                    degrees_of_freedom / (degrees_of_freedom + t0*t0));
}

float student_outside_t0(float t0, float degrees_of_freedom)
{
    return 1 -student_within_t0(t0, degrees_of_freedom);
}

float student_up_to_t0(float t0, float degrees_of_freedom)
{
    if (t0 >= 0)
        return 1.0 -student_outside_t0(t0, degrees_of_freedom) / 2;
    else
        return student_outside_t0(t0, degrees_of_freedom) / 2;
}

float student_from_t0(float t0, float degrees_of_freedom)
{
    return 1 -student_up_to_t0(t0, degrees_of_freedom);
}

float student_outside_probb(float p, float degrees_of_freedom)
{
    // Make a Bolzano search for the right value
    float p1, p2, pm, t1, t2, tm;
    t1 = 0;
    t2 = 100;
    do
    {
        tm = (t1 + t2) / 2;
        p1 = student_outside_t0(t1, degrees_of_freedom);
        p2 = student_outside_t0(t2, degrees_of_freedom);
        pm = student_outside_t0(tm, degrees_of_freedom);
        if (pm > p)
            t1 = tm;
        else
            t2 = tm;
    }
    while (ABS(pm - p) / p > 0.005);
    return tm;
}

float chi2_up_to_t0(float t0, float degrees_of_freedom)
{
    return gammp(degrees_of_freedom / 2, t0 / 2);
}

float chi2_from_t0(float t0, float degrees_of_freedom)
{
    return 1 -chi2_up_to_t0(t0, degrees_of_freedom);
}

// Log uniform distribution ................................................
float rnd_log(float a, float b)
{
    if (a == b)
        return a;
    else
        return exp(rnd_unif(log(a), log(b)));
}

/* Check if a file exists -------------------------------------------------- */
int exists(const FileName &fn)
{
    FILE *aux;
    if ((aux = fopen(fn.c_str(), "r")) == NULL)
        return 0;
    fclose(aux);
    return 1;
}

/* Wait until file has a stable size --------------------------------------- */
void wait_until_stable_size(const FileName &fn,
                            unsigned long time_step)
{
    if (!exists(fn))
        return;
    struct stat info1, info2;
    if (stat(fn.c_str(), &info1))
        REPORT_ERROR(1,
                     (string)"wait_until_stable_size: Cannot get size of file " + fn);
    off_t size1 = info1.st_size;
    do
    {
        usleep(time_step);
        if (stat(fn.c_str(), &info2))
            REPORT_ERROR(1,
                         (string)"wait_until_stable_size: Cannot get size of file " + fn);
        off_t size2 = info2.st_size;
        if (size1 == size2)
            break;
        size1 = size2;
    }
    while (true);
    return;
}

/* Create empty file ------------------------------------------------------- */
void create_empty_file(const FileName &fn, unsigned long long size,
                       unsigned long long block_size)
{
    unsigned char * buffer = (unsigned char*) calloc(sizeof(unsigned char),
                             block_size);
    if (buffer == NULL)
        REPORT_ERROR(1, "create_empty_file: No memory left");
    FILE * fd = fopen(fn.c_str(), "w");
    if (fd == NULL)
        REPORT_ERROR(1, (string)"create_empty_file: Cannot open file" + fn);
    for (unsigned long i = 0; i < size / block_size; i++)
        fwrite(buffer, sizeof(unsigned char), block_size, fd);
    fwrite(buffer, sizeof(unsigned char), size % block_size, fd);
    fclose(fd);
}

/* Get the Xmipp Base directory -------------------------------------------- */
FileName xmippBaseDir()
{
    string path = getenv("PATH");
    vector<string> directories;
    int number_directories = splitString(path, ":", directories);
    if (number_directories == 0)
        REPORT_ERROR(1, "xmippBaseDir::Cannot find Xmipp Base directory");
    bool found = false;
    int i;
    for (i = 0; i < number_directories; i++)
    {
        FileName fn = directories[i] + "/xmipp_art";
        FILE *aux;
        if ((aux = fopen(fn.c_str(), "r")) != NULL)
        {
            fclose(aux);
            found = true;
            break;
        }
    }
    if (found)
        return directories[i].substr(0, directories[i].length() - 4);
    else
        REPORT_ERROR(1, "xmippBaseDir::Cannot find Xmipp Base directory");
}

// Constructor with root, number and extension .............................
void FileName::compose(const string &str, int no, const string &ext)
{
    *this = (FileName) str;
    if (no != -1)
    {
        char aux_str[5];
        sprintf(aux_str, "%05d", no);
        *this += aux_str;
    }

    if (ext != "")
        *this += "." + ext;
}

// Get the root name of a filename .........................................
FileName FileName::get_root() const
{
    int skip_directories = find_last_of("/") + 1;
    int point = find_first_of(".", skip_directories);
    if (point == -1)
        point = length();
    int root_end = find_last_not_of("0123456789", point - 1);
    if (root_end + 1 != point)
        if (point - root_end > 5)
            root_end = point - 5 - 1;
    return (FileName) substr(0, root_end + 1);
}

// Get the base name of a filename .........................................
string FileName::get_baseName() const
{
    string basename = "";
    string myname = *this;
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
int FileName::get_number() const
{
    int skip_directories = find_last_of("/") + 1;
    int point = find_first_of(".", skip_directories);
    if (point == -1)
        point = length();
    int root_end = find_last_not_of("0123456789", point - 1);
    if (root_end + 1 != point)
    {
        if (point - root_end > 5)
            root_end = point - 5 - 1;
        string aux = substr(root_end + 1, point - root_end + 1);
        return atoi(aux.c_str());
    }
    else
        return -1;
}

// Get the extension of a filename .........................................
string FileName::get_extension() const
{
    int skip_directories = find_last_of("/") + 1;
    int first_point = find_first_of(".", skip_directories);
    if (first_point == -1)
        return "";
    else
        return substr(first_point + 1);
}

// Init random .............................................................
void FileName::init_random(int length)
{
    randomize_random_generator();
    *this = "";
    for (int i = 0; i < length; i++)
        *this += 'a' + FLOOR(rnd_unif(0, 26));
}

// Add at beginning ........................................................
FileName FileName::add_prefix(const string &prefix) const
{
    FileName retval = *this;
    int skip_directories = find_last_of("/") + 1;
    return retval.insert(skip_directories, prefix);
}

// Add at the end ..........................................................
FileName FileName::add_extension(const string &ext) const
{
    if (ext == "")
        return *this;
    else
    {
        FileName retval = *this;
        retval = retval.append((string)"." + ext);
        return retval;
    }
}

// Remove last extension ...................................................
FileName FileName::without_extension() const
{
    FileName retval = *this;
    return retval.substr(0, rfind("."));
}

// Remove root .............................................................
FileName FileName::without_root() const
{
    return without(get_root());
}

// Insert before extension .................................................
FileName FileName::insert_before_extension(const string &str) const
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
FileName FileName::remove_extension(const string &ext) const
{
    int first = find((string)"." + ext);
    if (first == -1)
        return *this;
    else
    {
        FileName retval = *this;
        return retval.erase(first, 1 + ext.length());
    }
}

// Remove all extensions....................................................
FileName FileName::remove_all_extensions() const
{
    int first = rfind("/");
    first = find(".", first + 1);
    if (first == -1)
        return *this;
    else
        return substr(0, first);
}

// Substitute one extension by other .......................................
FileName FileName::substitute_extension(const string &ext1,
                                        const string &ext2) const
{
    int first = find((string)"." + ext1);
    if (first == -1)
        return *this;
    else
    {
        FileName retval = *this;
        return retval.replace(first, 1 + ext1.length(), (string)"." + ext2);
    }
}

// Remove a substring ......................................................
FileName FileName::without(const string &str) const
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
FileName FileName::remove_until_prefix(const string &str) const
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
FileName FileName::remove_directories(int keep) const
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

/* Time managing ----------------------------------------------------------- */
#ifdef _NO_TIME
void time_config()
{}
void annotate_time(TimeStamp *time)
{}
void print_elapsed_time(TimeStamp &time)
{}
float elapsed_time(TimeStamp &time)
{}
float time_to_go(TimeStamp &time, float fraction_done)
{}
void TimeMessage(string message)
{}
void progress_bar(long rlen)
{}
#else
// A global ................................................................
int XmippTICKS;

// Time configuration ......................................................
// The clock frequency for each machine must be known
void time_config()
{
    XmippTICKS = sysconf(_SC_CLK_TCK);
}

// Annotate actual time ....................................................
void annotate_time(TimeStamp *time)
{
    times(time);
}

// Acumulative time
void acum_time(TimeStamp *orig, TimeStamp *dest)
{
    TimeStamp now;
    times(&now);
    (*dest).tms_utime += (*dest).tms_utime + (now.tms_utime - (*orig).tms_utime);
    (*dest).tms_stime += (*dest).tms_stime + (now.tms_utime - (*orig).tms_utime);
}

// Show elapsed time since last annotation .................................
void print_elapsed_time(TimeStamp &time, bool _IN_SECS)
{
    TimeStamp now;
    times(&now);
    float userTime = now.tms_utime - time.tms_utime;
    float sysTime = now.tms_stime - time.tms_stime;
    if (_IN_SECS)
    {
        userTime /= XmippTICKS;
        sysTime /= XmippTICKS;
    }
    cout << "Elapsed time: User(" << userTime << ") System(" << sysTime
    << ")\n";
}

// Calculate elapsed time since last annotation .............................
float elapsed_time(TimeStamp &time, bool _IN_SECS)
{
    TimeStamp now;
    times(&now);
    float userTime = now.tms_utime - time.tms_utime;
    float sysTime = now.tms_stime - time.tms_stime;
    if (_IN_SECS)
    {
        userTime /= XmippTICKS;
        sysTime /= XmippTICKS;
    }
    return userTime + sysTime;
}

// Compute the predicted time left .........................................
float time_to_go(TimeStamp &time, float fraction_done)
{
    TimeStamp now;
    times(&now);
    float totalTime = (now.tms_utime - time.tms_utime +
                       now.tms_stime - time.tms_stime) / XmippTICKS;
    return totalTime*(1 - fraction_done) / fraction_done;
}

// Show a message with the time it is produced .............................
void TimeMessage(string message)
{
    struct tm *T;
    time_t     seconds;

    if (time(&seconds) < 0)
        seconds = 0;
    T = localtime(&seconds);

    printf("%2d:%2d:%2d (day=%2d) =>%s ", T->tm_hour,
           T->tm_min, T->tm_sec, T->tm_mday, message.c_str());
}

// Show a bar with the progress in time ....................................
// When the input is negative then we are setting the progress bar, this
// will be the total of elements to process. Afterwards the call to this
// routine must be in ascending order, ie, 0, 1, 2, ... No. elements
void progress_bar(long rlen)
{
    static time_t startt, prevt;
    time_t currt;
    static long totlen;
    long t1, t2;
    int min, i, hour;
    float h1, h2, m1, m2;

    if (rlen == 0)
        return;
    currt = time(NULL);

    if (rlen < 0)
    {
        totlen = -rlen;
        prevt = startt = currt;
        fprintf(stderr, "0000/???? sec. ");
        for (i = 0; i < 10; i++)
            fprintf(stderr, "------");
        fflush(stderr);
    }
    else if (totlen > 0)
    {
        t1 = currt - startt;               // Elapsed time
        t2 = (long)(t1 * (float)totlen / rlen); // Total time

        hour = 0;
        min = 0;
        if (t2 > 60)
        {
            m1 = (float)t1 / 60.0;
            m2 = (float)t2 / 60.0;
            min = 1;
            if (m2 > 60)
            {
                h1 = (float)m1 / 60.0;
                h2 = (float)m2 / 60.0;
                hour = 1;
                min = 0;
            }
            else
                hour = 0;
        }
        else
            min = 0;

        if (hour)
            fprintf(stderr, "\r%3.2f/%3.2f %s ", h1, h2, "hours");
        else if (min)
            fprintf(stderr, "\r%3.2f/%3.2f %s ", m1, m2, "min");
        else
            fprintf(stderr, "\r%4u/%4u %4s ", (int)t1, (int)t2, "sec.");

        i = (int)(60 * (1 - (float)(totlen - rlen) / totlen));
        while (i--)
            fprintf(stderr, ".");
        if (rlen == totlen)
        {
            fprintf(stderr, "\n");
            totlen = 0;
        }
        fflush(stderr);
        prevt = currt;
    }
}


// Initialize progress bar.

void xmippTextualListener::OnInitOperation(unsigned long _est_it)
{
    progress_bar(-(_est_it));
}

// Show a bar with the progress in time ....................................
// When the input is negative then we are setting the progress bar, this
// will be the total of elements to process. Afterwards the call to this
// routine must be in ascending order, ie, 0, 1, 2, ... No. elements
// Almost identical to previous progress bar function, in fact, we call
// those functions inside.

void xmippTextualListener::OnProgress(unsigned long _it)
{
    progress_bar(_it);
}

// Shows a message indicating the operation in progress.
void xmippTextualListener::OnReportOperation(const string& _rsOp)
{
    fprintf(stderr, _rsOp.c_str());// cout << _rsOp;
}


#endif

/* Little/big endian ------------------------------------------------------- */
// Read in reverse/normal order --------------------------------------------
size_t FREAD(void *dest, size_t size, size_t nitems, FILE * &fp, bool reverse)
{
    size_t retval;
    if (!reverse)
        retval = fread(dest, size, nitems, fp);
    else
    {
        char *ptr = (char *)dest;
        bool end = false;
        retval = 0;
        for (int n = 0; n < nitems; n++)
        {
            for (int i = size - 1; i >= 0; i--)
            {
                if (fread(ptr + i, 1, 1, fp) != 1)
                {
                    end = true;
                    break;
                }
            }
            if (end)
                break;
            else
                retval++;
            ptr += size;
        }
    }
    return retval;
}

// Read in reverse/normal order --------------------------------------------
size_t FWRITE(const void *src, size_t size, size_t nitems, FILE * &fp,
              bool reverse)
{
    size_t retval;
    if (!reverse)
        retval = fwrite(src, size, nitems, fp);
    else
    {
        char *ptr = (char *)src;
        bool end = false;
        retval = 0;
        for (int n = 0; n < nitems; n++)
        {
            for (int i = size - 1; i >= 0; i--)
            {
                if (fwrite(ptr + i, 1, 1, fp) != 1)
                {
                    end = true;
                    break;
                }
            }
            if (end)
                break;
            else
                retval++;
            ptr += size;
        }
    }
    return retval;
}


/** Conversion little-big endian any size */
//#define little22bigendian(x) ByteSwap((unsigned char *) &(x),sizeof(x))
//the prevous define is in the hh file, I leave a commented copy here
//so the code can be more easily understood

void ByteSwap(unsigned char * b, int n)
{
    register int i = 0;
    register int j = n - 1;
    while (i < j)
    {
        std::swap(b[i], b[j]);
        i++, j--;
    }
}

/** Returns 1 if machine is little endian else 0*/

int IsLittleEndian(void)
{
    static const unsigned long ul = 0x00000001;

    return ((int)(*((unsigned char *) &ul)));
}

/** Returns 1 if machine is big endian else 0*/

int IsBigEndian(void)
{
    static const unsigned long ul = 0x01000000;

    return ((int)(*((unsigned char *) &ul)));
}
