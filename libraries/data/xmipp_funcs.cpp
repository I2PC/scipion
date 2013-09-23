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

#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <fcntl.h>
#include <string.h>
#ifndef __MINGW32__
 #include <sys/mman.h>
 #ifdef __MACH__
  #include <mach/clock.h>
  #include <mach/mach.h>
 #endif
#endif
#include "xmipp_funcs.h"
#include "numerical_recipes.h"


/* Numerical functions ----------------------------------------------------- */
// Kaiser-Bessel constructor
KaiserBessel::KaiserBessel(double alpha_, int K_, double r_, double v_,
                           int N_, double vtable_, int ntable_)
        : alpha(alpha_), v(v_), r(r_), N(N_), K(K_), vtable(vtable_),
        ntable(ntable_)
{
    // Default values are alpha=1.25, K=6, r=0.5, v = K/2
    if (0.f == v)
        v = double(K)/2;
    if (0.f == vtable)
        vtable = v;
    alphar = alpha*r;
    fac = static_cast<double>(2.*PI)*alphar*v;
    vadjust = 1.0f*v;
    facadj = static_cast<double>(2.*PI)*alphar*vadjust;
    build_I0table();
}

// Kaiser-Bessel I0 selfWindow function
double KaiserBessel::i0win(double x) const
{
    double val0 = double(bessi0(facadj));
    double absx = fabs(x);
    if (absx > vadjust)
        return 0.f;
    double rt = sqrt(1.f - pow(absx/vadjust, 2));
    double res = bessi0(facadj*rt)/val0;
    return res;
}

// Tabulate I0 selfWindow for speed
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
        if (s < vadjust)
        {
            double rt = sqrt(1.f - pow(s/vadjust, 2));
            //i0table[i] = gsl_sf_bessel_I0(facadj*rt)/val0;
            i0table[i] = bessi0(facadj*rt)/val0;
        }
        else
        {
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

// Kaiser-Bessel Sinh selfWindow function
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
int solve_2nd_degree_eq(double a, double b, double c, double &x1, double &x2,
                        double prec)
{
    // Degenerate case?
    if (fabs(a) < prec)
    {
        if (fabs(b) < prec)
            return -1;
        else
        {
            x1 = -c / b;
            return 1;
        }
    }

    // Normal case
    double d = b * b - 4 * a * c;
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

/* t-student value -------------------------------------------------------- */
double tstudent1D(double x, double df, double sigma, double mu)
{
    x -= mu;
    double norm = exp(gammln((df+1.)/2.)) / exp(gammln(df/2.));
    norm /= sqrt(df*PI*sigma*sigma);
    return norm * pow((1 + (x/sigma)*(x/sigma)/df),-((df+1.)/2.));

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

/* ICDF Gaussian ----------------------------------------------------------- */
double icdf_gauss(double p)
{
    const double c[] =
        {
            2.515517, 0.802853, 0.010328
        };
    const double d[] =
        {
            1.432788, 0.189269, 0.001308
        };
    if (p < 0.5)
    {
        // F^-1(p) = - G^-1(p)
        double t=sqrt(-2.0*log(p));
        double z=t - ((c[2]*t + c[1])*t + c[0]) /
                 (((d[2]*t + d[1])*t + d[0])*t + 1.0);
        return -z;
    }
    else
    {
        // F^-1(p) = G^-1(1-p)
        double t=sqrt(-2.0*log(1-p));
        double z=t - ((c[2]*t + c[1])*t + c[0]) /
                 (((d[2]*t + d[1])*t + d[0])*t + 1.0);
        return z;
    }
}

/* CDF Gaussian ------------------------------------------------------------ */
double cdf_gauss(double x)
{
    return 0.5 * (1. + erf(x/sqrt(2.)));
}

/*************************************************************************
Student's t distribution

Computes the integral from minus infinity to t of the Student
t distribution with integer k > 0 degrees of freedom:

                                     t
                                     -
                                    | |
             -                      |         2   -(k+1)/2
            | ( (k+1)/2 )           |  (     x   )
      ----------------------        |  ( 1 + --- )        dx
                    -               |  (      k  )
      sqrt( k pi ) | ( k/2 )        |
                                  | |
                                   -
                                  -inf.

Relation to incomplete beta integral:

       1 - stdtr(k,t) = 0.5 * incbet( k/2, 1/2, z )
where
       z = k/(k + t**2).

For t < -2, this is the method of computation.  For higher t,
a direct method is derived from integration by parts.
Since the function is symmetric about t=0, the area under the
right tail of the density is found by calling the function
with -t instead of t.

ACCURACY:

Tested at random 1 <= k <= 25.  The "domain" refers to t.
                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE     -100,-2      50000       5.9e-15     1.4e-15
   IEEE     -2,100      500000       2.7e-15     4.9e-17

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
double cdf_tstudent(int k, double t)
{
    double EPS=5E-16;
    double result;
    double x;
    double rk;
    double z;
    double f;
    double tz;
    double p;
    double xsqk;
    int j;

    if ( t==0 )
    {
        result = 0.5;
        return result;
    }
    if ( t<-2.0 )
    {
        rk = k;
        z = rk/(rk+t*t);
        result = 0.5*betai(0.5*rk, 0.5, z);
        return result;
    }
    if ( t<0 )
    {
        x = -t;
    }
    else
    {
        x = t;
    }
    rk = k;
    z = 1.0+x*x/rk;
    if ( k%2 != 0 )
    {
        xsqk = x/sqrt(rk);
        p = atan(xsqk);
        if ( k > 1 )
        {
            f = 1.0;
            tz = 1.0;
            j = 3;
            while ( j <= k-2 && tz/f > EPS )
            {
                tz = tz*((j-1)/(z*j));
                f = f+tz;
                j = j+2;
            }
            p = p+f*xsqk/z;
        }
        p = p*2.0/PI;
    }
    else
    {
        f = 1.0;
        tz = 1.0;
        j = 2;
        while ( j<= k-2 && tz/f > EPS)
        {
            tz = tz*((j-1)/(z*j));
            f = f+tz;
            j = j+2;
        }
        p = f*x/sqrt(z*rk);
    }
    if ( t<0 )
    {
        p = -p;
    }
    result = 0.5+0.5*p;
    return result;
}

/* Snedecor's F ------------------------------------------------------------ */
// http://en.wikipedia.org/wiki/F-distribution
double cdf_FSnedecor(int d1, int d2, double x)
{
    return betai(0.5*d1,0.5*d2,(d1*x)/(d1*x+d2));
}

double icdf_FSnedecor(int d1, int d2, double p)
{
    double xl=0, xr=1e6;
    double xm, pm;
    do
    {
        xm=(xl+xr)*0.5;
        pm=cdf_FSnedecor(d1,d2,xm);
        if (pm>p)
        {
            xr=xm;
        }
        else
        {
            xl=xm;
        }
    }
    while (fabs(pm-p)/p>0.001);
    return xm;
}

/* Random functions -------------------------------------------------------- */
int idum;

// Uniform distribution ....................................................
void  init_random_generator(int seed)
{
    idum = -1;
    ran1(&idum);
    if (seed != -1)
    {
        // Prevent seeds larger than 65000
        seed %=0xffff;
        for (int i = 0; i < seed; i++)
            ran1(&idum);
    }
}

unsigned int  randomize_random_generator()
{
    static  unsigned int seed;
    int rand_return;
    struct timespec highresTime;

#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    highresTime.tv_sec = mts.tv_sec;
    highresTime.tv_nsec = mts.tv_nsec;
#else
    clock_gettime(CLOCK_REALTIME, &highresTime);
#endif
    srand(rand()+clock()+time(NULL)+highresTime.tv_nsec);
    rand_return = rand();

    time_t t;
    time(&t);
    //rand_return = abs(rand_return);
    idum = (-(int)(t % 10000)
            - (int)(rand_return % 10000));
    ran1(&idum);
    seed = (unsigned int)rand_return;
    return seed;
}

double rnd_unif()
{
    return ran1(&idum);
}
double rnd_unif(double a, double b)
{
    if (a == b)
        return a;
    else
        return a + (b - a)*ran1(&idum);
}

// t-distribution
double rnd_student_t(double nu)
{
    return tdev(nu, &idum);
}
double rnd_student_t(double nu, double a, double b)
{
    if (b == 0)
        return a;
    else
        return b*tdev(nu, &idum) + a;
}

// Gaussian distribution ...................................................
double rnd_gaus()
{
    return gasdev(&idum);
}
double rnd_gaus(double a, double b)
{
    if (b == 0)
        return a;
    else
        return b*gasdev(&idum) + a;
}
double gaus_within_x0(double x0, double mean, double stddev)
{
    double z0 = (x0 - mean) / stddev;
    return erf(ABS(z0) / sqrt(2.0));
}

double gaus_outside_x0(double x0, double mean, double stddev)
{
    double z0 = (x0 - mean) / stddev;
    return erfc(ABS(z0) / sqrt(2.0));
}

double gaus_up_to_x0(double x0, double mean, double stddev)
{
    if (x0 > mean)
        return 1.0 -gaus_outside_x0(x0, mean, stddev) / 2;
    else if (x0 == mean)
        return 0.5;
    else
        return gaus_outside_x0(x0, mean, stddev) / 2;
}

double gaus_from_x0(double x0, double mean, double stddev)
{
    if (x0 > mean)
        return gaus_outside_x0(x0, mean, stddev) / 2;
    else if (x0 == mean)
        return 0.5;
    else
        return 1.0 -gaus_outside_x0(x0, mean, stddev) / 2;
}

double gaus_outside_probb(double p, double mean, double stddev)
{
    // Make a Bolzano search for the right value
    double pm, x1, x2, xm;
    x1 = mean;
    x2 = mean + 5 * stddev;
    do
    {
        xm = (x1 + x2) / 2;
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
double student_within_t0(double t0, double degrees_of_freedom)
{
    return 1 -betai(degrees_of_freedom / 2, 0.5,
                    degrees_of_freedom / (degrees_of_freedom + t0*t0));
}

double student_outside_t0(double t0, double degrees_of_freedom)
{
    return 1 -student_within_t0(t0, degrees_of_freedom);
}

double student_up_to_t0(double t0, double degrees_of_freedom)
{
    if (t0 >= 0)
        return 1.0 -student_outside_t0(t0, degrees_of_freedom) / 2;
    else
        return student_outside_t0(t0, degrees_of_freedom) / 2;
}

double student_from_t0(double t0, double degrees_of_freedom)
{
    return 1 -student_up_to_t0(t0, degrees_of_freedom);
}

double student_outside_probb(double p, double degrees_of_freedom)
{
    // Make a Bolzano search for the right value
    double pm, t1, t2, tm;
    t1 = 0;
    t2 = 100;
    do
    {
        tm = (t1 + t2) / 2;
        pm = student_outside_t0(tm, degrees_of_freedom);
        if (pm > p)
            t1 = tm;
        else
            t2 = tm;
    }
    while (fabs(pm - p) / p > 0.005);
    return tm;
}

double chi2_up_to_t0(double t0, double degrees_of_freedom)
{
    return gammp(degrees_of_freedom / 2, t0 / 2);
}

double chi2_from_t0(double t0, double degrees_of_freedom)
{
    return 1 -chi2_up_to_t0(t0, degrees_of_freedom);
}

// Log uniform distribution ................................................
double rnd_log(double a, double b)
{
    if (a == b)
        return a;
    else
        return exp(rnd_unif(log(a), log(b)));
}

/* Time managing ----------------------------------------------------------- */
#ifdef _NO_TIME
void time_config()
{}
void annotate_time(ProcessorTimeStamp *time)
{}
void print_elapsed_time(ProcessorTimeStamp &time)
{}
double elapsed_time(ProcessorTimeStamp &time)
{}
double time_to_go(ProcessorTimeStamp &time, double fraction_done)
{}
void TimeMessage(const std::string &message)
{}
void progress_bar(long rlen)
{}
#else
#if defined __MINGW32__ || defined __APPLE__
struct tm* localtime_r (const time_t *clock, struct tm *result)
{
    if (!clock || !result)
        return NULL;
    memcpy(result,localtime(clock),sizeof(*result));
    return result;
}
void sincos(double angle, double * sine, double * cosine)
{
    *sine = sin(angle);
    *cosine = cos(angle);
}
#endif

// A global ................................................................
int XmippTICKS;

// Time configuration ......................................................
// The clock frequency for each machine must be known
void time_config()
{
#ifndef __MINGW32__
    XmippTICKS = sysconf(_SC_CLK_TCK);
#else

    XmippTICKS = CLK_TCK;
#endif
}

#if !defined _NO_TIME && !defined __MINGW32__
// Annotate actual time ....................................................
void annotate_processor_time(ProcessorTimeStamp *time)
{
    times(time);
}
#endif

void annotate_time(TimeStamp *time)
{
    struct timeval  tv;
    gettimeofday(&tv, NULL);
    struct tm tm;
    localtime_r(&tv.tv_sec,&tm);
    *time = tm.tm_hour * 3600 * 1000 + tm.tm_min * 60 * 1000 + tm.tm_sec * 1000 +
            tv.tv_usec / 1000;
}

#if !defined _NO_TIME && !defined __MINGW32__
// Acumulative time
void acum_time(ProcessorTimeStamp *orig, ProcessorTimeStamp *dest)
{
    ProcessorTimeStamp now;
    times(&now);
    (*dest).tms_utime += (*dest).tms_utime + (now.tms_utime - (*orig).tms_utime);
    (*dest).tms_stime += (*dest).tms_stime + (now.tms_utime - (*orig).tms_utime);
}
#endif

#if !defined _NO_TIME && !defined __MINGW32__
// Show elapsed time since last annotation .................................
void print_elapsed_time(ProcessorTimeStamp &time, bool _IN_SECS)
{
    ProcessorTimeStamp now;
    times(&now);
    double userTime = now.tms_utime - time.tms_utime;
    double sysTime = now.tms_stime - time.tms_stime;
    if (_IN_SECS)
    {
        userTime /= XmippTICKS;
        sysTime /= XmippTICKS;
    }
    std::cout << "Elapsed time: User(" << userTime << ") System(" << sysTime
    << ")\n";
}
#endif

void print_elapsed_time(TimeStamp& time, bool _IN_SECS)
{
    struct timeval  tv;
    gettimeofday(&tv, NULL);
    struct tm tm;
    localtime_r(&tv.tv_sec,&tm);
    TimeStamp now = tm.tm_hour * 3600 * 1000 + tm.tm_min * 60 * 1000 + tm.tm_sec * 1000 +
                    tv.tv_usec / 1000;
    TimeStamp diff=now-time;

    std::cout << "Elapsed time: ";
    if (_IN_SECS)
        std::cout << diff/1000.0 << " secs." << std::endl;
    else
        std::cout << diff << " msecs." << std::endl;
}

size_t Timer::now()
{
  gettimeofday(&tv, NULL);
  localtime_r(&tv.tv_sec,&tm);
  return tm.tm_hour * 3600 * 1000 + tm.tm_min * 60 * 1000 + tm.tm_sec * 1000 +
                  tv.tv_usec / 1000;
}

void Timer::tic()
{
  tic_time = now();
}

void Timer::toc(const char * msg)
{
  size_t diff = now() - tic_time;
  if (msg != NULL)
      std::cout << msg;
  std::cout << "Elapsed time: ";
  std::cout << diff/1000.0 << " secs." << std::endl;
}


// Calculate elapsed time since last annotation .............................
double elapsed_time(ProcessorTimeStamp &time, bool _IN_SECS)
{
#if !defined _NO_TIME && !defined __MINGW32__
    ProcessorTimeStamp now;
    times(&now);
    double userTime = now.tms_utime - time.tms_utime;
    double sysTime = now.tms_stime - time.tms_stime;
    if (_IN_SECS)
    {
        userTime /= XmippTICKS;
        sysTime /= XmippTICKS;
    }
    return userTime + sysTime;
#endif
}

#if !defined _NO_TIME && !defined __MINGW32__
// Compute the predicted time left .........................................
double time_to_go(ProcessorTimeStamp &time, double fraction_done)
{
    ProcessorTimeStamp now;
    times(&now);
    double totalTime = (now.tms_utime - time.tms_utime +
                        now.tms_stime - time.tms_stime) / XmippTICKS;
    return totalTime*(1 - fraction_done) / fraction_done;
}
#endif

// Show a message with the time it is produced .............................
void TimeMessage(const std::string & message)
{
    struct tm *T;
    time_t     seconds;

    if (time(&seconds) < 0)
        seconds = 0;
    T = localtime(&seconds);

    printf("%2d:%2d:%2d (day=%2d) =>%s ", T->tm_hour,
           T->tm_min, T->tm_sec, T->tm_mday, message.c_str());
}

// Init progress bar
void init_progress_bar(long total)
{
    progress_bar(-(total));
}

// Show a bar with the progress in time ....................................
// When the input is negative then we are setting the progress bar, this
// will be the total of elements to process. Afterwards the call to this
// routine must be in ascending order, ie, 0, 1, 2, ... No. elements
void progress_bar(long rlen)
{
    static time_t startt;
    time_t currt;
    static long totlen;
    long t1, t2;
    int min, i, hour;
    double h1, h2, m1, m2;
    bool queue = (getenv("XMIPP_IN_QUEUE") != NULL);

    if (rlen == 0)
        return;
    currt = time(NULL);

    if (rlen < 0)
    {
        totlen = -rlen;
        startt = currt;
        fprintf(stderr, "0000/???? sec. ");
        if (!queue)
            for (i = 0; i < 10; i++)
                fprintf(stderr, "------");
        fflush(stderr);
    }
    else if (totlen > 0)
    {
        t1 = currt - startt;               // Elapsed time
        t2 = (long)(t1 * (double)totlen / rlen); // Total time

        hour = 0;
        min = 0;
        if (t2 > 60)
        {
            m1 = (double)t1 / 60.0;
            m2 = (double)t2 / 60.0;
            min = 1;
            if (m2 > 60)
            {
                h1 = (double)m1 / 60.0;
                h2 = (double)m2 / 60.0;
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

        if (!queue)
        {
            i = (int)(60 * (1 - (double)(totlen - rlen) / totlen));
            while (i--)
                fprintf(stderr, ".");
        }

        if (rlen == totlen)
        {
            fprintf(stderr, "\n");
            totlen = 0;
        }
        fflush(stderr);
    }
}

char * getCurrentTimeString()
{
    time_t rawtime;
    time ( &rawtime );
    char * str = ctime (&rawtime);
    char * pos = strrchr(str, '\n');
    pos[0] = '\0'; //Remove \n end character
    return str;

}

// Initialize progress bar.

void TextualListener::OnInitOperation(unsigned long _est_it)
{
    progress_bar(-(_est_it));
}

// Show a bar with the progress in time ....................................
// When the input is negative then we are setting the progress bar, this
// will be the total of elements to process. Afterwards the call to this
// routine must be in ascending order, ie, 0, 1, 2, ... No. elements
// Almost identical to previous progress bar function, in fact, we call
// those functions inside.

void TextualListener::OnProgress(unsigned long _it)
{
    progress_bar(_it);
}

// Shows a message indicating the operation in progress.
void TextualListener::OnReportOperation(const std::string& _rsOp)
{
    fprintf(stderr, "%s", _rsOp.c_str());// std::cout << _rsOp;
}


#endif

/* Little/big endian ------------------------------------------------------- */
// Read in reverse/normal order --------------------------------------------
size_t xmippFREAD(void *dest, size_t size, size_t nitems, FILE * &fp, bool reverse)
{
    size_t retval;
    if (!reverse)
        retval = fread(dest, size, nitems, fp);
    else
    {
        char *ptr = (char *)dest;
        bool end = false;
        retval = 0;
        for (size_t n = 0; n < nitems; n++)
        {
            char * ptrp = ptr + size - 1;
            for (size_t i = 0; i < size; ++i, --ptrp)
            {
                if (fread(ptrp, 1, 1, fp) != 1)
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
    if (retval != nitems)
        REPORT_ERROR(ERR_IO_NOREAD,"XmippFREAD: An error occurred or End of File was reached.");

    return retval;
}

// Read in reverse/normal order --------------------------------------------
size_t xmippFWRITE(const void *src, size_t size, size_t nitems, FILE * &fp,
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
        for (size_t n = 0; n < nitems; n++)
        {
        	char * ptrp = ptr + size - 1;
            for (size_t i = 0; i < size; ++i, --ptrp)
            {
                if (fwrite(ptrp, 1, 1, fp) != 1)
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

/* Map file */
void mapFile(const FileName &filename, char*&map, size_t &size, int &fileDescriptor, bool readOnly)
{
    if (size<0)
    {
        struct stat file_status;
        if(stat(filename.c_str(), &file_status) != 0)
            REPORT_ERROR(ERR_IO_NOPATH,"Cannot get filesize for file "+filename);
        size = file_status.st_size;
    }
#ifdef XMIPP_MMAP
    struct stat file_status;
    if(stat(filename.c_str(), &file_status) != 0)
        REPORT_ERROR(ERR_IO_NOPATH,(String)"Cannot get filesize for file "+filename);
    size = file_status.st_size;
    if(size==0)
        REPORT_ERROR(ERR_IO_NOPATH,(String)"File size=0, cannot read it ("+filename+")");

    if (readOnly)
        fileDescriptor = open(filename.c_str(),  O_RDONLY, S_IREAD);
    else
        fileDescriptor = open(filename.c_str(),  O_RDWR, S_IREAD | S_IWRITE);
    if (fileDescriptor == -1)
        REPORT_ERROR(ERR_IO_NOPATH,(String)"Cannot open file named "+filename);

    if (readOnly)
        map = (char *) mmap(0, size, PROT_READ, MAP_SHARED, fileDescriptor, 0);
    else
        map = (char *) mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fileDescriptor, 0);
    if (map == MAP_FAILED)
        REPORT_ERROR(ERR_MEM_BADREQUEST,"Write can not map memory ");
#else

    map = new char[size];
    fileDescriptor = open(filename.data(), O_RDONLY);
    if (fileDescriptor == -1)
        REPORT_ERROR(ERR_IO_NOPATH,(String)"Cannot open file named "+filename);
    int ok=read(fileDescriptor,map,size);
    if (ok==-1)
        REPORT_ERROR(ERR_IO_NOREAD,(String)"Cannot read from file named"+filename);
#endif
}

/* Unmap file*/
void unmapFile(char *&map, size_t &size, int& fileDescriptor)
{
#ifdef XMIPP_MMAP
    if (munmap(map, size) == -1)
        REPORT_ERROR(ERR_MEM_NOTDEALLOC,"Cannot unmap memory");
#else

    delete []map;
    map=NULL;
#endif

    close(fileDescriptor);
}

/* Conversion little-big endian any size */
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

// Bsoft function
void swapbytes(char* v, unsigned long n)
{
    char            t, t0, t1, t2, t3;
    switch (n)
    {
    case 4:
        t0 = v[0];
        t1 = v[1];
        v[0]=v[3];
        v[1]=v[2];
        v[2]=t1;
        v[3]=t0;
        break;
    case 2:
        t = v[0];
        v[0] = v[1];
        v[1] = t;
        break;
    case 1:
        break;
    case 8:
        t0 = v[0];
        t1 = v[1];
        t2 = v[2];
        t3 = v[4];
        v[0]=v[7];
        v[1]=v[6];
        v[2]=v[5];
        v[3]=v[4];
        v[4]=t3;
        v[5]=t2;
        v[6]=t1;
        v[7]=t0;
        break;
    default:
        for (size_t i=0; i<n/2; i++ )
        {
            t = v[i];
            v[i] = v[n-1-i];
            v[n-1-i] = t;
        }
    }
}

/** Returns true if machine is little endian else false */
bool IsLittleEndian(void)
{
    static const unsigned long ul = 0x00000001;
    return ((int)(*((unsigned char *) &ul)))!=0;
}

/** Returns true if machine is big endian else false*/
bool IsBigEndian(void)
{
    static const unsigned long ul = 0x01000000;
    return ((int)(*((unsigned char *) &ul)))!=0;
}

/** Divides a number into most equally groups */
size_t divide_equally(size_t N, size_t size, size_t rank, size_t &first, size_t &last)
{
    size_t jobs_per_worker = N / size;
    size_t jobs_resting = N % size;

    if (rank < jobs_resting)
    {
        first = rank * (jobs_per_worker + 1);
        last = first + jobs_per_worker;
    }
    else
    {
        first = rank * jobs_per_worker + jobs_resting;
        last = first + jobs_per_worker - 1;
    }

    return last - first + 1;
}
/** In which group from divide_equally is myself? */
size_t divide_equally_group(size_t N, size_t size, size_t myself)
{
    size_t first, last;
    for (size_t rank = 0; rank < size; rank++)
    {
        divide_equally(N, size, rank, first, last);
        if (myself >= first && myself <= last)
            return rank;
    }
    return -1;

}

/**Compare two files **/
bool compareTwoFiles(const FileName &fn1, const FileName &fn2, size_t offset)
{
    char *map1,*map2;
    int fd1, fd2;
    size_t size1=-1, size2=-1;
    mapFile(fn1,map1,size1,fd1,true);
    mapFile(fn2,map2,size2,fd2,true);
    int result=memcmp(map1+offset,map2+offset,size1-offset);
    unmapFile(map1,size1,fd1);
    unmapFile(map2,size2,fd2);
    return (result==0);
}
