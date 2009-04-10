/*******************************************/
/* libSSpectral: -lsrfftw -lsfftw -lIO -lm */
/*                                         */
/* Single-precision spectral analyses and  */
/* correlation functions.                  */
/*                                         */
/* Jun. 16, 2000 Ju Li <liju99@mit.edu>.   */
/*******************************************/

#include "SSpectral.h"

#include "../Spectral.c"

#ifdef _ssmile
int main (int argc, char *argv[])
{
    if (argc!=2)
    {
        printf ("Usage: %s 0 (conductivity), or 1 (viscosity)\n", argv[0]);
        return (1);
    }
    smile (atoi(argv[1]), stdout);
    return (0);
} /* end main() */
#endif  /* _ssmile */


#ifdef _ssmile1
int main (int argc, char *argv[])
{
    if (argc!=2)
    {
        printf ("Usage: %s 0 (conductivity), or 1 (viscosity)\n", argv[0]);
        return (1);
    }
    smile1 (atoi(argv[1]), stdout);
    return (0);
} /* end main() */
#endif  /* _ssmile1 */


#ifdef _spoke1
int main (int argc, char *argv[])
{
    char *fname;
    int correlation_step;
    if (argc==1)
    {
        printf ("Usage: %s <filename> correlation_step\n", argv[0]);
        return (1);
    }
    else if (argc==2)
    {
        fname = "j.out";
        correlation_step = atoi(argv[1]);
    }
    else
    {
        fname = argv[1];
        correlation_step = atoi(argv[2]);
    }
    poke1 (fname, correlation_step, stdout);
    return (0);
} /* end main() */
#endif  /* _spoke1 */
