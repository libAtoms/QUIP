/***************************************/
/* libIO:                              */
/*                                     */
/* Strings, parsers & files            */
/*                                     */
/* Dec.12, 1999 Ju Li <liju99@mit.edu> */
/***************************************/

#include "IO.h"

/* Miscellaneous IO routines */


/* check the consistency of our library types */
#ifdef _check_TEST
int main (int argc, char *argv[])
{/** common types **/
    
    if (sizeof(int2) != 2) pe ("int2 type undefined\n");
    else printf ("int2 type correct ... \n");
    
    if (sizeof(int4) != 4) pe ("int4 type undefined\n");
    else printf ("int4 type correct ... \n");
    
    if (sizeof(int8) != 8) pe ("int8 type undefined\n");
    else printf ("int8 type correct ... \n");
    
    if (sizeof(real4) != 4) pe ("real4 type undefined\n");
    else printf ("real4 type correct ... \n");
    
    if (sizeof(real8) != 8) pe ("real8 type undefined\n");
    else printf ("real8 type correct ... \n");
    
    /** specific features **/
#ifdef _IRIX64
    if (sizeof(real16) != 16) pe ("real16 type undefined\n");
    else printf ("real16 type correct ... \n");
#endif
    return (1);
}
#endif /* _check_TEST */


/* press return or ctrl-D to continue */
void press_return_to (FILE *in, FILE *out, char *do_what)
{
    int c;
    fprintf (out, "press return to %s ... ", do_what);
    fflush (out);
    while (((c=getc(in))!=EOF)&&(c!='\n'));
    if (c==EOF) fcr(out);
    return;
} /* end press_return_to() */

#ifdef _Press_return_TEST
void main()
{
    int i;
    printf ("Just a test:\n");
    Press_return();
    printf ("well done.\n");
    scanf("%d", &i); getchar();
    printf ("%d\n", i);
}
#endif /* end _Press_return_TEST */


/* Same as fprintf() except it does not print if "stream" is NULL */
void Fprintf (FILE *stream, char *format, ...)
{
    va_list ap;
    if (stream == NULL) return;
    va_start(ap, format);
    vfprintf(stream, format, ap);
    va_end(ap);
    return;
} /* end Fprintf() */


/* same as fgets() except newline is treated as EOF (not stored) */
/* and what is left in that line is read without storing.        */
char *Fgets (char *s, int size, FILE *stream)
{
    char *p, c;
    int fd = fileno (stream);
    for (p=s; p<s+size; p++)
        if ( (read(fd,p,1) == 0) ||
             (*p == '\n') ||
             (*p == 0x0d) )
        {
            if (*p == 0x0d) read(fd,&c,1);
            *p = EOS;
            return (s);
        }
    p--;
    if (p >= s) *p = EOS;
    while ((read(fd,&c,1) != 0) && (c != '\n') && (c != 0x0d));
    if (c == 0x0d) read(fd,&c,1);
    return (s);
} /* end Fgets() */


/* Same as fgets() except we replace ^M in the stream by '\n' */
char *FGETS(char *s, int size, FILE *stream)
{
    register int i,j;
    for (i=0; i<size-1; i++)
    {
        j = fgetc(stream);
        if (j==EOF)
        {
            if (i != 0)
            {
                s[i] = EOS;
                return(s);
            }
            else return(NULL);
        }
        s[i] = j;
        if (s[i] == 0x0d) s[i] = '\n';
        if (s[i] == '\n')
        {
            s[i+1] = EOS;
            return (s);
        }
    }
    s[i] = EOS;
    return (s);
} /* end FGETS() */


/********************************************************/
/* Scan a line from "in", and see if "keys" strings are */
/* at the front. If the line is longer than "linesize", */
/* claim that the file might be corrupted and exit. If  */
/* the line is a simple '\n', return -1. If the line is */
/* a simple EOF, return -4. If multiple "keys" match    */
/* the line, use the first match (therefore if you want */
/* to match "book" AND anything else starting with "bo" */
/* but not "book", you put "book" ahead of "bo" in      */
/* "keys"). Ending EOF or '\n' will always be stripped  */
/* from "linebuffer". If there is a match, the "key"    */
/* string will also be stripped. In the event of '\n'   */
/* termination, return keyword index; in the event of   */
/* EOF termination, return -(index+5). If the line has  */
/* content but does not match, return -2 if '\n'-       */
/* terminated, or -3 if EOF-terminated.                 */
/* -------------------- Summary ----------------------- */
/*  >=0: match by '\n'-terminated lines;                */
/*  -1:  empty '\n' line                                */
/*  -2: '\n'-terminated with content but does not match */
/*  -3:  EOF-terminated with content but does not match */
/*  -4:  empty EOF line                                 */
/*  <=-5: match by EOF-terminated lines                 */
/* "linebuffer" will be stripped of EOF or '\n', and    */
/* if there is a match, the matching string itself.     */
/********************************************************/
int freadline_matchheader
(FILE *in, int n_key, char *keys[], int linesize, char *linebuffer)
{
    register int i,j,k;
    if (!FGETS(linebuffer,linesize,in))
    {
        linebuffer[0] = EOS;
        return(-4);
    }
    if (linebuffer[0] == '\n')
    {
        linebuffer[0] = EOS;
        return(-1);
    }
    for (i=0; linebuffer[i]!=EOS; i++);
    if (i == linesize-1)
        pe ("freadline_matchheader: file has line\n\"%s...\"\n"
            "which is longer than linesize = %d characters;\n"
            "file may be corrupted.\n", linebuffer, linesize);
    for (i=0; i<n_key; i++)
    {
        for (j=0; (linebuffer[j]!=EOS)&&(keys[i][j]==linebuffer[j]); j++);
        if (keys[i][j] == EOS)
        { /* we have a match */
            for (k=0; (linebuffer[j]!='\n') && (linebuffer[j]!=EOS); k++,j++)
                linebuffer[k] = linebuffer[j];
            linebuffer[k] = EOS;
            if (linebuffer[j] == '\n') return(i);
            else return(-(i+5));
        }
    }
    /* does not match */
    for (j=0; linebuffer[j] != EOS; j++);
    if (linebuffer[j-1] == '\n')
    {
        linebuffer[j-1] = EOS;
        return(-2);
    }
    return(-3);
} /* end freadline_matchheader() */


/* return string "1st" for i=1, "2nd" for i=2, "100th" for i=100, etc. */
char *word_for_order (int i)
{
    static char BUFFER [INTEGER_DECIMALCHAR+1];
    switch (i)
    {
    case 1:
        strcpy (BUFFER, "1st");
        break;
    case 2:
        strcpy (BUFFER, "2nd");
        break;
    case 3:
        strcpy (BUFFER, "3rd");
        break;
    default:
        snprintf (BUFFER, INTEGER_DECIMALCHAR+1, "%dth", i);
    }
    return (BUFFER);
} /* end word_for_order() */


/* Convert a string to all uppercase and return back its pointer */
char *convert_to_uppercase (char *s)
{
    register char *p;
    for (p=s; *p != EOS; p++)
        *p = toupper(*p);
    return (s);
} /* end convert_to_uppercase() */


/* Convert a string to all lowercase and return back its pointer */
char *convert_to_lowercase (char *s)
{
    register char *p;
    for (p=s; *p!=EOS; p++)
        *p = tolower(*p);
    return (s);
} /* end convert_to_lowercase() */


/* Similar as UNIX tr: translate characters in s[] from set1 to set2 */
char *convert_tr (char *s, char *set1, char *set2)
{
    register int i, j;
    int len = strlen(set1);
    if (len != strlen(set2))
        pe ("set1 \"%s\" is not equal in length as set2 \"%s\"\n",
            set1, set2);
    for (i=strlen(s); i--;)
        for (j=len; j--;)
            if (s[i] == set1[j]) s[i] = set2[j];
    return (s);
} /* end convert_tr() */
    

/* Return nature of the file */
int Fstat (char filename[])
{
    struct stat buf;
    if ( lstat(filename, &buf) < 0 ) return( FSTAT_LSTAT_RUN_ERROR );
    else if ( S_ISREG(buf.st_mode) ) return( FSTAT_IS_REGULAR_FILE );
    else if ( S_ISDIR(buf.st_mode) ) return( FSTAT_IS_DIRECTORY );
    else if ( S_ISCHR(buf.st_mode) ) return( FSTAT_IS_CHARACTER_SPACIAL );
    else if ( S_ISBLK(buf.st_mode) ) return( FSTAT_IS_BLOCK_SPACIAL );
    else if ( S_ISFIFO(buf.st_mode) ) return( FSTAT_IS_FIFO );
    else if ( S_ISLNK(buf.st_mode) ) return( FSTAT_IS_SYMBOLIC_LINK );
    else if ( S_ISSOCK(buf.st_mode) ) return( FSTAT_IS_SOCKET );
    return ( FSTAT_IS_UNKNOWN );
} /* end Fstat() */


#ifdef _Fexists_TEST
int main (int argc, char *argv[])
{
    printf ("%d\n",Fexists("/etc/passwd"));
    printf ("%d\n",Fexists("/etc/nosuchfile"));
    return (0);
}
#endif /* _Fexists_TEST */


/* Return current umask */
mode_t Umask()
{
    mode_t result;
    result = umask((mode_t)0);
    umask (result);
    return(result);
} /* end Umask() */

#ifdef _Umask_TEST
int main (int argc, char *argv[])
{
    printf ("%03o\n", Umask());
    return (0);
}
#endif /* _Umask_TEST */


/* Create a directory if necessary. Return success or error code */
int dopen (char dirname[])
{
    int status;
    char *p, *buf;
    if (Finvalid(dirname)) return(DOPEN_DIR_NAME_INVALID);
    if (Fexists(dirname))
    {
        if (Fstat(dirname) != FSTAT_IS_DIRECTORY)
            return( DOPEN_NAME_EXISTS_BUT_IS_NOT_A_DIR );
    }
    else
    {
        buf = IOClone(dirname);
      restart:
        p = strrchr(buf, DIRECTORY_SEPARATOR);
        if (p == NULL)
        {
            free (buf);
            if ( mkdir (dirname, 0777 & ~Umask()) == -1)
                return( DOPEN_CANNOT_CREATE_DIR );
        }
        else if (*(p+1) == EOS)
        {
            *p = EOS;
            goto restart;
        }
        else
        {
            *p = EOS;
            status = dopen(buf);
            free(buf);
            if ( status != DOPEN_DIR_EXISTS ) return(status);
            if ( mkdir (dirname, 0777 & ~Umask()) == -1)
                return( DOPEN_CANNOT_CREATE_DIR );
        }
    }
    return(DOPEN_DIR_EXISTS);
} /* end dopen() */


/* Create a directory if necessary; print error message upon failure */
void Dopen (char dirname[])
{
    int status = dopen(dirname);
    /* if (status == DOPEN_NAME_EXISTS_BUT_IS_NOT_A_DIR) */
    /* pe("\"%s\" exists but is not a directory.\n", dirname); */
    if (status == DOPEN_CANNOT_CREATE_DIR)
        pe("Cannot create directory \"%s\".\n", dirname);
    return;
} /* end Dopen() */

#ifdef _Dopen_TEST
int main (int argc, char *argv[])
{
    char *dirname = "/etc///tmp/a//b///c////";
    Dopen (dirname);
    return (0);
}
#endif /* _Dopen_TEST */


/* obtain file handle for specified operations, else print error */
FILE *Fopen (char filename[], char operations[])
{
    char *p, *buf;
    FILE *file;
    if ( (strchr(filename, DIRECTORY_SEPARATOR) != NULL) &&
         (strchr(operations, 'w') != NULL) )
    {
        buf = IOClone(filename);
        p = strrchr(buf, DIRECTORY_SEPARATOR);
        *(p+1) = EOS;
        Dopen (buf);
        free(buf);
    }
    file = fopen (filename, operations);
    if (file == NULL)
        pe ("Fopen: cannot open file \"%s\"\n"
            "for \"%s\" operations.\n", filename, operations);
    return (file);
} /* end Fopen() */


/* Same as strstr() except there is no distinction in case */
char *strcasestr (const char *haystack, const char *needle)
{
    register char c, sc;
    register size_t len;
    if ((c = *needle++) != 0)
    {
        c = LETTER_LOWERCASE(c);
        len = strlen(needle);
        do {
            do {
                if ((sc = *haystack++) == 0)
                    return (NULL);
                sc = LETTER_LOWERCASE(sc);
            } while (sc != c);
        } while (strncasecmp(haystack, needle, len) != 0);
        haystack--;
    }
    return ((char *)haystack);
} /* end strcasestr() */


#ifdef _strcasestr_TEST
int main (int argc, char *argv[])
{
    const char *haystack="VmRSS:      2496 kB";
    const char *needle="rss";
    if ( strcasestr(haystack, needle) )
        printf ("Found \"%s\" in \"%s\".\n", needle, haystack);
    else printf ("Did NOT find \"%s\" in \"%s\".\n", needle, haystack);
    return (0);
}
#endif /* _strcasestr_TEST */


/* If "haystack" contains "needle" as its last word, return */
/* pointer at that last occurrence; otherwise return NULL.  */
char *str_end_with (char *haystack, char *needle)
{
    int i,j;
    i = strlen(haystack);
    j = strlen(needle);
    if (i<j) return(NULL);
    else if (!strcmp(haystack+i-j,needle)) return(haystack+i-j);
    return(NULL);
} /* end str_end_with() */


/* Same as str_end_with() except there is no distinction in case */
char *str_caseend_with (char *haystack, char *needle)
{
    int i,j;
    i = strlen(haystack);
    j = strlen(needle);
    if (i<j) return(NULL);
    else if (!strcasecmp(haystack+i-j,needle)) return(haystack+i-j);
    return(NULL);
} /* end str_caseend_with() */


/* COMPATIBILITY: strsep - extract token from string          */
/* The  strsep()  function  returns  the  next token from the */
/* string stringp which is delimited by delim.  The token  is */
/* terminated with a `\0' character and stringp is updated to */
/* point past the token.   RETURN VALUE:                      */
/* The strsep() function returns a pointer to the  token,  or */
/* NULL if delim is not found in stringp.                     */
char *STRSEP (char **stringp, char *delim)
{
    register char *s;
    register const char *spanp;
    register int c, sc;
    char *tok;
    if ((s = *stringp) == NULL)
        return (NULL);
    for (tok = s;;) {
        c = *s++;
        spanp = delim;
        do
        {
            if ((sc = *spanp++) == c)
            {
                if (c == 0)
                    s = NULL;
                else
                    s[-1] = 0;
                *stringp = s;
                return (tok);
            }
        } while (sc != 0);
    }
}  /* end STRSEP() */


/* Get the absolute pathname of a shell path */
char *absolute_pathname (char *pathname)
{
    static char pname[PATH_MAX+1];
    if (*pathname == '/') return (pathname);
    if (*pathname == '~')
        if (*(pathname+1) == '/')
        {
            if (snprintf (pname, PATH_MAX+1, "%s/%s", getenv("HOME"),
                          pathname+2) == -1)
                pe ("absolute_pathname: pathname too long.\n");
            return (pname);
        }
    if (snprintf (pname, PATH_MAX+1, "%s/%s", getenv("PWD"),
                  pathname) == -1)
        pe ("absolute_pathname: pathname too long.\n");
    return (pname);
} /* end absolute_pathname() */

#ifdef _absolute_pathname_TEST
int main (int argc, char *argv[])
{
    char *testpath[] = { "/tmp/a.x", "a.x", "~/a.x" };
    int i;
    for (i=0; i<sizeof(testpath)/sizeof(char *); i++)
        printf ("\"%s\" -> \"%s\" \n",
                testpath[i], absolute_pathname(testpath[i]) );
    return(0);
}
#endif /* _absolute_pathname_TEST */


/* Return the basename of a file (without suffix): "/tmp/a.b.c"->"a.b"  */
/* same as echo /tmp/a.b.c | sed -e "s/.*\///g" | sed -e 's/\.[^.]*$//' */
char *file_basename1 (char *pathname)
{
    static char pname[PATH_MAX+1];
    char *q;

    for (q=pathname+strlen(pathname); q>pathname; q--)
        if (*(q-1)=='/') break;
    strcpy(pname,q);
    for (q=pname+strlen(pname)-1; q>=pname; q--)
        if (*q=='.')
        {
            *q=EOS;
            break;
        }
    return (pname);
} /* end file_basename1() */

#ifdef _file_basename1_TEST
int main (int argc, char *argv[])
{
    char *testpath[] = { "/tmp/c.x/a.x.exe", "a.x", "~/a.x" };
    int i;
    for (i=0; i<sizeof(testpath)/sizeof(char *); i++)
        printf ("\"%s\" -> \"%s\" \n",
                testpath[i], file_basename1(testpath[i]) );
    return(0);
}
#endif /* _file_basename1_TEST */


/* Check whether an executable "command" exists in current PATH */
int command_exists (char *command)
{
    char *path, *p, *q, filename[PATH_MAX+1];
    int commandlen, plen;
    if (strchr(command, '/')) return(Fexecutable(command));
    if ((path=getenv("PATH")) == NULL)
        pe("command_exists: can't get $PATH from environment.\n");
    if ((path = strdup(path)) == NULL)
        pe("command_exists: can't allocate memory.");
    commandlen = strlen(command);
    for (q=path; (p=STRSEP(&q,":"))!=NULL;)
    {
        if (*p == EOS) p = ".";
        plen = strlen(p);
        while (p[plen-1] == '/') p[--plen] = EOS;
        if (plen+1+commandlen > PATH_MAX)
            pe("command_exists:\n%s/%s:\n%s", p, command,
               strerror(ENAMETOOLONG));
        strcpy(filename, p);
        filename[plen] = '/';
        strcpy(filename+plen+1, command);
        if (Fexecutable(filename))
        {
            free(path);
            return(TRUE);
        }
    }
    free(path);
    return(FALSE);
} /* end command_exists() */

#ifdef _command_exists_TEST
int main (int argc, char *argv[])
{
    printf ("%d\n", command_exists("em"));
    printf ("%d\n", command_exists("nosus"));
    return (0);
}
#endif /* _command_exists_TEST */


int commands_exists (int n_commands, char *commands[])
{
    int i;
    for (i=0; i<n_commands; i++)
        if (command_exists(commands[i])) return(TRUE);
    return(FALSE);
} /* end commands_exists() */


char *postscript_viewers[NUMBER_POSTSCRIPT_VIEWERS] =
{"gv", "ghostview", "xpsview", "imagetool", "pageview", "xv"};
char *raster_viewers[NUMBER_RASTER_VIEWERS] =
{"xv", "imagetool", "imgview", "netscape"};

#define BUFFER_SIZE 1024
int try_to_runfg (int n_commands, char *commands[], char *arg)
{
    int i;
    char *command, buffer[BUFFER_SIZE];
    buffer[0] = EOS;
    for (i=0; i<n_commands; i++)
    {
        command=commands[i];
        STR_append(buffer,command);
        STR_append(buffer,",");
        if ( command_exists(command) )
            return(system(str3(command," ",arg)));
    }
    if ((command=eos(buffer)) > buffer)
        if (*(command-1)==',') *(command-1) = EOS;
    fprintf(stderr, "try_to_runfg: \"%s\" not found in your run PATH.\n",
            buffer);
    return(-1);
} /* end try_to_runfg() */


/* run command in foreground; return -1 if command not found */
int TRY_to_runfg (char *arg, ...)
{
    char *command, buffer[BUFFER_SIZE];
    va_list ap;
    va_start (ap, arg);
    buffer[0] = EOS;
    while ( (command=va_arg(ap,char *)) != NULL )
    {
        STR_append(buffer,command);
        STR_append(buffer,",");
        if ( command_exists(command) )
        {
            va_end(ap);
            return(system(str3(command," ",arg)));
        }
    }
    va_end(ap);
    if ((command=eos(buffer)) > buffer)
        if (*(command-1)==',') *(command-1) = EOS;
    fprintf(stderr, "TRY_to_runfg: \"%s\" not found in your run PATH.\n",
            buffer);
    return(-1);
} /* end TRY_to_runfg() */

#ifdef _TRY_to_runfg_TEST
int main (int argc, char *argv[])
{
    TRY_to_runfg ("/asm/home/Moon/www/Archive/Photos/Baby/Later/DCP_1126.JPG",
                  "nosuchfile", "xv", NULL);
    printf ("so be it!\n");
    return(0);
}
#endif /* _TRY_to_runfg_TEST */


/* run command in background; return -1 if command not found */
int try_to_runbg (int n_commands, char *commands[], char *arg)
{
    int i;
    char *command, buffer[BUFFER_SIZE];
    buffer[0] = EOS;
    for (i=0; i<n_commands; i++)
    {
        command=commands[i];
        STR_append(buffer,command);
        STR_append(buffer,",");
        if ( command_exists(command) )
            return(system(str4(command," ",arg," &")));
    }
    if ((command=eos(buffer)) > buffer)
        if (*(command-1)==',') *(command-1) = EOS;
    fprintf(stderr, "try_to_runbg: \"%s\" not found in your run PATH.\n",
            buffer);
    return(-1);
} /* end try_to_runbg() */


int TRY_to_runbg (char *arg, ...)
{
    char *command, buffer[BUFFER_SIZE];
    va_list ap;
    va_start (ap, arg);
    buffer[0] = EOS;
    while ( (command=va_arg(ap,char *)) != NULL )
    {
        STR_append(buffer,command);
        STR_append(buffer,",");
        if ( command_exists(command) )
        {
            va_end(ap);
            return(system(str4(command," ",arg," &")));
        }
    }
    va_end(ap);
    if ((command=eos(buffer)) > buffer)
        if (*(command-1)==',') *(command-1) = EOS;
    fprintf(stderr, "TRY_to_runbg: \"%s\" not found in your run PATH.\n",
            buffer);
    return(-1);
} /* end TRY_to_runbg() */

#ifdef _TRY_to_runbg_TEST
int main (int argc, char *argv[])
{
    TRY_to_runbg ("/asm/home/Moon/www/Archive/Photos/Baby/Later/DCP_1126.JPG",
                  "nosuchfile", "xv", NULL);
    printf ("so be it!\n");
    return(0);
}
#endif /* _TRY_to_runbg_TEST */
#undef BUFFER_SIZE


/********************************************************************/
/* File format magic numbers:                                       */
/* http://cvs.gnome.org/lxr/source/gnome-libs/gnome-data/mime-magic */
/********************************************************************/

#define BZIP2_MAGIC_CHAR 3
static int BZIP2_MAGIC [BZIP2_MAGIC_CHAR] = {'B','Z','h'};

/*****************************************************************/
/* Check if the file is in BZIP2 format: if it is, decompress it */
/* using environmental 'bzip2' and return a temporary filename   */
/* that points to the decompressed file. Otherwise return NULL.  */
/*****************************************************************/
char *RBZIP2_fname (char *original_fname)
{
    register int i;
    char *fname;
    FILE *file = rOpen (original_fname);
    for (i=0; i<BZIP2_MAGIC_CHAR; i++)
        if (BZIP2_MAGIC[i] != fgetc(file))
        {
            fclose(file);
            return (NULL);
        }
    fclose(file);
    fname = TMPNAM(NULL);
    if (system(str4("bzip2 -cd ",original_fname," > ",fname)) == -1)
        pe("RBZIP2_fname: file \"%s\" is thought to be a .bz2 file but\n"
           "we encounter error when executing \"%s\".\n", original_fname,
           str4("bzip2 -cd ",original_fname," > ",fname));
    return (fname);
} /* end RBZIP2_fname() */


#define GZIP_MAGIC_CHAR 2
static int GZIP_MAGIC [GZIP_MAGIC_CHAR] = {0x1f,0x8b};

/****************************************************************/
/* Check if the file is in GZIP format: if it is, decompress it */
/* using environmental 'gzip' and return a temporary filename   */
/* that points to the decompressed file. Otherwise return NULL. */
/****************************************************************/
char *RGZIP_fname (char *original_fname)
{
    register int i;
    char *fname;
    FILE *file = rOpen (original_fname);
    for (i=0; i<GZIP_MAGIC_CHAR; i++)
        if (GZIP_MAGIC[i] != fgetc(file))
        {
            fclose(file);
            return (NULL);
        }
    fclose(file);
    fname = TMPNAM(NULL);
    if (system(str4("gzip -cd ",original_fname," > ",fname)) == -1)
        pe("RGZIP_fname: file \"%s\" is thought to be a .gz file but\n"
           "we encounter error when running \"%s\".\n", original_fname,
           str4("gzip -cd ",original_fname," > ",fname));
    return (fname);
} /* end RGZIP_fname() */


#define UNIXZ_MAGIC_CHAR 2
int UNIXZ_MAGIC [UNIXZ_MAGIC_CHAR] = {037,0235};

/******************************************************************/
/* Check if the file is in UNIX .z format: if it is, decompress   */
/* using environmental 'compress' and return a temporary filename */
/* that points to the decompressed file. Otherwise return NULL.   */
/******************************************************************/
char *RUNIXZ_fname (char *original_fname)
{
    register int i;
    char *fname;
    FILE *file = rOpen (original_fname);
    for (i=0; i<UNIXZ_MAGIC_CHAR; i++)
        if (UNIXZ_MAGIC[i] != fgetc(file))
        {
            fclose(file);
            return (NULL);
        }
    fclose(file);
    fname = TMPNAM(NULL);
    if (system(str4("compress -cd ",original_fname," > ",fname)) == -1)
        pe("RUNIXZ_fname: file \"%s\" is thought to be a .z file but\n"
           "we encounter error when running \"%s\".\n", original_fname,
           str4("compress -cd ",original_fname," > ",fname));
    return (fname);
} /* end RUNIXZ_fname() */


/***************************************************************/
/* If the file is .bz2, .gz or .z compressed, decompress using */
/* the environment 'bzip2', 'gzip' or 'compress' and return a  */
/* handle to file which will be automatically removed when the */
/* handle is closed. Otherwise return an ordinary read handle. */
/***************************************************************/
FILE *ROpen (char *original_fname)
{
    char *fname;
    int fd;
    if ( ! (fname = RBZIP2_fname(original_fname)) )
        if ( ! (fname = RGZIP_fname(original_fname)) )
            if ( ! (fname = RUNIXZ_fname(original_fname)) )
                return(rOpen(original_fname));
    fd = open(fname, O_RDWR);
    unlink(fname);
    return(fdopen(fd, "rw"));
} /* end ROpen() */


/*************************************************************/
/* Realize the compression of a file whose name ends with    */
/* ".bz2" or ".gz". Return after-processed filesize in bytes */
/* Also, any original_fname that ends with ".blitz" is going */
/* to have ".blitz" taken down in a blitzkrieg fashion.      */
/*************************************************************/
long Zrealize (char *original_fname)
{
    char *p, fname1[1024], fname2[1024];
    if (Fexists(p=BLITZKRIEG_FMASK(original_fname))) strcpy (fname1,p);
    else
    {
        if ( Finvalid(original_fname) ) return(0);
        strcpy (fname1, original_fname);
    }
    strcpy (fname2, fname1);
    if ( str_end_with(fname1, ".bz2") ||
         str_end_with(fname1, ".gz") )
    {
        STR_append (fname1, BLITZKRIEG_SUFFIX);
        rename( fname2, fname1 );
    }
    if ( (p = str_end_with(fname1, ".bz2"BLITZKRIEG_SUFFIX)) )
    {
        if ( system(str2("bzip2 -7 ", fname1)) != 0 )
            pr("Zrealize: sorry, file \"%s\" remains "
               "uncompressed.\n", fname1);
        STR_append (fname1, ".bz2");
    }
    else if ( (p = str_end_with(fname1, ".gz"BLITZKRIEG_SUFFIX)) )
    {
        if ( system(str2("gzip -7 ", fname1)) != 0 )
            pr("Zrealize: sorry, file \"%s\" remains "
               "uncompressed.\n", fname1);
        STR_append (fname1, ".gz");
    }
    if ( (p = strrstr(fname1, BLITZKRIEG_SUFFIX)) )
    {
        strcpy (fname2, fname1);
        *p = EOS;
        rename( fname2, fname1);
        /* 'rename' is sometimes implementation dependent */
        /* system(str4("mv ", fname2, " ", fname1)); */
    }
    return(Fsize(fname1));
} /* end Zrealize() */


#ifdef _Zrealize_TEST
/* #define FNAME "/tmp/try" */
#define FNAME "/tmp/try.bz2"
int main (int argc, char *argv[])
{
    FILE *fp = Wopen (FNAME);
    fprintf (fp, "This is a test.\n");
    fclose (fp);
    printf ("from %ld bytes.\n", Fsize(BLITZKRIEG_FMASK(FNAME)));
    printf ("to %ld bytes.\n", Zrealize(FNAME));
    return (0);
}
#endif /* _Zrealize_TEST */


/* Test if a filename is writable by opening and closing it */
int tested_to_be_writable (char filename[])
{
    char *p, *buf;
    FILE *fp;
    if ( strchr(filename, DIRECTORY_SEPARATOR) != NULL )
    {
        buf = IOClone(filename);
        p = strrchr(buf, DIRECTORY_SEPARATOR);
        *(p+1) = EOS;
        dopen (buf);
        free(buf);
    }
    fp = wopen(filename);
    if (fp)
    {
        fclose (fp);
        return (TRUE);
    }
    return (FALSE);
} /* end tested_to_be_writable() */


/* check if a file exists, and if so ask the user whether to overwrite it */
int Freetowrite (char filename[])
{
    char c;
    if (Fexists(filename)) 
    {
        printf ("File \"%s\" exists, overwrite (y/n)? ", filename);
        /* first char is answer, but keep at it until return is received */
        if ((c=getc(stdin))!='\n') while(getc(stdin)!='\n');
        /* free to write only when answer is affirmative */
        if (!((c=='y')||(c=='Y'))) return (FALSE);
        if (!Fwritable(filename))
        {
            printf ("Sorry, you do not have write permission.\n");
            return (FALSE);
        }
    }
    return (TRUE);
} /* end Freetowrite() */

#ifdef _Freetowrite_TEST
#define testfile "/etc/passwd"
void main()
{
    if (Freetowrite(testfile))
        printf("You've got user's consent, "
               "but go ahead and try to delete me :)\n");
    else printf("Even the user forbids you to do that.\n");
}
#endif


/* obtain size in bytes of a file with filename "fname" */
long Fsize (char *fname)
{
    struct stat buf[1];
    stat (fname, buf);
    return ((long)(buf->st_size));
} /* end Fsize() */


/* obtain size in bytes of an opened file with file number "fid" */
long fsize (int fid)
{
    struct stat buf[1];
    fstat (fid, buf);
    return ((long)(buf->st_size));
} /* end fsize() */

#ifdef _fpsize_TEST
int main (int argc, char *argv[])
{
    FILE *fp = rOpen ("/etc/passwd");
    printf ("%d\n", fpsize(fp));
    return (0);
}
#endif /* _fpsize_TEST */


/* end of string */
char *eos (char *s)
{
    char *p;
    for (p=s; *p!=EOS; p++);
    return (p);
} /* end eos() */


/* copy char to string at cell_numberth static cell and return its pointer */
char *c2s (int cell_number, char c)
{
    static char buffer[C2S_MAX_CELLS][2] = {{0}};
    if ( (cell_number < 0) || (cell_number >= C2S_MAX_CELLS) )
    {
        printf ("error: c2s: convert '%c' to string \n"
                "uses cell number = %d > C2S_MAX_CELLS = %d.\n",
                c, cell_number, C2S_MAX_CELLS);
        exit(1);
    }
    buffer[cell_number][0] = c;
    return (&buffer[cell_number][0]);
} /* end c2s() */


/* forward to char position in "s" that is not a white space, */
/* i.e., form-feed '\f', newline '\n', carriage return '\r',  */
/* horizontal tab '\t' and vertical tab ('\v').               */
char *space_advance(char *s)
{
    register char *p;
    for (p=s; ((*p)!=EOS) && isspace(*p); p++);
    return (p);
} /* end space_advance() */


/* forward to char position in "s" that is not a blank, i.e., space or tab */
char *blank_advance(char *s)
{
    register char *p;
    for (p=s; ((*p)!=EOS) && ISBLANK(*p); p++);
    return (p);
} /* end blank_advance() */


/* forward to char position in "s" that IS a space or tab, or stop at EOS */
char *nonblank_advance(char *s)
{
    register char *p;
    for (p=s; ((*p)!=EOS) && ISNOTBLANK(*p); p++);
    return (p);
} /* end nonblank_advance() */


/* Remove trailing blanks: "Ju Li  " -> "Ju Li" */
char *remove_trailing_blanks (char *p)
{
    register char *q;
    for (q=eos(p); (q>p)&&ISBLANK(*(q-1)); q--);
    *q = EOS;
    return (p);
} /* end remove_trailing_blanks() */


/***************************************************************/
/* Concatenate several strings together: at most STRS_MAX_CHAR */
/* characters in all. strs() is reusable STRS_CARTRIDGE times. */
/***************************************************************/
char *strs (int number_of_strings, ...)
{
    int i, j, cartridge_used[STRS_CARTRIDGE]={0};
    char *source, *dest, *r, buffer[STRS_MAX_CHAR+1];
    static char BUFFER[STRS_CARTRIDGE][STRS_MAX_CHAR+1];
    va_list ap;
    va_start(ap, number_of_strings);
    dest = buffer;
    for (i=0; i<number_of_strings; i++)
    {
        source = va_arg(ap, char *);
        if (source == NULL) continue;
        for (j=0; j<STRS_CARTRIDGE; j++)
            if ( (source >= BUFFER[j]) &&
                 (source <  BUFFER[j]+STRS_MAX_CHAR+1) )
                cartridge_used[j] = 1;
        for (r = dest;
             (r < buffer+STRS_MAX_CHAR) &&
                 (source[r-dest] != EOS); r++)
            *r = source[r-dest];
        dest = r;
        if (r == buffer+STRS_MAX_CHAR) break;
    }
    *dest = EOS;
    va_end(ap);
    for (j=0; j<STRS_CARTRIDGE; j++)
        if (!cartridge_used[j])
            return(strcpy(BUFFER[j],buffer));
    pe ("strs: run out of cartridges (%d)\n", STRS_CARTRIDGE);
    return(NULL);
} /* end strs() */

#ifdef _strs_TEST
void main()
{
    printf("\"%s\"\n",
           strs(3,strs(3,"Mary", " has a ", "little lamb."),
                "\n",
                strs(3,"Mary", " has a ", "little lamb.")));
    return;
}
#endif


/* static memory driver of vsnprintf() */
char *vstrf (char *format, va_list ap)
{
    static char buffer[STRF_MAX_CHAR+1];
    vsnprintf(buffer, STRF_MAX_CHAR+1, format, ap);
    return (buffer);
} /* end strf() */


char *strf (char *format, ...)
{
    char *s;
    va_list ap;
    va_start(ap, format);
    s = vstrf(format,ap);
    va_end(ap);
    return (s);
} /* end strf() */


#ifdef _strf_TEST
void main()
{
    printf("\"%s\"\n", strf("%s%s%s","Mary", " has a ", "little lamb."));
}
#endif


/* change stderr to a different stream */
void redirect_stderr_to (FILE *fp)
{
    fclose (stderr);
#if ( \
  defined(_SunOS)  ||\
  defined(_IRIX64) ||\
  defined(_IRIX)   ||\
  defined(_HPUX)   ||\
  defined(_OSF1)   ||\
  defined(_Darwin) )
    /* stderr is a macro: "&__iob[2]", not a valid left-value */
    *stderr = *fp;
#else
    /* stderr is a global variable */
    stderr = fp;
#endif
    return;
} /* end redirect_stderr_to() */


/* perror() driver using strf() */
void pe (char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    fprintf (stderr, "error: ");
    if (errno==0) vfprintf(stderr,format,ap);
    else perror(vstrf(format,ap));
    va_end(ap);
    exit(1);
    return;
} /* end pe() */

#ifdef _pe_TEST
int main()
{
    errno = E2BIG;
    pe ("test: errno = %d\n", errno);
    return(0);
}
#endif


/* print to stderr stream */
void pr (char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr,format,ap);
    va_end(ap);
    return;
} /* end pr() */

#ifdef _pr_TEST
int main()
{
    errno = E2BIG;
    pr ("test: errno = %d\n", errno);
    return(0);
}
#endif


#ifndef _HPUX
/* Find the last occurrence of "needle" in "haystack"     */
/* and return pointer to that location; else return NULL. */
char *strrstr (char *haystack, char *needle)
{
    char *last = NULL;
    while ( (haystack=strstr(haystack,needle)) )
        last = (haystack++);
    return (last);
} /* end strrstr() */
#endif


#ifdef _strrstr_TEST
int main (int argc, char *argv[])
{
    char *a = "rrrrr";
    char *b = "r";
    printf ("\"%s\"\n", strrstr(a,b));
    return (0);
}
#endif /* _strrstr_TEST */


/***********************************************************/
/* assume c[0] is the first byte of char rendering device  */
/* crd[rows][cols+1]. We snprintf in fmt to &crd[i][0] the */
/* ith row of A for at most cols characters (not including */
/* EOS), i=0..rows-1, then return the largest width (not   */
/* including EOS). crd[i][?..width-1] is padded with ' '.  */
/***********************************************************/
#define crd(i,j) (c[(i)*(cols+1)+(j)])
#define printtoken(TYPE) \
 (snprintf(r, &crd(i,cols)-r+1, inert, *((TYPE *)B)), B += sizeof(TYPE))

int Asprintf (char c[], int cols, int rows, char *fmt, void *A)
{
    int i, j, islong, width;
    char *inert, *start, *end=NULL, *r, cc, f;
    char *linefmt, *B, *s;

    linefmt = IOclone(fmt);
    B = (char *)A;
    width = 0;
    for (i=0; i<rows; i++)
    {
        r = &crd(i,0);
        for ( inert=linefmt; ; inert = end+1 )
        {
            for ( start=inert; *start!=EOS; start++ )
            {
                if ( *start == '%' )
                {
                    if ( *(start+1) == '%' )
                    { 
                        start++;
                        continue;
                    }
                    else
                    {
                        for ( islong=FALSE, end=start+1;
                              ISDIGIT(*end)
                                  || (*end=='.')
                                  || (*end=='#')
                                  || (*end=='+')
                                  || (*end=='-')
                                  || (*end=='l');
                              end++ )
                            if (*end == 'l') islong = TRUE;
                        if ( !ISALPHA(*end) )
                        {
                            printf ("error: Asprintf: \"%s\" "
                                    "contains single %%\n", linefmt);
                            exit(1);
                        }
                        cc = *(end+1);
                        *(end+1) = EOS;
                        switch ( *end )
                        {
                        case 'o':
                        case 'u':
                        case 'x':
                        case 'X':
                            f = *end;
                            *end = EOS;
                            if (end == start+1) j = 0;
                            else j = atoi(start+1);
                            *end = f;
                            if (j == 2)
                            {  /* %02x is a call sign to print bytes */
                                s = B;
                                j = *B;
                                B = (char *)(&j);
                                printtoken (int);
                                B = s + sizeof(char);
                            }
                            else if ( islong ) printtoken(unsigned int);
                            else printtoken(unsigned short);
                            break;
                        case 'd':
                        case 'i':
                            if ( islong ) printtoken(long int);
                            else printtoken(int);
                            break;
                        case 'f':
                            if ( islong ) printtoken(double);
                            else printtoken(float);
                            /* cosmetics */
                            for (j=FALSE; *r!=EOS; r++)
                                if (*r=='.') j=TRUE;
                            if (j)
                                for (r--;
                                     ((*r==' ')||(*r=='0')||(*r=='.'));
                                     r--)
                                {
                                    if (*r=='.')
                                    {
                                        *r = ' ';
                                        break;
                                    }
                                    else *r = ' ';
                                }
                            break;
                        case 'g':
                        case 'G':
                        case 'e':
                        case 'E':
                            if ( islong ) printtoken(double);
                            else printtoken(float);
                            break;
                        case 'c':
                            printtoken(char);
                            break;
                        case 's':
                            printtoken(char *);
                            break;
                        default:
                            printf ("error: Asprintf: \"%s\" descriptor in "
                                    "\"%s\" unsupported\n", start, linefmt);
                            exit(1);
                        }
                        *(end+1) = cc;
                        while (*r != EOS) r++;
                        break;
                    }
                }
            }  /* find escape sequence start and treat it */
            if (*start==EOS) break;  /* run out of commands */
        } /* loop over inert beginnings */
        snprintf(r, &crd(i,cols)-r+1, inert);
        while (*r != EOS) r++;
        if ( r - &crd(i,0) > width )  width = r - &crd(i,0);
    } /* loop over row number i */
    for (i=0; i<rows; i++)
    {
        for (j=0; crd(i,j)!=EOS; j++);
        for (; j<width; j++) crd(i,j) = ' ';
        crd(i,j) = EOS;
    }
    IOfree(linefmt);
    return (width);
} /* end Asprintf() */


/* fprintf to fp rows x fmt from A; return largest row strlen */
int Afprintf (FILE *fp, int cols, int rows, char *fmt, void *A)
{
    int i, j, width;
    char *crd;
    crd = IOmem(rows*(cols+1));
    width = Asprintf (crd, cols, rows, fmt, A);
    for (i=0; i<rows; i++)
    {
        for (j=0; j<width; j++) fputc (crd[i*(cols+1)+j], fp);
        fputc ('\n', fp);
    }
    return (width);
} /* end Afprintf() */


#ifdef _Aprint_TEST
#define ROWS 3
void main()
{
    int i, width;
    double A[ROWS][ROWS] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    int B[ROWS][ROWS] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    char C[ROWS][ROWS] = {{'a', 'b', 'c'}, {'d', 'e', 'f'},
                          {'g', 'h', 'i'}};
    char *D[ROWS][ROWS] = {{"a", "b", "c"}, {"d", "e", "f"},
                           {"", "h", "ii"}};
    width = Aprint (ROWS, "\"%4.2lf %4.2lf %4.2lf %%\"", A);
    printf ("width = %d\n\n", width);
    width = Aprint (ROWS, "\"%ld %ld %ld\"", B);
    printf ("width = %d\n\n", width);
    width = Aprint (ROWS, "\"%c %c %c\"", C);
    printf ("width = %d\n\n", width);
    width = Aprint (ROWS, "\"%s %s %s\"", D);
    printf ("width = %d\n\n", width);
    return;
} 
#endif


/* fork out a child process running slave_work_till_die() */
void spinoff( void (*slave_work_till_die)() )
{
    int pid;
    if ( (pid = fork()) < 0 )
    {
        printf ("error: spinoff: cannot fork\n");
        exit(1);
    }
    else if (pid == 0)
    { /* child process */
        slave_work_till_die();
        exit(0);
    }
    return;
} /* end spinoff() */


/* base64 encoding map */
char Base64Alphabet[64] =
{'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P',
 'Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f',
 'g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v',
 'w','x','y','z','0','1','2','3','4','5','6','7','8','9','+','/'};


/* Economically print floating-pt numbers by "0.123" -> ".123" reduction. */
#define VZFPRINTF_SIZE 1024
void vzfprintf (FILE *stream, char *format, va_list ap)
{
    char buf[VZFPRINTF_SIZE], *p, *q;
    if (stream == NULL) return;
    if (vsnprintf (buf, VZFPRINTF_SIZE, format, ap) == -1)
        pe ("vzfprintf: VZFPRINTF_SIZE = %d exceeded.\n", VZFPRINTF_SIZE);
    va_end(ap);
    q = buf;
    if (*q=='0') if (*(q+1)=='.') q++;
    while ((p=strstr(q," 0.")))
    {
        for (; q<p; q++) fputc (*q, stream);
        fputc (*(q++), stream);
        q++;
        fputc (*(q++), stream);
    }
    while ((*q)!=EOS) fputc (*(q++), stream);
    return;
} /* end vzfprintf() */
#undef VZFPRINTF_SIZE


void zfprintf (FILE *stream, char *format, ...)
{
    va_list ap;
    va_start (ap, format);
    vzfprintf (stream, format, ap);
    return;
} /* end zfprintf() */


void zprintf (char *format, ...)
{
    va_list ap;
    va_start (ap, format);
    vzfprintf (stdout, format, ap);
    return;
} /* end zprintf() */


#ifdef _zprintf_TEST
int main (int argc, char *argv[])
{
    printf  (" printf: %.16g %.16g\n", 1.5423543, 0.5224354325);
    zprintf ("zprintf: %.16g %.16g\n", 1.5423543, 0.5224354325);
    return (0);
}
#endif /* _zprintf_TEST */


/*************************************************************/
/* Starting at "series", there are "N" elements of "valtype" */
/* and separated in memory by regular "bytes_separation",    */
/* calculate the statistical properties of common interest.  */
/*************************************************************/
void CalculateSimpleStatistics
(int N, char *series, int bytes_separation, int valtype, SimpleStatistics *s)
{
    register int i;
    register char *p;
    double val;
    if (N < 1)
    {
        s->N = 0;
        s->idx_min = -1;
        s->min = 0;
        s->idx_max = -1;
        s->max = 0;
        s->average = 0;
        s->variance = 0;
        s->standard_deviation = 0;
        return;
    }
    else s->N = N;
    s->idx_min = 0;
    s->min = IOVAL(series,valtype);
    s->idx_max = 0;
    s->max = IOVAL(series,valtype);
    s->average = 0;
    s->variance = 0;
    for (i=0,p=series; i<N; i++,p+=bytes_separation)
    {
        val = IOVAL(p,valtype);
        if (val < s->min)
        {
            s->idx_min = i;
            s->min = val;
        }
        if (val > s->max)
        {
            s->idx_max = i;
            s->max = val;
        }
        s->average += val;
        s->variance += val*val;
    }
    s->average /= N;
    s->variance -= N*s->average*s->average;
    if (N>1) s->variance /= (N-1);
    else s->variance = 0;
    s->standard_deviation = sqrt(s->variance);
    return;
} /* end CalculateSimpleStatistics() */

#ifdef _CalculateSimpleStatistics_TEST
int main (int argc, char *argv[])
{
    typedef struct 
    {
        float val;
        char c;
    } CalculateSimpleStatisticsT;
    CalculateSimpleStatisticsT a[3] = {{1.,'a'},{2.,'b'},{3.,'d'}};
    SimpleStatistics s;
    CalculateSimpleStatistics
        (3, (char *)a, sizeof(CalculateSimpleStatisticsT), IOVAL_FLOAT, &s);
    printf ("N = %d\n", s.N);
    printf ("idx_min = %d\n", s.idx_min);
    printf ("min = %g\n", s.min);
    printf ("idx_max = %d\n", s.idx_max);
    printf ("max = %g\n", s.max);
    printf ("average = %g\n", s.average);
    printf ("variance = %g\n", s.variance);
    printf ("standard_deviation = %g\n", s.standard_deviation);
    return (0);
}
#endif /* _CalculateSimpleStatistics_TEST */
