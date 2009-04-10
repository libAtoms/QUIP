/***************************************/
/* libIO: -lreadline -ltermcap         */
/*                                     */
/* Strings, parsers & files            */
/*                                     */
/* Dec.12, 1999 Ju Li <liju99@mit.edu> */
/***************************************/

#include "IO.h"

/* Interface to GNU Readline Library */

static char *line_read = (char *)NULL;

/* Read a string, and return a pointer to it. Returns NULL on EOF. */
char *readline_gets (char *prompt, char *default_answer)
{
    char *real_prompt, *answer;
    if (line_read)
    {
        free (line_read);
        line_read = (char *)NULL;
    }
    if (default_answer)
    {
        MALLOC ( readline_gets, real_prompt,
                 strlen(prompt)+strlen(default_answer)+6, char);
    }
    else MALLOC ( readline_gets, real_prompt, strlen(prompt)+3, char);
    if (prompt)
    {
        if (default_answer)
            sprintf (real_prompt, "%s (%s): ", prompt, default_answer);
        else sprintf (real_prompt, "%s: ", prompt);
    }
    else
    {
        if (default_answer)
            sprintf (real_prompt, "(%s): ", default_answer);
        else sprintf (real_prompt, ": ");
    }
    line_read = readline(real_prompt);
    free (real_prompt);
    if (line_read)
    {
        if (*line_read)
        {
            add_history(answer=line_read);
            real_prompt = eos(line_read)-1;
            if (*real_prompt==' ') *real_prompt=EOS;
        }
        else if (default_answer) add_history(answer=default_answer);
        else answer=line_read;
    }
    else answer = NULL;
    return (answer);
} /* end readline_gets() */


#ifdef _readline_gets_TEST
int main (int argc, char *argv[])
{
    char *buf;
    buf = readline_gets("What are the files", "readline.c");
    buf = readline_gets("Again, what are the files", "readline.c");
    printf ("Your selection is: \"%s\"\n", buf);
    return (0);
}
#endif /* _readline_gets_TEST */


#ifdef _readlineSelectfile
int main (int argc, char *argv[])
{
    char *buf;
    if (argc == 1)
    {
        printf ("\nPurpose: Use GNU readline to select file.\n\n");
        printf ("Usage: %s 'Select file'\n", argv[0]);
        printf ("       %s 'Select file' 'Miranda02.pdf'\n\n", argv[0]);
        return (1);
    }
    else if (argc == 2) buf = readline_gets(argv[1], NULL);
    else buf = readline_gets(argv[1], argv[2]);
    printf ("%s\n", buf);
    return (0);
}
#endif /* _readlineSelectfile */
