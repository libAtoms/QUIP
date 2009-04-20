/* atomeye C api to be called from Python */

int atomeyelib_main(int argc, char* argv[], void (*onclick)(int atom), 
		    void (*onclose)(), void (*onadvance)(char *instr), int *done_init);
int atomeyelib_run_command(int iw, char *line, char **outstr);
int atomeyelib_help(int iw, char *line, char **outstr);
int atomeyelib_redraw(int iw);
int atomeyelib_close(int iw);
int atomeyelib_load_libatoms(int iw, Atoms *atoms, char *title, char **outstr);

/* int atomeyelib_set_output(int on_off); */
