/* H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */
/* H0 X                                                                            */
/* H0 X   libAtoms+QUIP: atomistic simulation library                              */
/* H0 X                                                                            */
/* H0 X   Portions of this code were written by                                    */
/* H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,      */
/* H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.          */
/* H0 X                                                                            */
/* H0 X   Copyright 2006-2010.                                                     */
/* H0 X                                                                            */
/* H0 X   These portions of the source code are released under the GNU General     */
/* H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html          */
/* H0 X                                                                            */
/* H0 X   If you would like to license the source code under different terms,      */
/* H0 X   please contact Gabor Csanyi, gabor@csanyi.net                            */
/* H0 X                                                                            */
/* H0 X   Portions of this code were written by Noam Bernstein as part of          */
/* H0 X   his employment for the U.S. Government, and are not subject              */
/* H0 X   to copyright in the USA.                                                 */
/* H0 X                                                                            */
/* H0 X                                                                            */
/* H0 X   When using this software, please cite the following reference:           */
/* H0 X                                                                            */
/* H0 X   http://www.libatoms.org                                                  */
/* H0 X                                                                            */
/* H0 X  Additional contributions by                                               */
/* H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras                 */
/* H0 X                                                                            */
/* H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */

/* required for isblank() */
#ifndef __USE_ISOC99
#define __USE_ISOC99
#endif

#include <stdio.h>
#include <stdint.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <libgen.h>

#include "libatoms.h"

#define LINESIZE 2048
#define MAX_ENTRY_COUNT 100
#define MAX_FIELD_COUNT 200
#define PROPERTY_STRING_LENGTH 10
#define PARAM_STRING_LENGTH 1024

/*
** Designation:  StriStr
**
** Call syntax:  char *stristr(char *String, char *Pattern)
**
** Description:  This function is an ANSI version of strstr() with
**               case insensitivity.
**
** Return item:  char *pointer if Pattern is found in String, else
**               pointer to 0
**
** Rev History:  16/07/97  Greg Thayer  Optimized
**               07/04/95  Bob Stout    ANSI-fy
**               02/03/94  Fred Cole    Original
**               09/01/03  Bob Stout    Bug fix (lines 40-41) per Fred Bulback
**
** Hereby donated to public domain.
*/

typedef unsigned int uint;

char *stristr(const char *String, const char *Pattern)
{
      char *pptr, *sptr, *start;

      for (start = (char *)String; *start != NULL; start++)
      {
            /* find start of pattern in string */
            for ( ; ((*start!=NULL) && (toupper(*start) != toupper(*Pattern))); start++)
                  ;
            if (NULL == *start)
                  return NULL;

            pptr = (char *)Pattern;
            sptr = (char *)start;

            while (toupper(*sptr) == toupper(*pptr))
            {
                  sptr++;
                  pptr++;

                  /* if end of pattern then pattern was found */

                  if (NULL == *pptr)
                        return (start);
            }
      }
      return NULL;
}

/* xyz_find_frames() 
 *
 * Find starting positions of xyz frames within a file
 * Uses a disk cache to save recomputing if xyz
 * file hasn't been modified. Returns number of frames.
 *
 * File format: text
 * First line:  number of frames (int)
 * Subsequent nframes+1 lines: offset (long), natoms (int)
 * Last offset is end of final frame scanned.
 */
void realloc_frames(long **frames, int **atoms, int *frames_array_size, int new_frames_array_size) {
  long *t_frames = NULL;
  int *t_atoms = NULL;
  int i;

  if (new_frames_array_size <= *frames_array_size) {
    return;
  }

  // if there's nothing, probably starting out, so assume a small array
  if (*frames_array_size == 0) {
    new_frames_array_size = (new_frames_array_size <= 0) ? 100 : new_frames_array_size;
  } else { // otherwise, need to expand
    if (new_frames_array_size <= 1.5*(*frames_array_size)) { // at least by factor of 1.5
      new_frames_array_size = 1.5*(*frames_array_size);
    }
  }

  // if there's data in current arrays
  if (*frames_array_size > 0) {
    // allocate temporaries
    t_frames = (long *) malloc (*frames_array_size * sizeof(long));
    t_atoms = (int *) malloc (*frames_array_size * sizeof(int));
    // copy data into temporaries
    for (i=0; i < *frames_array_size; i++) {
      t_frames[i] = (*frames)[i];
      t_atoms[i] = (*atoms)[i];
    }
    // free preexisting arrays
    free (*frames);
    free (*atoms);
  }
  // allocate new arrays
  *frames = (long *) malloc (new_frames_array_size * sizeof(long));
  *atoms = (int *) malloc (new_frames_array_size * sizeof(int));
  // copy data from temporaries
  for (i=0; i < *frames_array_size; i++) {
    (*frames)[i] = t_frames[i];
    (*atoms)[i] = t_atoms[i];
  }
  // update frames_array_size
  *frames_array_size = new_frames_array_size;

  // free temporaries
  if (t_frames != NULL) {
    free(t_frames);
    free(t_atoms);
  }
}

int xyz_find_frames(char *fname, long **frames, int **atoms, int *frames_array_size, int *error) {
  FILE *in, *index;
  char *bname;
  char indexname[LINESIZE], linebuffer[LINESIZE], buf1[LINESIZE], buf2[LINESIZE];
  int natoms, i, nframes;
  int from_scratch, do_update;
  struct stat xyz_stat, idx_stat;

  INIT_ERROR;

  natoms = 0;
  nframes = 0;
  strncpy(indexname, fname, LINESIZE);
  strncat(indexname, ".idx", LINESIZE-strlen(indexname)-1);

  if (stat(fname, &xyz_stat) != 0) {
    RAISE_ERROR_WITH_KIND(ERROR_IO, "Cannot stat xyz file %s\n", fname);
  }

  from_scratch = stat(indexname, &idx_stat) != 0;

  if (from_scratch) {
    // Try to read from current dir instead
    strncpy(buf1, indexname, LINESIZE);
    bname = basename(buf1);
    if (getcwd(buf2, LINESIZE) != NULL) {
      strncat(buf2, "/", LINESIZE-strlen(buf2)-1);
      strncat(buf2, bname, LINESIZE-strlen(buf2)-1);
      if (stat(buf2, &idx_stat) == 0) {
	fprintf(stderr,"Found index %s\n",buf2);
	strncpy(indexname,buf2, LINESIZE);
	from_scratch = 0;
      }
    }
  }

  do_update = xyz_stat.st_mtime > idx_stat.st_mtime;

  if (!from_scratch) {
    debug("xyz_find_frames: reading XYZ index from file %s\n", indexname);
    index = fopen(indexname, "r");
    if (index == NULL) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "Index file %s cannot be opened\n", indexname);
    }
    if (!fgets(linebuffer,LINESIZE,index)) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "Index file %s is empty\n",indexname);
    }
    sscanf(linebuffer, "%d", &nframes);
    realloc_frames(frames, atoms, frames_array_size, nframes+2);
    for (i=0; i<=nframes; i++) {
      if (!fgets(linebuffer,LINESIZE,index)) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "Premature end of indexfile %s\n",indexname);
      }
      sscanf(linebuffer, "%ld %d", &(*frames)[i], &(*atoms)[i]);
    }
    fclose(index);
  }

  if (from_scratch || do_update) {
    debug("xyz_find_frames: writing XYZ index to file %s\n", indexname);
    in = fopen(fname, "r");
    if (in == NULL) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "xyz_find_frames: cannot open %s for reading", fname);
    }
    if (from_scratch)
      nframes = 0;
    else {
      debug("xyz_find_frames: trying to update XYZ index... \n");

      // Try to seek past last frame, and check this looks
      // like the start of a new frame
      if (fseek(in,(*frames)[nframes],SEEK_SET) != 0 ||
	  !fgets(linebuffer,LINESIZE,in) ||
	  sscanf(linebuffer, "%d", &natoms) != 1 ||
	  natoms != (*atoms)[nframes]) {
	// Seek failed, we'll have to rebuild index from start
	fseek(in,0,SEEK_SET);
	nframes = 0;
	debug(" failed, rebuilding from scratch.\n");
      }
      else {
	// Rewind to start of frame
	fseek(in,(*frames)[nframes],SEEK_SET);
      }

      // TODO - improve check - fails if number of atoms has changed
    }

    debug("xyz_find_frames: starting to build index from file pos %ld nframes=%d\n", ftell(in), nframes);

    while (fgets(linebuffer,LINESIZE,in)) {
      realloc_frames(frames, atoms, frames_array_size, nframes+2);
      (*frames)[nframes] = ftell(in)-strlen(linebuffer);
      if (sscanf(linebuffer, "%d", &natoms) != 1) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "xyz_find_frames: malformed XYZ file %s at frame %d\n",fname,nframes);
      }

      (*atoms)[nframes] = natoms;
      nframes++;

      // Skip the whole frame, as quickly as possible
      for (i=0; i<natoms+1; i++)
	if (!fgets(linebuffer,LINESIZE,in)) {
	  fseek(in, (*frames)[nframes], SEEK_SET); // return file pointer to beginning of frame
	  nframes--; // incomplete last frame
	  goto XYZ_END;
	}
    }
  XYZ_END:
    if (nframes == 0) return 0;

    (*frames)[nframes] = ftell(in); // end of last frame in file
    (*atoms)[nframes] = natoms;
    index = fopen(indexname, "w");
    if (index == NULL) {
      // Try to write in current dir instead
      strncpy(buf1, indexname, LINESIZE);
      bname = basename(buf1);
      if (getcwd(buf2, LINESIZE) != NULL) {
	strncat(buf2, "/", LINESIZE-strlen(buf2)-1);
	strncat(buf2, bname, LINESIZE-strlen(buf2)-1);
	index = fopen(buf2, "w");
	debug("xyz_find_frames: writing index to %s\n", buf2);
      }
      if (index == NULL) {
	fclose(in);
	RAISE_ERROR_WITH_KIND(ERROR_IO, "Cannot write index file.\n");
      }
    } else
      debug("xyz_find_frames: writing index to %s\n", indexname);
    fprintf(index, "%d\n", nframes);
    for (i=0; i<=nframes; i++)
      fprintf(index, "%ld %d\n", (*frames)[i], (*atoms)[i]);
    fclose(in);
    fclose(index);
  }

  debug("xyz_find_frames: %s: found %d complete frames\n", fname, nframes);
  return nframes;
}


void query_xyz (char *filename, int compute_index, int frame, int *n_frame, int *n_atom, int *error)
{
  long *frames;
  int *atoms;
  int frames_array_size;

  INIT_ERROR;

  if (strcmp(filename, "stdout") == 0 || strcmp(filename, "stdin") == 0) {
    RAISE_ERROR_WITH_KIND(ERROR_IO, "query_xyz: cannot query special filename %s", filename);
  }

  *n_frame = 0;
  *n_atom = 0;

  if (!compute_index) return;

  frames_array_size = 0;
  *n_frame = xyz_find_frames(filename, &frames, &atoms, &frames_array_size, error);
  PASS_ERROR;

  if (*n_frame == 0) {
    RAISE_ERROR_WITH_KIND(ERROR_IO_EOF, "query_xyz: empty file", frame, *n_frame-1);
  }

  if (frame < 0 || frame >= *n_frame) {
    RAISE_ERROR_WITH_KIND(ERROR_IO, "query_xyz: frame %d out of range 0 <= frame < %d", frame, *n_frame-1);
  }
  *n_atom = atoms[frame];

  free(frames);
  free(atoms);
}


#define min(a,b) ((a) < (b) ? (a) : (b))


char* get_line(char *linebuffer, int string, int string_length, char *orig_stringp, char *stringp, char **prev_stringp,
		FILE *in, char *info, int *error)
{
  INIT_ERROR;

  if (string) {
    if (*stringp == '\0' || (string_length != 0 && (stringp-orig_stringp >= string_length))) {
      RAISE_ERROR_WITH_KIND(ERROR_IO_EOF, info);
    }
    *prev_stringp = stringp;
    while (*stringp != '\n' && *stringp != '\0' && (string_length == 0 || stringp-orig_stringp < string_length)) stringp++;
    strncpy(linebuffer, *prev_stringp, stringp-*prev_stringp);
    linebuffer[stringp-*prev_stringp] = '\0';
    //debug("line = <%s>\n", linebuffer);
    if (*stringp == '\n') stringp++;
    return stringp;
  } else {
    if (!fgets(linebuffer,LINESIZE,in)) {
      RAISE_ERROR_WITH_KIND(ERROR_IO_EOF, info);
    }
    linebuffer[strlen(linebuffer)-1] = '\0';
    //debug("line = <%s>\n", linebuffer);
    return NULL;
  }
}

#define GET_LINE(info) stringp = get_line(linebuffer, string, string_length, orig_stringp, stringp, &prev_stringp, in, info, error); PASS_ERROR;

void read_xyz (char *filename, fortran_t *params, fortran_t *properties, fortran_t *selected_properties, double lattice[3][3], int *n_atom,
	       int compute_index, int frame, int *range, int string, int string_length, int *error)
{
  FILE *in;
  int i,n, entry_count,j=0,k=0,ncols,m, atidx, at_start, at_end;
  char linebuffer[LINESIZE], tmpbuf[LINESIZE], param_key[LINESIZE], param_value[LINESIZE];
  char fields[MAX_FIELD_COUNT][LINESIZE], subfields[MAX_FIELD_COUNT][LINESIZE],
    finalfields[MAX_FIELD_COUNT][LINESIZE];
  char *p, *p1, tmp_logical, *orig_stringp, *prev_stringp, *stringp;
  int nxyz, nfields=0, offset, error_occured;
  double tmpd;
  int n_frame, n_selected;
  long *frames;
  int *atoms, type, shape[2], tmp_error, tmp_type, tmp_shape[2];
  int frames_array_size, got_index, n_buffer;
  void *data, *tmp_data;
  int property_type[MAX_ENTRY_COUNT], property_shape[MAX_ENTRY_COUNT][2], property_ncols[MAX_ENTRY_COUNT], n_property;
  void *property_data[MAX_ENTRY_COUNT];

  INIT_ERROR;

  if (strcmp(filename, "stdout") == 0) {
    RAISE_ERROR_WITH_KIND(ERROR_IO, "read_xyz: cannot open \"stdout\" for reading.");
  }

  n_buffer = 0;
  got_index = 0;
  orig_stringp = stringp = prev_stringp = NULL;
  in = NULL;
  if (string) {
    debug("read_xyz: reading from string\n");
    orig_stringp = stringp = filename;
  } else if (strcmp(filename, "stdin") == 0) {
    debug("read_xyz: reading from STDIN\n");
    in = stdin;
  } else if (compute_index) {
    // Not reading from stdin or from a string, so we can compute an index
    debug("read_xyz: computing index for file %s\n", filename);
    frames_array_size = 0;
    n_frame = xyz_find_frames(filename, &frames, &atoms, &frames_array_size, error);
    PASS_ERROR;

    if (frame < 0 || frame >= n_frame) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "read_xyz: frame %d out of range 0 <= frame < %d", frame, n_frame-1);
    }
    got_index = 1;
    in = fopen(filename, "r");
    if (in == NULL) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "read_xyz: cannot open file %s for reading", filename);
    }

    n_buffer = *n_atom;
    *n_atom = atoms[frame];
    if (fseek(in, frames[frame], SEEK_SET) == -1) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "cannot seek XYZ input file %s", filename);
    }
    free(frames);
    free(atoms);
  } else {
    // compute_index = 0, so we just open the file and start at the beginning
    in = fopen(filename, "r");
    if (in == NULL) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "read_xyz: cannot open file %s for reading", filename);
    }
  }

  if (!got_index && frame != 0 && in != stdin) {
    debug("read_xyz: skipping to frame %d\n", frame);

    for (i=0; i<frame-1; i++) {
      GET_LINE("read_xyz: premature end when skipping, expecting number of atoms");
      if (sscanf(linebuffer, "%d", &nxyz) != 1) {
  	RAISE_ERROR_WITH_KIND(ERROR_IO, "read_xyz: first line (%s) must be number of atoms when skipping frame %d", linebuffer, i);
      }
      GET_LINE("read_xyz: premature end when skipping, expecting comment line");
      for (j=0; j<nxyz; j++)
  	GET_LINE("read_xyz: premature end when skipping");
    }
  }

  GET_LINE("read_xyz: premature end, expecting number of atoms");

  if (sscanf(linebuffer, "%d", &nxyz) != 1) {
    RAISE_ERROR_WITH_KIND(ERROR_IO, "read_xyz: first line (%s) must be number of atoms", linebuffer);
  }

  if (got_index) {
    if (nxyz != *n_atom) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "read_xyz: mismatch in number of atoms - expecting %d but got %d\n", *n_atom, nxyz);
    }
  } else
    *n_atom = nxyz;

  // Have we been asked to read only a specific range of atom indices?
  if (range[0] != 0 && range[1] != 0) {
    if (range[0] < 1) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "read_xyz: lower limit of range (%d) must be >= 1", range[0]);
    }
    if (range[1] > *n_atom) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "read_xyz: upper limit of range (%d) must be <= %d", range[1], *n_atom);
    }
    if (range[1] <= range[0]) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "read_xyz: upper limit of range (%d) must be > lower limit (%d)", range[1], range[0]);
    }
    *n_atom = range[1] - range[0] + 1;
    at_start = range[0]-1;
    at_end = range[1]-1;
  }
  else {
    at_start = 0;
    at_end = *n_atom-1;
  }
  if (n_buffer < *n_atom)
      n_buffer = *n_atom;

  // Read comment line, which should contain 'Lattice=' and 'Properties=' keys
  GET_LINE("premature end - expecting comment line");

  if (!stristr(linebuffer, "lattice")) {
    // It's not an extended XYZ file. Try to guess what's going on.
    // If comment line contains nine or more fields, assume last nine are
    // lattice in cartesian coordinates.

    p = linebuffer;
    k = 0;
    while ((p1 = strsep(&p, " \t")) != NULL) {
      if (*p1 == '\0') continue;
      strncpy(fields[k++], p1, LINESIZE);
    }

    if (k >= 9) {
      offset = k-9;
      error_occured = 0;
      for (i=0; i<9; i++)
	if (sscanf(fields[offset+i], "%lf", &tmpd) != 1) {
	  error_occured = 1;
	  break;
	}

      if ((p = strstr(linebuffer, "\n")) != NULL) *p = '\0';
      if (!error_occured) {
	sprintf(tmpbuf, " Lattice=\"%s %s %s %s %s %s %s %s %s\"", fields[offset+0], fields[offset+1], fields[offset+2], fields[offset+3],
		fields[offset+4], fields[offset+5], fields[offset+6], fields[offset+7], fields[offset+8]);
	strncat(linebuffer, tmpbuf, LINESIZE-strlen(linebuffer)-1);

      } else {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "Cannot extract lattice from line %s\n", linebuffer);
      }
    } else {
      // Put in a bogus lattice
      sprintf(tmpbuf, " Lattice=\"0 0 0 0 0 0 0 0 0\"");
      strncat(linebuffer, tmpbuf, LINESIZE-strlen(linebuffer)-1);
    }
  }

  if (!stristr(linebuffer, "properties")) {
    // No Properties key. Add a default one.
    if ((p = strstr(linebuffer, "\n")) != NULL) *p = '\0';
    strncat(linebuffer, "Properties=species:S:1:pos:R:3",LINESIZE-strlen(linebuffer)-1);
  }

  // Parse parameters. First split on ", ', { or }
  p = linebuffer;
  k = 0;
  while ((p1 = strsep(&p, "\"'{}")) != NULL) {
    if (*p1 == '\0') continue;
    strncpy(fields[k++], p1, LINESIZE);
  }

  // Now split things outside quotes on whitespace
  nfields = 0;
  for (i=0; i<k; i++) {
    if (i % 2 == 0) {
      p = fields[i];
      j = 0;
      while ((p1 = strsep(&p, " \t")) != NULL) {
	if (*p1 == '\0') continue;
	strncpy(subfields[j++], p1, LINESIZE);
      }
      for (n=0; n<j; n++, nfields++) {
	strncpy(finalfields[nfields],subfields[n], LINESIZE);
      }

    } else {
      strncat(finalfields[nfields-1],fields[i],LINESIZE-strlen(finalfields[nfields-1])-1);
    }
  }

  // Finally, split on '=' to get key/value pairs
  for (i=0; i<nfields; i++) {
    strncpy(linebuffer, finalfields[i], LINESIZE);
    if ((p = strchr(linebuffer,'=')) == NULL) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "Badly formed key/value pair %s\n", linebuffer);
    }

    *p = '\0';
    strncpy(param_key, linebuffer, PARAM_STRING_LENGTH);
    strncpy(param_value, p+1, PARAM_STRING_LENGTH);
    // Now we have key in 'param_key' and value in 'param_value'

    //if (strcasecmp(param_key, "Lattice") == 0 ||
    //strcasecmp(param_key, "Properties") == 0) continue;

    strncpy(linebuffer, param_value, LINESIZE);
    k = 0;
    p = linebuffer;
    while ((p1 = strsep(&p, " ")) != NULL) {
      if (*p1 == '\0') continue;
      strncpy(fields[k++], p1, LINESIZE);
    }
    if (k == 0) {
      k = 1;
      strncpy(fields[0], linebuffer, LINESIZE);
    }

    debug("read_xyz: param key=%s value=%s k=%d\n", param_key, param_value, k);

    for (j=0; j<k; j++) {
      n=0;
      while(n < strlen(fields[j]) && isblank(fields[j][n])) n++;  // skip leading blanks
      if (n < strlen(fields[j]) && fields[j][n] == '-') n++;        // leading minus sign is OK
      for (; n<strlen(fields[j]); n++)
	if (!isdigit(fields[j][n])) goto NOT_INT;
    }

    if (k==1) {
      type = T_INTEGER;
      shape[0] = 1;
      shape[1] = 1;
    } else if (k == 3) {
      type = T_INTEGER_A;
      shape[0] = 3;
      shape[1] = 1;
    } else if (k == 9) {
      type = T_INTEGER_A2;
      shape[0] = 3;
      shape[1] = 3;
    } else {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "Bad number of fields %d in integer parameter %s\n", k, param_key);
    }
    goto FORMAT_DONE;

  NOT_INT:
    for (j=0; j<k; j++)
      if (strtod(fields[j], &p), strlen(p) != 0) goto NOT_REAL;

    if (k==1) {
      type = T_REAL;
      shape[0] = 1;
      shape[1] = 1;
    } else if (k == 3) {
      type = T_REAL_A;
      shape[0] = 3;
      shape[1] = 1;
    } else if (k == 9) {
      type = T_REAL_A2;
      shape[0] = 3;
      shape[1] = 3;
    } else {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "Bad number of fields %d in real parameter %s\n", k, param_key);
    }
    goto FORMAT_DONE;

  NOT_REAL:
    for (j=0; j<k; j++)
      if (strcmp(fields[j],"F") != 0 && strcmp(fields[j],"T") != 0)
	goto NOT_LOGICAL;

    if (k==1) {
      type = T_LOGICAL;
      shape[0] = 1;
      shape[1] = 1;
    } else if (k == 3) {
      type = T_LOGICAL_A;
      shape[0] = 3;
      shape[1] = 1;
    } else {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "Bad number of fields %d in logical parameter %s\n", k, param_key);
    }
    goto FORMAT_DONE;

  NOT_LOGICAL:
    // Fallback option: treat as a single string
    type = T_CHAR;
    shape[0] = strlen(param_value);
    goto FORMAT_DONE;

  FORMAT_DONE:
    if ((type == T_INTEGER_A || type == T_REAL_A || type == T_LOGICAL_A) && shape[0] != 3) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "Parameter %s must have size 1 or 3, but got %d\n", param_key, shape[0]);
    }

    // Create key in Fortran dictionary
    debug("read_xyz: adding key %s type=%d shape=[%d %d]\n", param_key, type, shape[0], shape[1]);
    dictionary_add_key(params, param_key, &type, shape, &data, error, strlen(param_key));
    PASS_ERROR;

    if (data == NULL) continue;

    switch(type) {
    case(T_INTEGER):
      INTEGER(data) = strtol(param_value, &p, 10);
      break;
    case(T_REAL):
      REAL(data) = strtod(param_value, &p);
      break;
    case(T_LOGICAL):
      if (param_value[0] == 'T')
	LOGICAL(data) = 1;
      else
	LOGICAL(data) = 0;
      break;
    case(T_INTEGER_A):
      for (m=0; m<k; m++)
	INTEGER_A(data,m) = strtol(fields[m], &p, 10);
      break;
    case(T_REAL_A):
      for (m=0; m<k; m++)
	REAL_A(data,m) = strtod(fields[m], &p);
      break;
    case(T_LOGICAL_A):
      for (m=0; m<k; m++) {
	if (fields[m][0] == 'T')
	  LOGICAL_A(data,m) = 1;
	else
	  LOGICAL_A(data,m) = 0;
      }
      break;
    case(T_INTEGER_A2):
      for (m=0; m<shape[1]; m++)
	for (n=0; n<shape[0]; n++)
	  INTEGER_A2(data,shape,n,m) = strtol(fields[3*m+n], &p, 10);
      break;
    case(T_REAL_A2):
      for (m=0; m<shape[1]; m++)
	for (n=0; n<shape[0]; n++)
	  REAL_A2(data,shape,n,m) = strtod(fields[3*m+n], &p);
      break;
    case(T_CHAR):
      strncpy(CHAR(data), param_value, strlen(param_value));
      break;
    default:
      RAISE_ERROR_WITH_KIND(ERROR_IO, "Unknown param type %d\n", type);
    }
  }

  // Now parse Properties string
  memset(param_key, ' ', LINESIZE);
  strncpy(param_key, "Properties", strlen("Properties"));
  dictionary_query_key(params, param_key, &type, shape, &data, error, strlen("Properties"));
  PASS_ERROR;
  strncpy(linebuffer, (char *)data, shape[0]);
  linebuffer[shape[0]] = '\0';

  debug("properties string %s\n", linebuffer);

  p = linebuffer;
  k = 0;
  while ((p1 = strsep(&p, ":")) != NULL) {
    if (k >= MAX_FIELD_COUNT) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "Maximum field count (%d) exceeded", MAX_FIELD_COUNT);
    }
    strncpy(fields[k++], p1, LINESIZE);
  }

  entry_count = 0;
  n_property = 0;
  n_selected = 0;
  if (selected_properties != NULL) dictionary_get_n(selected_properties, &n_selected);

  for (i=0; i<k/3; i++) {
    debug("read_xyz: got property %s:%s:%s\n", fields[3*i], fields[3*i+1], fields[3*i+2]);

    if (sscanf(fields[3*i+2], "%d", &ncols) != 1) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "Bad column count %s line=%s\n", fields[3*i+2], linebuffer);
    }

    entry_count += ncols;
    if (entry_count > MAX_ENTRY_COUNT) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "Maximum entry count(%d) exceeded\n", MAX_ENTRY_COUNT);
    }

    if (strcmp(fields[3*i+1],"I") == 0) {
      if (ncols == 1) {
	type = T_INTEGER_A;
	shape[0] = n_buffer;
      } else {
	type = T_INTEGER_A2;
	shape[0] = ncols;
	shape[1] = n_buffer;
      }
    } else if (strcmp(fields[3*i+1],"R") == 0) {
      if (ncols == 1) {
	type = T_REAL_A;
	shape[0] = n_buffer;
      } else {
	type = T_REAL_A2;
	shape[0] = ncols;
	shape[1] = n_buffer;
      }
    } else if (strcmp(fields[3*i+1],"S") == 0) {
      if (ncols != 1) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "String property %s with ncols != 1 no longer supported", fields[3*i]);
      }
      type = T_CHAR_A;
      shape[0] = PROPERTY_STRING_LENGTH; // FIXME: this could be variable
      shape[1] = n_buffer;
    } else if (strcmp(fields[3*i+1],"L") == 0) {
      if (ncols != 1) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "Logical property %s with ncols != 1 no longer supported", fields[3*i]);
      }
      type = T_LOGICAL_A;
      shape[0] = n_buffer;
    } else  {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "Bad property type %s\n", fields[3*i+1]);
    }

    property_type[n_property] = type;
    property_data[n_property] = NULL;
    property_ncols[n_property] = ncols;
    property_shape[n_property][0] = shape[0];
    property_shape[n_property][1] = shape[1];

    tmp_error = ERROR_NONE;
    if (n_selected != 0) {
      // Check selected_properties to see if should skip this property?
      dictionary_query_key(selected_properties, fields[3*i], &tmp_type, tmp_shape, &tmp_data, &tmp_error, strlen(fields[3*i]));
      CLEAR_ERROR;
    }
    if (tmp_error == ERROR_NONE) {
      dictionary_add_key(properties, fields[3*i], &type, shape, &data, error, strlen(fields[3*i]));
      if (data != NULL) {
	if (type == T_CHAR_A) memset(data, ' ', shape[0]*shape[1]); // Zero string data
	property_data[n_property] = data;
      }
    }
    n_property++;
  }

  // Read lattice
  memset(param_key, ' ', LINESIZE);
  strncpy(param_key, "Lattice", strlen("Lattice"));
  dictionary_query_key(params, param_key, &type, shape, &data, error, strlen("Lattice"));
  PASS_ERROR;

  if (type == T_REAL_A2) {
    for (m=0; m<shape[0]; m++)
      for (n=0; n<shape[1]; n++)
	lattice[m][n] = REAL_A2(data, shape, n, m);
  } 
   else if (type == T_INTEGER_A2) {
    for (m=0; m<shape[0]; m++)
      for (n=0; n<shape[1]; n++)
  	lattice[m][n] = INTEGER_A2(data, shape, n, m);
  } else {
    RAISE_ERROR_WITH_KIND(ERROR_IO, "read_xyz: bad type for Lattice (%d)", type);
  }

  // Now it's just one line per atom
  n = 0;
  for (atidx=0; atidx < nxyz; atidx++) {
    GET_LINE("premature file ending");

    if (atidx < at_start || atidx > at_end) continue;

    k = 0;
    p = linebuffer;
    while ((p1 = strsep(&p, " \t\n")) != NULL) {
      if (*p1 == '\0') continue;
      strncpy(fields[k++], p1, LINESIZE);
    }
    if (k != entry_count) {
      for (i=0;i<k;i++) fprintf(stderr, "fields[%d] = %s, length %lu\n", i, fields[i], (unsigned long)strlen(fields[i]));
      RAISE_ERROR_WITH_KIND(ERROR_IO, "incomplete row, frame %d atom %d - got %d/%d entries\n", frame, n, k, entry_count, frame, linebuffer);
    }

    k = 0;
    for (i=0; i<n_property; i++) {
      if (property_data[i] == NULL) {
	k += property_ncols[i];
	continue;
      }
      switch(property_type[i]) {
      case(T_INTEGER_A):
	if (sscanf(fields[k], "%d", &INTEGER_A(property_data[i], n)) != 1)  {
	  RAISE_ERROR_WITH_KIND(ERROR_IO, "Can't convert int value %s\n", fields[k]);
	}
	k++;
	break;

      case(T_INTEGER_A2):
	for (j=0; j < property_ncols[i]; j++)
	  if (sscanf(fields[k+j], "%d", &INTEGER_A2(property_data[i], property_shape[i], j, n)) != 1)  {
	    RAISE_ERROR_WITH_KIND(ERROR_IO, "Can't convert int value %s\n", fields[k+j]);
	  }
	k += property_shape[i][0];
	break;

      case(T_REAL_A):
	if (sscanf(fields[k], "%lf", &REAL_A(property_data[i], n)) != 1)  {
	  RAISE_ERROR_WITH_KIND(ERROR_IO, "Can't convert int value %s\n", fields[k+j]);
	}
	k++;
	break;

      case(T_REAL_A2):
	for (j=0; j < property_ncols[i]; j++)
	  if (sscanf(fields[k+j], "%lf", &REAL_A2(property_data[i], property_shape[i], j, n)) != 1)  {
	    RAISE_ERROR_WITH_KIND(ERROR_IO, "Can't convert real value %s\n", fields[k+j]);
	  }
	k += property_shape[i][0];
	break;

      case(T_LOGICAL_A):
	if (sscanf(fields[k], "%c", &tmp_logical) != 1)  {
	  RAISE_ERROR_WITH_KIND(ERROR_IO, "Can't convert logical value %s\n", fields[k+j]);
	}
	if (tmp_logical == 'T' || tmp_logical == '1')
	  LOGICAL_A(property_data[i],n) = 1;
	else
	  LOGICAL_A(property_data[i],n) = 0;
	k++;
	break;

      case(T_CHAR_A):
	if (sscanf(fields[k], "%10c", &CHAR_A(property_data[i], property_shape[i], 0, n)) != 1)  {
	  RAISE_ERROR_WITH_KIND(ERROR_IO, "Can't convert str value %s\n", fields[k]);
	}
	k++;
	break;

      default:
	RAISE_ERROR_WITH_KIND(ERROR_IO, "Bad property type %d", property_type[i]);
      }
    }
    n++;
  }
  if (!string && in != stdin) fclose(in);
}

#define PUT_LINE(line) { if (string) extendable_str_concat(estr, line, &tmp_zero, &tmp_one, strlen(line)-1); else fputs(line, out); }

void write_xyz (char *filename, fortran_t *params, fortran_t *properties, fortran_t *selected_properties, double lattice[3][3], int n_atom,
		int append, char *prefix, char *int_format, char *real_format, char *str_format, char *logical_format, 
		int string, fortran_t *estr, int *error) {
  FILE *out;
  char linebuffer[LINESIZE], tmpbuf[LINESIZE], param_key[LINESIZE], param_value[LINESIZE], property_name[C_KEY_LEN];
  int i, j, m, n, type, shape[2], tmp_type, tmp_shape[2];
  char *trimmed;
  int property_type[MAX_ENTRY_COUNT], property_shape[MAX_ENTRY_COUNT][2], n_property, ncols;
  void *property_data[MAX_ENTRY_COUNT];
  char property_code;
  void *data, *tmp_data;
  int tmp_zero = 0, tmp_one = 1;

  INIT_ERROR;

  if (strcmp(filename, "stdin") == 0) {
    RAISE_ERROR_WITH_KIND(ERROR_IO, "write_xyz: cannot open \"stdin\" for writing.");
  }

  debug("write_xyz: string=%d\n", string);
  if (!string) {
    if (strcmp(filename, "stdout") == 0) {
      out = stdout;
    } else {
      if (append)
	out = fopen(filename, "a");
      else
	out = fopen(filename, "w");
      if (out == NULL) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "write_xyz: cannot open %s for writing", filename);
      }
    }
  } else
    out = NULL;

  // Number of atoms
  if (prefix != NULL && prefix[0] != '\0')
    sprintf(linebuffer, "%s %d\n", prefix, n_atom);
  else
    sprintf(linebuffer, "%d\n", n_atom);
  PUT_LINE(linebuffer);

  memset(param_key, ' ', LINESIZE);
  strncpy(param_key, "Lattice", strlen("Lattice"));
  type = T_REAL_A2;
  shape[0] = 3;
  shape[1] = 3;
  dictionary_add_key(params, param_key, &type, shape, &data, error, strlen("Lattice"));
  PASS_ERROR;

  for (m=0; m<shape[0]; m++)
    for (n=0; n<shape[1]; n++)
      REAL_A2(data, shape, n, m) = lattice[m][n];

  dictionary_get_n(selected_properties, &n_property);
  memset(param_value, ' ', LINESIZE);
  param_value[0] = '\0';
  for (i=1; i<=n_property; i++) {
    dictionary_query_index(selected_properties, &i, property_name, &tmp_type, tmp_shape, &tmp_data, error, C_KEY_LEN);
    PASS_ERROR;

    dictionary_query_key(properties, property_name, &property_type[i], shape, &property_data[i], error, C_KEY_LEN);
    PASS_ERROR;

    // null-terminate the Fortran string
    property_name[C_KEY_LEN-1] = '\0';
    if (strchr(property_name, ' ') == NULL) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "write_xyz: property name %s not terminated with blank", property_name);
    }
    *strchr(property_name, ' ') = '\0';

    property_code = ' ';
    switch(property_type[i]) {
    case(T_INTEGER_A):
    case(T_INTEGER_A2):
      property_code = 'I';
      break;

    case(T_REAL_A):
    case(T_REAL_A2):
      property_code = 'R';
      break;

    case(T_CHAR_A):
      property_code = 'S';
      break;

    case(T_LOGICAL_A):
      property_code = 'L';
      break;

    default:
      RAISE_ERROR_WITH_KIND(ERROR_IO, "write_xyz: bad property type %d for key %s", type, property_name);
    }

    property_shape[i][0] = shape[0];
    property_shape[i][1] = shape[1];
    if (property_type[i] == T_INTEGER_A || property_type[i] == T_REAL_A ||
	property_type[i] == T_LOGICAL_A || property_type[i] == T_CHAR_A)
      ncols = 1;
    else {
      ncols = property_shape[i][0];
    }

    sprintf(tmpbuf, "%s:%c:%d:", property_name, property_code, ncols);
    strncat(param_value, tmpbuf, PARAM_STRING_LENGTH-strlen(param_value)-1);
  }
  param_value[strlen(param_value)-1] = ' '; // remove trailing :
  debug("properties string <%s>\n", param_value);

  memset(param_key, ' ', LINESIZE);
  strncpy(param_key, "Properties", strlen("Properties"));
  type = T_CHAR;
  shape[0] = strlen(param_value)-1;
  dictionary_add_key(params, param_key, &type, shape, &data, error, strlen("Properties"));
  PASS_ERROR;
  memset(data, ' ', shape[0]);
  strncpy(CHAR(data), param_value, strlen(param_value)-1);

  // Build parameter values

  dictionary_get_n(params, &n);
  param_value[0] = '\0';
  linebuffer[0] = '\0';
  if (prefix != NULL && prefix[0] != '\0') sprintf(linebuffer, "%s ", prefix);
  for (i=1; i<=n; i++) {
    dictionary_query_index(params, &i, param_key, &type, shape, &data, error, C_KEY_LEN);
    PASS_ERROR;

    // null-terminate the Fortran string
    param_key[C_KEY_LEN-1] = '\0';
    if (strchr(param_key, ' ') == NULL) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "write_xyz: key %s not terminated with blank", param_key);
    }
    *strchr(param_key, ' ') = '\0';

    if (type == T_INTEGER)
      sprintf(param_value, int_format, INTEGER(data));
    else if (type == T_REAL)
      sprintf(param_value, real_format, REAL(data));
    else if (type == T_LOGICAL) {
      if (LOGICAL(data))
	sprintf(param_value, "T");
      else
	sprintf(param_value, "F");
    }
    else if (type == T_INTEGER_A) {
      if (shape[0] != 3) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "Parameter %s integer array size=%d != 3", param_key, shape[0]);
      }
      param_value[0]='\0';
      for (j=0; j<shape[0]; j++) {
	sprintf(tmpbuf, int_format, INTEGER_A(data,j));
	strncat(param_value, tmpbuf, PARAM_STRING_LENGTH-strlen(param_value)-1);
      }
    } else if (type == T_REAL_A) {
      if (shape[0] != 3) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "Parameter %s real array size=%d != 3", param_key, shape[0]);
      }
      param_value[0]='\0';
      for (j=0; j<shape[0]; j++) {
	sprintf(tmpbuf, real_format, REAL_A(data,j));
	strncat(param_value, tmpbuf, PARAM_STRING_LENGTH-strlen(param_value)-1);
      }
    } else if (type == T_LOGICAL_A) {
      if (shape[0] != 3) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "Parameter %s logical array size=%d != 3", param_key, shape[0]);
      }
      sprintf(param_value, "%s %s %s",
	      LOGICAL_A(data,0) ? "T" : "F",
	      LOGICAL_A(data,1) ? "T" : "F",
	      LOGICAL_A(data,2) ? "T" : "F");
    } else if (type == T_INTEGER_A2) {
      if (shape[0] != 3 || shape[1] != 3) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "Parameter %s 2D integer array shape=[%d %d] != [3 3]", param_key, shape[0], shape[1]);
      }
      sprintf(tmpbuf, "%s %s %s %s %s %s %s %s %s", int_format, int_format, int_format,
	      int_format, int_format, int_format, int_format, int_format, int_format);
      sprintf(param_value, tmpbuf,
	      INTEGER_A2(data,shape,0,0), INTEGER_A2(data,shape,1,0), INTEGER_A2(data,shape,2,0),
	      INTEGER_A2(data,shape,0,1), INTEGER_A2(data,shape,1,1), INTEGER_A2(data,shape,2,1),
	      INTEGER_A2(data,shape,0,2), INTEGER_A2(data,shape,1,2), INTEGER_A2(data,shape,2,2));
    } else if (type == T_REAL_A2) {
      sprintf(tmpbuf, "%s %s %s %s %s %s %s %s %s", real_format, real_format, real_format,
	      real_format, real_format, real_format, real_format, real_format, real_format);
      sprintf(param_value, tmpbuf,
	      REAL_A2(data,shape,0,0), REAL_A2(data,shape,1,0), REAL_A2(data,shape,2,0),
	      REAL_A2(data,shape,0,1), REAL_A2(data,shape,1,1), REAL_A2(data,shape,2,1),
	      REAL_A2(data,shape,0,2), REAL_A2(data,shape,1,2), REAL_A2(data,shape,2,2));
    } else if (type == T_CHAR) {
      memset(param_value, ' ', LINESIZE);
      strncpy(param_value, CHAR(data), shape[0]);
      param_value[LINESIZE-1] = '\0';
      while (param_value[strlen(param_value)-1] == ' ')
	param_value[strlen(param_value)-1] = '\0';
    }

    trimmed = param_value;
    while (isblank(trimmed[0])) trimmed++;

    sprintf(tmpbuf, "%s=%s%s%s ", param_key,
	    strchr(trimmed,' ') != NULL ? "\"" : "",
	    trimmed,
	    strchr(trimmed,' ') != NULL ? "\"" : "");
    strncat(linebuffer, tmpbuf, LINESIZE-strlen(linebuffer)-1);
  }

  linebuffer[strlen(linebuffer)-1] = '\n';
  PUT_LINE(linebuffer);

  for (n=0; n<n_atom; n++) {
    linebuffer[0] = '\0';
    if (prefix != NULL && prefix[0] != '\0') strncat(linebuffer, prefix, LINESIZE-strlen(linebuffer)-1);

    for (i=1; i<=n_property; i++) {

      switch(property_type[i]) {
      case(T_INTEGER_A):
	sprintf(tmpbuf, int_format, INTEGER_A(property_data[i], n));
	if (strlen(linebuffer) != 0 && linebuffer[strlen(linebuffer)-1] != ' ' && tmpbuf[0] != ' ') 
	  strncat(linebuffer, " ", LINESIZE-strlen(linebuffer)-1);
	strncat(linebuffer, tmpbuf, LINESIZE-strlen(linebuffer)-1);
	break;

      case(T_INTEGER_A2):
	for (j=0; j < property_shape[i][0]; j++) {
	  sprintf(tmpbuf, int_format, INTEGER_A2(property_data[i], property_shape[i], j, n));
	  if (strlen(linebuffer) != 0 && linebuffer[strlen(linebuffer)-1] != ' ' && tmpbuf[0] != ' ') 
	    strncat(linebuffer, " ", LINESIZE-strlen(linebuffer)-1);
	  strncat(linebuffer, tmpbuf, LINESIZE-strlen(linebuffer)-1);
	}
	break;

      case(T_REAL_A):
	sprintf(tmpbuf, real_format, REAL_A(property_data[i], n));
	if (strlen(linebuffer) != 0 && linebuffer[strlen(linebuffer)-1] != ' ' && tmpbuf[0] != ' ') 
	  strncat(linebuffer, " ", LINESIZE-strlen(linebuffer)-1);
	strncat(linebuffer, tmpbuf, LINESIZE-strlen(linebuffer)-1);
	break;

      case(T_REAL_A2):
	for (j=0; j < property_shape[i][0]; j++) {
	  sprintf(tmpbuf, real_format, REAL_A2(property_data[i], property_shape[i], j, n));
	  if (strlen(linebuffer) != 0 && linebuffer[strlen(linebuffer)-1] != ' ' && tmpbuf[0] != ' ') 
	    strncat(linebuffer, " ", LINESIZE-strlen(linebuffer)-1);
	  strncat(linebuffer, tmpbuf, LINESIZE-strlen(linebuffer)-1);
	}
	break;

      case(T_CHAR_A):
	sprintf(tmpbuf, str_format, (char *)property_data[i] + property_shape[i][0]*n);
	if (strlen(linebuffer) != 0 && linebuffer[strlen(linebuffer)-1] != ' ' && tmpbuf[0] != ' ') 
	  strncat(linebuffer, " ", LINESIZE-strlen(linebuffer)-1);
	strncat(linebuffer, tmpbuf, LINESIZE-strlen(linebuffer)-1);
	break;

      case(T_LOGICAL_A):
	sprintf(tmpbuf, logical_format, LOGICAL_A(property_data[i], n) ? 'T' : 'F');
	if (strlen(linebuffer) != 0 && linebuffer[strlen(linebuffer)-1] != ' ' && tmpbuf[0] != ' ') 
	  strncat(linebuffer, " ", LINESIZE-strlen(linebuffer)-1);
	strncat(linebuffer, tmpbuf, LINESIZE-strlen(linebuffer)-1);
	break;
      }
    }
    strncat(linebuffer, "\n", LINESIZE-strlen(linebuffer)-1);
    PUT_LINE(linebuffer);
  }

  if (!string) {
    if (out != stdout) 
      fclose(out);
    else
      fflush(out);
  }
}

