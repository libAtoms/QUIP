/********************************************/
/* libAX: -lX11 -lXext -lpng -lz -ljpeg -lm */
/*        -lScalar -lIO -lTimer             */
/*                                          */
/* Accelerated low-level graphics library   */
/* with X Window Shared Memory Extension.   */
/*                                          */
/* Jan.13, 2000 Ju Li <liju99@mit.edu>      */
/********************************************/

int jpg2eps (char *filename_jpg, char *filename_eps, char *options);

/**********************************************/
/* Bundled jpeg2ps (C) 1994-1999 Thomas Merz  */
/* http://www.pdflib.com/jpeg2ps/. See also   */
/* /usr/share/ghostscript/6.0/lib/viewjpeg.ps */
/**********************************************/

#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <ctype.h>

/* -------------------------- psimage.h ------------------------- */
/*  (C) 1994-1999 Thomas Merz                                     */
/* -------------------------- psimage.h ------------------------- */

#ifndef TRUE
#define TRUE  1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#ifdef min
#undef min
#endif
#define min(a, b) ((a) < (b) ? (a) : (b))

typedef int BOOL;

/* Simple prototype macros for K&R and ANSI-C */

#ifdef KNR
#define P0(v)			()
#define P1(t1, p1)		(p1) t1 p1
#define P2(t1, p1, t2, p2)	(p1, p2) t1 p1; t2 p2
#define P3(t1, p1, t2, p2, t3, p3)	(p1, p2, p3) t1 p1; t2 p2; t3 p3
#else
#define P0(v)			(void)
#define P1(t1, p1)		(t1 p1)
#define P2(t1, p1, t2, p2)	(t1 p1, t2 p2)
#define P3(t1, p1, t2, p2, t3, p3)	(t1 p1, t2 p2, t3 p3)
#endif

/* data output mode: binary, ascii85, hex-ascii */
typedef enum { BINARY, ASCII85, ASCIIHEX } DATAMODE;

typedef struct {
  FILE     *fp;                   /* file pointer for jpeg file		 */
  char     *filename;             /* name of image file			 */
  int      width;                 /* pixels per line			 */
  int      height;                /* rows				 */
  int      components;            /* number of color components		 */
  int      bits_per_component;    /* bits per color component		 */
  float    dpi;                   /* image resolution in dots per inch   */
  DATAMODE mode;                  /* output mode: 8bit, ascii, ascii85	 */
  long     startpos;              /* offset to jpeg data		 */
  BOOL     landscape;             /* rotate image to landscape mode?	 */
  BOOL     adobe;                 /* image includes Adobe comment marker */
} imagedata;

#define	DPI_IGNORE	-1.0      /* dummy value for imagedata.dpi       */
#define DPI_USE_FILE     0.0      /* dummy value for imagedata.dpi       */

int Margin     	= 20;           /* safety margin */
/* Modified Fri Mar 10 15:22:29 2000 for library call */
/* BOOL quiet	= FALSE; */
/* suppress informational messages */
BOOL quiet	= TRUE;
BOOL autorotate = FALSE;	/* disable automatic rotation */

/* --------------------------- asc85ec.c --------------------------- */
/* ASCII85 and Hex encoding for PostScript Level 2 and PDF           */
/* (C) Thomas Merz 1994-99                                           */
/* --------------------------- asc85ec.c --------------------------- */

typedef unsigned char byte;
static unsigned char buf[4];
static unsigned long power85[5] = { 1L, 85L, 85L*85, 85L*85*85, 85L*85*85*85};
static int outbytes;		/* Number of characters in an output line */

/* read 0-4 Bytes. result: number of bytes read */
static int ReadSomeBytes P1(FILE *, in)
{
  register int count, i;

  for (count = 0; count < 4; count++) {
    if ((i = getc(in)) == EOF)
      break;
    else
      buf[count] = (byte) i;
  }
  return count;
}

/* Two percent characters at the start of a line will cause trouble
 * with some post-processing software. In order to avoid this, we
 * simply insert a line break if we encounter a percent character
 * at the start of the line. Of course, this rather simplistic
 * algorithm may lead to a large line count in pathological cases,
 * but the chance for hitting such a case is very small, and even
 * so it's only a cosmetic flaw and not a functional restriction.
 */
static void outbyte P2(byte, c, FILE *, out)
{    /* output one byte */

  if (fputc(c, out) == EOF) {
    fprintf(stderr, "jpeg2ps: write error - exit!\n");
    exit(1);
  }

  if (++outbytes > 63 ||		/* line limit reached */
    (outbytes == 1 && c == '%') )
  {	/* caution: percent character at start of line */
      fputc('\n', out);			/* insert line feed */
      outbytes = 0;
  }
}

int ASCII85Encode P2(FILE *, in, FILE *, out)
{
  register int i, count;
  unsigned long word, v;

  outbytes = 0;

  /* 4 bytes read ==> output 5 bytes */
  while ((count = ReadSomeBytes(in)) == 4) {
    word = ((unsigned long)(((unsigned int)buf[0] << 8) + buf[1]) << 16) +
                   (((unsigned int)buf[2] << 8) + buf[3]);
    if (word == 0)
      outbyte('z', out);       /* shortcut for 0 */
    else
      /* calculate 5 ASCII85 bytes and output them */
      for (i = 4; i >= 0; i--) {
	v = word / power85[i];
	outbyte((byte) (v + '!'), out);
	word -= v * power85[i];
      }
  }

  word = 0;

  if (count != 0) {   /* 1-3 bytes left */
    for (i = count-1; i >= 0; i--)   /* accumulate bytes */
      word += (unsigned long)buf[i] << 8 * (3-i);
    
    /* encoding as above, but output only count+1 bytes */
    for (i = 4; i >= 4-count; i--) {
      v = word / power85[i];
      outbyte((byte) (v + '!'), out);
      word -= v * power85[i];
    }
  }

  fputc('~', out);	/* EOD marker */
  fputc('>', out);
  return 0;
}

void ASCIIHexEncode P2(FILE *, in, FILE *, out)
{
  static char buffer[512];
  static char BinToHex[] = "0123456789ABCDEF";
  int CharsPerLine;
  size_t i, n;
  unsigned char *p;

  CharsPerLine = 0;
  fputc('\n', out);

  while ((n = fread(buffer, 1, sizeof(buffer), in)) != 0)
    for (i = 0, p = (unsigned char *) buffer; i < n; i++, p++) {
      fputc(BinToHex[*p>>4], out);           /* first nibble  */
      fputc(BinToHex[*p & 0x0F], out);       /* second nibble */
      if ((CharsPerLine += 2) >= 64) {
        fputc('\n', out);
        CharsPerLine = 0;
      }
    }

  fputc('>', out);         /* EOD marker for PostScript hex strings */
}


/* ------------------------- readjpeg.c ------------------------- */
/*  (C) Thomas Merz 1994-99                                       */
/* ------------------------- readjpeg.c ------------------------- */

/* The following enum is stolen from the IJG JPEG library
 * Comments added by tm
 * This table contains far too many names since jpeg2ps
 * is rather simple-minded about markers
 */

/* extern BOOL quiet; */

typedef enum {		/* JPEG marker codes			*/
  M_SOF0  = 0xc0,	/* baseline DCT				*/
  M_SOF1  = 0xc1,	/* extended sequential DCT		*/
  M_SOF2  = 0xc2,	/* progressive DCT			*/
  M_SOF3  = 0xc3,	/* lossless (sequential)		*/
  
  M_SOF5  = 0xc5,	/* differential sequential DCT		*/
  M_SOF6  = 0xc6,	/* differential progressive DCT		*/
  M_SOF7  = 0xc7,	/* differential lossless		*/
  
  M_JPG   = 0xc8,	/* JPEG extensions			*/
  M_SOF9  = 0xc9,	/* extended sequential DCT		*/
  M_SOF10 = 0xca,	/* progressive DCT			*/
  M_SOF11 = 0xcb,	/* lossless (sequential)		*/
  
  M_SOF13 = 0xcd,	/* differential sequential DCT		*/
  M_SOF14 = 0xce,	/* differential progressive DCT		*/
  M_SOF15 = 0xcf,	/* differential lossless		*/
  
  M_DHT   = 0xc4,	/* define Huffman tables		*/
  
  M_DAC   = 0xcc,	/* define arithmetic conditioning table	*/
  
  M_RST0  = 0xd0,	/* restart				*/
  M_RST1  = 0xd1,	/* restart				*/
  M_RST2  = 0xd2,	/* restart				*/
  M_RST3  = 0xd3,	/* restart				*/
  M_RST4  = 0xd4,	/* restart				*/
  M_RST5  = 0xd5,	/* restart				*/
  M_RST6  = 0xd6,	/* restart				*/
  M_RST7  = 0xd7,	/* restart				*/
  
  M_SOI   = 0xd8,	/* start of image			*/
  M_EOI   = 0xd9,	/* end of image				*/
  M_SOS   = 0xda,	/* start of scan			*/
  M_DQT   = 0xdb,	/* define quantization tables		*/
  M_DNL   = 0xdc,	/* define number of lines		*/
  M_DRI   = 0xdd,	/* define restart interval		*/
  M_DHP   = 0xde,	/* define hierarchical progression	*/
  M_EXP   = 0xdf,	/* expand reference image(s)		*/
  
  M_APP0  = 0xe0,	/* application marker, used for JFIF	*/
  M_APP1  = 0xe1,	/* application marker			*/
  M_APP2  = 0xe2,	/* application marker			*/
  M_APP3  = 0xe3,	/* application marker			*/
  M_APP4  = 0xe4,	/* application marker			*/
  M_APP5  = 0xe5,	/* application marker			*/
  M_APP6  = 0xe6,	/* application marker			*/
  M_APP7  = 0xe7,	/* application marker			*/
  M_APP8  = 0xe8,	/* application marker			*/
  M_APP9  = 0xe9,	/* application marker			*/
  M_APP10 = 0xea,	/* application marker			*/
  M_APP11 = 0xeb,	/* application marker			*/
  M_APP12 = 0xec,	/* application marker			*/
  M_APP13 = 0xed,	/* application marker			*/
  M_APP14 = 0xee,	/* application marker, used by Adobe	*/
  M_APP15 = 0xef,	/* application marker			*/
  
  M_JPG0  = 0xf0,	/* reserved for JPEG extensions		*/
  M_JPG13 = 0xfd,	/* reserved for JPEG extensions		*/
  M_COM   = 0xfe,	/* comment				*/
  
  M_TEM   = 0x01,	/* temporary use			*/

  M_ERROR = 0x100	/* dummy marker, internal use only	*/
} JPEG_MARKER;

/*
 * The following routine used to be a macro in its first incarnation:
 *  #define get_2bytes(fp) ((unsigned int) (getc(fp) << 8) + getc(fp))
 * However, this is bad programming since C doesn't guarantee
 * the evaluation order of the getc() calls! As suggested by
 * Murphy's law, there are indeed compilers which produce the wrong
 * order of the getc() calls, e.g. the Metrowerks C compilers for BeOS
 * and Macintosh.
 * Since there are only few calls we don't care about the performance 
 * penalty and use a simplistic C function.
 */

/* read two byte parameter, MSB first */
static unsigned int get_2bytes(FILE *fp)
{
    unsigned int val;
    val = getc(fp) << 8;
    val += getc(fp);
    return val;
}

static int next_marker P1(FILE *, fp)
{ /* look for next JPEG Marker  */
  int c, nbytes = 0;

  if (feof(fp))
    return M_ERROR;                 /* dummy marker               */

  do {
    do {                            /* skip to FF 		  */
      nbytes++;
      c = getc(fp);
    } while (c != 0xFF);
    do {                            /* skip repeated FFs  	  */
      c = getc(fp);
    } while (c == 0xFF);
  } while (c == 0);                 /* repeat if FF/00 	      	  */

  return c;
}

/* analyze JPEG marker */
BOOL AnalyzeJPEG P1(imagedata *, image)
{
  int b, c, unit;
  unsigned long i, length = 0;
#define APP_MAX 255
  unsigned char appstring[APP_MAX];
  BOOL SOF_done = FALSE;

  /* Tommy's special trick for Macintosh JPEGs: simply skip some  */
  /* hundred bytes at the beginning of the file!		  */
  do {
    do {                            /* skip if not FF 		  */
      c = getc(image->fp);
    } while (!feof(image->fp) && c != 0xFF);

    do {                            /* skip repeated FFs 	  */
      c = getc(image->fp);
    } while (c == 0xFF);

    /* remember start position */
    if ((image->startpos = ftell(image->fp)) < 0L) {
      fprintf(stderr, "Error: internal error in ftell()!\n");
      return FALSE;
    }
    image->startpos -= 2;           /* subtract marker length     */

    if (c == M_SOI) {
      fseek(image->fp, image->startpos, SEEK_SET);
      break;
    }
  } while (!feof(image->fp));

  if (feof(image->fp)) {
    fprintf(stderr, "Error: SOI marker not found!\n");
    return FALSE;
  }

  if (image->startpos > 0L && !quiet) {
    fprintf(stderr, "Note: skipped %ld bytes ", image->startpos);
    fprintf(stderr, "Probably Macintosh JPEG file?\n");
  }

  /* process JPEG markers */
  while (!SOF_done && (c = next_marker(image->fp)) != M_EOI) {
    switch (c) {
      case M_ERROR:
	fprintf(stderr, "Error: unexpected end of JPEG file!\n");
	return FALSE;

      /* The following are not officially supported in PostScript level 2 */
      case M_SOF2:
      case M_SOF3:
      case M_SOF5:
      case M_SOF6:
      case M_SOF7:
      case M_SOF9:
      case M_SOF10:
      case M_SOF11:
      case M_SOF13:
      case M_SOF14:
      case M_SOF15:
	fprintf(stderr, 
         "Warning: JPEG file uses compression method %X - proceeding anyway.\n",
		  c);
        fprintf(stderr, 
		"PostScript output does not work on all PS interpreters!\n");
	/* FALLTHROUGH */

      case M_SOF0:
      case M_SOF1:
	length = get_2bytes(image->fp);    /* read segment length  */

	image->bits_per_component = getc(image->fp);
	image->height             = get_2bytes(image->fp);
	image->width              = get_2bytes(image->fp);
	image->components         = getc(image->fp);

	SOF_done = TRUE;
	break;

      case M_APP0:		/* check for JFIF marker with resolution */
	length = get_2bytes(image->fp);

	for (i = 0; i < length-2; i++) {	/* get contents of marker */
	  b = getc(image->fp);
	  if (i < APP_MAX)			/* store marker in appstring */
	    appstring[i] = b;
	}

	/* Check for JFIF application marker and read density values
	 * per JFIF spec version 1.02.
	 * We only check X resolution, assuming X and Y resolution are equal.
	 * Use values only if resolution not preset by user or to be ignored.
	 */

#define ASPECT_RATIO	0	/* JFIF unit byte: aspect ratio only */
#define DOTS_PER_INCH	1	/* JFIF unit byte: dots per inch     */
#define DOTS_PER_CM	2	/* JFIF unit byte: dots per cm       */

	if (image->dpi == DPI_USE_FILE && length >= 14 &&
	    !strncmp((const char *)appstring, "JFIF", 4)) {
	  unit = appstring[7];		        /* resolution unit */
	  					/* resolution value */
	  image->dpi = (float) ((appstring[8]<<8) + appstring[9]);

	  if (image->dpi == 0.0) {
	    image->dpi = DPI_USE_FILE;
	    break;
	  }

	  switch (unit) {
	    /* tell the caller we didn't find a resolution value */
	    case ASPECT_RATIO:
	      image->dpi = DPI_USE_FILE;
	      break;

	    case DOTS_PER_INCH:
	      break;

	    case DOTS_PER_CM:
	      image->dpi *= (float) 2.54;
	      break;

	    default:				/* unknown ==> ignore */
	      fprintf(stderr, 
		"Warning: JPEG file contains unknown JFIF resolution unit - ignored!\n");
	      image->dpi = DPI_IGNORE;
	      break;
	  }
	}
        break;

      case M_APP14:				/* check for Adobe marker */
	length = get_2bytes(image->fp);

	for (i = 0; i < length-2; i++) {	/* get contents of marker */
	  b = getc(image->fp);
	  if (i < APP_MAX)			/* store marker in appstring */
	    appstring[i] = b;
	}

	/* Check for Adobe application marker. It is known (per Adobe's TN5116)
	 * to contain the string "Adobe" at the start of the APP14 marker.
	 */
	if (length >= 12 && !strncmp((const char *) appstring, "Adobe", 5))
	  image->adobe = TRUE;			/* set Adobe flag */

	break;

      case M_SOI:		/* ignore markers without parameters */
      case M_EOI:
      case M_TEM:
      case M_RST0:
      case M_RST1:
      case M_RST2:
      case M_RST3:
      case M_RST4:
      case M_RST5:
      case M_RST6:
      case M_RST7:
	break;

      default:			/* skip variable length markers */
	length = get_2bytes(image->fp);
	for (length -= 2; length > 0; length--)
	  (void) getc(image->fp);
	break;
    }
  }

  /* do some sanity checks with the parameters */
  if (image->height <= 0 || image->width <= 0 || image->components <= 0) {
    fprintf(stderr, "Error: DNL marker not supported in "
            "PostScript Level 2!\n");
    return FALSE;
  }

  /* some broken JPEG files have this but they print anyway... */
  if (length != (unsigned int) (image->components * 3 + 8))
    fprintf(stderr, "Warning: SOF marker has incorrect length - ignored!\n");

  if (image->bits_per_component != 8) {
    fprintf(stderr, "Error: %d bits per color component ",
	image->bits_per_component);
    fprintf(stderr, "not supported in PostScript level 2!\n");
    return FALSE;
  }

  if (image->components!=1 && image->components!=3 && image->components!=4)
  {
      fprintf(stderr, "Error: unknown color space (%d components)!\n",
              image->components);
      return FALSE;
  }

  return TRUE;
}


/* ------------------------- jpeg2ps.c ------------------------- */
/* convert JPEG files to compressed PostScript Level 2 EPS       */
/* (C) Thomas Merz 1994-99                                       */
/* ------------------------- jpeg2ps.c ------------------------- */
#define VERSION		"V1.8"
#define READMODE        "r"
#define WRITEMODE       "w"
/* write (some) PS files in binary mode */

extern BOOL	AnalyzeJPEG P1(imagedata *, image);
extern int	ASCII85Encode P2(FILE *, in, FILE *, out);
extern void     ASCIIHexEncode P2(FILE *, in, FILE *, out);
extern char *optarg;
extern int optind, opterr, optopt;

#define BUFFERSIZE 1024
static char buffer[BUFFERSIZE];
static char *ColorSpaceNames[] = {"", "Gray", "", "RGB", "CMYK" };

/* Array of known page sizes including name, width, and height */
typedef struct { const char *name; int width; int height; } PageSize_s;

PageSize_s PageSizes[] = {
    {"a0",	2380, 3368},
    {"a1",	1684, 2380},
    {"a2",	1190, 1684},
    {"a3",	842, 1190},
    {"a4",	595, 842},
    {"a5",	421, 595},
    {"a6",	297, 421},
    {"b5",	501, 709},
    {"letter",	612, 792},
    {"legal",	612, 1008},
    {"ledger",	1224, 792},
    {"p11x17",	792, 1224}
};

#define PAGESIZELIST	(sizeof(PageSizes)/sizeof(PageSizes[0]))

#ifdef A4
int PageWidth  = 595;           /* page width A4 */
int PageHeight = 842;           /* page height A4 */
#else
int PageWidth  = 612;           /* page width letter */
int PageHeight = 792;           /* page height letter */
#endif

static void JPEGtoPS P2(imagedata *, JPEG, FILE *, PSfile)
{
  int llx, lly, urx, ury;        /* Bounding box coordinates */
  size_t n;
  float scale, sx, sy;           /* scale factors            */
  time_t t;
  int i;
  /* read image parameters and fill JPEG struct*/
  if (!AnalyzeJPEG(JPEG)) {
    fprintf(stderr, "Error: '%s' is not a proper JPEG file!\n",
            JPEG->filename);
    return;
  }

  if (!quiet)
      fprintf(stderr, "Note on file '%s': %dx%d pixel, %d color component%s\n",
	JPEG->filename, JPEG->width, JPEG->height, JPEG->components,
	(JPEG->components == 1 ? "" : "s"));

  /* "Use resolution from file" was requested, but we couldn't find any */
  if (JPEG->dpi == DPI_USE_FILE && !quiet)
  { 
      fprintf(stderr, "Note: no resolution values found in JPEG file - "
              "using standard scaling.\n");
      JPEG->dpi = DPI_IGNORE;
  }

  if (JPEG->dpi == DPI_IGNORE) {
    if (JPEG->width > JPEG->height && autorotate)
    {	/* switch to landscape if needed */
      JPEG->landscape = TRUE;
      if (!quiet)
	  fprintf(stderr, "Note: image width exceeds height - "
                  "producing landscape output!\n");
    }
    if (!JPEG->landscape) {       /* calculate scaling factors */
      sx = (float) (PageWidth - 2*Margin) / JPEG->width;
      sy = (float) (PageHeight - 2*Margin) / JPEG->height;
    }else {
      sx = (float) (PageHeight - 2*Margin) / JPEG->width;
      sy = (float) (PageWidth - 2*Margin) / JPEG->height;
    }
    scale = min(sx, sy);	/* We use at least one edge of the page */
  }
  else
  {
      if (!quiet)
          fprintf(stderr, "Note: Using resolution %d dpi.\n",
                  (int) JPEG->dpi);
      scale = 72 / JPEG->dpi;     /* use given image resolution */
  }

  if (JPEG->landscape) {
    /* landscape: move to (urx, lly) */
    urx = PageWidth - Margin;
    lly = Margin;
    ury = (int) (Margin + scale*JPEG->width + 0.9);    /* ceiling */
    llx = (int) (urx - scale * JPEG->height);          /* floor  */
  }else {
    /* portrait: move to (llx, lly) */
    llx = lly = Margin;
    urx = (int) (llx + scale * JPEG->width + 0.9);     /* ceiling */
    ury = (int) (lly + scale * JPEG->height + 0.9);    /* ceiling */
  }

  time(&t);

  /* produce EPS header comments */
  fprintf(PSfile, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(PSfile, "%%%%Creator: jpeg2ps %s by Thomas Merz\n", VERSION);
  fprintf(PSfile, "%%%%Title: %s\n", JPEG->filename);
  fprintf(PSfile, "%%%%CreationDate: %s", ctime(&t));
  fprintf(PSfile, "%%%%BoundingBox: %d %d %d %d\n", 
                   llx, lly, urx, ury);
  fprintf(PSfile, "%%%%DocumentData: %s\n", 
                  JPEG->mode == BINARY ? "Binary" : "Clean7Bit");
  fprintf(PSfile, "%%%%LanguageLevel: 2\n");
  fprintf(PSfile, "%%%%EndComments\n");
  fprintf(PSfile, "%%%%BeginProlog\n");
  fprintf(PSfile, "%%%%EndProlog\n");
  fprintf(PSfile, "%%%%Page: 1 1\n");

  fprintf(PSfile, "/languagelevel where {pop languagelevel 2 lt}");
  fprintf(PSfile, "{true} ifelse {\n");
  fprintf(PSfile, "  (JPEG file '%s' needs PostScript Level 2!",
                  JPEG->filename);
  fprintf(PSfile, "\\n) dup print flush\n");
  fprintf(PSfile, "  /Helvetica findfont 20 scalefont setfont ");
  fprintf(PSfile, "100 100 moveto show showpage stop\n");
  fprintf(PSfile, "} if\n");

  fprintf(PSfile, "save\n");
  fprintf(PSfile, "/RawData currentfile ");

  if (JPEG->mode == ASCIIHEX)            /* hex representation... */
    fprintf(PSfile, "/ASCIIHexDecode filter ");
  else if (JPEG->mode == ASCII85)        /* ...or ASCII85         */
    fprintf(PSfile, "/ASCII85Decode filter ");
  /* else binary mode: don't use any additional filter! */

  fprintf(PSfile, "def\n");

  fprintf(PSfile, "/Data RawData << ");
  fprintf(PSfile, ">> /DCTDecode filter def\n");

  /* translate to lower left corner of image */
  fprintf(PSfile, "%d %d translate\n", (JPEG->landscape ? 
                   PageWidth - Margin : Margin), Margin);

  if (JPEG->landscape)                 /* rotation for landscape */
    fprintf(PSfile, "90 rotate\n");
      
  fprintf(PSfile, "%.2f %.2f scale\n", /* scaling */
                   JPEG->width * scale, JPEG->height * scale);
  fprintf(PSfile, "/Device%s setcolorspace\n", 
                  ColorSpaceNames[JPEG->components]);
  fprintf(PSfile, "{ << /ImageType 1\n");
  fprintf(PSfile, "     /Width %d\n", JPEG->width);
  fprintf(PSfile, "     /Height %d\n", JPEG->height);
  fprintf(PSfile, "     /ImageMatrix [ %d 0 0 %d 0 %d ]\n",
                  JPEG->width, -JPEG->height, JPEG->height);
  fprintf(PSfile, "     /DataSource Data\n");
  fprintf(PSfile, "     /BitsPerComponent %d\n", 
                  JPEG->bits_per_component);

  /* workaround for color-inverted CMYK files produced by Adobe Photoshop:
   * compensate for the color inversion in the PostScript code
   */
  if (JPEG->adobe && JPEG->components == 4) {
    if (!quiet)
	fprintf(stderr, "Note: Adobe-conforming CMYK file - "
                "applying workaround for color inversion.\n");
    fprintf(PSfile, "     /Decode [1 0 1 0 1 0 1 0]\n");
  }
  else
  {
      fprintf(PSfile, "     /Decode [0 1");
      for (i = 1; i < JPEG->components; i++) 
          fprintf(PSfile," 0 1");
      fprintf(PSfile, "]\n");
  }

  fprintf(PSfile, "  >> image\n");
  fprintf(PSfile, "  Data closefile\n");
  fprintf(PSfile, "  RawData flushfile\n");
  fprintf(PSfile, "  showpage\n");
  fprintf(PSfile, "  restore\n");
  fprintf(PSfile, "} exec");

  /* seek to start position of JPEG data */
  fseek(JPEG->fp, JPEG->startpos, SEEK_SET);

  switch (JPEG->mode) {
  case BINARY:
    /* important: ONE blank and NO newline */
    fprintf(PSfile, " ");
    /* copy data without change */
    while ((n = fread(buffer, 1, sizeof(buffer), JPEG->fp)) != 0)
      fwrite(buffer, 1, n, PSfile);
    break;

  case ASCII85:
    fprintf(PSfile, "\n");

    /* ASCII85 representation of image data */
    if (ASCII85Encode(JPEG->fp, PSfile)) {
      fprintf(stderr, "Error: internal problems with ASCII85Encode!\n");
      exit(1);
    }
    break;

  case ASCIIHEX:
    /* hex representation of image data (useful for buggy dvips) */
    ASCIIHexEncode(JPEG->fp, PSfile);
    break;
  }
  fprintf(PSfile, "\n%%%%EOF\n");
}

static void usage P0(void)
{
    fprintf(stderr, "jpeg2ps %s: convert JPEG files to PostScript Level 2.\n",
            VERSION);
    fprintf(stderr, "(C) Thomas Merz 1994-1999\n\n");
    fprintf(stderr,
            "usage: jpeg2ps [options] jpegfile > epsfile\n");
    fprintf(stderr, "-a        auto rotate: "
            "produce landscape output if width > height\n");
    fprintf(stderr, "-b        binary mode: "
            "output 8 bit data (default: 7 bit with ASCII85)\n");
    fprintf(stderr, "-h        hex mode: "
            "output 7 bit data in ASCIIHex encoding\n");
    fprintf(stderr, "-o <name> output file name\n");
    fprintf(stderr, "-p <size> page size name. Known names are:\n");
    fprintf(stderr, "          "
            "a0, a1, a2, a3, a4, a5, a6, b5, letter, legal, ledger, p11x17\n");
    fprintf(stderr, "-q        "
            "quiet mode: suppress all informational messages\n");
    fprintf(stderr, "-r <dpi>  resolution value (dots per inch)\n");
    fprintf(stderr, "          "
            "0 means use value given in file, if any (disables autorotate)\n");
    exit(1);
}

imagedata image;
FILE *outfile;
int opt, pagesizeindex = -1;

#ifndef _OLD_MAIN_TEST
int jpg2eps(char *filename_jpg, char *filename_eps, char *options)
{
    int argc;
    char *buffer, **argv;
    if (options==NULL) options="";
    /* count the number of text segments in options */
    for (argc=1,buffer=options;;)
    {
        while (isspace(*buffer)) buffer++;
        if (*buffer!=0) argc++;
        while ((*buffer!=0)&&(!isspace(*buffer))) buffer++;
        if (*buffer==0) break;
    }
    argc += 3;
    argv = (char **)malloc(argc*sizeof(char *));
    buffer = (char *)malloc(strlen(options)+1);
    strcpy(buffer,options);
    for (argc=1,argv[0]=buffer;;)
    {
        while (isspace(*argv[0])) argv[0]++;
        if (*argv[0]!=0) argv[argc++] = argv[0];
        while ((*argv[0]!=0)&&(!isspace(*argv[0]))) argv[0]++;
        if (*argv[0]==0) break; else *(argv[0]++)=0;
    }
    argv[0] = "jpg2eps";
    argv[argc++] = "-o";
    argv[argc++] = filename_eps;
    argv[argc++] = filename_jpg;
    /* filename_jpg should always be the last argument */
    optind = 1; opterr = 1; /* getoptreset(); */
#else
    int main P2(int, argc, char **, argv) {
#endif
    image.filename = NULL;
    image.mode     = ASCII85;
    image.startpos = 0L;
    image.landscape= FALSE;
    image.dpi      = DPI_IGNORE;
    image.adobe    = FALSE;
    outfile = stdout;
    if (argc == 1) usage();
    while ((opt = getopt(argc, argv, "abho:p:qr:")) != -1)
        switch (opt) {
            case 'a':
                autorotate = TRUE;
                break;
            case 'b':
                image.mode = BINARY;
                break;
            case 'h':
                image.mode = ASCIIHEX;
                break;
            case 'o':
                outfile = fopen(optarg, "w");
                if (outfile == NULL)
                {
                    fprintf(stderr, "Error: cannot open output file %s.\n",
                            optarg);
                    exit(-2);
                } 
                break;
            case 'p':
                for(pagesizeindex=0; pagesizeindex < PAGESIZELIST;
                    pagesizeindex++)
                    if (!strcmp((const char *) optarg,
                                PageSizes[pagesizeindex].name))
                    {
                        PageHeight = PageSizes[pagesizeindex].height;
                        PageWidth = PageSizes[pagesizeindex].width;
                        break;
                    }
                if (pagesizeindex == PAGESIZELIST)
                {	/* page size name not found */
                    fprintf(stderr, "Error: Unknown page size %s.\n", optarg);
                    exit(-3);
                }
                break;
            case 'q':
                quiet = (!quiet);
                break;
            case 'r':
                image.dpi = (float) atof(optarg);
                if (image.dpi < 0)
                {
                    fprintf(stderr, "Error: bad resolution value %f !\n",
                            image.dpi);
                    exit(1);
                }
                break;
            case '?':
                usage();
        }
    if (pagesizeindex != -1 && ! quiet)	/* page size user option given */
        fprintf(stderr, "Note: Using %s page size.\n",
                PageSizes[pagesizeindex].name);
    if (optind == argc)	/* filename missing */
        usage();
    else
        image.filename = argv[optind];
    if (!image.filename) usage();
    if ((image.fp = fopen(image.filename, READMODE)) == NULL)
    {
        fprintf(stderr, "Error: couldn't read JPEG file '%s'!\n", 
                image.filename),
            exit(1);
    }
    JPEGtoPS(&image, outfile);      /* convert JPEG data */
    argc = ftell(outfile);
    fclose(image.fp);
    fclose(outfile);
#ifndef _OLD_MAIN_TEST
    free(argv);
    free(buffer);
    return(argc);
#else
    return 0;
#endif 
}

#ifdef _TEST
int main (int argc, char *argv[])
{
    jpg2eps("/tmp/cylinders.jpg", "/tmp/a.eps", NULL);
    return (0);
}
#endif /* _TEST */
/* cc -D_TEST jpg2eps.c -o jpg2eps; jpg2eps */
