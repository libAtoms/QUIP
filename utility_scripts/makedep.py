#!/usr/bin/env python

import sys, os, cPickle, f90doc
from distutils.dep_util import newer

def warn(s):
    sys.stderr.write(s+'\n')

argv = sys.argv[1:]
args = []
while len(argv) > 0 and argv[0] != '--':
    args.append(argv.pop(0))

if argv ==[] or argv[0] != '--':
    print 'Usage: makedep.py [ -lc | -uc ] [ -suffix suffix ] -- filename [ filename2 ... ]'
    sys.exit(1)
else:
    argv.pop(0)

mod_define = {}
mod_uses = {}
mod_includes = {}
filename_from_mod = {}

f90doc.do_debug = True

for src in argv:
    base, ext = os.path.splitext(os.path.basename(src))
    if ext not in ['.f', '.F', '.f90', '.F90', '.f95', '.F95']:
        warn('warning: %s is not a fortran file' % src)
        continue
    if not os.path.exists(src):
        warn("warning: file %s doesn't exist" % src)
        continue
    
    fpp = base + '.fpp'
    if not os.path.exists(fpp):
        warn("preprocessed sources not found for file %s, skipping" % src)
        continue
    
    f90doc_file = base + '.f90doc'

    if os.path.exists(f90doc_file):
        (programs, modules, functs, subts) = cPickle.load(open(f90doc_file))
    else:
        programs = None
        modules = None
        functs = None
        subrts = None

    if newer(src, f90doc_file):
        new_programs, new_modules, new_functs, new_subts = f90doc.read_files([fpp])

        if (new_programs != programs or new_modules != modules or
            new_functs != functs or new_subts != subts):
            programs = new_programs
            modules  = new_modules
            functs   = new_functs
            subts    = new_subts
            cPickle.dump((programs, modules, functs, subts), open(f90doc_file, 'w'))

    mod_define[src] = []
    mod_uses[src] = []
    for mod,filename in modules:
        name = mod.name.lower()
        mod_define[src].append(name)
        filename_from_mod[name] = base+ext
        mod_uses[src].extend([u.lower() for u in mod.uses])

    mod_includes[src] = []
    include_files = [line.split()[2].replace('"','') for line in open(fpp).readlines() if line.startswith('# ') and len(line.split()) > 2 ]
    for inc in include_files:
        if os.path.isfile(inc) and os.path.basename(inc) != base+ext and inc not in mod_includes[src]:
            mod_includes[src].append(inc)

for src in argv:
    base, ext = os.path.splitext(os.path.basename(src))
    if src not in mod_uses:
        continue

    fobj = base + '.o'
    f90doc_file = base + '.f90doc'
    f90doc_uses = [ os.path.splitext(filename_from_mod[mod])[0]+'.f90doc' for mod in mod_uses[src] if mod in filename_from_mod ]

    print "%s: %s\n" % (f90doc_file, base+ext)
    print "%s: \\\n  %s %s %s\n" % (fobj, ' '.join(f90doc_uses), ' '.join(mod_includes[src]), base+ext)

