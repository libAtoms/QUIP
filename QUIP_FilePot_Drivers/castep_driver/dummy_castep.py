#!/usr/bin/env python
# HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HND X
# HND X   libAtoms+QUIP: atomistic simulation library
# HND X
# HND X   Portions of this code were written by
# HND X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
# HND X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
# HND X
# HND X   Copyright 2006-2010.
# HND X
# HND X   Not for distribution
# HND X
# HND X   Portions of this code were written by Noam Bernstein as part of
# HND X   his employment for the U.S. Government, and are not subject
# HND X   to copyright in the USA.
# HND X
# HND X   When using this software, please cite the following reference:
# HND X
# HND X   http://www.libatoms.org
# HND X
# HND X  Additional contributions by
# HND X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
# HND X
# HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Dummy castep that returns successive segments of castep log file
# Useful for debugging a failed supercomputer run. 

import sys, os, shutil

if len(sys.argv) != 2:
    print 'Usage: %s stem\n'
    sys.exit(1)

stem = sys.argv[1]

if not os.path.exists(stem+'.castep_log.dummy_castep'):
    if not os.path.exists(stem+'.castep_log'):
        print 'Log file %s.castep_log not found.' % stem
        sys.exit(1)
    else:
        shutil.copyfile(stem+'.castep_log', stem+'.castep_log.dummy_castep')

try:
    frame = int(open(stem+'.dummy_castep','r').read())
except IOError:
    frame = 0


lines = open(stem+'.castep_log.dummy_castep').readlines()

   
offsets = [i-2  for (i,line) in  enumerate(lines) if line == ' |      CCC   AA    SSS  TTTTT  EEEEE  PPPP        |\n']
n_frames = len(offsets)

if frame >= n_frames:
    if os.path.exists(stem+'.castep'): os.unlink(stem+'.castep')
    print 'dummy_castep: frame %d out of range' % frame
    sys.exit(1)

print 'dummy_castep: returning frame %d/%d of %s.castep_log.dummy_castep' % (frame, n_frames, stem)

offsets.append(len(lines)-1)

#for frame in range(n_frames):
#    print frame,offsets[frame],offsets[frame+1]

output = open(stem+'.castep','w')
output.writelines(lines[offsets[frame]:offsets[frame+1]])
output.close()

open(stem+'.dummy_castep','w').write(str(frame+1))

