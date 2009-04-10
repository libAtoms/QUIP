#!/usr/bin/env python

import os, sys, re, getopt
from difflib import ndiff

def loadtests(testfiles):
   tests = []
   for testfile in testfiles:
      data = open(testfile).read()

      endl = data.index('\n---')
      command = data[:endl]
      data = data[endl+4:]
      
      fields = data.split('---')

      if len(fields) % 2 != 0:
         raise ValueError('Badly formed testfile %s' % testfile)

      infiles = {}
      outfiles = {}
      for name, contents in zip(fields[::2], fields[1::2]):
         if name[0] == '<':
            infiles[name[1:].strip()]  = contents[1:].splitlines(1)
         elif name[0] == '>':
            outfiles[name[1:].strip()] = contents[1:].splitlines(1)
      
      if not 'stdout' in outfiles:
         raise ValueError('Missing stdout section in testfile %s' % testfile)
      
      tests.append((testfile, command, infiles, outfiles))

   return tests


def remove_ignored_stdout_lines(line):
   return not (line.startswith('System::Hello World') or
               line.startswith('System::Finalise') or
	       line.startswith('TIMER'))


def runtest(command, diff_method, infiles, outfiles):

   for name, contents in infiles.iteritems():
      if name != 'stdin':
         f = open(name,'w')
         f.writelines(contents)
         f.close()

   # escape "s from shell
   command = command.replace('"',r'\"')

   # Remove old output files
   for name, contents in outfiles.iteritems():
      if name != 'stdout':
         if os.path.exists(name): os.unlink(name)

   stdin, stdout, stderr = os.popen3('./build.%s/%s' % (ARCH, command))

   if 'stdin' in infiles:
      stdin.writelines(infiles['stdin'])
      stdin.close()

   cmpout = {}
   cmpout['stdout'] = (filter(remove_ignored_stdout_lines,stdout.readlines()),
                       filter(remove_ignored_stdout_lines,outfiles['stdout']))

   sys.stdout.writelines(stderr.readlines())

   for name, contents in outfiles.iteritems():
      if name != 'stdout':
         if not os.path.exists(name):
            print 'Output file %s not produced' % name
            return False
         
         f = open(name,'r')
         exp_contents = f.readlines()
         cmpout[name] = (contents, exp_contents)

   diffs = {}
   reduced_diffs = {}
   for name, (contents, exp_contents) in cmpout.iteritems():
      diffs[name] = list(do_diff(contents, exp_contents, diff_method))
   
   no_differences = True
   for name,rd in diffs.iteritems():
      if rd != []:
         no_differences = False
         print 'Mismatch in file %s' % name
         print ''.join(rd)
         #print ''.join(cmpout[name][0])
         #print ''.join(cmpout[name][1])

   # remove tmp files
   for fname in infiles.keys()+outfiles.keys():
      if not fname in ('stdin','stdout'):
         if os.path.exists(fname): os.unlink(fname)

   return no_differences

# Use ndiff(1) if found, otherwise use built in diff
def do_diff(a, b, diff_method):
   if diff_method != 'built in':
      fnamea, fnameb = os.tmpnam(), os.tmpnam()
      fa, fb  = open(fnamea, 'w'), open(fnameb, 'w')
      fa.writelines(a)
      fb.writelines(b)
      fa.close()
      fb.close()
      cin, cout, cerr = os.popen3("%s -abserr 1e-10 -separators '[ \t=\(\)\n]' %s %s" % (diff_method, fnamea, fnameb))
      result = cout.readlines()
      result.extend(cerr.readlines())
      result = filter(lambda s: not s.startswith('###'), result)
      os.remove(fnamea)
      os.remove(fnameb)
      return result
   else:
      return ndiff(a, b)


def runtests(tests):
   # Look for ndiff(1) executable somewhere on PATH
   diff_method = 'built in'
   for p in os.environ['PATH'].split(':'):
      candidate_path = os.path.join(p,'ndiff')
      if os.path.exists(candidate_path):
         diff_method = candidate_path
         break

   if diff_method != 'built in':
      # Disable RuntimeWarning printed by os.tmpnam
      def nowarning(message, category, filename, lineno): pass
      import warnings
      warnings.showwarning = nowarning

   print 'Running %d regression tests using diff method: %s' % (len(tests), diff_method)
   print
   n_fail = 0
   for name, command, infiles, outfiles in tests:

      print '  Running test : %s  ' % name, 
      if runtest(command, diff_method, infiles, outfiles):
         print 'OK'
      else:
         print 'FAIL'
         n_fail += 1

   print
   print 'Failed %d tests.' % n_fail
   return n_fail == 0


def mktest(name, command, read_stdin=False, infiles=[], outfiles=[]):

   testfile = open(name, 'w')
   testfile.write(command+'\n')

   if infiles != []:
      for fname in infiles:
         testfile.write('---<%s---\n' % fname)
         lines = open(fname,'r').readlines()
         testfile.writelines(lines)

   print 'running command %s' % command

   # escape "s from shell
   command = command.replace('"',r'\"')

   stdin, stdout = os.popen2('./build.%s/%s' % (ARCH, command))

   if read_stdin:
      testfile.write('---<stdin---\n')
      inputlines = sys.stdin.readlines()
      testfile.writelines(inputlines)
      stdin.writelines(inputlines)
      stdin.close()

   testfile.write('--->stdout---\n')
   testfile.writelines(stdout.readlines())

   if outfiles != []:
      for fname in outfiles:
         testfile.write('--->%s---\n' % fname)
         lines = open(fname, 'r').readlines()
         testfile.writelines(lines)

   testfile.close()
         

def print_usage():
   print 'Regression testing program'
   print 'James Kermode <jrk33@cam.ac.uk>'
   print 
   print 'Usage:  '
   print '  To run tests:       %s ARCH -r TESTFILE...' % sys.argv[0]
   print '  To make a new test: %s ARCH -m TESTFILE [-i INFILE]...'  % sys.argv[0]
   print '                         [-o OUTFILE]... COMMAND'
   print
   print 'where ARCH is the Makefile arch suffix, TESTFILE is a file containing a'
   print 'test case, COMMAND is the test case command, and INFILE and OUTFILE are'
   print 'input and output files for the new test case. If "-i stdin" is present'
   print 'then stdin will be read and passed along to the test program.'

      
if __name__ == '__main__':

   if len(sys.argv) < 4:
      print_usage()
      sys.exit(1)

   ARCH = sys.argv[1]

   opts, args = getopt.getopt(sys.argv[2:],'rm:i:o:')

   print opts, args

   if opts[0][0] == '-r':
      try:
         tests = loadtests(args)
      except ValueError,v:
         print 'Error loading testfile: %s' % v
         sys.exit(1)
      runtests(tests)
      
   elif opts[0][0] == '-m':

      testfile = opts[0][1]

      infiles = []
      outfiles = []

      for opt, value in opts:
         if opt == '-i': infiles.append(value)
         if opt == '-o': outfiles.append(value)

      command = ' '.join(args)

      print infiles, outfiles, command

      read_stdin = False
      if 'stdin' in infiles:
         read_stdin = True
         del infiles[infiles.index('stdin')]

      mktest(testfile, command, read_stdin, infiles, outfiles)
                         
      
   else:
      print_usage()
      sys.exit(1)
