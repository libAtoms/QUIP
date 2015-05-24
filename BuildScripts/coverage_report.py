#!/usr/bin/env python

# parse output of gcov and make a nicely formatted coverage report

import sys, os, operator

if len(sys.argv[1:]) == 0:
    print 'Usage: coverage_report.py GCDA_FILES'
    sys.exit(1)

lines = os.popen('gcov %s' % ' '.join(sys.argv[1:])).readlines()
if len(lines) % 4 != 0:
    print 'Number of lines should be multiple of four'
    sys.exit(1)

coverage = {}

quip_root = 'QUIP_ROOT' in os.environ and os.environ['QUIP_ROOT']

for file_line, coverage_line in zip(lines[::4], lines[1::4]):
    filename = file_line.split()[1].replace("'","")
    if filename.startswith('/usr'): continue
    if quip_root and filename.startswith(quip_root):
        filename = filename[len(quip_root)+1:]
    coverage[filename] = float(coverage_line[coverage_line.index(':')+1:coverage_line.index('%')])

def print_coverage(coverage):
    print '%-50s %-10s' % ('Filename', 'Coverage')
    print '-'*50 + ' '+'-'*10
    for filename, cov in sorted(coverage.iteritems(), key=operator.itemgetter(1), reverse=True):
        print '%-50s %10.1f' % (filename, cov)

modules = set(sorted([k.split('/')[0] for k in coverage.keys()]))

for module in modules:
    print '='*61
    print module
    print '='*61
    module_coverage = dict([(k[len(module)+1:],v) for (k,v) in coverage.iteritems() if k.startswith(module)])
    print_coverage(module_coverage)
    print

