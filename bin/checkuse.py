#! /usr/bin/env python

#
# Check order of use statements
#

import sys
import string

###

def parse_use_statements(fn):
    f = open(fn, 'r')
    l = f.readline()
    start_parse = False
    stop_parse = False
    in_interface = False
    module_name = None
    use_statements = [ ]
    while l and not stop_parse:
        s = map(string.lower, map(string.strip, l.split()))
        if len(s) > 0:
            if s[0] == 'module' and not start_parse:
                start_parse = True
                module_name = s[1]
            elif start_parse:
                if s[0] == 'interface':
                    in_interface = True
                elif s[0] == 'endinterface' or ( s[0] == 'end' and s[1] == 'interface' ):
                    in_interface = False
                if s[0] == 'use' and not in_interface:
                    use_statements += [ s[1] ]
                elif s[0] == 'contains':
                    stop_parse = True
        l = f.readline()
    f.close()
    return module_name, use_statements

###

def check_deps(deplist, name, used=None):
    if used is None:
        used = set()

    for mod in deplist[name]:
        if mod in used:
            print 'In module %s: Module %s should be included earlier.' % ( name, mod )
        if mod in deplist:
            for depmod in deplist[mod]:
                used.add(depmod)


###

deplist = { }
for fn in sys.argv[1:]:
    name, use = parse_use_statements(fn)
    if name is not None:
        deplist[name] = use

for mod in deplist.keys():
    check_deps(deplist, mod)
