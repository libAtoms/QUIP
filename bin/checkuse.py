#! /usr/bin/env python

"""
Check order of use statements
"""

from __future__ import print_function, unicode_literals

import sys


def parse_use_statements(filename):
    f = open(filename, 'r')
    line = f.readline()
    start_parse = False
    stop_parse = False
    in_interface = False
    module_name = None
    use_statements = [ ]
    while line and not stop_parse:
        split_line = [x.strip().lower() for x in line.split()]
        if len(split_line) > 0:
            if split_line[0] == 'module' and not start_parse:
                start_parse = True
                module_name = split_line[1]
            elif start_parse:
                if split_line[0] == 'interface':
                    in_interface = True
                elif split_line[0] == 'endinterface' or (split_line[0] == 'end' and split_line[1] == 'interface'):
                    in_interface = False
                if split_line[0] == 'use' and not in_interface:
                    use_statements += [split_line[1]]
                elif split_line[0] == 'contains':
                    stop_parse = True
        line = f.readline()
    f.close()
    return module_name, use_statements


def check_deps(deplist, name, used=None):
    if used is None:
        used = set()

    for mod in deplist[name]:
        if mod in used:
            print('In module {0}: Module {1} should be included earlier.'
                  ''.format(name, mod))
        if mod in deplist:
            for depmod in deplist[mod]:
                used.add(depmod)


deplist = {}
for fn in sys.argv[1:]:
    name, use = parse_use_statements(fn)
    if name is not None:
        deplist[name] = use

for mod in deplist.keys():
    check_deps(deplist, mod)
