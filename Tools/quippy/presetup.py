#!/usr/bin/env python

import os
from numpy.distutils.ccompiler import new_compiler

# Compile simple C program to find sizeof(void *) on this arch
print 'Calculating sizeof(void *)...', 
if not os.path.exists('./sizeof_void_ptr'):
    cc = new_compiler()
    cc.link_executable(cc.compile(['sizeof_void_ptr.c']),'sizeof_void_ptr')
sizeof_void_ptr = int(os.popen('./sizeof_void_ptr').read())
print sizeof_void_ptr, 'bytes.'

# Generate .f2py_f2cmap
print 'Generating .f2py_f2cmap...',
size_t_lookup = {4:'int', 8:'long_long'}
if sizeof_void_ptr not in size_t_lookup:
    raise ValueError("Can't guess C type for size_t with sizeof(void *) = %d" % sizeof_void_ptr)
cmapf = open('.f2py_f2cmap','w')
cmapf.write("""{
 'real':{'dp':'double'},
 'complex':{'dp':'complex_double'},
 'integer':{'C_SIZE_T':'%s'}
}
""" % size_t_lookup[sizeof_void_ptr])
cmapf.close()
print 'done.'
