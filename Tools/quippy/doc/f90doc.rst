:mod:`f90doc` -- Fortran 90 source code scanner
===============================================

.. module:: f90doc
   :platform: Unix
   :synopsis: Read Fortan 90 sources and extract documenataion strings

.. moduleauthor:: Ian Rutt <i.c.rutt@swansea.ac.uk>
.. moduleauthor:: James Kermode <james.kermode@kcl.ac.uk>

:mod:`f90doc` can be used either as a command line tool or as a module.

Command line usage
------------------

When used as a command line program it reads from the Fortran source
files given as arguments and writes LaTeX source to :file:`stdout`. The
program accepts the following options:

.. program:: f90doc

.. cmdoption:: -t <title>

   specify a title.


.. cmdoption::  -a <author> 

   specify an author.

.. cmdoption::  -i <intro> 

   specify introductory LaTeX file to be included.

.. cmdoption::  -n 

   don't use an introductory LaTeX file.

.. cmdoption::  -s 

   use short documentation format.

.. cmdoption::  -l 

   write latex output to :file:`stdout`.


Module usage
------------

The only function that is designed to be used when :mod:`f90doc` is
imported as a module is :func:`read_files`.  The classes
:class:`f90doc.C_prog`, :class:`f90doc.C_module`,
:class:`f90doc.C_subt`, :class:`f90doc.C_funct`,
:class:`f90doc.C_decl` and :class:`f90doc.C_interface` contain the
implementation of LaTeX output code, with each class represnting
structures in the Fortran code.

.. function:: read_files(in_files)
   
   Read Fortran 90 sources from the list of filenames `in_files`, and return
   a tuple `(programs, modules, functs, subts)`. 

   `programs` is a list of pairs of instances of :class:`f90doc.C_prog` and program names,
   `modules` is a list of pairs of instances of :class:`f90doc.C_module` and module names,
   `functs` is a list of pairs of instances of :class:`f90doc.C_funct` and function names
   `subts` is a list of pairs of instances of :class:`f90doc.C_subt` and subroutine names

   For example, to print the names of all modules defined in a list of filenames `in_files`::

       programs, modules, functs, subts = f90doc.read_files(in_files)
       for mod, name in modules:
          print name


