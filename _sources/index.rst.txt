.. HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
.. HQ X
.. HQ X   quippy: Python interface to QUIP atomistic simulation library
.. HQ X
.. HQ X   Copyright James Kermode 2010
.. HQ X
.. HQ X   These portions of the source code are released under the GNU General
.. HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
.. HQ X
.. HQ X   If you would like to license the source code under different terms,
.. HQ X   please contact James Kermode, james.kermode@gmail.com
.. HQ X
.. HQ X   When using this software, please cite the following reference:
.. HQ X
.. HQ X   http://www.jrkermode.co.uk/quippy
.. HQ X
.. HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

.. quippy documentation master file, created by
   sphinx-quickstart on Wed Sep  2 14:17:01 2009.

QUIP and quippy documentation
=============================

.. module:: quippy

The ``QUIP`` package (`GitHub <https://github.com/libAtoms/QUIP>`_) is a
collection of software tools to carry out molecular dynamics
simulations. It implements a variety of interatomic potentials and
tight binding quantum mechanics, and is also able to call external
packages, and serve as plugins to other software such as `LAMMPS
<http://lammps.sandia.gov>`_, `CP2K <http://www.cp2k.org>`_ and also
the python framework `ASE <https://wiki.fysik.dtu.dk/ase>`_.  Various
hybrid combinations are also supported in the style of QM/MM, with a
particular focus on materials systems such as metals and
semiconductors.

``quippy`` is a Python interface to ``QUIP`` that provides deep access to
most of the Fortran types and routines. The quippy interface is principally
maintained by `James Kermode <http://www.warwick.ac.uk/jrkermode>`_.

Long term support of the package is ensured by:
 - Noam Bernstein (Naval Research Laboratory)
 - Gabor Csanyi (University of Cambridge)
 - James Kermode (University of Warwick)

Portions of this code were written by: Albert Bartok-Partay, Livia
Bartok-Partay, Federico Bianchini, Anke Butenuth, Marco Caccin,
Silvia Cereda, Gabor Csanyi, Alessio Comisso, Tom Daff, ST John,
Chiara Gattinoni, Gianpietro Moras, James Kermode, Letif Mones,
Alan Nichol, David Packwood, Lars Pastewka, Giovanni Peralta, Ivan
Solt, Oliver Strickson, Wojciech Szlachta, Csilla Varnai, Steven
Winfield.

Copyright 2006-2020

Most of the publicly available version is released under the GNU
General Public license, version 2, with some portions in the public
domain.

Features
========

The following interatomic potentials are presently coded or linked in QUIP:

* EAM (fcc metals)
* Fanourgakis-Xantheas (water)
* Finnis-Sinclair (bcc metals)
* Flikkema-Bromley
* GAP (Gaussian Approximation Potentials: general many-body)
* Guggenheim-McGlashan
* Brenner (carbon)
* OpenKIM (general interface)
* Lennard-Jones
* Morse
* Partridge-Schwenke (water monomer)
* Stillinger-Weber (carbon, silicon, germanium)
* SiMEAM (silicon)
* Sutton-Chen
* Tangney-Scandolo (silica, titania etc)
* Tersoff (silicon, carbon)

The following tight-binding functional forms and parametrisations are implemented:

* Bowler
* DFTB
* GSP
* NRL-TB

Contents
========

.. toctree::
   :maxdepth: 2

   install.rst
   Tutorials/index.rst
   quippy.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
* `Module Source Code <_modules/index.html>`_
