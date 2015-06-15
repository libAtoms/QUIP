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

The QUIP package is a software library written in Fortran 95+ for the
purposes of carrying out molecular dynamics simulations. The code
implements a wide variety of interatomic potentials and tight binding
quantum mechanics, and is also able to call external packages. Various
hybrid combinations are also supported in the style of QM/MM.

quippy is a Python interface to the `libAtoms/QUIP <https://github.com/libAtoms/QUIP>`_
molecular dynamics framework.  

The main libAtoms/QUIP contributors are:

* University of Cambridge: Albert P. Bartók, Gábor Csányi, Alan Nichol, Letif Mones, Wojciech Szlachta
* University of Warwick: James Kermode
* King's College London: Alessandro De Vita
* Naval Research Laboratory, Washington DC: Noam Bernstein
* Fraunhofer IWM, Freiburg: Lars Pastewka

The quippy interface is principally maintained by `James Kermode <http://www.jrkermode.co.uk>`_.  

Features
========

The following interatomic potentials are presently coded or linked in QUIP:

* EAM (fcc metals)
* Fanourgakis-Xantheas (water)
* Finnis-Sinclair (bcc metals)
* Flikkema-Bromley
* GAP (Gaussian Approximation Potentials: general many-body)
* Guggenheim-!McGlashan
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

The following external packages can be called:

* CASTEP
* VASP
* CP2K
* ASAP
* ASE (latest svn trunk recommended)
* Molpro

Contents
========

.. toctree::
   :maxdepth: 2

   intro.rst
   install.rst
   Examples/index.rst
   tutorials.rst
   quippy.rst
   visualisation.rst
   fortran_wrapper.rst
   

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
* `Module source code <_modules/index.html>`_
