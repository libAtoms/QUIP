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

Atoms objects
=============

.. automodule:: quippy.atoms
   :synopsis: Representation of atomic configurations

.. autoclass:: Atoms
   :members:
   :inherited-members:
   :show-inheritance:

.. autofunction::     make_lattice
.. autofunction::     get_lattice_params
.. autofunction::     bond_length
.. autofunction::     termination_bond_rescale
.. autofunction::     cell_volume
.. autofunction::     map_into_cell
.. autofunction::     parse_atom_mask
