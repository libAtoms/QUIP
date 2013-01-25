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

:mod:`quippy.atoms` --- representation of atomic configurations
===============================================================

.. automodule:: quippy.atoms
   :synopsis: Representation of atomic configurations

Atoms objects
-------------

.. autoclass:: Atoms
    :members:
    :inherited-members:
    :show-inheritance:

.. autosummary::
    make_lattice
    get_lattice_params
    bond_length
    termination_bond_rescale
    cell_volume
    map_into_cell
    parse_atom_mask

.. autofunction::     make_lattice
.. autofunction::     get_lattice_params
.. autofunction::     bond_length
.. autofunction::     termination_bond_rescale
.. autofunction::     cell_volume
.. autofunction::     map_into_cell
.. autofunction::     parse_atom_mask


Connection objects
------------------

.. autoclass:: Connection
    :members:
    :inherited-members:
    :show-inheritance:

.. autosummary::
    fit_box_in_cell
    divide_cell
    get_min_max_images
    max_cutoff

.. autofunction::       fit_box_in_cell
.. autofunction::       divide_cell
.. autofunction::       get_min_max_images
.. autofunction::       max_cutoff


DomainDecomposition objects
---------------------------

.. autoclass:: DomainDecomposition
    :members:
    :inherited-members:
    :show-inheritance:

