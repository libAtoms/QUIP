Tutorial
********

.. currentmodule:: quippy

This tutorial assumes you've successfully installed quippy and that all the
tests pass. If not, see :ref:`installation`. 

If you're new to Python, you might like to start with the `offical
Python tutorial
<http://docs.python.org/dev/tutorial/index.html>`_. There's also a
good `introductory tutorial
<https://wiki.fysik.dtu.dk/ase/python.html>`_ as part of the ASE
documentation. The numpy project also provides an (unfinished) `tutorial
<http://www.scipy.org/Tentative_NumPy_Tutorial>`_. If you follow this
tutorial, be aware that quippy's :class:`FortranArray` class uses one-
rather than zero-based indexing, in order to fit in better with
Fortran.

There's a couple of different ways you can use quippy: either by
typing commands interactively in Python (or, better, `ipython
<http://ipython.scipy.org>`_), or by writing a script. It's easier to
play around with things interactively, which is what we'll do in this
tutorial, but when you want to do something more than once it's worth
saving the commands in a script file. Just create a text file with the
extension ``.py``.


Getting Started
---------------

Let's start by firing up Python and importing quippy::

   $ python
   >>> from quippy import *

The ``>>>`` is the Python prompt, at which you type commands. Here
we've asked Python to import all the names defined in quippy.

Now, let's create an :class:`Atoms` object. This is the fundamental
class in quippy, and it represents a collection of atoms which together
make up a crystal or molecule. We'll make an 8-atom silicon unit cell::

   >>> dia = diamond(5.44, 14)

The first argument to :func:`diamond` is the lattice constant in Angstrom,
and the second is the atomic number. Let's poke around a little inside
our new Atoms object. Here are two equivalent ways to determine the
number of atoms::

    >>> print len(at)
    8
    >>> print at.n
    8

If you have the optional :mod:`atomeye` module installed, you can visualise
your new Atoms object::

     >>> dia.show()

.. image:: si8.png
   :align: center

An AtomEye window should pop up with a 3D representation of the
silicon crystal. Right clicking on an atom prints information about it::

   frame 0, atom 3 clicked
   Z = 14
   pos = [ 2.72  2.72  0.  ] (norm 3.846661)
   species = Si
    
The positions and atomic numbers of each atom are available in the ``pos``
and ``z`` :attr:`~Atoms.properties`::

    >>> print dia.z
    [14 14 14 14 14 14 14 14]
    >>> print dia.pos
    [[ 0.    0.    0.  ]
     [ 1.36  1.36  1.36]
     [ 2.72  2.72  0.  ]
     [ 4.08  4.08  1.36]
     [ 2.72  0.    2.72]
     [ 4.08  1.36  4.08]
     [ 0.    2.72  2.72]
     [ 1.36  4.08  4.08]]

Atomic properties are stored in arrays (actually, in a special type of array
called a :class:`FortranArray`), so it's easy to access parts of the data.
Array indexing and slicing works in the same way as in Fortran, except that
square brackets are used::

   >>> print dia.pos[1,1]   # x coordinate of atom 1
   0.0
   >>> print dia.pos[:,1]   # position of atom 1
   [ 0.  0.  0.]
   >>> print dia.pos[1]     # alternative syntax for position of atom 1
   [ 0.  0.  0.]
   >>> print dia.pos[1,:]   # all the x coordinates
   [ 0.    1.36  2.72  4.08  2.72  4.08  0.    1.36]
   >>> print dia.z[1:3]        # atomic numbers of atoms 1 to 3 inclusive
   [14 14 14]

You can also do fancier indexing as we'll see below.


Alpha Quartz
------------

Let's make a more complex structure to play with, a :func:`supercell`
of 3 x 3 x 3 alpha quartz unit cells::

   >>> unit = alpha_quartz(4.92, 5.40, 0.4697, 0.4135, 0.2669, 0.1191)
   >>> aq = supercell(unit, 3, 3, 3)

.. image:: alphaquartz.png
   :align: center

This cell contains 243 atoms. Let's look at the lattice vectors, which
are stored in a 3 x 3 matrix with ``a = aq.lattice[:,1]``, ``b =
aq.lattice[:,2]`` and ``c = aq.lattice[:,3]``::

   >>> print aq.lattice
   [[  7.38         7.38         0.        ]
    [-12.78253496  12.78253496   0.        ]
    [  0.           0.          16.2       ]]

We can convert from a cartesian representation to cell lengths and
angles to confirm that the cell has trigonal symmetry, using a
function called :func:`get_lattice_params`. To get online help on any
anything from within Python, you type ``help(name)``. In ipython, it's
even easier, just postfix the name with a question mark (or two
question marks for even more information). Here's the help for
:func:`get_lattice_params`::

   >>> help(get_lattice_params)
   Help on function get_lattice_params in module quippy.xyz_netcdf:

   get_lattice_params(lattice)
       Return 6-tuple of cell lengths and angles a,b,c,alpha,beta,gamma

From this we can see that this function takes a 3x3 ``lattice`` matrix
and returns six numbers: the three cell lengths and three cell
angles. Let's call it for the lattice associated with our new
alpha-quartz Atoms object::

   >>> a, b, c, alpha, beta, gamma =  get_lattice_params(aq.lattice)
   >>> print a, b, c
   14.76 14.76 16.2
   >>> print alpha*180.0/PI, beta*180.0/PI, gamma*180.0/PI
   90.0 90.0 120.0

So we have :math:`a = b \ne c`, and :math:`\alpha = \beta \ne \gamma`,
i.e. trigonal symmetry. Let's try out some of the more advanced indexing
operations. Firstly, logical operations on arrays return arrays, which
can be used for indexing::

   >>> (aq.z == 8).count()      # number of Oxygen atoms
   162
   >>> (aq.z == 14).count()     # number of Silicon atoms
   81
   >>> aq.pos[aq.z == 8]        # positions of the Oxygen atoms
   ...

Secondly, it's also possible to specify a list of indices, and to use
negative numbers to count backwards fromt the end of an array::

   >>> print aq.pos[:, [1,-1]]   # positions of first and last atoms
   [[ 1.155462   -2.00131889  3.6       ]
    [ 9.544062   -1.7618594   8.35686   ]]

Here's an example of making a new Atoms object containing only
the atoms with positive y coordinates. We first ensure all atoms
are within the first unit cell using :meth:`~Atoms.map_into_cell`::

   >>> aq.map_into_cell()
   >>> aq2 = aq.select(aq.pos[2,:] > 0)

.. image:: alphaquartz2.png
   :align: center

.. note:: 

   AtomEye has a different idea of the first unit-cell than quippy.
   I've shifted the cell by lattice coordinages of ``(0.5, 0.5, 0.5)``
   to make this image.  With the extended version of AtomEye included
   in quippy, you can do this by pressing :kbd:`Shift+z`.


A simple frame tool
-------------------

One common use of quippy is for post-processing the results of
simulations conducted with QUIP and libAtoms. It can read and write
Atoms objects to and from a variety of file formats, with extended XYZ
and NetCDF being the natively supported formats. Let's see how we
might construct a simple script to load analyse some simulation data.

First since this is an example let's generate some random input data::

   a0 = supercell(diamond(5.44, 14), 2, 2, 2)   # 64 atom Si cell
   L = []
   for i in range(100):
       a = a0.copy()
       a.pos += numpy.random.uniform(-0.1, 0.1, 3*a.n).reshape((3,a.n))
       L.append(a)
   al = AtomsList(L)
   al.write('test-data.xyz')
       
In this code we save 100 slightly distored silicon cells to the file
:file:`test-data.xyz`. We build up a Python list of the
configurations, which we then convert to an :class:`AtomsList`
object. This class knows how to write itself to an XYZ file.
To read the file back in and loop over the frames, we would do
something like::

   for at in AtomsList('test-data.xyz'):
      pass



Error Handling
--------------

A :exc:`RuntimeError` is raised if anything goes wrong within the
Fortran routine. You can catch this exception if necessary, with a
``try...except`` statement, for example::

   try:
      n = a.n_neighbours(i) # get number of neighbours of atom i
   except RuntimeError, e:
      print 'Something went wrong:', str(e)


