Tutorial
********

Other Tutorials
---------------

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


Error Handling
--------------

A :exc:`RuntimeError` is raised if anything goes wrong within the
Fortran routine.


Overview
--------

Important classes:

#. :class:`~quippy.Atoms`
#. :class:`~quippy.FortranArray`
#. :class:`~quippy.atomeye.AtomEyeView`
#. :class:`~quippy.AtomsList`
#. :class:`~quippy.Dictionary`
#. :class:`~quippy.Potential`
#. :class:`~quippy.DynamicalSystem`
#. :class:`~quippy.MetaPotential`

