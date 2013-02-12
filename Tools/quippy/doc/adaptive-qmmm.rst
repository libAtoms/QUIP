.. currentmodule:: quippy

************************************************
Adaptive QM/MM Simulation of Fracture in Silicon
************************************************

**DRAFT -- not yet ready for public use**

This tutorial has been prepared for use at a hands-on session at the
`ADGLASS Winter School on Advanced Molecular Dynamics Simulations
<http://www.adglass.eu/adglass_news.html>`_, ICTP, Trieste, February
2013.

:Authors: James Kermode and Gianpietro Moras
:Date: February 2013

Scope of the tutorial
=====================

In this tutorial, we will carry out classical and hybrid multiscale
QM/MM (quantum mechanics/molecular mechanics) molecular dynamics
simulations of fracture in silicon. For the classical simulations we
will use the Stillinger-Weber [SW]_ interatomic potential, which
provides a generally good description of many properties of silicon,
but, not, however of brittle fracture as we will see. Since the focus
of this tutorial is on embedding schemes, we will use an approximate
quantum mechanical method which demonstrates the technique while
remaining fast enough to carry out calculations during the time
available, namely Density Functional Tight Binding [DFTB]_.

The tutorial is divided into three sections. At the end of each
section there are extension tasks which you can complete if you have
time. In the first section we will prepare the model fracture system
for our simulations. In the second part, we will carry out classical
molecular dynamics, and in the third part of the tutorial, the 'Learn
on the Fly' (LOTF) [Csanyi2004]_ embedding scheme will be used to
carry out coupled classical and quantum calculations. One of the main
advantages of the LOTF approach is that the QM region can be moved
during the simulation, which is the origin of the term *adaptive*
QM/MM. You can read more details about multiscale embedding methods
applied to materials systems in [Bernstein2009]_, and more about brittle
fracture in silicon investigated with the 'Learn on the Fly' technique
in [Kermode2008]_.

Practical considerations
========================

In this tutorial we will carry out simulations using a combination of
two packages: `quippy <http://www.jrkermode.co.uk/quippy>`_, a Python
interface to the `libAtoms/QUIP <http://www.libatoms.org>`_ MD code
developed at King's College London, Cambridge University, the Naval
Research Lab in Washington and the Fraunhofer IWM in Freiberg, and the
Atomic Simulation Environment, `ASE <https://wiki.fysik.dtu.dk/ase>`_,
a Python framework developed at the Centre for Atomistic Materials
Design (`CAMd <http://www.camd.dtu.dk/>`_) at DTU, Copenhagen.

If you're not familiar with Python don't worry, it's quite easy to
pick up and the syntax is very similar to Fortran. We have provided
template Python code which you should copy and paste into your text
editor. Save your script with a ``.py`` file extension. Whenever you
see ``...`` it means there is something you need to add to the code in
order to complete it; usually there will be a comment to give you a
hint. The source code listings are all in boxes like this::

   print 'Example code listing'     # Print something out
   print ...                        # You need to add something here
                                    # e.g. replace the ... with a number

To run code, we will be using `ipython <http://ipython.org>`_, an interactive Python
shell. You can start `ipython` with a simple shell command::
   
   $ ipython

Start by importing everything from the `quippy` package with the
command::

   In [1]: from qlab import *

You should prepare your Python scripts in a text editor and then run
then from within `ipython` with the `run
<http://ipython.org/ipython-doc/stable/interactive/tutorial.html#running-and-editing>`_
command. For example::

   In [2]: run make_crack.py

The tutorial is structured so that you will build up your script as
you go along. You should keep the same `ipython` session open, and
simply run the script again each time you add something new. You will
find it easier to follow the tutorial if you use the same variable
names etc. as given those here.

You can also execute individual Python statements by typing them
directly into the `ipython` shell. This is very useful for debugging
and for interactive visualisation (described in more detail
below). Finally, it's possible to copy and paste code into ipython
using the `paste` and `cpaste` commands.

You can follow the links below to read more documentation for
functions and classes we will be using, and for `quippy` routines you
can also browse the source code by clicking on the `[source]` link
which you will find to the right of class and function definitions in
the online documentations. You can also explore functions
interactively from within `ipython` by following the name of an object
or function with a ``?``, e.g. ::

   In [3]: Atoms ?

displays information about the :class:`~.Atoms` class, or ::

   In [4]: write ?

displays information about the :func:`~quippy.io.write` function. See
the `ipython tutorial
<http://ipython.org/ipython-doc/stable/interactive/tutorial.html>`_
for more information.

Each subsection indicates the approximate amount of time you should
spend on it. At the end of each subsection there is a **Milestone**
where an script for that stage of the tutorial is provided. If you run
out of time, just skip ahead to the next milestone, download the
script and then continue with the next section.

Tutorial contents
=================

.. toctree::
   :maxdepth: 2

   adaptive-qmmm-step1.rst
   adaptive-qmmm-step2.rst
   adaptive-qmmm-step3.rst


References
==========

.. [SW] Stillinger, F. H., & Weber, T. A. Computer simulation
   of local order in condensed phases of silicon. Physical Review B,
   31(8),
   5262–5271. (1985). http://link.aps.org/doi/10.1103/PhysRevB.31.5262

.. [DFTB] Elsterner, M., Porezag, D., Jungnickel, G., Elsner, J.,
   Haugk, M., Frauenheim, T., Suhai, S., et
   al. Self-consistent-charge density-functional tight-binding
   method for simulations of complex materials
   properties. Phys. Rev. B. 58 7260 (1998).
   http://prola.aps.org/abstract/PRB/v58/i11/p7260_1

.. [Csanyi2004] Csányi, G., Albaret, T., Payne, M., & De Vita,
   A. 'Learn on the Fly': A Hybrid Classical and Quantum-Mechanical
   Molecular Dynamics Simulation. Physical Review Letters,
   93(17), 175503. (2004)
   http://prl.aps.org/abstract/PRL/v93/i17/e175503>

.. [Bernstein2009] Bernstein, N., Kermode, J. R., & Csányi,
   G. Hybrid atomistic simulation methods for materials
   systems. Reports on Progress in Physics,
   72(2), 026501 (2009). http://dx.doi.org/10.1088/0034-4885/72/2/026501

.. [Kermode2008] Kermode, J. R., Albaret, T., Sherman, D., Bernstein,
   N., Gumbsch, P., Payne, M. C., Csányi, G., and A. De Vita. Low-speed
   fracture instabilities in a brittle crystal. Nature, 455,
   1224–1227 (2008). http://www.nature.com/doifinder/10.1038/nature07297

