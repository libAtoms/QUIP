Introduction
************

Scope of the tutorial
=====================

In this tutorial, we will carry out classical and hybrid multiscale
QM/MM (quantum mechanics/molecular mechanics) molecular dynamics
simulations using the Learn on the Fly (LOTF) schmee for fracture in
silicon. For the classical simulations we will use the
Stillinger-Weber [Stillinger1985]_ interatomic potential, which
provides a generally good description of many properties of silicon,
but, not, however of brittle fracture as we will see. The focus of the
tutorial is on embedding techniques, so we will use an approximate
quantum mechanical method which demonstrates the approach while
remaining fast enough to carry out calculations during the time
available, namely Density Functional Tight Binding [Elstner1998]_.

The tutorial is divided into three sections. There are also some
extension tasks at the end which you can complete if you have time. In
the first section we will prepare the model fracture system for our
simulations. In the second part, we will carry out classical molecular
dynamics, and in the third part of the tutorial, the 'Learn on the
Fly' ([DeVita1998]_, [Csanyi2004]_) embedding scheme will be used to
carry out coupled classical and quantum calculations. One of the main
advantages of the LOTF approach is that the QM region can be moved
during the simulation (*adaptive* QM/MM). You can read more details
about multiscale embedding methods applied to materials systems in
[Bernstein2009]_, and more about brittle fracture in silicon
investigated with the 'Learn on the Fly' technique in [Kermode2008]_.

.. _practical:

Practical considerations
========================

In this tutorial we will carry out simulations using a combination of
two packages: `quippy <http://www.jrkermode.co.uk/quippy>`_, a Python
interface to the `libAtoms/QUIP <http://www.libatoms.org>`_ MD code
developed at King's College London, Cambridge University, the Naval
Research Lab in Washington and the Fraunhofer IWM in Freiburg, and the
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

To run code, we will be using `ipython <http://ipython.org>`_, an
interactive Python shell. You can start `ipython` with a simple shell
command::
   
   $ ipython

Start by importing everything from the `quippy` package with the
command::

   In [1]: from qlab import *

.. note::

   If you have not installed the `AtomEye` plugin then the command
   above will give an :exc:`ImportError`. You can use::

     from quippy import *

   instead, which only imports `quippy` and not `AtomEye`, but you
   will then not be able to visualise interactively using the
   :func:`~qlab.view` function - you can always save configurations to
   external files and visualise them with other tools such as `VMD
   <http://www.ks.uiuc.edu/Research/vmd/>`_ or `OVITO
   <http://www.ovito.org>`_ instead. See the :mod:`qlab` and
   :mod:`atomeye` documentation for further details.

You should prepare your Python scripts in a text editor and then run
them from within `ipython` with the `%run
<http://ipython.org/ipython-doc/stable/interactive/tutorial.html#running-and-editing>`_
command. For example::

   In [2]: run make_crack.py

If your script is taking a long time to run and you realise something
is wrong you can interrupt it with `Ctrl+C`.

The tutorial is structured so that you will build up your script as
you go along. You should keep the same `ipython` session open, and
simply run the script again each time you add something new. You will
find it easier to follow the tutorial if you use the same variable
names as we have used here.

You can also execute individual Python statements by typing them
directly into the `ipython` shell. This is very useful for debugging
and for interactive visualisation (described in more detail in the
:ref:`visualisation2` section). Finally, it's possible to copy and
paste code into `ipython` using the `%paste` or `%cpaste` commands.

You can follow the links in the tutorial to read more documentation for
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
spend on it. At the end of each subsection there is a *Milestone*
where a script for that stage of the tutorial is provided. If you run
out of time, just skip ahead to the next milestone, download the
script and then continue with the next section.

Continue with :ref:`step1`.
