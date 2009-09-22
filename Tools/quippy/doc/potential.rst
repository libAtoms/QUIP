:class:`quippy.Potential` -- evaluate interatomic potentials
============================================================

.. currentmodule:: quippy

.. class:: Potential(args_str, param_str)

   A :class:`Potential` object represents an interatomic potential, a
   tight binding model or an interface to an external code used to
   perform calculations. It is initialised from an `args_str` describing
   the type of potential, and an XML formatted string `param_str`
   giving the parameters.

   Types of Potential:

     ====================  ==========================================================
     `arg_str` prefix      Description
     ====================  ==========================================================   
     ``IP``                Interatomic Potential
     ``TB``                Tight Binding Model
     ``FilePot``           File potential, used to communicate with external program
     ====================  ==========================================================

   
   Types of interatomic potential available:

     ======================= ==========================================================
     ``arg_str`` prefix      Description
     ======================= ==========================================================   
     ``IP ASAP``             ASAP quartz potential
     ``IP BOP``              Bond order potential for metails
     ``IP Brenenr``          Brenner (1990) potential for carbon
     ``IP Brenner_2002``     Brenner (2002) reactive potential for carbon
     ``IP Brenner_Screened`` Pastewka modified Brenner reactive potential for carbon
     ``IP EAM_ErcolAd``      Embedded atom potential of Ercolessi and Adams
     ``IP FB``               Flikkema and Bromley potential for SiO2
     ``IP FS``	  	     Finnis=Sinclair potential for metals
     ``IP GAP``              Gaussian approximation potential
     ``IP LJ``               Lennard-Jones potential
     ``IP Si_MEAM``          Silicon modified embedded attom potential
     ``IP SW``               Stillinger-Weber potential for silicon
     ``IP Tersoff``          Tersoff potential for silicon
     ======================= ==========================================================   

   Types of tight binding potential available:

     ======================= ==========================================================
     ``arg_str`` prefix      Description
     ======================= ==========================================================   
     ``TB Bowler``           Bowler tight binding model 
     ``TB DFTB``             Density functional tight binding
     ``TB GSP``              Goodwin-Skinner-Pettifor tight binding model
     ``TB NRL_TB``           Naval Research Laboratory tight binding model
     ======================= ==========================================================   

   
   Examples of the XML parameters for each of these potential can be
   found in the `QUIP_Core/parameters
   <http://src.tcm.phy.cam.ac.uk/viewvc/jrk33/repo/trunk/QUIP/QUIP_Core/parameters>`_
   directory of the QUIP svn repository.

   Important methods of Potential objects are listed below.

   .. method:: cutoff

      Return the cutoff of this Potential, in Angstrom. This is the
      minimum neighbour connectivity cutoff that should be used: if
      you're doing MD you'll want to use a slightly larger cutoff so
      that new neighbours don't drift in between connectivity updates.

   .. method:: calc(at[, e, local_e, f, df, virial, calc_energy, calc_force, calc_local_e, calc_df, calc_virial, args_str])

      Apply this Potential to the Atoms object `at`, which must have
      connectivity information (i.e. :meth:`Atoms.calc_connect` should
      have been called). The optional arguments determine what should be
      calculated and how it will be returned. Each physical quantity
      has two correspond optional arguments. The first should be an array of
      the correct shape to receive the quantity in question, as set out
      in the table below. 

	================ =================== ============= ===========
	Array argument   Logical argument    Quantity      Shape
	================ =================== ============= ===========
	`e`              `calc_energy`	     Energy        ``()``
	`local_e`        `calc_local_e`	     Local energy  ``(at.n,)``
	`f`              `calc_force`	     Force         ``(3,at.n)``
	`virial`         `calc_virial`       Virial tensor ``(3,3)``
	================ =================== ============= ===========

      Note that for the energy `e`, a rank-0 array is required to
      store the results of this calculation. This is necessary since
      the underlying Fortran argument is both ``optional`` and
      ``intent(out)``. We can avoid this inconvenience using the
      `calc_energy` argument which returns the result in the `energy`
      parameter of the Atoms object `at`. The same is true for the
      other `calc` logical arguments listed in the table above.

      The `args_str` argument can be a string or dictionary containing
      additional arguments which depend on the particular Potential
      being used. As mentioned in the documentation of
      :meth:`oo_fortran.FortranDerivedType._runroutine`, as a special
      case unexpected keyword arguments will be converted to
      additional `args_str` options (in fact this is how the
      `calc_energy` etc.  arguments are implemented).

      Not all Potentials support all of these quantities: a
      :exc:`RuntimeError` will be raised if you ask for something that
      is not supported.
 