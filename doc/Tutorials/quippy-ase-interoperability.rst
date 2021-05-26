Interoperability with Atomic Simulation Environment
===================================================

-  quippy uses the standard ``ase.atoms.Atoms`` class to represent
   Atomic structures
-  quippy ``Potential`` objects can be used as ASE calculators, and
   vice-versa
-  Can use standard ASE tools, plus communicate with other packages
   using ASE as *lingua franca*

Example: vacancy formation energy
---------------------------------

-  Generate structure with ``ASE`` lattice tools
-  Stillinger-Weber potential implementation from ``QUIP``
-  Elastic constant fitting routine from ``matscipy``, internal
   relaxations with ``ASE`` FIRE minimiser

.. code:: ipython3

    %pylab inline
    from ase.build import bulk
    from ase.optimize import BFGS
    from ase.optimize.precon import PreconLBFGS
    from quippy.potential import Potential
    
    si = bulk('Si', a=5.44, cubic=True)
    sw_pot = Potential('IP SW', 
                       param_filename='../../share/Parameters/ip.parms.SW.xml')
    si.set_calculator(sw_pot)
    e_bulk_per_atom = si.get_potential_energy()/len(si)
    
    # call general purpose elastic constants calculator 
    #   using ASE Atoms and QUIP Potential
    from matscipy.elasticity import fit_elastic_constants
    Cij = fit_elastic_constants(si, optimizer=BFGS,
                                symmetry='cubic', logfile='-')
    vac1 = si.copy()
    vac1 *= (3, 3, 3)
    half_cell = np.diag(vac1.cell)/2.
    vac_atom = ((vac1.positions - half_cell)**2).sum(axis=1).argmin()
    del vac1[vac_atom]
    
    vac1.set_calculator(sw_pot)
    vac1.rattle(0.01)
    opt = PreconLBFGS(vac1) # big cell, use preconditioned minimiser
    opt.run(fmax=1e-6)
    e_vac = vac1.get_potential_energy() - e_bulk_per_atom*len(vac1)
    print('SW vacancy formation energy', e_vac, 'eV')


.. parsed-literal::

    Populating the interactive namespace from numpy and matplotlib
          Step     Time          Energy         fmax
    BFGS:    0 21:22:10      -34.635777        0.3247
    BFGS:    1 21:22:10      -34.644590        0.1504
    BFGS:    2 21:22:10      -34.646997        0.0001
          Step     Time          Energy         fmax
    BFGS:    0 21:22:10      -34.670667        0.1584
    BFGS:    1 21:22:10      -34.672777        0.0749
    BFGS:    2 21:22:10      -34.673385        0.0000
          Step     Time          Energy         fmax
    BFGS:    0 21:22:10      -34.678737        0.0000
          Step     Time          Energy         fmax
    BFGS:    0 21:22:10      -34.660845        0.1508
    BFGS:    1 21:22:10      -34.662784        0.0742
    BFGS:    2 21:22:10      -34.663403        0.0000
          Step     Time          Energy         fmax
    BFGS:    0 21:22:10      -34.617822        0.2945
    BFGS:    1 21:22:10      -34.625262        0.1477
    BFGS:    2 21:22:10      -34.627763        0.0001
    Fitting C_11
    Strain array([-0.02, -0.01,  0.  ,  0.01,  0.02])
    Stress array([-2.56044687, -1.01247671,  0.5027424 ,  1.98366491,  3.42893711]) GPa
    Cij (gradient) / GPa    :     149.74909575669915
    Error in Cij / GPa      :     1.1696996170085603
    Correlation coefficient :     0.9999084935045536
    Setting C11 (1) to 0.934660 +/- 0.007301
    
    
    Fitting C_21
    Strain array([-0.02, -0.01,  0.  ,  0.01,  0.02])
    Stress array([-1.07663577, -0.26643655,  0.5027424 ,  1.23345414,  1.92818198]) GPa
    Cij (gradient) / GPa    :     75.0952617697293
    Error in Cij / GPa      :     1.3149075235415504
    Correlation coefficient :     0.9995404242109733
    Setting C21 (7) to 0.468708 +/- 0.008207
    
    
    Fitting C_31
    Strain array([-0.02, -0.01,  0.  ,  0.01,  0.02])
    Stress array([-1.07663577, -0.26643655,  0.5027424 ,  1.23345414,  1.92818198]) GPa
    Cij (gradient) / GPa    :     75.09526176972929
    Error in Cij / GPa      :     1.31490752354155
    Correlation coefficient :     0.9995404242109733
    Updating C31 (7) with value 0.468708 +/- 0.008207
    
    
    Fitting C_44
    Strain array([-0.02, -0.01,  0.  ,  0.01,  0.02])
    Stress array([-1.13572340e+00, -5.65842409e-01, -9.46072689e-15,  5.60142655e-01,
            1.11304586e+00]) GPa
    Cij (gradient) / GPa    :     56.23523568430933
    Error in Cij / GPa      :     0.19437884854805132
    Correlation coefficient :     0.9999820790695022
    Setting C44 (4) to 0.350993 +/- 0.001213
    
    
    [[b C11 b C12 b C12 b     b     b    ]
     [b C12 b C11 b C12 b     b     b    ]
     [b C12 b C12 b C11 b     b     b    ]
     [b     b     b     b C44 b     b    ]
     [b     b     b     b     b C44 b    ]
     [b     b     b     b     b     b C44]]
    
     = 
    
    [[149.75  75.1   75.1    0.     0.     0.  ]
     [ 75.1  149.75  75.1    0.     0.     0.  ]
     [ 75.1   75.1  149.75   0.     0.     0.  ]
     [  0.     0.     0.    56.24   0.     0.  ]
     [  0.     0.     0.     0.    56.24   0.  ]
     [  0.     0.     0.     0.     0.    56.24]]
    C_11 = 149.75 +/- 1.17 GPa
    C_12 = 75.10 +/- 1.31 GPa
    C_44 = 56.24 +/- 0.19 GPa
    PreconLBFGS:   0  21:22:10     -927.087471       0.8332
    estimate_mu(): mu=2.3315998829549316, mu_c=1.0
    PreconLBFGS:   1  21:22:10     -927.611735       0.1487
    PreconLBFGS:   2  21:22:11     -927.643943       0.0642
    PreconLBFGS:   3  21:22:11     -927.652380       0.0393
    PreconLBFGS:   4  21:22:11     -927.655827       0.0210
    PreconLBFGS:   5  21:22:11     -927.656988       0.0092
    PreconLBFGS:   6  21:22:11     -927.657131       0.0052
    PreconLBFGS:   7  21:22:11     -927.657193       0.0030
    PreconLBFGS:   8  21:22:11     -927.657220       0.0011
    PreconLBFGS:   9  21:22:11     -927.657224       0.0009
    PreconLBFGS:  10  21:22:11     -927.657226       0.0004
    PreconLBFGS:  11  21:22:11     -927.657226       0.0002
    PreconLBFGS:  12  21:22:11     -927.657226       0.0001
    PreconLBFGS:  13  21:22:11     -927.657226       0.0001
    PreconLBFGS:  14  21:22:11     -927.657226       0.0000
    PreconLBFGS:  15  21:22:11     -927.657226       0.0000
    PreconLBFGS:  16  21:22:11     -927.657226       0.0000
    PreconLBFGS:  17  21:22:11     -927.657226       0.0000
    PreconLBFGS:  18  21:22:12     -927.657226       0.0000
    PreconLBFGS:  19  21:22:12     -927.657226       0.0000
    PreconLBFGS:  20  21:22:12     -927.657226       0.0000
    PreconLBFGS:  21  21:22:12     -927.657226       0.0000
    SW vacancy formation energy 4.33384020178346 eV


