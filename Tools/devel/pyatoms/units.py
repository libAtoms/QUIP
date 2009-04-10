# Some constants for the eV/A/fs system of units

ELECTRONMASS_GPERMOL =  5.48579903e-4 # grams/mol
ELEM_CHARGE = 1.60217653e-19 # coulombs
HARTREE = 27.2113961 # eV
RYDBERG = 0.5*HARTREE # eV                                                
BOHR = 0.529177249 # Angstrom                                             
HBAR_EVSEC = 6.5821220e-16 # hbar in eV seconds                           
HBAR_AU = 1.0              # hbar in a.u.                                 
HBAR = (HBAR_EVSEC*1e-15)    # hbar in eV fs                              
ONESECOND = 1e15           # 1 second in fs                               
ONESECOND_AU = (1.0/(HBAR_EVSEC/(HBAR_AU*HARTREE))) # 1 second in a.u.    
AU_FS = (1.0/ONESECOND_AU*ONESECOND) # a.u. time in fs                    
MASSCONVERT = (1.0/ELECTRONMASS_GPERMOL*HARTREE*AU_FS*AU_FS/(BOHR*BOHR))   
BOLTZMANN_K = 8.617385e-5 # eV/Kelvin                                     
PI = 3.141592653589793238462643383                                         
N_A = 6.0221479e23 # Avogadro's number                                    
KCAL_MOL = 4.3383e-2 # eV                                                 
DEGREES_PER_RADIAN = 180.0 / PI                                            
RADIANS_PER_DEGREE = PI / 180.0                                            
GPA = 1.6022e-19*1.0e30/1.0e9 
