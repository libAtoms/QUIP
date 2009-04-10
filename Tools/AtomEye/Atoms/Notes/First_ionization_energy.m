%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first_ionization_energy.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

number;

% The first ionisation energy is the energy which is necessary to remove one %
% electron from an atom in the ground state. In general metals have a low %
% ionisation energy and non-metals have a high ionisation energy. The noble %
% gases Helium and Neon have the highest ionisation energies. Caesium and %
% Francium have the lowest ionisation energies. The ionisation energies in %
% the table below are given in units of electron volts (eV).   %

% http://klbproductions.com/yogi/periodic/ %

% [eV] %
first_ionization_energy = [
13.598 
24.587 
5.392 
9.322 
8.298 
11.260 
14.534 
13.618 
17.422 
21.564 
5.139 
7.646 
5.986 
8.151 
10.486 
10.360 
12.967 
15.759 
4.341 
6.113 
6.54 
6.82 
6.74 
6.766 
7.435 
7.870 
7.86 
7.635 
7.726 
9.394 
5.999 
7.899 
9.81 
9.752 
11.814 
13.999 
4.177 
5.695 
6.38 
6.84 
6.88 
7.099 
7.28 
7.37 
7.46 
8.34 
7.576 
8.993 
5.786 
7.344 
8.641 
9.009 
10.451 
12.130 
3.894 
5.212 
5.577 
5.47 
5.42 
5.49 
5.55 
5.63 
5.67 
6.14 
5.85 
5.93 
6.02 
6.10 
6.18 
6.254 
5.426 
7.0 
7.89 
7.98 
7.88 
8.7 
9.1 
9.0 
9.225 
10.437 
6.108 
7.416 
7.289 
8.42 
9.5 
10.748 
4.0 
5.279 
6.9 
6.95 
NVL 
6.08 
NVL 
5.8 
6.0 
NVL 
NVL 
NVL 
NVL 
NVL 
NVL 
NVL 
NVL 
NVL 
NVL 
NVL 
NVL 
NVL 
NVL
];
% [eV] %

if size(first_ionization_energy) ~= ZMAX
  fprintf (1, 'first_ionization_energy has %d elements\n', ...
      size(first_ionization_energy));
  return;
else
  fprintf (1, 'first_ionization_energy added ..\n');
end;
