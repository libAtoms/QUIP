%%%%%%%%%%%%%%%%%%%%%%
% empirical_radius.m %
%%%%%%%%%%%%%%%%%%%%%%

number;

% http://www.webelements.com/webelements/properties/text %
% /definitions/atomic-radius-emp.html %

% The idea is that for a bond A-B, the atomic radius of A added to the atomic %
% radius of B will give a reasonable estimate for the A-B separation in %
% whatever environment. A single set of radii is very useful for most %
% purposes, however, for very accurate work adjustments would have to be made %
% to the values quoted to reflect the specific environment of the element %
% (such as coordination number). Values are given to an accuracy of about 5 %
% pm. %

% J.C. Slater, J. Chem. Phys. 1964, 39, 3199. %

% These values derived by J.C. Slater are an empirical set of atomic radii %
% derived by the careful comparison of bond lengths in over 1200 bond types %
% in ionic, metallic, and covalent crystals and molecules. %

% [Angstrom] %
empirical_radius = [
35e-2   % Slater's original 25e-2 is too small
1.785   % He: inferred from lattice constant
145e-2
105e-2
85e-2
72e-2
65e-2
60e-2
50e-2
1.5662  % Ne: inferred from lattice constant
180e-2
150e-2
    1.4255 % Al: inferred from Al-Al bond distance (1st and 2nd shell)
107e-2  % Si: inferred from lattice constant
100e-2
100e-2
100e-2
    1.8597 % Ar: inferred from lattice constant
220e-2
180e-2
160e-2
140e-2
    1.51995 % V: inferred from lattice constant (2nd and 3rd shell)
    1.44225 % Cr: inferred from lattice constant (2nd and 3rd shell)
140e-2
    1.43325 % Fe: inferred from Fe-Fe lattice constant (2nd and 3rd shell)
135e-2
135e-2
    1.2780 % Cu: inferred from Cu-Cu bond distance (1st and 2nd shell)
135e-2
130e-2
125e-2
115e-2
115e-2
115e-2
    2.0223  % Kr: inferred from lattice constant
235e-2
200e-2
180e-2
155e-2
    1.6504  % Nb:  inferred from lattice constant (2nd and 3rd shell)
    1.3872  % Mo:  for nice splitting coloring
135e-2
130e-2
135e-2
    140e-2 % Palladium:
160e-2
155e-2
155e-2
145e-2
145e-2
140e-2
140e-2
    2.1920 % Xe: inferred from lattice constant
260e-2
215e-2
195e-2
185e-2
185e-2
185e-2
185e-2
185e-2
185e-2
180e-2
175e-2
175e-2
175e-2
175e-2
175e-2
175e-2
175e-2
155e-2
    1.6529 % Ta: inferred from lattice constant (2nd and 3rd shell)
    1.5826 % W:  inferred from lattice constant (2nd and 3rd shell)
135e-2
130e-2
135e-2
135e-2
135e-2
150e-2
190e-2
180e-2
160e-2 % Bi: 
190e-2
160e-2 % At: set the same as Bi
NVL
NVL
215e-2
195e-2
180e-2
180e-2
175e-2
175e-2
175e-2
175e-2
NVL
NVL
160e-2 % Cf: set the same as Bi
160e-2 % Es: set the same as Bi
NVL
NVL
NVL
NVL
NVL
160e-2 % Db: set the same as Bi
NVL
NVL
NVL
NVL
];
% [Angstrom] %

if size(empirical_radius) ~= ZMAX
  fprintf (1, 'empirical_radius has %d elements\n', ...
      size(empirical_radius));
  return;
else
  fprintf (1, 'empirical_radius added ..\n');
end;
