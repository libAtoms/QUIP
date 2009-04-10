%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% electronegativity_Pauling.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

number;

% The electronegativity is a measure for the attractive strength of an %
% atom on the binding electron pair of an atomic binding. The common %
% abreviation is EN. The EN is a quantity which cannot be measured %
% directly.  It is derived form other quantities like the ionisation %
% energy, the binding energy, etc. In general the EN increases from the %
% left to the right and from the botton to the top in the periodic %
% table. Fluorine is the element with the highest EN. %

% The alcali metals Potassium, Rubidium, Caseium, and Francium have the %
% lowest EN.  The electronegativity according to Bohr, Allred and %
% Rochow, respectively is a dimensionless quantity which uses some %
% quantity as a reference value. The different values of the EN %
% according to Bohr, Allred and Rochow, and according to Pearson result %
% from the different methods of calculation. The dimension of the %
% abslolute values of the EN according to Pearson is the electron volt %
% (eV). %

% http://klbproductions.com/yogi/periodic/ %

electronegativity_Pauling = [
2.20 
NVL 
0.98 
1.57 
2.04 
2.55 
3.04 
3.44 
3.98 
NVL 
0.93 
1.31 
1.61 
1.90 
2.19 
2.58 
3.16 
NVL 
0.82 
1.00 
1.36 
1.54 
1.63 
1.66 
1.55 
1.83 
1.88 
1.91 
1.90 
1.65 
1.81 
2.01 
2.18 
2.55 
2.96 
NVL 
0.82 
0.95 
1.22 
1.33 
1.6 
2.16 
1.9 
2.2 
2.28 
2.20 
1.93 
1.69 
1.78 
1.96 
2.05 
2.1 
2.66 
2.6 
0.79 
0.89 
1.10 
1.12 
1.13 
1.14 
NVL 
1.17 
NVL 
1.20 
NVL 
1.22 
1.23 
1.24 
1.25 
NVL 
1.27 
1.3 
1.5 
2.36 
1.9 
2.2 
2.20 
2.28 
2.54 
2.00 
1.62 
2.33 
2.02 
2.0 
2.2 
NVL 
0.7 
0.89 
1.1 
1.3 
1.5 
1.38 
1.36 
1.28 
1.3 
1.3 
1.3 
1.3 
1.3 
1.3 
1.3 
1.3 
1.3 
NVL 
NVL 
NVL 
NVL 
NVL 
NVL 
];

if size(electronegativity_Pauling) ~= ZMAX then
  fprintf (1, 'electronegativity_Pauling has %d elements\n', size(electronegativity_Pauling));
  break;
else
  fprintf (1, 'electronegativity_Pauling added ..\n');
end;
