%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% electronegativity_Pearson.m %
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

% [eV] %
electronegativity_Pearson = [
7.18 
12.3 
3.01 
4.9 
4.29 
6.27 
7.30 
7.54 
10.41 
10.6 
2.85 
3.75 
3.23 
4.77 
5.62 
6.22 
8.30 
7.70 
2.42 
2.2 
3.34 
3.45 
3.6 
3.72 
3.72 
4.06 
4.3 
4.40 
4.48 
4.45 
3.2 
4.6 
5.3 
5.89 
7.59 
6.8 
2.34 
2.0 
3.19 
3.64 
4.0 
3.9 
3.91 
4.5 
4.30 
4.45 
4.44 
4.33 
3.1 
4.30 
4.85 
5.49 
6.76 
5.85 
2.18 
2.4 
3.1 
3.0 
3.0 
3.0 
3.0 
3.1 
3.1 
3.3 
3.2 
NVL 
3.3 
3.3 
3.4 
3.5 
3.0 
3.8 
4.11 
4.40 
4.02 
4.9 
5.4 
5.6 
5.77 
4.91 
3.2 
3.90 
4.69 
5.16 
6.2 
5.1 
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
3.5 
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

if size(electronegativity_Pearson) ~= ZMAX then
  fprintf (1, 'electronegativity_Pearson has %d elements\n', size(electronegativity_Pearson));
  break;
else
  fprintf (1, 'electronegativity_Pearson added ..\n');
end;
