%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% electronegativity_Allred.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

electronegativity_Allred = [
2.2 
NVL 
1.0 
1.5 
2.0 
2.5 
3.1 
3.5 
4.1 
NVL 
1.0 
1.2 
1.5 
1.7 
2.1 
2.4 
2.8 
NVL 
0.9 
1.0 
1.2 
1.3 
1.5 
1.6 
1.6 
1.6 
1.7 
1.8 
1.8 
1.7 
1.8 
2.0 
2.2 
2.5 
2.7 
NVL 
0.9 
1.0 
1.1 
1.2 
1.2 
1.3 
1.4 
1.4 
1.5 
1.4 
1.4 
1.5 
1.5 
1.7 
1.8 
2.0 
2.2 
NVL 
0.9 
1.0 
1.1 
1.1 
1.1 
1.1 
1.1 
1.1 
1.0 
1.1 
1.1 
1.1 
1.1 
1.1 
1.1 
1.1 
1.1 
1.2 
1.4 
1.4 
1.5 
1.5 
1.6 
1.4 
1.4 
1.5 
1.4 
1.6 
1.7 
1.8 
2.0 
NVL 
0.9 
1.0 
1.0 
1.1 
1.1 
1.2 
1.2 
1.2 
1.2 
1.2 
1.2 
1.2 
1.2 
1.2 
1.2 
NVL 
NVL
NVL 
NVL 
NVL 
NVL 
NVL 
NVL
];

if size(electronegativity_Allred) ~= ZMAX then
  fprintf (1, 'electronegativity_Allred has %d elements\n', size(electronegativity_Allred));
  break;
else
  fprintf (1, 'electronegativity_Allred added ..\n');
end;
