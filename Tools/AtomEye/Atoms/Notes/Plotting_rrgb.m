%%%%%%%%%%%%%%%%%%%
% plotting_rrgb.m %
%%%%%%%%%%%%%%%%%%%

number;

% http://www.chemicalgraphics.com/paul/Distribution.html periodic.tab %
% http://www.chemicalgraphics.com/paul/Manual.html#CPK-Atom-Colors %

% + %

% http://www.thch.uni-bonn.de/tc/doc/viewmol/viewmol.html */
% http://www.theochem.kth.se/html/pemac/links/VIEWMOL.en.html %
% http://www.ccl.net/chemistry/resources/software/SOURCES/C/viewmol/index.shtml
% viewmol-2.2.1-1.i386.rpm + viewmol-2.2.1.bin.lesstif.tgz %
% viewmol-2.2.1/viewmolrc, /usr/local/lib/viewmol/viewmolrc %

% radius[Angstrom]  colour from r    g    b     to r    g    b  %

copper = [184 115 51] / 184;
gold = [255 215 0] / 255;
steelblue = [176 196 222] / 255;
silver = [ hex2dec('A4') hex2dec('AA') hex2dec('AD') ] / 255;

plotting_rrgb = [
% H:
      0.46   0.8    0.8    0.8
% He:
      1.60   silver
% Li:
      0.68   0.7    0.7    0.7
% Be:
      0.35   silver
% B:
      0.82   0.9    0.4    0.0
% C:
      0.85   0.35   0.35   0.35
% N:
      0.75   0.2    0.2    0.8
% O:
      0.74   0.8    0.2    0.2
% F:
      0.72   0.7    0.85   0.45
% Ne:
      1.12   silver
% Na:
      1.54   0.6    0.6    0.6
% Mg:
      1.37   0.6    0.6    0.7
% Al:
      1.33   silver
% Si:
    1.17   steelblue
% P:
      1.06   0.1    0.7    0.3
% S:
      1.02   0.95   0.9    0.2
% Cl:
      0.99   0.15   0.5    0.1
% Ar:
      1.54   silver
% K:
      1.33   0.5    0.5    0.5
% Ca:
      1.74   0.8    0.8    0.7
% Sc:
      1.91   silver
% Ti:
      1.63   silver
% V:
      1.61   silver
% Cr:
      1.55   0.00   0.80   0.00
% Mn:
      1.42   silver
% Fe:
    1.37   steelblue * 0.75
% Co:
      1.55   silver
% Ni:
    1.55   silver * 0.4
% Cu:
    1.17   0.05*copper+0.9*gold
% Zn:
      1.47   silver
% Ga:
      1.52   0.90   0.00   1.00
% Ge:
      1.53   silver
% As:
      1.55   1.00   1.     0.30
% Se:
      1.46   silver
% Br:
      1.44   0.5    0.08   0.12
% Kr:
      1.60   silver
% Rb:
      2.78   silver
% Sr:
      2.45   silver
% Y:
      2.08   silver
% Zr:
      1.89   silver
% Nb:
      1.73   silver
% Mo:
      1.56   silver
% Tc:
      1.65   silver
% Ru:
      1.63   silver
% Rh:
      1.65   silver
% Pd:
      1.68   silver
% Ag:
      1.75   silver
% Cd:
      1.79   silver
% In:
      1.93   silver
% Sn:
      1.71   silver
% Sb:
      1.75   silver
% Te:
      1.73   silver
% I:
      1.63   0.5    0.1    0.5
% Xe:
      1.70   silver
% Cs:
      2.95   silver
% Ba:
      2.47   silver
% La:
      2.17   silver
% Ce:
      1.82   0.80   0.80   0.00
% Pr:
      2.12   silver
% Nd:
      2.11   silver
% Pm:
      2.11   silver
% Sm:
      2.10   silver
% Eu:
      2.29   silver
% Gd:
      2.09   gold
% Tb:
      2.06   silver
% Dy:
      2.05   silver
% Ho:
      2.04   silver
% Er:
      2.03   silver
% Tm:
      2.03   silver
% Yb:
      2.24   silver
% Lu:
      2.02   silver
% Hf:
      1.86   silver
% Ta:
      1.73   silver
% W:
      1.58   silver
% Re:
      1.67   silver
% Os:
      1.64   silver
% Ir:
      1.66   silver
% Pt:
      1.69   silver
% Au:
      1.74   0.90   0.80   0.00
% Hg:
      1.80   silver
% Tl:
      2.00   silver
% Pb:
      2.05   silver
% Bi:
      1.85   silver
% Po:
      1.48   silver
% At:
      5.43   0.8    0.2    0.2
% Rn:
      1.90   silver
% Fr:
      3.15   silver
% Ra:
      2.52   silver
% Ac:
      2.18   silver
% Th:
      2.10   silver
% Pa:
      1.91   silver
% U:
      1.69   silver
% Np:
      0.95   silver
% Pu:
      1.81   silver
% Am:
      1.61   silver
% Cm:
      0.91   silver
% Bk:
      0.90   silver
% Cf:
      1.85   0.1    0.7    0.3
% Es:
      2.22   0.1    0.3    0.7 
% Fm:
      0.87   silver
% Md:
      0.86   silver
% No:
      0.85   silver
% Lr:
      0.84   silver
% Rf:
      0.84   silver
% Db:
      1.85   0.90   0.80   0.00
% Sg:
      0.84   silver
% Bh:
      0.84   silver
% Hs:
      0.84   silver
% Mt:
      0.84   silver 
];

if size(plotting_rrgb) ~= ZMAX
  fprintf (1, 'plotting_rrgb has %d elements\n', ...
        size(plotting_rrgb));
  return;
else
  fprintf (1, 'plotting_rrgb added ..\n');
end;
