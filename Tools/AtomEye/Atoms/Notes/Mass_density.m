%%%%%%%%%%%%%%%%%%
% mass_density.m %
%%%%%%%%%%%%%%%%%%

number;

% The density is given by the ratio of mass and volume. The density is a %
% material constant that depends on the temperature and the pressure,  %
% especially for gases. Thus the temperature and the pressure have to be  %
% given along with the density. %

% Note, that in this table the densities for solids are given in  %
% units of g/cm3 at a temperature of 20 C and 1013 hPa. Hydrogen,  %
% Helium, and Lithium have the lowest densities, while the noble %
% metals Iridium, Osmium, and Platinum have the highest densities.  %
 
% [g/cm3] %
mass_density = [
0.084e-3
0.17e-3
0.53
1.85
2.46
3.51
1.17e-3
1.33e-3
1.58e-3
0.84e-3
0.97
1.74
2.70
2.33
1.82
2.06
2.95e-3
1.66e-3
0.86
1.54
2.99
4.51
6.09
7.14
7.44
7.87
8.89
8.91
8.92
7.14
5.91
5.32
5.72
4.82
3.14
3.48e-3
1.53
2.63
4.47
6.51
8.58
10.28
11.49
12.45
12.41
12.02
10.49
8.64
7.31
7.29
6.69
6.25
4.94
4.49e-3
1.90
3.65
6.16
6.77
6.48
7.00
7.22
7.54
5.25
7.89
8.25
8.56
8.78
9.05
9.32
6.97
9.84
13.31
16.68
19.26
21.03
22.61
22.65
21.45
19.32
13.55
11.85
11.34
9.80
9.20
NVL 
9.23e-3
NVL 
5.50
10.07
11.72
15.37
18.97
20.48
19.74
13.67
13.51
13.25
15.1
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
% [g/cm3] %

if size(mass_density) ~= ZMAX
  fprintf (1, 'mass_density has %d elements\n', ...
      size(mass_density));
  return;
else
  fprintf (1, 'mass_density added ..\n');
end;
