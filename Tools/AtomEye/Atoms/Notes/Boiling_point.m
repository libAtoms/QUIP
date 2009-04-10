%%%%%%%%%%%%%%%%%%%
% boiling_point.m %
%%%%%%%%%%%%%%%%%%%

number;

% The boiling point of an element is the temperature at which there is the %
% transition from the liquid to the gas phase. The element with the lowest %
% boiling point is Helium (4.216 K). Tungsten has the highest boiling point %
% (6200 K). %

% [K] %
boiling_point = [
20.28 
4.216 
1590 
3243 
2823 
5100 
77.4 
90.188 
85.01 
27.1 
1165 
1380 
2740 
2628 
553 
717.824 
238.55 
87.29 
1047 
1760 
3105 
3533 
3653 
2755 
2370 
3023 
3143 
3005 
2868 
1180 
2676 
3103 
NVL
958.1 
331.93 
120.85 
961 
1657 
3610 
4650 
5200 
5833 
5303 
4173 
4000 
3413 
2485 
1038 
2353 
2543 
2023 
1263 
457.55 
166.1 
963 
1913 
3727 
3530 
3485 
3400 
3000 
2051 
1870 
3506 
3314 
2608 
2993 
2783 
2000 
1466 
3588 
5673 
5698 
6200 
5900 
5300 
4403 
4100 
3213 
629.73 
1730 
2013 
1833 
1235 
610 
211.4 
950 
1413 
3470 
5060 
4300 
4091 
4175 
3600 
2880 
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
% [K] %

if size(boiling_point) ~= ZMAX
  fprintf (1, 'boiling_point has %d elements\n', ...
      size(boiling_point));
  return;
else
  fprintf (1, 'boiling_point added ..\n');
end;
