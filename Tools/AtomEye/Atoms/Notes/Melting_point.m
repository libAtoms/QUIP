%%%%%%%%%%%%%%%%%%%
% melting_point.m %
%%%%%%%%%%%%%%%%%%%

number;

% The melting point of an element is the %
% temperature at which the transition from the %
% solid phase to the liquid phase occurs. The %
% element with the lowest melting point is %
% Helium (0,95 K). Arsenic sublimates at a %
% temperature of about 886 K, i.e. makes a %
% direct transition from the solid to the gas %
% phase. Carbon has the highest melting point %
% (3823 K). %
 
% [K] %
melting_point = [
14.01 
0.95 
453.69 
1551 
2573 
3823 
63.29 
54.75 
53.53 
24.48 
370.95 
921.95 
933.52 
1683 
317.3 
386 
172.17 
83.78 
336.8 
1112 
1812 
1933 
2163 
2130 
1517 
1808 
1768 
1726 
1356.6 
692.73 
302.93 
1210.55
886 
490 
265.9 
116.55 
312.2 
1042 
1796 
2125 
2741 
2890 
2445 
2583 
2239 
1825 
1235.08
594.1 
429.32 
505.118
903.89 
722.7 
386.65 
161.3 
301.55 
998 
1193 
1071 
1204 
1283 
1353 
1345 
1095 
1584 
1633 
1682 
1743 
1795 
1818 
1097 
1929 
2423 
3269 
3680 
3453 
3318 
2683 
2045 
1337.58
234.28 
576.7 
600.65 
544.5 
527 
575 
202 
300 
973 
1320 
2023 
1827 
1405.5 
913 
914 
1267 
1613 
1259 
1173 
1133 
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

if size(melting_point) ~= ZMAX
  fprintf (1, 'melting_point has %d elements\n', ...
      size(melting_point));
  return;
else
  fprintf (1, 'melting_point added ..\n');
end;
