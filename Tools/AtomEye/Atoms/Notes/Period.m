%%%%%%%%%%%%
% Period.m %
%%%%%%%%%%%%

number;

% Periods: While groups are characterized by the number of electrons present %
% in the outer shell, periods are characterized by the number of energy %
% levels (shells) of electrons surrounding the nucleus. Elements in Period %
% 1 has only one shell. As you probably recall, the elements in the first %
% period have a 2 electrons maximum (hydrogen has 1 electron and helium has %
% 2 electrons. As we move to the first group of the second period, we find %
% that lithium, which has the two electrons in the first shell and one in %
% the second. Neon is in Group 18 of Period 2 and therefore has the two %
% electrons in the first shell and eight electrons in the second shell. %
% Sodium starts Period 3 with 11 electrons, two in the first shell, eight %
% in the second shell and one in the third shell. In other words, the %
% element in Group 1 always has one more electron (in a new shell) than the %
% Group 18 element in the previous period. %

period = [
1
1

2
2
2
2
2
2
2
2

3
3
3
3
3
3
3
3

4
4
4
4
4
4
4
4
4
4
4
4
4
4
4
4
4
4

5
5
5
5
5
5
5
5
5
5
5
5
5
5
5
5
5
5

6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6

7
7
7
7
7
7
7
7
7
7
7
7
7
7
7
7
7
7
7
7
7
7
7
];

if size(period) ~= ZMAX then
  fprintf (1, 'period has %d elements\n', size(period));
  break;
else
  fprintf (1, 'period added ..\n');
end;
