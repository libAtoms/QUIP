%%%%%%%%
% ZZ.m %
%%%%%%%%

number;

% Atomic charge

Z = (1:ZMAX)';

if size(Z) ~= ZMAX then
  fprintf (1, 'Z has %d elements\n', size(Z));
  break;
else
  fprintf (1, 'Z added ..\n');
end;
