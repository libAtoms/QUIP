%%%%%%%%%%%
% runme.m %
%%%%%%%%%%%

clear;
more off;

fname = '../Mendeleyev.c';

fprintf (1, '\n');

Symbol;
Name;

fprintf (1, '\n');

ZZ;
AA;
Period;
Iupac;

fprintf (1, '\n');

IUPAC_to_old_IUPAC;
IUPAC_to_CAS;
IUPAC_to_classification;
IUPAC_to_groupmates;

fprintf (1, '\n');

Electron_config;
First_ionization_energy;

fprintf (1, '\n');

Electronegativity_Allred;
Electronegativity_Pauling;
Electronegativity_Pearson;

fprintf (1, '\n');

Aggregate_state;
Mass_density;
Melting_point;
Boiling_point;

fprintf (1, '\n');

Xtal;

fprintf (1, '\n');

Empirical_radius; 
Charge_radius;
Plotting_rrgb;

fprintf (1, '\nwriting to file ...\n');

fid = fopen (fname, 'w');

for i = 0 : ZMAX,
  
  fprintf (fid, '{');
  
  Unknown = 0;
  if i == 0,
    fprintf (fid, '" A",  ');
    fprintf (fid, '"Unknown",  ');
    Unknown = 1;
    i = 1;  
    % hydrogen character %
  else
    symboll = [' ' symbol{i}];
    j = size(symboll,2);
    fprintf (fid, '"%2s",  ', symboll(j-1:j));
    fprintf (fid, '"%s",  ', name{i});
  end;
  
  fprintf (fid, '\n ');
  
  fprintf (fid, '%d,  ', Z(i));
  fprintf (fid, '%.10g,  ', A(i));
  fprintf (fid, '%d,  ', period(i));
  if IUPAC(i) ~= 3.5
    fprintf (fid, '%d,  ', IUPAC(i));
  else
    fprintf (fid, 'NVL, ');
  end;
  
  fprintf (fid, '\n ');
  
  fprintf (fid, '"%s",  ', old_IUPAC{i});
  fprintf (fid, '"%s",  ', CAS{i});
  fprintf (fid, '"%s",  ', classification{i});
  fprintf (fid, '"%s",  ', groupmates{i});
  
  fprintf (fid, '\n ');
  
  fprintf (fid, '"%s",  ', electron_config{i});
  if first_ionization_energy(i) ~= NVL
    fprintf (fid, '%.10g,  ', first_ionization_energy(i));
  else
    fprintf (fid, 'NVL, ');
  end;
  
  fprintf (fid, '\n ');
  
  if electronegativity_Allred(i) ~= NVL
    fprintf (fid, '%.10g,  ', electronegativity_Allred(i));
  else
    fprintf (fid, 'NVL, ');
  end;
  if electronegativity_Pauling(i) ~= NVL
    fprintf (fid, '%.10g,  ', electronegativity_Pauling(i));
  else
    fprintf (fid, 'NVL, ');
  end;
  if electronegativity_Pearson(i) ~= NVL
    fprintf (fid, '%.10g,  ', electronegativity_Pearson(i));
  else
    fprintf (fid, 'NVL, ');
  end;
  
  fprintf (fid, '\n ');
  
  fprintf (fid, '"%s",  ', aggregate_state{i});
  if mass_density(i) ~= NVL
    fprintf (fid, '%.10g,  ', mass_density(i));
  else
    fprintf (fid, 'NVL, ');
  end;
  if melting_point(i) ~= NVL
    fprintf (fid, '%.10g,  ', melting_point(i));
  else
    fprintf (fid, 'NVL, ');
  end;
  if boiling_point(i) ~= NVL
    fprintf (fid, '%.10g,  ', boiling_point(i));
  else
    fprintf (fid, 'NVL, ');
  end;
  
  fprintf (fid, '\n ');
  
  fprintf (fid, '"%s",  ', xtal{i,1});
  if xtal{i,2} ~= NVL
    fprintf (fid, '%.10g,  ', xtal{i,2});
  else
    fprintf (fid, 'NVL,  ');
  end;
  fprintf (fid, '"%s",  ', xtal{i,3});
  if xtal{i,4} ~= NVL
    fprintf (fid, '%.10g, ', xtal{i,4}/100);
  else
    fprintf (fid, 'NVL, ');
  end;
  if xtal{i,5} ~= NVL
    fprintf (fid, '%.10g, ', xtal{i,5}/100);
  else
    fprintf (fid, 'NVL, ');
  end;
  if xtal{i,6} ~= NVL
    fprintf (fid, '%.10g,  ', xtal{i,6}/100);
  else
    fprintf (fid, 'NVL,  ');
  end;
  if xtal{i,7} ~= NVL
    fprintf (fid, '%.10g, ', xtal{i,7});
  else
    fprintf (fid, 'NVL,');
  end;
  if xtal{i,8} ~= NVL
    fprintf (fid, '%.10g, ', xtal{i,8});
  else
    fprintf (fid, 'NVL,');
  end;
  if xtal{i,9} ~= NVL
    fprintf (fid, '%.10g, ', xtal{i,9});
  else
    fprintf (fid, 'NVL, ');
  end;
  
  fprintf (fid, '\n ');
  
  if empirical_radius(i) ~= NVL
    fprintf (fid, '%.10g,  ', empirical_radius(i));
  else
    fprintf (fid, 'NVL,  ');
  end;
  fprintf (fid, '%.10g,  ', charge_radius(i));
  
  if Unknown == 1,  
%  dark hydrogen %
   fprintf (fid, '%.10g, %.10g, %.10g', 0.5, 0.5, 0.5);
  else
    fprintf (fid, '%.10g, %.10g, %.10g', ...
             plotting_rrgb(i,2), plotting_rrgb(i,3), plotting_rrgb(i,4));
  end;
  
  fprintf (fid, '},\n\n');
end;

fprintf (1, 'saved on file "%s"\n\n', fname);
