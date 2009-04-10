%%%%%%%%%%%%%%%%%%%
% charge_radius.m %
%%%%%%%%%%%%%%%%%%%

number;

% http://www.fhi-berlin.mpg.de/th/balsac/balm.47.html %
% For clusters and general lattices (codes 10,-10) values of atomic radii may %
% be determined as renormalized default radii according to respective nuclear %
% charges such that maximum packing without overlapping spheres is achieved, %
% see Secs. 6.2.1, 6.3.2. The internal table of default atomic radii 
% (given in Angstroms) is shown in the following. %
 
% [Angstrom] %
charge_radius = [
 0.4350
 1.4000
 1.5199
 1.1430
 0.9750
 0.6550
 0.7500
 0.7300 % O
 0.7200
 1.6000
 1.8579
 1.6047
 1.4318
 1.1758
 1.0600
 1.0200
 0.9900
 1.9000
 2.2620
 1.9758
 1.6545
 1.4755
 1.3090
 1.2490 % Cr
 1.3500
 1.2411
 1.2535
 1.2460 % Ni
 1.2780
 1.3325
 1.3501
 1.2248
 1.2000
 1.1600
 1.1400
 2.0000
 2.4700
 2.1513
 1.8237
 1.6156
 1.4318
 1.3626 % Mo
 1.3675
 1.3529
 1.3450
 1.3755
 1.4447 % Ag
 1.4894
 1.6662
 1.5375
 1.4000
 1.3600
 1.3300
 2.2000
 2.6325
 2.1705
 1.8725
 1.8243
 1.8362
 1.8295
 1.8090
 1.8040
 1.9840
 1.8180
 1.8005
 1.7951
 1.7886
 1.7794
 1.7687
 1.9396
 1.7515
 1.5973
 1.4280
 1.3705
 1.3800
 1.3676
 1.3573
 1.3873
 1.4419
 1.5025
 1.7283
 1.7501
 1.4600  % Bi
 1.4600
 4.3500  % At
 1.4300
 2.5000
 2.1400
 1.8775
 1.7975
 1.6086
 1.5683
 1.0000 
 1.0000 
 1.0000 
 1.0000 
 1.0000 
 1.4600  % Cf
 1.752   % Es
 1.0000 
 1.0000 
 1.0000 
 1.0000 
 1.0000 
 1.4600 % Db
 1.0000 
 1.0000 
 1.0000 
 1.0000 
];
% [Angstrom] %

if size(charge_radius) ~= ZMAX
  fprintf (1, 'charge_radius has %d elements\n', ...
      size(charge_radius));
  return;
else
  fprintf (1, 'charge_radius added ..\n');
end;
