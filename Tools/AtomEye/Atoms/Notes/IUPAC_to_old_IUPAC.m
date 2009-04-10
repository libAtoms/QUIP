%%%%%%%%%%%%%%%%%%%%%%%%
% IUPAC_to_old_IUPAC.m %
%%%%%%%%%%%%%%%%%%%%%%%%


% http://klbproductions.com/yogi/periodic/ %

% Group: the vertical columns (major classes or divisions) into which %
% elements are arranged in the periodic table of elements. There are three %
% common numbering systems for these groups:  %

% The new IUPAC system, which numbers each column with Arabic %
% numbers from 1 (one) through 18 (eighteen). To reduce confusion %
% caused by the other two systems, this is the system that is used in %
% articles on this web site.  %

% The old IUPAC system, which labeled columns with Roman numerals %
% followed by either the letter 'A' or 'B'. Columns were numbed such that %
% columns one through seven were numbered 'IA' through 'VIIA', %
% columns 8 through 10 were labeled 'VIIIA', columns 11 through 17 %
% were numbered 'IB' through 'VIIB' and column 18 was numbered %
% 'VIII'.  %

% The CAS system, which also used Roman numerals followed by an 'A' %
% or 'B'. This method, however, labeled columns 1 and 2 as 'IA' and 'IIA', %
% columns 3 through 7 as 'IIIB' through 'VIB', column 8 through 10 as %
% 'VIII', columns 11 and 12 as 'IB' and 'IIB' and columns 13 through 18 as %
% 'IIIA' through 'VIIIA'.  %

% Because of the confusion the old IUPAC and the CAS system created, %
% the IUPAC adopted their new system.  %

% Elements are arranged in these groups according to whose proprieties %
% are similar. All elements in Group 1 for instance are alkali metals. They %
% have only one electron in the outer shell (valence electron) and as a %
% result are highly reactive. Elements in Group 17 are the halogens. They %
% all have seven electrons in the outer orbital (two in level s and five in %
% level p). They are also very reactive because they have seven electrons %
% in the outer shell and will readily accept an electron in order to reach
% the ion configuration with the ideal number of eight electrons in the outer %
% shell. Elements Group 18 have a complete outer shell with eight %
% electrons. These noble gases are highly stable and do not react to form %
% compounds under normal conditions. %

old_IUPAC = cell(ZMAX,1);

for i = 1:ZMAX,
  switch IUPAC(i)
    case 1,
      old_IUPAC{i} = 'IA';
    case 2,
      old_IUPAC{i} = 'IIA';
    case 3,
      old_IUPAC{i} = 'IIIA';
   case 3.5,
    old_IUPAC{i} = 'lanthanoids';
   case 3.6,
    old_IUPAC{i} = 'actanoids';
   case 4,
      old_IUPAC{i} = 'IVA';
    case 5,
      old_IUPAC{i} = 'VA';
    case 6,
      old_IUPAC{i} = 'VIA';
    case 7,
      old_IUPAC{i} = 'VIIA';
    case 8,
      old_IUPAC{i} = 'VIIIA';
    case 9,
      old_IUPAC{i} = 'VIIIA';
    case 10,
      old_IUPAC{i} = 'VIIIA';
    case 11,
      old_IUPAC{i} = 'IB';
    case 12,
      old_IUPAC{i} = 'IIB';
    case 13,
      old_IUPAC{i} = 'IIIB';
    case 14,
      old_IUPAC{i} = 'IVB';
    case 15,
      old_IUPAC{i} = 'VB';
    case 16,
      old_IUPAC{i} = 'VIB';
    case 17,
      old_IUPAC{i} = 'VIIB';
    case 18,
      old_IUPAC{i} = 'VIIIB';
    otherwise, 
      fprintf (1, 'group %f does not exist\n', IUPAC(i));
      break;
  end
end

fprintf (1, 'IUPAC converted to old_IUPAC ..\n');
