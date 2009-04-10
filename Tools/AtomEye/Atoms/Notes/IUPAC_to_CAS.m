%%%%%%%%%%%%%%%%%%%%%%%%
% group_IUPAC_to_CAS.m %
%%%%%%%%%%%%%%%%%%%%%%%%

% CAS: Chicago Academy of Science, also known as Chemical %
% Abstracts Registry %

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

CAS = cell(ZMAX,1);

for i = 1:ZMAX,
  switch IUPAC(i)
    case 1,
      CAS{i} = 'IA';
    case 2,
      CAS{i} = 'IIA';
    case 3,
      CAS{i} = 'IIIB';
   case 3.5,
    CAS{i} = 'lanthanoids';
   case 3.6,
    CAS{i} = 'actanoids';
    case 4,
      CAS{i} = 'IVB';
    case 5,
      CAS{i} = 'VB';
    case 6,
      CAS{i} = 'VIB';
    case 7,
      CAS{i} = 'VIIB';
    case 8,
      CAS{i} = 'VIII';
    case 9,
      CAS{i} = 'VIII';
    case 10,
      CAS{i} = 'VIII';
    case 11,
      CAS{i} = 'IB';
    case 12,
      CAS{i} = 'IIB';
    case 13,
      CAS{i} = 'IIIA';
    case 14,
      CAS{i} = 'IVA';
    case 15,
      CAS{i} = 'VBA';
    case 16,
      CAS{i} = 'VIA';
    case 17,
      CAS{i} = 'VIIA';
    case 18,
      CAS{i} = 'VIIIA';
    otherwise, 
      fprintf (1, 'group %f does not exist\n', IUPAC(i));
      break;
  end
end

fprintf (1, 'IUPAC converted to CAS ..\n');
