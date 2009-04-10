%%%%%%%%%%%%%%%%%%%%%%%%%
% IUPAC_to_groupmates.m %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Group: the vertical columns (major classes or divisions) into which %
% elements are arranged in the periodic table of elements. There are three %
% common numbering systems for these groups:  %

% Groups 1-2 (except hydrogen) and 13-18 are termed main group elements. %

% Groups 3-11 are termed transition elements. Transition elements are those %
% whose atoms have an incomplete d-subshell or whose cations have an %
% incomplete d-subshell. %

% Main group elements in the first two rows of the table are called typical %
% elements. %

groupmates = cell(ZMAX,1);

for i = 1:ZMAX,
  switch IUPAC(i)
    case 1,
      groupmates{i} = 'H,Li,Na,K,Rb,Cs,Fr';
    case 2,
      groupmates{i} = 'Be,Mg,Ca,Sr,Ba,Ra';
   case 3,
    groupmates{i} = 'Sc,Y,La,Ac';
   case 3.5,
    groupmates{i} = [symbol{i} ',' symbol{i+32}];
   case 3.6,
    groupmates{i} = [symbol{i-32} ',' symbol{i}];
   case 4,
      groupmates{i} = 'Ti,Zr,Hf';
    case 5,
      groupmates{i} = 'V,Nb,Ta';
    case 6,
      groupmates{i} = 'Cr,Mo,W';
    case 7,
      groupmates{i} = 'Mn,Tc,Re';
    case 8,
      groupmates{i} = 'Fe,Ru,Os';
    case 9,
      groupmates{i} = 'Co,Rh,Ir';
    case 10,
      groupmates{i} = 'Ni,Pd,Pt';
    case 11,
      groupmates{i} = 'Cu,Ag,Au';
    case 12,
      groupmates{i} = 'Zn,Cd,Hg';
    case 13,
      groupmates{i} = 'B,Al,Ga,In,Tl';
    case 14,
      groupmates{i} = 'C,Si,Ge,Sn,Pb';
    case 15,
      groupmates{i} = 'N,P,As,Sb,Bi';
    case 16,
      groupmates{i} = 'O,S,Se,Te,Po';
    case 17,
      groupmates{i} = 'F,Cl,Br,I,At';
    case 18,
      groupmates{i} = 'He,Ne,Ar,Kr,Xe,Rn';
    otherwise,
      fprintf (1, 'group %f does not exist\n', IUPAC(i));
      break;
  end
end

fprintf (1, 'IUPAC converted to groupmates ..\n');
