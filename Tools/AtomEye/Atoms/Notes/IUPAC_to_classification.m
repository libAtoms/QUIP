%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IUPAC_to_classification.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Group: the vertical columns (major classes or divisions) into which %
% elements are arranged in the periodic table of elements. There are three %
% common numbering systems for these groups:  %

% Groups 1-2 (except hydrogen) and 13-18 are termed main group elements. %

% Groups 3-11 are termed transition elements. Transition elements are those %
% whose atoms have an incomplete d-subshell or whose cations have an %
% incomplete d-subshell. %

% Main group elements in the first two rows of the table are called typical %
% elements. %

% http://www.mcs.net/~ars/spectro/elements.htm %

classification = cell(ZMAX,1);

for i = 1:ZMAX,
  switch IUPAC(i)
   case 1,
    if (i==1) classification{i} = 'hydrogen group non-metal';
    else classification{i} = 'hydrogen group alkali metal'; end;
   case 2,
    classification{i} = 'alkali earth';
   case 3,
    classification{i} = 'scandium group transition metal';
   case 3.5,
    classification{i} = 'lanthanide series rare earth';
   case 3.6,
    classification{i} = 'actinide series rare earth';
   case 4,
    classification{i} = 'titanium group transition metal';
   case 5,
    classification{i} = 'vanadium group transition metal';
   case 6,
    classification{i} = 'chromium group transition metal';
   case 7,
    classification{i} = 'manganese group transition metal';
   case 8,
    classification{i} = 'iron group transition metal';
   case 9,
    classification{i} = 'cobalt group transition metal';
   case 10,
    classification{i} = 'nickel group transition metal';
   case 11,
    classification{i} = 'copper group transition metal';
   case 12,
    classification{i} = 'zinc group transition metal';
   case 13,
    if (i==5) classification{i} = 'boron group non-metal';
    else classification{i} = 'boron group metal'; end;
   case 14,
    if (i<=14) classification{i} = 'carbon group non-metal';
    else classification{i} = 'carbon group metal'; end;
   case 15,
    if (i<=33) classification{i} = ['nitrogen group (pnictogen) non-metal'];
    else classification{i} = 'nitrogen group (pnictogen) metal'; end;
   case 16,
    if (i<=52) classification{i} = 'oxygen group (chalcogen) non-metal';
    else classification{i} = 'oxygen group (chalcogen) metal'; end;
   case 17,
    classification{i} = 'halogen';
   case 18,
    classification{i} = 'noble gas';
   otherwise,
    fprintf (1, 'group %f does not exist\n', IUPAC(i));
    break;
  end
end

fprintf (1, 'IUPAC converted to classification ..\n');
