#!/usr/bin/perl
# HP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HP X
# HP X   pyatoms: atomistic simulations tools
# HP X
# HP X   Copyright James Kermode 2010
# HP X
# HP X   These portions of the source code are released under the GNU General
# HP X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HP X
# HP X   If you would like to license the source code under different terms,
# HP X   please contact James Kermode, james.kermode@gmail.com
# HP X
# HP X   When using this software, please cite the following reference:
# HP X
# HP X   http://www.jrkermode.co.uk/PyAtoms
# HP X
# HP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

while (<>) {
    if (/AuxProps=([^ \n]*)/) {
	@fields = split(/:/,$1);
	$props = "Properties=species:S:1:pos:R:3:";
	foreach $f (@fields) { $props .= $f.":R:1:"; }
	chop($props);
	s/AuxProps=([^ \n]*)/$props/;
	print;
    }
    else {
	print;
    }
}
