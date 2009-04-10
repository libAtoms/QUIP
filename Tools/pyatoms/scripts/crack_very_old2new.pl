#!/usr/bin/perl

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
