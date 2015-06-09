#!/usr/bin/perl -w
# HJ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HJ X
# HJ X   libAtoms+QUIP: atomistic simulation library
# HJ X
# HJ X   Portions of this code were written by
# HJ X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
# HJ X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
# HJ X
# HJ X   Copyright 2006-2010.
# HJ X
# HJ X   These portions of the source code are released under the GNU General
# HJ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HJ X
# HJ X   If you would like to license the source code under different terms,
# HJ X   please contact Gabor Csanyi, gabor@csanyi.net
# HJ X
# HJ X   Portions of this code were written by Noam Bernstein as part of
# HJ X   his employment for the U.S. Government, and are not subject
# HJ X   to copyright in the USA.
# HJ X
# HJ X
# HJ X   When using this software, please cite the following reference:
# HJ X
# HJ X   http://www.libatoms.org
# HJ X
# HJ X  Additional contributions by
# HJ X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
# HJ X
# HJ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Perl script to generate quip.dtd Document Type Definition
# file, which can be used to validate an input XML file. 

#----------------------------------------------------------------------

# Attributes that match this regexp will be marked #REQUIRED rather than 
# #IMPLIED in the DTD.
my $required_atts = qr/(n_types)/;

# Elements that match this regexp can contain textual data
my $pcdata = qr/(point|comment|orb_set_type|E)/;

# Elements that don't get picked up in the automagic parsing
@extra_elements = qw(comment orb_set_type E);

# Describe nesting of tags - this is not picked up in the parsing code
my %subelements_add = (
		       qr/^per_pair_data$/ => ["point", "H_spline", "S_spline", "Vrep_spline"],
		       qr/^per_type_data$/  => ["orb_set_type", "E"],
		       qr/^H_spline$/ => ["point"],
		       qr/^S_spline$/ => ["point"],
		       qr/^Vrep_spline$/ => ["point"]
		       );

my %subelements_remove = (
			  qr/_params$/ => "point"
			  );


#----------------------------------------------------------------------

if ($#ARGV+1 <= 0) {
    die("\nUsage: mkdtd.pl [-t] IPModel*.f95 TBModel*.f95 > quip.dtd\n\nIf -t is present, document level element is included.\n\n");
}

if ($ARGV[0] =~ /-t/) {
    $print_toplevel = 1;
    shift(@ARGV);
} else {
    $print_toplevel = 0;
}

my $insub = 0;
my $inentity = 0;
my $inatts = 0;
my $gotname = 0;
my ($name, $topname);

my %elements = ();
my %validsubelements = ();

my @topnames = ();
my %subnames = ();

foreach $ex (@extra_elements) {
    $elements{$ex} = {};
    $subnames{$ex} = 1;
}

while (<>) {
    if (!$insub) {
	next until (/^\W*subroutine (\w*)_startElement_handler/);
	$insub = 1;
    } else {
	
	if (/^\W*end subroutine (\w*)_startElement_handler/) {
	    $gotname = 0;
	    $insub = 0;
	    $inentity = 0;
	    next
	}

	if (!$inentity) {
	    if (!$gotname) {
		next until (($name) = ($_ =~ /name == [\'\"](.*)[\'\"]/));
		$topname = $name;
		push @topnames, $topname;
		$validsubelements{$topname} = [];
		push @{$validsubelements{$topname}}, "comment";
	    }
	    $inentity = 1;
	    $elements{$name} = {} if (!defined($elements{$name}));

	} 

	if ($inentity) {
	    if (/name == [\'\"](.*)[\'\"]/) {
		$name = $1;
		$gotname = 1;
		$inentity = 0;
		if ($name ne $topname) { 
		    push (@{$validsubelements{$topname}},$name); 
		    $subnames{$name} = 1;
		    $validsubelements{$name} = [];
		}
		next;
	    }

	    next until ($att) = ($_ =~ /^\W*call QUIP_FoX_get_value\(attributes,\W*[\'\"]([a-zA-Z0-9_]*)[\'\"]/);
	    $elements{$name}{$att} = 1;
	}
    }
}

# Add sub elements required to describe nesting
while (($re, $subel) = each %subelements_add) {
    foreach $key (keys %validsubelements) {
	if ($key =~ $re) {
	    push @{$validsubelements{$key}}, @$subel;
	}
    }
}

# Remove subelements that shouldn't be include in top level elements
# (those that match /_params$/). 
while (($re, $subel) = each %subelements_remove) {
    foreach $key (keys %validsubelements) {
	if ($key =~ $re) {
	    @{$validsubelements{$key}} = grep (!/$subel/, @{$validsubelements{$key}});
	}
    }
}

# Print elements and entity definitions
foreach $el (@topnames, keys %subnames) {

    if (defined(${validsubelements{$el}}) && @{$validsubelements{$el}}) {
	$subelements = "(".join("|", @{$validsubelements{$el}}).")*";
    } elsif ($el =~ $pcdata) {
	$subelements = "(#PCDATA)";
    } else {
	$subelements = "EMPTY";
    }

    print "<!ELEMENT $el $subelements>\n";
    if (%{$elements{$el}}) {
	$attributes = join("\n  ", map({$_." CDATA ". ($_ =~ $required_atts ? "#REQUIRED" : "#IMPLIED")} (sort keys %{$elements{$el}})));
	print "<!ATTLIST $el\n  $attributes\n>\n\n";
    }
}


# Entity parameter
$joined = join("|", @topnames);
print qq|<!ENTITY % QUIP_params "$joined">\n\n|;

# Document level element
print "<!ELEMENT QUIP_params (%QUIP_params;)*>\n" if $print_toplevel;
