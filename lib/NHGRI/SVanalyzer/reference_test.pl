#!/usr/local/gtb/vendor/perlbrew/perls/perl-5.22.1/bin/perl

use strict;

my $ra_array = get_ref_to_string();
my $referent = $ra_array->[0] || "Nothing there!\n";

print $referent;

sub get_ref_to_string {

    my $string = "Hello there!\n";
    my @array = ();
    push @array, $string;

    return \@array;
}
