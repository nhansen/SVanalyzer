#!/usr/local/gtb/vendor/perlbrew/perls/perl-5.22.1/bin/perl -w
# $Id:$

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);

our %Opt;

=head1 NAME

SVcluster.pl - cluster structural variants based on distance metrics printed by SVmerge.pl

=head1 SYNOPSIS

  SVcluster.pl --dist <SVmerge distance file> --ids <optional file of ids to allow output of single "clusters">

=head1 DESCRIPTION

The program steps through a distance file, creating an input file for a simple clustering 
program based on user-specified criteria.

=cut

#------------
# Begin MAIN 
#------------

$|=1;

process_commandline();

my $workingdir = $Opt{workdir}; # good to allow use of a temporary file

my $dist_file = $Opt{dist};
my $id_file = $Opt{ids};
my $max_posdiff = $Opt{posdiff};
my $max_relshift = $Opt{relshift};
my $max_relsizediff = $Opt{relsizediff};
my $max_reldist = $Opt{reldist};

my $dist_fh = Open($dist_file);

my %node_index = (); # stores the index of each node id
my @nodes = (); # array of nodes

if ($id_file) { # read in list of node ids and assign indices
    my $ids_fh = Open($id_file);
    while (<$ids_fh>) {
        if (/^(\S+)$/) {
            push @nodes, $1;
            $node_index{$1} = $#nodes + 1;
        }
        else {
            die "Invalid line in id file $id_file--lines must contain only one entry with no spaces!\n";
        }
    }
    close $ids_fh;
}

my @node_edges = (); # array of nodes

while (<$dist_fh>) {

    chomp;
    my ($id1, $id2, $posdiff, $relshift, $relsizediff, $reldist) = split /\t/, $_;
    next if ($id1 eq 'ID1');

    if (!(defined($reldist))) {
        print "Skipping line--not enough fields:\n$_\n";
        next;
    }
    if ($posdiff <= $max_posdiff && $relshift <= $max_relshift &&
        $relsizediff <= $max_relsizediff && $reldist <= $max_reldist) {
        my $index1 = $node_index{$id1};
        if (!(defined($index1))) {
            push @nodes, $id1;
            $index1 = $#nodes + 1;
            $node_index{$id1} = $index1;
        }
        my $index2 = $node_index{$id2};
        if (!(defined($index2))) {
            push @nodes, $id2;
            $index2 = $#nodes + 1;
            $node_index{$id2} = $index2;
        }
        push @{$node_edges[$index1 - 1]}, $index2 - 1;
        push @{$node_edges[$index2 - 1]}, $index1 - 1;
    }
}

close $dist_fh;

# Now do DFS to find connected components:

my @visited = ();
for (my $i=0; $i<=$#nodes; $i++) {
    $visited[$i] = 0;
}

for (my $i=0; $i<=$#nodes; $i++) {
    next if ($visited[$i]);
    my @connectednodes = dfs($i, \@node_edges, \@visited);
    my @nodenames = map { $nodes[$_] } @connectednodes;
    my $nodestring = join ':', @nodenames;
    print "$nodestring\n";
}


#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( workdir => '.', 'posdiff' => 100000, 'relshift' => 1.0, 'relsizediff' => 1.0, 'reldist' => 1.0 );
    GetOptions(\%Opt, qw( dist=s workdir=s ids=s posdiff=i relshift=f relsizediff=f reldist=f cleanup manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "SVcluster.pl, ", q$Revision:$, "\n"; }

    if ($Opt{workdir} ne '.') {
        mkdir $Opt{workdir}; # don't stress if it was already there, so not checking return value
    }

}

sub dfs {
    my $i = shift;
    my $ra_nodeedges = shift;
    my $ra_visited = shift;

    my @stack = ();
    push @stack, $i;

    my @connectednodes = ();
    while (@stack) {
        my $j = pop @stack;
        if (!$ra_visited->[$j]) {
            $ra_visited->[$j] = 1;
            push @connectednodes, $j;
            foreach my $k (@{$ra_nodeedges->[$j]}) {
                next if ($ra_visited->[$k]);
                push @stack, $k;
            }
        }
    }

    return @connectednodes;
}

__END__

=head1 OPTIONS

=over 4

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=back

=head1 AUTHOR

 Nancy Hansen - nhansen@mail.nih.gov

=head1 LEGAL

This software/database is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software/database is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use. 

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S.  Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose. 

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
Genome Technology Branch" as the citation. 

=cut
