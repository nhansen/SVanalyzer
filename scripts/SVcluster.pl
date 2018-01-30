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
my $rh_distances = {}; # hash of distance arrays for pairs of ids

while (<$dist_fh>) {

    chomp;
    s/^DIST\s//; # some versions of distance file have "DIST\t" at beginning of each line
    my $line = $_;
    my @fields = split /\t/, $line; # first two fields must be ids, last four dist measures
    my $id1 = $fields[0];
    my $id2 = $fields[1];
    my $posdiff = $fields[$#fields-3];
    my $relshift = $fields[$#fields-2];
    my $relsizediff = $fields[$#fields-1];
    my $reldist = $fields[$#fields];
    next if ($id1 eq 'ID1');
    next if ($id1 =~ /Read \d+ SVs/);
    #next if (!$id1);

    if (!(defined($reldist))) {
        print "Skipping line--not enough fields:\n$line\n";
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
    my $idlesser = (($id1 cmp $id2) > 0) ? $id2 : $id1;
    my $idgreater = (($id1 cmp $id2) > 0) ? $id1 : $id2;
    $rh_distances->{"$idlesser:$idgreater"} = [$posdiff, $relshift, $relsizediff, $reldist];
}

close $dist_fh;

# Now do DFS to find connected components:

my @visited = ();
for (my $i=0; $i<=$#nodes; $i++) {
    $visited[$i] = 0;
}

my $rh_cluster_info = {};
for (my $i=0; $i<=$#nodes; $i++) {
    next if ($visited[$i]);
    my @connectednodes = dfs($i, \@node_edges, \@visited);
    my @nodenames = map { $nodes[$_] } @connectednodes;
    my $nodestring = join ':', @nodenames;
    my ($clustercall, $no_cluster_calls, $dummynodestring, $no_exact_matches, $exactsubcluster, $maxdist1, $maxdist2, $maxdist3) = 
                 analyze_cluster($nodestring, $rh_distances);
    print "$nodestring\t$clustercall\t$no_cluster_calls\t$nodestring\t$no_exact_matches\t$exactsubcluster\t$maxdist1\t$maxdist2\t$maxdist3\n";
    $rh_cluster_info->{$clustercall} = {'nodestring' => $nodestring, 'exactstring' => $exactsubcluster, 'maxdists' => [$maxdist1, $maxdist2, $maxdist3]};
}

if ($Opt{'vcf'}) { # if user specified a VCF file, write out a "clustered" VCF file with only annotated representative variants
    my $vcf = $Opt{'vcf'};
    my $newvcf = $vcf;
    $newvcf =~ s/\.vcf/.clustered.$max_relshift.$max_relsizediff.$max_reldist.vcf/;

    if ($vcf eq $newvcf) {
        print STDERR "VCF file name passed with --vcf option must contain the string \".vcf\". Skipping VCF file creation.\n";
        exit;
    }

    my $vcf_fh = Open($vcf);
    my $newvcf_fh = Open($newvcf, "w");

    my $ra_allowed_info_fields = [];
    while (<$vcf_fh>) {
        if (/^##/) {
            print $newvcf_fh $_; # include original VCF header in new file

            if (/^##INFO=<ID=([^,]+)/) {
                push @{$ra_allowed_info_fields}, $1;
            }
        }
        elsif (/^#CHROM/) { # include new INFO lines
            my $chromline = $_;

            # print lots of INFO lines
            print $newvcf_fh "##INFO=<ID=ClusterIDs,Number=1,Type=String,Description=\"IDs of SVs that cluster with this SV\">\n";
            print $newvcf_fh "##INFO=<ID=NumClusterSVs,Number=1,Type=Integer,Description=\"Total number of SVs in this cluster\">\n";
            print $newvcf_fh "##INFO=<ID=ExactMatchIDs,Number=1,Type=String,Description=\"IDs of SVs that are exactly the same call as this SV\">\n";
            print $newvcf_fh "##INFO=<ID=NumExactMatchSVs,Number=1,Type=Integer,Description=\"Total number of SVs in this exact cluster\">\n";
            print $newvcf_fh "##INFO=<ID=ClusterMaxShiftDist,Number=1,Type=Float,Description=\"Maximum relative shift distance between two SVs in this cluster\">\n";
            print $newvcf_fh "##INFO=<ID=ClusterMaxSizeDiff,Number=1,Type=Float,Description=\"Maximum relative size difference between two SVs in this cluster\">\n";
            print $newvcf_fh "##INFO=<ID=ClusterMaxEditDist,Number=1,Type=Float,Description=\"Maximum relative edit distance between two SVs in this cluster\">\n";
            print $newvcf_fh $chromline; # include original VCF header in new file
        }
        else {
            chomp;
            my @fields = split /\t/, $_;
            my $id_field = $fields[2];

            if ($rh_cluster_info->{$id_field}) { # this is a cluster rep!
                write_vcf_line_with_info($newvcf_fh, \@fields, $rh_cluster_info->{$id_field}, $ra_allowed_info_fields);
            }
        }
    }
}

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( workdir => '.', 'posdiff' => 100000, 'relshift' => 1.0, 'relsizediff' => 1.0, 'reldist' => 1.0 );
    GetOptions(\%Opt, qw( dist=s workdir=s ids=s posdiff=i relshift=f relsizediff=f reldist=f vcf=s cleanup manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "SVcluster.pl, ", q$Revision:$, "\n"; }

    if ($Opt{workdir} ne '.') {
        mkdir $Opt{workdir}; # don't stress if it was already there, so not checking return value
    }

    if (!$Opt{dist}) {
        die "SVcluster.pl requires distance file as input with --dist parameter!\n";
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

sub analyze_cluster {
    my $nodestring = shift;
    my $rh_dists = shift;

    my @nodes = split /:/, $nodestring;
    my $no_cluster_calls = @nodes;

    if ($no_cluster_calls == 1) {
        return ($nodestring, 1, $nodestring, 1, $nodestring, 0, 0, 0);
    }

    # at least two calls--analyze exact subclusters
    my ($maxdist1, $maxdist2, $maxdist3) = (0, 0, 0);
    my @exact_clusters = @nodes;
    my %node_cluster = map { $exact_clusters[$_] => $_ } ( 0 .. $#nodes ); # all start out in their own cluster
    for (my $i=0; $i<=$#nodes; $i++) {
        for (my $j=$i + 1; $j<=$#nodes; $j++) {
            my $idlesser = (($nodes[$i] cmp $nodes[$j]) > 0) ? $nodes[$j] : $nodes[$i];
            my $idgreater = (($nodes[$i] cmp $nodes[$j]) > 0) ? $nodes[$i] : $nodes[$j];
            if (!$rh_dists->{"$idlesser:$idgreater"}) {
                #print STDERR "Missing distance information for $idlesser:$idgreater!\n";
                next;
            }
            my ($thisposdist, $this_dist1, $this_dist2, $this_dist3) = @{$rh_dists->{"$idlesser:$idgreater"}};
            $maxdist1 = ($this_dist1 > $maxdist1) ? $this_dist1 : $maxdist1;
            $maxdist2 = ($this_dist2 > $maxdist2) ? $this_dist2 : $maxdist2;
            $maxdist3 = ($this_dist3 > $maxdist3) ? $this_dist3 : $maxdist3;
            if (!$this_dist1 && !$this_dist2 && !$this_dist3) { # exact match--merge exact clusters
                if ($node_cluster{$idlesser} != $node_cluster{$idgreater}) { # different clusters merge greater into lesser
                    #print "Merging $idgreater ($node_cluster{$idgreater}) into $idlesser ($node_cluster{$idlesser})!\n";
                    if (!($exact_clusters[$node_cluster{$idgreater}])) {
                        print "Empty node cluster for $idgreater being merged into $idlesser!\n";
                    }
                    my $old_cluster_index = $node_cluster{$idgreater};
                    my $new_cluster_index = $node_cluster{$idlesser};
                    $exact_clusters[$new_cluster_index] .= ":".$exact_clusters[$old_cluster_index];
                    $exact_clusters[$old_cluster_index] = '';
                    foreach my $key (keys %node_cluster) {
                        if ($node_cluster{$key} == $old_cluster_index) {
                            $node_cluster{$key} = $new_cluster_index;
                        }
                    }
                }
                else {
                    #print "$idgreater and $idlesser are already in the same exact cluster!\n";
                } 
            }

            my $clusterstring = join '/', @exact_clusters;
            #print "i$i j$j $nodes[$i] $nodes[$j] $clusterstring\n";
        }
    }

    my @unique_clusters = grep { $_ ne '' } @exact_clusters;
    my @largest_exact_clusters = ();
    my $max_exact_clustersize = 0;

    foreach my $unique_cluster (@unique_clusters) {
        my $no_vars = split /:/, $unique_cluster;
        if ($no_vars == $max_exact_clustersize) { # it's a tie
            push @largest_exact_clusters, $unique_cluster;
            $max_exact_clustersize = $no_vars;
        }
        elsif ($no_vars > $max_exact_clustersize) { # replace
            @largest_exact_clusters = ($unique_cluster);
            $max_exact_clustersize = $no_vars;
        }
    }

    my $exactsubcluster = $largest_exact_clusters[rand()*$#largest_exact_clusters];
    my @exact_vars = split /:/, $exactsubcluster;
    my $clustercall = $exact_vars[rand()*$#exact_vars];

    return ($clustercall, $no_cluster_calls, $nodestring, $max_exact_clustersize, $exactsubcluster, $maxdist1, $maxdist2, $maxdist3);
}

sub write_vcf_line_with_info {
    my $vcf_fh = shift;
    my $ra_vcf_fields = shift;
    my $rh_cluster_info = shift;
    my $ra_allowed_info_fields = shift;

    my $clustersvcalls = $rh_cluster_info->{'nodestring'};
    my $no_clustersvcalls = split /:/, $clustersvcalls;
    my $exactsvcalls = $rh_cluster_info->{'exactstring'};
    my $no_exactsvcalls = split /:/, $exactsvcalls;
    my $ra_maxdists = $rh_cluster_info->{'maxdists'};

    my $info_field = $ra_vcf_fields->[7];
    my @info_entries = ($info_field eq '.') ? () : split /;/, $info_field;
    my @allowed_info_entries = ();
    foreach my $info_entry (@info_entries) {
        my $info_id = $info_entry;
        $info_id =~ s/=.*$//;
        if (grep {$_ eq $info_id} @{$ra_allowed_info_fields}) {
            push @allowed_info_entries, $info_entry;
        }
    }
    @info_entries = @allowed_info_entries;

    my ($maxshiftdist, $maxsizediff, $maxeditdist) = @{$ra_maxdists};

    push @info_entries, "ClusterIDs=$clustersvcalls";
    push @info_entries, "NumClusterSVs=$no_clustersvcalls";
    push @info_entries, "ExactMatchIDs=$exactsvcalls";
    push @info_entries, "NumExactMatchSVs=$no_exactsvcalls";
    push @info_entries, "ClusterMaxShiftDist=$maxshiftdist";
    push @info_entries, "ClusterMaxSizeDiff=$maxsizediff";
    push @info_entries, "ClusterMaxEditDist=$maxeditdist";

    $ra_vcf_fields->[7] = join ';', @info_entries;
    $ra_vcf_fields->[8] = '.'; # no format/genotype for now
    $ra_vcf_fields->[9] = '.'; # no format/genotype for now

    my $vcf_record = join "\t", @{$ra_vcf_fields};

    print $vcf_fh "$vcf_record\n";
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
