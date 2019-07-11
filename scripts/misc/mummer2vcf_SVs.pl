#!/usr/bin/perl -w
# $Id:$

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);
use GTB::FASTA;
use NISC::Sequencing::Date;
use NHGRI::MUMmer::AlignSet;
our %Opt;

=head1 NAME

mummer2vcf_SVs.pl - Read MUMmer output and write structural variants to VCF.

=head1 SYNOPSIS

  mummer2vcf_SVs.pl --delta <path to delta file of alignments> --diffs <path to output of show-diff> --ref_fasta <path to reference multi-FASTA file> --query_fasta <path to query multi-FASTA file> --outvcf <path to output VCF file>

For complete documentation, run C<mummer2vcf_SVs.pl -man>

=head1 DESCRIPTION

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $query_fasta = $Opt{query_fasta}; # will use values from delta file if not supplied as arguments
my $ref_fasta = $Opt{ref_fasta}; # will use values from delta file if not supplied as arguments

my $delta_file = $Opt{delta};
my $diffs_file = $Opt{diffs};

my $outvcf = $Opt{outvcf};
my $vcf_fh = Open($outvcf, "w");
write_header($vcf_fh) if (!$Opt{noheader});

my $delta_obj = NHGRI::MUMmer::AlignSet->new(-delta_file => $delta_file);

# set up assembly FASTA objects:
$ref_fasta = $delta_obj->{reference_file} if (!$ref_fasta);
$query_fasta = $delta_obj->{query_file} if (!$query_fasta);

my $ref_db = GTB::FASTA->new($ref_fasta);
my $query_db = GTB::FASTA->new($query_fasta);

my $diffs_fh = Open($diffs_file);

write_variants($delta_obj, $diffs_fh, $vcf_fh);

close $vcf_fh;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( );
    GetOptions(\%Opt, qw( delta=s diffs=s ref_fasta=s query_fasta=s outvcf=s maxsize=i includeseqs
                                 refname=s noheader manual help+ verbose version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "mummer2vcf_SVs.pl, ", q$Revision:$, "\n"; }

    if (!($Opt{delta})) {
        print STDERR "Must specify a delta file path with --delta option!\n"; 
        pod2usage(0);
    }

    if (!($Opt{diffs})) {
        print STDERR "Must specify a diffs file path with --diffs option!\n"; 
        pod2usage(0);
    }

    if (!($Opt{outvcf})) {
        print STDERR "Must specify a VCF file path to output confirmed variants!\n"; 
        pod2usage(0);
    }
}

sub write_header {
    my $fh = shift;

    print $fh "##fileformat=VCFv4.2\n";
    my $date_obj = NISC::Sequencing::Date->new(-plain_language => 'today');
    my $year = $date_obj->year();
    my $month = $date_obj->month();
    $month =~ s/^(\d)$/0$1/;
    my $day = $date_obj->day();
    $day =~ s/^(\d)$/0$1/;
    print $fh "##fileDate=$year$month$day\n";
    print $fh "##source=mummer2vcf_SVs.pl\n";
    print $fh "##reference=$Opt{refname}\n" if ($Opt{refname});
    print $fh "##ALT=<ID=DEL,Description=\"Deletion\">\n";
    print $fh "##ALT=<ID=INS,Description=\"Insertion\">\n";
    #print $fh "##ALT=<ID=INV,Description=\"Inversion\">\n";
    print $fh "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Left end coordinate of SV\">\n";
    print $fh "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV:DEL=Deletion, CON=Contraction, INS=Insertion, DUP=Duplication\">\n";
    print $fh "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between ALT and REF alleles (negative for deletions from reference)\">\n";
    #print $fh "##INFO=<ID=PREDREGION,Number=1,Type=String,Description=\"Region of SV Prediction that was refined to create this record\">\n";
    print $fh "##INFO=<ID=HOMAPPLEN,Number=.,Type=Integer,Description=\"Length of alignable homology at event breakpoints as determined by MUMmer\">\n";
    #print $fh "##INFO=<ID=HOMAPPSEQ,Number=.,Type=String,Description=\"Sequence of alignable homology at event breakpoints as determined by MUMmer\">\n";
    #print $fh "##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype\">\n";
    #print $fh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tQUERY";
    print $fh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    print $fh "\n";
}

sub write_variants {
    my $delta_obj = shift;
    my $diff_fh = shift;
    my $outvcf_fh = shift;

    # examine SV calls for reference against contig:

    while (<$diff_fh>) {
        next if (m:^/:);
        next if (m/^NUCMER$/);
        next if (m/^\s*$/);

        if (/^(\S+)\tGAP\t(\d+)\t(\d+)\t(\-{0,1}\d+)\t(\-{0,1}\d+)\t(\-{0,1}\d+)/) {
            my ($ref_entry, $refstart, $refend, $refjump, $queryjump, $svsize) = ($1, $2, $3, $4, $5, $6);

            my $ref1 = $refstart - 1; # last aligned base in ref entry coords
            my $ref2 = $refend + 1; # first aligned base in second alignment in ref entry coords

            # jump values are the number of unaligned bases in ref/query
            if ($refjump != $ref2 - $ref1 - 1) {
                die "Refjump not as expected: $ref_entry, $refstart, $refend, $refjump, $ref1, $ref2\n";
            }

            # find query positions of alignment endpoints in delta file:
            print "Looking for $ref_entry, ref1=$ref1, ref2=$ref2\n" if ($Opt{verbose});
            my ($varcontig, $query1, $query2, $comp) = $delta_obj->find_gap_query_coords($ref_entry, $ref1, $ref2);

            my $calc_queryjump = ($comp) ? $query1 - $query2 - 1 : $query2 - $query1 - 1; # from show-diff code--number of unaligned bases
            if ($queryjump != $calc_queryjump) {
                print "Comp $comp, var contig $varcontig, query1 $query1, query2 $query2, query jump $queryjump, svsize $svsize\n" if ($Opt{verbose});
                print "Calculated queryjump = $calc_queryjump (reported $queryjump)--SKIPPING THIS VARIANT!!!\n";
                next;
            }

            # SVsize is difference betweeen refjump and queryjump, and helps determine SVTYPE
            if ($svsize != $refjump - $queryjump) {
                die "SVsize not as expected: $svsize, $queryjump, $refjump\n";
            }

            my $type = ($refjump <= 0 && $queryjump <= 0) ? (($svsize < 0) ? 'DUP' : 'CON') :
                           (($svsize < 0 ) ? 'INS' : 'DEL'); # we reverse this later on so that deletions have negative SVLEN
            if (($svsize < 0) && ($refjump > 0)) {
                $type = 'SUBSINS';
            }
            elsif (($svsize > 0) && ($queryjump > 0)) {
                $type = 'SUBSDEL';
            }
            elsif ($svsize == 0) { # is this even possible?
                $type = 'SUBS';
            }
            print "Comp $comp, ref jump $refjump, query jump $queryjump, svsize $svsize has type $type\n" if ($Opt{verbose});
            $svsize = abs($svsize);

            my $repeat_bases = ($type eq 'INS' || $type eq 'DUP') ? -1*$refjump : 
                                (($type eq 'DEL' || $type eq 'CON') ? -1*$queryjump : 0);

            if ($type eq 'INS' || $type eq 'DUP') {

                my ($query1p, $query2p);
                foreach my $rh_entry (@{$delta_obj->{entry_pairs}}) {
                    next if (($rh_entry->{ref_entry} ne $ref_entry) || ($rh_entry->{query_entry} ne $varcontig));
                    foreach my $rh_align (@{$rh_entry->{aligns}}) {
                        if (defined($rh_align->{query_matches}->{$ref1})) { # potential query1p
                            my ($rstart, $rend, $qstart, $qend) = ($rh_align->{ref_start},
                                $rh_align->{ref_end}, $rh_align->{query_start}, $rh_align->{query_end});
                            if ($rend == $ref1) { # found $query1
                                if ($query1 != $qend) {
                                    print "Found different query1 value $qend for q1 (prior $query1)---THIS IS BAD\n";
                                }
                            }
                            elsif ($rstart == $ref2) { # found $query1p
                                $query1p = $rh_align->{query_matches}->{$ref1};
                            }
                        }
                        if (defined($rh_align->{query_matches}->{$ref2})) { # potential query2p
                            my ($rstart, $rend, $qstart, $qend) = ($rh_align->{ref_start},
                                $rh_align->{ref_end}, $rh_align->{query_start}, $rh_align->{query_end});
                            if ($rstart == $ref2) { # found $query2
                                if ($query2 != $qstart) {
                                    print "Found different query2 value $qstart for q2 (prior $query2)---THIS IS BAD\n";
                                }
                            }
                            elsif ($rend == $ref1) { # found $query2p
                                $query2p = $rh_align->{query_matches}->{$ref2};
                            }
                        }
                    }
                }

                if (!defined($query1p)) {
                    print "NONREPETITIVE INSERTION--need to check\n";
                    $query1p = $query2 - 1;
                }
                if (!defined($query2p)) {
                    print "NONREPETITIVE INSERTION--need to check\n";
                    $query2p = $query1 - 1;
                }
                print "$ref_entry:$ref1-$ref2:\n" if ($Opt{verbose});
                print "Assigned query1=$query1, query1p=$query1p\n" if ($Opt{verbose});
                print "Assigned query2=$query2, query2p=$query2p\n" if ($Opt{verbose});

                my $rh_var = {'type' => $type,
                              'svsize' => -1.0*$svsize,
                              'repbases' => $repeat_bases,
                              'chrom' => $ref_entry,
                              'pos' => $ref2,
                              'end' => $ref2,
                              'contig' => $varcontig,
                              'altpos' => $query2p,
                              'altend' => $query2,
                              'ref1' => $ref1,
                              'ref2' => $ref2,
                              'refjump' => $refjump,
                              'query1' => $query1,
                              'query2' => $query2,
                              'queryjump' => $queryjump,
                              'comp' => $comp,
                              'query1p' => $query1p,
                              'query2p' => $query2p,
                             };

                write_simple_variant($outvcf_fh, $rh_var);

            }
            elsif ($type eq 'DEL' || $type eq 'CON') {

                print "Type $type comp $comp, varcontig $varcontig query1 $query1, query2 $query2, ref $ref_entry ref1 $ref1, ref2 $ref2\n" if ($Opt{verbose});
                print "Refjump $refjump, query jump $queryjump, svsize $svsize\n" if ($Opt{verbose});
                $delta_obj->find_ref_coords_from_query_coord($query1, $varcontig);
                $delta_obj->find_ref_coords_from_query_coord($query2, $varcontig);

                my ($ref1p, $ref2p);
                foreach my $rh_entry (@{$delta_obj->{entry_pairs}}) {
                    next if (($rh_entry->{ref_entry} ne $ref_entry) || ($rh_entry->{query_entry} ne $varcontig));
                    foreach my $rh_align (@{$rh_entry->{aligns}}) {
                        if (defined($rh_align->{ref_matches}->{$query1})) { # potential ref1p
                            my ($rstart, $rend, $qstart, $qend) = ($rh_align->{ref_start},
                                $rh_align->{ref_end}, $rh_align->{query_start}, $rh_align->{query_end});
                            if ($qend == $query1) { # found $ref1
                                if ($ref1 != $rend) {
                                    print "Found different ref1 value $rend for r1 (prior $ref1)--THIS IS BAD\n";
                                }
                            }
                            elsif ($qstart == $query2) { # found $ref1p
                                $ref1p = $rh_align->{ref_matches}->{$query1};
                            }
                        }
                        if (defined($rh_align->{ref_matches}->{$query2})) { # potential ref2p
                            my ($rstart, $rend, $qstart, $qend) = ($rh_align->{ref_start},
                                $rh_align->{ref_end}, $rh_align->{query_start}, $rh_align->{query_end});
                            if ($qstart == $query2) { # found $ref2
                                if ($ref2 != $rstart) {
                                    print "Found different ref2 value $rstart for r2 (prior $ref2)--THIS IS BAD\n";
                                }
                                else {
                                    print "Confirmed alignment start $ref2\n";
                                }
                            }
                            elsif ($qend == $query1) { # found $ref2p
                                $ref2p = $rh_align->{ref_matches}->{$query2};
                            }
                        }
                    }
                }
                if (!defined($ref1p)) {
                    print "NONREPETITIVE DELETION--need to check\n";
                    $ref1p = $ref2 - 1;
                }
                if (!defined($ref2p)) {
                    print "NONREPETITIVE DELETION--need to check\n";
                    $ref2p = $ref1 - 1;
                }
                print "Assigned ref1=$ref1, ref1p=$ref1p\n" if ($Opt{verbose});
                print "Assigned ref2=$ref2, ref2p=$ref2p\n" if ($Opt{verbose});

                my $rh_var = {'type' => $type,
                              'svsize' => -1.0*$svsize,
                              'repbases' => $repeat_bases,
                              'chrom' => $ref_entry,
                              'pos' => $ref2p,
                              'end' => $ref2,
                              'contig' => $varcontig,
                              'altpos' => $query2,
                              'altend' => $query2,
                              'ref1' => $ref1,
                              'ref2' => $ref2,
                              'refjump' => $refjump,
                              'query1' => $query1,
                              'query2' => $query2,
                              'queryjump' => $queryjump,
                              'comp' => $comp,
                              'ref1p' => $ref1p,
                              'ref2p' => $ref2p,
                             };

                write_simple_variant($outvcf_fh, $rh_var);

            }
            else {

                $svsize = ($type =~ /DEL/) ? -1.0*$svsize : $svsize;
                my $altend = ($comp) ? $query2 + 1 : $query2 - 1;
                my $rh_var = {'type' => $type,
                              'svsize' => $svsize,
                              'repbases' => 0,
                              'chrom' => $ref_entry,
                              'pos' => $ref1,
                              'end' => $ref2-1,
                              'contig' => $varcontig,
                              'altpos' => $query1,
                              'altend' => $altend,
                              'ref1' => $ref1,
                              'ref2' => $ref2,
                              'refjump' => $refjump,
                              'query1' => $query1,
                              'query2' => $query2,
                              'queryjump' => $queryjump,
                              'comp' => $comp,
                             };

                write_simple_variant($outvcf_fh, $rh_var);

            }
        }
    }

    close $diff_fh;
}

sub write_simple_variant {
    my $fh = shift;
    my $rh_var = shift;

    my $vartype = $rh_var->{type};
    my $svsize = $rh_var->{svsize};
    my $repbases = $rh_var->{repbases};
    my $chrom = $rh_var->{chrom};
    my $pos = $rh_var->{pos};
    my $end = $rh_var->{end};
    my $varcontig = $rh_var->{contig};
    my $altpos = $rh_var->{altpos};
    my $altend = $rh_var->{altend};
    my $ref1 = $rh_var->{ref1};
    my $ref2 = $rh_var->{ref2};
    my $refjump = $rh_var->{refjump};
    my $query1 = $rh_var->{query1};
    my $query2 = $rh_var->{query2};
    my $queryjump = $rh_var->{queryjump};
    my $comp = $rh_var->{comp};
    my $query1p = $rh_var->{query1p};
    my $query2p = $rh_var->{query2p};
    my $ref1p = $rh_var->{ref1p};
    my $ref2p = $rh_var->{ref2p};

    if ($vartype eq 'INS' || $vartype eq 'DUP') {
        my ($refseq, $altseq) = ('N', '<INS>');
        $svsize = abs($svsize); # insertions always positive
        if ($Opt{includeseqs}) {
            $refseq = uc($ref_db->seq($chrom, $pos, $pos));
            $altseq = uc($query_db->seq($varcontig, $altpos, $altend)); # GTB::FASTA will reverse complement if necessary unless altpos = altend
            if (($comp) && ($altpos == $altend)) { # need to complement altseq
                $altseq =~ tr/ATGCatgc/TACGtacg/;
            }
        }

        my $compstring = ($comp) ? '_comp' : '';
        my $varstring = "$chrom\t$pos\t.\t$refseq\t$altseq\t.\tPASS\tEND=$end;SVTYPE=$vartype;SVLEN=$svsize;HOMAPPLEN=$repbases;REFWIDENED=$chrom:$ref2-$ref1;CONTIGALTPOS=$varcontig:$altpos-$altend;CONTIGWIDENED=$varcontig:$query2p-$query1p$compstring";
        print $fh "$varstring\n";
    }
    elsif ($vartype eq 'DEL' || $vartype eq 'CON') {
        my ($refseq, $altseq) = ('N', '<DEL>');
        $svsize = -1.0*abs($svsize); # deletions always negative
        if (!$repbases) { # kludgy for now
            $pos++;
            $end--;
            $ref2p++;
            if ($comp) {
                $query2++;
            }
            else {
                $query2--;
            }
        }
        if ($Opt{includeseqs}) {
            $refseq = uc($ref_db->seq($chrom, $pos, $end));
            $altseq = uc($ref_db->seq($chrom, $pos, $pos));
        }
        my $compstring = ($comp) ? '_comp' : '';
        my $varstring = "$chrom\t$pos\t.\t$refseq\t$altseq\t.\tPASS\tEND=$end;SVTYPE=$vartype;SVLEN=$svsize;HOMAPPLEN=$repbases;REFWIDENED=$chrom:$ref2p-$ref1p;CONTIGALTPOS=$varcontig:$query2;CONTIGWIDENED=$varcontig:$query2-$query1$compstring";
        print $fh "$varstring\n";
    }
    else {
        my ($refseq, $altseq) = ('N', 'N');
        $svsize = ($vartype =~ /DEL/) ? -1.0*abs($svsize) : abs($svsize); # insertions always positive
        if ($Opt{includeseqs}) {
            $refseq = uc($ref_db->seq($chrom, $pos, $end));
            $altseq = uc($query_db->seq($varcontig, $altpos, $altend));
        }
        my $compstring = ($comp) ? '_comp' : '';
        my $varstring = "$chrom\t$pos\t.\t$refseq\t$altseq\t.\tPASS\tEND=$end;SVTYPE=$vartype;SVLEN=$svsize;HOMAPPLEN=$repbases;REFWIDENED=$chrom:$ref1-$ref2;CONTIGALTPOS=$varcontig:$altpos;CONTIGWIDENED=$varcontig:$query1-$query2$compstring";
        print $fh "$varstring\n";
    }
}

__END__

=head1 OPTIONS

=over 4

=item B<--delta <path to delta file>>

Specify a delta file produced by MUMmer with the alignments to be used for
retrieving SV sequence information.  Generally, one would use the same
filtered delta file that was used to create the "diff" file (see below).
(Required).

=item B<--diffs <path to diff file>>

Specify a diff file produced by MUMmer's "show-diff" command run with the
-r option on the delta file passed with the --delta option to this program.
(Required).

=item B<--ref_fasta <path to reference multi-fasta file>>

Specify the path to a multi-fasta file containing the sequences used as 
reference in the MUMmer alignment.  If not specified on the command line,
the script uses the reference path obtained by parsing the delta file's
first line.

=item B<--query_fasta <path to query multi-fasta file>>

Specify the path to a multi-fasta file containing the sequences used as 
the query in the MUMmer alignment.  If not specified on the command line,
the script uses the query path obtained by parsing the delta file's
first line.

=item B<--outvcf <path to which to write a new VCF-formatted file>>

Specify the path to which to write a new VCF file containing the structural
variants discovered in this comparison.  BEWARE: if this file already 
exists, it will be overwritten!

=item B<--refname <string to include as the reference name in the VCF header>>

Specify a string to be written as the reference name in the ##reference line 
of the VCF header.

=item B<--noheader>

Flag option to suppress printout of the VCF header.

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=back

=head1 AUTHOR

 Nancy F. Hansen - nhansen@mail.nih.gov

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
