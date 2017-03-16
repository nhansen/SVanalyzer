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

SVrefine.pl - Read regions from a BED file and use MUMmer alignments of an assembly to the reference to refine structural variants in those regions and print them out in VCF format.

=head1 SYNOPSIS

  SVrefine.pl --delta <path to delta file of alignments> --regions <path to BED-formatted file of regions> --ref_fasta <path to reference multi-FASTA file> --query_fasta <path to query multi-FASTA file> --outvcf <path to output VCF file>

For complete documentation, run C<SVrefine.pl -man>

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $query_fasta = $Opt{query_fasta}; # will use values from delta file if not supplied as arguments
my $ref_fasta = $Opt{ref_fasta}; # will use values from delta file if not supplied as arguments

my $delta_file = $Opt{delta};
my $regions_file = $Opt{regions};

my $samplename = $Opt{samplename};

my $outvcf = $Opt{outvcf};
my $vcf_fh = Open($outvcf, "w");

my $refbedfile = $outvcf; # file to write no coverage and coverage that agrees with reference
$refbedfile =~ s/\.vcf(.*)$//;
$refbedfile .= ".bed";
my $refregions_fh = Open($refbedfile, "w"); # will write regions with support for reference sequence across inquiry regions to bed formatted file 

write_header($vcf_fh) if (!$Opt{noheader});

my $delta_obj = NHGRI::MUMmer::AlignSet->new(
                  -delta_file => $delta_file,
                  -storerefentrypairs => 1, # object stores hash of each ref/query entry's align pairs (or "edges")
                  -storequeryentrypairs => 1);

# set up assembly FASTA objects:
$ref_fasta = $delta_obj->{reference_file} if (!$ref_fasta);
$query_fasta = $delta_obj->{query_file} if (!$query_fasta);

my $ref_db = GTB::FASTA->new($ref_fasta);
my $query_db = GTB::FASTA->new($query_fasta);

my $regions_fh = Open($regions_file); # regions to refine

process_regions($delta_obj, $regions_fh, $vcf_fh, $refregions_fh, $ref_db, $query_db);

close $vcf_fh;
close $regions_fh;
close $refregions_fh;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( samplename => 'SAMPLE', buffer => 1000, bufseg => 50 );
    GetOptions(\%Opt, qw( delta=s regions=s ref_fasta=s query_fasta=s outvcf=s buffer=i 
                           bufseg=i verbose includeseqs samplename=s noheader manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "SVrefine.pl, ", q$Revision: 7771 $, "\n"; }

    if (!($Opt{delta})) {
        print STDERR "Must specify a delta file path with --delta option!\n"; 
        pod2usage(0);
    }

    if (!($Opt{regions})) {
        print STDERR "Must specify a regions BED file path with --regions option!\n"; 
        pod2usage(0);
    }

    if (!($Opt{outvcf})) {
        print STDERR "Must specify a VCF file path to output confirmed variants!\n"; 
        pod2usage(0);
    }
}

sub write_header {
    my $fh = shift;

    print $fh "##fileformat=VCFv4.3\n";
    my $date_obj = NISC::Sequencing::Date->new(-plain_language => 'today');
    my $year = $date_obj->year();
    my $month = $date_obj->month();
    $month =~ s/^(\d)$/0$1/;
    my $day = $date_obj->day();
    $day =~ s/^(\d)$/0$1/;
    print $fh "##fileDate=$year$month$day\n";
    print $fh "##source=SVrefine.pl\n";
    print $fh "##reference=$Opt{refname}\n" if ($Opt{refname});
    print $fh "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of SV\">\n";
    print $fh "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV:DEL=Deletion, INS=Insertion, INV=Inversion\">\n";
    print $fh "##INFO=<ID=REPTYPE,Number=1,Type=String,Description=\"Type of SV, with designation of uniqueness of new or deleted sequence:SIMPLEDEL=Deletion of at least some unique sequence, SIMPLEINS=Insertion of at least some unique sequence, CONTRAC=Contraction, or deletion of sequence entirely similar to remaining sequence, DUP=Duplication, or insertion of sequence entirely similar to pre-existing sequence, INV=Inversion, SUBSINS=Insertion of new sequence with alteration of some pre-existing sequence, SUBSDEL=Deletion of sequence with alteration of some remaining sequence\">\n";
    print $fh "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between ALT and REF alleles (negative for deletions from reference, positive for insertions to reference)\">\n";
    print $fh "##INFO=<ID=BREAKSIMLENGTH,Number=1,Type=Integer,Description=\"Length of alignable similarity at event breakpoints as determined by the aligner\">\n";
    print $fh "##INFO=<ID=REFWIDENED,Number=1,Type=String,Description=\"Widened boundaries of the event in the reference allele\">\n";
    print $fh "##INFO=<ID=ALTWIDENED,Number=1,Type=String,Description=\"Widened boundaries of the event in the alternate allele\">\n";
    print $fh "##INFO=<ID=ALTPOS,Number=1,Type=String,Description=\"Position (CHROM:POS-END) of the event in the sample assembly\">\n";
    print $fh "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    print $fh "##FORMAT=<ID=CONTIG,Number=1,Type=String,Description=\"Supporting contigs, in same order as alleles reported in genotype\">\n";
    print $fh "##FORMAT=<ID=GTMATCH,Number=1,Type=Character,Description=\"Genotype match type: E=Exact match, H=Position boundary match, L=Inexact match\">\n";
    print $fh "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
}

sub process_regions {
    my $delta_obj = shift;
    my $regions_fh = shift;
    my $outvcf_fh = shift;
    my $outref_fh = shift;
    my $ref_db = shift;
    my $query_db = shift;

    # examine regions for broken alignments:

    my $rh_ref_entrypairs = $delta_obj->{refentrypairs}; # store entry pairs by ref entry for quick retrieval
    my $rh_query_entrypairs = $delta_obj->{queryentrypairs}; # store entry pairs by query entry for quick retrieval

    while (<$regions_fh>) {
        chomp;
        my ($chr, $start, $end, $remainder) = split /\s/, $_;
        $start++; # 0-based to 1-based
        if ($end < $start) {
            die "End less than start illegal bed format: $_";
        }

        if (!check_region($ref_db, $chr, $start, $end)) {
            print "SKIPPING REGION $chr\t$start\t$end because it is too close to chromosome end!\n" if ($Opt{verbose});
            next;
        }

        my $ra_rentry_pairs = $rh_ref_entrypairs->{$chr}; # everything aligning to this chromosome

        # check small regions to the left and right of region of interest for aligning contigs:
        my $refleftstart = $start - $Opt{buffer};
        my $refleftend = $start - $Opt{buffer} + $Opt{bufseg} - 1;
        my $refrightstart = $end + $Opt{buffer} - $Opt{bufseg} + 1;
        my $refrightend = $end + $Opt{buffer};

        # find entries with at least one align spanning left and right contig align, with same comp value (can be a single align, or more than one)

        my @valid_entry_pairs = find_valid_entry_pairs($ra_rentry_pairs, $refleftstart, $refleftend, $refrightstart, $refrightend);

        if (!@valid_entry_pairs) { # no consistent alignments across region
            print $outref_fh "NOCOV\t$chr\t$start\t$end\n";
            next;
        }
        
        my @aligned_contigs = map { $_->{query_entry} } @valid_entry_pairs;
        my $contig_string = join ':', @aligned_contigs;
        my $no_entry_pairs = @valid_entry_pairs;
        print "VALID PAIRS $chr\t$start\t$end $no_entry_pairs aligning contigs: $contig_string\n" if (($no_entry_pairs > 1) && ($Opt{verbose}));

        my @refaligns = ();
        foreach my $rh_rentry_pair (@valid_entry_pairs) { # for each contig aligned across region
            my $contig = $rh_rentry_pair->{query_entry};
            my $ra_qentry_pairs = $rh_query_entrypairs->{$contig};
            my @contig_aligns = ();
            my $comp = $rh_rentry_pair->{comp}; # determines whether to sort contig matches in decreasing contig coordinate order
            foreach my $rh_entrypair (@{$ra_qentry_pairs}) {
                push @contig_aligns, @{$rh_entrypair->{aligns}};
            }
            my @sorted_contigaligns = ($comp) ? 
                      sort {$b->{query_start} <=> $a->{query_start}} @contig_aligns : 
                      sort {$a->{query_start} <=> $b->{query_start}} @contig_aligns;
            my @region_aligns = (); # will contain only contig alignments in region aligned to reference region
            my $in_region = 0;
            foreach my $rh_contigalign (@sorted_contigaligns) {
                if (!($in_region) && ($rh_contigalign->{ref_entry} eq $chr && $rh_contigalign->{ref_end} > $refleftend)) {
                    $in_region = 1;
                }
                elsif (($in_region) && ($rh_contigalign->{ref_entry} eq $chr && $rh_contigalign->{ref_start} > $refrightend)) {
                    last;
                }

                if ($in_region) {
                    push @region_aligns, $rh_contigalign;
                }
            }
            my @contigrefentries = map { $_->{ref_entry} } @region_aligns;
            my $refstring = join ':', @contigrefentries;
            print "CONTIG $contig matches to ref entries $refstring\n" if ($Opt{verbose});

            # is there only one alignment to $chr and does it span the region of interest?  If so, output reference coverage:

            my @ref_aligns = grep {$_->{ref_entry} eq $chr} @region_aligns;
            if (@ref_aligns == 1 && $ref_aligns[0]->{ref_start} < $refleftstart && 
                      $ref_aligns[0]->{ref_end} > $refrightend) {
                my $ref_align = $ref_aligns[0];
                print $outref_fh "HOMREF\t$chr\t$start\t$end\t$chr\t$ref_align->{ref_start}\t$ref_align->{ref_end}\t$ref_align->{query_entry}\t$ref_align->{query_start}\t$ref_align->{query_end}\n";
                next;
            }

            if (@region_aligns > 1) { # potential SV--print out alignments if verbose option
                foreach my $rh_align (@region_aligns) {
                    my $ref_entry = $rh_align->{ref_entry};
                    my $query_entry = $rh_align->{query_entry};
                    my $ref_start = $rh_align->{ref_start};
                    my $ref_end = $rh_align->{ref_end};
                    my $query_start = $rh_align->{query_start};
                    my $query_end = $rh_align->{query_end};
    
                    print "ALIGN\t$chr\t$start\t$end\t$ref_entry\t$ref_start\t$ref_end\t$query_entry\t$query_start\t$query_end\n" if ($Opt{verbose});
                }
            }

            # find left and right matches among aligns:
           
            if ((@region_aligns <= 3) && !(grep {$_->{ref_entry} ne $chr} @region_aligns)) { # simple insertion, deletion, or inversion
                my @simple_breaks = ();
                my @inversion_aligns = ();
                if (@region_aligns == 2 && $region_aligns[0]->{comp} == $region_aligns[1]->{comp}) {
                    # order them, if possible:
                    @simple_breaks = find_breaks($region_aligns[0], $region_aligns[1]);
                    print "ONE SIMPLE BREAK\n" if (@simple_breaks==1);
                }
                elsif (($region_aligns[0]->{comp} == $region_aligns[1]->{comp}) && 
                    ($region_aligns[1]->{comp} == $region_aligns[2]->{comp})) {
                    @simple_breaks = find_breaks($region_aligns[0], $region_aligns[1], $region_aligns[2]);
                    print "TWO SIMPLE BREAKS\n" if (@simple_breaks==2);
                }
                elsif (($region_aligns[0]->{comp} != $region_aligns[1]->{comp}) &&
                       ($region_aligns[1]->{comp} != $region_aligns[2]->{comp})) {
                    @inversion_aligns = ($region_aligns[0], $region_aligns[1], $region_aligns[2]);
                    print "INVERSION!\n";
                }

                foreach my $ra_simple_break (@simple_breaks) {
                    my $left_align = $ra_simple_break->[0];
                    my $right_align = $ra_simple_break->[1];
    
                    my $ref1 = $left_align->{ref_end};
                    my $ref2 = $right_align->{ref_start};
                    my $refjump = $ref2 - $ref1 - 1;
               
                    my $query1 = $left_align->{query_end}; 
                    my $query2 = $right_align->{query_start}; 

                    my $queryjump = ($comp) ? $query1 - $query2 - 1 : $query2 - $query1 - 1; # from show-diff code--number of unaligned bases
                    my $svsize = $refjump - $queryjump;
                 
                    my $type = ($refjump <= 0 && $queryjump <= 0) ? (($svsize < 0) ? 'DUP' : 'CONTRAC') :
                               (($svsize < 0 ) ? 'SIMPLEINS' : 'SIMPLEDEL'); # we reverse this later on so that deletions have negative SVLEN

                    if (($svsize < 0) && ($refjump > 0)) {
                        $type = 'SUBSINS';
                    }
                    elsif (($svsize > 0) && ($queryjump > 0)) {
                        $type = 'SUBSDEL';
                    }
                    elsif ($svsize == 0) { # is this even possible?
                        $type = 'SUBS';
                    }
                    $svsize = abs($svsize);
                    my $repeat_bases = ($type eq 'SIMPLEINS' || $type eq 'DUP') ? -1*$refjump : 
                                        (($type eq 'SIMPLEDEL' || $type eq 'CONTRAC') ? -1*$queryjump : 0);
  
                    if ($repeat_bases > 100*$svsize) { # likely alignment artifact?
                        next;
                    }
                    print "TYPE $type size $svsize (REFJUMP $refjump QUERYJUMP $queryjump)\n"; 
    
                    if ($type eq 'SIMPLEINS' || $type eq 'DUP') {
                        process_insertion($delta_obj, $left_align, $right_align, $chr, $contig,
                              $ref1, $query1, $ref2, $query2, $comp, $type, $svsize, 
                              $refjump, $queryjump, $repeat_bases, $outvcf_fh);
                    }
                    elsif ($type eq 'SIMPLEDEL' || $type eq 'CONTRAC') {
                        process_deletion($delta_obj, $left_align, $right_align, $chr, $contig,
                              $ref1, $query1, $ref2, $query2, $comp, $type, $svsize, 
                              $refjump, $queryjump, $repeat_bases, $outvcf_fh);
                    }
                    else {
                        my $altpos = ($comp) ? $query1 - 1 : $query1 + 1;
                        my $altend = ($comp) ? $query2 + 1 : $query2 - 1;
                        my $rh_var = {'type' => $type,
                                      'svsize' => $svsize,
                                      'repbases' => 0,
                                      'chrom' => $chr,
                                      'pos' => $ref1 + 1,
                                      'end' => $ref2 - 1,
                                      'contig' => $contig,
                                      'altpos' => $altpos,
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

                # is it an inversion?
                if (@inversion_aligns) {

                }
            }
        }
    }

    close $regions_fh;
}

sub check_region {
    my $ref_db = shift;
    my $chr = shift;
    my $start = shift;
    my $end = shift;

    my $chrlength = $ref_db->len($chr);
    if ((0 >= $start - $Opt{buffer}) || ($chrlength < $end + $Opt{buffer})) {
        print "REGION $chr:$start-$end is too close to an end of the chromosome $chr!\n";
        return 0;
    }
    else {
        return 1;
    }
}

sub find_valid_entry_pairs {
    my $ra_rentry_pairs = shift;
    my $ref_left_start = shift;
    my $ref_left_end = shift;
    my $ref_right_start = shift;
    my $ref_right_end = shift;

    my @valid_entry_pairs = ();
    foreach my $rh_entrypair (@{$ra_rentry_pairs}) {
        my @left_span_aligns = grep {$_->{ref_start} < $ref_left_start && $_->{ref_end} > $ref_left_end} @{$rh_entrypair->{aligns}};
        my @right_span_aligns = grep {$_->{ref_start} < $ref_right_start && $_->{ref_end} > $ref_right_end} @{$rh_entrypair->{aligns}};
        if (!(@left_span_aligns) || !(@right_span_aligns)) { # need to span both sides
            next;
        }
        else {
            my $no_left_aligns = @left_span_aligns;
            my $no_right_aligns = @right_span_aligns;
            if ($no_left_aligns != 1 || $no_right_aligns != 1) {
                print "Found $no_left_aligns left aligns and $no_right_aligns right aligns!\n";
            }
        }

        # want to find the same comp value on left and right side:
        if (grep {$_->{comp} == 0} @left_span_aligns && grep {$_->{comp} == 0} @right_span_aligns) {
            $rh_entrypair->{comp} = 0;
            push @valid_entry_pairs, $rh_entrypair;
        }
        elsif (grep {$_->{comp} == 1} @left_span_aligns && grep {$_->{comp} == 1} @right_span_aligns) {
            $rh_entrypair->{comp} = 1;
            push @valid_entry_pairs, $rh_entrypair;
        }
        else {
            print "Odd line up of alignments for contig $rh_entrypair->{query_entry}\n";
            next;
        }
    }

    return @valid_entry_pairs;
}

sub find_breaks {
    my @aligns = @_;

    # right now this routine accepts 2 or 3 alignments, all of which have the same value of "comp".
    # it attempts to order them left to right wrt the reference, and create pairs of alignments
    # that represent breaks corresponding to structural variants

    my @valid_breaks = ();
    my @sorted_aligns = sort {$a->{ref_start} <=> $b->{ref_start}} @aligns;
    for (my $left_index = 0; $left_index <= $#sorted_aligns - 1; $left_index++) {
        print "LEFTINDEX=$left_index!\n";
        my $left_align = $sorted_aligns[$left_index];
        my $right_align = $sorted_aligns[$left_index + 1];

        # is it a valid pair of aligns?

        if ($left_align->{ref_end} < $right_align->{ref_end}) {
            if ((!$left_align->{comp} && # should line up left to right
                   $left_align->{query_start} < $right_align->{query_start} && 
                   $left_align->{query_end} < $right_align->{query_end}) ||
                ($left_align->{comp} && # query decreasing
                   $left_align->{query_start} > $right_align->{query_start} && 
                   $left_align->{query_end} > $right_align->{query_end})) {

                push @valid_breaks, [$left_align, $right_align];
            }
            else {
                print "INVALID ALIGNS $left_align->{query_start}-$left_align->{query_end}, then $right_align->{query_start}-$right_align->{query_end}\n";
            }
        }
        else {
            print "REF ENDS OUT OF ORDER $left_align->{ref_start}-$left_align->{ref_end}, then $right_align->{ref_start}-$right_align->{ref_end}\n";
        }
    }

    return @valid_breaks;
}

sub process_insertion {
    my $delta_obj = shift;
    my $left_align = shift;
    my $right_align = shift;
    my $chr = shift;
    my $contig = shift;
    my $ref1 = shift;
    my $query1 = shift;
    my $ref2 = shift;
    my $query2 = shift;
    my $comp = shift;
    my $type = shift;
    my $svsize = shift;
    my $refjump = shift;
    my $queryjump = shift;
    my $repeat_bases = shift;
    my $outvcf_fh = shift;

    print "Type $type comp $comp, varcontig $contig query1 $query1, query2 $query2, ref $chr ref1 $ref1, ref2 $ref2\n" if ($Opt{verbose});
    print "Refjump $refjump, query jump $queryjump, svsize $svsize\n" if ($Opt{verbose});

    $delta_obj->find_query_coords_from_ref_coord($ref1, $chr);
    $delta_obj->find_query_coords_from_ref_coord($ref2, $chr);

    if (!$left_align->{query_matches}->{$ref1} || $query1 != $left_align->{query_matches}->{$ref1}) {
        die "QUERY1 $query1 doesn\'t match $left_align->{query_matches}->{$ref1}\n"; 
    }
    
    if (!$right_align->{query_matches}->{$ref2} || $query2 != $right_align->{query_matches}->{$ref2}) {
        die "QUERY2 $query2 doesn\'t match $right_align->{query_matches}->{$ref2}\n"; 
    }
    
    my $query1p = $right_align->{query_matches}->{$ref1};
    my $query2p = $left_align->{query_matches}->{$ref2};

    print "query1p $query1p corresponds to ref $ref1 in second align, query2p $query2p corresponds to ref2 $ref2 in first align\n" if ($Opt{verbose});

    if (!defined($query1p)) {
        print "NONREPETITIVE INSERTION--need to check\n";
        $query1p = $query2 - 1;
    }
    if (!defined($query2p)) {
        print "NONREPETITIVE INSERTION--need to check\n";
        $query2p = $query1 - 1;
    }

    my $altpos = ($comp) ? $query2p + 1 : $query2p - 1;
    my $altend = ($comp) ? $query2 + 1 : $query2 - 1;
    my $rh_var = {'type' => $type,
                  'svsize' => -1.0*$svsize,
                  'repbases' => $repeat_bases,
                  'chrom' => $chr,
                  'pos' => $ref2 - 1,
                  'end' => $ref2 - 1,
                  'contig' => $contig,
                  'altpos' => $altpos,
                  'altend' => $altend,
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

sub process_deletion {
    my $delta_obj = shift;
    my $left_align = shift;
    my $right_align = shift;
    my $chr = shift;
    my $contig = shift;
    my $ref1 = shift;
    my $query1 = shift;
    my $ref2 = shift;
    my $query2 = shift;
    my $comp = shift;
    my $type = shift;
    my $svsize = shift;
    my $refjump = shift;
    my $queryjump = shift;
    my $repeat_bases = shift;
    my $outvcf_fh = shift;

    print "Type $type comp $comp, varcontig $contig query1 $query1, query2 $query2, ref $chr ref1 $ref1, ref2 $ref2\n" if ($Opt{verbose});
    print "Refjump $refjump, query jump $queryjump, svsize $svsize\n" if ($Opt{verbose});

    $delta_obj->find_ref_coords_from_query_coord($query1, $contig);
    $delta_obj->find_ref_coords_from_query_coord($query2, $contig);

    if (!$left_align->{ref_matches}->{$query1} || $ref1 != $left_align->{ref_matches}->{$query1}) {
        die "REF1 $ref1 doesn\'t match $left_align->{ref_matches}->{$query1}\n"; 
    }
    
    if (!$right_align->{ref_matches}->{$query2} || $ref2 != $right_align->{ref_matches}->{$query2}) {
        die "REF2 $ref2 doesn\'t match $right_align->{ref_matches}->{$query2}\n"; 
    }
    
    my $ref1p = $right_align->{ref_matches}->{$query1};
    my $ref2p = $left_align->{ref_matches}->{$query2};

    print "ref1p $ref1p corresponds to query $query1 in second align, ref2p $ref2p corresponds to query2 $query2 in first align\n" if ($Opt{verbose});

    if (!defined($ref1p)) {
        print "NONREPETITIVE DELETION--need to check\n";
        $ref1p = $ref2 - 1;
    }
    if (!defined($ref2p)) {
        print "NONREPETITIVE DELETION--need to check\n";
        $ref2p = $ref1 - 1;
    }

    my $altpos = ($comp) ? $query2 + 1 : $query2 - 1;
    my $altend = ($comp) ? $query2 + 1 : $query2 - 1;
    my $rh_var = {'type' => $type,
                  'svsize' => -1.0*$svsize,
                  'repbases' => $repeat_bases,
                  'chrom' => $chr,
                  'pos' => $ref2p - 1,
                  'end' => $ref2 - 1,
                  'contig' => $contig,
                  'altpos' => $altpos,
                  'altend' => $altend,
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

    if ($vartype eq 'SIMPLEINS' || $vartype eq 'DUP') {
        my ($refseq, $altseq) = ('N', '<INS>');
        $svsize = abs($svsize); # insertions always positive
        if ($Opt{includeseqs}) {
            $refseq = uc($ref_db->seq($chrom, $pos, $pos));
            if ($altpos > $altend) { # extract and check:
                $altseq = uc($query_db->seq($varcontig, $altend, $altpos)); # GTB::FASTA will reverse complement if necessary unless altpos = altend
                if ($altseq =~ /[^ATGCatgcNn]/) { # foundit!
                    die "Seq $varcontig:$altend-$altpos has non ATGC char!\n";
                }
            }
            print "Retrieving $varcontig:$altpos-$altend\n";
            $altseq = uc($query_db->seq($varcontig, $altpos, $altend)); # GTB::FASTA will reverse complement if necessary unless altpos = altend
            if (($comp) && ($altpos == $altend)) { # need to complement altseq
                $altseq =~ tr/ATGCatgc/TACGtacg/;
            }
        }

        my $compstring = ($comp) ? '_comp' : '';
        my $svtype = 'INS';

        my $varstring = "$chrom\t$pos\t.\t$refseq\t$altseq\t.\tPASS\tEND=$end;SVTYPE=$svtype;REPTYPE=$vartype;SVLEN=$svsize;BREAKSIMLENGTH=$repbases;REFWIDENED=$chrom:$ref2-$ref1;ALTPOS=$varcontig:$altpos-$altend$compstring;ALTWIDENED=$varcontig:$query2p-$query1p$compstring";
        print $fh "$varstring\n";
    }
    elsif ($vartype eq 'SIMPLEDEL' || $vartype eq 'CONTRAC') {
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
            $altseq = uc($query_db->seq($varcontig, $altpos, $altpos));
            if ($comp) {
                $altseq =~ tr/ATGC/TACG/;
            }
        }
        my $svtype = 'DEL';

        my $compstring = ($comp) ? '_comp' : '';
        my $varstring = "$chrom\t$pos\t.\t$refseq\t$altseq\t.\tPASS\tEND=$end;SVTYPE=$svtype;REPTYPE=$vartype;SVLEN=$svsize;BREAKSIMLENGTH=$repbases;REFWIDENED=$chrom:$ref2p-$ref1p;ALTPOS=$varcontig:$altpos-$altpos$compstring;ALTWIDENED=$varcontig:$query2-$query1$compstring";
        print $fh "$varstring\n";
    }
    else { # complex "SUBS" types--note, this code is NOT appropriate for inversions!
        my ($refseq, $altseq) = ('N', 'N');
        $svsize = ($vartype =~ /DEL/) ? -1.0*abs($svsize) : abs($svsize); # insertions always positive
        my $svtype = ($vartype =~ /DEL/) ? 'DEL' : 'INS';

        if ($Opt{includeseqs}) {
            $refseq = uc($ref_db->seq($chrom, $pos, $end));
            $altseq = uc($query_db->seq($varcontig, $altpos, $altend));
            if (($altpos == $altend) && ($comp)) {
                $altseq =~ tr/ATGC/TACG/;
            }
        }
        my $compstring = ($comp) ? '_comp' : '';
        my $varstring = "$chrom\t$pos\t.\t$refseq\t$altseq\t.\tPASS\tEND=$end;SVTYPE=$svtype;REPTYPE=$vartype;SVLEN=$svsize;BREAKSIMLENGTH=$repbases;REFWIDENED=$chrom:$ref1-$ref2;ALTPOS=$varcontig:$altpos-$altend$compstring;ALTWIDENED=$varcontig:$query1-$query2$compstring";
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

=item B<--regions <path to a BED file of regions>>

Specify a BED file of regions to be investigated for structural variants
in the assembly (i.e., the query fasta file).
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

=item B<--samplename <string to include as the sample name in the "CHROM" line>>

Specify a string to be written as the sample name in the header specifying a 
genotype column in the VCF line beginning with "CHROM".

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
