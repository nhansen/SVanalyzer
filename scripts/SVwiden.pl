#!/usr/bin/perl -w
# $Id:$

use strict;

use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);
use GTB::FASTA;
use NHGRI::MUMmer::AlignSet;

our %Opt;

# Custom tags added by this script to describe widened variants:
our %TAGSTOADD = ('INFO' => [ 
     '<ID=BREAKSIMLENGTH,Number=1,Type=Integer,Description="Length of alignable similarity at event breakpoints as determined by the aligner">',
     '<ID=REFWIDENED,Number=1,Type=String,Description="Widened boundaries of the event in the reference allele">',
     '<ID=ALTWIDENED,Number=1,Type=String,Description="Widened boundaries of the event in the alternate allele">' 
                   ] );

=head1 NAME

SVwiden.pl - Read a VCF file and use MUMmer to determine widened coordinates for 
structural variants, adding custom tags to the VCF record.


=head1 SYNOPSIS

  SVwiden.pl --invcf <path to input VCF file> --ref <path to reference multi-FASTA file> --outvcf <path to output VCF file>

For complete documentation, run C<SVwiden.pl -man>

=cut

#------------
# Begin MAIN 
#------------

process_commandline();

my $ref_fasta = $Opt{ref};
my $ref_db = GTB::FASTA->new($ref_fasta);

my $invcf = $Opt{invcf};
my $invcf_fh = Open($invcf, "r");

my $outvcf = $Opt{outvcf};
my $outvcf_fh = Open($outvcf, "w");

write_header($invcf_fh, $outvcf_fh) if (!$Opt{noheader});

while (<$invcf_fh>) {
    next if (/^#/);

    my $vcf_line = $_;
    my ($chr, $pos, $end, $id, $ref, $alt) = parse_vcf_line($vcf_line);

    my $rh_tags = process_variant($chrom, $pos, $end, $ref, $alt, $ref_db);
    write_vcf($outvcf_fh, $vcf_line, $rh_tags);
}

close $outvcf_fh;

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = ( );
    GetOptions(\%Opt, qw( invcf=s ref=s outvcf=s verbose noheader manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "SVwiden.pl, ", q$Revision:$, "\n"; }

    if (!($Opt{invcf})) {
        print STDERR "Must specify a VCF file path of variants to widen with --invcf option!\n"; 
        pod2usage(0);
    }

    if (!($Opt{outvcf})) {
        print STDERR "Must specify a VCF file path to output variants with custom widened tags with --outvcf option!\n"; 
        pod2usage(0);
    }

    if (!($Opt{ref})) {
        print STDERR "Must specify a reference FASTA file path with --ref option!\n"; 
        pod2usage(0);
    }

}

sub parse_vcf_line {
    my $vcf_line = shift;

    if ($vcf_line =~ /^(\S+)\t(\d+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/) {
        my ($chr, $start, $id, $ref, $alt, $info) = ($1, $2, $3, $4, $5, $8);
        $ref = uc($ref);
        $alt = uc($alt);

        my $end; # if we have sequence, end will be determined as last base of REF, otherwise, from END=
        if (($ref =~ /^([ATGC]+)$/) && ($alt =~ /^([ATGC]+)$/)) {
            $end = $start + length($ref) - 1;
        }
        else { # check for END= INFO tag
            if ($info =~ /END=\s*(\d+)/) {
                $end = $1;
            }
            else {
                die "Variants without END in INFO field must have sequence alleles in REF and ALT fields!\n";
            }
        }

        if (!$Opt{'ignore_length'} && ($ref =~ /^([ATGC]+)$/) && ($end - $start + 1 != length($ref))) {
            die "Length of reference allele does not match provided POS, END!  Use --ignore_length option to ignore this discrepancy.\n";
        }

        return ($chr, $start, $end, $id, $ref, $alt);
    }
    else {
        die "Unexpected VCF line:\n$vcf_line";
    }
}

sub write_header {
    my $infh = shift;
    my $outfh = shift;

    my $lasttag;

    while (<$infh>) {
        if (/^##([^=]+)=/) { # transfer all header lines from the input VCF to the output VCF
            my $thistag = $1; 
            my $thisline = $_; 
            if (($lasttag) && ($lasttag != $thistag)) { # add our tags immediately after similar tags, if they exist
               if ($TAGSTOADD{$thistag}) {
                   foreach my $tag (@{$TAGSTOADD{$thistag}}) {
                       print $outfh "##$thistag=$tag\n";
                   }
                   delete $TAGSTOADD{$thistag};
               }
            }
            print $outfh $thisline;
            $lasttag = $thistag;
        }
        elsif (/^#CHROM/) {
            my $thisline = $_; 
            # any of our tags left to print?
            foreach my $thistag (keys %TAGSTOADD) {
                foreach my $tag (@{$TAGSTOADD{$thistag}}) {
                    print $outfh "##$thistag=$tag\n";
                }
                delete $TAGSTOADD{$thistag};
            }
            print $outfh $thisline;
        }
    }
    return 1;
}

sub process_variant {
    my $chrom = shift;
    my $pos = shift;
    my $end = shift;
    my $ref = shift;
    my $alt = shift;
    my $ref_db = shift;

    # examine regions for broken alignments:

}

sub construct_alt_hap {
    my $fai_obj = shift;
    my $chr = shift;
    my $left = shift;
    my $right = shift;
    my $start = shift;
    my $end = shift;
    my $r_ref = shift;
    my $r_alt = shift;

    #my $alt_allele = uc($gtb_ref->seq("$chr:$left-$right"));
    #print "Extracting sequence $chr:$left-$right\n";
    my $alt_allele = uc($fai_obj->fetch("$chr:$left-$right"));
    my $offset = $start - $left;
    my $length = $end - $start + 1;
    my $extracted_ref = substr($alt_allele, $offset, $length); 
    if (($Opt{check_ref}) && ($extracted_ref ne ${$r_ref})) {
        die "Extracted reference for $chr:$start-$end does not match provided REF!\n";
    }

    substr($alt_allele, $offset, $length) = ${$r_alt}; 

    return $alt_allele;
}

sub run_mummer {
    my $fasta1 = shift;
    my $fasta2 = shift;
    my $nucmer_name = shift;

    my $mummer_cmd_file = "$workingdir/run_mummer.$nucmer_name.sh";
    my $mummer_fh = Open("$mummer_cmd_file", "w");
    print $mummer_fh "#!/bin/bash\n/home/nhansen/projects/MUMmer/MUMmer3.23/nucmer -o -p $workingdir/nucmer.$nucmer_name -maxmatch $fasta1 $fasta2\n";
    close $mummer_fh;
    chmod 0755, $mummer_cmd_file;
    
    my $cmd = "$mummer_cmd_file > $workingdir/nucmer.$nucmer_name.out 2> $workingdir/nucmer.$nucmer_name.err";
    system($cmd) == 0
        or print STDERR "Something went wrong running $mummer_cmd_file!\n";

    #if ($Opt{'cleanup'}) {
        unlink $fasta1;
        unlink $fasta2;
        unlink "$workingdir/nucmer.$nucmer_name.err";
        unlink "$workingdir/nucmer.$nucmer_name.coords";
        unlink "$workingdir/nucmer.$nucmer_name.out";
        unlink "$workingdir/run_mummer.$nucmer_name.sh";
    #}

    return "$workingdir/nucmer.$nucmer_name.delta";
}

sub write_vcf {
    my $vcf_fh = shift;
    my $ref_fh = shift;
    my $ra_variant_lines = shift;
    my $ra_ref_cov = shift;

    # write VCF lines in order, adding genotypes and avoiding redundancy:
    my @sorted_vcf_lines = sort byposthenend @{$ra_variant_lines};

    my %written = ();
    foreach my $vcf_line (@sorted_vcf_lines) {
        my ($pos, $end) = ($vcf_line =~ /^(\S+)\s(\d+).*END=(\d+)/) ? ($2, $3) : (0, 0);
        my $gt = covered($ra_ref_cov, $pos, $end) ? '0/1' : '1';
        if (!$written{"$pos:$end"}) {
            print $vcf_fh "$vcf_line\tGT\t$gt\n";
            $written{"$pos:$end"} = 1;
        }
    }

    # write reference regions to BED formatted "ref coverage" file
    foreach my $ra_refregion (sort {$a->[1] <=> $b->[1]} @{$ra_ref_cov}) {
        my $region_line = join "\t", @{$ra_refregion};
        print $ref_fh "$region_line\n";
    }
}

sub covered {
    my $ra_ref_cov = shift;
    my $pos = shift;
    my $end = shift;

    foreach my $ra_hr (@{$ra_ref_cov}) {
        if ($ra_hr->[1] <= $pos && $ra_hr->[2] >= $end) {
            return 1;
        }
    }
    return 0;
}

sub byposthenend {
    my ($pos_a, $end_a) = ($a =~ /^(\S+)\s(\d+).*END=(\d+)/) ? ($2, $3) : (0, 0);
    my ($pos_b, $end_b) = ($b =~ /^(\S+)\s(\d+).*END=(\d+)/) ? ($2, $3) : (0, 0);

    if ($pos_a != $pos_b) {
        return $pos_a <=> $pos_b;
    }
    else { # equal starts, sort ends:
        return $end_a <=> $end_b;
    }
}

__END__

=head1 OPTIONS

=over 4

=item B<--ref <path to reference multi-fasta file>>

Specify the path to the multi-fasta file that serves as a reference
for the structural variants in the VCF file.

=item B<--outvcf <path to which to write a new VCF-formatted file>>

Specify the path to which to write a new VCF file containing the structural
variants from the input VCF file, but now with tags specifying widened
coordinates

=item B<--refname <string to include as the reference name in the VCF header>>

Specify a string to be written as the reference name in the ##reference line 
of the VCF header.

=item B<--samplename <string to include as the sample name in the "CHROM" line>>

Specify a string to be written as the sample name in the header specifying a 
genotype column in the VCF line beginning with "CHROM".

=item B<--maxsize <maximum size of SV to report>>

Specify an integer for the maximum size of SV to report. 

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
