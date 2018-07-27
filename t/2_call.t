# $Id$
# t/02_call.t - check for accurate results from programs.

use strict;
use Test::More;
use Module::Build;

# Direct useless output to STDERR, to avoid confusing Test::Harness
#my $stdin = select STDERR;
# Restore STDOUT as default filehandle
#select $stdin;

my $has_samtools = `which samtools 2>/dev/null`;
my $has_edlib = `which edlib-aligner 2>/dev/null`;
my $has_nucmer = `which nucmer 2>/dev/null`;
my $has_delta_filter = `which delta-filter 2>/dev/null`;

plan tests => 8;

my $out;
# Test SVrefine: (2 tests--requires samtools)
my $script = 'blib/script/SVrefine.pl';

SKIP: {
    ($has_samtools) or skip "Skipping SVrefine tests because no samtools in path", 3;

    system("perl -w -I lib $script --delta t/refine.qdelta --regions t/regions.bed --outvcf t/refined.vcf --ref_fasta t/hs37d5_1start.fa --query_fasta t/utg7180000002239.fa --includeseqs --maxsize 1000000 > t/refine.out 2>&1");
    $out = `grep '#' t/refined.vcf | wc -l`;
    like $out, qr/^\s*13\s*$/, "$script headerlines";
    $out = `awk -F"\t" '\$2==1016054 {print \$8}' t/refined.vcf`;
    like $out, qr/REPTYPE=SIMPLEINS/, "$script vartype";
    system("perl -w -I lib $script --delta t/refine.qdelta --regions t/regions.bed --outvcf t/refined.vcf --ref_fasta t/hs37d5_1start.fa --query_fasta t/utg7180000002239.fa --svregions t/refine.sv.bed --includeseqs > t/refine.out 2>&1");
    $out = `wc -l t/refine.sv.bed`;
    like $out, qr/^\s*3\s/, "$script svregions";
    $out = `awk -F"\t" '\$2==1074450 {print \$6}' t/refine.sv.bed`;
    like $out, qr/REPTYPE=CONTRAC/, "$script svbedregions";
}

SKIP: {
    # Test SVwiden:
    ($has_samtools && $has_nucmer && $has_delta_filter) or skip "Skipping SVwiden tests because one of samtools, nucmer or delta-filter in path", 2;
    $script = 'blib/script/SVwiden.pl';
    
    mkdir "t/test";
    system("perl -w -I lib $script --variants t/widen.vcf --prefix t/widened --ref t/hs37d5_1start.fa --workdir t/test > t/test2.out 2>&1");
    $out = `awk -F"\t" '\$2==821604 {print \$8}' t/widened.vcf`;
    like $out, qr/REPTYPE=DUP/, "$script vartype";
    $out = `awk -F"\t" '\$2==842057 {print \$8}' t/widened.vcf`;
    like $out, qr/REFWIDENED=1:842056-842090/, "$script widensmall";
    #system("rm t/widened.vcf");
    #system("rm t/test2.out");
}

SKIP: {
    # Test SVmerge:
    ($has_samtools && $has_edlib) or skip "Skipping SVmerge tests because one of samtools or edlib-aligner is missing from path", 2;
    $script = 'blib/script/SVmerge.pl';
    
    mkdir "t/test";
    system("perl -w -I lib $script --variants t/merge.vcf --outvcf t/merged.vcf --ref t/hs37d5_1start.fa --workdir t/test > t/test3.out 2>&1");
    $out = `awk -F"\t" '\$2==821604 {print \$8}' t/widened.vcf`;
    like $out, qr/REPTYPE=DUP/, "$script vartype";
    $out = `awk -F"\t" '\$2==842057 {print \$8}' t/widened.vcf`;
    like $out, qr/REFWIDENED=1:842056-842090/, "$script widensmall";
    #system("rm t/merged.vcf");
    #system("rm t/test3.out");
}
