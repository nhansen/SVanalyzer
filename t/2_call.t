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

my $no_tests = 0;
if ($has_samtools) {
    $no_tests += 2;
}
if ($has_samtools && $has_nucmer && $has_delta_filter) {
    $no_tests += 2;
}
plan tests => $no_tests;

my $out;
# Test SVrefine: (2 tests--requires samtools)
my $script = 'blib/script/SVrefine.pl';
if ($has_samtools) {
    system("perl -w -I lib $script --delta t/refine.qdelta --regions t/regions.bed --outvcf t/refined.vcf --ref_fasta t/hs37d5_1start.fa --query_fasta t/utg7180000002239.fa --includeseqs --maxsize 1000000 > t/refine.out 2>&1");
    $out = `grep '#' t/refined.vcf | wc -l`;
    like $out, qr/^\s*13\s*$/, "$script headerlines";
    $out = `awk -F"\t" '\$2==1016054 {print \$8}' t/refined.vcf`;
    like $out, qr/REPTYPE=SIMPLEINS/, "$script vartype";
    #system("rm t/test1.out");
    #system("rm t/test1.vcf");
}
else {
    print "Skipping SVrefine.pl tests--no samtools in path!\n";
}

# Test SVwiden:
$script = 'blib/script/SVwiden.pl';
if ($has_samtools && $has_nucmer && $has_delta_filter) {
    mkdir "t/test";
    system("perl -w -I lib $script --invcf t/widen.vcf --outvcf t/widened.vcf --ref t/hs37d5_1start.fa --workdir t/test > t/test2.out 2>&1");
    $out = `awk -F"\t" '\$2==821604 {print \$8}' t/widened.vcf`;
    like $out, qr/REPTYPE=DUP/, "$script vartype";
    $out = `awk -F"\t" '\$2==842057 {print \$8}' t/widened.vcf`;
    like $out, qr/REFWIDENED=1:842056-842090/, "$script widensmall";
    #system("rm t/widened.vcf");
    #system("rm t/test2.out");
}
else {
    print "Skipping SVwiden tests--must have samtools, nucmer, and delta-filter in path!\n";
}
