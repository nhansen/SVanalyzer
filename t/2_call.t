# $Id$
# t/02_call.t - check for accurate results from programs.

use strict;
use Test::More;
use Module::Build;

my $script = 'blib/script/SVrefine.pl';

# Direct useless output to STDERR, to avoid confusing Test::Harness
my $stdin = select STDERR;
# Restore STDOUT as default filehandle
select $stdin;

plan tests => 2;

my $out;
system("perl -w -I lib $script --delta t/test.qdelta --regions t/regions.bed --outvcf t/test1.vcf --ref_fasta t/hs37d5_1start.fa --query_fasta t/utg7180000002239.fa --includeseqs --maxsize 1000000 > t/test1.out 2>&1");
$out = `grep '#' t/test1.vcf | wc -l`;
like $out, qr/^\s*15\s*$/, "$script headerlines";
$out = `awk -F"\t" '\$2==1016054 {print \$8}' t/test1.vcf`;
like $out, qr/REPTYPE=SIMPLEINS/, "$script vartype";
#system("rm t/test1.out");
#system("rm t/test1.vcf");
