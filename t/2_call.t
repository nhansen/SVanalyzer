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
system("perl -w -I lib $script --delta t/test.delta --regions t/regions.bed --outvcf t/test1.vcf --ref_fasta t/hs37d5_1start.fa --query_fasta t/utg7180000002239.fa --includeseqs --maxsize 1000000 > t/test1.out 2>&1");
$out = `awk '\$2==11589022 {print \$9}' t/test1.vcf`;
like $out, qr/43/, "$script count";
$out = `awk '\$1==1 {print \$3}' t/test1out/somatic_diffs.vs`;
like $out, qr/11589021/, "$script variant";
system("rm t/test1.out");
system("rm t/test1.vcf");
