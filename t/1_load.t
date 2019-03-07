# $Id$
# t/01_load.t - check module loading and program syntax

use strict;
use Test::More;
use Module::Build;

# Four scripts to test at this point:
my @scripts = qw(blib/script/SVrefine blib/script/SVcomp blib/script/SVwiden blib/script/SVmerge blib/script/SVbenchmark);
my @modules = qw(blib/lib/GTB/FASTA.pm blib/lib/GTB/File.pm blib/lib/NHGRI/MUMmer/AlignSet.pm blib/lib/NHGRI/MiniMap2/AlignSet.pm blib/lib/NHGRI/SVanalyzer/Comp.pm blib/lib/NISC/Sequencing/Date.pm );

# Direct useless output to STDERR, to avoid confusing Test::Harness
my $stdin = select STDERR;
# Restore STDOUT as default filehandle
select $stdin;

plan tests => scalar(@modules) * 2 + scalar(@scripts) * 2;

my $out;
for my $script (@scripts) {
    $out = `perl -cw -I lib $script 2>&1`;
    print $out;
    like $out, qr/syntax OK/, "$script syntax";
    $out = `podchecker $script 2>&1`;
    print $out;
    like $out, qr/pod syntax OK|does not contain any pod/, "$script POD";
}

for my $mod (@modules) {
    $out = `perl -I lib -cw $mod 2>&1`;
    print $out;
    like $out, qr/syntax OK/, "$mod syntax";
    $out = `podchecker $mod 2>&1`;
    print $out;
    like $out, qr/pod syntax OK/, "$mod POD";
}

