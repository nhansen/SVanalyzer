# $Id$
# t/01_load.t - check module loading and program syntax

use strict;
use Test::More;
use Module::Build;

my @scripts = qw(blib/script/SVrefine.pl);
my @modules = qw();

# Direct useless output to STDERR, to avoid confusing Test::Harness
my $stdin = select STDERR;
# Restore STDOUT as default filehandle
select $stdin;

plan tests => scalar(@scripts) * 2;

my $out;
for my $script (@scripts) {
    $out = `perl -cw -I lib $script 2>&1`;
    print $out;
    like $out, qr/syntax OK/, "$script syntax";
    $out = `podchecker $script 2>&1`;
    print $out;
    like $out, qr/pod syntax OK|does not contain any pod/, "$script POD";
}

