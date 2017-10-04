# $Id:$
# To create build script, run, "perl Build.PL"
use Module::Build;

my $build = Module::Build->new(
        module_name   => 'SVanalyzer',
        dist_version  => 0.10,
        dist_author   => 'NHGRI Perl Programmers',
        dist_abstract => 'Package for analyzing structural variation in sequence data',
        license       => 'perl',
        scripts       => [qw(scripts/SVrefine.pl scripts/SVcomp.pl scripts/SVwiden.pl)],
        requires      => { },
        );
if (!`which samtools`) {
    warn "ERROR: samtools executable is not found.\nThis script requires "
        . "the samtools binary to be in the user's PATH\n";
}
$build->create_build_script;

