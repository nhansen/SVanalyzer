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
        configure_requires      => { 'Module::Build' => '0.36' },
        requires      => { },
        );
if (!(`which samtools`)) {
    warn "ERROR: samtools executable is not found.\nSVanalyzer scripts requires "
        . "the samtools executable to be in the user's PATH\n";
}
if (!(`which edlib-aligner`)) {
    warn "ERROR: Please install the edlib aligner (http://martinsosic.com/edlib/).\nThis script requires "
        . "the edlib-aligner binary to be in the user's PATH\n";
}
if (!(`which nucmer`) || !(`which delta-filter`)) {
    warn "ERROR: Please install the edlib aligner (http://martinsosic.com/edlib/).\nThis script requires "
        . "the edlib-aligner binary to be in the user's PATH\n";
}
$build->create_build_script;

