# $Id:$
# To create build script, run, "perl Build.PL"
use Module::Build;

my $build = Module::Build->new(
        module_name   => 'SVanalyzer',
        dist_version  => 0.36,
        dist_author   => 'NHGRI Perl Programmers',
        dist_abstract => 'Package for analyzing structural variation in sequence data',
        license       => 'perl',
        scripts       => [qw(scripts/SVrefine scripts/SVcomp scripts/SVwiden scripts/SVmerge scripts/SVbenchmark scripts/svanalyzer)],
        configure_requires      => { 'Module::Build' => '0.36' },
        requires      => { 'Log::Log4perl' => '0.21' },
        );
if (!(`which samtools 2>/dev/null`)) {
    warn "ERROR: samtools executable is not found.\nSVanalyzer scripts requires "
        . "the samtools (http://www.htslib.org) executable to be in the user's PATH\n";
}
if (!(`which edlib-aligner 2>/dev/null`)) {
    warn "WARNING: You do not have the edlib aligner (http://martinsosic.com/edlib/) in your path, so will be unable to run some SVanalyzer commands.\n";
}
if (!(`which nucmer 2>/dev/null`) || !(`which delta-filter 2>/dev/null`)) {
    warn "WARNING: You do not seem to have the MUMmer alignment tools (available at https://github.com/mummer4/mummer) installed.\nSVanalyzer requires "
        . "the executables for nucmer and delta-filter to be in the user's PATH to run some commands.\n";
}
if (!(`which bedtools 2>/dev/null`)) {
    warn "WARNING: You do not have the bedtools (http://bedtools.readthedocs.io/en/latest/) executable in your path, so will be unable to run some SVbenchmark with the --includebed option.\n";
}
$build->create_build_script;

