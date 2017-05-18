# $Id$
package GTB::File;

use strict;
use Carp qw(croak);
use GTB::Unix qw(which);
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(Open first_is_older wait_for_file);
our %EXPORT_TAGS = ('use_bgzip' => [] );
our $Gzip = 'gzip';
our $GzipIn;
our $GzipOut;

our $VERSION = '0.23';

sub import {
    my ($self, @in) = @_;
    my @in_pragmas = grep { /^:/ } @in;
    for my $pragma (@in_pragmas) {
        if ($pragma eq ':use_bgzip') {
            use_bgzip();
        }
        else {
            warn "GTB::File pragma '$pragma' is unrecognized\n";
        }
    }
    GTB::File->export_to_level(1,@_);
}

sub use_bgzip {
    my ($bgzip) = @_;
    if (!$bgzip) {
        $bgzip = which("bgzip");
      
    }
    if ($bgzip) {
        $GzipOut = $bgzip;
    }
}

sub Open {
    my ($file, $mode) = @_;
    if (!$file && $file ne '0') {
        croak "Open: no filename provided";
    }
    if ($mode) {
        $mode = lc $mode;
    }
    elsif ($file =~ /^\s*\|/) {
        $mode = 'w';
    }
    else {
        $mode = 'r';
    }
    my $fh;
    if ($file =~ /\|/) {
        if ($mode eq 'r') {
            if ($file =~ /\|\s*$/) {
                open $fh, $file or die "Can't open pipe '$file', $!\n";
            }
            else {
                croak "To open pipe for reading, pipe character must "
                    . "appear at end of command";
            }
        }
        elsif ($mode eq 'w') {
            if ($file =~ /^\s*\|/) {
                open $fh, $file or die "Can't open pipe '$file', $!\n";
            }
            else {
                croak "To open pipe for writing, pipe character must "
                    . "appear at beginning of command";
            }
        }
        else { # pipe, but not first or last in sequence
            croak << "END_MSG";
If a pipe character is present in the open string, there must be a pipe at
the beginning or end of the string, depending upon whether you plan to
write or read to the filehandle; '$file' is not valid.  If you need to read
and write to a program, try IPC::Open2 or IPC::Open3.
END_MSG
        }
    }
    elsif ($file =~ /\.(b?gz|bz2|zip|Z)$/) {
        if ($mode eq 'r') {
            my $prog = $1 eq 'bz2' ? 'bzip2' : ($GzipIn || $Gzip);
            croak "File ($file) not found" unless (-e $file);
            croak "File ($file) was not readable" unless (-r $file);
            open $fh, "$prog -dc $file |"
                or die "Can't read $file with $prog, $!\n";
        }
        elsif ($mode eq 'w') {
            my $prog = $1 eq 'bz2' ? 'bzip2' : ($GzipOut || $Gzip);
            open $fh, "| $prog > $file"
                or die "Can't create $prog file $file, $!\n";
        }
        elsif ($mode eq 'a') {
            if ($1 eq 'bz2') {
                croak "Open: mode 'a' not supported for bzip2 file $file";
            }
            my $prog = $GzipOut || $Gzip;
            open $fh, "| $prog >> $file"
                or die "Can't append $prog output to $file, $!\n";
        }
        else {
            croak "Open: mode '$mode' not supported; use 'r', 'w' or 'a'";
        }
    }
    elsif ($file eq '-') {
        if ($mode eq 'r') {
            open $fh, '-' or die "Can't read from STDIN, $!\n";
        }
        elsif ($mode eq 'w' || $mode eq 'a') {
            open $fh, '>-' or die "Can't write to STDOUT, $!\n";
        }
        else {
            croak "Open: mode '$mode' not supported; use 'r', 'w' or 'a'";
        }
    }
    elsif ($mode eq 'r') {
        open $fh, '<', $file or die "Can't open $file, $!\n";
    }
    elsif ($mode eq 'w') {
        open $fh, '>', $file or die "Can't create $file, $!\n";
    }
    elsif ($mode eq 'a') {
        open $fh, '>>', $file or die "Can't append to $file, $!\n";
    }
    else {
        croak "Open: mode '$mode' not supported; use 'r', 'w' or 'a'";
    }
    return $fh;
}

sub first_is_older {
    my ($file1, @files) = @_;
    if (!-e $file1) {
        return -1;
    }
    my $m1 = -M _;
    for (my $i = 0; $i < @files; ++$i) {
        if (-e $files[$i]) {
            return $i+1 if $m1 > -M _;
        }
    }
    return 0;
}

sub wait_for_file {
    my ($file, $num_sec_to_wait) = @_;
    if (!$num_sec_to_wait) {
        $num_sec_to_wait = 60;
    }
    if ($num_sec_to_wait < 1 || $num_sec_to_wait !~ /^[0-9]+$/) {
        croak "You should specify a positive number of seconds to wait for the file.\n";
    }
    my $sleep_count = 0;
    my $min_size = $file =~ /gz|bz2$/ ? 20 : 0;
    my $size = -s $file || 0;
    while ($size <= $min_size) {
        if ($sleep_count < $num_sec_to_wait){
            sleep 1;
            ++$sleep_count;
        }
        else {
            die "File $file not found in $num_sec_to_wait seconds.\n";
        }
        $size = -s $file || 0;
    }
    return 1;
}

1;
__END__

=head1 NAME

GTB::File - common file-related operations

=head1 SYNOPSIS

    use GTB::File qw(Open first_is_older);
    my $infile  = '/path/to/infile.gz';
    my $outfile = '/path/to/outfile.gz';
    if (first_is_older($outfile, $infile)) {
        my $fh_in  = Open($infile);
        my $fh_out = Open($outfile, 'w');
        while (<$fh_in>) {
            print $fh_out $_;
        }
    }

=head1 DESCRIPTION

This is a library of exportable functions.

=head1 EXPORTABLE FUNCTIONS

=head2 Open($filename [,'w'])

Returns a file handle to the given file, or dies if file can't be opened.
If the filename ends with .gz, it is assumed to be a gzipped file.  if the
filename ends with .bz2, it is assumed to be a bzip2 compressed file.

If the optional 'w' mode is specified, the file is created/overwritten, and an
output filehandle is returned.  This allows gzipped files to be created
directly, rather than writing them uncompressed first. Append 'a' mode also
works for uncompressed and gzipped files.

  my $fh = Open("/path/to/file.gz");
  while (<$fh>) {
      # do stuff
  }

  my $ofh = Open("/path/to/output.bz2", "w");
  print $ofh "Something long\n";
  close $ofh or die "Error closing file, $!";

As of version 0.11, Open() handles one-way pipe commands as well as ordinary
files and STDIN/STDOUT (the latter two represented by '-').
Does not use the '|-' or '-|' syntax.

 my $stdin = Open('-', 'r'); # opens stdin (buffered)
 my $sh = Open('|/bin/sh', 'w');  # open a stream to stdin of program sh
 while(<$stdin>) {
     print $sh "echo $_";
 }
 close($stdin);
 close($sh);

As of version 0.12, Open() will use the program specified in
$GTB::File::GzipIn to decompress and $GTB::File::GzipOut to compress files
with extensions .gz, .Z, .zip.  $GTB::File::Gzip sets the default for both
compression and decompression.

As of version 0.23, there is now a helper method to switch the output
compression method for gzip files to use bgzip, if it is available.  Before
opening the file, call

    GTB::File::use_bgzip([optional path to bgzip]).

As a bit of syntactic sugar, one can also specify this in the import
statement:

    use GTB::File qw(:use_bgzip Open);

=head2 first_is_older

Returns true when the first file is at least one second older than any of the
others.  Frequently used when updating files (a la make), but be careful in
cases where files are changed rapidly or filesystem timestamps are
unreliable.  The value returned is the index number of the first file that is
newer, starting with the value 1.

  if (first_is_older("my.out", "m1.in", "m2.in")) {
      # recreate my.out from m1.in and m2.in
  }

Also returns true (-1) when the first file doesn't exist at all, because it
clearly needs to be recreated in this case.

This method was added in version 0.21.

=head2 wait_for_file

Checks once per second to see if file exists (and is not empty). User can
specify the number of seconds to wait for the file (default is 60 seconds).          Dies if non-empty file is not found.
usage: wait_for_file($file_name, $num_seconds_to_wait)

This method was added in version .22.

=head1 AUTHORS

 Peter Chines <pchines@mail.nih.gov>

=head1 LEGAL

This software is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use. 

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S. Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose. 

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
Genome Technology Branch" as the citation. 

=cut
