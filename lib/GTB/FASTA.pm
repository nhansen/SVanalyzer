# $Id$
package GTB::FASTA;

use strict;
use Carp qw(croak);
use GTB::File qw(Open);
use File::Spec;

# revcomp(), alleles_to_iupac(), and iupac_to_alleles() are useful
# independently of a FASTA file
use Exporter; 
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(revcomp alleles_to_iupac iupac_to_alleles);


our $VERSION = '0.05';
our $SAMTOOLS = 'samtools';
our %I2A = (
        N => 'ACGT',
        B => 'CGT',
        D => 'AGT',
        H => 'ACT',
        V => 'ACG',
        K => 'GT',
        M => 'AC',
        R => 'AG',
        S => 'CG',
        W => 'AT',
        Y => 'CT',
        A => 'A',
        C => 'C',
        G => 'G',
        T => 'T',
        U => 'U',
);
our %A2I = reverse %I2A;
our $EMPTY = q{};

sub new {
    my ($pkg, $param, $dir) = @_;
    if (!$param) {
        croak "new: assembly name or .mfa file expected";
    }
    my $self = bless { _cached_chr => $EMPTY }, ref $pkg || $pkg;
    my $mfa = $self->find_reference_fasta($param, $dir);
    return $self;
}

sub find_reference_fasta {
    my ($self, $param, $dir) = @_;
    my $mfa;
    my @dirs;
    if (!$param) {
        $mfa = $self->{_mfa};
    }
    elsif ($param =~ m{[/\\]}) {
        if (-f $param) {
            $mfa = $param;
        }
        elsif (-f "$param.mfa") {
            $mfa = $param;
        }
        else {
            croak "Assembly '$param' is not a multi-fasta file";
        }
    }
    elsif ($dir && ref $dir) {
        @dirs = @$dir;
    }
    elsif ($dir) {
        push @dirs, $dir;
    }
    if (!$mfa && $ENV{GENOME_FASTA} && -d $ENV{GENOME_FASTA}) {
        push @dirs, File::Spec->catdir($ENV{GENOME_FASTA},$param),
             $ENV{GENOME_FASTA};
    }
    if (!$mfa) {
        for my $d (@dirs) {
            next if !-d $d;
            my $file = File::Spec->join($d, "$param.mfa");
            for my $suf ('', '.gz', '.rz') {
                if (-f "$file$suf") {
                   $mfa = "$file$suf";
                   last;
                }
            }
        }
    }
    if (!$mfa) {
        if (@dirs) {
            croak "Assembly '$param.mfa' not found in GENOME_FASTA "
                . "directory path\ni.e. @dirs\n";
        }
        else {
            croak "Environment variable GENOME_FASTA is not set, and"
                . " '$param' does not appear\nto be the full path to a"
                . " multi-FASTA file";
        }
    }
    if (!-r $mfa) {
        croak "Assembly '$mfa' is not readable";
    }
    if (ref $self) {
        $self->{_mfa} = $mfa;
        delete $self->{_index};
        delete $self->{_fh};
        $self->{_use_samtools} = ($mfa =~ /\.[gr]z$/);
    }
    return $mfa;
}

sub _get_index {
    my ($self) = @_;
    if (!$self->{_index}) {
        my $fai = "$self->{_mfa}.fai";
        if (!-f $fai) {
            if (!$self->reindex()) {
                die "Failed trying to create FASTA index $fai, $!\n";
            }
        }
        my $fh = Open("$self->{_mfa}.fai");
        $self->{_order} = [];
        while (<$fh>) {
            chomp;
            my @d = split /\t/;
            $self->{_index}{$d[0]} = {
                len        => $d[1],
                offset     => $d[2],
                bases_line => $d[3],
                chars_line => $d[4],
            };
            push @{ $self->{_order} }, $d[0];
        }
    }
    #return
    $self->{_index};
}

sub ids {
    my ($self) = @_;
    $self->_get_index();
    return @{ $self->{_order} };
}

sub len {
    my ($self, $chr) = @_;
    croak "len: chromosome/contig name required" if !defined($chr);
    my $rh_len = $self->_get_index();
    return $rh_len->{$chr}{len};
}

# use length() as a synonym for len(), to maintain compatibility with
# Bio::DB::Fasta
{
    no strict 'refs';
    *length = \&len;
}

sub seq {
    my ($self, $chr, $start, $end) = @_;
    my $rev_comp;
    if ($chr =~ /^([^:\s]+):(\d+)(?:(-)|-(\d+))?$/) {
        if ($start || $end) {
            croak "seq: supply start/end either as separate arguments or in"
                . " single 'chr:start-end'\nregion spec, but not both";
        }
        ($chr, $start) = ($1, $2);
        if ($3 && !$4) {
            $end = "END";
        }
        else {
            $end = $4 || $2;
        }
    }
    my $rh_seq = $self->_get_index();
    if (!$rh_seq->{$chr}) {
        # if can't find with chr prefix, try without, and vice-versa
        $chr =~ s/^chr// or $chr = "chr$chr";
        if (!$rh_seq->{$chr}) {
            $chr =~ s/^chr// or $chr = "chr$chr";
            croak "seq: chromosome '$chr' is not found in reference "
                . "'$self->{_mfa}'";
        }
    }
    if ($end && ($end eq "END" || $end > $rh_seq->{$chr}{len})) {
        $end = $rh_seq->{$chr}{len};
    }
    if ($start && $end && $end < $start) {
        ($start, $end) =  ($end, $start);
        $rev_comp = 1;
    }
    elsif ($start && !$end) {
        $end = $start;
    }
    $start ||= 1;
    $end   ||= $rh_seq->{$chr}{len};
    if ($start > $rh_seq->{$chr}{len}) {
        warn "Attempt to access position $start beyond end of chrom $chr "
            . "($rh_seq->{$chr}{len})\n";
        return $EMPTY;
    }
    my $rs_seq = $self->_read_seq($chr, $start, $end);
    if ($rev_comp) {
        my $seq = $self->revcomp($$rs_seq);
        $rs_seq = \$seq;
    }
    # ASSERT: sequence length should match request
    if (CORE::length($$rs_seq) != $end - $start + 1) {
        die "Assert failure: sequence $chr:$start-$end is wrong length";
    }
    #return
    $$rs_seq;
}

sub reindex {
    my ($self) = @_;
    delete $self->{_index};
    delete $self->{_order};
    my $rc = system $SAMTOOLS, 'faidx', $self->{_mfa};
    if ($rc) {
        return 0;
    }
    return 1;
}

sub revcomp {
    my($self, $seq);
    if ($_[0]->isa('GTB::FASTA')) {
        ($self, $seq) = @_;
    } else {
        $seq = $_[0];
    }
    $seq = reverse $seq;
    my $n = ( $seq =~ tr{ACGTUBDHVMKYRNXSWacgtubdhvmkyrnxsw[]/-}
                        {TGCAAVHDBKMRYNXSWtgcaavhdbkmrynxsw][/-} );
    if ($n != CORE::length($seq)) {
        croak("revcomp: non-sequence characters in sequence string, now '$seq'");
    }
    return $seq;
}

sub alleles_to_iupac {
    my ($self, $all);
    if ($_[0]->isa('GTB::FASTA')) {
        ($self, $all) = @_;
    } else { # handle being called directly
        $all = $_[0];
    }
    
    my %all;
    if (!ref $all) {
        my @a;
        if ($all =~ m{/}) {
            @a = split '/', uc $all;
        }
        else {
            @a = split $EMPTY, uc $all;
        }
        $all = \@a;
    }
    for my $a (@$all) {
        if ($a !~ /^[ACGT]$/i) {
            return 'N';
        }
        else {
            $all{uc $a} = undef;
        }
    }
    my $lookup = join $EMPTY, sort keys %all;
    return $A2I{$lookup};
}

sub iupac_to_alleles {
    my($self, $iupac);
    if ($_[0]->isa('GTB::FASTA')) {
        ($self, $iupac) = @_;
    } else {
        $iupac = $_[0];
    }

    if ($iupac !~ /^[ACGTUBDHVMKYRNXSW]$/i) {
        croak "iupac_to_alleles: '$iupac' is not a valid DNA IUPAC code";
    }
    $iupac = uc $iupac;
    my $a = $I2A{$iupac} || $iupac;
    if (wantarray) {
        my @a = split $EMPTY, $a;
        return @a;
    }
    return $a;
}

sub use_buffer {
    my ($self) = shift;
    warn "use_buffer() is deprecated, and unnecessary\n";
    if ($self->{_use_samtools}) {
        warn "for fastest access, use an uncompressed .mfa file\n";
    }
}

sub _read_seq {
    my ($self, $chr, $start, $end) = @_;
    return if $chr =~ /^\*/;
    my $seq;
    if ($self->{_cached_chr} eq $chr
            && $self->{_cached_start} < $start
            && $self->{_cached_end}   > $end) {
        --$start;
        $seq = substr($self->{_cached_seq},
                      $start - $self->{_cached_start},
                      $end - $start);
    }
    elsif ($self->{_use_samtools}) {
        my $fh = Open("$SAMTOOLS faidx $self->{_mfa} $chr:$start-$end |");
        $_ = <$fh>; # defline
        while (<$fh>) {
            chomp;
            if (/[^ACGTUBDHVMKYRWSNXacgtubdhvmkyrwsnx]/) {
                warn "Sequence '$chr:$start-$end' in $self->{_mfa} "
                    . "contains non-sequence character: $_\n";
            }
            $seq .= $_;
        }
        $self->{_cached_seq}   = $seq;
        $self->{_cached_chr}   = $chr;
        $self->{_cached_start} = $start - 1;
        $self->{_cached_end}   = $end;
    }
    else {  # not compressed, can use normal file ops
        $self->{_fh} ||= Open($self->{_mfa});
        my $rh_index = $self->_get_index();
        my $rh = $rh_index->{$chr};
        my $fh = $self->{_fh};
        --$start;   # to make position zero-based;
        my $start_line = int($start / $rh->{bases_line});
        my $end_line   = int($end / $rh->{bases_line});
        seek $fh, $start_line * $rh->{chars_line} + $rh->{offset}, 0;
        my $s = <$fh>;
        $s =~ s/\s+//g;
        my $pos = $start - $start_line * $rh->{bases_line};
        if ($end_line > $start_line) {
            $seq = substr $s, $pos;
            for (my $line = $start_line + 1; $line < $end_line; ++$line) {
                $s = <$fh>;
                $s =~ s/\s+//g;
                $seq .= $s;
            }
            $s = <$fh>;
            $s =~ s/\s+//g;
            $self->{_cached_seq}   = $seq . $s;
            $self->{_cached_chr}   = $chr;
            $self->{_cached_start} = $start;
            $self->{_cached_end}   = ($end_line+1) * $rh->{bases_line};
            $seq .= substr $s, 0, $end - $end_line * $rh->{bases_line};
        }
        else {
            $self->{_cached_seq}   = $s;
            $self->{_cached_chr}   = $chr;
            $self->{_cached_start} = $start_line * $rh->{bases_line};
            $self->{_cached_end}   = ($end_line+1) * $rh->{bases_line};
            $seq = substr $s, $pos, $end - $start;
        }
    }
    #return
    \$seq;
}

1;
__END__

=head1 NAME

GTB::FASTA - extract sequence from FASTA files

=head1 SYNOPSIS

    use GTB::FASTA;
    my $fdb = GTB::FASTA->new('hg18');
    my $seq = $fdb->seq('chr10', 123456, 243567);
    my $len = $fdb->len('chrX');

=head1 DESCRIPTION

This is an object for extracting sequence from multi-sequence FASTA files
(.mfa).  Requires that samtools be installed in the users path, and that FASTA
files be indexed, or the directories where files exist are writable, allowing
such indexes to be created.

To use a non-standard version of samtools, set $GTB::FASTA::SAMTOOLS to point
to the executable.

=head1 METHODS

=head2 new

Create new FASTA object.

  Usage: my $fdb = GTB::FASTA->new('hg18');
     or: my $fdb = GTB::FASTA->new('/full/path/to/hg19.mfa');
 Return: object

If only assembly name is given, the .mfa file is assumed to be
"$ENV{GENOME_FASTA}/$assembly/$assembly.mfa".  To specify the full path to a
multi-FASTA file, be sure to include a directory separator, e.g.
"./genome.mfa".

=head2 find_reference_fasta

Returns full path to reference multi-FASTA file, identified by searching the
path specified by environment variable GENOME_FASTA, or an optional path.

  Usage: my $mfa = GTB::FASTA->find_reference_fasta("hg18");
     or: my $mfa = GTB::FASTA->find_reference_fasta("hg18", "/my/path");
     or: my $mfa = GTB::FASTA->find_reference_fasta("hg18", \@dirs);
     or: my $mfa = $fdb->find_reference_fasta();
 Return: full path to MFA file
 Except: croaks if no such file is found

=head2 ids

Return names of chromosomes/contigs.

  Usage: my @chroms = $fdb->ids();
 Return: array of sequence names
 Except: dies if .fai index does not exist and can't be created

=head2 len

Return length of chromosome.

  Usage: my $len = $fdb->len('chr1');
 Return: length of sequence, in bases
 Except: dies if .fai index does not exist and can't be created

=head2 seq

Extract sequence from FASTA file.  Coordinates are assumed to be 1-based.  If
end is less than start position, sequence will be reverse-complemented.  If end
position is missing, the single base at the start position is returned.

  Usage: my $seq = $fdb->seq('chr1');   # whole chromosome
     or: my $seq = $fdb->seq('chr1:123456-234567');
     or: my $seq = $fdb->seq('chr1', 123456, 234567);
     or: my $rev = $fdb->seq('chr1', 234567, 123456);
     or: my $base = $fdb->seq('chr1', 123456);
 Return: sequence string

=head2 reindex

Create or re-create .fai index file

  Usage: $fdb->reindex();
 Return: true if success

=head2 revcomp

Helper method to reverse complement a sequence, including handling of IUPAC
ambiguity codes.  Will die if non-sequence characters are encountered (IUPAC
codes as well as additional characters "[]/-" are allowed, since they
are often found in sequences of variants).

  Usage: my $rc   = GTB::FASTA->revcomp($seq);
 Return: reverse-complemented sequence
 Except: dies if non-sequence characters are encountered

=head2 alleles_to_iupac

Helper method to convert alleles to IUPAC code.

  Usage: my $iupac = GTB::FASTA->alleles_to_iupac("TC"); # returns "Y"
 Return: IUPAC code for DNA [ABCDGHKMNRSWY]

=head2 iupac_to_alleles

Helper method to convert IUPAC code to alleles.

  Usage: my $all = GTB::FASTA->iupac_to_alleles("W"); # returns "AT"
         my @all = GTB::FASTA->iupac_to_alleles("W"); # returns ("A","T")
 Return: string or array of alleles, depending on calling context

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
