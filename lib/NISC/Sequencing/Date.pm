#!/usr/local/bin/perl -w 
############################################################
# Date.pm: Module for parsing and outputting dates
#
# Author:       Nancy F. Hansen
# Version: $Id: Date.pm,v 1.9 2005/03/10 17:42:40 nhansen Exp $
############################################################

package NISC::Sequencing::Date;
use strict;
use Carp;
use Time::Local;

###############################################################################
# Date constructor
#
# INPUT: parameters:
#    One of the following must be specified:
#    -consed - the date in 6 digit "consed" format: i.e., '010427'.
#    -oracle - the date in Oracle format: i.e., '04-FEB-2001'.
#    -month_day_slash - the date in slashed format: i.e., '4/6/00',or '4/6/2000'
#    -plain_language - 'today', 'tomorrow', 'yesterday', etc.
#    -sixdigit - six digits in order mmddyy
#    -time_local - seconds since the "Epoch"
#    -barcode_style - style on NISC barcodes
#    -trace_archive - style reported by NCBI (e.g., Feb 27 2004 12:00AM)
#    -formal - full name of month day, year (e.g., February 27, 2003)
#
# OUTPUT: a new Date object
###############################################################################

sub new {

    my $this = shift;
    my %params = @_;
 
    my ($day, $month, $year, $hour, $minute, $second);
    my $consed = $params{-consed};
    my $oracle = $params{-oracle};
    my $month_day_slash = $params{-month_day_slash};
    my $plain_language = $params{-plain_language};
    my $barcode_style = $params{-barcode_style};
    my $sixdigit = $params{-sixdigit}; # mmddyy
    my $time_local = $params{-time_local};
    my $trace_archive = $params{-trace_archive};
    my $formal = $params{-formal};
 
    if (defined($consed) ) # consed format in
    {
        if ($consed =~ /(\d{2})(\d{2})(\d{2})(:\d{6}){0,1}/) # date with or without time
        {
            ($year, $month, $day) = ($1, $2, $3);
            $consed = "$year$month$day";
            if ($year < 80)
            {
                $year += 2000;
            }
            else
            {
                $year += 1900;
            }
        }
        else
        {
            croak "Unexpected consed date formatting in $consed!\n";
        }
    }
    elsif (defined ($oracle) ) # oracle format in
    {
        if ($oracle =~ /(\d{1,2})\-([A-Z]{3})\-(\d+)(,.*){0,1}$/)
        {
            $day = $1;
            $month = month_name_number($2);
            $year = $3;

            if ($year < 100)
            {
                if ($year > 80)
                {
                    $year += 1900;
                }
                else
                {
                    $year += 2000;
                }
            }
        }
        else
        {
            croak "Unexpected oracle date formatting in $oracle!\n";
        }
    } 
    elsif (defined ($month_day_slash) ) # month/day slash format in
    {
        if ($month_day_slash =~ m:(\d{1,2})/(\d{1,2})/(\d+):)
        {
            $month = $1;
            $day = $2;
            $year = $3;
 
            if ($year < 80) # two-digit after 2000
            {
                 $year += 2000;
            }
            elsif ($year < 100) # two-digit before 2000
            {
                $year += 1900;
            }
        }
        else
        {
            croak "Unexpected month/day slash date formatting in $month_day_slash!\n";
        }
    }
    elsif (defined ($barcode_style) ) # month/day slash format in
    {
        if ($barcode_style =~ m:(\d{2})/([a-zA-Z]+)/(\d+):)
        {
            $day = $1;
            my $month_name = $2;
            $month_name =~ tr/a-z/A-Z/;
            $month = month_name_number($month_name);
            $year = $3;
 
            if ($year < 80) # two-digit after 2000
            {
                 $year += 2000;
            }
            elsif ($year < 100) # two-digit before 2000
            {
                $year += 1900;
            }
        }
        else
        {
            croak "Unexpected barcode style date formatting in $barcode_style!\n";
        }
    }
    elsif (defined ($plain_language))
    {
        my $seconds;
        if ($plain_language eq 'today')
        {
            $seconds = time;
        }
        elsif ($plain_language eq 'tomorrow')
        {
            $seconds = time + 24*60*60;
        }
        elsif ($plain_language eq 'yesterday')
        {
            $seconds = time - 24*60*60;
        }
 
        my ($thissec, $thismin, $thishour, $mday, $mon, $yr, $wday, $yday, $idsdst) = localtime($seconds);
        $month = $mon + 1;
        $day = $mday;
        $year = $yr + 1900;

        $hour = $thishour;
        $minute = $thismin;
        $second = $thissec;
    }
    elsif (defined ($sixdigit))
    {
        if ($sixdigit =~ /^(\d{2})(\d{2})(\d{2})$/)
        {
            $month = $1;
            $day = $2;
            $year = $3;
    
            if ($year > 50)
            {
                $year += 1900;
            }
            else
            {
                $year += 2000;
            }
        }
        else
        {
            croak "Date passed in -sixdigit to Date constructor is not in six digit format!\n";
        }
    }
    elsif (defined ($time_local))
    {
        my ($sec, $min, $hour, $mday, $mon, $yr, $wday, $yday, $idsdst) = localtime($time_local);
        $month = $mon + 1;
        $day = $mday;
        $year = $yr + 1900;
    }
    elsif (defined ($trace_archive))
    {
        if ($trace_archive =~ /^(\S+)\s+(\d+)\s+(\d{2,4})\s+(\S+)$/)
        {
            $month = month_name_number($1)
                or croak "Don\'t recognize month $1 in NISC::Sequencing::Date!\n";
            $day = $2;
            $year = $3;
        }
        else
        {
            croak "Date passed in -trace_archive to Date constructor is not in trace_archive format!\n";
        }
    }
    elsif (defined ($formal))
    {
        if ($formal =~ /^(\S+)\s+(\d+),\s*(\d{4})$/)
        {
            $month = formal_month_name_number($1)
                or croak "Don\'t recognize month $1 in NISC::Sequencing::Date!\n";
            $day = $2;
            $year = $3;
        }
        else
        {
            croak "Date passed in -formal to Date constructor is not in formal format!\n";
        }
    }
 
    my $self = { consed => $consed, oracle => $oracle, month_day_slash => $month_day_slash, sixdigit => $sixdigit, day => $day, month => $month, year => $year, time_local => $time_local, formal => $formal, hour => $hour, minute => $minute, second => $second };
    my $class = ref($this) || $this;
    bless $self, $class;
    return($self);
}

###############################################################################
# Subroutine to output date in trace_archive format
#
# INPUT: Date object, 
# OUTPUT: trace-archive-formatted date (e.g., 'Feb 21 1996')
###############################################################################

sub trace_archive{

    my $self = shift;
    my %params = @_;
    my $day = $self->day();
    my $month = $self->three_letter_month();
    my $year = $self->year();

    return "$month $day $year"; 

} ## trace_archive

###############################################################################
# Subroutine to output date in a formal format
#
# INPUT: Date object, -one_digit_day => 1 will leave out "0" in the date.
# OUTPUT: formal-formatted date (e.g., 'February 21, 1996')
###############################################################################

sub formal{

    my $self = shift;
    my %params = @_;
    my $one_digit_day = $params{-one_digit_day};
    my $day = $self->day();
    my $month = $self->month_formal();
    my $year = $self->year();

    if ($one_digit_day)
    {
        $day =~ s:^0+::;
    }

    return "$month $day, $year";

} ## end formal
###############################################################################
# Subroutine to output date in consed format
#
# INPUT: Date object.
# OUTPUT: consed-formatted date (e.g., '010908')
###############################################################################

sub consed{

    my $self = shift;

    if (!defined ($self->{consed}) )
    {
        my $day = $self->day();
        $day = make_two_digits($day);
        my $month = $self->month();
        $month = make_two_digits($month);
        my $year = $self->year();

        # take only last two digits of year:
        $year =~ s:\d{2}(\d{2}):$1:;
        $self->{consed} = "$year$month$day";
    } 

    return $self->{consed};

}

###############################################################################
# Subroutine to output date in sixdigit format
#
# INPUT: Date object, -yearfirst => 1 will print the year in the first
#     two spots.
# OUTPUT: sixdigit-formatted date (e.g., '090901')
###############################################################################

sub sixdigit{

    my $self = shift;
    my %params = @_;
    my $yearfirst = $params{'-yearfirst'};

    if (!defined ($self->{sixdigit}) )
    {
        my $day = $self->day();
        $day = make_two_digits($day);
        my $month = $self->month();
        $month = make_two_digits($month);
        my $year = $self->year();

        # take only last two digits of year:
        $year =~ s:\d{2}(\d{2}):$1:;
        $self->{sixdigit} = "$month$day$year";
        $self->{sixdigityearfirst} = "$year$month$day";
    } 

    return ($yearfirst) ? $self->{sixdigityearfirst} : $self->{sixdigit};
}

###############################################################################
# Subroutine to output date in eightdigit format
#
# INPUT: Date object.
# OUTPUT: eightdigit-formatted date (e.g., '20030410')
###############################################################################

sub eightdigit{

    my $self = shift;

    if (!defined ($self->{eightdigit}) )
    {
        my $day = $self->day();
        $day = make_two_digits($day);
        my $month = $self->month();
        $month = make_two_digits($month);
        my $year = $self->year();

        $self->{eightdigit} = "$year$month$day";
    } 

    return $self->{eightdigit};

} ## end eightdigit

###############################################################################
# Subroutine to output date in format used by Stephen Granite in his
#    "TimeSlug" routine.
#
# INPUT: Date object, -no_time => 1 only returns 'YYYY-MM-DD'.
# OUTPUT: date/time (e.g., '2003-04-10_23.33.22')
###############################################################################

sub time_slug{

    my $self = shift;
    my %params = @_;
    my $no_time = (defined ($params{-no_time})) ? $params{-no_time} : 0;

    my $day = $self->day();
    $day = make_two_digits($day);
    my $month = $self->month();
    $month = make_two_digits($month);
    my $year = $self->year();
    my $hour = $self->hour() || '00';
    $hour = make_two_digits($hour);
    my $min = $self->minute() || '00';
    $min = make_two_digits($min);
    my $sec = $self->second() || '00';
    $sec = make_two_digits($sec);

    $self->{time_slug} = ($no_time) ? 
                         "$year-$month-$day" :
                         "$year-$month-$day\_$hour.$min.$sec";

    return $self->{time_slug};

} ## end time_slug

###############################################################################
# Subroutine to output date in oracle format
#
# INPUT: Date object.
# OUTPUT: oracle-formatted date (e.g., '04-FEB-2001')
###############################################################################

sub oracle{

    my $self = shift;

    if (!defined ($self->{oracle}) )
    {
        my $day = $self->day();
        my $month = $self->month();
        my $year = $self->year();

        # convert month to "named" month:
        $month = make_two_digits($month);
        $month = month_name_number($month);
        $self->{oracle} = "$day-$month-$year";
    } 

    return $self->{oracle};

}

###############################################################################
# Subroutine to computer the seconds since the Epoch (1/1/1970) at noon
#     on a given date.
#
# INPUT: Date object.
# OUTPUT: scalar number (seconds)
###############################################################################

sub time_local{

    my $self = shift;

    unless (defined ($self->{time_local}))
    {
        my $day = $self->day();
        my $month = $self->month();
        my $year = $self->year();
    
        $month--; # to get into localtime format
    
        $self->{time_local} = timelocal(0,0,12,$day,$month,$year);
    }

    return $self->{time_local};

}

###############################################################################
# Subroutine to return a Date object translated from the given date by a
#     fixed period.
#
# INPUT: Date object, -days => n translates n days from the given Date.
# OUTPUT: new Date object
###############################################################################

sub translate {

    my $self = shift;
    my %params = @_;
    my $date;

    if (defined $params{-days})
    {
        my $time_local = $self->time_local();
        my $new_time_local = $time_local + 24*60*60*$params{-days};

        $date = $self->new(-time_local => $new_time_local);
    }

    return $date;
}

###############################################################################
# Subroutine to output date in month/day slash format
#
# INPUT: Date object.
# OUTPUT: month_day_slash-formatted date (e.g., '4/6/2000')
###############################################################################

sub month_day_slash{

    my $self = shift;

    if (!defined ($self->{month_day_slash}) )
    {
        my $day = $self->day();
        my $month = $self->month();
        my $year = $self->year();

        $self->{month_day_slash} = "$month/$day/$year";
    } 

    return $self->{month_day_slash};

}

###############################################################################
# Subroutine to output date in barcode month/day slash format
#
# INPUT: Date object.
# OUTPUT: barcode-style date (e.g., '04/May/2003')
###############################################################################

sub barcode_style{

    my $self = shift;

    if (!defined ($self->{barcode_style}) )
    {
        my $day = $self->day();
        $day =~ s:^(\d)$:0$1:;
        my $month = $self->month();
        my $month_name = month_name_number($month);
        $month_name =~ tr/A-Z/a-z/;
        if ($month_name =~ m:^([a-z]):)
        {
            my ($lc, $uc) = ($1, uc($1));
            $month_name =~ s:^$lc:$uc:;
        }
        my $year = $self->year();

        $self->{barcode_style} = "$day/$month_name/$year";
    } 

    return $self->{barcode_style};

} ## end barcode_style

###############################################################################
# Subroutine to return the year for a date object
#
# INPUT: Date object.
# OUTPUT: 4-digit year
###############################################################################
sub year{

    my $self = shift;
   
    return $self->{year};

}

###############################################################################
# Subroutine to return the month for a date object
#
# INPUT: Date object.
# OUTPUT: month value
###############################################################################
sub month{

    my $self = shift;
   
    return $self->{month};

}

###############################################################################
# Subroutine to return the month name for a date object
#
# INPUT: Date object.
# OUTPUT: full month name (e.g., July, February)
###############################################################################
sub month_formal{

    my $self = shift;
   
    my %formal_month = ('01' => 'January', '02' => 'February', 
                        '03' => 'March', '04' => 'April',
                        '05' => 'May', '06' => 'June',
                        '07' => 'July', '08' => 'August',
                        '09' => 'September', '10' => 'October',
                        '11' => 'November', '12' => 'December');

    my $month = $self->month();
    $month = make_two_digits($month);

    return $formal_month{$month};

} ## end month_formal

###############################################################################
# Subroutine to return the three letter version of the month name for a 
#    date object
#
# INPUT: Date object.
# OUTPUT: three letter month name (e.g., Jul, Feb)
###############################################################################
sub three_letter_month{

    my $self = shift;
   
    my %three_letter_month = ('01' => 'Jan', '02' => 'Feb', 
                        '03' => 'Mar', '04' => 'Apr',
                        '05' => 'May', '06' => 'Jun',
                        '07' => 'Jul', '08' => 'Aug',
                        '09' => 'Sep', '10' => 'Oct',
                        '11' => 'Nov', '12' => 'Dec');

    my $month = $self->month();
    $month = make_two_digits($month);

    return $three_letter_month{$month};

} ## end three_letter_month

###############################################################################
# Subroutine to return the day for a date object
#
# INPUT: Date object.
# OUTPUT: day
###############################################################################
sub day{

    my $self = shift;
   
    return $self->{day};

}

###############################################################################
# Subroutine to return the hour for a date object
#
# INPUT: Date object.
# OUTPUT: hour
###############################################################################
sub hour{

    my $self = shift;
   
    return $self->{hour};

}

###############################################################################
# Subroutine to return the minute for a date object
#
# INPUT: Date object.
# OUTPUT: minute
###############################################################################
sub minute{

    my $self = shift;
   
    return $self->{minute};

}

###############################################################################
# Subroutine to return the second for a date object
#
# INPUT: Date object.
# OUTPUT: second
###############################################################################
sub second{

    my $self = shift;
   
    return $self->{second};

}

###############################################################################
# Subroutine to return a 2-digit month number from a 3-character month
#     name, and vice-versa.
#
# INPUT: Month in one of the two formats (e.g., "MAR" or "03")
# OUTPUT: Month in the other of the two formats (e.g., "03" or "MAR")
###############################################################################

sub month_name_number
{
    my $thing_to_convert = shift @_;

    my %month_number = ('UNK' => '00', 'JAN' => '01', 'FEB' => '02', 'MAR' => '03', 'APR' => '04', 'MAY' => '05', 'JUN' => '06', 'JUL' => '07', 'AUG' => '08', 'SEP' => '09', 'OCT' => '10', 'NOV' => '11', 'DEC' => '12');

    my %month_name = reverse %month_number;

    if ($thing_to_convert =~ m/\d/) # we have a number
    {
        $thing_to_convert =~ s:^(\d)$:0$1:;

        return $month_name{$thing_to_convert};
    }
    else 
    {
        $thing_to_convert =~ tr/a-z/A-Z/;
        return $month_number{$thing_to_convert};
    }
}

###############################################################################
# Subroutine to return a 2-digit month number from a complete month
#     name, and vice-versa.
#
# INPUT: Month in one of the two formats (e.g., "March" or "03")
# OUTPUT: Month in the other of the two formats (e.g., "03" or "March")
###############################################################################

sub formal_month_name_number
{
    my $thing_to_convert = shift @_;

    my %month_number = ('January' => '01', 'February' => '02', 'March' => '03', 'April' => '04', 'May' => '05', 'June' => '06', 'July' => '07', 'August' => '08', 'September' => '09', 'October' => '10', 'November' => '11', 'December' => '12');

    my %month_name = reverse %month_number;

    if ($thing_to_convert =~ m/\d{2}/) # we have a number
    {
        return $month_name{$thing_to_convert};
    }
    else 
    {
        return $month_number{$thing_to_convert};
    }
}

###############################################################################
# Subroutine to return the number of days ago this date was
#     (will be negative if date is in the future)
#
# INPUT: Date object
# OUTPUT: integer number of full days
###############################################################################

sub days_ago
{
    my $self = shift @_;
    my $seconds = $self->time_local();

    my $today = $self->new(-plain_language => 'today');
    my $today_seconds = $today->time_local();

    my $number_of_days = ($today_seconds - $seconds)/(24*60*60);

    return int($number_of_days);

} ## end days_ago

###############################################################################
# Subroutine to return the day number on the NISC calendar (whose
#     first day was 1/7/97, inaccurate though that may be!)
#
# INPUT: Date object
# OUTPUT: integer day number
###############################################################################

sub nisc_day
{
    my $self = shift @_;
    my $end_seconds = $self->time_local();

    my $start_date = NISC::Sequencing::Date->new(-oracle => '07-JAN-1997');
    my $start_seconds = $start_date->time_local();

    my $number_of_days = ($end_seconds - $start_seconds)/(24*60*60);

    return int($number_of_days) + 1;

} ## end nisc_day

###############################################################################
# Subroutine to return the label for a "NISC day" number (which is
#     just the date in month_day_slash format)
#
# INPUT: day number
# OUTPUT: scalar string (formatted date)
###############################################################################

sub nisc_day_label
{
    my $day_number = shift @_;

    my $start_date = NISC::Sequencing::Date->new(-oracle => '07-JAN-1997');
    my $start_seconds = $start_date->time_local();

    my $date_seconds = $start_seconds + ($day_number - 1)*24*60*60;
    my $end_date = NISC::Sequencing::Date->new(-time_local => $date_seconds);

    return $end_date->month_day_slash();

} ## end nisc_day_label

###############################################################################
# Subroutine to return the week number on the NISC calendar (whose
#     first week was 1/7/97-1/13/97, inaccurate though that may be!)
#
# INPUT: Date object, -offset => -2 will move the start of the week
#     back two days, i.e., from Sunday through Saturday
# OUTPUT: integer week number
###############################################################################

sub nisc_week
{
    my $self = shift @_;
    my %params = @_;
    my $end_seconds = $self->time_local();

    my $start_date = NISC::Sequencing::Date->new(-oracle => '07-JAN-1997');

    if (my $offset = $params{-offset})
    {
        $start_date = $start_date->translate(-days => $offset);
    }
 
    my $start_seconds = $start_date->time_local();

    my $number_of_weeks = ($end_seconds - $start_seconds)/(7*24*60*60);

    return int($number_of_weeks) + 1;

} ## end nisc_week

###############################################################################
# Subroutine to return the label for a "NISC week" number (which is
#     just the date of the last day of the week in month_day_slash format
#
# INPUT: week number, -offset => -2 will start each week 2 days earlier,
#     i.e., on a Sunday
# OUTPUT: scalar string (formatted date)
###############################################################################

sub nisc_week_label
{
    my $week_number = shift @_;
    my %params = @_;

    my $start_date = NISC::Sequencing::Date->new(-oracle => '07-JAN-1997');

    if (my $offset = $params{-offset})
    {
        $start_date = $start_date->translate(-days => $offset);
    }

    my $start_seconds = $start_date->time_local();

    my $date_seconds = $start_seconds + $week_number*7*24*60*60;
    my $end_date = NISC::Sequencing::Date->new(-time_local => $date_seconds);

    return $end_date->month_day_slash();

} ## end nisc_week_label

###############################################################################
# Subroutine to return the month number on the NISC calendar (whose
#     first month is 1/97, inaccurate though that may be!)
#
# INPUT: Date objects
# OUTPUT: integer month number
###############################################################################

sub nisc_month
{
    my $self = shift @_;

    my $month = $self->month();
    my $year = $self->year();

    my $no_months = 12*($year - 1997);
    $no_months += $month;

} ## end nisc_month

###############################################################################
# Subroutine to return the label for a "NISC month" number.
#
# INPUT: Month number
# OUTPUT: scalar string (label in 1/97) format
###############################################################################

sub nisc_month_label
{
    my $nisc_month_number = shift @_;

    if (!defined ($nisc_month_number))
    {
        return '';
    }
    my $year = int(($nisc_month_number - 1)/12) + 1997;
    $year =~ s:^(\d{2})::; # strip first two digits
    my $month = $nisc_month_number - 12*int(($nisc_month_number - 1)/12);

    return "$month/$year";

} ## end nisc_month_label

###############################################################################
# Subroutine to return a list of month/year strings for all months from
#     the start_date to present.
#
# INPUT: -start_date specifies the first month to include (default is 12/00)
#        -end_date specifies the last month to include (default is present month)
# OUTPUT: list of strings in format 
###############################################################################

sub get_all_past_months
{

} ## end get_all_past_months

###############################################################################
# Subroutine to sort dates in chronological order
#
# INPUT: Two date objects
# OUTPUT: 1, 0, or -1
###############################################################################

sub sort_by_date
{
    my ($date_a, $date_b) = @_;
    my $year_a = $date_a->year();
    my $year_b = $date_b->year();

    if ($year_a != $year_b)
    {
        return $year_a <=> $year_b;
    }
   
    my $month_a = $date_a->month();
    my $month_b = $date_b->month();

    if ($month_a != $month_b) 
    {
        return $month_a <=> $month_b;
    }

    my $day_a = $date_a->day();
    my $day_b = $date_b->day();
  
    return $day_a <=> $day_b; 
}

sub make_two_digits
{
    my $arg = shift @_;
    if ( length $arg == 1 )
    {
        $arg = "0".$arg;
    }
    return $arg;
}

1;

__END__

=head1 NAME

NISC::Sequencing::Date - Perl extension for parsing and printing dates

=head1 SYNOPSIS

use NISC::Sequencing::Date;

$date_obj = NISC::Sequencing::Date->new( -oracle => '2001-FEB-03');

$consed_date = $date_obj->consed();

=head1 DESCRIPTION

This module is used to parse print, and sort dates in Oracle and consed formats.

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=head2 new

 Title   : new
 Usage   : $seq_chem_obj = NISC::Sequencing::Chemistry->new( -name => 
                  'rhodamine', -label => 'primer');
           $id = $seq_chem_obj->id();

 Function: Returns a new Chemistry object from the niscprod database.
           (Doesn't actually do a database retrieval unless a method
           calling for a new value is called).
 Returns : a new Chemistry object
 Args	 : Either:
           -id          => primary key value from seq_chem table
           Or:
           -name	=> the name of the chemistry     AND
	   -label	=> location of the dye label

=head2 name

 Title   : name
 Usage   : $seq_chem_name = $seq_chem_obj->name();
 
 Function: Returns the name of the sequencing chemistry (e.g., "big-dye", 
           "d-rhoadmine", or "energy-transfer".)
 Returns : a scalar string (name)
 Args	 : none.

=head2 label

 Title   : label
 Usage   : $seq_chem_label = $seq_chem_obj->label();
 
 Function: Returns the location of the dye-label for this chem (i.e., 
           "primer" or "terminator".)
 Returns : a scalar string (label)
 Args	 : none.

=head2 company

 Title   : company
 Usage   : $company = $seq_chem_obj->company();
 
 Function: Returns the company from the "company" field of the 
           seq_chem table.
 Returns : a scalar string (company name)
 Args	 : none.

=head2 id

 Title   : id
 Usage   : $seq_chem_id = $seq_chem_obj->id();
 
 Function: Returns the primary key value for the sequencing chemistry 
           from the seq_chem table.
 Returns : a scalar number (id value)
 Args	 : none.

=head2 _retrieve_seq_chem_info

 Title   : _retrieve_seq_chem_info
 Usage   : $seq_chem_obj = $seq_chem_obj->_retrieve_seq_chem_info();
 
 Function: Populates all fields of the Chemistry object with values
           from the seq_chem table, or, if no record exists for the
           chemistry, carps and assigns 0 or '' to all fields.
 Returns : The new object
 Args	 : none.

=head1 AUTHOR

Nancy F. Hansen <nhansen@nhgri.nih.gov>

=head1 SEE ALSO

NISC::Finishing::Experiment.

=cut

