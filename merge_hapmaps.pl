#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $hmp_file_pattern;
my $help;

parse_args();
parse_hmp();

exit(0);


sub parse_hmp {
	my $print_header = 1;

	foreach my $hmp_file (glob($hmp_file_pattern)) {
		open(HMP, '<', $hmp_file) or error("can't read $hmp_file: $!");

		while (my $line = <HMP>) {
			if ($line =~ /^rs#\t/) {
				if ($print_header == 1) {
					print(STDOUT "$line");

					$print_header = 0;
				}

				next();
			}

			print(STDOUT "$line");
		}

		close(HMP);
	}

	return(0);
}


sub error {
	my $msg = shift();

	if (defined($msg)) {
		print(STDERR "error: $msg\n");
	}

	exit(0);
}


sub arg_error {
	my $msg = shift();

	if (defined($msg)) {
		print(STDERR "error: $msg\n\n");
	}

	print("use option -h or --help to display help menu\n");

	exit(0);
}


sub parse_args {
	if ($#ARGV == -1) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	GetOptions ('p|pattern=s' => \$hmp_file_pattern,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	if (! defined($hmp_file_pattern)) {
		arg_error('hapmap file pattern required');
	}

	return(0);
}


__END__

=head1 NAME

merge_hapmaps.pl

=head1 SYNOPSIS

merge_hapmaps.pl -p "pop.chr*.hmp.txt"

=head1 DESCRIPTION

merge_hapmaps.pl is a simple tool to combine individual hapmaps into a merged file

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -p --pattern  hapmap file pattern (required)

 -h --help     display help menu

=cut
