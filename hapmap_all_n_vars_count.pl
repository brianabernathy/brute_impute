#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $hmp_file;
my $help;

parse_args();
parse_hmp();

exit(0);


sub parse_hmp {
	my $all_n_count = 0;

	open(HMP, '<', $hmp_file) or error("can't read $hmp_file: $!");

	while (my $line = <HMP>) {
		chomp($line);

		my ($rs, $alleles, $chrom, $pos, $strand, $assy, $center, $prot_lsid, $array_lsid, $panel_lsid, $qc_code, @samples) = split(/\t/, $line);

		if ($rs eq 'rs#') {
			next();
		}

		my $all_n_rec = 1;

		foreach my $sample (@samples) {
			if ($sample ne 'N') {
				$all_n_rec = 0;

				last();
			}
		}

		if ($all_n_rec == 1) {
			$all_n_count++;
		}
	}

	close(HMP);

	print("variant records containing all Ns: $all_n_count\n");

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

	GetOptions ('hmp=s' => \$hmp_file,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	if (! defined($hmp_file)) {
		arg_error('hapmap file required');
	}

	return(0);
}


__END__

=head1 NAME

hapmap_all_n_vars_count.pl

=head1 SYNOPSIS

hapmap_all_n_vars_count.pl --hmp variants.hmp.txt

=head1 DESCRIPTION

hapmap_all_n_vars_count.pl counts the number of variant records with all 'N' alleles

TASSEL may silently write all 'N' records when it encounters issues.

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 --hmp        hapmap file (required)

 -h --help    display help menu

=cut
