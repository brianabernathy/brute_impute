#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $hmp_file;
my $min_dist = 300;
my $help;

parse_args();
parse_hmp();

exit(0);


sub parse_hmp {
	my $prev_pos;

	open(HMP, '<', $hmp_file) or error("can't read $hmp_file: $!");

	while (my $line = <HMP>) {
		chomp($line);

		my ($rs, $alleles, $chrom, $pos, $strand, $assy, $center, $prot_lsid, $array_lsid, $panel_lsid, $qc_code, @samples) = split(/\t/, $line);

		if ($rs eq 'rs#') {
			print(STDOUT "$line\n");

			next();
		}

		my $print_rec = 0;
		my $dist;

		if (! defined($prev_pos)) {
			$print_rec = 1;
		}

		else {
			$dist = $pos - $prev_pos;

			if ($dist >= $min_dist) {
				$print_rec = 1;
			}
		}


		if ($print_rec == 1) {
			print(STDOUT "$line\n");

			$prev_pos = $pos;
		}

		else {
			print(STDERR "skipping pos: $pos\tdist: $dist\n");
		}
	}

	close(HMP);

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
				'd|dist=i' => \$min_dist,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	if (! defined($hmp_file)) {
		arg_error('hapmap file required');
	}

	if ($min_dist < 1) {
		arg_error('minimum distance must be > 0');
	}

	return(0);
}


__END__

=head1 NAME

downsample_hapmap.pl

=head1 SYNOPSIS

downsample_hapmap.pl [options] --hmp variants.hmp.txt

=head1 DESCRIPTION

downsample_hapmap.pl generates a downsampled hapmap using the specified minimum distance between variants

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 --hmp        hapmap file (required)

 -d --dist    minimum distance between variants
                default: 200

 -h --help    display help menu

=cut
