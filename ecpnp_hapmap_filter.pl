#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $hmp_file;
my $help;

my $parent_a_index = 0;
my $parent_c_index = 1;

parse_args();
parse_hmp();

exit(0);


sub parse_hmp {
	my %filtered = ();

	open(HMP, '<', $hmp_file) or error("can't read $hmp_file: $!");

	while (my $line = <HMP>) {
		chomp($line);

		my ($rs, $alleles, $chrom, $pos, $strand, $assy, $center, $prot_lsid, $array_lsid, $panel_lsid, $qc_code, @samples) = split(/\t/, $line);

		if ($line =~ /^rs#\t/) {
			print(STDOUT "$line\n");

			next();
		}

		my $parent_a_val = $samples[$parent_a_index];
		my $parent_c_val = $samples[$parent_c_index];
		my $parents_explicitly_called = 1;

		if ($parent_a_val !~ /[ACGT-]/) {
			$parents_explicitly_called = 0;
			$filtered{'parent A not explicitly called'}++;
		}

		if ($parent_c_val !~ /[ACGT-]/) {
			$parents_explicitly_called = 0;
			$filtered{'parent C not explicitly called'}++;
		}

		if ($parent_a_val eq $parent_c_val) {
			if ($parents_explicitly_called == 1) {
				$filtered{'parents explicitly called and non-poly'}++;
			}

			else {
				$filtered{'parents not explicitly called and non-poly'}++;
			}
		}

		else {
			if ($parents_explicitly_called == 1) {
				print(STDOUT "$line\n");
			}

			else { 
				$filtered{'parents not explicitely called and poly'}++;
			}
		}
	}

	close(HMP);

	foreach my $filter_type (sort keys %filtered) {
		print(STDERR "#$filter_type: $filtered{$filter_type}\n");
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

ecpnp_hapmap_filter.pl

=head1 SYNOPSIS

ecpnp_hapmap_filter.pl --hmp variants.hmp.txt

=head1 DESCRIPTION

ecpnp_hapmap_filter.pl removes hapmap records where parents are not explicitely called and polymorphic

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 --hmp        hapmap file (required)

 -h --help    display help menu

=cut
