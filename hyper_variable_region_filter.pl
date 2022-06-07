#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $hmp_file;
my $min_div_pct = 20;
my $win_size = 100;
my $help;

my %chr_pos = ();
my %chr_pos_filter = ();

parse_args();
parse_hmp('get_pos');
win_proc();
parse_hmp('print');

exit(0);


sub win_proc {
	foreach my $chr (sort keys %chr_pos) {
		foreach my $win_start_index (0..$#{$chr_pos{$chr}} - 1) {
			my $win_start_pos = $chr_pos{$chr}[$win_start_index];
			my $win_stop_pos = $win_start_pos + $win_size - 1;
			my @win_snps = ();

			foreach my $test_index ($win_start_index..$#{$chr_pos{$chr}}) {
				my $test_pos = $chr_pos{$chr}[$test_index];

				if ($test_pos <= $win_stop_pos) {
					push(@win_snps, $test_pos);

					next();
				}

				my $div_pct = 100 * scalar(@win_snps) / $win_size;

				if ($div_pct > $min_div_pct) {
					foreach my $snp_pos (@win_snps) {
						$chr_pos_filter{$chr}{$snp_pos}++;
					}
				}

				undef(@win_snps);
				
				last();
			}
		}
	}

	foreach my $chr (sort keys %chr_pos_filter) {
		my $filtered_count = keys %{$chr_pos_filter{$chr}};

		print(STDERR "$filtered_count variants filtered from $chr\n");
	}

	return(0);
}


sub parse_hmp {
	my $mode = shift();

	open(HMP, '<', $hmp_file) or error("can't read $hmp_file: $!");

	while (my $line = <HMP>) {
		chomp($line);

		my ($rs, $alleles, $chrom, $pos, $strand, $assy, $center, $prot_lsid, $array_lsid, $panel_lsid, $qc_code, @samples) = split(/\t/, $line);

		if ($mode eq 'get_pos') {
			if ($rs eq 'rs#') {
				next();
			}

			push(@{$chr_pos{$chrom}}, $pos);
		}

		elsif ($mode eq 'print') {
			if (exists($chr_pos_filter{$chrom}{$pos})) {
				next();
			}

			print(STDOUT "$line\n");
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
				'd|div=i' => \$min_div_pct,
				'w|win=i' => \$win_size,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	if (! defined($hmp_file)) {
		arg_error('hapmap file required');
	}

	if (defined($min_div_pct) && ($min_div_pct < 0 || $min_div_pct > 100)) {
		arg_error('minimum divergence % must be between 0 and 100');
	}

	if (defined($win_size) && $win_size < 0) {
		arg_error('window size must be > 0');
	}

	return(0);
}


__END__

=head1 NAME

hyper_variable_region_filter.pl

=head1 SYNOPSIS

hyper_variable_region_filter.pl [options] --hmp variants.hmp.txt

=head1 DESCRIPTION

hyper_variable_region_filter.pl remove variants from hyper-variable regions

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 --hmp        hapmap file (required)

 -d --div     minimum tolerable window divergence %, filter variants above threshold
                default: 20

 -w --win     window size (bps)
                default: 100

 -h --help    display help menu

=cut
