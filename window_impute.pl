#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $hapmap_file;
my $het_pct = 0.04;
my $min_acm_pct = 0.3;
my $win_size = 31;
my $win_step = 1;
my $help;

my %chr_pos = ();
my %chr_pos_recs = ();
my %hets = ('A' => 'M','C' => 'M', 'T' => 'W', 'G' => 'W', 'N' => 'N');
my $offset = int($win_size / 2);
my $sample_count = 0;

parse_args();
parse_hmp();
win_proc();

exit(0);


sub win_proc {
	print(STDERR "warning: $0 will remove the first and last $offset variant records\n");

	foreach my $chr (sort keys %chr_pos) {
		for (my $win_start_index = 0; $win_start_index <= ($#{$chr_pos{$chr}} - $win_size); $win_start_index += $win_step) {
			my %imputed_calls = ();
			my %odd_haps = ();
			my %win_haps = ();
			my $win_stop_index = $win_start_index + $win_size - 1;

			foreach my $test_index ($win_start_index..$win_stop_index) {
				my $test_pos = $chr_pos{$chr}[$test_index];
				my ($rs, $alleles, $chrom, $pos, $strand, $assy, $center, $prot_lsid, $array_lsid, $panel_lsid, $qc_code, @samples) = split(/\t/, $chr_pos_recs{$chr}{$test_pos});

				foreach my $sample_index (0..$#samples) {
					push(@{$win_haps{$sample_index}}, $samples[$sample_index]);
				}
			}

			foreach my $sample_index (keys %win_haps) {
				my $imputed_call = 'NA';
				my $ac_count = 0;
				my %orig_calls = ();

				foreach my $sample_call (@{$win_haps{$sample_index}}) {
					$orig_calls{$sample_call}++;

					if ($sample_call eq 'A' || $sample_call eq 'C') {
						$ac_count++;
					}
				}

				my $a_count = 0;
				my $c_count = 0;

				if (exists($orig_calls{'A'})) {
					$a_count = $orig_calls{'A'};
				}

				if (exists($orig_calls{'C'})) {
					$c_count = $orig_calls{'C'};
				}

				if ($ac_count > $win_size * $min_acm_pct) {
					if ($a_count / $ac_count > (1 - $het_pct * 2)) {
						$imputed_call = 'A';
					}

					elsif ($c_count / $ac_count > (1 - $het_pct * 2)) {
						$imputed_call = 'C';
					}

					else {
						my $match;

						if (join('', @{$win_haps{$sample_index}}) =~ /([AN]+[CN]+|[CN]+[AN]+)/) {
							$match = $1;
						}

						if (defined($match) && length($match) == $win_size) {
							$imputed_call = 'N';

							if ($a_count > $c_count) {
								$imputed_call = 'A';
							}

							elsif ($c_count > $a_count) {
								$imputed_call = 'C';
							}
						}

						else {
							$odd_haps{$sample_index}++;
						}
					}
				}

				$imputed_calls{$sample_index} = $imputed_call;
			}

			foreach my $sample_index (keys %imputed_calls) {
				my $imputed_call = $imputed_calls{$sample_index};

				if ($imputed_call eq 'NA') {
					my $sample_offset_call = $win_haps{$sample_index}[$offset];

					if (exists($odd_haps{$sample_index}) && scalar(keys %odd_haps) < $sample_count * $het_pct) {
						# leaves as N if no explicit call, else homogenizes to het call
						$imputed_calls{$sample_index} = $hets{$sample_offset_call};
					}

					else {
						$imputed_calls{$sample_index} = $sample_offset_call;
					}
				}
			}

			# process offset record
			my $offset_pos = $chr_pos{$chr}[$win_start_index + $offset];
			my ($rs, $alleles, $chrom, $pos, $strand, $assy, $center, $prot_lsid, $array_lsid, $panel_lsid, $qc_code, @samples) = split(/\t/, $chr_pos_recs{$chr}{$offset_pos});

			if ($samples[1] eq 'G/T' || $samples[1] eq 'T/G') {
				my $offset_pos = $chr_pos{$chr}[$win_start_index + $offset];

				print(STDOUT "$chr_pos_recs{$chr}{$offset_pos}\n");
			}

			else {
				my @imputed_calls = ();

				foreach my $imputed_call_index (sort { $a <=> $b } keys %imputed_calls) {
					push(@imputed_calls, $imputed_calls{$imputed_call_index});
				}

				print(STDOUT join("\t", $rs, $alleles, $chrom, $pos, $strand, $assy, $center, $prot_lsid, $array_lsid, $panel_lsid, $qc_code, @imputed_calls), "\n");
			}
		}
	}

	return(0);
}


sub parse_hmp {
	open(HMP, '<', $hapmap_file) or error("can't read $hapmap_file: $!");

	while (my $line = <HMP>) {
		chomp($line);

		my ($rs, $alleles, $chrom, $pos, $strand, $assy, $center, $prot_lsid, $array_lsid, $panel_lsid, $qc_code, @samples) = split(/\t/, $line);

		if ($rs eq 'rs#') {
			$sample_count = scalar(@samples);

			print(STDOUT "$line\n");

			next();
		}

		push(@{$chr_pos{$chrom}}, $pos);
		$chr_pos_recs{$chrom}{$pos} = $line;
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

	GetOptions ('hmp=s' => \$hapmap_file,
				'het=f' => \$het_pct,
				'a|acm=f' => \$min_acm_pct,
				'w|win=i' => \$win_size,
				's|step=i' => \$win_step,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	if (! defined($hapmap_file)) {
		arg_error('hapmap file required');
	}

	if (defined($het_pct) && ($het_pct < 0 || $het_pct > 1)) {
		arg_error('expected heterozygosity must be between 0 and 1');
	}

	if (defined($min_acm_pct) && ($min_acm_pct < 0 || $min_acm_pct > 1)) {
		arg_error('minimum ACM % must be between 0 and 1');
	}

	if (defined($win_size) && $win_size < 0) {
		arg_error('window size must be > 0');
	}

	if (defined($win_step) && $win_step < 0) {
		arg_error('window step must be > 0');
	}

	return(0);
}


__END__

=head1 NAME

window_impute.pl

=head1 SYNOPSIS

window_impute.pl [options] --hmp variants.hmp.txt

=head1 DESCRIPTION

window_impute.pl imputes hapmap variants

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 --hmp        hapmap file (required)

 --het        expected heterozygosity
                default: 0.04

 -a --acm     min variant AC%
                default: 0.3

 -w --win     window size
                default: 31

 -s --step    window step
                default: 1 

 -h --help    display help menu

=cut
