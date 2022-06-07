#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $hapmap_file;
my $parent_a_geno;
my $parent_c_geno;
my $tassel_path;
my $ld_win_size = 20;
my $min_ld = 0.93;
my $help;

my %low_lds = ();

parse_args();
proc_tassel_lds();
parse_hapmap();

exit(0);


sub proc_tassel_lds {
	my $rand_int = int(rand(999999));
	my $ld_tmp_file = "$rand_int.ld.tmp.txt";
	my $tassel_cmd = "$tassel_path/run_pipeline.pl -Xmx32g -fork1 -h $hapmap_file -ld -ldWinSize $ld_win_size -td_tab $ld_tmp_file";

	print(STDERR "tassel command: $tassel_cmd\n");

	`$tassel_cmd`;

	my %lds = ();

	open(LD, '<', $ld_tmp_file) or error("can't read $ld_tmp_file: $!");

	while (my $line = <LD>) {
		if ($line =~ /^Locus1/) {
			next();
		}

		chomp($line);

		my ($locus1, $pos1, $site1, $num_states1, $states1, $freq1, $locus2, $pos2, $site2, $num_states2, $states2, $freq2, $dist_bp, $r_squared, $d_prime, $p_diseq, $n) = split(/\t/, $line);

		if ($p_diseq ne 'NaN') {
			push(@{$lds{$locus1}{$pos1}}, $d_prime);
		}
	}

	close(LD);

	unlink($ld_tmp_file);

	foreach my $locus (keys %lds) {
		foreach my $pos (keys %{$lds{$locus}}) {
			my $sum_ld = 0;

			foreach my $d_prime (@{$lds{$locus}{$pos}}) {
				$sum_ld += $d_prime;
			}

			my $mean_ld = $sum_ld / scalar(@{$lds{$locus}{$pos}});

			if ($mean_ld < $min_ld) {
				$low_lds{$locus}{$pos}++;
			}
		}
	}

	return(0);
}


sub parse_hapmap {
	my $parent_a_index;
	my $parent_c_index;
	my %rec_counts = ();

	open(HAPMAP, '<', $hapmap_file) or error("can't read $hapmap_file: $!");

	while (my $line = <HAPMAP>) {
		chomp($line);
		
		my ($rs, $alleles, $chrom, $pos, $strand, $assy, $center, $prot_lsid, $array_lsid, $panel_lsid, $qc_code, @samples) = split(/\t/, $line);
		
		if ($line =~ /^rs#/) {
			foreach my $index (0..$#samples) {
				my $sample = $samples[$index];

				if ($sample eq $parent_a_geno) {
					$parent_a_index = $index;
				}

				elsif ($sample eq $parent_c_geno) {
					$parent_c_index = $index;
				}
			}

			if (! defined($parent_a_index)) {
				error("parent A genotype: $parent_a_geno not found in $hapmap_file");
			}

			if (! defined($parent_c_index)) {
				error("parent C genotype: $parent_c_geno not found in $hapmap_file");
			}

			print(STDOUT "$line\n");

			next();
		}


		my $parent_a_allele = $samples[$parent_a_index];
		my $parent_c_allele = $samples[$parent_c_index];
		my $type = 'seg';

		if ($parent_a_allele =~ /[ACGT]/ && $parent_a_allele eq $parent_c_allele) {
			$type = 'pnp';
		}

		elsif ($parent_a_allele !~ /[ACGT]/ || $parent_c_allele !~ /[ACGT]/) {
			$type = 'amb';
		}


		if (exists($low_lds{$chrom}{$pos})) {
			$rec_counts{'lowLD'}{$type}++;
		}

		else {
			$rec_counts{'inLD'}{$type}++;

			print(STDOUT "$line\n");
		}
	}

	close(HAPMAP);

	foreach my $ld (sort keys %rec_counts) {
		foreach my $type (sort keys %{$rec_counts{$ld}}) {
			print(STDERR "$ld $type count: $rec_counts{$ld}{$type}\n");
		}
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

	GetOptions ('hapmap=s' => \$hapmap_file,
				'a=s' => \$parent_a_geno,
				'c=s' => \$parent_c_geno,
				't|tassel=s' => \$tassel_path,
				'w|win=i' => \$ld_win_size,
				'm|min=f' => \$min_ld,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	if (! defined($hapmap_file)) {
		arg_error('hapmap file required');
	}

	if (! defined($parent_a_geno)) {
		arg_error('parent A genotype required');
	}

	if (! defined($parent_c_geno)) {
		arg_error('parent C genotype required');
	}

	if ($parent_a_geno eq $parent_c_geno) {
		arg_error('parent A and C genotypes must be different');
	}

	if (defined($tassel_path)) {
		if (! -e $tassel_path) {
			error("TASSEL path: $tassel_path does not exist");
		}
	}

	else {
		$tassel_path = qx(which run_pipeline.pl);

		chomp($tassel_path);

		if (! defined($tassel_path) || $tassel_path eq '') {
			error('TASSEL not found in $PATH, specify using -t/--tassel option');
		}

		$tassel_path =~ s/\/*run_pipeline\.pl//;

		print(STDERR "using TASSEL found at $tassel_path\n");
	}

	$tassel_path =~ s/\/$//;

	if ($ld_win_size < 1) {
		arg_error('TASSEL LD window size must be > 0');
	}

	if ($min_ld < 0 || $min_ld > 1) {
		arg_error('minimum window LD must be between 0 and 1');
	}

	return(0);
}


__END__

=head1 NAME

tassel_ld_filter.pl

=head1 SYNOPSIS

tassel_ld_filter.pl --hapmap variants.hap.txt -a 'parentA' -c 'parentC'

=head1 DESCRIPTION

tassel_ld_filter.pl filters hapmap files using TASSEL's windowed linkage disequilibrium (LD) calculations

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 --hapmap     hapmap file (required)

 -a           parent 'A' genotype name from hapmap header (required)

 -c           parent 'C' genotype name from hapmap header required)

 -t --tassel  path to TASSEL script directory
                default: attempt to locate in $PATH

 -w --win     TASSEL LD window size
                default: 20

 -m --min     minimum window LD
                default: 0.93

 -h --help    display help menu

=cut
