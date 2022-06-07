#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $orig_hapmap_file;
my $imp_hapmap_file;
my $parent_a_geno;
my $parent_c_geno;
my $help;

my %orig_samples = ();
my $parent_a_index = 0;
my $parent_c_index = 1;
my $total = 0;
my $mismatch = 0;

parse_args();
parse_hapmap($orig_hapmap_file);
parse_hapmap($imp_hapmap_file);

my $freq = 0;

if ($total > 0) {
	$freq = sprintf("%.4f", $mismatch / $total);
}

print(STDOUT "total: $total\nmismatch: $mismatch\nfreq: $freq\n");

exit(0);


sub parse_hapmap {
	my $hapmap_file = shift();
	my $mode;

	if ($hapmap_file eq $orig_hapmap_file) {
		$mode = 'orig';
	}

	elsif ($hapmap_file eq $imp_hapmap_file) {
		$mode = 'imp';
	}

	else {
		error("hapmap parsing mode cannot be determined: $hapmap_file");
	}

	open(HAPMAP, '<', $hapmap_file) or error("can't read $hapmap_file: $!");

	while (my $line = <HAPMAP>) {
		chomp($line);
		
		my ($rs, $alleles, $chrom, $pos, $strand, $assy, $center, $prot_lsid, $array_lsid, $panel_lsid, $qc_code, @samples) = split(/\t/, $line);
		
		if ($line =~ /^rs#/) {
			next();
		}

		if ($mode eq 'orig') {
			push(@{$orig_samples{$chrom}{$pos}}, @samples);
		}

		elsif ($mode eq 'imp') {
			if (! exists($orig_samples{$chrom}{$pos})) {
				print(STDERR "warning: original hapmap record not found for corresponding imputed record at chr: $chrom, pos: $pos\n");

				next();
			}

			my @orig_samples = @{$orig_samples{$chrom}{$pos}};
			my @imp_samples = @samples;
			my $imp_parent_a_allele = $imp_samples[$parent_a_index];
			my $imp_parent_c_allele = $imp_samples[$parent_c_index];

			if ($imp_parent_a_allele eq 'A' && $imp_parent_c_allele eq 'C') {
				foreach my $sample_index (2..$#samples) {
					my $orig_sample_allele = $orig_samples[$sample_index];
					my $imp_sample_allele = $imp_samples[$sample_index];

					if ($orig_sample_allele =~ /[ACGT-]/ && $imp_sample_allele =~ /[AC]/) {
						$total++;

						my $orig_parent_allele;

						if ($imp_sample_allele eq 'A') {
							$orig_parent_allele = $orig_samples[$parent_a_index];
						}

						elsif ($imp_sample_allele eq 'C') {
							$orig_parent_allele = $orig_samples[$parent_c_index];
						}

						else {
							next();
						}

						if ($orig_sample_allele ne $orig_parent_allele) {
							print(STDERR "chr: $chrom\tpos: $pos\tsample_index: $sample_index\torig: $orig_sample_allele\timp: $imp_sample_allele\torig_parent: $orig_parent_allele\n");

							$mismatch++;
						}
					}
				}
			}
		}
	}

	close(HAPMAP);

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

	GetOptions ('o|orig=s' => \$orig_hapmap_file,
				'i|imp=s' => \$imp_hapmap_file,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	if (! defined($orig_hapmap_file)) {
		arg_error('original hapmap file required');
	}

	if (! defined($imp_hapmap_file)) {
		arg_error('imputed hapmap file required');
	}

	return(0);
}


__END__

=head1 NAME

imputed_hapmap_comp.pl

=head1 SYNOPSIS

imputed_hapmap_comp.pl -o orig.hmp.txt -i imputed.hmp.txt

=head1 DESCRIPTION

imputed_hapmap_comp.pl compares genotypes between original and imputed hapmap files

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -o --orig    original hapmap file (required)

 -i --imp     imputed hapmap file (required)

 -h --help    display help menu

=cut
