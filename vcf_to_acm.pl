#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $vcf_file;
my $parent_a_geno;
my $parent_c_geno;
my $help;

parse_args();
parse_vcf();

exit(0);


sub parse_vcf {
	my $parent_a_index;
	my $parent_c_index;
	my %processed = ();

	open(VCF, '<', $vcf_file) or error("can't read $vcf_file: $!");

	while (my $line = <VCF>) {
		chomp($line);
		
		my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @samples) = split(/\t/, $line);

		if ($line =~ /^#/) {
			if ($line =~ /^#CHROM/) {
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
					error("parent A genotype: $parent_a_geno not found in $vcf_file");
				}

				if (! defined($parent_c_index)) {
					error("parent C genotype: $parent_c_geno not found in $vcf_file");
				}
			}

			print(STDOUT "$line\n");

			next();
		}


		if ($alt =~ /,/) {
			next();
		}

		if (exists($processed{$chrom}{$pos})) {
			print(STDERR "duplicate chrom/pos found: $line\n");

			next();
		}


		my $parent_a_allele = $samples[$parent_a_index];
		my $parent_c_allele = $samples[$parent_c_index];

		if ($parent_a_allele eq $parent_c_allele) {
			next();
		}

		my $ref_acm_allele;
		my $alt_acm_allele;

		if ($parent_a_allele eq '0/0' && $parent_c_allele eq '1/1') {
			$ref_acm_allele = 'A';
			$alt_acm_allele = 'C';
		}

		elsif ($parent_a_allele eq '1/1' && $parent_c_allele eq '0/0') {
			$ref_acm_allele = 'C';
			$alt_acm_allele = 'A';
		}

		else {
			$ref_acm_allele = 'T';
			$alt_acm_allele = 'G';
		}

		print(STDOUT join("\t", $chrom, $pos, $id, $ref_acm_allele, $alt_acm_allele, $qual, $filter, $info, $format, @samples), "\n");

		$processed{$chrom}{$pos}++;
	}

	close(VCF);

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

	GetOptions ('v|vcf=s' => \$vcf_file,
				'a=s' => \$parent_a_geno,
				'c=s' => \$parent_c_geno,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	if (! defined($vcf_file)) {
		arg_error('vcf file required');
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

	return(0);
}


__END__

=head1 NAME

vcf_to_acm.pl

=head1 SYNOPSIS

vcf_to_acm.pl --vcf variants.vcf -a 'parentA' -c 'parentC'

=head1 DESCRIPTION

vcf_to_acm.pl converts vcf files to ACM format
also removes non-biallelic and monomorphic variants

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -v --vcf     vcf file (required)

 -a           parent 'A' genotype name (from vcf header)

 -c           parent 'C' genotype name (from vcf header)

 -h --help    display help menu

=cut
