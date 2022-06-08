#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $vcf_file;
my $prefix;
my $postfix;
my $help;

my @headers = ();
my %chr_recs = ();

parse_args();
parse_vcf();
print_vcfs();

exit(0);


sub parse_vcf {
	open(VCF, '<', $vcf_file) or error("can't read $vcf_file: $!");

	while (my $line = <VCF>) {
		chomp($line);
		
		if ($line =~ /^#/) {
			push(@headers, $line);

			next();
		}

		my @fields = split(/\t/, $line);
		my $chrom = $fields[0];

		push @{$chr_recs{$chrom}}, join("\t", @fields);
	}

	close(VCF);

	return(0);
}


sub print_vcfs {
	foreach my $chr (sort keys %chr_recs) {
		my $chr_vcf_file;

		if (defined($prefix)) {
			$chr_vcf_file .= "$prefix";
		}

		$chr_vcf_file .= "$chr";

		if (defined($postfix)) {
			$chr_vcf_file .= "$postfix";
		}

		$chr_vcf_file .= '.vcf';

		open(CHR_VCF, '>', $chr_vcf_file) or error("can't write $chr_vcf_file: $!");

		foreach my $header (@headers) {
			print(CHR_VCF "$header\n");
		}

		foreach my $record (@{$chr_recs{$chr}}) {
			print(CHR_VCF "$record\n");
		}

		close(CHR_VCF);
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

	GetOptions ('v|vcf=s' => \$vcf_file,
				'pre=s' => \$prefix,
				'post=s' => \$postfix,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	if (! defined($vcf_file)) {
		arg_error('vcf file required');
	}

	return(0);
}


__END__

=head1 NAME

split_vcf.pl

=head1 SYNOPSIS

split_vcf.pl [options] --vcf variants.vcf

=head1 DESCRIPTION

split_vcf.pl writes individual vcf files for each chromosome. Output vcf files
will be written to 'chromosome.vcf' by default. Text can be pre/postfixed to
chromosome names using the --pre/--post options.

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -v --vcf     vcf file (required)

 --pre        prefix text added before chromosome in output filename

 --post       postfix text added after chromosome in output filename

 -h --help    display help menu

=cut
