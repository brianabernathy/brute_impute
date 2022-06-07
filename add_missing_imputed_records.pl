#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $orig_hapmap_file;
my $imp_hapmap_file;
my $max_dist = 100;
my $help;

my %imp_pos = ();
my %imp_recs = ();
my $parent_a_index = 11;
my $parent_c_index = 12;

parse_args();
parse_hapmap($imp_hapmap_file);
parse_hapmap($orig_hapmap_file);

exit(0);


sub parse_hapmap {
	my $hapmap_file = shift();
	my $mode;
	my $recs_added = 0;

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
		
		#my ($rs, $alleles, $chrom, $pos, $strand, $assy, $center, $prot_lsid, $array_lsid, $panel_lsid, $qc_code, @samples) = split(/\t/, $line);
		my @fields = split(/\t/, $line);
		my $chr = $fields[2];
		my $pos = $fields[3];
		
		if ($line =~ /^rs#/) {
			if ($mode eq 'orig') {
				print(STDOUT "$line\n");
			}

			next();
		}

		if ($mode eq 'imp') {
			push(@{$imp_pos{$chr}}, $pos);
			push(@{$imp_recs{$chr}{$pos}}, @fields);
		}

		elsif ($mode eq 'orig') {
			if (exists($imp_recs{$chr}{$pos})) {
				print(STDOUT join("\t", @{$imp_recs{$chr}{$pos}}), "\n");

				next();
			}

			my $nearest_imp_index = num_bin_search($pos, \@{$imp_pos{$chr}});
			my $nearest_imp_pos = $imp_pos{$chr}[$nearest_imp_index];
			my $nearest_imp_dist = abs($nearest_imp_pos - $pos);
			my $pass_dist_check = 0;

			if ($nearest_imp_dist <= $max_dist) {
				print(STDERR join("\t", @fields[0..4], "nearest pos: $nearest_imp_pos\tdist: $nearest_imp_dist"), "\n");

				$pass_dist_check = 1;
			}

			else {
				print(STDERR join("\t", @fields[0..4], "nearest pos: $nearest_imp_pos\tdist: $nearest_imp_dist"), "\n");
			}

			my $imp_parent_a_val = $imp_recs{$chr}{$nearest_imp_pos}[$parent_a_index];
			my $imp_parent_c_val = $imp_recs{$chr}{$nearest_imp_pos}[$parent_c_index];
			my $invert_parents = 0;
	
			if ($imp_parent_a_val eq 'C' && ($imp_parent_c_val eq 'A' || $imp_parent_c_val eq 'N')) {
				$invert_parents = 1;
			}

			elsif ($imp_parent_c_val eq 'A' && ($imp_parent_a_val eq 'C' || $imp_parent_a_val eq 'N')) {
				$invert_parents = 1;
			}

			foreach my $index ($parent_a_index..$#fields) {
				my $mod_val = 'N';

				if ($pass_dist_check == 1) {
					$mod_val = $imp_recs{$chr}{$nearest_imp_pos}[$index];
				}

				if ($invert_parents == 1) {
					if ($mod_val eq 'A') {
						$mod_val = 'C';
					}

					elsif ($mod_val eq 'C') {
						$mod_val = 'A';
					}
				}

				$fields[$index] = $mod_val;
			}

			print(STDOUT join("\t", @fields), "\n");

			$recs_added++;
		}
	}

	close(HAPMAP);

	if ($mode eq 'orig') {
		print(STDERR "missing imputed variants added: $recs_added\n");
	}

	return(0);
}


# return exact or nearest match array index
sub num_bin_search {
	my $target_val = shift();
	my $ints_ref = shift();

	if (! defined($target_val) || ! defined($ints_ref)) {
		return(undef);
	}

	my $min = 0;
	my $max = $#$ints_ref;

	if ($target_val > $$ints_ref[$max]) {
		return($max);
	}

	elsif ($target_val <= $$ints_ref[$min]) {
		return(0);
	}

	while ($min <= $max) {
		if ($min == ($max - 1)) {
			my $dist_min = abs($target_val - $$ints_ref[$min]);
			my $dist_max = abs($target_val - $$ints_ref[$max]);

			if ($dist_min <= $dist_max) {
				return($min);
			}

			else {
				return($max);
			}
		}

		my $middle = int(($max + $min) / 2);

		if ($target_val < $$ints_ref[$middle]) {
			$max = $middle;

			next();
		}

		elsif ($target_val > $$ints_ref[$middle]) {
			$min = $middle;

			next();
		}

		return($middle);
	}

	return(undef);
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
				'd|dist=i' => \$max_dist,
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

	if ($max_dist < 1) {
		arg_error('maximum distance must be > 1');
	}

	return(0);
}


__END__

=head1 NAME

add_missing_imputed_records.pl

=head1 SYNOPSIS

add_missing_imputed_records.pl -o orig.hmp.txt -i imputed.hmp.txt

=head1 DESCRIPTION

add_missing_imputed_records.pl adds the nearest imputed variant alleles (within specified maximum distance, otherwise 'N' is used) for original hapmap variantpositions that are missing from the imputed hapmap.

window_impute.pl removes the first and last few (based on window size) variants from each chomosome. TASSEL will occasionally omit variants for various reasons. add_missing_imputed_records.pl is useful if you want to retain all original hapmap variant positions in the imputed hapmap.

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -o --orig    original hapmap file (required)

 -i --imp     imputed hapmap file (required)

 -d --dist    maximum distance for nearest imputed variant
                default: 100

              if no imputed variant is found within the maximum distance,
              the original hapmap position is retained with alleles set to 'N'

 -h --help    display help menu

=cut
