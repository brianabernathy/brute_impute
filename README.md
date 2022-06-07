# brute_impute
basic imputation for skim seq genotying

---

### hyper_variable_region_filter.pl

Description:

    hyper_variable_region_filter.pl remove variants from hyper-variable regions

Usage:

    hyper_variable_region_filter.pl [options] --hmp variants.hmp.txt

Options:

     --hmp        hapmap file (required)

     -d --div     minimum tolerable window divergence %, filter variants above threshold
                    default: 20

     -w --win     window size (bps)
                    default: 100

     -h --help    display help menu

---

### vcf_to_acm.pl

Description:

    vcf_to_acm.pl converts vcf files to ACM format, also removes non-biallelic and monomorphic variants

Usage:

    vcf_to_acm.pl --vcf variants.vcf -a 'parentA' -c 'parentC'

Options:

     -v --vcf     vcf file (required)

     -a           parent 'A' genotype name from vcf header (required)

     -c           parent 'C' genotype name from vcf header (required)

     -h --help    display help menu

---

### tassel_ld_filter.pl

Description:

    tassel_ld_filter.pl filters hapmap files using TASSEL's windowed linkage disequilibrium (LD) calculations

Usage:

    tassel_ld_filter.pl --hapmap variants.hap.txt -a 'parentA' -c 'parentC'

Options:

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

---

### window_impute.pl

Description

    window_impute.pl imputes hapmap variants

Usage:

    window_impute.pl [options] --hmp variants.hmp.txt

Options:

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

---

### imputed_hapmap_comp.pl

Description:
    imputed_hapmap_comp.pl compares genotypes between original and imputed hapmap files

Usage:

    imputed_hapmap_comp.pl -o orig.hmp.txt -i imputed.hmp.txt

Options:

     -o --orig    original hapmap file (required)

     -i --imp     imputed hapmap file (required)

     -h --help    display help menu

---

### hapmap_all_n_vars_count.pl

Description

    hapmap_all_n_vars_count.pl counts the number of variant records with all 'N' alleles. TASSEL may silently write all 'N' records when it encounters issues.

Usage:

    hapmap_all_n_vars_count.pl --hmp variants.hmp.txt

Options:

     --hmp        hapmap file (required)

     -h --help    display help menu

---

### add_missing_imputed_records.pl

Description:

    add_missing_imputed_records.pl adds original hapmap variants that are missing from the imputed hapmap.

    window_impute.pl removes the first and last few (based on window size) variants from each chomosome. TASSEL will occasionally omit variants for various reasons. add_missing_imputed_records.pl is useful if you want to retain all original hapmap variants in the imputed hapmap.

Usage:

    add_missing_imputed_records.pl -o orig.hmp.txt -i imputed.hmp.txt

Options:

     -o --orig    original hapmap file (required)

     -i --imp     imputed hapmap file (required)

     -d --dist    maximum distance for nearest imputed variant
                    default: 100

                  if no imputed variant is found within the maximum distance,
                  the original hapmap position is retained with alleles set to 'N'

     -h --help    display help menu

---

### downsample_hapmap.pl

Description: 

    downsample_hapmap.pl generates a downsampled hapmap using the specified minimum distance between variants

Usage:

    downsample_hapmap.pl [options] --hmp variants.hmp.txt

Options:

     --hmp        hapmap file (required)

     -d --dist    minimum distance between variants
                    default: 200

     -h --help    display help menu

