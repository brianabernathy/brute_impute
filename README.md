# brute_impute

basic imputation for skim seq genotying

---

## pipeline usage overview

Test data sets, scripts and additional documentation specific to skim-seq processing are also available at the [PanPipes](https://github.com/USDA-ARS-GBRU/PanPipes#demos "PanPipes") repository.

The following assumes a biparental mapping population VCF file is available. VCF files may be split into individual chromosomes and processed separately if desired or to reduce the computational resources required for processing, but is not required. See `split_vcf.pl` and `merge_hapmaps.pl` for help splitting VCFs and merging imputed hapmaps.

TASSEL [(https://tassel.bitbucket.io/)](https://tassel.bitbucket.io/ "TASSEL") is required for sorting, format conversion and MAF filtering operations. A population 'taxa' file is required for TASSEL MAF filtering, see TASSEL documentation for details.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sort vcf
- `${TASSEL_DIR}/run_pipeline.pl -SortGenotypeFilePlugin -inputFile pop.vcf -outputFile pop.sort.vcf -fileType VCF`

  convert vcf to hapmap
- `${TASSEL_DIR}/run_pipeline.pl -vcf pop.sort.vcf -export pop.sort.hmp.txt -exportType Hapmap`

  filter hyper-variable regions
- `${BRUTE_IMPUTE_DIR}/hyper_variable_region_filter.pl --hmp pop.sort.hmp.txt > pop.sort.hypervar.hmp.txt`

  minor allele frequency (MAF) filter
  sample parameters shown, adjust as appropriate
- `${TASSEL_DIR}/run_pipeline.pl -fork1 -h pop.sort.hmp.txt -includeTaxaInFile pop.taxas.txt -filterAlign -filterAlignMinCount 10 -filterAlignMinFreq 0.2 -filterAlignMaxFreq 0.8 -export pop.sort.hypervar.maf.hmp.txt,pop.sort.hypervar.maf.json.gz -exportType Hapmap -runfork1`

  convert hapmap to vcf
- `${TASSEL_DIR}/run_pipeline.pl -h pop.sort.hypervar.maf.hmp.txt -export pop.sort.hypervar.maf.vcf -exportType VCF`

  convert vcf to acm format
- `${BRUTE_IMPUTE_DIR}/vcf_to_acm.pl -v pop.sort.hypervar.maf.vcf -a 'parentA' -c 'parentC' > pop.sort.hypervar.maf.acm.vcf`

  convert acm-formatted vcf to hapmap
- `${TASSEL_DIR}/run_pipeline.pl -vcf pop.sort.hypervar.maf.acm.vcf -export pop.sort.hypervar.maf.acm.hmp.txt -exportType Hapmap`

  filter variants in low linkage disequilibrium (LD)
- `${BRUTE_IMPUTE_DIR}/tassel_ld_filter.pl --hapmap pop.sort.hypervar.maf.acm.hmp.txt -a 'parentA' -c 'parentC' --tassel $TASSEL_DIR > pop.sort.hypervar.maf.acm.ld.hmp.txt`

  impute (removes first and last few variants from each chromosome)
- `${BRUTE_IMPUTE_DIR}/window_impute.pl --hmp pop.sort.hypervar.maf.acm.ld.hmp.txt > pop.window_impute.hmp.txt`

  add missing variants removed by imputation
- `${BRUTE_IMPUTE_DIR}/add_missing_imputed_records.pl -o pop.sort.hypervar.maf.acm.ld.hmp.txt -i pop.window_impute.hmp.txt > pop.window_impute.mod.hmp.txt`

  compare pre and post-imputation hapmaps
- `${BRUTE_IMPUTE_DIR}/imputed_hapmap_comp.pl -o pop.tassel.sort.hypervar.maf.acm.ld.hmp.txt -i pop.window_impute.mod.hmp.txt 1> pop.window_impute.mod.mismatch.stats 2> pop.window_impute.mod.mismatch.genos`

  count variants with Ns for all population genotypes
- `${BRUTE_IMPUTE_DIR}/hapmap_all_n_vars_count.pl --hmp pop.window_impute.mod.hmp.txt > pop.window_impute.mod.hmp.all.n.var.count`


More complete documentation for each tool is available below. Tools can also be ran without options or -h/--help to display the help menu.


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

    vcf_to_acm.pl converts vcf files to ACM format, also removes non-biallelic and
    monomorphic variants

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

    tassel_ld_filter.pl filters hapmap files using TASSEL's windowed linkage
    disequilibrium (LD) calculations

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

    hapmap_all_n_vars_count.pl counts the number of variant records with all 'N' alleles.
    TASSEL may silently write all 'N' records when it encounters issues.

Usage:

    hapmap_all_n_vars_count.pl --hmp variants.hmp.txt

Options:

     --hmp        hapmap file (required)

     -h --help    display help menu

---

### add_missing_imputed_records.pl

Description:

    add_missing_imputed_records.pl adds original hapmap variants that are missing from
    the imputed hapmap.

    window_impute.pl removes the first and last few (based on window size) variants from 
    each chomosome. TASSEL will occasionally omit variants for various reasons.
    add_missing_imputed_records.pl is useful if you want to retain all original hapmap
    variants in the imputed hapmap.

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

    downsample_hapmap.pl generates a downsampled hapmap using the specified minimum
    distance between variants

Usage:

    downsample_hapmap.pl [options] --hmp variants.hmp.txt

Options:

     --hmp        hapmap file (required)

     -d --dist    minimum distance between variants
                    default: 300

     -h --help    display help menu

---

### split_vcf.pl

Description:

    split_vcf.pl writes individual vcf files for each chromosome. Output vcf files
    will be written to 'chromosome.vcf' by default. Text can be pre/postfixed to
    chromosome names using the --pre/--post options.

Usage:

    split_vcf.pl [options] --vcf variants.vcf

Options:

     -v --vcf     vcf file (required)

     --pre        prefix text added before chromosome in output filename

     --post       postfix text added after chromosome in output filename

     -h --help    display help menu

---

### merge_hapmaps.pl

Description:

merge_hapmaps.pl is a simple tool to combine individual hapmaps into a merged file

Usage:

    merge_hapmaps.pl -p "pop.chr*.hmp.txt"

Options:

     -p --pattern  hapmap file pattern (required)

     -h --help     display help menu

