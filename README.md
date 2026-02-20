# mstatspop V.1.1.2

## Variability Analyses of multiple populations: Calculation and estimation of statistics and neutrality tests.

#### Sebastian E. Ramos-Onsins, Luca Ferretti, Emanuele Raineri, Giacomo Marmorini, William Burgos, Joan Jene and Gonzalo Vera

mstatspop is designed for users interested in performing detailed analyses of nucleotide variability at the genomic scale in multiple populations. The program allows the analysis of genomic variability using various types of input files, including gVCF files, transposed FASTA (tFASTA), conventional FASTA files, and simulated data in ms format.

mstatspop supports the simultaneous analysis of multiple populations, the inclusion of one or more outgroups for the study of polarized variants (ancestral and derived), and allows users to define the order of the input sequences according to their preferences.

The program is available at:

	`https://github.com/CRAGENOMICA/mstatspop`

It also includes complementary tools, available at:

	`https://github.com/CRAGENOMICA/fastaconvtr`

	`https://github.com/sramosonsins/gVCF2tFasta`

### Annotation File Inclusion

mstatspop can integrate information contained in annotation files in GTF format, allowing for the analysis of variability in different functional categories, such as synonymous sites, non-synonymous sites, 0-fold, 4-fold coding, non-coding, and other features indicated in the feature column (e.g., introns, UTRs, etc.).

### Sliding Windows

The program performs analyses using a sliding window approach, with user-defined size and slide parameters. The maximum window size corresponds to the full length of the scaffold or chromosome. Windows can be defined in two ways: (i) by their physical position (based on genomic coordinates) or (ii) by their effective length (based on the actual number of positions analyzed). Additionally, a BED file with specific coordinates (e.g., genes or regions of interest) can be provided and used as analysis windows.

### Statistics and Neutrality and Differentiation Tests

mstatspop includes a comprehensive set of genetic variability statistics and neutrality tests, as well as population differentiation statistics and Fst significance estimates using permutations. A particularly relevant feature is the calculation of frequency-based statistics considering positions with missing data. The program also implements optimal neutrality tests, in which the user can specify an expected alternative hypothesis. Furthermore, the null hypothesis can be modified based on a user-defined demographic model.

### Output Formats

mstatspop offers several output formats, including (i) an extended format for viewing scaffolds or entire chromosomes, (ii) tabular formats optimized for further processing, where each line corresponds to an individual window, (iii) formats compatible with programs such as dadi and SweepFinder, and (iv) a multi-population frequency spectrum (SFS) in tabular format.

### Included Complementary Tools

The mstatspop package includes several programs and scripts for manipulating input and output files: (i) *VCF2Tfasta*: converts VCF files to Tfasta format for further analysis, (ii) *tfa_merge*: combines multiple Tfasta files (compressed and indexed) for joint analysis, (iii) *fastaconvtr*: transforms GTF files into index files with the information needed for analysis (e.g., silent positions), (iv) *collect_data_columns.pl*: allows you to extract specific columns from output files in windows.

The options available in mstatspop are described below:

##Flags:

      -f [input format file: ms, fasta OR tfa (gz file indexed)]
      -i [path and name of the input file]
      -o [output format file: 0 (extended),
                              1 (single line/window),
                              2 (single line SFS/window),
                              3 (dadi-like format),
                              4 (single line pairwise distribution)
                              5 (single line freq. variant per line/window)
                              6 (SNP genotype matrix)
                              7 (SweepFinder-like format -only first pop-)
                              8 (single line/window: Frequency of each haplotype in the populations)
                              9 (single line/window: Frequency of variants per line and population)
                             10 (full extended)]
      -N [#_pops] [#samples_pop1] ... [#samples_popN]
      -n [name of the file containing the name(s) of scaffold(s) and their length (separated by a tab), one per line (ex. fai file)]
      -T [path and name of the output file]. DEFAULT stdout.
      
  #####  OPTIONAL GENERAL PARAMETERS:
   
      -G [outgroup (0/1)] (last population). DEFAULT 0.
      -u [include unknown positions (0/1)].  DEFAULT 0.
      -R [performs analysis using only rSFS (0/1)].  DEFAULT 0.
      -A [Alternative Spectrum File (Only for Optimal Test): alternative_spectrum for each population (except outg)
          File format: (average absolute values) header plus fr(0,1) fr(0,2) ... fr(0,n-1) theta(0)/nt,
          fr(1,1) fr(1,2) ... fr(1,n-1) theta(1)/nt...; -u 1 not allowed yet]
      -S [Null Spectrum File (only if -A is defined): null_spectrum for each population (except outg).
          (average absolute values) header plus fr(0,1) fr(0,2) ... fr(0,n-1) theta(0)/nt,
          fr(1,1) fr(1,2) ... fr(1,n-1) theta(1)/nt...]. DEFAULT SNM.
    Optional Parameters for fasta and tfa input files:
      -O [#_nsam] [number order of first sample, number 0 is the first sample] [second sample] ...etc. up to nsamples.
         DEFAULT current order.
      -t [# permutations per window (H0: Fst=0). Only available with option -u 0]. DEFAULT 0.
      -s [seed]. DEFAULT 123456.
      
  #####  PARAMETERS FOR TFASTA INPUT (-f tfa): 'SLIDING WINDOW ANALYSIS OF EMPIRICAL DATA'
   
      -w [window size].
        OR
      -W [file with the coordinates of each window [scaffold init end] (instead options -w and -z).
        DEFAULT one whole window.
    Optional:
      -z [slide size (must be a positive value)]. DEFAULT window size.
      -Z [first window size displacement [for comparing overlapped windows])]. DEFAULT 0.
      -Y [define window lengths in 'physical' positions (1) or in 'effective' positions (0)]. DEFAULT 1.
      -E [input file with weights for positions:
         include three columns with a header,
         first the physical positions (1...end),
         second the weight for positions and
         third a boolean weight for the variant (eg. syn variant in nsyn counts is 0.000)].
         DEFAULT all 1.000
         
  #####  PARAMETERS FOR MS INPUT (-f ms):'SIMULATION ANALYSIS OF A SINGLE REGION'
    Optional:
       
      -r [# ms iterations]. DEFAULT 1.
      -m [include mask_filename] DEFAULT -1 (all positions included).
         [mask_file format: 1st row with 'length' weights, next sample rows x lengths: missing 0, sequenced 1)].
         DEFAULT no mask.
      -v [ratio transitions/transversions]. DEFAULT 0.5.
      -F [force analysis to include outgroup (0/1) (0 in ms means ancestral)]. DEFAULT 0.
      -q [frequency of reverted mutation] (only with -F 1). DEFAULT 0.
      
   ##### PARAMETERS FOR FASTA INPUT (-f fasta): 'WHOLE REGION ANALYSIS'
   
    Optional:
      -p [Number of lineages per sequence (1/2)]. DEFAULT 1.
      -g [GFF_file]
         [add also: coding,noncoding,synonymous,nonsynonymous,silent, others (or whatever annotated)]
         [if 'synonymous', 'nonsynonymous', 'silent' add: Genetic_Code: Nuclear_Universal,mtDNA_Drosophila,mtDNA_Mammals,Other]
         [if 'Other', introduce the code for the 64 triplets in the order UUU UUC UUA UUG ... etc.].
         DEFAULT no annotation.
      -c [in case use coding regions, criteria to consider transcripts (max/min/first/long)]. DEFAULT long.
      -K [make a MASK file with the valid positions for this fasta. Useful for running ms simulations (1/0)]. DEFAULT 0.
   
   ##### HELP:
   
      -h [help and exit]

### The statistics and tests described in mstatspop are the following:

#### Effective length for each population:

      - Eff_length1_pop_outg[0]: Effective length1 for population 0. That is, considering the outgroup, how many positions have at least one sequence per population 0 and the outgroup exist. Useful in case considering missing positions (-u 1). In case -u 0 all lengths are equal. Useful for calculating divergence vs outgroup or versus populations.
      - Eff_length2_pop_outg[0]: Effective length2 for population 0. That is, considering the outgroup, how many positions have at least two sequences per population 0 and the outgroup exist. Useful in case considering missing positions (-u 1). In case -u 0 all lengths are equal. Useful for calculating levels of variability.
      - Eff_length3_pop_outg[0]: Effective length3 for population 0. That is, considering the outgroup, how many positions have at least three sequences per population 0 and the outgroup exist.Useful in case considering missing positions (-u 1). In case -u 0 all lengths are equal. Useful for calculating tests of neutrality (variances).
      - Eff_length1_pop[0]: Effective length1 for population 0. That is, how many positions have at least one sequence per population 0. Useful in case considering missing positions (-u 1). In case -u 0 all lengths are equal. Useful for calculating divergence versus populations.
      - Eff_length2_pop[0]: Effective length2 for population 0. That is, how many positions have at least two sequences per population 0. Useful in case considering missing positions (-u 1). In case -u 0 all lengths are equal. Useful for calculating levels of variability.
      - Eff_length3_pop[0]: Effective length3 for population 0. That is, how many positions have at least three sequences per population 0. Useful in case considering missing positions (-u 1). In case -u 0 all lengths are equal. Useful for calculating tests of neutrality (variances).

In principle, if the analysis contains an outgroup, all analysis are conditioned to the presence of the outgroup variant.

#### Estimates of variability for each population and divergence:

      - S[0]: Number of segregating sites at population 0
      - Theta(Wat)[0]:	Watterson (1975) estimation of variability.
      - Theta(Taj)[0]:	Tajima (1983) estimation of variability.
      - Theta(Fu&Li)[0]:	Fu and Li’s (1991) estimation of variability.
      - Theta(Fay&Wu)[0]: Fay and Wu’s (2000) estimation of variability.
      - Theta(Zeng)[0]:	Zeng’s (2006) estimation of variability.
      - Theta(Achaz,Wat)[0]: Watterson estimator without considering singletons (Achaz 2008)
      - Theta(Achaz,Taj)[0]: Tajima estimator without considering singletons (Achaz 2008)

      - Divergence[0]: Divergence between population 0 and the outgroup.

      - an_x[0]: ($\sum_{i=1}^{n-1}1/i$, that is,  sumatory of (1/i) from i=1 to n-1, being n the number of samples. In case missing data (-u 1), sum of an for each position and divided by the eff_length2_pop. Used for Watterson estimator and neutrality tests.
      - bn_x[0]: ($\sum_{i=1}^{n-1}1/i^2$, that is,  sumatory of (1/i^2) from i=1 to n-1, being n the number of samples. In case missing data (-u 1), sum of an for each position and divided by the eff_length2_pop. Used for neutrality tests.
      - an_xo[0]: Same but considering only those positions were the outgroup is present.
      - bn_xo[0]: Same but considering only those positions were the outgroup is present.

      - Haplotype diversity and number of haplotypes for each population:
      - HapW[0]: Haplotype diversity at population 0.
      - nHap[0]: number of haplotypes at population 0.

#### Neutrality tests for each population:

      - Tajima D[0]: Test of Tajima (1989)
      - Fu&Li D[0]: 	Test of Fu and Li D (1993)
      - Fu&Li F[0]: 	Test of Fu and Li D (1993)
      - Fay&Wu norm H[0]: Test of Fay and Wu (2000) normalized (Achaz 2009, Zeng 2006)
      - Fay&WuH[0]: Test of Fay and Wu (2000)
      - Ferretti L[0]: Test of Ferretti similar to Fay and Wu but using theta  Watt instead Theta Taj (2017)
      - Zeng E[0]: 	Test of Zeng (2006)
      - Achaz Y[0]: 	Test of Achaz (2008)
      - Fs[0]: Test Fs (Fu 1997)
      - R2[0]: Test of R2 (2002)

#### Variants assigned to exclusive, fixed, polymorphic but fixed in rest of pops, and shared.

      - Ss[rest] are shared variants between populations but fixed within:
      - Sx[0]:  Exclusive variants at population 0
      - Sf[0]: Fixed variants at population 0
      - Ss: Shared variants among populations
      - Sxf[0]: Exclusive variants at population 0 that are fixed at other populations.

#### Mismatch distribution statistics:

      - SDev[0]: 	Standard deviation of Tajima’s estimator.
      - Skewness[0]: Third standardised moment of Tajima’s estimator.
      - Kurtosis[0]: 	Fourth standardised moment of Tajima’s estimator.

#### Differentiation with nucleotide sequences (even with missing data):

      - Fstall: 1 - mean(pi_within/nt)/mean(pi_among/nt)
      - Fst1all[pop1]: 1- mean(pi_within/nt[rest],pi_within_rest/nt[pop1])/mean(pi_among/nt)
      - Fst[pop1][pop2]: 1 - mean(pi_within/nt[pop1],pi_within/nt[pop2])/ pi_among/nt[pop1][pop2] (Hudson et al. 1992)

#### Differentiation with haplotype data (phasing is necessary):

      - Fsthall: 1 - mean(pih_within)/mean(pih_among)
      - Fsth1all[pop1]: 1- mean(pih_within[pop1]-pih_within[rest])/mean(pih_among)
      - Fsth[pop1][pop2]: 1 - mean(pih_within[pop1],pih_within[pop2])/pih_among[pop1][pop2] 

#### Frequency of variants for each population:

      - fr[0,1]:	Site frequency Spectrum per population (here population 0)

#### Frequency of each haplotype in the populations:

      - frH[0,hap00]:	frequency of each haplotype per population

#### Joint frequency distribution for each variant and population (No included variants that are missing or polymorphic in the outgroup):

      - SNP[xx]	frequency of each SNP at each population (in tabs).


### Create index file and tfasta convertion

To create an index file for tfasta file use `tfa_index` program. 
Note : if the tfasta file is in version 1, it will convert it to version 2 and create the index.
`tfa_index` program also accepts weight files. In this case, it will create an index for the weight file and convert it to the correct format if needed.

```bash
Usage:  tfa_index  [options] <input.tfa|input.tfa.bgz|input.tfa.gz|weights.txt|weights.txt.gz>
Options:
  --version
  --help
  --threads <int>
  --force
  --weight
   the input file is a weight file. In this case create and index for it. and convert it to the correct format
  --output FILE           set custom name to the output file, only used when converting from TFAv1 to TFAv2 or compressing the input file
  ```

Convertion example 
```bash
  tfa_index  ./Examples/V0.1.0/100Kchr10.tfa.gz -o ./Examples/V1.0.0//100Kchr10.tfa.gz
```

Create an index example
```bash
  tfa_index  ./Examples/V1.0.0/100Kchr10.tfa.gz
  # use -f to force create the index file if it already exists
  tfa_index  ./Examples/V1.0.0/100Kchr10.tfa.gz -f
```


### Merge multiple tfa files

use `tfa_merge` to merge two tfasta files into one.
### tfa_merge command line usage
```bash
tfa_merge -i file1.tfa.gz -i file2.tfa.gz -o merged.tfa.gz
###
### Options:
- `-i, --input FILE`: Input file (specify twice for two files)
- `-o, --output FILE`: Output file name
- `-f, --force`: Force overwrite of output file
- `-h, --help`: Show help message
```
