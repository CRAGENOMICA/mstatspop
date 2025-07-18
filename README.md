# mstatspop v.1.0.0

## Variability Analyses of multiple populations: Calculation and estimation of statistics and neutrality tests. 

#### Sebastian E. Ramos-Onsins, Luca Ferretti, Emanuele Raineri, Giacomo Marmorini, William Burgos, Joan Jene and Gonzalo Vera

### Flags:
      -f [input format file: ms, fasta OR tfa (gz file indexed)]
      -i [path and name of the input file]
      -o [output format file: 0 (extended),
                              1 (single line/window),
                              2 (single line SFS/window),
                              3 (dadi-like format),
                              4 (single line pairwise distribution)
                              5 (single line freq. variant per line/window)
                              6 (SNP genotype matrix)
                              7 (SweepFiinder format -only first pop-)
                              8 (single line/window: Frequency of each haplotype in the populations)
                              9 (single line/window: Frequency of variants per line and population)
                              92 (single line/window: -rSFS- Frequency of variants per population relative to all)
                             10 (full extended)]
      -N [#_pops] [#samples_pop1] ... [#samples_popN]
      -n [name of the file containing the name(s) of scaffold(s) and their length (separated by a tab), one per line (ex. fai file)]
      -T [path and name of the output file]. DEFAULT stdout.
##### OPTIONAL GENERAL PARAMETERS:
      -G [outgroup (0/1)] (last population). DEFAULT 0.
      -u [include unknown positions (0/1)].  DEFAULT 0.
      -A [Alternative Spectrum File (Only for Optimal Test): alternative_spectrum for each population (except outg)
          File format: (average absolute values) header plus fr(0,1) fr(0,2) ... fr(0,n-1) theta(0)/nt,
          fr(1,1) fr(1,2) ... fr(1,n-1) theta(1)/nt...]
      -S [Null Spectrum File (only if -A is defined): null_spectrum for each population (except outg).
          (average absolute values) header plus fr(0,1) fr(0,2) ... fr(0,n-1) theta(0)/nt,
          fr(1,1) fr(1,2) ... fr(1,n-1) theta(1)/nt...]. DEFAULT SNM.
    Optional Parameters for fasta and tfa input files:
      -O [#_nsam] [number order of first sample, number 0 is the first sample] [second sample] ...etc. up to nsamples.
         DEFAULT current order.
      -t [# permutations per window (H0: Fst=0). Only available with option -u 0]. DEFAULT 0.
      -s [seed]. DEFAULT 123456.
##### PARAMETERS FOR TFASTA INPUT (-f tfa): 'SLIDING WINDOW ANALYSIS OF EMPIRICAL DATA'
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
##### PARAMETERS FOR MS INPUT (-f ms):'SIMULATION ANALYSIS OF A SINGLE REGION'
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
  
  Conversion example 

```bash
  tfa_index  ./Examples/V0.1.0/100Kchr10.tfa.gz -o ./Examples/V1.0.0//100Kchr10.tfa.gz
```

Create an index example
```bash
  tfa_index  ./Examples/V1.0.0/100Kchr10.tfa.gz
  # use -f to force create the index file if it already exists
  tfa_index  ./Examples/V1.0.0/100Kchr10.tfa.gz -f
```