#mstatspop v.0.1beta (20170627)
#Sebastian E. Ramos-Onsins, Luca Ferretti, Emanuele Raineri, Giacomo Marmorini, William Burgos and Gonzalo Vera
#Variability Analyses of multiple populations: Calculation and estimation of statistics and neutrality tests.
#Flags:
#      -f [input format file: fasta, tfa, ms]
#      -i [path and name of the input file]
#      -o [output format file: 0 (extended),
#                              1 (single line/window),
#                              2 (single line SFS/window),
#                              3 (dadi-like format),
#                              4 (single line pairwise distribution)
#                              5 (single line freq. variant per line/window)
#                              6 (SNP genotype matrix)
#                             10 (Full extended)]
#      -N [#_pops] [#samples_pop1] ... [#samples_popN]
#      -n [name of a single scaffold to analyze. For tfa can be a list separated by commas(ex. -n chr1,chr2,chr3]
#   OPTIONAL GENERAL PARAMETERS:
#      -G [outgroup (0/1)] (last population). DEFAULT 0.
#      -u [include unknown positions (0/1)].  DEFAULT 0.
#      -T [path and name of the output file]. DEFAULT stdout.
#      -A [Alternative Spectrum File (Only for Optimal Test): alternative_spectrum for each population (except outg)
#          File format: (average absolute values) header plus fr(0,1) fr(0,2) ... fr(0,n-1) theta(0)/nt,
#          fr(1,1) fr(1,2) ... fr(1,n-1) theta(1)/nt...]
#      -S [Null Spectrum File (only if -a is defined): null_spectrum for each population (except outg).
#          (average absolute values) header plus fr(0,1) fr(0,2) ... fr(0,n-1) theta(0)/nt,
#          fr(1,1) fr(1,2) ... fr(1,n-1) theta(1)/nt...]. DEFAULT SNM.
#      -P [Only for Calculation of R2_p: first value is the number of values to include, next are the number of lines to consider. ex: -P 6 1 2 4 8 16 64]
#    Optional Parameters for fasta and tfa input files:
#      -O [#_nsam] [number order of first sample, number 0 is the first sample] [second sample] ...etc. up to nsamples.
#         DEFAULT current order.
#      -t [# permutations per window (H0: Fst=0). Only available with option -u 0]. DEFAULT 0.
#      -s [seed]. DEFAULT 123456.
#   PARAMETERS FOR TFASTA INPUT (-f tfa): "SLIDING WINDOW ANALYSIS OF EMPIRICAL DATA"
#      -w [window size].
#    Optional:
#      -z [slide size (must be equal or larger than window size)]. DEFAULT window size.
#      -Z [first window size displacement [for comparing overlapped windows])]. DEFAULT 0.
#      -Y [define window lengths in 'physical' positions (1) or in 'effective' positions (0)]. DEFAULT 1.
#      -W [file with the coordinates of each window [init end] (overwrite options -w and -z).
#         DEFAULT one whole window.
#      -E [input file with weights for positions:
#         include three columns with a header,
#         first the physical positions (1...end),
#         second the weight for positions and
#         third a boolean weight for the variant (eg. syn variant in nsyn counts is 0.000)].
#         DEFAULT all 1.000
#   PARAMETERS FOR MS INPUT (-f ms): "SIMULATION ANALYSIS OF A SINGLE REGION"
#      -l [length]
#    Optional:
#      -r [# ms iterations]. DEFAULT 1.
#      -m [include mask_filename] DEFAULT -1 (all positions included).
#         [mask_file format: 1st row with 'length' weights, next sample rows x lengths: missing 0, sequenced 1)].
#         DEFAULT no mask.
#      -v [ratio transitions/transversions]. DEFAULT 0.5.
#      -F [force analysis to include outgroup (0/1) (0 in ms means ancestral)]. DEFAULT 0.
#      -q [frequency of reverted mutation] (only with -F 1). DEFAULT 0.
#   PARAMETERS FOR FASTA INPUT (-f fasta): "WHOLE REGION ANALYSIS"
#    Optional:
#      -p [Number of lineages per sequence (1/2)]. DEFAULT 1.
#      -g [GFF_file]
#         [add also: coding,noncoding,synonymous,nonsynonymous,silent,others(whatever annotated)]
#         [if 'synonymous' or 'nonsynonymous' add: Genetic_Code: Nuclear_Universal,mtDNA_Drosophila,mtDNA_Mammals,Other]
#         [if 'Other', introduce the code for the 64 triplets in the order UUU UUC UUA UUG ... etc.].
#         DEFAULT no annotation.
#      -K [make a MASK file with the valid positions for this fasta. Useful for running ms simulations (1/0)]. DEFAULT 0   HELP:
#      -c [in case use coding regions, criteria to consider transcripts (max/min/first/long)]. DEFAULT long.
#   HELP:
#      -h [help and exit]
#      
#Validation: test a number of combinations.
#-f fasta tfa ms
#-i ./100Kchr10.fa
#-i ./100Kchr10.tfa
#-i ./100Kchr10_fa2ms_03.ms.txt
#-g ./100Kchr10.gtf nonsynonymous Nuclear_Universal
#-W ./coord_100Kb.txt 
#-E ./100Kchr10_fa2tfa_02.tfa_nonsynonymous_max_IncludeMissing_NOoutg_ploidy1_WEIGHTS.txt
#-o 1
#-N 2 40 2 
#-N 3 20 20 2
#options: [-G -u -T]
#-G 1/0 -u 1/0

# Check file "compile_mstatspop.sh" to compile and generate "mstatspop" binary.
# To compile "ms" program, just run:
gcc -o ./bin/ms ./sources_msHudson/ms.c ./sources_msHudson/streec.c ./sources_msHudson/rand2.c -lm -O3

cd ./Examples

echo
echo flags included in mstatspop_help.txt
../bin/mstatspop -h > ../mstatspop_help.txt

echo --------------------------------------------------------------------------------------------------
echo PARAMETERS FOR FASTA INPUT -f fasta: [-p -g -c -O -t -s] WHOLE REGION ANALYSIS
echo --------------------------------------------------------------------------------------------------
echo 
echo Example fa.01.txt
echo ../bin/mstatspop -f fasta -i ./100Kchr10.fa -o 0 -N 1 42   -T ./mstatspop_100chr10.fa.01.txt -K 1 -n chr10
../bin/mstatspop -f fasta -i ./100Kchr10.fa -o 0 -N 1 42   -T ./mstatspop_100chr10.fa.01.txt -K 1 -n chr10
echo 
echo Example fa.01b.txt: Same result but in a single line
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 1 -N 1 42   -T ./mstatspop_100chr10.fa.01b.txt -K 1 -n chr10
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 1 -N 1 42   -T ./mstatspop_100chr10.fa.01b.txt -K 1 -n chr10
echo 
echo Example fa.02.txt: Two populations
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ./mstatspop_100chr10.fa.02.txt -K 1 -n chr10
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ./mstatspop_100chr10.fa.02.txt -K 1 -n chr10
echo 
echo Example fa.03.txt: Two pops + outgroup
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ./mstatspop_100chr10.fa.03.txt -G 1 -K 1 -n chr10
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ./mstatspop_100chr10.fa.03.txt -G 1 -K 1 -n chr10
echo 
echo Example fa.04.txt: Two pops + outg + missing
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ./mstatspop_100chr10.fa.04.txt -G 1 -u 1 -K 1 -n chr10
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ./mstatspop_100chr10.fa.04.txt -G 1 -u 1 -K 1 -n chr10
echo 
echo Example fa.05.txt: Three pops + outg + missing
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.fa.05.txt -G 1 -u 1 -K 1 -n chr10
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.fa.05.txt -G 1 -u 1 -K 1 -n chr10
echo 
echo Example fa.06.txt: Two pops + outg + missing + IUPAC diploid
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ./mstatspop_100chr10.fa.06.txt -G 1 -u 1 -p 2 -K 1 -n chr10
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ./mstatspop_100chr10.fa.06.txt -G 1 -u 1 -p 2 -K 1 -n chr10
echo 
echo Example fa.07.txt: Two pops + outg + missing + GTF2 for nonsyn
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ./mstatspop_100chr10.fa.07.txt -G 1 -u 1 -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -K 1 -n chr10
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ./mstatspop_100chr10.fa.07.txt -G 1 -u 1 -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -K 1 -n chr10
echo 
echo Example fa.08.txt: Three pops + outg + missing + permutation test
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.fa.08.txt -G 1 -u 1 -t 1000 -s 1684 -K 1 -n chr10
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.fa.08.txt -G 1 -u 1 -t 1000 -s 1684 -K 1 -n chr10
echo 
echo Example fa.09.txt: Three pops + outg + missing + reordering the two pops in inverse order
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.fa.09.txt -G 1 -K 1-u 1 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41 -n chr10
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.fa.09.txt -G 1 -u 1 -K 1 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41 -n chr10
echo 
echo Example fa.10.txt: Three pops + outg + missing + GTF2 for nonsyn
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.fa.10.txt -G 1 -u 1  -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -K 1 -n chr10
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.fa.10.txt -G 1 -u 1  -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -K 1 -n chr10
echo 
echo Example fa.11.txt: Three pops + outg + missing + reordering the two pops in inverse order + GTF2 for nonsyn
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.fa.11.txt -G 1 -K 0 -u 1 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41  -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -n chr10
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.fa.11.txt -G 1 -u 1 -K 0 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41  -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -n chr10
echo 
echo --------------------------------------------------------------------------------------------------
echo PARAMETERS FOR MS INPUT -f ms: -l [-r -m -F -v -q] SIMULATION ANALYSIS OF A SINGLE REGION
echo --------------------------------------------------------------------------------------------------
echo 
echo Example ms.00.txt
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_02.ms.txt -o 0 -N 1 42 -T ./mstatspop_100chr10.ms.00.txt -l 100000 -n chr10
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_02.ms.txt -o 0 -N 1 42 -T ./mstatspop_100chr10.ms.00.txt -l 100000 -n chr10
echo 
echo Example ms.01.txt
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_02.ms.txt -o 0 -N 1 42 -T ./mstatspop_100chr10.ms.01.txt -l 100000 -m ./100Kchr10_fa2ms_02.ms.txt_npops1_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_02.ms.txt -o 0 -N 1 42 -T ./mstatspop_100chr10.ms.01.txt -l 100000 -m ./100Kchr10_fa2ms_02.ms.txt_npops1_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10
echo 
echo Example ms.01b.txt: Same result but in a single line
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_02.ms.txt -o 1 -N 1 42 -T ./mstatspop_100chr10.ms.01b.txt -l 100000 -m ./100Kchr10_fa2ms_02.ms.txt_npops1_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_02.ms.txt -o 1 -N 1 42 -T ./mstatspop_100chr10.ms.01b.txt -l 100000 -m ./100Kchr10_fa2ms_02.ms.txt_npops1_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10
echo 
echo Example ms.02.txt: Two populations
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_03.ms.txt -o 0 -N 2 40 2 -T ./mstatspop_100chr10.ms.02.txt -l 100000 -m ./100Kchr10_fa2ms_07.ms.txt_npops2_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_03.ms.txt -o 0 -N 2 40 2 -T ./mstatspop_100chr10.ms.02.txt -l 100000 -m ./100Kchr10_fa2ms_07.ms.txt_npops2_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10
echo 
echo Example ms.03.txt: Two pops + outgroup
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_03.ms.txt -o 0 -N 2 40 2 -T ./mstatspop_100chr10.ms.03.txt -l 100000 -G 1 -m ./100Kchr10_fa2ms_07.ms.txt_npops2_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_03.ms.txt -o 0 -N 2 40 2 -T ./mstatspop_100chr10.ms.03.txt -l 100000 -G 1 -m ./100Kchr10_fa2ms_07.ms.txt_npops2_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10
echo 
echo Example ms.04.txt: One pop + several ms iterations
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_07.ms.txt -o 1 -N 1 42 -T ./mstatspop_100chr10.ms.04.txt -l 10000 -r 5 -n chr10
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_07.ms.txt -o 1 -N 1 42 -T ./mstatspop_100chr10.ms.04.txt -l 10000 -r 5 -n chr10
echo 
echo Example ms.05.txt: Two pops + several ms iterations
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ./mstatspop_100chr10.ms.05.txt -l 10000 -r 5 -n chr10
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ./mstatspop_100chr10.ms.05.txt -l 10000 -r 5 -n chr10
echo 
echo Example ms.06.txt: Two pops + several ms iterations + forcing outgroup the value '0'
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ./mstatspop_100chr10.ms.06.txt -l 10000 -r 5 -F 1 -n chr10
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ./mstatspop_100chr10.ms.06.txt -l 10000 -r 5 -F 1 -n chr10
echo 
echo Example ms.07.txt: Two pops + several ms iterations  + forcing outgroup the value '0' + revert mutation 
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ./mstatspop_100chr10.ms.07.txt -l 10000 -r 5 -F 1 -q 0.01 -n chr10
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ./mstatspop_100chr10.ms.07.txt -l 10000 -r 5 -F 1 -q 0.01 -n chr10
echo 
echo Example ms.08.txt: Two pops + outg + several ms iterations + s/v ratio
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ./mstatspop_100chr10.ms.08.txt -l 10000 -r 5 -G 1 -v 2.0 -n chr10
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ./mstatspop_100chr10.ms.08.txt -l 10000 -r 5 -G 1 -v 2.0 -n chr10
echo 
echo Example ms.09.txt: Two pops +  missing
../bin/ms 42 1 -t 1000 -I 2 40 2 -ej 2.0 2 1 -p 9 > ./100Kchr10_simulation.ms.txt
echo ../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -T ./mstatspop_100chr10.ms.09.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_07.ms.txt_npops2_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10
../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -T ./mstatspop_100chr10.ms.09.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_07.ms.txt_npops2_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10
echo 
echo Example ms.10.txt: Two pops + missing + weight for nonsyn
echo ../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -T ./mstatspop_100chr10.ms.10.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_04.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10
../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -T ./mstatspop_100chr10.ms.10.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_04.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10
echo 
echo Example ms.11.txt: Two pops + missing + weight for nonsyn + permutation test
echo ../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -T ./mstatspop_100chr10.ms.11.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_04.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -t 1000 -s 1684 -n chr10
../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -T ./mstatspop_100chr10.ms.11.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_04.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -t 1000 -s 1684 -n chr10
echo 
echo Example ms.12.txt: Two pops + outg + missing + weight for nonsyn + permutation test + s/v ratio
echo ../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -G 1 -T ./mstatspop_100chr10.ms.12.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_04b.ms.txt_npops2_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_outg_ploidy1_MASK.txt -t 1000 -s 1684 -v 2.0 -n chr10
../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -G 1 -T ./mstatspop_100chr10.ms.12.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_04b.ms.txt_npops2_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_outg_ploidy1_MASK.txt -t 1000 -s 1684 -v 2.0 -n chr10
echo 
echo Example ms.13.txt: Two pops + outg + missing + weight for nonsyn + permutation test + s/v ratio + 1000 replicates
../bin/ms 42 1000 -t 1000 -I 2 40 2 -ej 2.0 2 1 -p 9 > ./100Kchr10_simulation1000.ms.txt
echo ../bin/mstatspop  -f ms -i ./100Kchr10_simulation1000.ms.txt -o 1 -N 2 40 2 -G 1 -T ./mstatspop_100chr10.ms.13.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_04b.ms.txt_npops2_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_outg_ploidy1_MASK.txt -t 1000 -s 1684 -v 2.0 -r 1000 -n chr10
../bin/mstatspop  -f ms -i ./100Kchr10_simulation1000.ms.txt -o 1 -N 2 40 2 -G 1 -T ./mstatspop_100chr10.ms.13.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_04b.ms.txt_npops2_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_outg_ploidy1_MASK.txt -t 1000 -s 1684 -v 2.0 -r 1000 -n chr10
echo 
echo --------------------------------------------------------------------------------------------------
echo PARAMETERS FOR TFASTA INPUT -f tfa: -w [-z -Y -W -E -O -t -s] SLIDING WINDOW ANALYSIS OF EMPIRICAL DATA
echo --------------------------------------------------------------------------------------------------
echo 
echo Example tfa.01.txt
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 1 42   -T ./mstatspop_100chr10.tfa.01.txt -w 100000 -n chr10
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 1 42   -T ./mstatspop_100chr10.tfa.01.txt -w 100000 -n chr10
echo 
echo Example tfa.01b.txt: Same result but in a single line
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 1 42   -T ./mstatspop_100chr10.tfa.01b.txt -w 100000 -n chr10
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 1 42   -T ./mstatspop_100chr10.tfa.01b.txt -w 100000 -n chr10
echo 
echo Example tfa.01c.txt: Sliding window of 1000 in lines 
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 1 42   -T ./mstatspop_100chr10.tfa.01c.txt -w 1000 -n chr10
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 1 42   -T ./mstatspop_100chr10.tfa.01c.txt -w 1000 -n chr10
echo 
echo Example tfa.02.txt: Two populations
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100chr10.tfa.02.txt -w 100000 -n chr10
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100chr10.tfa.02.txt -w 100000 -n chr10
echo 
echo Example tfa.03.txt: Two pops + outgroup
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100chr10.tfa.03.txt -G 1 -w 100000 -n chr10
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100chr10.tfa.03.txt -G 1 -w 100000 -n chr10
echo 
echo Example tfa.04.txt: Two pops + outg + missing
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100chr10.tfa.04.txt -G 1 -u 1 -w 100000 -n chr10
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100chr10.tfa.04.txt -G 1 -u 1 -w 100000 -n chr10
echo
echo Example tfa.06.txt: Two pops + outg + missing + window 10000
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ./mstatspop_100chr10.tfa.06.txt -G 1 -u 1 -w 10000 -n chr10
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ./mstatspop_100chr10.tfa.06.txt -G 1 -u 1 -w 10000 -n chr10
echo 
echo Example tfa.07.txt: Two pops + outg + missing + window 10000 + slide 20000
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ./mstatspop_100chr10.tfa.07.txt -G 1 -u 1 -w 10000 -z 20000 -n chr10
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ./mstatspop_100chr10.tfa.07.txt -G 1 -u 1 -w 10000 -z 20000 -n chr10
echo 
echo Example tfa.08.txt: Two pops + outg + missing + coordinates. Equal than the previous example
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100chr10.tfa.08.txt -G 1 -u 1 -W ./coord_100Kb.txt -n chr10
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100chr10.tfa.08.txt -G 1 -u 1 -W ./coord_100Kb.txt -n chr10
echo 
echo Example tfa.09.txt: Two pops + outg + missing + window 10000 + weight for nonsyn
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100chr10.tfa.09.txt -G 1 -u 1 -w 10000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100chr10.tfa.09.txt -G 1 -u 1 -w 10000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10
echo 
echo Example tfa.10.txt: Two pops + outg + missing + window 1000 + weight for nonsyn + Effective_positions
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100chr10.tfa.10.txt -G 1 -u 1 -w 1000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100chr10.tfa.10.txt -G 1 -u 1 -w 1000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10
echo 
echo Example tfa.10b.txt: Two pops + outg + missing + window 1000 + weight for nonsyn + Effective_positions + single line results
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ./mstatspop_100chr10.tfa.10b.txt -G 1 -u 1 -w 1000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ./mstatspop_100chr10.tfa.10b.txt -G 1 -u 1 -w 1000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10
echo 
echo Example tfa.10c.txt: Two pops + outg + missing + window 1000 + 500 displacement + weight for nonsyn + Effective_positions + single line results
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ./mstatspop_100chr10.tfa.10c.txt -G 1 -u 1 -w 1000 -z 500 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ./mstatspop_100chr10.tfa.10c.txt -G 1 -u 1 -w 1000 -z 500 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10
echo 
echo Example tfa.10d.txt: Two pops + outg + missing + window 1000 + 2000 displacement + weight for nonsyn + Effective_positions + single line results
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ./mstatspop_100chr10.tfa.10d.txt -G 1 -u 1 -w 1000 -z 2000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ./mstatspop_100chr10.tfa.10d.txt -G 1 -u 1 -w 1000 -z 2000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10
echo 
echo Example tfa.11.txt: Two pops + outg + missing + weight for nonsyn + Effective_positions
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100chr10.tfa.11.txt -G 1 -u 1 -w 100000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100chr10.tfa.11.txt -G 1 -u 1 -w 100000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10
echo 
echo Example tfa.12.txt: Three pops + outg + missing + reordering the two pops in inverse order
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.tfa.12.txt -G 1 -u 1 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41 -w 100000 -n chr10
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.tfa.12.txt -G 1 -u 1 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41 -w 100000 -n chr10
echo 
echo Example tfa.13.txt: Three pops + outg + missing + weight for nonsyn
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.tfa.13.txt -G 1 -u 1  -w 100000 -E ./100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.tfa.13.txt -G 1 -u 1  -w 100000 -E ./100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10
echo 
echo Example tfa.14.txt: Three pops + outg + missing + reordering the two pops in inverse order + weight for nonsyn
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.tfa.14.txt -G 1 -u 1  -w 100000 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41  -E ./100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.tfa.14.txt -G 1 -u 1  -w 100000 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41  -E ./100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10
echo 
echo Example tfa.14b.txt: Three pops + outg + missing + reordering the two pops in inverse order + weight for nonsyn + Effective positions
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.tfa.14b.txt -G 1 -u 1  -w 100000 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41  -E ./100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ./mstatspop_100chr10.tfa.14b.txt -G 1 -u 1  -w 100000 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41  -E ./100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10
echo 
echo Example tfa.15.txt: Two pops + outg + missing + window 10000 + slide 20000 + permutations: like example 07
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ./mstatspop_100chr10.tfa.15.txt -G 1 -u 1 -w 10000 -z 20000 -t 1000 -s 1684 -n chr10
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ./mstatspop_100chr10.tfa.15.txt -G 1 -u 1 -w 10000 -z 20000 -t 1000 -s 1684 -n chr10
echo
echo Example tfa.16.txt: Two pops + outg + missing + window 100 + slide 100 + permutations:
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ./mstatspop_100chr10.tfa.16.txt -G 1 -u 1 -w 100 -z 100 -t 1000 -s 1684 -n chr10
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ./mstatspop_100chr10.tfa.16.txt -G 1 -u 1 -w 100 -z 100 -t 1000 -s 1684 -n chr10
echo 
echo Example tfa.17.txt: Three pops + outg + missing + window 100 + slide 100 + permutations:
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 3 20 20 2 -T ./mstatspop_100chr10.tfa.17.txt -G 1 -u 1 -w 100 -z 100 -t 1000 -s 1684 -n chr10
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 3 20 20 2 -T ./mstatspop_100chr10.tfa.17.txt -G 1 -u 1 -w 100 -z 100 -t 1000 -s 1684 -n chr10
echo
echo Example tfa.18.txt: Three pops LESS SAMPLES + outg + missing + window 100 + slide 100 :
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 3 10 10 2 -T ./mstatspop_100chr10.tfa.18.txt -G 1 -u 1 -w 100 -z 100 -n chr10
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 3 10 10 2 -T ./mstatspop_100chr10.tfa.18.txt -G 1 -u 1 -w 100 -z 100 -n chr10
echo 
echo Example tfa.19.txt: Three pops LESS SAMPLES + reordering + outg + missing + window 100 + slide 100 :
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 3 10 10 2 -T ./mstatspop_100chr10.tfa.19.txt -G 1 -u 1 -w 100 -z 100 -O 22 20 21 22 23 24 25 26 27 28 29 0 1 2 3 4 5 6 7 8 9 40 41 -n chr10
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 3 10 10 2 -T ./mstatspop_100chr10.tfa.19.txt -G 1 -u 1 -w 100 -z 100 -O 22 20 21 22 23 24 25 26 27 28 29 0 1 2 3 4 5 6 7 8 9 40 41 -n chr10
echo 
echo --------------------------------------------------------------------------------------------------
echo Examples using MULTICHROMOSOME TFA FILE
echo --------------------------------------------------------------------------------------------------
echo
echo Three different chromosme analyses
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 1 42   -T ./mstatspop_100allchr.tfa.01.txt   -w 100000 -n chr10,chr12,chr14
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 1 42   -T ./mstatspop_100allchr.tfa.01.txt   -w 100000 -n chr10,chr12,chr14
echo
echo Same as before but single line output
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 1 42   -T ./mstatspop_100allchr.tfa.01b.txt   -w 100000 -n chr10,chr12,chr14
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 1 42   -T ./mstatspop_100allchr.tfa.01b.txt   -w 100000 -n chr10,chr12,chr14
echo
echo Same as before but using 1000bp sliding window
 ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 1 42   -T ./mstatspop_100allchr.tfa.01c.txt -w 1000 -n chr10,chr12,chr14
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 1 42   -T ./mstatspop_100allchr.tfa.01c.txt -w 1000 -n chr10,chr12,chr14
echo
echo Same but using two populations
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100allchr.tfa.02.txt   -w 100000 -n chr10,chr12,chr14
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100allchr.tfa.02.txt   -w 100000 -n chr10,chr12,chr14
echo
echo Same but using two populations, last is outgroup
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -G 1 -T ./mstatspop_100allchr.tfa.03.txt   -w 100000 -n chr10,chr12,chr14
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -G 1 -T ./mstatspop_100allchr.tfa.03.txt   -w 100000 -n chr10,chr12,chr14
echo
echo Same but including missing data
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -G 1 -T ./mstatspop_100allchr.tfa.04.txt  -u 1 -w 100000 -n chr10,chr12,chr14
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -G 1 -T ./mstatspop_100allchr.tfa.04.txt  -u 1 -w 100000 -n chr10,chr12,chr14
echo
echo Same but only including chr10 and chr12
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -G 1 -T ./mstatspop_100allchr.tfa.05.txt  -u 1 -w 100000 -n chr10,chr12
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -G 1 -T ./mstatspop_100allchr.tfa.05.txt  -u 1 -w 100000 -n chr10,chr12
echo
echo Same but in single line output
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 2 40 2 -G 1 -T ./mstatspop_100allchr.tfa.06.txt  -u 1 -w 10000 -n chr10,chr12,chr14
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 2 40 2 -G 1 -T ./mstatspop_100allchr.tfa.06.txt  -u 1 -w 10000 -n chr10,chr12,chr14
echo
echo Same but using displacement of 20000bp
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 2 40 2 -G 1 -T ./mstatspop_100allchr.tfa.07.txt  -u 1 -z 20000 -w 10000 -n chr10,chr12,chr14
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 2 40 2 -G 1 -T ./mstatspop_100allchr.tfa.07.txt  -u 1 -z 20000 -w 10000 -n chr10,chr12,chr14
echo
echo Same but using coordinates file instead sliding window
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 2 40 2 -G 1 -T ./mstatspop_100allchr.tfa.08c1.txt  -u 1 -W ./coord_100Kb_allchr.txt -n chr10,chr12,chr14
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 2 40 2 -G 1 -T ./mstatspop_100allchr.tfa.08c1.txt  -u 1 -W ./coord_100Kb_allchr.txt -n chr10,chr12,chr14
echo
echo Use the weight file to calculate nonsuynonymous variability
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100allchr.tfa.09.txt -u 1 -w 10000 -E ./100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10,chr12,chr14
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100allchr.tfa.09.txt -u 1 -w 10000 -E ./100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10,chr12,chr14
echo
echo same but using outgroup weights
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100allchr.tfa.10.txt -G 1 -u 1 -w 10000 -E ./100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10,chr12,chr14
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100allchr.tfa.10.txt -G 1 -u 1 -w 10000 -E ./100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10,chr12,chr14
echo
echo same but with 100bp sliding window
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100allchr.tfa.11.txt -G 1 -u 1 -w 100 -E ./100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10,chr12,chr14
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ./mstatspop_100allchr.tfa.11.txt -G 1 -u 1 -w 100 -E ./100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10,chr12,chr14
echo
echo three pops and different sample order
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 3 10 10 2 -T ./mstatspop_100allchr.tfa.12.txt -G 1 -u 1 -w 100 -z 100 -O 22 20 21 22 23 24 25 26 27 28 29 0 1 2 3 4 5 6 7 8 9 40 41 -n chr10,chr12,chr14
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 3 10 10 2 -T ./mstatspop_100allchr.tfa.12.txt -G 1 -u 1 -w 100 -z 100 -O 22 20 21 22 23 24 25 26 27 28 29 0 1 2 3 4 5 6 7 8 9 40 41 -n chr10,chr12,chr14
echo
echo --------------------------------------------------------------------------------------------------
echo Collect some columns from the mstatspop output file:
echo --------------------------------------------------------------------------------------------------
echo 
echo 'perl ./bin/collect_data_columns.pl -in ./mstatspop_100chr10.tfa.15.txt -fc ./required_columns1.txt > ./mstatspop_100chr10.tfa.15_COLUMNS.txt'
perl ../bin/collect_data_columns.pl -in ./mstatspop_100chr10.tfa.15.txt -fc ./required_columns1.txt > ./mstatspop_100chr10.tfa.15_COLUMNS.txt
echo
echo --------------------------------------------------------------------------------------------------
echo Optimal tests using non-missing data:
echo --------------------------------------------------------------------------------------------------
echo 
echo  ../bin/mstatspop -f fasta -i ./MC1R_PigsOutg_aligned.fas -o 0 -p 1 -u 0 -t 1000 -s 123456 -G 1 -N 3 48 46 1 -A ./MC1R_H1frq.txt -g ./MC1R.gff nonsynonymous Nuclear_Universal -T ./MC1R_PigsOutg_NSyn_Opttest.txt -n MC1R-1
 ../bin/mstatspop -f fasta -i ./MC1R_PigsOutg_aligned.fas -o 0 -p 1 -u 0 -t 1000 -s 123456 -G 1 -N 3 48 46 1 -A ./MC1R_H1frq.txt -g ./MC1R.gff nonsynonymous Nuclear_Universal -T ./MC1R_PigsOutg_NSyn_Opttest.txt -n MC1R-1
echo
echo
