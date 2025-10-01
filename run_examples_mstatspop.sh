#To compile:
#for osx: look at library versions and modify
brew install gsl
brew install htslib
brew install zlib
sh build.sh
#OR
#gcc  -lgsl -lgslcblas -lm -lhts -lz -o ./bin/mstatspop -Wall -DinGSL=1 -O3 ./mstatspop/zutil.c ./mstatspop/files_util.c ./mstatspop/freq_stats.c ./mstatspop/fstcalc.c ./mstatspop/get_msdata.c ./mstatspop/get_obsdata.c ./mstatspop/get_obsdatastats.c ./mstatspop/get_tfadata.c ./mstatspop/jointfreqdist.c ./mstatspop/log.c ./mstatspop/mismatch.c ./mstatspop/missing_freqs.c ./mstatspop/mstatspop.c ./mstatspop/optimal-RinfR0.c ./mstatspop/permute.c ./mstatspop/print_output.c ./mstatspop/ran1.c ./mstatspop/sancestral.c ./mstatspop/tfasta.c ./mstatspop/tn93.c ./mstatspop/usegff.c ./mstatspop/util.c ./mstatspop/zindex.c ./mstatspop/calc_Toptimal_tests.c ./mstatspop/calcFs.c ./mstatspop/calcR2.c -I/opt/homebrew/Cellar/gsl/2.8/include  -I/opt/homebrew/Cellar/htslib/1.21/include -L/opt/homebrew/Cellar/gsl/2.8/lib -L/opt/homebrew/Cellar/htslib/1.21/lib
#for linux:
#sh build.sh
#gcc  -lgsl -lgslcblas -lm -lhts -lz -o ./bin/mstatspop -Wall -DinGSL=1 -O3 ./mstatspop/zutil.c ./mstatspop/files_util.c ./mstatspop/freq_stats.c ./mstatspop/fstcalc.c ./mstatspop/get_msdata.c ./mstatspop/get_obsdata.c ./mstatspop/get_obsdatastats.c ./mstatspop/get_tfadata.c ./mstatspop/jointfreqdist.c ./mstatspop/log.c ./mstatspop/mismatch.c ./mstatspop/missing_freqs.c ./mstatspop/mstatspop.c ./mstatspop/optimal-RinfR0.c ./mstatspop/permute.c ./mstatspop/print_output.c ./mstatspop/ran1.c ./mstatspop/sancestral.c ./mstatspop/tfasta.c ./mstatspop/tn93.c ./mstatspop/usegff.c ./mstatspop/util.c ./mstatspop/zindex.c ./mstatspop/calc_Toptimal_tests.c ./mstatspop/calcFs.c ./mstatspop/calcR2.c
#to compile Hudson's ms
#gcc -o ./bin/ms ./sources_msHudson/ms.c ./sources_msHudson/streec.c ./sources_msHudson/rand2.c -lm -O3

cd ./Examples

echo
echo flags included in mstatspop_help.txt
../bin/mstatspop -h > ../mstatspop_help.txt

echo --------------------------------------------------------------------------------------------------
echo PARAMETERS FOR FASTA INPUT -f fasta: [-p -g -c -O  -T -s] WHOLE REGION ANALYSIS
echo --------------------------------------------------------------------------------------------------
echo 
echo Example fa.01.txt
echo ../bin/mstatspop -f fasta -i ./100Kchr10.fa -o 0 -N 1 42    -T ../Results/mstatspop_100chr10.fa.01.txt -K 1 -n chr10.txt
../bin/mstatspop -f fasta -i ./100Kchr10.fa -o 0 -N 1 42   -T ../Results/mstatspop_100chr10.fa.01.txt -K 1 -n chr10.txt
echo 
echo Example fa.01b.txt: Same result but in a single line
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 1 -N 1 42   -T ../Results/mstatspop_100chr10.fa.01b.txt -K 1 -n chr10.txt
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 1 -N 1 42   -T ../Results/mstatspop_100chr10.fa.01b.txt -K 1 -n chr10.txt
echo 
echo Example fa.02.txt: Two populations
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.fa.02.txt -K 1 -n chr10.txt
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.fa.02.txt -K 1 -n chr10.txt
echo 
echo Example fa.03.txt: Two pops + outgroup
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.fa.03.txt -G 1 -K 1 -n chr10.txt
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.fa.03.txt -G 1 -K 1 -n chr10.txt
echo 
echo Example fa.04.txt: Two pops + outg + missing
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.fa.04.txt -G 1 -u 1 -K 1 -n chr10.txt
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.fa.04.txt -G 1 -u 1 -K 1 -n chr10.txt
echo 
echo Example fa.05.txt: Three pops + outg + missing
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.fa.05.txt -G 1 -u 1 -K 1 -n chr10.txt
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.fa.05.txt -G 1 -u 1 -K 1 -n chr10.txt
echo 
echo Example fa.06.txt: Two pops + outg + missing + IUPAC diploid
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.fa.06.txt -G 1 -u 1 -p 2 -K 1 -n chr10.txt
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.fa.06.txt -G 1 -u 1 -p 2 -K 1 -n chr10.txt
echo 
echo Example fa.07.txt: Two pops + outg + missing + GTF2 for nonsyn
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.fa.07.txt -G 1 -u 1 -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -K 1 -n chr10.txt
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.fa.07.txt -G 1 -u 1 -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -K 1 -n chr10.txt
echo 
echo Example fa.08.txt: Three pops + outg + permutation test
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.fa.08.txt -G 1 -u 0 -t 1000 -s 1684 -K 1 -n chr10.txt
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.fa.08.txt -G 1 -u 0 -t 1000 -s 1684 -K 1 -n chr10.txt
echo
echo Example fa.09.txt: Three pops + outg + missing + reordering the two pops in inverse order
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.fa.09.txt -G 1 -K 1-u 1 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41 -n chr10.txt
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.fa.09.txt -G 1 -u 1 -K 1 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41 -n chr10.txt
echo 
echo Example fa.10.txt: Three pops + outg + missing + GTF2 for nonsyn
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.fa.10.txt -G 1 -u 1  -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -K 1 -n chr10.txt
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.fa.10.txt -G 1 -u 1  -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -K 1 -n chr10.txt
echo 
echo Example fa.11.txt: Three pops + outg + missing + reordering the two pops in inverse order + GTF2 for nonsyn
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.fa.11.txt -G 1 -K 0 -u 1 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41  -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -n chr10.txt
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.fa.11.txt -G 1 -u 1 -K 0 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41  -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -n chr10.txt
echo
echo Example fa.12.txt: One pop + rSFS + missing + GTF2 for nonsyn
echo ../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -R 1 -N 20 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -T ../Results/mstatspop_100chr10.fa.12.txt -G 0 -u 1  -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -K 1 -n chr10.txt
../bin/mstatspop  -f fasta -i ./100Kchr10.fa -o 0 -R 1 -N 20 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -T ../Results/mstatspop_100chr10.fa.12.txt -G 0 -u 1  -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -K 1 -n chr10.txt
echo
echo --------------------------------------------------------------------------------------------------
echo PARAMETERS FOR MS INPUT -f ms: [-r -m -F -v -q] SIMULATION ANALYSIS OF A SINGLE REGION
echo --------------------------------------------------------------------------------------------------
echo 
echo Example ms.00.txt
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_02.ms.txt -o 0 -N 1 42 -T ../Results/mstatspop_100chr10.ms.00.txt -l 100000 -n chr10.txt
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_02.ms.txt -o 0 -N 1 42 -T ../Results/mstatspop_100chr10.ms.00.txt -l 100000 -n chr10.txt
echo 
echo Example ms.01.txt
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_02.ms.txt -o 0 -N 1 42 -T ../Results/mstatspop_100chr10.ms.01.txt -l 100000 -m ./100Kchr10_npop1_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10.txt
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_02.ms.txt -o 0 -N 1 42 -T ../Results/mstatspop_100chr10.ms.01.txt -l 100000 -m ./100Kchr10_npop1_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10.txt
echo
echo Example ms.01b.txt: Same result but in a single line
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_02.ms.txt -o 1 -N 1 42 -T ../Results/mstatspop_100chr10.ms.01b.txt -l 100000 -m ./100Kchr10_npop1_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10.txt
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_02.ms.txt -o 1 -N 1 42 -T ../Results/mstatspop_100chr10.ms.01b.txt -l 100000 -m ./100Kchr10_npop1_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10.txt
echo 
echo Example ms.02.txt: Two populations
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_03.ms.txt -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.02.txt -l 100000 -m ./100Kchr10_npop2_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10.txt
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_03.ms.txt -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.02.txt -l 100000 -m ./100Kchr10_npop2_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10.txt
echo 
echo Example ms.03.txt: Two pops + outgroup
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_03.ms.txt -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.03.txt -l 100000 -G 1 -m ./100Kchr10_npop2_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10.txt
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_03.ms.txt -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.03.txt -l 100000 -G 1 -m ./100Kchr10_npop2_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10.txt
echo 
echo Example ms.04.txt: One pop + several ms iterations
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_07.ms.txt -o 1 -N 1 42 -T ../Results/mstatspop_100chr10.ms.04.txt -l 10000 -r 5 -n chr10.txt
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_07.ms.txt -o 1 -N 1 42 -T ../Results/mstatspop_100chr10.ms.04.txt -l 10000 -r 5 -n chr10.txt
echo 
echo Example ms.05.txt: Two pops + several ms iterations
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.05.txt -l 10000 -r 5 -n chr10.txt
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.05.txt -l 10000 -r 5 -n chr10.txt
echo 
echo Example ms.06.txt: Two pops + several ms iterations + forcing outgroup the value '0'
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.06.txt -l 10000 -r 5 -F 1 -n chr10.txt
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.06.txt -l 10000 -r 5 -F 1 -n chr10.txt
echo 
echo Example ms.07.txt: Two pops + several ms iterations  + forcing outgroup the value '0' + revert mutation 
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.07.txt -l 10000 -r 5 -F 1 -q 0.01 -n chr10.txt
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.07.txt -l 10000 -r 5 -F 1 -q 0.01 -n chr10.txt
echo 
echo Example ms.08.txt: Two pops + outg + several ms iterations + s/v ratio
echo ../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.08.txt -l 10000 -r 5 -G 1 -v 2.0 -n chr10.txt
../bin/mstatspop  -f ms -i ./100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.08.txt -l 10000 -r 5 -G 1 -v 2.0 -n chr10.txt
echo 
echo Example ms.09.txt: Two pops +  missing
../bin/ms 42 1 -t 1000 -I 2 40 2 -ej 2.0 2 1 -p 9 > ./100Kchr10_simulation.ms.txt
echo ../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.09.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_npop2_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10.txt
../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.09.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_npop2_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10.txt
echo
echo Example ms.10.txt: Two pops + missing + weight for nonsyn
echo ../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.10.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_04.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10.txt
../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.10.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_04.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10.txt
echo
echo Example ms.11.txt: Two pops + missing + weight for nonsyn + permutation test
echo ../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.11.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_04.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -t 1000 -s 1684 -n chr10.txt
../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.ms.11.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_04.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -t 1000 -s 1684 -n chr10.txt
echo
echo Example ms.12.txt: Two pops + outg + missing + weight for nonsyn + permutation test + s/v ratio
echo ../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -G 1 -T ../Results/mstatspop_100chr10.ms.12.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_04b.ms.txt_npops2_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_outg_ploidy1_MASK.txt -t 1000 -s 1684 -v 2.0 -n chr10.txt
../bin/mstatspop  -f ms -i ./100Kchr10_simulation.ms.txt -o 0 -N 2 40 2 -G 1 -T ../Results/mstatspop_100chr10.ms.12.txt -l 100000 -G 1 -u 1 -m ./100Kchr10_fa2ms_04b.ms.txt_npops2_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_outg_ploidy1_MASK.txt -t 1000 -s 1684 -v 2.0 -n chr10.txt
echo
echo Example ms.14.txt: One pop + rSFS + missing + weight for nonsyn
echo ../bin/mstatspop  -f ms -R 1 -i ./100Kchr10_simulation.ms.txt -o 0 -N 20 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -T ../Results/mstatspop_100chr10.ms.14.txt -l 100000 -G 0 -u 1 -m ./100Kchr10_fa2ms_04.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10.txt
../bin/mstatspop  -f ms -R 1 -i ./100Kchr10_simulation.ms.txt -o 0 -N 20 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -T ../Results/mstatspop_100chr10.ms.14.txt -l 100000 -G 0 -u 1 -m ./100Kchr10_fa2ms_04.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n chr10.txt
echo
echo --------------------------------------------------------------------------------------------------
echo PARAMETERS FOR TFASTA INPUT -f tfa: -w [-z -Y -W -E -O  -T -s] SLIDING WINDOW ANALYSIS OF EMPIRICAL DATA
echo --------------------------------------------------------------------------------------------------
echo 
echo Example tfa.01.txt
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 1 42   -T ../Results/mstatspop_100chr10.tfa.01.txt -w 100000 -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 1 42   -T ../Results/mstatspop_100chr10.tfa.01.txt -w 100000 -n chr10.txt
echo 
echo Example tfa.01b.txt: Same result but in a single line
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 1 42   -T ../Results/mstatspop_100chr10.tfa.01b.txt -w 100000 -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 1 42   -T ../Results/mstatspop_100chr10.tfa.01b.txt -w 100000 -n chr10.txt
echo 
echo Example tfa.01c.txt: Sliding window of 1000 in lines 
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 1 42   -T ../Results/mstatspop_100chr10.tfa.01c.txt -w 1000 -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 1 42   -T ../Results/mstatspop_100chr10.tfa.01c.txt -w 1000 -n chr10.txt
echo 
echo Example tfa.02.txt: Two populations
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.02.txt -w 100000 -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.02.txt -w 100000 -n chr10.txt
echo 
echo Example tfa.03.txt: Three pops + outgroup
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 10 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.tfa.03.txt -G 1 -w 100000 -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 10 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.tfa.03.txt -G 1 -w 100000 -n chr10.txt
echo
echo Example tfa.04.txt: Two pops + outg + missing
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.04.txt -G 1 -u 1 -w 100000 -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.04.txt -G 1 -u 1 -w 100000 -n chr10.txt
echo
echo Example tfa.06.txt: Two pops + outg + missing + window 10000
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.06.txt -G 1 -u 1 -w 10000 -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.06.txt -G 1 -u 1 -w 10000 -n chr10.txt
echo 
echo Example tfa.07.txt: Two pops + outg + missing + window 10000 + slide 20000
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.07.txt -G 1 -u 1 -w 10000 -z 20000 -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.07.txt -G 1 -u 1 -w 10000 -z 20000 -n chr10.txt
echo 
echo Example tfa.08.txt: Two pops + outg + missing + coordinates. Equal than the previous example
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.08.txt -G 1 -u 1 -W ./coord_100Kb.txt -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.08.txt -G 1 -u 1 -W ./coord_100Kb.txt -n chr10.txt
echo 
echo Example tfa.09.txt: Two pops + outg + missing + window 10000 + weight for nonsyn
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.09.txt -G 1 -u 1 -w 10000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.09.txt -G 1 -u 1 -w 10000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10.txt
echo 
echo Example tfa.10.txt: Two pops + outg + missing + window 1000 + weight for nonsyn + Effective_positions
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.10.txt -G 1 -u 1 -w 1000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.10.txt -G 1 -u 1 -w 1000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10.txt
echo 
echo Example tfa.10b.txt: Two pops + outg + missing + window 1000 + weight for nonsyn + Effective_positions + single line results
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.10b.txt -G 1 -u 1 -w 1000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.10b.txt -G 1 -u 1 -w 1000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10.txt
echo 
echo Example tfa.10c.txt: Two pops + outg + missing + window 1000 + 500 displacement + weight for nonsyn + Effective_positions + single line results
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.10c.txt -G 1 -u 1 -w 1000 -z 500 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.10c.txt -G 1 -u 1 -w 1000 -z 500 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10.txt
echo 
echo Example tfa.10d.txt: Two pops + outg + missing + window 1000 + 2000 displacement + weight for nonsyn + Effective_positions + single line results
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.10d.txt -G 1 -u 1 -w 1000 -z 2000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.10d.txt -G 1 -u 1 -w 1000 -z 2000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10.txt
echo 
echo Example tfa.11.txt: Two pops + outg + missing + weight for nonsyn + Effective_positions
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.11.txt -G 1 -u 1 -w 100000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.11.txt -G 1 -u 1 -w 100000 -E ./100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10.txt
echo 
echo Example tfa.12.txt: Three pops + outg + missing + reordering the two pops in inverse order
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.tfa.12.txt -G 1 -u 1 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41 -w 100000 -n chr10.txt
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.tfa.12.txt -G 1 -u 1 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41 -w 100000 -n chr10.txt
echo 
echo Example tfa.13.txt: Three pops + outg + missing + weight for nonsyn
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.tfa.13.txt -G 1 -u 1  -w 100000 -E ./100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10.txt
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.tfa.13.txt -G 1 -u 1  -w 100000 -E ./100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10.txt
echo 
echo Example tfa.14.txt: Three pops + outg + missing + reordering the two pops in inverse order + weight for nonsyn
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.tfa.14.txt -G 1 -u 1  -w 100000 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41  -E ./100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10.txt
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.tfa.14.txt -G 1 -u 1  -w 100000 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41  -E ./100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10.txt
echo 
echo Example tfa.14b.txt: Three pops + outg + missing + reordering the two pops in inverse order + weight for nonsyn + Effective positions
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.tfa.14b.txt -G 1 -u 1  -w 100000 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41  -E ./100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10.txt
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.tfa.14b.txt -G 1 -u 1  -w 100000 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41  -E ./100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10.txt
echo 
echo Example tfa.15.txt: Two pops + outg + missing + window 10000 + slide 20000 + permutations: like example 07
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.15.txt -G 1 -u 1 -w 10000 -z 20000 -t 1000 -s 1684 -n chr10.txt
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.15.txt -G 1 -u 1 -w 10000 -z 20000 -t 1000 -s 1684 -n chr10.txt
echo
echo Example tfa.16.txt: Two pops + outg + missing + window 100 + slide 100 + permutations:
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.16.txt -G 1 -u 1 -w 100 -z 100 -t 1000 -s 1684 -n chr10.txt
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 2 40 2 -T ../Results/mstatspop_100chr10.tfa.16.txt -G 1 -u 1 -w 100 -z 100 -t 1000 -s 1684 -n chr10.txt
echo 
echo Example tfa.17.txt: Three pops + outg + missing + window 100 + slide 100 + permutations:
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.tfa.17.txt -G 1 -u 1 -w 100 -z 100 -t 1000 -s 1684 -n chr10.txt
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 3 20 20 2 -T ../Results/mstatspop_100chr10.tfa.17.txt -G 1 -u 1 -w 100 -z 100 -t 1000 -s 1684 -n chr10.txt
echo
echo Example tfa.18.txt: Three pops LESS SAMPLES + outg + missing + window 100 + slide 100 :
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 3 10 10 2 -T ../Results/mstatspop_100chr10.tfa.18.txt -G 1 -u 1 -w 100 -z 100 -n chr10.txt
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 3 10 10 2 -T ../Results/mstatspop_100chr10.tfa.18.txt -G 1 -u 1 -w 100 -z 100 -n chr10.txt
echo 
echo Example tfa.19.txt: Three pops LESS SAMPLES + reordering + outg + window 100 + slide 100 :
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 3 10 10 2 -T ../Results/mstatspop_100chr10.tfa.19a.txt -G 1 -u 0 -w 100 -z 100 -O 22 20 21 22 23 24 25 26 27 28 29 0 1 2 3 4 5 6 7 8 9 40 41 -n chr10.txt
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 3 10 10 2 -T ../Results/mstatspop_100chr10.tfa.19a.txt -G 1 -u 0 -w 100 -z 100 -O 22 20 21 22 23 24 25 26 27 28 29 0 1 2 3 4 5 6 7 8 9 40 41 -n chr10.txt
echo
echo Example tfa.19.txt: Three pops LESS SAMPLES + reordering + outg + missing + window 100 + slide 100 :
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 3 10 10 2 -T ../Results/mstatspop_100chr10.tfa.19.txt -G 1 -u 1 -w 100 -z 100 -O 22 20 21 22 23 24 25 26 27 28 29 0 1 2 3 4 5 6 7 8 9 40 41 -n chr10.txt
 ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -N 3 10 10 2 -T ../Results/mstatspop_100chr10.tfa.19.txt -G 1 -u 1 -w 100 -z 100 -O 22 20 21 22 23 24 25 26 27 28 29 0 1 2 3 4 5 6 7 8 9 40 41 -n chr10.txt
echo
echo Example tfa.20.txt: rSFS in a single line
echo ../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -R 1 -N 20 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1   -T ../Results/mstatspop_100chr10.tfa.20.txt -w 100000 -n chr10.txt
../bin/mstatspop -f tfa -i ./100Kchr10.tfa.gz -o 1 -R 1 -N 20 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1   -T ../Results/mstatspop_100chr10.tfa.20.txt -w 100000 -n chr10.txt
echo
echo --------------------------------------------------------------------------------------------------
echo Examples using MULTICHROMOSOME TFA FILE
echo --------------------------------------------------------------------------------------------------
echo
echo Three different chromosme analyses
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 1 42   -T ../Results/mstatspop_100allchr.tfa.01.txt   -w 100000 -n chr10-12-14.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 1 42   -T ../Results/mstatspop_100allchr.tfa.01.txt   -w 100000 -n chr10-12-14.txt
echo
echo Same as before but single line output
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 1 42   -T ../Results/mstatspop_100allchr.tfa.01b.txt   -w 100000 -n chr10-12-14.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 1 42   -T ../Results/mstatspop_100allchr.tfa.01b.txt   -w 100000 -n chr10-12-14.txt
echo
echo Same as before but using 10000bp sliding window
 ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 1 42   -T ../Results/mstatspop_100allchr.tfa.01c.txt -w 10000 -n chr10-12-14.txt
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 1 42   -T ../Results/mstatspop_100allchr.tfa.01c.txt -w 10000 -n chr10-12-14.txt
echo
echo Same but using two populations
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100allchr.tfa.02.txt   -w 100000 -n chr10-12-14.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100allchr.tfa.02.txt   -w 100000 -n chr10-12-14.txt
echo
echo Same but using two populations, last is outgroup
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -G 1 -T ../Results/mstatspop_100allchr.tfa.03.txt   -w 100000 -n chr10-12-14.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -G 1 -T ../Results/mstatspop_100allchr.tfa.03.txt   -w 100000 -n chr10-12-14.txt
echo
echo Same but including missing data
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -G 1 -T ../Results/mstatspop_100allchr.tfa.04.txt  -u 1 -w 100000 -n chr10-12-14.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -G 1 -T ../Results/mstatspop_100allchr.tfa.04.txt  -u 1 -w 100000 -n chr10-12-14.txt
echo
echo Same but only including chr10 and chr12
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -G 1 -T ../Results/mstatspop_100allchr.tfa.05.txt  -u 1 -w 100000 -n chr10-12.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -G 1 -T ../Results/mstatspop_100allchr.tfa.05.txt  -u 1 -w 100000 -n chr10-12.txt
echo
echo Same but in single line output
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 2 40 2 -G 1 -T ../Results/mstatspop_100allchr.tfa.06.txt  -u 1 -w 10000 -n chr10-12-14.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 2 40 2 -G 1 -T ../Results/mstatspop_100allchr.tfa.06.txt  -u 1 -w 10000 -n chr10-12-14.txt
echo
echo Same but in single line output, no missing and three pops
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 3 20 20 2 -G 1 -T ../Results/mstatspop_100allchr.tfa.06b.txt  -u 0 -w 10000 -n chr10-12-14.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 3 20 20 2 -G 1 -T ../Results/mstatspop_100allchr.tfa.06b.txt  -u 0 -w 10000 -n chr10-12-14.txt
echo
echo Same but in single line output, no missing and three pops and NO outgroup
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 3 20 20 2 -G 0 -T ../Results/mstatspop_100allchr.tfa.06c.txt  -u 0 -w 10000 -n chr10-12-14.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 3 20 20 2 -G 0 -T ../Results/mstatspop_100allchr.tfa.06c.txt  -u 0 -w 10000 -n chr10-12-14.txt
echo
echo Same but using displacement of 20000bp
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 2 40 2 -G 1 -T ../Results/mstatspop_100allchr.tfa.07.txt  -u 1 -z 20000 -w 10000 -n chr10-12-14.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 2 40 2 -G 1 -T ../Results/mstatspop_100allchr.tfa.07.txt  -u 1 -z 20000 -w 10000 -n chr10-12-14.txt
echo
echo Same but using coordinates file instead sliding window
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 2 40 2 -G 1 -T ../Results/mstatspop_100allchr.tfa.08c1.txt  -u 1 -W ./coord_100Kb_allchr.txt -n chr10-12-14.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 2 40 2 -G 1 -T ../Results/mstatspop_100allchr.tfa.08c1.txt  -u 1 -W ./coord_100Kb_allchr.txt -n chr10-12-14.txt
echo
echo Use the weight file to calculate nonsuynonymous variability
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100allchr.tfa.09.txt -u 1 -w 10000 -E ./100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10-12-14.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100allchr.tfa.09.txt -u 1 -w 10000 -E ./100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10-12-14.txt
echo
echo same but using outgroup weights
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100allchr.tfa.10.txt -G 1 -u 1 -w 10000 -E ./100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10-12-14.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100allchr.tfa.10.txt -G 1 -u 1 -w 10000 -E ./100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n chr10-12-14.txt
echo
echo same but with 100bp sliding window
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100allchr.tfa.11.txt -G 1 -u 1 -w 100 -E ./100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10-12-14.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 0 -N 2 40 2 -T ../Results/mstatspop_100allchr.tfa.11.txt -G 1 -u 1 -w 100 -E ./100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n chr10-12-14.txt
echo
echo three pops and different sample order
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 3 10 10 2 -T ../Results/mstatspop_100allchr.tfa.12.txt -G 1 -u 1 -w 1000 -z 100 -O 22 20 21 22 23 24 25 26 27 28 29 0 1 2 3 4 5 6 7 8 9 40 41 -n chr10-12-14.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -N 3 10 10 2 -T ../Results/mstatspop_100allchr.tfa.12.txt -G 1 -u 1 -w 1000 -z 100 -O 22 20 21 22 23 24 25 26 27 28 29 0 1 2 3 4 5 6 7 8 9 40 41 -n chr10-12-14.txt
echo
echo rSFS in single line output
echo ../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -R 1 -N 20 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1   -T ../Results/mstatspop_100allchr.tfa.20.txt   -w 100000 -n chr10-12-14.txt
../bin/mstatspop -f tfa -i ./100Kallchr.tfa.gz -o 1 -R 1 -N 20 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1   -T ../Results/mstatspop_100allchr.tfa.20.txt   -w 100000 -n chr10-12-14.txt
echo
echo --------------------------------------------------------------------------------------------------
echo Collect some columns from the mstatspop output file:
echo --------------------------------------------------------------------------------------------------
echo 
echo 'perl ../bin/collect_data_columns.pl -in ./mstatspop_100chr10.tfa.15.txt -fc ./required_columns1.txt > ./mstatspop_100chr10.tfa.15_COLUMNS.txt'
perl ../bin/collect_data_columns.pl -in ../Results/mstatspop_100chr10.tfa.15.txt -fc ./required_columns1.txt > ./mstatspop_100chr10.tfa.15_COLUMNS.txt
echo
echo --------------------------------------------------------------------------------------------------
echo Optimal tests using non-missing data:
echo --------------------------------------------------------------------------------------------------
echo 
echo  ../bin/mstatspop -f fasta -i ./MC1R_PigsOutg_aligned.fas -o 0 -p 1 -u 0 -t 1000 -s 123456 -G 1 -N 3 48 46 1 -A ./MC1R_H1frq.txt -g ./MC1R.gff nonsynonymous Nuclear_Universal -T ../Results/MC1R_PigsOutg_NSyn_Opttest.txt -n MC1R.txt
 ../bin/mstatspop -f fasta -i ./MC1R_PigsOutg_aligned.fas -o 0 -p 1 -u 0 -t 1000 -s 123456 -G 1 -N 3 48 46 1 -A ./MC1R_H1frq.txt -g ./MC1R.gff nonsynonymous Nuclear_Universal -T ../Results/MC1R_PigsOutg_NSyn_Opttest.txt -n MC1R.txt
 echo ../bin/mstatspop -f fasta -i ./MC1R_PigsOutg_aligned.fas -o 0 -p 1 -u 0 -s 123456 -G 1 -N 3 48 46 1 -A ./MC1R_H1frq.txt -S ./MC1R_H0frq.txt -T ./MC1R_PigsOutg_NSyn_Opttest_H0.txt -n MC1R.txt
  ../bin/mstatspop -f fasta -i ./MC1R_PigsOutg_aligned.fas -o 0 -p 1 -u 0 -s 123456 -G 1 -N 3 48 46 1 -A ./MC1R_H1frq.txt -S ./MC1R_H0frq.txt -T ./MC1R_PigsOutg_NSyn_Opttest_H0.txt -n MC1R.txt
echo
echo --------------------------------------------------------------------------------------------------
echo "indexing tfasta files:"
echo --------------------------------------------------------------------------------------------------
echo
echo #Usage: ./tfa_index  [options] <input.tfa|input.tfa.gz|input.tfa.gz|weights.txt|weights.txt.gz>
echo #Options:
echo # --version
echo # --help
echo # --threads <int>
echo # --force      "force to overwrite the file"
echo # --weight     "the input file is a weight file. In this case create and index for it. and convert it to the correct format"
echo # --output FILE "set custom name to the output file, only used when converting from TFAv1 to TFAv2 or compressing the input file"
echo # --min-shift <int> "set minimal interval size for CSI indices to 2^INT recommended value is 14, 0 means build tbi index. Default [0]"
echo
echo "Converting files to tfa.gz v2"
echo # list of files
# Examples/V0.1.0/100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz
echo "Converting Examples/V0.1.0/100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz"
../bin/tfa_index  ./V0.1.0/100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -f -o ../Results/V1.0.0/100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz
echo # Examples/V0.1.0/100Kallchr.tfa.gz
echo "Converting Examples/V0.1.0/100Kallchr.tfa.gz"
../bin/tfa_index  ./V0.1.0/100Kallchr.tfa.gz -f -o ../Results/V1.0.0/100Kallchr.tfa.gz
echo
echo # Examples/V0.1.0/100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz
echo "Converting Examples/V0.1.0/100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz"
../bin/tfa_index  ./V0.1.0/100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -f -o ../Results/V1.0.0/100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz
echo # Examples/V0.1.0/100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz
echo "Converting Examples/V0.1.0/100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz"
../bin/tfa_index  ./V0.1.0/100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -f -o ../Results/V1.0.0/100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz
echo # Examples/V0.1.0/100Kchr10.tfa.gz
echo "Converting Examples/V0.1.0/100Kchr10.tfa.gz"
../bin/tfa_index  ./V0.1.0/100Kchr10.tfa.gz -f -o ../Results/V1.0.0/100Kchr10.tfa.gz
echo
echo --------------------------------------------------------------------------------------------------
echo "Merge tfasta files (same assembly reference but different individuals):"
echo --------------------------------------------------------------------------------------------------
echo
echo #Usage: ./tfa_merge [options]
echo #Options:
echo # -i, --input FILE      Input file name (can be specified multiple times)
echo # -o, --output FILE     Output file name
echo # -f, --force           Force overwrite of output file
echo # -h, --help            This help message
echo #     --version         Print version
echo
