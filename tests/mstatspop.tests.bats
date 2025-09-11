load ./base.tests.bats
load '../node_modules/bats-support/load'
load '../node_modules/bats-assert/load'
load '../node_modules/bats-file/load'
# bats test_tags=tag:bin
@test "mstatspop ok" {
  run mstatspop -h
  assert_success
  [ "${lines[2]}" = "#mstatspop " ]
  # [ "${lines[0]}" = "mstatspop v.0.1beta (20231105)" ]
  # [ "${lines[2]}" = "Usage: ./mstatspop -v input.vcf(.gz) -r reference.fa(.gz) -o outputname -n chromosomes.txt" ]
}


# echo --------------------------------------------------------------------------------------------------
# echo PARAMETERS FOR FASTA INPUT -f fasta: [-p -g -c -O  -T -s] WHOLE REGION ANALYSIS
# echo --------------------------------------------------------------------------------------------------

# bats test_tags=tag:mstatspop
@test "fa.01" {
  run mstatspop -f fasta -i $TEST_FILES_DIR/100Kchr10.fa -o 0 -N 1 42   -T $TEST_OUTPUT/mstatspop_100chr10.fa.01.txt -K 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.fa.01.txt
}

# bats test_tags=tag:mstatspop
@test "fa.01b" {
  # ...
  run mstatspop -f fasta -i $TEST_FILES_DIR/100Kchr10.fa -o 1 -N 1 42   -T $TEST_OUTPUT/mstatspop_100chr10.fa.01b.txt -K 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.fa.01b.txt
}

# bats test_tags=tag:mstatspop
@test "fa.02" {
  # ...
  run mstatspop -f fasta -i $TEST_FILES_DIR/100Kchr10.fa -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.fa.02.txt -K 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.fa.02.txt
}

# bats test_tags=tag:mstatspop
@test "fa.03" {
  # ...
  run mstatspop -f fasta -i $TEST_FILES_DIR/100Kchr10.fa -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.fa.03.txt -G 1 -K 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.fa.03.txt
}

# bats test_tags=tag:mstatspop
@test "fa.04" {
  # ...
  run mstatspop -f fasta -i $TEST_FILES_DIR/100Kchr10.fa -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.fa.04.txt -G 1 -u 1 -K 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.fa.04.txt
}

# bats test_tags=tag:mstatspop
@test "fa.05" {
  # ...
  run mstatspop -f fasta -i $TEST_FILES_DIR/100Kchr10.fa -o 0 -N 3 20 20 2 -T $TEST_OUTPUT/mstatspop_100chr10.fa.05.txt -G 1 -u 1 -K 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.fa.05.txt
}

# bats test_tags=tag:mstatspop
@test "fa.06" {
  # ...
  run mstatspop -f fasta -i $TEST_FILES_DIR/100Kchr10.fa -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.fa.06.txt -G 1 -u 1 -p 2 -K 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.fa.06.txt
}

# bats test_tags=tag:mstatspop
@test "fa.07" {
  # ...
  run mstatspop -f fasta -i $TEST_FILES_DIR/100Kchr10.fa -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.fa.07.txt -G 1 -u 1 -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -K 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.fa.07.txt
}

# bats test_tags=tag:mstatspop
@test "fa.08" {
  # ...
  run mstatspop -f fasta -i $TEST_FILES_DIR/100Kchr10.fa -o 0 -N 3 20 20 2 -T $TEST_OUTPUT/mstatspop_100chr10.fa.08.txt -G 1 -u 0 -t 1000 -s 1684 -K 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.fa.08.txt
}

# bats test_tags=tag:mstatspop
@test "fa.09" {
  # ...
  run mstatspop -f fasta -i $TEST_FILES_DIR/100Kchr10.fa -o 0 -N 3 20 20 2 -T $TEST_OUTPUT/mstatspop_100chr10.fa.09.txt -G 1 -u 1 -K 1 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.fa.09.txt
}

# bats test_tags=tag:mstatspop
@test "fa.10" {
  # ...
  run mstatspop -f fasta -i $TEST_FILES_DIR/100Kchr10.fa -o 0 -N 3 20 20 2 -T $TEST_OUTPUT/mstatspop_100chr10.fa.10.txt -G 1 -u 1  -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -K 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.fa.10.txt
}

# bats test_tags=tag:mstatspop
@test "fa.11" {
  # ...
  run mstatspop -f fasta -i $TEST_FILES_DIR/100Kchr10.fa -o 0 -N 3 20 20 2 -T $TEST_OUTPUT/mstatspop_100chr10.fa.11.txt -G 1 -u 1 -K 0 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41  -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.fa.11.txt
}

# bats test_tags=tag:mstatspop
@test "fa.12" {
  # ...
  run mstatspop -f fasta -i $TEST_FILES_DIR/100Kchr10.fa -o 0 -R 1 -N 20 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -T $TEST_OUTPUT/mstatspop_100chr10.fa.12.txt -G 0 -u 1  -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -K 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.fa.12.txt
}

# echo --------------------------------------------------------------------------------------------------
# echo PARAMETERS FOR MS INPUT -f ms: [-r -m -F -v -q] SIMULATION ANALYSIS OF A SINGLE REGION
# echo --------------------------------------------------------------------------------------------------

# bats test_tags=tag:mstatspop
@test "ms.00" {
  # ...
  run mstatspop -f ms -i $TEST_FILES_DIR/100Kchr10_tfa2ms_02.ms.txt -o 0 -N 1 42 -T $TEST_OUTPUT/mstatspop_100chr10.ms.00.txt -l 100000 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.ms.00.txt
}

# bats test_tags=tag:mstatspop
@test "ms.01" {
  # ...
  run mstatspop -f ms -i $TEST_FILES_DIR/100Kchr10_tfa2ms_02.ms.txt -o 0 -N 1 42 -T $TEST_OUTPUT/mstatspop_100chr10.ms.01.txt -l 100000 -m $TEST_FILES_DIR/100Kchr10_npop1_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.ms.01.txt
}

# bats test_tags=tag:mstatspop
@test "ms.01b" {
  # ...
  run mstatspop -f ms -i $TEST_FILES_DIR/100Kchr10_tfa2ms_02.ms.txt -o 1 -N 1 42 -T $TEST_OUTPUT/mstatspop_100chr10.ms.01b.txt -l 100000 -m $TEST_FILES_DIR/100Kchr10_npop1_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.ms.01b.txt
}

# bats test_tags=tag:mstatspop
@test "ms.02" {
  # ...
  run mstatspop -f ms -i $TEST_FILES_DIR/100Kchr10_tfa2ms_02.ms.txt -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.ms.02.txt -l 100000 -m $TEST_FILES_DIR/100Kchr10_npop2_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.ms.02.txt
}

# bats test_tags=tag:mstatspop
@test "ms.03" {
  # ...
  run mstatspop -f ms -i $TEST_FILES_DIR/100Kchr10_tfa2ms_03.ms.txt -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.ms.03.txt -l 100000 -G 1 -m $TEST_FILES_DIR/100Kchr10_npop2_nsam42_ExcludeMissingVariantsmhits_NOoutg_ploidy1_MASK.txt -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.ms.03.txt
}

# bats test_tags=tag:mstatspop
@test "ms.04" {
  # ...
  run mstatspop -f ms -i $TEST_FILES_DIR/100Kchr10_tfa2ms_07.ms.txt -o 1 -N 1 42 -T $TEST_OUTPUT/mstatspop_100chr10.ms.04.txt -l 10000 -r 5 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.ms.04.txt
}

# bats test_tags=tag:mstatspop
@test "ms.05" {
  # ...
  run mstatspop -f ms -i $TEST_FILES_DIR/100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.ms.05.txt -l 10000 -r 5 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.ms.05.txt
}

# bats test_tags=tag:mstatspop
@test "ms.06" {
  # ...
  run mstatspop -f ms -i $TEST_FILES_DIR/100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.ms.06.txt -l 10000 -r 5 -F 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.ms.06.txt
}

# bats test_tags=tag:mstatspop
@test "ms.07" {
  # ...
  run mstatspop -f ms -i $TEST_FILES_DIR/100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.ms.07.txt -l 10000 -r 5 -F 1 -q 0.01 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.ms.07.txt
}

# bats test_tags=tag:mstatspop
@test "ms.08" {
  # ...
  run mstatspop -f ms -i $TEST_FILES_DIR/100Kchr10_tfa2ms_06.ms.txt -o 1 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.ms.08.txt -l 10000 -r 5 -G 1 -v 2.0 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.ms.08.txt
}


# echo --------------------------------------------------------------------------------------------------
# echo PARAMETERS FOR TFASTA INPUT -f tfa: -w [-z -Y -W -E -O  -T -s] SLIDING WINDOW ANALYSIS OF EMPIRICAL DATA
# echo --------------------------------------------------------------------------------------------------

# bats test_tags=tag:mstatspop
@test "tfa.01" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 0 -N 1 42   -T $TEST_OUTPUT/mstatspop_100chr10.tfa.01.txt -w 100000 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.01.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.01b" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 1 -N 1 42   -T $TEST_OUTPUT/mstatspop_100chr10.tfa.01b.txt -w 100000 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.01b.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.01c" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 1 -N 1 42   -T $TEST_OUTPUT/mstatspop_100chr10.tfa.01c.txt -w 1000 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.01c.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.02" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.02.txt -w 100000 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.02.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.03" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 10 -N 3 20 20 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.03.txt -G 1 -w 100000 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.03.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.04" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.04.txt -G 1 -u 1 -w 100000 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.04.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.06" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 1 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.06.txt -G 1 -u 1 -w 10000 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.06.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.07" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 1 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.07.txt -G 1 -u 1 -w 10000 -z 20000 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.07.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.08" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.08.txt -G 1 -u 1 -W $TEST_FILES_DIR/coord_100Kb.txt -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.08.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.09" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.09.txt -G 1 -u 1 -w 10000 -E $TEST_FILES_DIR/100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.09.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.10" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.10.txt -G 1 -u 1 -w 1000 -E $TEST_FILES_DIR/100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.10.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.10b" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 1 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.10b.txt -G 1 -u 1 -w 1000 -E $TEST_FILES_DIR/100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.10b.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.10c" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 1 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.10c.txt -G 1 -u 1 -w 1000 -z 500 -E $TEST_FILES_DIR/100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.10c.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.10d" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 1 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.10d.txt -G 1 -u 1 -w 1000 -z 2000 -E $TEST_FILES_DIR/100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.10d.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.11" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.11.txt -G 1 -u 1 -w 100000 -E $TEST_FILES_DIR/100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.11.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.12" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.12.txt -G 1 -u 1 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41 -w 100000 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.12.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.13" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.13.txt -G 1 -u 1  -w 100000 -E $TEST_FILES_DIR/100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.13.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.14" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.14.txt -G 1 -u 1  -w 100000 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41  -E $TEST_FILES_DIR/100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.14.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.14b" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 0 -N 3 20 20 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.14b.txt -G 1 -u 1  -w 100000 -O 42 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 40 41  -E $TEST_FILES_DIR/100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.14b.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.15" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 1 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.15.txt -G 1 -u 1 -w 10000 -z 20000 -t 1000 -s 1684 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.15.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.16" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 1 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.16.txt -G 1 -u 1 -w 100 -z 100 -t 1000 -s 1684 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.16.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.17" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 1 -N 3 20 20 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.17.txt -G 1 -u 1 -w 100 -z 100 -t 1000 -s 1684 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.17.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.18" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 1 -N 3 10 10 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.18.txt -G 1 -u 1 -w 100 -z 100 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.18.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.19" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 1 -N 3 10 10 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.19a.txt -G 1 -u 0 -w 100 -z 100 -O 22 20 21 22 23 24 25 26 27 28 29 0 1 2 3 4 5 6 7 8 9 40 41 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.19a.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.19b" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 1 -N 3 10 10 2 -T $TEST_OUTPUT/mstatspop_100chr10.tfa.19b.txt -G 1 -u 1 -w 100 -z 100 -O 22 20 21 22 23 24 25 26 27 28 29 0 1 2 3 4 5 6 7 8 9 40 41 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.19b.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.20" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kchr10.tfa.gz -o 1 -R 1 -N 20 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1   -T $TEST_OUTPUT/mstatspop_100chr10.tfa.20.txt -w 100000 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.tfa.20.txt
}

# echo --------------------------------------------------------------------------------------------------
# echo Examples using MULTICHROMOSOME TFA FILE
# echo --------------------------------------------------------------------------------------------------

# bats test_tags=tag:mstatspop
@test "tfa.01m" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 0 -N 1 42   -T $TEST_OUTPUT/mstatspop_100allchr.tfa.01m.txt   -w 100000 -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.01m.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.01bm" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 1 -N 1 42   -T $TEST_OUTPUT/mstatspop_100allchr.tfa.01bm.txt   -w 100000 -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.01bm.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.01cm" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 1 -N 1 42   -T $TEST_OUTPUT/mstatspop_100allchr.tfa.01cm.txt -w 10000 -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.01cm.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.02m" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100allchr.tfa.02m.txt   -w 100000 -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.02m.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.03m" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 0 -N 2 40 2 -G 1 -T $TEST_OUTPUT/mstatspop_100allchr.tfa.03m.txt   -w 100000 -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.03m.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.04m" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 0 -N 2 40 2 -G 1 -T $TEST_OUTPUT/mstatspop_100allchr.tfa.04m.txt  -u 1 -w 100000 -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.04m.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.05m" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 0 -N 2 40 2 -G 1 -T $TEST_OUTPUT/mstatspop_100allchr.tfa.05m.txt  -u 1 -w 100000 -n $TEST_FILES_DIR/chr10-12.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.05m.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.06m" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 1 -N 2 40 2 -G 1 -T $TEST_OUTPUT/mstatspop_100allchr.tfa.06m.txt  -u 1 -w 10000 -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.06m.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.06bm" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 1 -N 3 20 20 2 -G 1 -T $TEST_OUTPUT/mstatspop_100allchr.tfa.06bm.txt  -u 0 -w 10000 -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.06bm.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.06cm" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 1 -N 3 20 20 2 -G 0 -T $TEST_OUTPUT/mstatspop_100allchr.tfa.06cm.txt  -u 0 -w 10000 -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.06cm.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.07m" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 1 -N 2 40 2 -G 1 -T $TEST_OUTPUT/mstatspop_100allchr.tfa.07m.txt  -u 1 -z 20000 -w 10000 -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.07m.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.08m" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 1 -N 2 40 2 -G 1 -T $TEST_OUTPUT/mstatspop_100allchr.tfa.08m.txt  -u 1 -W $TEST_FILES_DIR/coord_100Kb_allchr.txt -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.08m.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.09m" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100allchr.tfa.09m.txt -u 1 -w 10000 -E $TEST_FILES_DIR/100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.09m.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.10m" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100allchr.tfa.10m.txt -G 1 -u 1 -w 10000 -E $TEST_FILES_DIR/100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.10m.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.11m" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 0 -N 2 40 2 -T $TEST_OUTPUT/mstatspop_100allchr.tfa.11m.txt -G 1 -u 1 -w 100 -E $TEST_FILES_DIR/100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -Y 0 -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.11m.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.12m" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 1 -N 3 10 10 2 -T $TEST_OUTPUT/mstatspop_100allchr.tfa.12m.txt -G 1 -u 1 -w 1000 -z 100 -O 22 20 21 22 23 24 25 26 27 28 29 0 1 2 3 4 5 6 7 8 9 40 41 -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.12m.txt
}

# bats test_tags=tag:mstatspop
@test "tfa.20m" {
  # ...
  run mstatspop -f tfa -i $TEST_FILES_DIR/100Kallchr.tfa.gz -o 1 -R 1 -N 20 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1   -T $TEST_OUTPUT/mstatspop_100allchr.tfa.20m.txt   -w 100000 -n $TEST_FILES_DIR/chr10-12-14.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100allchr.tfa.20m.txt
}

# echo --------------------------------------------------------------------------------------------------
# echo Optimal tests using non-missing data:
# echo --------------------------------------------------------------------------------------------------

# bats test_tags=tag:mstatspop
@test "opt.1" {
  # ...
  run mstatspop -f fasta -i $TEST_FILES_DIR/MC1R_PigsOutg_aligned.fas -o 0 -p 1 -u 0 -t 1000 -s 123456 -G 1 -N 3 48 46 1 -A $TEST_FILES_DIR/MC1R_H1frq.txt -g $TEST_FILES_DIR/MC1R.gff nonsynonymous Nuclear_Universal -T $TEST_OUTPUT/MC1R_PigsOutg_NSyn_Opttest.txt -n $TEST_FILES_DIR/MC1R.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/MC1R_PigsOutg_NSyn_Opttest.txt 
}

# bats test_tags=tag:mstatspop
@test "opt.2" {
  # ...
  run mstatspop -f fasta -i $TEST_FILES_DIR/MC1R_PigsOutg_aligned.fas -o 0 -p 1 -u 0 -s 123456 -G 1 -N 3 48 46 1 -A $TEST_FILES_DIR/MC1R_H1frq.txt -S $TEST_FILES_DIR/MC1R_H0frq.txt -T $TEST_OUTPUT/MC1R_PigsOutg_NSyn_Opttest_H0.txt -n $TEST_FILES_DIR/MC1R.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/MC1R_PigsOutg_NSyn_Opttest_H0.txt
}

# echo --------------------------------------------------------------------------------------------------
# echo "indexing tfasta files:"
# echo --------------------------------------------------------------------------------------------------
# echo #Usage: ./tfa_index  [options] <input.tfa|input.tfa.gz|input.tfa.gz|weights.txt|weights.txt.gz>
# echo #Options:
# echo # --version
# echo # --help
# echo # --threads <int>
# echo # --force      "force to overwrite the file"
# echo # --weight     "the input file is a weight file. In this case create and index for it. and convert it to the correct format"
# echo # --output FILE "set custom name to the output file, only used when converting from TFAv1 to TFAv2 or compressing the input file"
# echo # --min-shift <int> "set minimal interval size for CSI indices to 2^INT recommended value is 14, 0 means build tbi index. Default [0]"


# bats test_tags=tag:tfa_index
@test "ind.01" {
  # ...
  run tfa_index  $TEST_FILES_DIR/V0.1.0/100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -f -o $TEST_OUTPUT/100Kallchr07.ind.01_WEIGHTS.gz
  assert_success
  assert_file_exist $TEST_OUTPUT/100Kallchr07.ind.01_WEIGHTS.gz.tfa.gz
}

# bats test_tags=tag:tfa_index
@test "ind.02" {
  # ...
  run tfa_index  $TEST_FILES_DIR/V0.1.0/100Kallchr.tfa.gz -f -o $TEST_OUTPUT/100Kallchr_ind.02.tfa.gz
  assert_success
  assert_file_exist $TEST_OUTPUT/100Kallchr_ind.02.tfa.gz
}

# bats test_tags=tag:tfa_index
@test "ind.03" {
  # ...
  run tfa_index  $TEST_FILES_DIR/V0.1.0/100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -f -o $TEST_OUTPUT/100Kchr10_fa2tfa_08.ind.03_WEIGHTS.gz
  assert_success
  assert_file_exist $TEST_OUTPUT/100Kchr10_fa2tfa_08.ind.03_WEIGHTS.gz.tfa.gz
}

# bats test_tags=tag:tfa_index
@test "ind.04" {
  # ...
  run tfa_index  $TEST_FILES_DIR/V0.1.0/100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -f -o $TEST_OUTPUT/100Kchr10_tfa2tfa_03c.ind.04_WEIGHTS.gz
  assert_success
  assert_file_exist $TEST_OUTPUT/100Kchr10_tfa2tfa_03c.ind.04_WEIGHTS.gz.tfa.gz
}

# bats test_tags=tag:tfa_index
@test "ind.05" {
  # ...
  run tfa_index  $TEST_FILES_DIR/V0.1.0/100Kchr10.tfa.gz -f -o $TEST_OUTPUT/100Kchr10.ind.05.tfa.gz
  assert_success
  assert_file_exist $TEST_OUTPUT/100Kchr10.ind.05.tfa.gz
}

# @test 'refute_output' {
#   run echo 'want'
#   refute_output 'want'
#   echo 'want' | assert_output
# }