load ./base.tests.bats
load '../node_modules/bats-support/load'
load '../node_modules/bats-assert/load'
load '../node_modules/bats-file/load'
# bats test_tags=tag:bin
@test "mstatspop ok" {
  run mstatspop -h
  assert_success
  [ "${lines[0]}" = "mstatspop v.0.1beta (20231105)" ]
  # [ "${lines[2]}" = "Usage: ./mstatspop -v input.vcf(.gz) -r reference.fa(.gz) -o outputname -n chromosomes.txt" ]
}


# echo --------------------------------------------------------------------------------------------------
# echo PARAMETERS FOR FASTA INPUT -f fasta: [-p -g -c -O  -T -s] WHOLE REGION ANALYSIS
# echo --------------------------------------------------------------------------------------------------
# echo 
# echo Example fa.01.txt
# echo ../bin/mstatspop -f fasta -i ./100Kchr10.fa -o 0 -N 1 42    -T ../Results/mstatspop_100chr10.fa.01.txt -K 1 -n chr10.txt
# ../bin/mstatspop -f fasta -i ./100Kchr10.fa -o 0 -N 1 42   -T ../Results/mstatspop_100chr10.fa.01.txt -K 1 -n chr10.txt
# echo 
@test "fa.01" {
  run mstatspop -f fasta -i $TEST_FILES_DIR/100Kchr10.fa -o 0 -N 1 42   -T $TEST_OUTPUT/mstatspop_100chr10.fa.01.txt -K 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success
  assert_file_exist $TEST_OUTPUT/mstatspop_100chr10.fa.01.txt

}