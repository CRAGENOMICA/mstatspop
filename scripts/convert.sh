#!/bin/bash


## echo "Converting files to tfa.gz v2"
echo "Converting files to tfa.gz v2"

# list of files 
# Examples/V0.1.0/100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz
echo "Converting Examples/V0.1.0/100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz"
./build/tfa_index  ./Examples/V0.1.0/100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -o ./Examples/V1.0.0//100Kallchr07.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz
# Examples/V0.1.0/100Kallchr.tfa.gz
echo "Converting Examples/V0.1.0/100Kallchr.tfa.gz"
./build/tfa_index  ./Examples/V0.1.0/100Kallchr.tfa.gz -o ./Examples/V1.0.0//100Kallchr.tfa.gz


# Examples/V0.1.0/100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz
echo "Converting Examples/V0.1.0/100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz"
./build/tfa_index  ./Examples/V0.1.0/100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -o ./Examples/V1.0.0//100Kchr10_fa2tfa_08.tfa.gz_npops2_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz
# Examples/V0.1.0/100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz
echo "Converting Examples/V0.1.0/100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz"
./build/tfa_index  ./Examples/V0.1.0/100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz -o ./Examples/V1.0.0//100Kchr10_tfa2tfa_03c.tfa.gz_npops3_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_outg_ploidy1_WEIGHTS.gz
# Examples/V0.1.0/100Kchr10.tfa.gz
echo "Converting Examples/V0.1.0/100Kchr10.tfa.gz"
./build/tfa_index  ./Examples/V0.1.0/100Kchr10.tfa.gz -o ./Examples/V1.0.0//100Kchr10.tfa.gz
