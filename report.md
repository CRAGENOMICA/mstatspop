# Update Version 1.0.0 

mstatspop now accept tfasta file  in version 2  compressed with htslib (bgzip) and indexed format (.tbi).

### Implementation
mstatspop now uses htslib to read and write Files,the input tfasta file/weight files is compressed (bgzip) and indexed. The tfasta index file is created with htslib tabix format.

### TFasta file format
The tfasta file format has the following format:
```ini
##fileformat=TFAv2.0
# command1 
# command2
#NAMES: >chr1_a >chr1_b .....
#CHR	POSITION	GENOTYPES
chr10	31338	GGGGNNGGGGGGGGGGGGGGGGGGGGGGGGGGNNGGGGGGGG
chr10	31339	GGGGNNGGGGGGGGGGGGGGGGGGGGGGGGGGNNGGGGGGGG
chr10	31340	AAAANNAAAAAAAAAAAAAAAAAAAAAAAAAANNAAAAAAAA
chr10	31341	TTTTNNTTTTTTTTTTTTTTTTTTTTTTTTTTNNTTTTTTTT
chr10	31342	AAAANNAAAAAAAAAAAAAAAAAAAAAAAAAANNAAAAAAAA
chr10	31343	CCCCNNCCCCCCCCCCCCCCCCCCCCCCCCCCNNCCCCCCCC
chr10	31344	TTTTNNTTTTTTTTTTTTTTTTTTTTTTTTTTNNTTTTTTTT
chr10	31345	GGGGNNGGGGGGGGGGGGGGGGGGGGGGGGGGNNGGGGGGGG
```

Header section:
- `##fileformat=TFAv2.0` : The version of the tfasta file format.
- `#command1`: The command used to generate the tfasta file.
- `#command2`: other command used to generate the tfasta file.
- `#NAMES` : The names of the sequences in the tfasta file. The names are separated by a space and start with the "`>`" character.

Data section is a tab-separated file with the following columns:
- `CHR`: Chromosome name
- `POSITION`: Position in the chromosome
- `GENOTYPES`: Genotypes in the position. The genotypes are represented by the IUPAC nucleotide code. Missing data is represented by N.

### Weight file format
The weight file format aslo follow the same fomat with extra columns

