#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <htslib/tbx.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h> // Include the header file that defines tbx_intv_t.
#include <getopt.h>
#include <assert.h>
#include "log.h"
#include "tfasta.h"
#include "files_util.h"
#include <unistd.h>

// Move the function prototype here
static void print_usage(const char *program_name, int exit_code);

#define MAX_INPUT_FILES 2

/**
 * @brief Structure representing the arguments for tfa_merge.
 * 
 * This structure holds the arguments required for the tfa_merge function.
 * It includes the input file name and the output file name.
 */
typedef struct
{
  char **input_fnames;  // Array of input filenames
  int input_count;      // Number of input files
  char *output_fname;
  int force;           // Add force flag
} tfa_merge_args_t;

int main(int argc, char *argv[])
{

  static const struct option loptions[] =
      {
          {"version", no_argument, NULL, 1},
          {"help", no_argument, NULL, 2},
          {"input", required_argument, NULL, 'i'},
          {"output", required_argument, NULL, 'o'},
          {"force", no_argument, NULL, 'f'},  // Add force option
      };
  int c, ret = 0;

  char *reheader = NULL;
  tfa_merge_args_t args;
  memset(&args, 0, sizeof(tfa_merge_args_t));
  args.input_fnames = malloc(MAX_INPUT_FILES * sizeof(char*));
  args.input_count = 0;

  while ((c = getopt_long(argc, argv, "i:o:f@:", loptions, NULL)) >= 0)
  {
    switch (c)
    {
    case 'h':
    case '?':
      print_usage(argv[0], EXIT_SUCCESS);
      return 0;
    case 2:
      print_usage(argv[0], EXIT_SUCCESS);
      return 0;
    case 1:
      print_usage(argv[0], EXIT_SUCCESS);
      return 0;
    case 'i':
      if (args.input_count >= MAX_INPUT_FILES) {
        log_error("Currently only merging of two files is supported. You provided more than 2 input files.");
        return EXIT_FAILURE;
      }
      args.input_fnames[args.input_count++] = optarg;
      break;
    case 'o':
      args.output_fname = optarg;
      break;
    case 'f':
      args.force = 1;
      break;
    default:
      fprintf(stderr, "Unknown option: %s\n", optarg);
      return 1;
    }
  }

  // Modify the input check
  if (args.input_count == 0)
  {
    log_error("At least one input file is required");
    print_usage(argv[0], EXIT_FAILURE);
  }
  else if (args.input_count != 2)
  {
    log_error("Exactly two input files are required for merging");
    print_usage(argv[0], EXIT_FAILURE);
  }

   #ifdef DEBUG
	log_set_level(LOG_TRACE);
	#else
	log_set_level(LOG_INFO);
	#endif
	// set program name as tfa_merge
  const char *program_name = "tfa_merge";
	log_start(program_name, argc, argv);

  // Echo filenames
  for (int i = 0; i < args.input_count; i++) {
    log_info("Input file %d: %s", i + 1, args.input_fnames[i]);
  }

  if (args.output_fname != NULL)
  {
    if (ends_with(args.output_fname, ".tfa") || ends_with(args.output_fname, ".tfasta"))
    {
      args.output_fname = strcat(args.output_fname, ".gz");
    }
    else
    {
      args.output_fname = strcat(args.output_fname, ".tfa.gz");
    }
  }
  else
  {
   // use merged_tfa.gz as output file name
   args.output_fname = "merged_tfa.gz";

  }

  // Modify the output file check to respect force flag
  if (file_exists(args.output_fname) && !args.force)
  {
    log_error("Output file %s already exists. Use --force to overwrite", args.output_fname);
    return EXIT_FAILURE;
  }


  log_info("Merging files into %s", args.output_fname);

  tfasta_file *tfasta1;
  tfasta_file *tfasta2;

  log_info("Validating tfasta file %s", args.input_fnames[0]);
  // validate tfasta file exists and is readable and in TFAv2.0 format
  if (access(args.input_fnames[0], R_OK) != 0)
  {
    // fprintf(file_logerr, "\nError: the file %s does not exist or is not readable.\n", file_in);
    log_error("Error: the file %s does not exist or is not readable.", args.input_fnames[0]);
    return EXIT_FAILURE;
  }

  if ((tfasta1 = (tfasta_file *)malloc(sizeof(tfasta_file))) == NULL)
  {

    log_fatal("Error: can not allocate memory for tfasta_file structure.");
    return EXIT_FAILURE;
  }

  // validate weights file is in TFAv2.0 format
  // memset the structure to 0
  memset(tfasta1, 0, sizeof(tfasta_file));

  int ret_status = init_tfasta_file(tfasta1, args.input_fnames[0]);
  if (ret_status == TFA_ERROR)
  {
    log_error("Can not initialize tfasta file %s", args.input_fnames[0]);
    return EXIT_FAILURE;
  }
  if (ret_status == TFA_INVALID_FORMAT)
  {
    log_error("The file %s is not in TFAv2.0 format", args.input_fnames[0]);
    return EXIT_FAILURE;
  }
  if (ret_status == TFA_INVALID_INDEX)
  {
    log_error("The index file for %s is invalid", args.input_fnames[0]);
    return EXIT_FAILURE;
  }

  // read the second file
  log_info("Validating tfasta file %s", args.input_fnames[1]);
  // validate tfasta file exists and is readable and in TFAv2.0 format
  if (access(args.input_fnames[1], R_OK) != 0)
  {
    log_error("Error: the file %s does not exist or is not readable.", args.input_fnames[1]);
    return EXIT_FAILURE;
  }

  if ((tfasta2 = (tfasta_file *)malloc(sizeof(tfasta_file))) == NULL)
  {

    log_fatal("Error: can not allocate memory for tfasta_file structure.");
    return EXIT_FAILURE;
  }

  // validate weights file is in TFAv2.0 format
  // memset the structure to 0
  memset(tfasta2, 0, sizeof(tfasta_file));

  ret_status = init_tfasta_file(tfasta2, args.input_fnames[1]);
  if (ret_status == TFA_ERROR)
  {
    log_error("Can not initialize tfasta file %s", args.input_fnames[1]);
    return EXIT_FAILURE;
  }
  if (ret_status == TFA_INVALID_FORMAT)
  {
    log_error("The file %s is not in TFAv2.0 format", args.input_fnames[1]);
    return EXIT_FAILURE;
  }
  if (ret_status == TFA_INVALID_INDEX)
  {
    log_error("The index file for %s is invalid", args.input_fnames[1]);
    return EXIT_FAILURE;
  }


  // so far so good, now we need to merge the two files
  // before merging we need to check if the two files have the same sequences
  log_info("Checking if both files have the same sequences...");
  
  if (tfasta1->nseq != tfasta2->nseq) {
    log_error("Files have different number of sequences (%d vs %d)", 
              tfasta1->nseq, tfasta2->nseq);
    return EXIT_FAILURE;
  }

  // Create a lookup array for matching sequence indices
  int *seq_index_map = (int *)calloc(tfasta1->nseq, sizeof(int));
  if (!seq_index_map) {
    log_error("Failed to allocate memory for sequence mapping");
    return EXIT_FAILURE;
  }

  // For each sequence in file1, find its matching sequence in file2
  for (int i = 0; i < tfasta1->nseq; i++) {
    seq_index_map[i] = -1;  // Initialize to "not found"
    for (int j = 0; j < tfasta2->nseq; j++) {
      if (strcmp(tfasta1->seq_names[i], tfasta2->seq_names[j]) == 0) {
        seq_index_map[i] = j;
        // Check record counts match
        if (tfasta1->n_records[i] != tfasta2->n_records[j]) {
          log_error("Number of records mismatch for sequence %s: %lu vs %lu",
                    tfasta1->seq_names[i], 
                    tfasta1->n_records[i], 
                    tfasta2->n_records[j]);
          free(seq_index_map);
          return EXIT_FAILURE;
        }
        break;
      }
    }
    if (seq_index_map[i] == -1) {
      log_error("Sequence %s from first file not found in second file", 
                tfasta1->seq_names[i]);
      free(seq_index_map);
      return EXIT_FAILURE;
    }
  }

  log_info("Sequence validation passed - both files have matching sequences and records");


  // start writing the output file
  BGZF *output_fp;
  output_fp = bgzf_open(args.output_fname, "w");
  if (!output_fp) {
    log_error("Failed to open output file %s", args.output_fname);
    return EXIT_FAILURE;
  }

  // write the TFAv2.0 header to the output file
  if (bgzf_write(output_fp, v2_header, strlen(v2_header)) < 0) {
    log_error("Failed to write TFAv2.0 header to output file %s", args.output_fname);
    return EXIT_FAILURE;
  }
  // write command line to the output file
  kstring_t str = {0, 0, NULL};
  ksprintf(&str, "#tfa_merge ");
  for (int i = 1; i < argc; i++) {
    kputsn(argv[i], strlen(argv[i]), &str);
    kputc(' ', &str);
  }
  // add a newline at the end
  kputc('\n', &str);
  if (bgzf_write(output_fp, str.s, str.l) < 0) {
    log_error("Failed to write command line to output file %s", args.output_fname);
    return EXIT_FAILURE;
  }

  // reset the kstring_t structure
  str.l = 0;
  

  // write the merged sample names to the output file 
  ksprintf(&str, "#NAMES: ");
  // tfasta1->names is an array of strings separated by spaces  
  for (int i = 0; i < tfasta1->n_sam; i++) {
    kputsn(tfasta1->names[i], strlen(tfasta1->names[i]), &str);
    kputc(' ', &str);
  }
  // tfasta2->names is an array of strings separated by spaces
  for (int i = 0; i < tfasta2->n_sam; i++) {
    kputsn(tfasta2->names[i], strlen(tfasta2->names[i]), &str);
    // do not add a space after the last sample name
    if (i < tfasta2->n_sam - 1) {
      kputc(' ', &str);
    }
  }
  kputc('\n', &str);
  if (bgzf_write(output_fp, str.s, str.l) < 0) {
    log_error("Failed to write sample names to output file %s", args.output_fname);
    return EXIT_FAILURE;
  }

  // reset the kstring_t structure
  str.l = 0;
  // write #CHR	POSITION	GENOTYPES
  ksprintf(&str, "#CHR\tPOSITION\tGENOTYPES\n");
  if (bgzf_write(output_fp, str.s, str.l) < 0) {
    log_error("Failed to write header to output file %s", args.output_fname);
    return EXIT_FAILURE;
  }
  // reset the kstring_t structure
  str.l = 0;
  // free the kstring_t structure
  free(str.s);

  // Use seq_index_map for merging here...
  // let's start by merging the DNA sequences
  // loop through all sequences and merge the DNA sequences
  for (int i = 0; i < tfasta1->nseq; i++) {
    // log the sequence name
    long int position = 1;
    log_info("Merging sequence %s", tfasta1->seq_names[i]);
    char *chr_name = tfasta1->seq_names[i];
    hts_itr_t *iter1 = tbx_itr_querys(tfasta1->tbx, chr_name);
    hts_itr_t *iter2 = tbx_itr_querys(tfasta2->tbx, chr_name);
    if (iter1 == NULL) {
      log_error("Failed to parse region: %s", chr_name);
      return EXIT_FAILURE;
    }
    if (iter2 == NULL) {
      log_error("Failed to parse region: %s", chr_name);
      return EXIT_FAILURE;
    }
    kstring_t str1 = {0, 0, NULL};
    kstring_t str2 = {0, 0, NULL};
    // str output is a string with the following format:
    //  seq_name1\tpos1\tdna_seq1dna_seq2
    kstring_t output_str = {0, 0, NULL};

    const char *delim = ":\t\n";
   
    // keep track DNA_matr2 size
    int count = 0;
    while (tbx_itr_next(tfasta1->fp, tfasta1->tbx, iter1, &str1) >= 0 &&
           tbx_itr_next(tfasta2->fp, tfasta2->tbx, iter2, &str2) >= 0) {
      // get the DNA sequence from tfasta1
      char *seq_name1 = strtok(str1.s, delim);
      char *str_pos1 = strtok(NULL, delim);
      long int pos1 = atol(str_pos1);
      char *dna_seq1 = strtok(NULL, delim);

      // get the DNA sequence from tfasta2
      char *seq_name2 = strtok(str2.s, delim);
      char *str_pos2 = strtok(NULL, delim);
      long int pos2 = atol(str_pos2);
      char *dna_seq2 = strtok(NULL, delim);

       // Verify positions match
       if (pos1 != pos2) {
        log_error("Positions do not match for sequence %s: %ld vs %ld", 
                  tfasta1->seq_names[i], pos1, pos2);
        return EXIT_FAILURE;
       }

      ksprintf(&output_str, "%s\t%ld\t%s%s\n", 
               seq_name1, pos1, dna_seq1, dna_seq2);
      if (bgzf_write(output_fp, output_str.s, output_str.l) < 0) {
        log_error("Failed to write merged sequence to output file %s", args.output_fname);
        return EXIT_FAILURE;
      }
      // reset the kstring_t structure
      output_str.l = 0;
      // log the sequence name and position and DNA sequence
      // log_info("Sequence1 %s at position %ld: %s", seq_name1, pos1, dna_seq1);
      //log_info("Sequence2 %s at position %ld: %s", seq_name2, pos2, dna_seq2);
    }
  }

  

  // close the output file
  bgzf_close(output_fp);


  log_info("Creating index for the output file: %s", args.output_fname);
  // create index for the output file
  create_index(args.output_fname, 2, 0, tfasta_conf);
  log_info("Index was created successfully");
  // Free the mapping when done
  free(seq_index_map);

  // end of merging
  // close the two files
  close_tfasta_file(tfasta1);
  close_tfasta_file(tfasta2);


  // Don't forget to free memory before exit
  free(args.input_fnames);

  log_info("Program finished successfully");

  return EXIT_SUCCESS;
  
}

static void print_usage(const char *program_name, int exit_code)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s [options]\n", program_name);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -i, --input FILE      Input file name (can be specified multiple times)\n");
    fprintf(stderr, "  -o, --output FILE     Output file name\n");
    fprintf(stderr, "  -f, --force           Force overwrite of output file\n");  // Add force description
    fprintf(stderr, "  -h, --help            This help message\n");
    fprintf(stderr, "      --version         Print version\n");
    fprintf(stderr, "\n");
    exit(exit_code);
}
