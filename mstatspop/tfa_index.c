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

/**
 * @brief Structure representing the arguments for tfa_index.
 * 
 * This structure holds the arguments required for the tfa_index function.
 * It includes the output file name, the minimum shift value, the number of threads,
 * and a flag indicating whether to force the operation.
 */
typedef struct
{
  char *output_fname; /**< The output file name. */
  int min_shift; /**< The minimum shift value. */
  int threads; /**< The number of threads. */
  int is_force; /**< Flag indicating whether to force the operation. */
  int is_weight; /**< Flag indicating whether the input file is a weight file. */
} tfa_index_args_t;


const tbx_conf_t tfasta_conf = {TBX_GENERIC, 1, 2, 2, '#', 0};


/**
 * Compresses a file.
 *
 * This function takes an input file and compresses it, saving the compressed
 * data to an output file.
 *
 * @param input_filename The path to the input file.
 * @param output_filename The path to the output file.
 * @return Returns 0 if the file was successfully compressed, or a negative
 *         value if an error occurred.
 */
int compress_file(const char *input_filename, const char *output_filename)
{
  char buffer[4096];
  int read_len;
  BGZF *bgzf_output;

  // Open the input file in binary read mode
  // FILE *input_file = fopen(input_filename, "rb");
  BGZF *input_file = bgzf_open(input_filename, "r");
  if (!input_file)
  {
    perror("Failed to open input file");
    return EXIT_FAILURE;
  }

  // Open the output file in binary write mode using BGZF
  bgzf_output = bgzf_open(output_filename, "wb");
  if (!bgzf_output)
  {
    // fprintf(stderr, "Failed to open output file with BGZF\n");
    log_error("Failed to open output file with BGZF");
    // fclose(input_file);
    bgzf_close(input_file);
    return EXIT_FAILURE;
  }

  // Read from the input file and write compressed data to the output file
  // while ((read_len = fread(buffer, 1, sizeof(buffer), input_file)) > 0)
  while ((read_len = bgzf_read(input_file, buffer, sizeof(buffer))) > 0)
  {
    if (bgzf_write(bgzf_output, buffer, read_len) < 0)
    {
      // fprintf(stderr, "Failed to write to BGZF file\n");
      log_error("Failed to write to BGZF file");
      bgzf_close(bgzf_output);
      // fclose(input_file);
      bgzf_close(input_file);
      return EXIT_FAILURE;
    }
  }

  // Close files
  bgzf_close(bgzf_output);
  // fclose(input_file);
  bgzf_close(input_file);

  return EXIT_SUCCESS;
}


// This program takes an input tfasta filename as an argument
// And do one of the following:
// if the file is not compressed, then it compresses the file using BGZF and creates an index file
// if the file is compressed, then it only creates an index file
// if the tfasta file is v1 then it converts it to v2 and compress the new file and create an index file
// if the tfasta file is v2 then it only compress the file and create an index file








/**
 * Converts a TFA v1 file to a TFA v2 file.
 *
 * @param input_filename The path to the input TFA v1 file.
 * @param v2_filename The path to the output TFA v2 file.
 * @return Returns 0 on success, -1 on failure.
 */
int convert_tfa_v1_to_v2(const char *input_filename, const char *v2_filename)
{
  // Open the input file for reading using BGZF
  BGZF *input_fp = bgzf_open(input_filename, "r");
  if (!input_fp)
  {
    log_error("Failed to open input file %s\n", input_filename);
    return 1;
  }

  // Open the output file for writing using BGZF
  BGZF *output_fp;
  output_fp = bgzf_open(v2_filename, "w");
  if (!output_fp)
  {
    log_error("Failed to open output file %s\n", v2_filename);
    bgzf_close(input_fp);
    return 1;
  }

  // Prepare to read lines using kstring_t
  kstring_t str = {0, 0, NULL};
  int ret;

  bool in_header = true;

  // Write the TFAv2.0 header to the output file
  const char *v2_header = "##fileformat=TFAv2.0\n";
  if (bgzf_write(output_fp, v2_header, strlen(v2_header)) < 0)
  {
    log_error("Error writing the TFAv2.0 header to the output file\n");
    bgzf_close(input_fp);
    bgzf_close(output_fp);
    return 1;
  }

  // copy the rest of the lines from the input file to the output file
  while ((ret = bgzf_getline(input_fp, '\n', &str)) >= 0)
  {
    // if lins start with # then it is a header line
    if (str.s[0] == '#')
    {
      in_header = true;
    }
    else
    {
      in_header = false;
    }
    // Check if the line is a header line
    if (!in_header)
    {
      // modify the line to replace ':' with '\t'
      for (int i = 0; i < str.l; i++)
      {
        if (str.s[i] == ':')
        {
          str.s[i] = '\t';
        }
      }
    }
    else
    {
      // Check if the line is the end of the header start with "#CHR:POSITION"
      if (strncmp(str.s, "#CHR:POSITION", 13) == 0 || strncmp(str.s, "#CHROMOSOME", 11 )   == 0)
      {
        // modify the line to replace ':' with '\t'
        for (int i = 0; i < str.l; i++)
        {
          if (str.s[i] == ':')
          {
            str.s[i] = '\t';
          }
        }
        in_header = false; // Set in_header to false to indicate the end of the header
      }
      else if (strncmp(str.s, "#NAMES", 6) == 0) {
        // skip and do nothing
      }
      else 
      {
        // modify the line to replace ':' with '\t'
        for (int i = 0; i < str.l; i++)
        {
          if (str.s[i] == ':')
          {
            str.s[i] = '\t';
          }
        }
        in_header = false; // Set in_header to false to indicate the end of the header
      }
    }

    // append new line character to the line
    kputsn("\n", 1, &str);

    if (bgzf_write(output_fp, str.s, str.l) < 0)
    {
      log_error("Error writing to the output file\n");
      bgzf_close(input_fp);
      bgzf_close(output_fp);
      return 1;
    }
  }
  //
  // Check for a read error

  if (ret < -1)
  { // bgzf_getline returns -1 for EOF, less than -1 for errors
    log_error("Error reading input file\n");
    free(str.s);
    bgzf_close(input_fp);
    bgzf_close(output_fp);
    return 1;
  }

  // Free the allocated kstring
  free(str.s);

  // Close the input and output files
  bgzf_close(input_fp);
  bgzf_close(output_fp);

  return 0;
}


int convert_wtfa_v1_to_v2(const char *input_filename, const char *v2_filename)
{
  // Open the input file for reading using BGZF
  BGZF *input_fp = bgzf_open(input_filename, "r");
  if (!input_fp)
  {
    log_error("Failed to open input file %s\n", input_filename);
    return 1;
  }

  // Open the output file for writing using BGZF
  BGZF *output_fp;
  output_fp = bgzf_open(v2_filename, "w");
  if (!output_fp)
  {
    log_error("Failed to open output file %s\n", v2_filename);
    bgzf_close(input_fp);
    return 1;
  }

  // Prepare to read lines using kstring_t
  kstring_t str = {0, 0, NULL};
  int ret;


  // Write the TFAv2.0 header to the output file
  const char *v2_header = "##fileformat=wTFAv2.0\n";
  if (bgzf_write(output_fp, v2_header, strlen(v2_header)) < 0)
  {
    log_error("Error writing the TFAv2.0 header to the output file\n");
    bgzf_close(input_fp);
    bgzf_close(output_fp);
    return 1;
  }

  // copy the rest of the lines from the input file to the output file
  while ((ret = bgzf_getline(input_fp, '\n', &str)) >= 0)
  {
    // modify the line to replace ':' with '\t'
    for (int i = 0; i < str.l; i++)
    {
      if (str.s[i] == ':')
      {
        str.s[i] = '\t';
      }
    }

    // append new line character to the line
    kputsn("\n", 1, &str);

    if (bgzf_write(output_fp, str.s, str.l) < 0)
    {
      log_error("Error writing to the output file\n");
      bgzf_close(input_fp);
      bgzf_close(output_fp);
      return 1;
    }
  }
  //
  // Check for a read error

  if (ret < -1)
  { // bgzf_getline returns -1 for EOF, less than -1 for errors
    log_error("Error reading input file\n");
    free(str.s);
    bgzf_close(input_fp);
    bgzf_close(output_fp);
    return 1;
  }

  // Free the allocated kstring
  free(str.s);

  // Close the input and output files
  bgzf_close(input_fp);
  bgzf_close(output_fp);

  return 0;
}

// function to create an index file
/**
 * Creates an index for the given input file and saves it to the specified output file.
 *
 * @param input_filename The path to the input file.
 * @param output_filename The path to the output file where the index will be saved.
 * @param n_threads The number of threads to use for creating the index.
 * @param min_shift The minimum shift value for the index, 0 for TBI and use 14 for CSI. CSI not implemented yet.
 * @param conf The configuration settings for creating the index.
 * @return Returns 0 on success, or a negative value on failure.
 */
int create_index(
    const char *input_filename,
    const char *output_filename,
    int n_threads,
    int min_shift,
    const tbx_conf_t conf)
{
  tbx_t *tbx;
  BGZF *fp;
  int ret;
  if ((fp = bgzf_open(input_filename, "r")) == 0)
    return -1;
  if (n_threads)
    bgzf_mt(fp, n_threads, 256);
  if (bgzf_compression(fp) != bgzf)
  {
    log_error("%s is not compressed with bgzf format", input_filename);
    bgzf_close(fp);
    return -2;
  }
  tbx = tbx_index(fp, min_shift, &conf);
  bgzf_close(fp);
  if (!tbx)
    return -1;
  ret = hts_idx_save_as(tbx->idx, input_filename, NULL, HTS_FMT_TBI);
  tbx_destroy(tbx);
  return ret;
}

// set program name as tfa_index
const char *program_name = "tfa_index";

int print_usage(char *program_name, int status)
{
  fprintf(stderr, "Usage: %s  [options] <input.tfa|input.tfa.bgz|input.tfa.gz|weights.txt|weights.txt.gz>\n", program_name);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  --version\n");
  fprintf(stderr, "  --help\n");
  fprintf(stderr, "  --threads <int>\n");
  fprintf(stderr, "  --force\n");
  fprintf(stderr, "  --weight\n   the input file is a weight file. In this case create and index for it. and convert it to the correct format\n");

  fprintf(stderr, "  --output FILE           set custom name to the output file, only used when converting from TFAv1 to TFAv2 or compressing the input file\n");
  // NOTE :: tabix require the input file to be compressed with bgzf
  // fprintf(stderr, "  --no-compression     Do not compress, Create Index only \n");
  fprintf(stderr, "  --min-shift <int>       set minimal interval size for CSI indices to 2^INT recommended value is 14, 0 means build tbi index. Default [0]\n");

  return status;
}

int tfa_main(char *input_fname, tfa_index_args_t *args){
  // make sure the input file is not a weight file, using assert here as this should not happen
  assert(!args->is_weight);



  int ret = 0;
  // ##########
  const char * input_basename = get_basename(input_fname, false);
  // get file name without extension
  const char *input_extension = get_extension(input_fname);

  // if we are compressing the file we should use the output filename as the basename
  kstring_t output_fname = KS_INITIALIZE; // {0, 0, NULL};
  char *output_basename = NULL;
  char *output_extension = NULL;

  // check index file
  kstring_t idx_fname = KS_INITIALIZE;

  tfa_info_t tfa_info;
  memset(&tfa_info, 0, sizeof(tfa_info_t));

  // check if input file exists
  struct stat input_sbuf;
  if (  !(stat(input_fname, &input_sbuf) == 0 && S_ISREG(input_sbuf.st_mode) ) )
  {
    log_error("[%s] File '%s' does not exist. Abort.", __func__, input_fname);
    return 1;
  }

  const char *format = check_tfa_format(input_fname, &tfa_info);
  if (format == NULL)
  {
    log_error("Failed to check tfasta format");
    return 1;
  }
  log_info("TFA format: %s", format);

  


  // char *output_path = get_filepath(args->output_fname);

  if (args->output_fname != NULL)
  {
    
    // if output_basename is empty then use the input basename
    
    // args.output_fname end with tfa.bgz or tfa.gz then use it as it is
    // else add .tfa.bgz to the output_fname
    if (ends_with(args->output_fname, ".tfa.bgz")  || ends_with(args->output_fname, ".tfa.gz") )
    {
      ksprintf(&output_fname, "%s", args->output_fname);
    }
    else
    {
      output_basename = get_basename(args->output_fname, true);
      
      // if args.output_fname end with tfa or tfasta then add .bgz to the output_fname
      if (ends_with(args->output_fname, ".tfa")  || ends_with(args->output_fname, ".tfasta") )
      {
        ksprintf(&output_fname, "%s.gz", args->output_fname);
      }
      else
      {
        if (strlen(output_basename) == 0)
        {
          ksprintf(&output_fname, "%s%s.tfa.gz", args->output_fname, input_basename);
        }
        else
        {
          ksprintf(&output_fname, "%s.tfa.gz", args->output_fname);
        }
      }
    }
  }
  else
  {
    // make a copy of the input basename and extension
    // here we use the input basename and extension
    output_basename = strdup(tfa_info.basename);
    output_extension = strdup(tfa_info.extension);
    if (tfa_info.version == 1)
    {
      // to avoid overwriting the input file
      // construct the output filename with tfav2 as suffix before the extension
      ksprintf(&output_fname, "%s_TFAv2.%s.gz", output_basename, output_extension);
    }
    else
    {

      if (!tfa_info.is_compressed)
      {
        // need to compress the file
        ksprintf(&output_fname, "%s.%s.gz", output_basename, output_extension);
      }
      else if (tfa_info.is_compressed && !tfa_info.is_bgzf)
      {
        // need to change the extension to bgz
        ksprintf(&output_fname, "%s.%s.bgz", output_basename, output_extension);
      }
      else
      {
        // file is compressed and in bgzf format
        // output_fname should be the same as the input file
        ksprintf(&output_fname, "%s", input_fname);
      }
    }
  }

  // version 1
  struct stat output_sbuf;
  if (tfa_info.version == 1)
  {
    // input_fname = output_fname.s add a suffix to the output filename before the extension

    if (stat(output_fname.s, &output_sbuf) == 0 && S_ISREG(output_sbuf.st_mode))
    {
      // if the two files are the same then change the output_fname and add TFAv2 to the file name
      if (output_sbuf.st_ino == input_sbuf.st_ino)
      {
        log_warn("[%s] File '%s' exists. adding TFAv2 to the file name.", __func__, output_fname.s);
        ksprintf(&output_fname, "%s_TFAv2.%s.bgz", output_basename, output_extension);
      }
    }
  }

  if (!args->is_force && stat(output_fname.s, &output_sbuf) == 0 && S_ISREG(output_sbuf.st_mode))
  {
    log_error("[%s] File '%s' exists. Please apply '-f' to overwrite. Abort.", __func__, output_fname.s);
    ret = 1;
    goto tfa_main_end;
  }

  if (strcmp(format, TFAv1.version) == 0)
  {
    // log_info("File is in TFAv1.0 format");
  }
  else if (strcmp(format, TFAv2.version) == 0)
  {
    // log_info("File is in TFAv2.0 format");
  }
  else
  {
    log_error("Unknown TFA version format");
    ret = 2;
    goto tfa_main_end;
  }
  
  if (!args->min_shift)
    ksprintf(&idx_fname, "%s.tbi", output_fname.s);
  else
    ksprintf(&idx_fname, "%s.csi", output_fname.s);

  if (!args->is_force && stat(idx_fname.s, &output_sbuf) == 0 && S_ISREG(output_sbuf.st_mode))
  {
    log_warn("[%s] Index file '%s' exists. Please apply '-f' to overwrite. Abort.", __func__, idx_fname.s);
    ret = 1;
    goto tfa_main_end;
  }

  // check the format version
  if (strcmp(format, TFAv1.version) == 0)
  {
    log_info("%s File will be converted to TFAv2.0 with bgzf format", input_fname);
    log_info("Output file: %s", output_fname.s);
    convert_tfa_v1_to_v2(input_fname, output_fname.s);
    create_index(output_fname.s, output_fname.s, 2, 0, tfasta_conf);
    log_info("File converted and index was created successfully");
  }
  else if (strcmp(format, TFAv2.version) == 0)
  {
    // log_info("File is in TFAv2.0 format");

    // if the file is not compressed then compress it
    if (!tfa_info.is_compressed || !tfa_info.is_bgzf)
    {
      log_info("Output compressed file: %s", output_fname.s);

      // if log_info("File is compressed but not in bgzf format");
      // TODO :: make sure file names are different

      int status = compress_file(input_fname, output_fname.s);
      if (status != 0)
      {
        log_error("Failed to compress file");
        return 1;
      }
      log_info("File compressed successfully : %s", output_fname.s);
    }

    int status = create_index(output_fname.s, output_fname.s, 2, 0, tfasta_conf);
    if (status != 0)
    {
      log_error("Failed to create index file");
      return 1;
    }
    log_info("Index file created successfully : %s", idx_fname.s);
  }
  else
  {
    log_error("Unknown TFA version format");
    return 1;
  }

tfa_main_end:
  if (output_fname.s)
    ks_free(&output_fname);
  if (ret != 0)
  {
    log_error("%s failed with error code %d", program_name, ret);
  }
  // free resources
  free(tfa_info.basename);
  // free(tfa_info.extension);
  // free(tfa_info.compressed_extension);
  if(output_basename)
    free(output_basename);
  if(output_extension)
    free(output_extension);
  if(idx_fname.s)
    free(idx_fname.s);

  return ret;
}
int wtfa_main(char *input_fname, tfa_index_args_t *args){
  // make sure the input file is a weight file
  assert(args->is_weight);


  int ret = 0;
  tfa_info_t tfa_info;


  memset(&tfa_info, 0, sizeof(tfa_info_t));

   // check index file
  kstring_t idx_fname = KS_INITIALIZE;
  // ##########
  const char * basename = get_basename(input_fname, false);
  // get file name without extension
  const char * extension = get_extension(input_fname);

  // if we are compressing the file we should use the output filename as the basename
  kstring_t output_fname = KS_INITIALIZE; // {0, 0, NULL};
  char *output_basename = NULL;
  char *output_extension = NULL;
  // check if input file exists
  struct stat input_sbuf;
  if (  !(stat(input_fname, &input_sbuf) == 0 && S_ISREG(input_sbuf.st_mode)) )
  {
    log_error("[%s] File '%s' does not exist. Abort.", __func__, input_fname);
    return 1;
  }

  const char *format = check_wtfa_format(input_fname, &tfa_info);
  if (format == NULL)
  {
    log_error("Failed to check tfasta weight format");
    return 1;
  }
  log_info("weight TFA format: %s", format);

  


  // char *output_path = get_filepath(args->output_fname);

    if (args->output_fname != NULL)
    {
      //char * output_path = get_filepath(args->output_fname);
      output_basename = get_basename(args->output_fname, true);
      // if output_basename is empty then use the input basename
      
      // args.output_fname end with .bgz or .gz then use it as it is
      // else add .tfa.bgz to the output_fname
      if (ends_with(args->output_fname, ".bgz")  || ends_with(args->output_fname, ".gz") )
      {
        ksprintf(&output_fname, "%s", args->output_fname);
      }
      else
      {

        if (strlen(output_basename) == 0)
        {
          ksprintf(&output_fname, "%s%s.gz", args->output_fname,basename);
        } 
        else {
          ksprintf(&output_fname, "%s.gz", args->output_fname);
        }
        
      }
      
    }
    else
    {
      // make a copy of the input basename and extension
      // here we use the input basename and extension
      output_basename = strdup(tfa_info.basename);
      if (tfa_info.version == 1)
      {
        // if input_fname end with gz then change the name 
        if(ends_with(input_fname, ".gz") )
        {
          ksprintf(&output_fname, "%s_wTFAv2.%s.gz", output_basename, extension   );
        }
        else
        {
          ksprintf(&output_fname, "%s.gz", input_fname);
        }
        
      }
      else
      {

        if (!tfa_info.is_compressed)
        {
          // need to compress the file
          ksprintf(&output_fname, "%s.%s.gz", output_basename, output_extension);
        }
        else if (tfa_info.is_compressed && !tfa_info.is_bgzf)
        {
          // need to change the extension to bgz
          ksprintf(&output_fname, "%s.%s.bgz", output_basename, output_extension);
        }
        else
        {
          // file is compressed and in bgzf format
          // output_fname should be the same as the input file
          ksprintf(&output_fname, "%s", input_fname);
        }
      }
    }
  

  // version 1
  struct stat output_sbuf;
  if (tfa_info.version == 1)
  {
    // input_fname = output_fname.s add a suffix to the output filename before the extension

    if (stat(output_fname.s, &output_sbuf) == 0 && S_ISREG(output_sbuf.st_mode))
    {
      // if the two files are the same then change the output_fname and add TFAv2 to the file name
      if (output_sbuf.st_ino == input_sbuf.st_ino)
      {
        log_warn("[%s] File '%s' exists. adding wTFAv2 to the file name.", __func__, output_fname.s);
        ksprintf(&output_fname, "%s_wTFAv2.%s.gz", output_basename, output_extension);
      }
    }
  }

  if (!args->is_force && stat(output_fname.s, &output_sbuf) == 0 && S_ISREG(output_sbuf.st_mode))
  {
    log_error("[%s] File '%s' exists. Please apply '-f' to overwrite. Abort.", __func__, output_fname.s);
    ret = 1;
    goto wtfa_main_end;
  }

  if (strcmp(format, wTFAv1.version) == 0)
  {
    // log_info("File is in TFAv1.0 format");
  }
  else if (strcmp(format, wTFAv2.version) == 0)
  {
    // log_info("File is in TFAv2.0 format");
  }
  else
  {
    log_error("Unknown wTFA version format");
    ret = 2;
    goto wtfa_main_end;
  }
 
  if (!args->min_shift)
    ksprintf(&idx_fname, "%s.tbi", output_fname.s);
  else
    ksprintf(&idx_fname, "%s.csi", output_fname.s);

  if (!args->is_force && stat(idx_fname.s, &output_sbuf) == 0 && S_ISREG(output_sbuf.st_mode))
  {
    log_warn("[%s] Index file '%s' exists. Please apply '-f' to overwrite. Abort.", __func__, idx_fname.s);
    ret = 1;
    goto wtfa_main_end;
  }

  // check the format version
  if (strcmp(format, wTFAv1.version) == 0)
  {
    log_info("%s File will be converted to wTFAv2.0 with bgzf format", input_fname);
    log_info("Output file: %s", output_fname.s);
    convert_wtfa_v1_to_v2(input_fname, output_fname.s);
    create_index(output_fname.s, output_fname.s, 2, 0, tfasta_conf);
    log_info("File converted and index was created successfully");
  }
  else if (strcmp(format, wTFAv2.version) == 0)
  {
    // log_info("File is in TFAv2.0 format");

    // if the file is not compressed then compress it
    if (!tfa_info.is_compressed || !tfa_info.is_bgzf)
    {
      log_info("Output compressed file: %s", output_fname.s);

      // if log_info("File is compressed but not in bgzf format");
      // TODO :: make sure file names are different

      int status = compress_file(input_fname, output_fname.s);
      if (status != 0)
      {
        log_error("Failed to compress file");
        return 1;
      }
      log_info("File compressed successfully : %s", output_fname.s);
    }

    int status = create_index(output_fname.s, output_fname.s, 2, 0, tfasta_conf);
    if (status != 0)
    {
      log_error("Failed to create index file");
      return 1;
    }
    log_info("Index file created successfully : %s", idx_fname.s);
  }
  else
  {
    log_error("Unknown wTFA version format");
    return 1;
  }

wtfa_main_end:
  if (output_fname.s)
    ks_free(&output_fname);
  if (ret != 0)
  {
    log_error("%s failed with error code %d", program_name, ret);
  }
  // free resources
  free(tfa_info.basename);
  // free(tfa_info.extension);
  // free(tfa_info.compressed_extension);
  if(output_basename)
    free(output_basename);
  if(output_extension)
    free(output_extension);
  if(idx_fname.s)
    free(idx_fname.s);

  return ret;
}


int main(int argc, char *argv[])
{

  static const struct option loptions[] =
      {
          {"version", no_argument, NULL, 1},
          {"help", no_argument, NULL, 2},
          {"threads", required_argument, NULL, '@'},
          {"force", no_argument, NULL, 'f'},
          {"weight", no_argument, NULL, 'w'},
          {"output", required_argument, NULL, 'o'},
          // {"no-compression", no_argument, NULL, 'n'},
          {"min-shift", required_argument, NULL, 'm'},

      };
  int c, ret = 0;

  char *reheader = NULL;
  tfa_index_args_t args;
  memset(&args, 0, sizeof(tfa_index_args_t));

  char *tmp;
  while ((c = getopt_long(argc, argv, "wfm:o:@:", loptions, NULL)) >= 0)
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
    case 'm':
      args.min_shift = atoi(optarg);
      break;
    case '@':
      args.threads = atoi(optarg);
      break;
    case 'f':
      args.is_force = 1;
      break;
    case 'w':
      args.is_weight = 1;
      break;
    // case 'n':
    //   args.no_compression = 1;
    //   break;
    case 'o':
      args.output_fname = optarg;
      break;
    default:
      fprintf(stderr, "Unknown option: %s\n", optarg);
      return 1;
    }
  }

  // if ( optind==argc ) return print_usage(argv[0] ,EXIT_FAILURE);

  // input file is required
  if (optind == argc)
  {
    // print_usage(argv[0]);
    log_error("Input file is required");
    return print_usage(argv[0], EXIT_FAILURE);

    // return EXIT_FAILURE;
  }
  // get the input filename
  char *input_fname = argv[optind];
  // check file as required
  if (input_fname == NULL)
  {
    log_error("Input file is required");
    // fprintf(stderr, "Input tfasta file is required\n");
    return print_usage(argv[0], EXIT_FAILURE);
  }

  
  // echo program started
  log_info("%s started", program_name);
  // echo command line arguments as received
  // construct a string of the command line arguments
  kstring_t cmd_args = KS_INITIALIZE;
  // cat program name to the command line arguments
  // cmd_args = strcat(cmd_args, program_name);
  for (int i = 0; i < argc; i++)
  {
    ksprintf(&cmd_args, "%s ", argv[i]);
  }
  log_info("%s", ks_c_str(&cmd_args));
  // free resources
  ks_free(&cmd_args);

  // echo filename
  log_info("Filename: %s", input_fname);


  if (args.is_weight) {
    return wtfa_main(input_fname, &args);
  } else {
    return tfa_main(input_fname, &args);
  }
  

  
}
