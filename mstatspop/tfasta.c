#include "tfasta.h"
#include "log.h"

// Path: mstatspop/tfasta.c

// Configuration for the tbx index

static inline char get_DNA_char (char *c) {
  switch (*c)
  {
  case 'T':
    return '1';
  case 't':
    return '1';
  case 'U':
    return '1';
  case 'u':
    return '1';
  case 'C':
    return '2';
  case 'c':
    return '2';
  case 'G':
    return '3';
  case 'g':
    return '3';
  case 'A':
    return '4';
  case 'a':
    return '4';
  case 0:
    return 0;
  case -1:
    return 0;
  case 10:
    return 0;
  case 13:
    return 0;
  case 32:
    return 0;
  case 9:
    return 0;
  case 'N':
    return '5';
  case 'n':
    return '5';
  case '?':
    return '5';
  case '-':
    return '6';
  case 'W':
    return 'w';
  case 'w':
    return 'w';
  case 'M':
    return 'm';
  case 'm':
    return 'm';
  case 'R':
    return 'r';
  case 'r':
    return 'r';
  case 'Y':
    return 'y';
  case 'y':
    return 'y';
  case 'K':
    return 'k';
  case 'k':
    return 'k';
  case 'S':
    return 's';
  case 's':
    return 's';
  case 'b':
    return '5';
  case 'B':
    return '5';
  case 'd':
    return '5';
  case 'D':
    return '5';
  case 'h':
    return '5';
  case 'H':
    return '5';
  case 'v':
    return '5';
  case 'V':
    return '5';
  default:
    return -1;
  }
}




// Function to get the name of the file format
static const char* get_format_name(int format) {
    switch (format) {
        case TBX_UCSC: return "UCSC (BED)";
        case TBX_SAM: return "SAM";
        case TBX_VCF: return "VCF";
        case TBX_GENERIC: return "Generic";
        default: return "Unknown";
    }
}

const char *check_tfa_format(const char *filename, tfa_info_t *tfa_info)
{
  // Open the BGZF file for reading
  BGZF *bgzf_fp = bgzf_open(filename, "r");
  if (!bgzf_fp)
  {
    fprintf(stderr, "Failed to open file %s\n", filename);
    return NULL;
  }
  tfa_info->is_compressed = bgzf_fp->is_compressed;

  if (bgzf_compression(bgzf_fp) == bgzf)
  {
    tfa_info->is_bgzf = 1;
  }
  else
  {
    tfa_info->is_bgzf = 0;
  }

  // allocate memory for the basename to copy the filename
  tfa_info->basename =  get_basename(filename, false);
  if (tfa_info->is_compressed)
  {
    // detect compression format
    // detect compression extension
    // remove the extension to get file name without the compression extension
    
    // get the extension
  char *extension = strrchr(tfa_info->basename, '.');
   if (extension != NULL && (strcmp(extension, ".gz") == 0 || strcmp(extension, ".bgz") == 0)) {
      tfa_info->compressed_extension = extension+1;
      *extension = '\0';  // Remove the extension by setting a null character
      // get the the actual file extension without the compression extension
      extension = strrchr(tfa_info->basename, '.');
      if (extension != NULL) {
        *extension = '\0';  // Remove the extension by setting a null character
        tfa_info->extension = extension + 1;
      }
    }
  }
  else {
    // get the extension
    char *extension = strrchr(tfa_info->basename, '.');
    if (extension != NULL) {
      *extension = '\0';  // Remove the extension by setting a null character
      tfa_info->extension = extension + 1;
    }
  }
  
  // Prepare to read lines using kstring_t
  kstring_t str = {0, 0, NULL};
  int ret = bgzf_getline(bgzf_fp, '\n', &str); // Read the first line

  if (ret < 0)
  {
    // If there is an error reading the first line, assume TFAv1.0
    bgzf_close(bgzf_fp);
    return "TFAv1.0";
  }
  // assume TFAv1.0 by default
  tfa_info->version = 1;
  // Check if the first line matches the TFA format version 2.0
  const char *expected_header = "##fileformat=TFAv2.0";
  if (strncmp(str.s, expected_header, strlen(expected_header)) == 0)
  {
    bgzf_close(bgzf_fp);
    free(str.s);
    tfa_info->version = 2;
    return "TFAv2.0"; // Return the format version
  }

  // If the first line does not match and is not empty, assume TFAv1.0
  bgzf_close(bgzf_fp);
  free(str.s);
  return "TFAv1.0";
}

const char *check_wtfa_format(const char *filename, tfa_info_t *tfa_info)
{
  // Open the BGZF file for reading
  BGZF *bgzf_fp = bgzf_open(filename, "r");
  if (!bgzf_fp)
  {
    fprintf(stderr, "Failed to open file %s\n", filename);
    return NULL;
  }
  tfa_info->is_compressed = bgzf_fp->is_compressed;

  if (bgzf_compression(bgzf_fp) == bgzf)
  {
    tfa_info->is_bgzf = 1;
  }
  else
  {
    tfa_info->is_bgzf = 0;
  }

  // allocate memory for the basename to copy the filename
  tfa_info->basename = strdup(filename);
  if (tfa_info->is_compressed)
  {
    // detect compression format
    // detect compression extension
    // remove the extension to get file name without the compression extension
    
    // get the extension
  char *extension = strrchr(tfa_info->basename, '.');
   if (extension != NULL && (strcmp(extension, ".gz") == 0 || strcmp(extension, ".bgz") == 0)) {
      tfa_info->compressed_extension = extension+1;
      *extension = '\0';  // Remove the extension by setting a null character
      // get the the actual file extension without the compression extension
      extension = strrchr(tfa_info->basename, '.');
      if (extension != NULL) {
        *extension = '\0';  // Remove the extension by setting a null character
        tfa_info->extension = extension + 1;
      }
    }

  }
  tfa_info->is_weight_file = 1;
  // Prepare to read lines using kstring_t
  kstring_t str = {0, 0, NULL};
  int ret = bgzf_getline(bgzf_fp, '\n', &str); // Read the first line

  if (ret < 0)
  {
    // If there is an error reading the first line, assume TFAv1.0
    bgzf_close(bgzf_fp);
    return "wTFAv1.0";
  }
  // assume TFAv1.0 by default
  tfa_info->version = 1;
  // Check if the first line matches the TFA format version 2.0
  const char *expected_header = "##fileformat=wTFAv2.0";
  if (strncmp(str.s, expected_header, strlen(expected_header)) == 0)
  {
    bgzf_close(bgzf_fp);
    free(str.s);
    tfa_info->version = 2;
    return "wTFAv2.0"; // Return the format version
  }

  // If the first line does not match and is not empty, assume TFAv1.0
  bgzf_close(bgzf_fp);
  free(str.s);
  return "wTFAv1.0";
}

// error codes

/**
 * Initializes a tfasta_file structure with the given tfasta_fname.
 *
 * @param tfasta A pointer to the tfasta_file structure to be initialized.
 * @param tfasta_fname The name of the tfasta file.
 */
int init_tfasta_file(tfasta_file *tfasta, char *tfasta_fname)
{
  if (tfasta->is_initialized)
  {
    log_debug("tfasta_file is already initialized");
    return TFA_ERROR;
  }

  const char *tfa_format = check_tfa_format(tfasta_fname, &tfasta->tfa_info);
  if (tfa_format == NULL)
  {
    log_error("Failed to check tfasta format for file: %s", tfasta_fname);
    return TFA_ERROR;
  }

  if (strcmp(tfa_format, TFAv2.version) != 0) {
    log_error("Unsupported TFA format: %s", tfa_format);
    log_error("Current Supported TFA format is: %s", TFAv2.version);
    return TFA_INVALID_FORMAT;
  }

  // return 0 for error, 1 for success
  tfasta->tfasta_fname = tfasta_fname;
  tfasta->n_sam = 0;
  tfasta->names = NULL;
  tfasta->tbx = tbx_index_load(tfasta_fname);
  if (tfasta->tbx == NULL)
  {
    log_error("Failed to load index for: %s", tfasta_fname);
    return TFA_INVALID_INDEX;
  }
  tfasta->fp = hts_open(tfasta_fname, "r");
  if (tfasta->fp == NULL)
  {
    log_error("Failed to open file: %s", tfasta_fname);
    tbx_destroy(tfasta->tbx);
    return TFA_ERROR;
  }

  // allocate memory for samples names
  if ((tfasta->names = (char **)calloc(maxsam, sizeof(char *))) == 0)
  {
    log_fatal("Error: memory not allocated, samples names init_tfasta_file");
    return TFA_ERROR;
  }

  // read and fill samples names
  kstring_t str = {0, 0, NULL};
  int len;
  while ((len = hts_getline(tfasta->fp, KS_SEP_LINE, &str)) >= 0)
  {
    if (str.s[0] != '#')
      break; // Stop at the first non-header line
    if (strstr(str.s, "#NAMES:") != 0)
    {
      // collect names
      int nseq = 0;
      char *cc = strtok(str.s, ">\n\r ");
      while (cc != NULL)
      {
        if (strstr(cc, "#NAMES:") == 0)
        {
          log_debug("Adding Sample : %s", cc);
          // allocate memory for the sample name with the size of cc
          // restrict the size of the sample name to 50 characters:: possible to change this later
          // TODO :: fixme, the 50 characters limit, affect downstream code, sorting and replacing names
          // if ((tfasta->names[nseq] = (char *)calloc(strlen(cc) + 1, sizeof(char))) == 0)
          if ((tfasta->names[nseq] = (char *)calloc(50, sizeof(char))) == 0)
          {
            log_fatal("Error: memory not allocated, samples names init_tfasta_file");
            free(str.s);
            return TFA_ERROR;
          }

          strcpy(tfasta->names[nseq], cc);
          nseq++;
          if (nseq == maxsam)
          {
            maxsam += NSAM_INC;
            if (maxsam > 32767)
            {
              log_error("Sorry, no more than 32767 samples are allowed.");
              free(str.s);
              return TFA_ERROR;
            }
            if ((tfasta->names = (char **)realloc(tfasta->names, maxsam * sizeof(char *))) == 0)
            {
              log_fatal("Error: memory not reallocated, names assigna.1");
              free(str.s);
              return TFA_ERROR;
            }
            // for (int x = nseq; x < maxsam; x++)
            // {
            //   if ((tfasta->names[x] = (char *)calloc(50, sizeof(char))) == 0)
            //   {
            //     log_fatal("Error: memory not reallocated, names assigna.2");
            //     return 0;
            //   }
            // }
          }
        }
        cc = strtok(NULL, ">\n\r ");
      }
      tfasta->n_sam = nseq;
    }
    
  }
  // free kstring_t str
  free(str.s);
  // if n_sam is 0, then the file does not contain any samples 
  if (tfasta->n_sam == 0)
  {
    log_error("No samples or names found in the tfasta: %s", tfasta_fname);
    return TFA_ERROR;
  }

  int nseq;
  const char **seqnames = tbx_seqnames(tfasta->tbx, &nseq);
  if (seqnames == NULL)
  {
    fprintf(stderr, "Failed to retrieve sequence names.\n");
    tbx_destroy(tfasta->tbx);
    return TFA_ERROR;
  }
  // Print sequence names
  for (int i = 0; i < nseq; i++)
  {
    // printf("%s\n", seqnames[i]);
    log_debug("Sequence Name : %s", seqnames[i]);
  }

  // Print sequence names and their interval counts
  log_debug("Sequence names and interval counts:");
  uint64_t n_records;
  uint64_t unmaped;
  for (int i = 0; i < nseq; i++)
  {
    hts_idx_get_stat(tfasta->tbx->idx, i, &n_records, &unmaped);
    printf("\t\t%s: %lu intervals\n", seqnames[i], n_records);
  }

  // Print the format of the indexed file
  int format = tfasta->tbx->conf.preset;
  // log_debug("Format of the indexed file: %s\n", get_format_name(format));

  tfasta->is_initialized = true;
  free(seqnames);
  return TFA_OK;
}

/**
 * Closes the specified tfasta file and free memory.
 *
 * @param tfasta A pointer to the tfasta_file structure representing the tfasta file to be closed.
 */
void close_tfasta_file(tfasta_file *tfasta)
{
  if (tfasta->fp != NULL)
  {
    hts_close(tfasta->fp);
  }
  if (tfasta->tbx != NULL)
  {
    tbx_destroy(tfasta->tbx);
  }
  if (tfasta->names != NULL)
  {
    for (int x = 0; x < tfasta->n_sam; x++)
    {
      free(tfasta->names[x]);
    }
    free(tfasta->names);
  }
  tfasta->is_initialized = false;
}



int init_wtfasta_file(wtfasta_file *wtfasta, char *wtfasta_fname)
{
  if (wtfasta->is_initialized)
  {
    log_debug("wtfasta_file is already initialized");
    return TFA_ERROR;
  }

  const char *tfa_format = check_wtfa_format(wtfasta_fname, &wtfasta->tfa_info);
  if (tfa_format == NULL)
  {
    log_error("Failed to check tfasta format for file: %s", wtfasta_fname);
    return TFA_ERROR;
  }

  if (strcmp(tfa_format, wTFAv2.version) != 0) {
    log_error("Unsupported wTFA format: %s", tfa_format);
    log_error("Current Supported wTFA format is: %s", wTFAv2.version);
    return TFA_INVALID_FORMAT;
  }

  // return 0 for error, 1 for success
  wtfasta->wtfasta_fname = wtfasta_fname;
  wtfasta->tbx = tbx_index_load(wtfasta_fname);
  if (wtfasta->tbx == NULL)
  {
    log_error("Failed to load index for: %s", wtfasta_fname);
    return TFA_INVALID_INDEX;
  }
  wtfasta->fp = hts_open(wtfasta_fname, "r");
  if (wtfasta->fp == NULL)
  {
    log_error("Failed to open file: %s", wtfasta_fname);
    tbx_destroy(wtfasta->tbx);
    return TFA_ERROR;
  }

  wtfasta->is_initialized = true;
  return TFA_OK;
}


/**
 * Closes the wtfasta file.
 *
 * @param wtfasta A pointer to the wtfasta_file structure representing the file to be closed.
 */
void close_wtfasta_file(wtfasta_file *wtfasta)
{
  if (wtfasta->fp != NULL)
  {
    hts_close(wtfasta->fp);
  }
  if (wtfasta->tbx != NULL)
  {
    tbx_destroy(wtfasta->tbx);
  }
  wtfasta->is_initialized = false;
}

int get_interval_length(const char* chr_name,tbx_t *tbx) {
  // get the length of the chromosome / scaffold
  // get the length of the sequence
  int nseq;
  const char **seqnames = tbx_seqnames(tbx, &nseq);
  if (seqnames == NULL)
  {
    log_error( "Failed to retrieve sequence names.");
    return -2;
  }
  // find the index of 
  for (int i = 0; i < nseq; i++)
  {
    if (strcmp(seqnames[i], chr_name) == 0)
    {
      uint64_t n_records;
      uint64_t unmaped;
      hts_idx_get_stat(tbx->idx, i, &n_records, &unmaped);
      free(seqnames);
      return n_records;
    }
  }

  free(seqnames);
  return 0;
}


// n_sam initialized once when first is 0
// names initialized once when first is 0
// the actual output of this function is the DNA_matr, n_site
// return 0 if error
// return 1 if success 
// return -1 if end of file : since we are using an index and we are reading a specific region of the file, we should not reach the end of the file
// we can not reach the end of the file, we should always read the region we are interested in
// or retun n_site = 0 if the region we are interested in is not in the file
int read_tfasta_DNA(
    tfasta_file *tfasta,
    char *chr_name,
    long int init_site,
    long int end_site,
    int *n_sam,
    long long *n_site,
    //char ***names,
    char **DNA_matr    
    )
{
  // n_sam : number of samples
  // n_site : number of sites, sites correspond to positions in the DNA sequences
  // names : names of the samples
  // DNA_matr : matrix of DNA sequences
  // matrix_pol_tcga :  ?? not used here ??
  // chr_name : name of the chromosome / scaffold to get data for
  // first : First time we call this function
  // length : length of the chromosome / scaffold

  // this get the data for chr_name sequence , start from init_site to end_site
  // could end earlier if the length of the sequence is smaller than end_site

  // return 0 for error, 1 for success but did not reach end of sequences  and -1 for end of sequence

  // DNA_matr structure :
  // DNA_matr is a matrix of DNA sequences,
  // each row is a site and each column is a variant value for a sample
  // example :
  // Nucleotides : A, C, G, T, N, -
  // Sample 1 : A, C, G, T, N, -
  // ACTGAAAAAAAAAACCA
  // ACCCCCCCCCCCCCCCC
  // means that the first site is A for sample 1 and C for sample 2 , and T for sample 3 and G for sample 4 , etc
  // then the second site is A for sample 1 and C for sample 2 , and C for sample 3 and C for sample 4 , etc
  // the matrix is transposed later, so each row is a sample and each column is a site

  // some static variables to keep track of the position in the file
  // TODO : check if we can remove them
  // sample increment to reallocate memory for samples names

  // static tfasta_file tfasta;
  // memset(&tfasta, 0, sizeof(tfasta_file));
  // char *tfasta_fname = file_input_gz->file_name;
  if(end_site == -1) {
    // reading all 
    log_debug("Reading DNA data for %s from %ld to end of sequence", chr_name, init_site);
    end_site =  get_interval_length(chr_name, tfasta->tbx);
    if (end_site <=0)
    {
      log_error("Failed to get the length of the sequence for %s", chr_name);
      return 0;
    }
  }
  else {
    log_debug("Reading DNA data for %s from %ld to %ld", chr_name, init_site, end_site);
  }
  
  
  static long int position = 1;

  // names :: is pre allocated to 128 samples
  // each sample name is allocated to 50 characters

  char *c;
  
  // // return 0 for error, 1 for success
  // // Open the file with tabix index
  // tbx_t *tbx = tbx_index_load(filename);
  // if (tbx == NULL)
  // {
  //   // fprintf(stderr, "Failed to load index for: %s\n", filename);
  //   log_error("Failed to load index for: %s", filename);
  //   return 0;
  // }
  // htsFile *fp = hts_open(filename, "r");
  // if (fp == NULL)
  // {
  //   // fprintf(stderr, "Failed to open file: %s\n", filename);
  //   log_error("Failed to open file: %s", filename);
  //   tbx_destroy(tbx);
  //   return 0;
  // }

  // if first time we call this function
  // TODO :: this should be moved to the main function and not here

  //init_tfasta_file(&tfasta, tfasta_fname);
  
  // copy the number of samples and names to the output variables
  *n_sam = tfasta->n_sam;
  // *names = tfasta->names;
  *n_site = 0;
  

  // construct and char string with chr_name and init_site and end_site
  // construct a region string in the format chr_name:init_site-end_site
  // allocate memory for the region string dynamically
  // 25 is the length of the string ":1234567890-1234567890"
  // leading to a total of 12 chars for the number which is enough for the length of a chromosome for up to 10^12 == 1 trillion bases
  // 999999999999-999999999999
  char *region = (char *)malloc(strlen(chr_name) + 25);

  sprintf(region, "%s:%ld-%ld", chr_name, init_site, end_site);

  hts_itr_t *iter = tbx_itr_querys(tfasta->tbx, region);

  // no need for the region string anymore
  free(region);

  if (iter == NULL)
  {
    // it is possible that the region is not found in the index
    //fprintf(stderr, "Failed to parse region: %s\n", chr_name);
    log_error("Failed to parse region: %s", chr_name);

    // reset the position to 0
    // free iter
    hts_itr_destroy(iter);
    
    position = 1;
    return 0;
  }

  // allocate
  // DNA_matr2 is a temporary matrix to store the data
  char *DNA_matr2;
  // allocate memory for the matrix by muliple of tfasta.n_sam by init n positions
  // we expect that the number of sites will be in the range of init_site and end_site
  int expected_sites = end_site - init_site + 1;
  long long DNA_matr2_size = tfasta->n_sam * expected_sites;
  if ((DNA_matr2 = (char *)calloc(DNA_matr2_size, sizeof(char))) == 0)
  {
    // fprintf(file_logerr,"\nError: memory not reallocated. get_tfadata.23d \n");
    log_fatal("Error: memory not reallocated, DNA_matr2 get_tfadata.23d");
    free(DNA_matr2);
    // close the file and free memory
    // close_tfasta_file(tfasta);
    return (0);
  }
  


  // tbx_itr_querys will return an iterator to the region in the file
  // it only return what we asked for, so we need to iterate over the iterator to get the data

  kstring_t str = {0, 0, NULL};
  const char *delim = ":\t\n";
  // keep track DNA_matr2 size
  int count = 0;
  while (tbx_itr_next(tfasta->fp, tfasta->tbx, iter, &str) >= 0)
  {
    // if line start with # then it is a comment
    if (str.s[0] == '#')
    {
      // skip comments :: not going to happen Here
      // log_debug("Comment : %s", str.s);
      continue;
    }

    // log_debug("TFA at site %d  : %s", *n_site ,  str.s);
    // tokenize the line to get the data
    char *cc = strtok(str.s, delim);
    int col = 0;
    while (cc != NULL)
    {
      // log_debug("col %d  : %s", col,  cc);

      // col 0 is the name of the sequence
      if(col == 0) {

      }

      // col 1 is the position
      if(col == 1) {
        position = atol(cc);
      }
      
      // col 2 is the nucleotides per sample 
      if(col == 2) {
        // TODO : fill the DNA_matr matrix

        // DOES we require that strlen(cc) == n_sam ??

        // iterate over the nucleotides and fill the matrix
        for (int i = 0; i < strlen(cc); i++)
        {
          // log_debug("Nucleotide %d  : %c", i,  cc[i]);
          // fill the matrix
          // DNA_matr[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '1';
          int DNA_matr2_index = (((long long)tfasta->n_sam * (unsigned long)*n_site) + (unsigned long)i);
          char dna_char =  get_DNA_char(&cc[i]);
          if(dna_char == -1) {
            log_error("Unexpected value in tfa file: position %ld, sample %d \n%c", position, i, cc[i]);
            free(DNA_matr2);
            hts_itr_destroy(iter);
            return (-1);
          }
          if(dna_char > 0) {
            DNA_matr2[DNA_matr2_index] = dna_char;
            count++;
          }
        }

      } 

      // in case we have more than 3 columns
      // we can ignore them for now
      if(col > 2) {
        log_debug("Ignoring unsupported data column %d  : %s", col,  cc);
      }

      // increment the column
      col++;
      // get the next token
      cc = strtok(NULL, delim);
    }
    // increment the site
    *n_site += 1;
    if(*n_site > expected_sites) {
      // need to reallocate memory for the matrix, if our calculations are correct, this should not happen
      log_debug("Reallocating memory for DNA_matr2");
    }
  }
  // if *n_site == 0, then the region is not found in the file
  if(*n_site == 0) {
   
    hts_itr_destroy(iter);
    free(DNA_matr2);
    free(str.s);
    return 1;
  }

  
  free(str.s);

  // possible side effect, memory was allocated for DNA_matr2 before
  free(*DNA_matr);
  // actual size of the matrix as n_site can be less than expected_sites
  // but it should not be more than expected_sites
  DNA_matr2_size = tfasta->n_sam * *n_site;
  if ((*DNA_matr = (char *)calloc(DNA_matr2_size, sizeof(char))) == 0)
  {
    // fprintf(file_logerr,"\nError: memory not reallocated. get_tfadata.23d2 \n");
    log_fatal("Error: memory not reallocated, DNA_matr get_tfadata.23d2");
    free(DNA_matr2);
    hts_itr_destroy(iter);
    return (0);
  }

  for (unsigned long x = 0; x < tfasta->n_sam; x++) {
    for (unsigned long xx = 0; xx < *n_site; xx++) { /*transpose */
      (*DNA_matr)[((*n_site * x) + xx)] =
          DNA_matr2[((tfasta->n_sam * xx) + x)];
    }
  }
  free(DNA_matr2);
  hts_itr_destroy(iter);

  return (1);
}


/**
 * Reads the DNA data from the specified tfasta file for the specified chromosome and region.
 * The data is stored in the DNA_matr matrix.
 * is a matrix of DNA sequences, each row is a site and each column is a variant value for a sample
 * the only difference between this function and read_tfasta_DNA is that this function does not transpose the matrix
 */
int read_tfasta_DNA_lite(
    tfasta_file *tfasta,
    char *chr_name,
    long int init_site,
    long int end_site,
    int *n_sam,
    long long *n_site,
    //char ***names,
    char **DNA_matr    
    )
{
  // n_sam : number of samples
  // n_site : number of sites, sites correspond to positions in the DNA sequences
  // names : names of the samples
  // DNA_matr : matrix of DNA sequences
  // matrix_pol_tcga :  ?? not used here ??
  // chr_name : name of the chromosome / scaffold to get data for
  // first : First time we call this function
  // length : length of the chromosome / scaffold

  // this get the data for chr_name sequence , start from init_site to end_site
  // could end earlier if the length of the sequence is smaller than end_site

  // return 0 for error, 1 for success but did not reach end of sequences  and -1 for end of sequence

  // DNA_matr structure :
  // DNA_matr is a matrix of DNA sequences,
  // each row is a site and each column is a variant value for a sample
  // example :
  // Nucleotides : A, C, G, T, N, -
  // Sample 1 : A, C, G, T, N, -
  // ACTGAAAAAAAAAACCA
  // ACCCCCCCCCCCCCCCC
  // means that the first site is A for sample 1 and C for sample 2 , and T for sample 3 and G for sample 4 , etc
  // then the second site is A for sample 1 and C for sample 2 , and C for sample 3 and C for sample 4 , etc
  // the matrix is transposed later, so each row is a sample and each column is a site

  // some static variables to keep track of the position in the file
  // TODO : check if we can remove them
  // sample increment to reallocate memory for samples names

  // static tfasta_file tfasta;
  // memset(&tfasta, 0, sizeof(tfasta_file));
  // char *tfasta_fname = file_input_gz->file_name;
  if(end_site == -1) {
    // reading all 
    log_debug("Reading DNA data for %s from %ld to end of sequence", chr_name, init_site);
    end_site =  get_interval_length(chr_name, tfasta->tbx);
    if (end_site <=0)
    {
      log_error("Failed to get the length of the sequence for %s", chr_name);
      return 0;
    }
  }
  else {
    log_debug("Reading DNA data for %s from %ld to %ld", chr_name, init_site, end_site);
  }
  
  
  static long int position = 1;

  // names :: is pre allocated to 128 samples
  // each sample name is allocated to 50 characters

  char *c;
  
  // // return 0 for error, 1 for success
  // // Open the file with tabix index
  // tbx_t *tbx = tbx_index_load(filename);
  // if (tbx == NULL)
  // {
  //   // fprintf(stderr, "Failed to load index for: %s\n", filename);
  //   log_error("Failed to load index for: %s", filename);
  //   return 0;
  // }
  // htsFile *fp = hts_open(filename, "r");
  // if (fp == NULL)
  // {
  //   // fprintf(stderr, "Failed to open file: %s\n", filename);
  //   log_error("Failed to open file: %s", filename);
  //   tbx_destroy(tbx);
  //   return 0;
  // }

  // if first time we call this function
  // TODO :: this should be moved to the main function and not here

  //init_tfasta_file(&tfasta, tfasta_fname);
  
  // copy the number of samples and names to the output variables
  *n_sam = tfasta->n_sam;
  // *names = tfasta->names;
  *n_site = 0;
  

  // construct and char string with chr_name and init_site and end_site
  // construct a region string in the format chr_name:init_site-end_site
  // allocate memory for the region string dynamically
  // 25 is the length of the string ":1234567890-1234567890"
  // leading to a total of 12 chars for the number which is enough for the length of a chromosome for up to 10^12 == 1 trillion bases
  // 999999999999-999999999999
  char *region = (char *)malloc(strlen(chr_name) + 25);

  sprintf(region, "%s:%ld-%ld", chr_name, init_site, end_site);

  hts_itr_t *iter = tbx_itr_querys(tfasta->tbx, region);

  // no need for the region string anymore
  free(region);

  if (iter == NULL)
  {
    // it is possible that the region is not found in the index
    //fprintf(stderr, "Failed to parse region: %s\n", chr_name);
    log_error("Failed to parse region: %s", chr_name);

    // reset the position to 0
    // free iter
    hts_itr_destroy(iter);
    
    position = 1;
    return 0;
  }

  // allocate
  // DNA_matr2 is a temporary matrix to store the data
  char *DNA_matr2;
  // allocate memory for the matrix by muliple of tfasta.n_sam by init n positions
  // we expect that the number of sites will be in the range of init_site and end_site
  int expected_sites = end_site - init_site + 1;
  long long DNA_matr2_size = tfasta->n_sam * expected_sites;
  if ((DNA_matr2 = (char *)calloc(DNA_matr2_size, sizeof(char))) == 0)
  {
    // fprintf(file_logerr,"\nError: memory not reallocated. get_tfadata.23d \n");
    log_fatal("Error: memory not reallocated, DNA_matr2 get_tfadata.23d");
    free(DNA_matr2);
    // close the file and free memory
    // close_tfasta_file(tfasta);
    return (0);
  }
  


  // tbx_itr_querys will return an iterator to the region in the file
  // it only return what we asked for, so we need to iterate over the iterator to get the data

  kstring_t str = {0, 0, NULL};
  const char *delim = ":\t\n";
  // keep track DNA_matr2 size
  int count = 0;
  while (tbx_itr_next(tfasta->fp, tfasta->tbx, iter, &str) >= 0)
  {
    // if line start with # then it is a comment
    if (str.s[0] == '#')
    {
      // skip comments :: not going to happen Here
      // log_debug("Comment : %s", str.s);
      continue;
    }

    // log_debug("TFA at site %d  : %s", *n_site ,  str.s);
    // tokenize the line to get the data
    char *cc = strtok(str.s, delim);
    int col = 0;
    while (cc != NULL)
    {
      // log_debug("col %d  : %s", col,  cc);

      // col 0 is the name of the sequence
      if(col == 0) {

      }

      // col 1 is the position
      if(col == 1) {
        position = atol(cc);
      }
      
      // col 2 is the nucleotides per sample 
      if(col == 2) {
        // TODO : fill the DNA_matr matrix

        // DOES we require that strlen(cc) == n_sam ??

        // iterate over the nucleotides and fill the matrix
        for (int i = 0; i < strlen(cc); i++)
        {
          // log_debug("Nucleotide %d  : %c", i,  cc[i]);
          // fill the matrix
          // DNA_matr[(((long long)nseq * (unsigned long)*n_site) + (unsigned long)col)] = '1';
          int DNA_matr2_index = (((long long)tfasta->n_sam * (unsigned long)*n_site) + (unsigned long)i);
          char dna_char =  get_DNA_char(&cc[i]);
          if(dna_char == -1) {
            log_error("Unexpected value in tfa file: position %ld, sample %d \n%c", position, i, cc[i]);
            free(DNA_matr2);
            hts_itr_destroy(iter);
            return (-1);
          }
          if(dna_char > 0) {
            DNA_matr2[DNA_matr2_index] = dna_char;
            count++;
          }
        }

      } 

      // in case we have more than 3 columns
      // we can ignore them for now
      if(col > 2) {
        log_debug("Ignoring unsupported data column %d  : %s", col,  cc);
      }

      // increment the column
      col++;
      // get the next token
      cc = strtok(NULL, delim);
    }
    // increment the site
    *n_site += 1;
    if(*n_site > expected_sites) {
      // need to reallocate memory for the matrix, if our calculations are correct, this should not happen
      log_debug("Reallocating memory for DNA_matr2");
    }
  }
  // if *n_site == 0, then the region is not found in the file
  if(*n_site == 0) {
   
    hts_itr_destroy(iter);
    free(DNA_matr2);
    free(str.s);
    return 1;
  }

  
  free(str.s);

  // possible side effect, memory was allocated for DNA_matr2 before
  free(*DNA_matr);
  *DNA_matr = DNA_matr2;
  hts_itr_destroy(iter);

  return (1);
}


/**
 * Creates an index for the given input file and saves it to the specified output file.
 *
 * @param input_filename The path to the input file.
 * @param n_threads The number of threads to use for creating the index.
 * @param min_shift The minimum shift value for the index, 0 for TBI and use 14 for CSI. CSI not implemented yet.
 * @param conf The configuration settings for creating the index.
 * @return Returns 0 on success, or a negative value on failure.
 */
int create_index(
    const char *input_filename,
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
