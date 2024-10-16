/*

*/
#ifndef MSTATSPOP_TFASTA_H
#define MSTATSPOP_TFASTA_H

// #include "common.h"
// #include <stdio.h>

#include <htslib/tbx.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h> 
#include <htslib/kseq.h>
#include "files_util.h"
#ifdef __cplusplus
extern "C"
{
#endif



// Define a structure to hold version information
typedef struct
{
  const char *version;
  const char *description;
} VersionInfo;

typedef struct
{
  int version;
  char *description;
  char *basename;
  char *extension;
  char *compressed_extension;
  short is_compressed;
  short is_bgzf;
  short is_weight_file; /* this is not a tfa, it is a weight file */
} tfa_info_t;



static const VersionInfo TFAv3 = {"TFAv3.0", "Handling TFA version 3.0"};
static const VersionInfo TFAv2 = {"TFAv2.0", "Handling TFA version 2.0"};
static const VersionInfo TFAv1 = {"TFAv1.0", "Handling TFA version 1.0"};
// weights for a TFA file
static const VersionInfo wTFAv1 = {"wTFAv1.0", "Handling wTFA version 1.0"};
static const VersionInfo wTFAv2 = {"wTFAv2.0", "Handling wTFA version 2.0"};


// error codes for init_tfasta_file
static const int TFA_OK = 1;
static const int TFA_ERROR = 0;
static const int TFA_INVALID_FORMAT = -1;
static const int TFA_INVALID_INDEX = -2;


/*
  *  TFA structure  for  a TFA file
*/
typedef struct 
{
  tfa_info_t tfa_info;
  char *tfasta_fname;
  // samples names
  char **names;
  // number of samples
  long long n_sam;
  tbx_t *tbx;
  htsFile *fp;
  int is_initialized;
} tfasta_file;

/*
  *  Weights TFA file assiciated with a TFA file
*/
typedef struct 
{
  tfa_info_t tfa_info;
  char *wtfasta_fname;
  tbx_t *tbx;
  htsFile *fp;
  int is_initialized;
} wtfasta_file;


/**
 * @brief The number of samples to increment by when resizing an array.
 *
 * This constant represents the number of samples to increment by when resizing an array.
 * It is used in the context of the `tfasta` module.
 */
static const int NSAM_INC = 5;
/**
 * @brief The maximum number of samples.
 * 
 * This variable represents the maximum number of samples and is used in the program for calculations and data storage.
 * The initial value is set to NSAM_INC.
 */
static int maxsam = NSAM_INC;

static const tbx_conf_t tfasta_conf = {TBX_GENERIC, 1, 2, 2, '#', 0};




/**
 * Returns the name of the index format based on the given format code.
 * Mostly used for debugging purposes as we are only using the TFA index format.
 *
 * @param format The format code.
 * @return The name of the format as a string.
 */
static const char* get_format_name(int format);



/**
 * Checks the format of a TFA file.
 *
 * This function takes a filename and a pointer to a `tfa_info_t` structure as input.
 * It checks the format of the TFA file specified by the filename and populates the `tfa_info` structure with relevant information.
 *
 * @param filename The name of the TFA file to check.
 * @param tfa_info A pointer to a `tfa_info_t` structure to store the TFA file information.
 * @return A pointer to a string indicating the format of the TFA file. 
 * Possible values are "TFAv2.0", "TFAv3.0", "TFAv1.0", or NULL if the format is not recognized.
 */
const char *check_tfa_format(const char *filename, tfa_info_t *tfa_info);


const char *check_wtfa_format(const char *filename, tfa_info_t *tfa_info);


/**
 * Initializes a tfasta_file structure with the given tfasta_fname.
 *
 * @param tfasta A pointer to the tfasta_file structure to be initialized.
 * @param tfasta_fname The name of the tfasta file.
 */
int init_tfasta_file(tfasta_file *tfasta, char *tfasta_fname);

/**
 * Closes the specified tfasta file and free memory.
 *
 * @param tfasta A pointer to the tfasta_file structure representing the tfasta file to be closed.
 */
void close_tfasta_file(tfasta_file *tfasta);

int init_wtfasta_file(wtfasta_file *wtfasta, char *wtfasta_fname);


void close_wtfasta_file(wtfasta_file *wtfasta);




/**
 * Reads DNA sequences from a tfasta file within the specified range.
 *
 * @param tfasta The tfasta file to read from.
 * @param chr_name The name of the chromosome.
 * @param init_site The initial site to start reading from.
 * @param end_site The end site to stop reading at.
 * @param n_sam Pointer to the variable that will store the number of samples read.
 * @param n_site Pointer to the variable that will store the number of sites read.
 * @param names Pointer to the array that will store the names of the samples.
 * @param DNA_matr Pointer to the matrix that will store the DNA sequences.
 * @return Returns an integer indicating the success or failure of the operation.
 *        0: Error.
 *        1: Success.
 *       -1 : end of file
 */
int read_tfasta_DNA(
    tfasta_file *tfasta,
    char *chr_name,
    long int init_site,
    long int end_site,
    int *n_sam,
    long long *n_site,
    //char ***names,
    char **DNA_matr    
    );



int read_tfasta_DNA_lite(
    tfasta_file *tfasta,
    char *chr_name,
    long int init_site,
    long int end_site,
    int *n_sam,
    long long *n_site,
    //char ***names,
    char **DNA_matr    
    );

int create_index(
    const char *input_filename,
    int n_threads,
    int min_shift,
    const tbx_conf_t conf);

#ifdef __cplusplus
}
#endif

#endif