/*
 * mstatspop, Statistical Analysis using Multiple Populations for Genomic Data.
 * Created by Sebastian E. Ramos Onsins.
 */
/**
 *  \brief     zutil.h
 *  \details
 *  \author    Joan Jené
 *  \version   1.5
 *  \date      Nov 23, 2016
 *  \pre
 *  \bug
 *  \warning
 *  \copyright
 *  \usage     1) Replace:    FILE *file;
 *
 *                By:         FILE *file;
 *                            SGZip file_gz;
 *
 *             2) Replace:    function(FILE *handle)
 *                            {
 *                               ...
 *                            }
 *
 *                By:         function(FILE *handle, SGZip *handle_gz)
 *                            {
 *                               ...
 *                            }
 *
 *             3) Replace:    fopen    -> fzopen
 *                            fclose   -> fzclose
 *                            feof     -> fzeof
 *                            fprintf  -> fzprintf
 *                            fgetc    -> fzgetc
 *                            fgets    -> fzgets
 *                            getc     -> fzgets
 *
 *             4) Always close files. Use fzclose()
 *
 *             5) If string length is smaller than MAX_FZPRINTF_MESSAGE, you can use '%':
 *
 *                            fzprintf(file, &file_gz, "%s%d\n", string, num);
 *
 *                Else, do not use '%' for the string:
 *
 *                            fzprintf(file, &file_gz, string);
 *                            fzprintf(file, &file_gz, "%d\n", num);
  */

#ifndef ZUTIL_H
#define ZUTIL_H

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif

	/* CHUNK is the size of every compressed piece of data */
	#define CHUNK             0x4000 /* 16384 bytes */

	/* The following constants are defined for standard GZ compression */
	#define windowBits        15
	#define ENABLE_ZLIB_GZIP  32
	#define GZIP_ENCODING     16

	/* All zlib functions are called through this macro. It catches and shows errors if they exist */
	#define CALL_ZLIB(x) {                                                  \
        int status;                                                     \
        status = x;                                                     \
        if (status < 0) {                                               \
            fprintf (stderr,                                            \
                     "%s:%d: %s returned a bad status of %d.\n",        \
                     __FILE__, __LINE__, #x, status);                   \
            exit (EXIT_FAILURE);                                        \
        }                                                               \
    }

	/* This is the struct used for storing gz file information */
    typedef struct {
        unsigned char in[CHUNK];  /* Data to be inflated or deflated */
        unsigned char out[CHUNK]; /* Processed data */
        int bytes_read;
        unsigned have;
        unsigned pointer;
        z_stream strm;
        int first_time;           /* 1:yes / 0:no */
        int file_compressed;      /* 1:yes / 0:no */
        int reading;              /* 1:reading / 0:writing / -1:file closed */
    } SGZip;

    /**
     * This function does the following tasks:
     *   - Initializes the SGZip structure
     *   - Decides if the file is a compressed file or it must be compressed.
     *     (if ".gz" extension then Yes)
     *   - Opens the file.
     */
    FILE * fzopen(const char *filename, const char *opentype, SGZip *z);

    /**
     * Use this function for getting characters from a file.
     * 
     * It uses an intermediate memory buffer. Only when it is empty
     * the function reads again from the file.
     *
     * @param file_handle is the open file.
     * @param z is the initialized SGZip structure.
     * 
     * @return the character or -1 if the end of file is reached.
     */
    int fzgetc(FILE * file_handle, SGZip *z);
    
    /**
     * Use this function for getting count number of characters from a file.
     * 
     * @param row is the output data vector. It must be allocated.
     * @param count is the maximum size of the row.
     * @param file_handle is the open file.
     * @param z is the initialized SGZip structure.
     * 
     * For example:
     *
     * char row[1024];
     * fzgets(row, 1024, file_handle, &file_gz);
     */
    char *fzgets(char *row, int count, FILE * file_handle, SGZip *z);


    /**
     * Use this function for writing strings to a file.
     * 
     * It uses an intermediate memory buffer. Only when this buffer is full the
     * data is written into the file.
     *
     * @param file_handle is the open file.
     * @param message is the string to be sent to the file.
     * @param z is the initialized SGZip structure.
     * 
     * For example:
     *
     * fzprintf(file_handle, &file_gz, "Number of columns = %d .\n", num_columns);
     *
     * @return 0 if OK.
     */
    #define MAX_FZPRINTF_MESSAGE 0x4000 /* 16384 bytes */

    /*
     * Do not used this function directly. This function is used by the fzprintf function.
     */
    int private_fzprintf(FILE * file_handle, SGZip *z, char *message);

    /*
     * ATTENTION! This function won't replace arguments (like %s, %d, ...) if the size of the message is larger than  MAX_FZPRINTF_MESSAGE
     * Usage:
     *     1) fzprintf(handle, gz, data); <---- very long data without %s, %d, ... (no length limit)
     *     2) fzprintf(handle, gz, "%s\n%d"); <---- use it as printf function but do not send messages longer than MAX_FZPRINTF_MESSAGE.
     */
    int fzprintf(FILE * file_handle, SGZip *z, char *message, ...);

    /**
     * This function checks for the end of file.
     * 
     * @param file_handle is the open file.
     * @param z is the initialized SGZip structure.
     * 
     * @return 1 if the end of file is reached. Else 0.
     */
    int fzeof(FILE * file_handle, SGZip *z);
    
    /**
     * This function closes the file.
     * 
     * @param file_handle is the open file.
     * @param z is the initialized SGZip structure.
     */
    void fzclose(FILE * file_handle, SGZip *z);

    /**
     * This function initializes the SGZip structure with default values.
     *
     * @param z is the SGZip structure to be initialized.
     */
    void init_gzip_structure(SGZip *z);
    
    /**
     * This function compresses a block of data in memory to a file.
     * 
     * @param file_name is the output file name.
     * @param start_address is the starting position of the data in memory.
     * @param end_address is the ending position of the data in memory.
     * 
     * @return 0 if OK.
     */
    int memory_deflate(char *file_name, char *start_address, char *end_address);
    
    /**
     * This function uncompresses a disk file into new disk file.
     * 
     * @param compressed_file_name is the file name of the compressed input file.
     * @param uncompressed_file_name is the file name of the uncompressed output file.
     * 
     * @return 0 if OK. 
     */
    int uncompress_file(const char * compressed_file_name, const char * uncompressed_file_name);
    
    /**
     * 
     * @param uncompressed_file_name is the file name of the uncompressed input file.
     * @param compressed_file_name is the file name of the compressed output file.
     * 
     * @return  0 if OK.
     */
    int compress_file(const char * uncompressed_file_name, const char * compressed_file_name);
    

#ifdef __cplusplus
}
#endif

#endif /* ZUTIL_H */
