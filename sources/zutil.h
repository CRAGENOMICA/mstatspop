/**
 * mstatspop, Statistical Analysis using Multiple Populations for Genomic Data.
 * Created by Sebastian E. Ramos Onsins.
 */
/**
 *  \brief     zutil.h
 *  \details
 *  \author    Joan Jené
 *  \version   1.15
 *  \date      May 22, 2017
 *  \history   - March 29, 2017 : Added stdarg.h include.
 *             - April 7, 2017  : fprintf generates box GZ and Index files.
 *             - April 12, 2017 : different size constants for "in" and "out" buffers & discard Z_BUF_ERROR message.
 *             - April 13, 2017 : Documentation added & some fixes to the fzgetc  & fzeof functions.
 *             - April 18, 2017 : Documentation updated & fzprintf, fzgetc, fzeof, fzclose functions updated, too.
 *             - April 20, 2017 : Library updated for Mac Os X.
 *             - May 19, 2017   : working on zindex::fzseekNearest.
 *             - May 22, 2017   : C++ things (comments, declararions, ...) changed to C things.
 *  \pre
 *  \bug
 *  \warning
 *  \copyright
 *  \usage
 *
 *  	Compilation Options: -lz
 *
 *      Read GZ Example:
 *      ================
 *
 *		FILE *h = 0;
 *		SGZip gz;    														 <---- Define one SGZip structure for every FILE variable.
 *
 * 		h = fzopen("input.tfa.gz", "r", &gz);								 <---- Always pass the SGZip structure to every function that uses a FILE variable.
 *
 *		if (h != NULL) {
 *
 * 			char ch = ' ';
 *			while((!fzeof(h, &gz)) && (ch != '\n') && (ch != '\x0')) {
 *				ch = fzgetc(h, &gz);
 *				printf("%c", ch);
 *			}
 *
 *			fzclose(h, &gz);												 <---- Close your files, always.
 *		}
 *
 *      Write GZ Example:
 *      =================
 *
 *		FILE *h = 0;
 *		SGZip gz;    														 <---- Define one SGZip structure for every FILE variable.
 *
 * 		h = fzopen("input.tfa.gz", "wb+", &gz);								 <---- Always pass the SGZip structure to every function that uses a FILE variable.
 *
 *		if (h != NULL) {
 *
 * 			fzprintf(h, &gz, "Hello World!\n");
 * 			fzprintf(h, &gz, "My name is %s. I am %d years old.\n", name, age);
 * 			fzprintf(h, &gz, data);                                          <------ if the size of the data is larger than MAX_FZPRINTF_MESSAGE do not use '%'
 *
 *			fzclose(h, &gz);												 <---- Close your files, always.
 *		}
 *
 *
 *		Equivalences:
 *		=============
 *
 *          --------    --------
 *			C			ZUtil
 *			--------    --------
 *          fopen    -> fzopen
 *          fclose   -> fzclose
 *          feof     -> fzeof
 *          fprintf  -> fzprintf
 *          fgetc    -> fzgetc
 *          fgets    -> fzgets
 *          getc     -> fzgets
 *
 */

#ifndef ZUTIL_H
#define ZUTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdarg.h>


#ifndef WIN32
    #define F_OK 0
    int access(const char *path, int mode); /* Always use F_OK as mode. Returns 0 if the file exists otherwise it returns -1 */
#endif

#ifdef __cplusplus
extern "C" {
#endif

	typedef int gz_return;			/* It is the return type of a function. It could be some of the following codes: */
	#define GZ_OK 					0		/* The function did it OK */
	#define GZ_UNKNOWN_ERROR 		1		/* Unidentified error. Sorry... */
	#define GZ_ERROR_CLOSE_FILE  	-1		/* The file can not be closed. Disk space? */
	#define GZ_ERROR_CREATE_FILE 	-2		/* The file can not be created. Does the file already exist? Correct Name? Disk space? */
	#define GZ_ERROR_DATA_FILE   	-3		/* The file does not contain the expected data. Maybe expected TXT or GZ file? */
    #define GZ_ERROR_OPEN_FILE		-4		/* The file can not be opened. Does it exist?  */
	#define GZ_EOF                  0      /* The end of file has been reached. It could be OK. Ask for feof() for knowning if the EOF has been reached. */
	#define GZ_PARAMS_ERROR         -6      /* Function parameters are null */

	/* CHUNK is the size of every compressed piece of data */
	#define CHUNK                16384  /* bytes 0x4000 */

	/* The following constants are defined for standard GZ compression */
	#define windowBits        15
	#define ENABLE_ZLIB_GZIP  32
	#define GZIP_ENCODING     16

	/* All zlib functions are called through this macro. It catches and shows errors if they exist */
    /* Z_BUF_ERROR is just an indication that there was nothing for inflate() to do on that call. Simply continue and provide more input data and more output space for the next inflate() call. */
    /* && (status != Z_BUF_ERROR)) */
	#define CALL_ZLIB(x) {                                              \
        int status;                                                     \
        status = x;                                                     \
        if (status < 0)  {                  \
            fprintf (stderr,                                            \
                     "%s:%d: %s returned a bad status of %d.\n",        \
                     __FILE__, __LINE__, #x, status);                   \
        }                                                               \
    }

	#define MAX_FILE_NAME_LEN 2000

	/* This is the struct used for storing gz file information */
    typedef struct {
        char file_name[MAX_FILE_NAME_LEN];
        char index_file_name[MAX_FILE_NAME_LEN];
        struct SGZIndex *index;
		long int compressed_bytes_written; /* Temporary variable used by index when creating the index file */
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
     *
     * @param filename is the file name to be open.
     * @param opentype is "r" for reading a text file; "rb" for reading from a gz file; "wb+" for creating a GZ file, ...
     * @param z is the initialized SGZip structure.
     *
     * @return the file open handle.
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
     * @return the character or GZ_EOF if the end of file is reached.
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


    #define MAX_FZPRINTF_MESSAGE 0x4000 /* 16384 bytes */

    /*
     * Do not used this function directly. This function is used by the fzprintf function.
     *
     * @return: GZ_OK
     * 			GZ_PARAMS_ERROR
     */
    gz_return private_fzprintf(FILE * file_handle, SGZip *z, char *message);

    /**
     * Use this function for writing strings to a file.
     * 
     * It uses an intermediate memory buffer. Only when this buffer is full the
     * data is written into the file.
     *
     * @param file_handle is the open file.
     * @param z is the initialized SGZip structure.
     * @param message is the string to be sent to the file.
     * 
     * For example:
     *
     * ATTENTION! This function won't replace arguments (like %s, %d, ...) if the size of the message is larger than  MAX_FZPRINTF_MESSAGE
     * Usage:
     *     1) fzprintf(handle, gz, data); <---- very long data without %s, %d, ... (no length limit)
     *     2) fzprintf(handle, gz, "%s\n%d"); <---- use it as printf function but do not send messages longer than MAX_FZPRINTF_MESSAGE.
     *
     * @return: GZ_OK
     * 			GZ_PARAMS_ERROR
     */
    gz_return fzprintf(FILE * file_handle, SGZip *z, char *message, ...);

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
     * @return 	GZ_OK
     * 			GZ_ERROR_CREATE_FILE
     */
    gz_return memory_deflate(char *file_name, char *start_address, char *end_address);
    
    /**
     * This function uncompresses a disk file into new disk file (only if it does not exist yet).
     * 
     * @param compressed_file_name is the file name of the compressed input file.
     * @param uncompressed_file_name is the file name of the uncompressed output file.
     * 
     * @return 	GZ_OK
     * 			GZ_ERROR_CREATE_FILE
     * 			GZ_ERROR_OPEN_FILE
     */
    gz_return uncompress_file(const char * compressed_file_name, const char * uncompressed_file_name);
    
    /**
     * This function compresses a file (only if it does not exist yet).
     * 
     * @param uncompressed_file_name is the file name of the uncompressed input file.
     * @param compressed_file_name is the file name of the compressed output file.
     * 
     * @return 	GZ_OK
     * 			GZ_ERROR_CREATE_FILE
     * 			GZ_ERROR_OPEN_FILE
     */
    gz_return compress_file(const char * uncompressed_file_name, const char * compressed_file_name);
    

#ifdef __cplusplus
}
#endif

#endif /* ZUTIL_H */
