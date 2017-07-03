/*
 * mstatspop, Statistical Analysis using Multiple Populations for Genomic Data.
 * Created by Sebastian E. Ramos Onsins.
 */
/**
 *  \brief     zindex.h
 *  \details
 *  \author    Joan Jené
 *  \version   1.9.3
 *  \date      June 9, 2017
 *  \history   - April 7, 2007 :  The SGZIndex struct has its own SScanState struct.
 *             - April 10, 2007 : Comments updated & items++ added.
 *             - April 11, 2007 : Search from last position found.
 *             - April 12, 2007 : define ';' removed.
 *             - April 18, 2017 : Documentation updated & fzprintf, fzgetc, fzeof, fzclose functions updated, too.
 *             - April 20, 2017 : Library updated for Mac Os X.
 *             - May 19, 2017   : fzseekNearest added.
 *             - May 22, 2017   : C++ things (comments, declararions, ...) changed to C things.
 *             - May 22, 2017   : fgetc returns 0 (GZ_EOF) when EOF instead of returning the -5 error code.
 *  \pre
 *  \bug
 *  \warning
 *  \copyright
 *  \usage
 *
 *  	1) Compress "data.tfa" and create its index:
 *
 *            compress_file_and_create_index("data.tfa", "data.tfa.gz", "data.tfa.index");
 *
 *      2) GZ Random Access to "data.tfa.gz":
 *
 *            struct SGZIndex idx;
 *
 *            load_index_from_file("data.tfa.index", &idx);
 *
 *            fzseek(handle, &gz, &idx, "000064681", false, DO_NOT_GET_NEAREST);
 *
 *            unload_all_index_positions(&idx);
 *
 *      Example:
 *
 *		FILE *h = 0;
 *		SGZip gz;
 *
 * 		h = fzopen("input.tfa.gz", "r", &gz);
 *
 *		if (h != NULL) {
 *
 *			struct SGZIndex idx;
 *			load_index_from_file(gz.index_file_name, &idx);
 *
 *          long int row_num = -1;                                                                    Set to -1 if you want to search by ID.
 *			if (fzseek(h, &gz, &idx, "00000001", &row_num, false, DO_NOT_GET_NEAREST) == GZ_OK) {     Or set NULL to the ID if you want to seach by position.
 *                                                                                                    After the call to fzseek, the row_num variable has the reached sequence number (0-based).
 *                                                                                                    Set the last parameter to true (1) if you want to search from the last position found. It is faster.
 *	 			char ch = ' ';
 * 				while((!fzeof(h, &gz)) && (ch != '\n') && (ch != '\x0')) {
 *					ch = fzgetc(h, &gz);
 *					printf("%c", ch);
 *				}
 *			}
 *
 *			unload_all_index_positions(&idx);
 *
 *			fzclose(h, &gz);
 *		}
 *
 *      Example of using fzseekNearest     
 *
 *      void test(void) {
 *      	FILE *h = 0;
 *      	SGZip gz;
 *      
 *      	h = fzopen("./examples/input.tfa.gz", "r", &gz);
 *      
 *      	if (h != NULL) {
 *      
 *      		struct SGZIndex idx;
 *      		load_index_from_file(gz.index_file_name, &idx);
 *      
 *              char search_id[10];
 *              strcpy(search_id, "1:2");
 *              printf("Positionate at existing %s\n", search_id);
 *      
 *              long int max = 10;
 *              long int seq_id = 0;
 *      
 *              if (fzseekNearest(h, &gz, &idx, search_id, max, &seq_id) == GZ_OK) {
 *                  printLine(h, &gz, &idx);
 *      
 *                  strcpy(search_id, "1:3");
 *                  printf("Searching for %s\n", search_id);
 *      
 *                  if (fzseekNearest(h, &gz, &idx, search_id, max, &seq_id) == GZ_OK) {
 *      
 *                      printLine(h, &gz, &idx);
 *      
 *                      strcpy(search_id, "1:7");
 *                      printf("Searching for %s\n", search_id);
 *      
 *                      if (fzseekNearest(h, &gz, &idx, search_id, max, &seq_id) == GZ_OK) {
 *      
 *                          printLine(h, &gz, &idx);
 *      
 *                          strcpy(search_id, "1:9");
 *                          printf("Searching for %s\n", search_id);
 *      
 *                          if (fzseekNearest(h, &gz, &idx, search_id, max, &seq_id) == GZ_OK) {
 *      
 *                              printLine(h, &gz, &idx);
 *      		            } else {
 *      			            printf ("Not found %s.\n", search_id);
 *      		            }
 *      		        } else {
 *      			        printf ("Not found %s.\n", search_id);
 *      		        }
 *      		    } else {
 *      			    printf ("Not found %s.\n", search_id);
 *      		    }
 *              }
 *      
 *      		unload_all_index_positions(&idx);
 *      
 *      		fzclose(h, &gz);
 *      	}
 *      }
 *      
 */

#ifndef ZINDEX_H
#define ZINDEX_H

#include "zutil.h"

#undef NodeJS /* Define this constant if you want to use this library in a NodeJS Addon */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NodeJS
	typedef int bool;					/* This is the definition of the boolean type. It could be used stdbool.h instead but */
	#define true 1              		/* this library is only for C99 or upper. */
	#define false 0
#endif

	/* More gz_return codes related wirh Index Files */

	#define GZ_INDEX_KEY_NOT_FOUND  -100    /* The ID (key) searched is not found inside the index file */
	#define GZ_INDEX_FILE_EXISTS    -101    /* The index file exists and can not be created */
	#define GZ_DEFLATED_FILE_EXISTS	-102    /* The compressed file already exists and can not be created */

	#define MAX_ID_LEN 1024				/* This is the maximum length of an ID (The first column of an uncompressed indexable file) */
	                            		/* It should not be longer than 50 characters */

	struct SScanState {					/* The scanner function uses some temporary variables. They are all stored in this structure */
		bool first_char_file;			/* It is true if the scanner starts scanning and it is the first char of the file */
		bool EOL_found;					/* It is set to true if the scanner founds an end of line */
		bool discarding_comments;		/* It is set to true if the scanner is discarding comments from the stream (those that start with '#') */
		bool getting_id;				/* It is set to true if the scanner is getting an ID. IDs start with a new line and end with a tabulator character */
		int id_byte;					/* While getting the ID, this variable points to every char of the ID word */
		int id_start_byte;				/* It is the position of the ID in an uncompressed block of data */
		long int id_start_block;        /* It is the starting position of the block where is the current ID */
		char id[MAX_ID_LEN];			/* It is the current processed ID */
	};

	/**
	 * This functions initializes the values of the structure
	 *
	 * @param ss is the structure to be initialized.
	 */
	void init_scanstate_structure(struct SScanState *ss);

	/* Index File Example on Disk */
	/* 0 position key offset key offset key offset key offset ... 0 position key offset key offset key offset key offset ... */
	/* 0 indicates that the following value is the position of all next keys */

	struct SIndexPosDisk {       		/* This is an entry of the index file. This is the structure that is stored on disk */
		unsigned long key;      		/* It is the hash(ID) */
		long int pos_offset;			/* It is the position in bytes of the compressed file for reaching the key */
	};

	/* Index File Example in Memory */
	/* [key position offset] -> [key position offset] -> ... */

	struct SIndexPos {                 	/* This is an entry of the index file in memory */
		unsigned long key;      		/* It is the hash(ID) */
		long int position;				/* It is the position in bytes of the compressed file for reaching the key */
		long int offset;				/* It is the offset in bytes of the uncompressed block of data for reaching the key */
        long int row_num;               /* It is a 0-based sequence count. From 0 to n. */
		struct SIndexPos *next;			/* This is the pointer to the next SIndexPos inside a list of these objects */
	};

	struct SGZIndex	{					/* This is the IndexFile object */
		struct SIndexPos *first;		/* It is the pointer to the first index entry */
		struct SIndexPos *last;			/* It is the pointer to the last index entry */
		long int items;					/* It has the length of the list */
		struct SIndexPos *last_search;  /* It is a pointer to the last searched position in order to perform next searches from the last position found */
		long int last_search_count;     /* The same as last_search but for position */
		struct SScanState state;
	};

	/**
	 * This functions initializes the values of the structure
	 *
	 * @param idx is the structure to be initialized.
	 */
	void init_gzindex_structure(struct SGZIndex *idx);

    /**
     * This is a private function. Do not call it.
     * It scans every received character and it creates index entries.
     *
     * @param ch is the received current character.
     * @param i is the position of the current character inside the uncompressed block of data.
     * @param state is an structure with some temporary variables. See SScanState comments for more information.
     * @param idx is the index file object that is being created in every call to this function.
     * @param gz_block_starting_pos is the starting position for all IDs found in the call of this function.
     *
     * @return 	GZ_OK
     * 			GZ_PARAMS_ERROR
     */
	gz_return private_scan_deflated_for_create_index_positions(char ch, int i, struct SGZIndex *idx, long int gz_block_starting_pos);

    /**
     * This function compresses a file and it creates its index.
     *
     * @param uncompressed_file_name is the name of the file to be compressed.
     * @param compressed_file_name is the name of the new compressed file.
     * @param index_file_name is the name of the new index file.
     *
     * @return  GZ_OK
     * 			GZ_PARAMS_ERROR
     * 			GZ_ERROR_OPEN_FILE
     * 			GZ_INDEX_FILE_EXISTS
     * 			GZ_DEFLATED_FILE_EXISTS
     */
    gz_return compress_file_and_create_index(const char *uncompressed_file_name, const char *compressed_file_name, const char *index_file_name);

    /**
     * Use this function for random access to a desired sequence.
     * You can access by setting the sequence ID (first column) or by setting the sequence number (row count).
     * It only works with compressed files with index file.
     *
     * @param file_handle is the compressed file handle.
     * @param z is the GZ structure of the compressed file.
     * @param idx is the index file structure.
     * @search_id is the ID to be found.
     * @row_num is the sequence number to be reached (0-based). It is also a return parameter that indicates the reached sequence number.
     * @param from_last_search is 1 if the search must be done from the last found position.
     * @param get_nearest is DO_NOT_GET_NEAREST if False and GET_R_NEAREST if get nearest from right.
     *
     * @return 	GZ_OK
     * 			GZ_PARAMS_ERROR
     * 			GZ_ERROR_DATA_FILE
     * 			GZ_INDEX_KEY_NOT_FOUND
     */
    gz_return fzseek(FILE *file_handle, SGZip *z, struct SGZIndex *idx, const char *search_id, long int *row_num, int from_last_search);


    /**
     * Use this function for random access to a desired sequence ID or the next available ID of the same chromosome.
     * You can access by setting the sequence ID (first column).
     * It only works with compressed files with index file.
     *
     * @param file_handle is the compressed file handle.
     * @param z is the GZ structure of the compressed file.
     * @param idx is the index file structure.
     * @param search_id is the ID to be found.
     * @param max is maximum sequence id to search in. For example, imagine that you are looking for "1:3" but it does not exist. So, the function will try "1:4", "1:5", ... until this max..
     * @param seq_id_found If the function has to search for the nearest sequence inside the same chromosome. Otherwise, -1.
     *        If the exact search_id is not found but it is found a nearest sequence within the same chromosome the value of seq_id_found will be the id of the found sequence (and GZ_OK). 
     *        For example: - you are searching for "1:2", It exists, then the file pointer will point to the start of "1:2". The function will return GZ_OK and this parameter will be -1.
     *                     - you are searching for "1:3", It does not exist, then the file pointer will point to the nearest sequence inside the same chromosome "1:6". The function will return GZ_OK and this parameter will be 6.
     *                     - you are searching for "1:9", It does not exist and it is the last sequence of this chromosome, so it does not have a nearest sequence. Then, the file pointer won't be moved. The function will return GZ_INDEX_KEY_NOT_FOUND and this parameter will be -1.
     *
     * @return 	GZ_OK
     *          GZ_PARAMS_ERROR
     * 			GZ_ERROR_DATA_FILE
     *          GZ_INDEX_KEY_NOT_FOUND
     */
    gz_return fzseekNearest(FILE *file_handle, SGZip *z, struct SGZIndex *idx, const char *search_id, long int max, long int *seq_id_found);

	/**
	 * Dan Bernstein djb2 has function.
	 * This function gets an string and it generated an unique hash value.
	 *
	 * @param str is the string to be converted into a key value.
	 *
	 * @return the hash.
	 */
	unsigned long hash(unsigned char *str);

	/**
	 * This functions creates a new index entry.
	 *
	 * @param idx is the index file structure.
	 *
	 * @return	GZ_OK
	 * 			GZ_PARAMS_ERROR
	 */
	gz_return add_index_position(struct SGZIndex *idx);

	/**
	 * This function gets the ID position inside a compressed file.
	 *
	 * @param search is the ID to be found (first column of the file).
	 * @param idx is the index file structure. The function returns the index entry position of the searched ID, here.
	 * @param from_last_search is 1 if the search must be done from the last found position.
	 *
	 * @return 	GZ_OK
	 * 			GZ_PARAMS_ERROR
	 */
    struct SIndexPos *get_index_position_by_id(const char *search, struct SGZIndex *idx, int from_last_search);

    /**
	 * This function gets the ID position inside a compressed file.
	 *
	 * @param pos is the sequence number (row count).
	 * @param idx is the index file structure. The function returns the index entry position of the searched ID, here.
	 * @param from_last_search is 1 if the search must be done from the last found position.
	 *
	 * @return 	GZ_OK
	 * 			GZ_PARAMS_ERROR
	 */
    struct SIndexPos *get_index_position_by_pos(long int pos, struct SGZIndex *idx, int from_last_search);

	/**
	 * This function load all index entries from file.
	 *
	 * @param file_name is the name of the index file.
	 * @param idx is the index file structure.
	 *
	 * @return 	GZ_OK
	 * 			GZ_PARAMS_ERROR
	 * 			GZ_ERROR_OPEN_FILE
	 */
    gz_return load_index_from_file(const char *file_name, struct SGZIndex *idx);

	/**
	 * This function save all index entries to a file.
	 *
	 * @param file_name is the name of the index file.
	 * @param idx is the index file structure.
	 *
	 * @return 	GZ_OK
	 * 			GZ_PARAMS_ERROR
	 * 			GZ_ERROR_CREATE_FILE
	 */
    gz_return save_index_to_file(const char *file_name, struct SGZIndex *idx);

	/**
	 * This function removes all index entries from memory.
	 *
	 * @param idx is the index file structure.
	 *
	 * @return 	GZ_OK
	 * 			GZ_PARAMS_ERROR
	 */
    gz_return unload_all_index_positions(struct SGZIndex *idx);


#ifdef __cplusplus
}
#endif

#endif /* ZINDEX_H */
