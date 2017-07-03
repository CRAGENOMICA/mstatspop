/**
 * mstatspop, Statistical Analysis using Multiple Populations for Genomic Data.
 * Created by Sebastian E. Ramos Onsins.
 */
/**
 *  \brief     zutil.c
 *  \details
 *  \author    Joan Jené
 *  \version   1.15
 *  \date      June 22, 2017
 *  \pre
 *  \bug
 *  \warning
 *  \copyright
 */

#include "zutil.h"

#ifndef WIN32
    #define F_OK 0
    int access(const char *path, int mode) {
        int ret = 0;

        FILE *handle = fopen(path, "r");
        if (handle != NULL) {
            fclose(handle);
        } else {
            ret = -1;
        }

        return ret;
    }
#endif

FILE * fzopen(const char *filename, const char *opentype, SGZip *z) {
	init_gzip_structure(z);

	strcpy((*z).file_name, filename);

	(*z).file_compressed = 0;
	if (strlen(filename) > 3) {
		if ((filename[strlen(filename) - 3] == '.')
				&& (filename[strlen(filename) - 2] == 'g')
				&& (filename[strlen(filename) - 1] == 'z')) {
			(*z).file_compressed = 1;

		}
	}

	/* If the file is not compressed, it has not an index file but the name should be the following one in both cases: */

	if ((*z).file_compressed == 1) {
		strncpy((*z).index_file_name, filename, strlen(filename) - 3); /* remove last ".gz" */
	}

	strcat((*z).index_file_name, ".index");

	return fopen(filename, opentype);
}

void init_gzip_structure(SGZip *z) {
	(*z).bytes_read = 0;
	(*z).have = 0;
	(*z).pointer = 0;
	(*z).first_time = 1;
	(*z).file_compressed = 1;
	(*z).index = NULL;
	(*z).strm.avail_out = CHUNK;
	(*z).reading = 1;
	(*z).strm.next_in = 0;
	(*z).strm.avail_in = 0;
	(*z).strm.total_in = 0;
	(*z).strm.next_out = 0;
	(*z).strm.avail_out = 0;
	(*z).strm.total_out = 0;
	(*z).strm.msg = 0;
	(*z).strm.state = 0;
	(*z).strm.zalloc = 0;
	(*z).strm.zfree = 0;
	(*z).strm.opaque = 0;
	(*z).strm.data_type = 0;
	(*z).strm.adler = 0;
	(*z).strm.reserved = 0;
	(*z).compressed_bytes_written = 0;
	memset((*z).file_name, '\x0', MAX_FILE_NAME_LEN);
	memset((*z).index_file_name, '\x0', MAX_FILE_NAME_LEN);
}

/*
 - initialize data (Only first call)
 - read data
 - uncompress a piece of data into a buffer
 - return one char of the buffer
 Next calls to this funcion:
 - If there are chars inside the buffer, return one
 - Else uncompress another piece of data
 - But if there are not more piece of data to uncompress then read from file
 */
/* return unsigned char = int */
long int num_chars = 0;

int fzgetc(FILE * file_handle, SGZip *z) {
	int ret = 0;

	(*z).reading = 1;
	if ((*z).file_compressed == 0) {
		/* The file is not compressed */
		ret = fgetc(file_handle);
		num_chars++;
	} else {
		/* The file is compressed */
		if ((*z).first_time == 1) {
			/* This is the first read of the input file. Let's initialize some variables */
			(*z).strm.zalloc = Z_NULL;
			(*z).strm.zfree = Z_NULL;
			(*z).strm.opaque = Z_NULL;
			(*z).strm.next_in = (*z).in;
			(*z).strm.avail_in = 0;

			CALL_ZLIB(inflateInit2(&((*z).strm), windowBits | ENABLE_ZLIB_GZIP));
		}

		if (((*z).have > 0) && ((*z).pointer < (*z).have)) {
			/* The buffer has uncompressed chars */
		} else {
			if (((*z).strm.avail_in > 0) && ((*z).first_time == 0)) { /*!!<--- Before avail_out == 0 */
				/* The buffer of uncompressed chars is empty but
				 It exists read data that can be uncompressed */
			} else {
				/* The buffer of uncompressed chars is empty and
				 It does not exist read data that could be uncompressed, so
				 Read data from the compressed input file */
				(*z).strm.next_in = (*z).in;
				(*z).strm.avail_in = 0; /*!!<--- New*/
				(*z).bytes_read = fread((*z).in, sizeof(unsigned char),
						sizeof((*z).in), file_handle);
				(*z).strm.avail_in = (*z).bytes_read;
			}

			/*	if ( ((*z).strm.avail_in == 0) && There is no data to be uncompressed
						     ((*z).strm.avail_out == 0)) There is no data uncompressed
						{
							Let's read again
							(*z).strm.next_in = (*z).in;
							(*z).bytes_read = fread((*z).in, sizeof(unsigned char),
						sizeof((*z).in), file_handle);
				(*z).strm.avail_in = (*z).bytes_read;
			
			}*/

			if ((*z).bytes_read > 0) {
				/* Uncompress the data */
				(*z).strm.avail_out = CHUNK;
				(*z).strm.next_out = (*z).out;
				(*z).have = 0; /*!!<--- New*/
                CALL_ZLIB(inflate(&((*z).strm), Z_NO_FLUSH));
				(*z).have = CHUNK - (*z).strm.avail_out;
				(*z).first_time = 0;

				num_chars += (*z).have;
			}

			/* The next char to be returned will be the first one in the buffer */
			(*z).pointer = 0;
		}

		if ((*z).bytes_read == 0) {
			inflateEnd(&((*z).strm));
			ret = GZ_EOF;
		} else {
            /*int aux=99;*/
            if ((*z).have > 0) { /* If the last inflated compressed block has uncompressed chars, then: */
                /* Return the current char */
                ret = (*z).out[(*z).pointer];
                /* Increment the pointer for the next call to this function */
                (*z).pointer++;
            } else {
                /* the last inflated uncompressed block has 0 uncompressed chars */
                ret = GZ_EOF; /*assumed is the last block of uncompressed data*/
                (*z).strm.avail_in = 0;
                /*aux = fzeof(file_handle, z);*/
            }
                
			if ((*z).pointer >= (*z).have) {
				(*z).pointer = 0;
				(*z).strm.avail_out = CHUNK;
				(*z).have = 0;
			}
		}
	}

	return ret;
}

char *fzgets(char *row, int count, FILE * file_handle, SGZip *z) {
	char *ret = row;
	int i = 0;

	(*z).reading = 1;
	if ((*z).file_compressed == 0) {
		ret = fgets(row, count, file_handle);
	} else {
		while ((i < count) && (!fzeof(file_handle, z))) {
			row[i] = fzgetc(file_handle, z);
			i++;
		}
	}

	return ret;
}

gz_return private_scan_deflated_for_create_index_positions(char ch, int i, struct SGZIndex *idx, long int gz_block_starting_pos);

/*
 * Do not used this function directly. This function is used by the fzprintf function.
 */
gz_return private_fzprintf(FILE * file_handle, SGZip *z, char *message) {
	gz_return ret = GZ_OK;
	unsigned int i = 0;
	int last_i = 0;

	if ((file_handle == NULL) || (z == NULL) || (message == NULL)) {
		ret = GZ_PARAMS_ERROR;
	} else {
		(*z).reading = 0;

		if ((*z).first_time == 1) {
			/* This is the first write of the output file. Let's initialize some variables */
			(*z).first_time = 0;
			(*z).strm.zalloc = Z_NULL;
			(*z).strm.zfree = Z_NULL;
			(*z).strm.opaque = Z_NULL;
			memset((*z).in, '\x0', CHUNK);
			(*z).strm.next_in = (*z).in;
			(*z).strm.avail_in = 0;

			CALL_ZLIB(deflateInit2(&((*z).strm), Z_DEFAULT_COMPRESSION, Z_DEFLATED, windowBits | GZIP_ENCODING, 8, Z_DEFAULT_STRATEGY));

			/*!! Removed: (*z).have = CHUNK;*/ /* the next_in buffer has all bytes free */
		}

		if ((*z).strm.avail_in < CHUNK) {

			/* There is still space inside the "(*z).in" buffer of uncompressed chars */
			for (i = 0;
					((i < strlen(message)) && ((*z).strm.avail_in + i < CHUNK));
					i++) {
				(*z).in[(*z).strm.avail_in + i] = message[i];
			}
			(*z).strm.avail_in += i;
			/* if (i > 0) { last_i = i - 1; } else { last_i = 0; } */
			last_i = i;

			if ((*z).strm.avail_in < CHUNK) {
				/* There is still space inside the "(*z).in" buffer of uncompressed chars. Deflate them the next time... */
			} else {
				if (z->index == NULL) {
					/* ************************************************************************* */
					/* * This is the fzprintf WITHOUT Index Version                            * */
					/* * It is faster and it generates smaller GZ files.                       * */
					/* ************************************************************************* */
					/* Deflate now */
					do {
						(*z).strm.avail_out = CHUNK;
						(*z).strm.next_out = (*z).out;
						CALL_ZLIB(deflate(&((*z).strm), Z_NO_FLUSH));
						(*z).have = CHUNK - (*z).strm.avail_out;
						fwrite((*z).out, sizeof(char), (*z).have, file_handle);
					} while ((*z).strm.avail_out == 0);

					/* Reset */
					memset((*z).in, '\x0', CHUNK);
					(*z).strm.next_in = (*z).in;
					(*z).strm.avail_in = 0;


				} else {
					/* ************************************************************************* */
					/* * This is the fzprintf WITH Index Version                               * */
					/* * It generates Index files and it enables applications to access        * */
					/* * GZ files randomly (faster)                                            * */
					/* ************************************************************************* */

					/* int isEOF = feof(file_handle); */
					int flush = Z_NO_FLUSH;
					long int bytes_write = 0;
					long int gz_block_starting_pos = 0;

					/* Scan the data for IDs */
					gz_block_starting_pos = (*z).compressed_bytes_written;

					for (i = 0; i < (*z).strm.avail_in; i++) {
						private_scan_deflated_for_create_index_positions((*z).in[i], i, (*z).index, gz_block_starting_pos);
					}

					/* Set the initial mark */
					/*
					 * If flush is set to Z_FULL_FLUSH, all output is flushed as with Z_SYNC_FLUSH, and the compression state is reset
					 * so that decompression can restart from this point if previous compressed data has been damaged or if random
					 * access is desired. Using Z_FULL_FLUSH too often can seriously degrade compression.
					 */

					/* THIS LOOP WILL ONLY ITERATE ONCE: (HERE ONLY FOR SECURITY PURPOSES) */
					do {
						/* compress the data */

						(*z).strm.next_out = (*z).out;
						(*z).strm.avail_out = CHUNK;

/*
						if (isEOF == 0) {
							flush = Z_FULL_FLUSH;
						}
						else {
							flush = Z_FINISH;
						}
*/
						flush = Z_FULL_FLUSH;

						CALL_ZLIB(deflate(&((*z).strm), flush));

						(*z).have = CHUNK - (*z).strm.avail_out;

						/* write compressed data to the output file */

						bytes_write = fwrite((*z).out, sizeof(char), (*z).have, file_handle);
						(*z).compressed_bytes_written += bytes_write;

					} while ((*z).strm.avail_in != 0);
					/*!! Removed } while ((*z).strm.avail_out == 0);*/ /* continue compressing until it does not exist more data to be compressed */


					/* Reset */
					memset((*z).in, '\x0', CHUNK);
					(*z).strm.next_in = (*z).in;
					(*z).strm.avail_in = 0;
				}

				/* if there are more chars in the current message, then insert them in the z.in buffer: */
				/* how can it be possible? because of this: (for example)
				 - previous message size : 1 byte
				 - current message size  : CHUCK bytes
				 - total space in (*z).strm.in is full because it is of CHUNK size but still remains 1 byte to be written into the disk.
				 So,  this byte is stored in the empty (*z).in for the next call.
				 */
				for (i = last_i;
						((i < strlen(message))
								&& ((*z).strm.avail_in + i - last_i < CHUNK));
						i++) {
					(*z).in[(*z).strm.avail_in] = message[i];
					(*z).strm.avail_in++;
				}
			}
		}



	}

	return ret;
}

/*http://stackoverflow.com/questions/1056411/how-to-pass-variable-number-of-arguments-to-printf-sprintf*/
/*
 * message -> buffer
 * buffer append into z.in
 * if buffer is full:
 *   deflate z.in into z.out
 *   write z.out
 */
gz_return fzprintf(FILE * file_handle, SGZip *z, char *message, ...) {
	gz_return ret = GZ_OK;

	if ((file_handle == NULL) || (z == NULL) || (message == NULL)) {
		ret = GZ_PARAMS_ERROR;
	} else {
		char *buffer = NULL;
		size_t message_len = 0;
		size_t sent_chars = 0;
		int replace_args = 0;
		char buffer_block[CHUNK + 1]; /* +1 for the \x0 */
    	va_list args;
        int rest = 0;

		message_len = strlen(message);
		replace_args = (message_len < MAX_FZPRINTF_MESSAGE);

		if (replace_args) {
			buffer = (char *) malloc(2 * MAX_FZPRINTF_MESSAGE * sizeof(char));
			va_start(args, message);
			vsnprintf(buffer, 2 * MAX_FZPRINTF_MESSAGE, message, args);
			va_end(args);
		} else {
			buffer = message;
		}

		if ((z != 0) && ((*z).file_compressed == 0)) {
			/* The file is not compressed */
			ret = fprintf(file_handle, "%s", buffer);
		} else {
			/* If the message is smaller than CHUNK then no problem */
			message_len = strlen(buffer);

			if (message_len < CHUNK) {
				ret += private_fzprintf(file_handle, z, buffer);
			} else {
				/* But if the message is longer than CHUNK then it must be divided into small blocks of CHUNK size */

				/* Send all blocks of CHUNK size */
				sent_chars = 0;
				while (message_len - sent_chars >= CHUNK) {
					memcpy(buffer_block, buffer + sent_chars, CHUNK);
					buffer_block[CHUNK] = '\x0';
					ret += private_fzprintf(file_handle, z, buffer_block);
					sent_chars += CHUNK;
				}

				/* Send the last block with size less than CHUNK */
				rest = message_len - sent_chars;
				if (rest < CHUNK) {
					memcpy(buffer_block, buffer + sent_chars, rest);
					buffer_block[rest] = '\x0';
					ret += private_fzprintf(file_handle, z, buffer_block);
					sent_chars += rest;
				}

				/*
				 memset(buffer_block, '\x0', CHUNK);

				 sent_chars = 0;
				 while (sent_chars < message_len) {
				 memcpy(buffer_block, buffer + sent_chars, CHUNK);
				 sent_chars += CHUNK;
				 ret += private_fzprintf(file_handle, z, buffer_block);
				 }*/
			}
		}

		if ((replace_args) && (buffer != NULL)) {
			free(buffer);
			buffer = NULL;
		}
	}

	return ret;
}

int fzeof(FILE * file_handle, SGZip *z) {
	int ret = feof(file_handle);

	if ((*z).file_compressed == 1) {

		ret = ((ret > 0) &&							/* End of file reached, AND */
			  ((*z).strm.avail_in == 0) &&			/* It does not exist compressed data that can be uncompressed, AND */
			  ((*z).pointer >= (*z).have));

		if (ret) { 					/* It does not exist available uncompressed data */
			 /*if ((ret = ((ret) &&  Reached the end of the file and
			  ((*z).strm.avail_out != 0) &&  there is no more input data to uncompress and
			  ((*z).pointer >= (*z).have))) == 1) { there is no more chars in the buffer*/

			inflateEnd(&((*z).strm));

			/* Initialize first time to true */
			(*z).first_time = 1;
		}
	}

	return ret;
}

gz_return save_index_to_file(const char *file_name, struct SGZIndex *idx);
gz_return unload_all_index_positions(struct SGZIndex *idx);

void fzclose(FILE * file_handle, SGZip *z) {
	if ((file_handle != 0) && (z != 0) && ((*z).reading != -1)) {
		if (((*z).file_compressed == 0) || ((*z).reading == 1)) {
		} else {
			if (((*z).file_compressed == 1) && ((*z).reading == 0)) {

				if ((*z).index != NULL) {
					/* ************************************************************************* */
					/* * This is the fzclose WITH Index Version                                * */
					/* * It generates Index files and it enables applications to access        * */
					/* * GZ files randomly (faster)                                            * */
					/* * Let's create the index positions for this last piece of output data   * */
					/* ************************************************************************* */

					long int gz_block_starting_pos = 0;
					long int i = 0;

					/* Scan the data for IDs */
					gz_block_starting_pos = (*z).compressed_bytes_written;

					for (i = 0; i < (*z).strm.avail_in; i++) {
						private_scan_deflated_for_create_index_positions((*z).in[i], i, (*z).index, gz_block_starting_pos);
					}
				}

				/* Deflate the z.in buffer into the z.out buffer (if there are chars) */
				do {
					(*z).strm.avail_out = CHUNK;
					(*z).strm.next_out = (*z).out;
					CALL_ZLIB(deflate(&((*z).strm), Z_FINISH));
					(*z).have = CHUNK - (*z).strm.avail_out;
					fwrite((*z).out, sizeof(char), (*z).have, file_handle);
				} while ((*z).strm.avail_out == 0);

				CALL_ZLIB(deflateEnd(&((*z).strm)));
			}
		}
		fclose(file_handle);
		file_handle = 0;
		(*z).reading = -1; /* file closed */
	}

	if ((*z).index != NULL) {
		/* ************************************************************************* */
		/* * This is the fzclose WITH Index Version                               * */
		/* * It generates Index files and it enables applications to access        * */
		/* * GZ files randomly (faster)                                            * */
		/* ************************************************************************* */

		save_index_to_file((*z).index_file_name, (*z).index);
		unload_all_index_positions((*z).index);
	}
}

gz_return memory_deflate(char *file_name, char *start_address, char *end_address) {
	gz_return ret = GZ_OK;

	if ((file_name == NULL) || (start_address == NULL)
			|| (end_address == NULL)) {
		ret = GZ_PARAMS_ERROR;
	} else {
		FILE *file_handle = NULL;

		if ((file_handle = fopen(file_name, "wb+")) != NULL) {
			SGZip z;
			z.first_time = 0;
			z.strm.zalloc = Z_NULL;
			z.strm.zfree = Z_NULL;
			z.strm.opaque = Z_NULL;
			CALL_ZLIB(deflateInit2(&(z.strm), Z_DEFAULT_COMPRESSION, Z_DEFLATED, windowBits | GZIP_ENCODING, 8, Z_DEFAULT_STRATEGY));

			z.strm.next_in = (unsigned char *) start_address;
			z.strm.avail_in = end_address - start_address;
			do {
				z.strm.avail_out = CHUNK;
				z.strm.next_out = z.out;
				CALL_ZLIB(deflate(&(z.strm), Z_FINISH));
				z.have = CHUNK - z.strm.avail_out;
				fwrite(z.out, sizeof(char), z.have, file_handle);
			} while (z.strm.avail_out == 0);
		} else {
			ret = GZ_ERROR_CREATE_FILE;
		}
	}

	return ret;
}

gz_return uncompress_file(const char * compressed_file_name, const char * uncompressed_file_name) {
	gz_return ret = GZ_OK;

	if ((compressed_file_name == NULL) || (uncompressed_file_name == NULL)) {
		ret = GZ_PARAMS_ERROR;
	} else {

		if (access(uncompressed_file_name, F_OK) != -1) {
			/* uncompressed file exists */
		} else {
			/* uncompressed file doesn't exist */

			gzFile compressed_file;
			FILE * uncompressed_file;

			unsigned char buffer[CHUNK];

			compressed_file = gzopen(compressed_file_name, "rb");
			if (compressed_file != NULL) {
				uncompressed_file = fopen(uncompressed_file_name, "wb+");
				if (uncompressed_file != NULL) {
					/* ................................................................. */
					int bytes_read = 0;

					while (!gzeof(compressed_file)) {
						memset(buffer, '\x0', sizeof(buffer));
						bytes_read = gzread(compressed_file, buffer, CHUNK - 1);
						fwrite(buffer, sizeof(unsigned char), bytes_read,
								uncompressed_file);
					}
					/* ................................................................. */
					fclose(uncompressed_file);
				} else {
					ret = GZ_ERROR_CREATE_FILE;
				}

				gzclose(compressed_file);
			} else {
				ret = GZ_ERROR_OPEN_FILE;
			}
		}
	}
	return ret;
}

gz_return compress_file(const char * uncompressed_file_name, const char * compressed_file_name) {
	gz_return ret = GZ_OK;

	if ((uncompressed_file_name == NULL) || (compressed_file_name == NULL)) {
		ret = GZ_PARAMS_ERROR;
	} else {
		if (access(compressed_file_name, F_OK) != -1) {
			/* compressed file exists */
		} else {
			/* compressed file doesn't exist */
			FILE * uncompressed_file;
			gzFile compressed_file;

			char buffer[CHUNK];

			uncompressed_file = fopen(uncompressed_file_name, "rb");
			if (uncompressed_file != NULL) {
				compressed_file = gzopen(compressed_file_name, "wb");
				if (compressed_file != NULL) {
					/* ................................................................. */
					while (!feof(uncompressed_file)) {
						unsigned int bytes_read = 0;
						memset(buffer, '\x0', sizeof(buffer));
						bytes_read = fread(buffer, sizeof(char), sizeof(buffer),
								uncompressed_file);

						if (bytes_read > 0) {
							gzwrite(compressed_file, buffer, bytes_read); /*strlen(buffer));*/
						}
					}
					/* ................................................................. */
					gzclose(compressed_file);
				} else {
					ret = GZ_ERROR_CREATE_FILE;
				}

				fclose(uncompressed_file);
			} else {
				ret = GZ_ERROR_OPEN_FILE;
			}
		}
	}

	return ret;
}
