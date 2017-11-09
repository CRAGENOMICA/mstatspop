/**
 * mstatspop, Statistical Analysis using Multiple Populations for Genomic Data.
 * Created by Sebastian E. Ramos Onsins.
 */
 /**
 *  \brief     zindex.c
 *  \details
 *  \author    Joan Jené
 *  \version   1.9.3
 *  \date      May 22, 2017
 *  \pre
 *  \bug
 *  \warning
 *  \copyright
 *  \usage
 */

#include "zindex.h"

void init_gzindex_structure(struct SGZIndex *idx) {
	idx->first = NULL;
	idx->last = NULL;
	idx->items = 0;
	idx->last_search = NULL;
	idx->last_search_count = 0;
	init_scanstate_structure(&idx->state);
}

gz_return private_scan_deflated_for_create_index_positions(char ch, int i, struct SGZIndex *idx, long int gz_block_starting_pos) {
	gz_return ret = GZ_OK;

	if (idx == NULL) {
		ret = GZ_PARAMS_ERROR;
	} else {

		/* Decide States */

		if (((idx->state.first_char_file == true) || (idx->state.EOL_found == true))
				&& (ch == '#')) {
			/* If the line starts with '#' it means that its a comment. So, the whole line is discarded.
			 * But the line could be long enough to be in more than one CHUNK. So, we set the discarding
			 * state variable to true.
			 */
			idx->state.discarding_comments = true;
		}

		if ((ch != '#')
				&& ((idx->state.first_char_file == true) || (idx->state.EOL_found == true))) {
			/* If it is the first char of the file and it does not start with '#' or
			 * It is the start of a new line and it does not start with '#', then it is the start of an Id
			 */
			idx->state.getting_id = true;
			idx->state.id_byte = 0;
		}

		/* Do States */

		if (idx->state.discarding_comments == true) {
			/* Discard the entire comment line until an End of Line is reached */
			if (ch == '\n') {
				idx->state.discarding_comments = false;
			}
		} else {
			if (idx->state.getting_id == true) {
				if (ch == '\t') {
					/* Getting Id until TAB is reached */
					idx->state.getting_id = false;

					/* Create index position */
					add_index_position(idx);
					idx->last->key = hash((unsigned char *)idx->state.id);
					idx->last->position = idx->state.id_start_block; /*gz_block_starting_pos;*/
					idx->last->offset = idx->state.id_start_byte;

				} else {
					if (idx->state.id_byte == 0) {
						/* here is where the ID begins */
						memset(idx->state.id, '\x0', MAX_ID_LEN);
						idx->state.id_start_byte = i;
						idx->state.id_start_block = gz_block_starting_pos;
					}
					if (idx->state.id_byte < MAX_ID_LEN) {
						/*  Getting Id without exceeding the size of the id array */
						idx->state.id[idx->state.id_byte] = ch;
						idx->state.id_byte++;
					}
				}
			}
		}

		idx->state.EOL_found = (ch == '\n') ? true : false;
		idx->state.first_char_file = false;
	}

	return ret;
}

void init_scanstate_structure(struct SScanState *ss) {
	if (ss != NULL) {
		ss->first_char_file = true;
		ss->EOL_found = false;
		ss->discarding_comments = false;
		ss->getting_id = false;
		ss->id_byte = 0;
		ss->id_start_byte = 0;
		ss->id_start_block = 0;
		memset(ss->id, '\x0', MAX_ID_LEN);
	}
}

gz_return compress_file_and_create_index(const char *uncompressed_file_name, const char *compressed_file_name, const char *index_file_name) {
	gz_return ret = GZ_OK;

	unsigned int i = 0;
	FILE *uncompressed_handle = NULL;
    FILE *compressed_handle = NULL;
	long int compressed_bytes_written = 0;
	long int gz_block_starting_pos = 0;
	long int bytes_write = 0;
	int isEOF = false;
	int flush = Z_FULL_FLUSH;
	struct SGZIndex idx;

	if ((uncompressed_file_name == NULL) || (compressed_file_name == NULL)
			|| (index_file_name == NULL)) {
		ret = GZ_PARAMS_ERROR;
	} else {
		init_gzindex_structure(&idx);


		/* Files check */

		if (access((char *)uncompressed_file_name, F_OK) == -1) {
			/*printf("The file '%s' does not exist.", uncompressed_file_name);*/
			return GZ_ERROR_OPEN_FILE;
		}

		if (access((char *)index_file_name, F_OK) != -1) {
			/*printf("The index file '%s' already exists. If you want to update it, remove it first.", index_file_name);*/
			return GZ_INDEX_FILE_EXISTS;
		}

		if (access((char *)compressed_file_name, F_OK) != -1) {
			/*printf("The compressed file '%s' already exists. If you want to update it, remove it first.", compressed_file_name);*/
			return GZ_DEFLATED_FILE_EXISTS;
		}

		/* Process */

		/* Loop the uncompressed file
		 * - Get X Kb of uncompressed input data
		 * - Create an index position for every ID:
		 *   - Use the hash(ID) as key
		 *   - Use the starting position of the current GZ block as the position
		 *   - Use the starting ID position as the offset
		 * - Compress the data
		 *
		 * Save the index file
		 */

        uncompressed_handle = fopen((char *)uncompressed_file_name, "r");

		if (uncompressed_handle != NULL) {

			compressed_handle = fopen((char *)compressed_file_name, "wb+");

			if (compressed_handle != NULL) {

				SGZip z;
				init_gzip_structure(&z);

				z.reading = false;
				z.first_time = false;
				z.strm.zalloc = Z_NULL;
				z.strm.zfree = Z_NULL;
				z.strm.opaque = Z_NULL;

				memset(z.in, '\x0', CHUNK);
				z.strm.next_in = z.in;
				z.strm.avail_in = 0;

				CALL_ZLIB(deflateInit2(&(z.strm), Z_DEFAULT_COMPRESSION, Z_DEFLATED, windowBits | GZIP_ENCODING, 8, Z_DEFAULT_STRATEGY));
				z.have = CHUNK;

				compressed_bytes_written = 0;
				gz_block_starting_pos = 0;
				bytes_write = 0;
				isEOF = false;
				flush = Z_FULL_FLUSH;

				do {
					/* read uncompressed input data */

					z.strm.next_in = z.in;
					z.strm.avail_in = 0;

					z.bytes_read = fread(z.in, sizeof(unsigned char),
							sizeof(z.in), uncompressed_handle);
					z.strm.avail_in = z.bytes_read;
					isEOF = feof(uncompressed_handle);

					if (z.bytes_read > 0) {
						/* Scan the data for IDs */

						gz_block_starting_pos = compressed_bytes_written;

						for (i = 0; i < z.strm.avail_in; i++) {
							private_scan_deflated_for_create_index_positions(z.in[i], i, &idx, gz_block_starting_pos);
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

							z.strm.next_out = z.out;
							z.strm.avail_out = CHUNK;

/*
							if (isEOF == false) {
								flush = Z_FULL_FLUSH;
							} else {
								flush = Z_FULL_FLUSH;
							}
*/

							CALL_ZLIB(deflate(&(z.strm), flush));

							z.have = CHUNK - z.strm.avail_out;

							/* write compressed data to the output file */

							bytes_write = fwrite(z.out, sizeof(char), z.have,
									compressed_handle);
							compressed_bytes_written += bytes_write;

						/*!! Removed } while (z.strm.avail_out == 0);*/ /* continue compressing until it does not exist more data to be compressed */
						} while (z.strm.avail_in != 0); /* continue compressing until it does not exist more data to be compressed */
					}
				} while (isEOF == false);

				/* Finish */
				CALL_ZLIB(deflate(&(z.strm), Z_FINISH));
				CALL_ZLIB(deflateEnd(&(z.strm)));

				fclose(compressed_handle);
			}

			fclose(uncompressed_handle);
		}

		save_index_to_file(index_file_name, &idx);

		unload_all_index_positions(&idx);
	}

	return ret;
}

gz_return fzseek(FILE *file_handle, SGZip *z, struct SGZIndex *idx, const char *search_id, long int *row_num, int from_last_search) {
	gz_return ret = GZ_OK;
	long int uncompressed_bytes = 0;
	int flush = 0;

    if (((((search_id != NULL) && (*row_num == -1))) || (((search_id == NULL) && (*row_num != -1)))) && 
        (file_handle != NULL) &&
        (z != NULL)) {
		struct SIndexPos *search = NULL;

		if ((search_id != NULL) && (*row_num == -1)) {
			search = get_index_position_by_id(search_id, idx, from_last_search);
		}

		if ((search_id == NULL) && (*row_num != -1)) {
			search = get_index_position_by_pos(*row_num, idx, from_last_search);
		}

		if (search != NULL) {
            *row_num = search->row_num;

			/*printf("Found: position for ID '%s' is (%ld, %ld)\n", search_id, search->position, search->offset);*/

			uncompressed_bytes = 0;
			flush = 0;

			(*z).reading = 1;
			if ((*z).file_compressed == 1) {

				memset((*z).in, '\x0', CHUNK);
				memset((*z).out, '\x0', CHUNK);

				inflateEnd(&((*z).strm));

				/* This is the first read of the input file. Let's initialize some variables */
				(*z).strm.zalloc = Z_NULL;
				(*z).strm.zfree = Z_NULL;
				(*z).strm.opaque = Z_NULL;
				(*z).strm.next_in = (*z).in;
				(*z).strm.avail_in = 0;
				(*z).first_time = false;
				(*z).strm.data_type = 0;
				(*z).strm.next_out = (*z).out;

				/* Go to the desired position */

				/* There is a difference between accessing the first block of data or accessing the others */
				if (search->position == 0) {
					/* Accessing the first block of data */
					CALL_ZLIB(inflateInit2(&((*z).strm), windowBits | ENABLE_ZLIB_GZIP));
					flush = Z_FULL_FLUSH;
				} else {
					/* Accessing the second block of data and others */
					CALL_ZLIB(inflateInit2(&((*z).strm), -windowBits));
					flush = Z_BLOCK;
				}

				fseek(file_handle, search->position, SEEK_SET);

				while (search->offset >= uncompressed_bytes) {
					(*z).strm.next_in = (*z).in;
					(*z).bytes_read = fread((*z).in, sizeof(unsigned char),
							sizeof((*z).in), file_handle);
					(*z).strm.avail_in = (*z).bytes_read;

					do {
						if ((*z).bytes_read > 0) { /* && ((*z).strm.avail_in > 0)) { */
							/* Uncompress the data */
							(*z).strm.avail_out = CHUNK;
							(*z).strm.next_out = (*z).out;

							CALL_ZLIB(inflate(&((*z).strm), flush));
							flush = Z_NO_FLUSH;

							(*z).have = CHUNK - (*z).strm.avail_out;
							uncompressed_bytes += (*z).have;
						}
					} while ((search->offset > uncompressed_bytes)
							&& ((*z).strm.avail_out == 0));
				}

				/* The next char to be returned will be the pointer by the offset */
				(*z).pointer = search->offset % CHUNK;

				/*
				printf("pointer=%ld, offset=%ld, CHUNK=%d\n", (*z).pointer, search->offset, CHUNK);
				printf("%s\n", (*z).out);*/

				if ((*z).bytes_read == 0) {
					CALL_ZLIB(inflateEnd(&((*z).strm)));
				}

			} else {
				ret = GZ_ERROR_DATA_FILE;
			}
		} else {
			ret = GZ_INDEX_KEY_NOT_FOUND;
		}
	} else {
		ret = GZ_PARAMS_ERROR;
    }

	return ret;
}

gz_return fzseekNearest(FILE *file_handle, SGZip *z, struct SGZIndex *idx, const char *search_id, long int max_search, long int *seq_id_found) {
    gz_return ret = GZ_OK;

    struct SIndexPos *last_search = NULL;
    long int row_num = -1;
    long int i = 0;

    /*
     Example:
     1:1 TTTTTTTTTTTTTTTT
     1:2 TTTTTTTTTTTTTTTT
     1:5 TTTTTTTTTTTTTTTT
     1:6 TTTTTTTTTTTTTTTT
     2:1 TTTTTTTTTTTTTTTT
     We are looking for "1:3"
    */

    *seq_id_found = -1;
    last_search = idx->last_search; /* store the current last_search position */

    if ((ret = fzseek(file_handle, z, idx, search_id, &row_num, false /*true*/)) != GZ_OK) {

        /* "1,3" does not exist. So, this function must return the "1,5" but the index does not store (x:y) in clear text
           so, the function must try "1:3", "1:4", "1:5" until it gets the sequence. */

        /* Let's separate "1:3" into two chars first="1" and second="3" */
        const char SEPARATOR[2] = ":";
        char *first;
        char *second;
        long int searched_sequence_id = 0;
        char new_search_id[200];
        char search_coords[200];
        char buffer[100];
        memset(search_coords, '\x0', 200);
        strcpy(search_coords, search_id);

        /* get the first token */
        first = strtok(search_coords, SEPARATOR);
        if (first != NULL) {
            second = strtok(NULL, SEPARATOR);
            if (second != NULL) {

                char *ptr;
                searched_sequence_id = strtol(second, &ptr, 10);

                i = 0;
                for (i = searched_sequence_id + 1; ((i < max_search) && (ret != GZ_OK)); i++) {
                    idx->last_search = last_search; /* restore the last_search position in order to start searching from this position (through the index) */

                    strcpy(new_search_id, first);
                    strcat(new_search_id, ":");
                    sprintf(buffer, "%ld", i);
                    strcat(new_search_id, buffer);

                    row_num = -1;
                    ret = fzseek(file_handle, z, idx, new_search_id, &row_num, false /*true*/);

                    if (ret == GZ_OK) {
                        *seq_id_found = i;
                    }
                }
            }
            else {
                ret = GZ_ERROR_DATA_FILE;
            }
        } else {
            ret = GZ_ERROR_DATA_FILE;
        }
    }

    if (ret != GZ_OK) {
        idx->last_search = last_search; /* restore the last_search position in order to start searching from this position (through the index) */
    }

    return ret;
}

unsigned long hash(unsigned char *str) {
	unsigned long hash = 5381;
	int c;

	while ((c = *str++) != 0) {
		hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
	}

	return hash;
}

gz_return add_index_position(struct SGZIndex *idx) {
	gz_return ret = GZ_OK;

	if (idx == NULL) {
		ret = GZ_PARAMS_ERROR;
	} else {
		struct SIndexPos *new_node = NULL;

		new_node = (struct SIndexPos *) malloc(sizeof(struct SIndexPos) * 1);

		new_node->key = 0;
		new_node->position = 0;
		new_node->offset = 0;
		new_node->next = NULL;

		if (idx->first == NULL) {
			idx->first = new_node;
		}
		if (idx->last != NULL) {
			idx->last->next = new_node;
		}
		idx->last = new_node;
		/* idx->items++; */
	}

	return ret;
}

struct SIndexPos *get_index_position_by_id(const char *search, struct SGZIndex *idx, int from_last_search) {
	struct SIndexPos *ret = NULL;
	unsigned long search_hash = hash((unsigned char *)search);
	struct SIndexPos *iterator = ((from_last_search == 1) && (idx->last_search != NULL))?idx->last_search:idx->first;

    while ((iterator != NULL) && (ret == NULL)) {
        if (iterator->key == search_hash) {
	        ret = iterator;
	        idx->last_search = iterator;
        }

        iterator = iterator->next;
    }

	return ret;
}

struct SIndexPos *get_index_position_by_pos(long int pos, struct SGZIndex *idx, int from_last_search) {
	struct SIndexPos *ret = NULL;
	long int current = ((from_last_search == 1) && (idx->last_search != NULL))?idx->last_search_count:0;

	struct SIndexPos *iterator = ((from_last_search == 1) && (idx->last_search != NULL))?idx->last_search:idx->first;

	while ((iterator != NULL) && (ret == NULL)) {
		if (pos == current) {
			ret = iterator;
			idx->last_search = iterator;
			idx->last_search_count = current;
		}

		iterator = iterator->next;
		current++;
	}

	return ret;
}

gz_return load_index_from_file(const char *file_name, struct SGZIndex *idx) {
	gz_return ret = GZ_OK;
    FILE *h = NULL;

	if ((file_name == NULL) || (idx == NULL)) {
		ret = GZ_PARAMS_ERROR;
	} else {

		init_gzindex_structure(idx);

		h = fopen((char *)file_name, "rb");

		if (h != NULL) {
			struct SIndexPosDisk obj;

			long int current_position = 0;

			while (!feof(h)) {
				if (fread(&obj, 1, sizeof(struct SIndexPosDisk), h) != 0) {
					if (obj.key == 0) { /* 0 means that it is not a real key. It means that pos_offset has the position for all the following keys */
						current_position = obj.pos_offset;
					} else {
                        add_index_position(idx);
						idx->last->key = obj.key;
						idx->last->position = current_position;
						idx->last->offset = obj.pos_offset;
                        idx->last->row_num = idx->items;
						idx->items++;
					}
				}
			}

			fclose(h);

		} else {
			return GZ_ERROR_OPEN_FILE;
		}
	}

	return ret;
}

gz_return save_index_to_file(const char *file_name, struct SGZIndex *idx) {
	gz_return ret = GZ_OK;

	if ((file_name == NULL) || (idx == NULL)) {
		ret = GZ_PARAMS_ERROR;
	} else {

		FILE *h = fopen((char *)file_name, "wb+");

		if (h != NULL) {
			struct SIndexPosDisk obj;
			struct SIndexPos *iterator = idx->first;
			long int current_position = -1;

			while (iterator != NULL) {
				if (iterator->position != current_position) {
					obj.key = 0; /* 0 means that it is not a real key. It means that pos_offset has the position for all the following keys */
					obj.pos_offset = iterator->position;
					fwrite(&obj, 1, sizeof(struct SIndexPosDisk), h);
					current_position = iterator->position;
				}
				obj.key = iterator->key;
				obj.pos_offset = iterator->offset;
				fwrite(&obj, 1, sizeof(struct SIndexPosDisk), h);
				iterator = iterator->next;
			}
			fclose(h);
		} else {
			ret = GZ_ERROR_CREATE_FILE;
		}
	}

	return ret;
}

gz_return unload_all_index_positions(struct SGZIndex *idx) {
	gz_return ret = GZ_OK;

	if (idx == NULL) {
		ret = GZ_PARAMS_ERROR;
	} else {
		struct SIndexPos *iterator = idx->first;
		struct SIndexPos *previous = NULL;

		while (iterator != NULL) {
			previous = iterator;
			iterator = iterator->next;

			if (previous != NULL) {
				free(previous);
				previous = NULL;
			}
		}

		idx->first = NULL;
		idx->last = NULL;
        idx->items = 0;
	}
	return ret;
}

