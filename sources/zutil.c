/*
 * mstatspop, Statistical Analysis using Multiple Populations for Genomic Data.
 * Created by Sebastian E. Ramos Onsins.
 */
/**
 *  \brief     zutil.c
 *  \details
 *  \author    Joan Jené
 *  \version   1.5
 *  \date      Nov 23, 2016
 *  \pre
 *  \bug
 *  \warning
 *  \copyright
 */

#include "zutil.h"

FILE * fzopen(const char *filename, const char *opentype, SGZip *z) {
	init_gzip_structure(z);

	(*z).file_compressed = 0;

	if (strlen(filename) > 3) {
		if((filename[strlen(filename)-3] == '.') &&
		   (filename[strlen(filename)-2] == 'g') &&
		   (filename[strlen(filename)-1] == 'z')) {
			(*z).file_compressed = 1;
		}
	}

	return fopen(filename, opentype);
}

void init_gzip_structure(SGZip *z) {
  (*z).bytes_read = 0;
  (*z).have = 0;
  (*z).pointer = 0;
  (*z).first_time = 1;
  (*z).file_compressed = 1;
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
  int ret;

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
      if (((*z).strm.avail_out == 0) && ((*z).first_time == 0)) {
        /* The buffer of uncompressed chars is empty but
           It exists read data that can be uncompressed */
      } else {
        /* The buffer of uncompressed chars is empty and
           It does not exist read data that could be uncompressed, so
           Read data from the compressed input file */
        (*z).strm.next_in = (*z).in;
        (*z).bytes_read = fread((*z).in, sizeof (unsigned char), sizeof ((*z).in), file_handle);
        (*z).strm.avail_in = (*z).bytes_read;
      }

      if ((*z).bytes_read > 0) {
        /* Uncompress the data */
        (*z).strm.avail_out = CHUNK;
        (*z).strm.next_out = (*z).out;
        CALL_ZLIB(inflate(&((*z).strm), Z_NO_FLUSH));
        (*z).have = CHUNK - (*z).strm.avail_out;
        (*z).first_time = 0;

        num_chars += (*z).have;
      }

      /* The next char to be returned will be the first one in the buffer */
      (*z).pointer = 0;
    }


    if ((*z).bytes_read == 0) {
      inflateEnd(& ((*z).strm));
      ret = -1;
    } else {
      /* Return the current char */
      ret = (*z).out[(*z).pointer];

      /* Increment the pointer for the next call to this function */
      (*z).pointer++;
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
    while((i < count) && (!fzeof(file_handle, z))) {
    	row[i] = fzgetc(file_handle, z);
    	i++;
    }
  }

  return ret;
}

/*
 * Do not used this function directly. This function is used by the fzprintf function.
 */
int private_fzprintf(FILE * file_handle, SGZip *z, char *message) {
  int ret = 0;
  int i = 0;
  int last_i = 0;

  if (z == 0) {
	  printf("z is 0");
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

			CALL_ZLIB(deflateInit2(&((*z).strm), Z_DEFAULT_COMPRESSION, Z_DEFLATED,
					windowBits | GZIP_ENCODING, 8,
					Z_DEFAULT_STRATEGY));

			(*z).have = CHUNK; /* the next_in buffer has all bytes free */
		}

		if ((*z).strm.avail_in < CHUNK) {
			/* There is still space inside the "(*z).in" buffer of uncompressed chars */
			for (i = 0; ((i < strlen(message)) && ((*z).strm.avail_in + i < CHUNK)); i++) {
				(*z).in[(*z).strm.avail_in + i] = message[i];
			}
			(*z).strm.avail_in += i;
			/* if (i > 0) { last_i = i - 1; } else { last_i = 0; } */
			last_i = i;


			if ((*z).strm.avail_in < CHUNK) {
				/* There is still space inside the "(*z).in" buffer of uncompressed chars. Deflate them the next time... */
			} else {
				/* Deflate now */
				do {
				  (*z).strm.avail_out = CHUNK;
				  (*z).strm.next_out = (*z).out;
				  CALL_ZLIB(deflate(&((*z).strm), Z_NO_FLUSH));
				  (*z).have = CHUNK - (*z).strm.avail_out;
				  fwrite((*z).out, sizeof (char), (*z).have, file_handle);
				} while ((*z).strm.avail_out == 0);

				/* Reset */
				memset((*z).in, '\x0', CHUNK);
				(*z).strm.next_in = (*z).in;
				(*z).strm.avail_in = 0;

				/* if there are more chars in the current message, then insert them in the z.in buffer: */
				/* how can it be possible? because of this: (for example)
					 - previous message size : 1 byte
					 - current message size  : CHUCK bytes
					 - total space in (*z).strm.in is full because it is of CHUNK size but still remains 1 byte to be written into the disk.
					   So,  this byte is stored in the empty (*z).in for the next call.
				*/
				for (i = last_i; ((i < strlen(message)) && ((*z).strm.avail_in + i - last_i < CHUNK)); i++) {
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

int fzprintf(FILE * file_handle, SGZip *z, char *message, ...) {
  int ret = 0;


  char *buffer       = NULL;
  size_t message_len = 0;
  size_t sent_chars  = 0;
  int replace_args   = 0;
  char buffer_block[CHUNK + 1]; // +1 for the \x0

  message_len  = strlen(message);
  replace_args = (message_len < MAX_FZPRINTF_MESSAGE);

  if (replace_args) {
	  buffer = (char *)malloc(2 * MAX_FZPRINTF_MESSAGE * sizeof(char));
	  va_list args;
	  va_start (args, message);
	  vsnprintf(buffer,2 * MAX_FZPRINTF_MESSAGE,message, args);
	  va_end (args);
  }
  else {
	  buffer = message;
  }

  if ((z != 0) && ((*z).file_compressed == 0)) {
	  /* The file is not compressed */
	  ret = fprintf(file_handle, buffer);
  } else {
	  /* If the message is smaller than CHUNK then no problem */
	  message_len  = strlen(buffer);

	  if (message_len < CHUNK) {
		  ret += private_fzprintf(file_handle, z, buffer);
	  } else {
		  /* But if the message is longer than CHUNK then it must be divided into small blocks of CHUNK size */


		  /* Send all blocks of CHUNK size */
		  sent_chars = 0;
		  while (message_len - sent_chars >= CHUNK) {
			  memcpy(buffer_block, buffer + sent_chars, CHUNK);
			  buffer_block[CHUNK]='\x0';
			  ret += private_fzprintf(file_handle, z, buffer_block);
			  sent_chars += CHUNK;
		  }

		  /* Send the last block with size less than CHUNK */
		  int rest = message_len - sent_chars;
		  if (rest < CHUNK) {
			  memcpy(buffer_block, buffer + sent_chars, rest);
			  buffer_block[rest]='\x0';
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

  return ret;
}


int fzeof(FILE * file_handle, SGZip *z) {
  int ret = feof(file_handle);

  if ((*z).file_compressed == 1) {
    if ((ret = ((ret) && /* Reached the end of the file and */
            ((*z).strm.avail_out != 0) && /* there is no more input data to uncompress and */
            ((*z).pointer >= (*z).have))) == 1) { /* there is no more chars in the buffer */
      inflateEnd(&((*z).strm));

      /* Initialize first time to true */
      (*z).first_time = 1;
    }
  }

  return ret;
}

void fzclose(FILE * file_handle, SGZip *z) {
  if ((file_handle != 0) && (z != 0) && ((*z).reading != -1)) {
	  if (((*z).file_compressed == 0) || ((*z).reading == 1)) {
	  } else {
		  if (((*z).file_compressed == 1) && ((*z).reading == 0)) {
			/* Deflate the z.in buffer into the z.out buffer (if there are chars) */
			do {
			  (*z).strm.avail_out = CHUNK;
			  (*z).strm.next_out = (*z).out;
			  CALL_ZLIB(deflate(&((*z).strm), Z_FINISH));
			  (*z).have = CHUNK - (*z).strm.avail_out;
			  fwrite((*z).out, sizeof (char), (*z).have, file_handle);
			} while ((*z).strm.avail_out == 0);

			CALL_ZLIB(deflateEnd(&((*z).strm)));
		  }
	  }
	  fclose(file_handle);
	  file_handle = 0;
	  (*z).reading = -1; /* file closed */
  }
}

int memory_deflate(char *file_name, char *start_address, char *end_address) {
  int ret = 0;
  FILE *file_handle = NULL;

  if ((file_handle = fopen(file_name, "wb+")) != NULL) {
    SGZip z;
    z.first_time = 0;
    z.strm.zalloc = Z_NULL;
    z.strm.zfree = Z_NULL;
    z.strm.opaque = Z_NULL;
    CALL_ZLIB(deflateInit2(&(z.strm), Z_DEFAULT_COMPRESSION, Z_DEFLATED,
            windowBits | GZIP_ENCODING, 8,
            Z_DEFAULT_STRATEGY));


    z.strm.next_in = (unsigned char *) start_address;
    z.strm.avail_in = end_address - start_address;
    do {
      z.strm.avail_out = CHUNK;
      z.strm.next_out = z.out;
      CALL_ZLIB(deflate(&(z.strm), Z_FINISH));
      z.have = CHUNK - z.strm.avail_out;
      fwrite(z.out, sizeof (char), z.have, file_handle);
    } while (z.strm.avail_out == 0);
  } else {
    ret = 1;
  }

  return ret;
}

int uncompress_file(const char * compressed_file_name, const char * uncompressed_file_name) {
  int ret = 0;

  gzFile compressed_file;
  FILE * uncompressed_file;

  unsigned char buffer[CHUNK];

  compressed_file = gzopen(compressed_file_name, "rb");
  if (compressed_file != NULL) {
    uncompressed_file = fopen(uncompressed_file_name, "wb+");
    if (uncompressed_file != NULL) {
      /* ................................................................. */
      int bytes_read;

      while (!gzeof(compressed_file)) {
        memset(buffer, '\x0', sizeof (buffer));
        bytes_read = gzread(compressed_file, buffer, CHUNK - 1);
        fwrite(buffer, sizeof (unsigned char), bytes_read, uncompressed_file);
      }

      /* ................................................................. */
      fclose(uncompressed_file);
    } else {
      ret = 1;
    }

    gzclose(compressed_file);
  } else {
    ret = 1;
  }

  return ret;
}

int compress_file(const char * uncompressed_file_name, const char * compressed_file_name) {
  int ret = 0;

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
        memset(buffer, '\x0', sizeof (buffer));
        bytes_read = fread(buffer, sizeof (char), sizeof (buffer), uncompressed_file);

        if (bytes_read > 0) {
          gzwrite(compressed_file, buffer, strlen(buffer));
        }
      }
      /* ................................................................. */
      gzclose(compressed_file);
    } else {
      ret = 1;
    }

    fclose(uncompressed_file);
  } else {
    ret = 1;
  }

  return ret;
}


