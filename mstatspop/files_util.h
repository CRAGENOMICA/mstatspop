

#ifndef FILES__UTIL_H_
#define FILES__UTIL_H_

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <sys/stat.h>
#ifdef __cplusplus
extern "C"
{
#endif

    char *get_filepath(const char *fname);
    /**
     * Retrieves the base name of a file path.
     *
     * @param filepath The file path from which to extract the base name.
     * @param remove_extension Specifies whether to remove the file extension from the base name.
     * @return The base name of the file path with or without extension.
     */
    char *get_basename(const char *filepath, bool remove_extension);

    char *get_extension(const char *filepath);
    bool ends_with(const char *str, const char *suffix);

    /**
     * @brief Check if a file exists
     * @param filename The path to the file to check
     * @return 1 if the file exists, 0 otherwise
     */
    int file_exists(const char *filename);

#ifdef __cplusplus
}
#endif

#endif /* FILES__UTIL_H_ */
