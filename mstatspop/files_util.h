

#ifndef FILES__UTIL_H_
#define FILES__UTIL_H_

#include <stdio.h>
#include <stdbool.h>
#include <string.h>

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

#ifdef __cplusplus
}
#endif

#endif /* FILES__UTIL_H_ */
