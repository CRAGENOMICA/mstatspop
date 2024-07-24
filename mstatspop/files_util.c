

#include "files_util.h"


/**
 * Retrieves the file path for the given file name.
 *
 * @param fname The name of the file.
 * @return The file path corresponding to the given file name.
 */
char *get_filepath(const char *fname)
{
  char *path = strdup(fname);
  char *last_slash = strrchr(path, '/');
  if (last_slash)
  {
    *last_slash = '\0'; // Remove the file name by setting a null character
  }
  else
  {
    path[0] = '.';  // Set the current directory
    path[1] = '\0'; // Null-terminate the string
  }
  return path;
}


char *get_basename(const char *filepath, bool remove_extension)
{
  const char *last_slash = strrchr(filepath, '/');
  // copy the file path to a new buffer
  const char *basename = last_slash ? strdup(last_slash + 1) : strdup(filepath);   

  if (!remove_extension)
  {
    return (char *)basename; // Return full basename with extension
  }
  else
  {
    //static char base[256];                     // Static buffer to hold the modified basename
    //strncpy(base, basename, sizeof(base) - 1); // Copy the basename to a local buffer
    //base[sizeof(base) - 1] = '\0';             // Ensure null-termination

    char *last_dot = strrchr(basename, '.');
    if (last_dot && last_dot != basename)
    {                   // Ensure the dot is not the first character (hidden files)
      *last_dot = '\0'; // Remove the extension
    }
    return (char *)basename;
  }
}

// Function to extract the file extension
char *get_extension(const char *filepath)
{
  const char *dot = strrchr(filepath, '.');
  if (dot && dot > strrchr(filepath, '/'))
  {
    // return a copy of the extension
    return strdup(dot + 1); // Return the substring after the last '.'
    // return dot + 1; // Return the substring after the last '.'
  }
  return ""; // No extension found
}
// Function to check if str ends with suffix
bool ends_with(const char *str, const char *suffix) {
    if (!str || !suffix) {
        return false;
    }
    size_t str_len = strlen(str);
    size_t suffix_len = strlen(suffix);
    if (suffix_len > str_len) {
        return false;
    }
    return strncmp(str + str_len - suffix_len, suffix, suffix_len) == 0;
}