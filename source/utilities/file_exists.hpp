#ifndef FILE_EXISTS_HPP
#define FILE_EXISTS_HPP

#include "../general_definitions.hpp"

namespace Utilities {
/**
 * Determine whether a file exists.
 * Uses fopen to check if a file exists. This method of checking if a file exists was based on
 * performance measurements done by
 * <a href="https://stackoverflow.com/a/12774387">PherricOxide's answer</a>.
 *
 * @param file_name name of file
 * @return Whether or not the file exists
 */
inline bool file_exists(const std::string& file_name) {
    if (FILE* file = fopen(file_name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}
}
#endif