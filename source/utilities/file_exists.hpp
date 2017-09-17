#ifndef FILE_EXISTS_HPP
#define FILE_EXISTS_HPP

namespace Utilities {
// This method of checking if a file exists was based on
// https://stackoverflow.com/a/12774387
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