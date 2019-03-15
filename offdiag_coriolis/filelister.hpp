#pragma once

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

class FileLister
{
public:
    explicit FileLister( std::string dir ) : dir(std::move(dir))
    {
        files = listdir();
        std::sort(files.begin(), files.end(), numeric_string_compare);
    }

    std::string full_path( std::string const& file );
    std::vector<std::string> listdir( );
    std::vector<std::string> const& get_files() const { return files; }

    static bool numeric_string_compare( const std::string& str1, const std::string& str2 );
    static int get_number( std::string str );
    static std::vector<std::string> split( const std::string& s, char delimiter );
private:
    std::string dir;
    std::vector<std::string> files;
};
