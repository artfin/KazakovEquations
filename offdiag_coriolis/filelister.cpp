//
// Created by artfin on 01.03.19.
//

#include "filelister.hpp"

int FileLister::get_number( std::string str )
{
    // For atoi, the input string has to start with a digit, so lets search for the first digit
    size_t i = 0;
    for (; i < str.length(); i++) { if (isdigit(str[i])) break; }

    // remove the first chars, which aren't digits
    str = str.substr(i, str.length() - i);

    for (; i < str.length(); i++) { if (!isdigit(str[i])) break; }
    // remove the last chars, which aren't digits
    str = str.substr(0, i);

    // convert the remaining text to an integer
    int n = std::stoi(str);
    return n;
}

std::string FileLister::full_path( std::string const& file )
{
    return dir + file;
}

std::vector<std::string> FileLister::split(const std::string& s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<std::string> FileLister::listdir( )
{
    DIR *dp;
    struct dirent *dirp;

    if( (dp  = opendir(dir.c_str())) == nullptr ) {
        std::cerr << "Error(" << errno << ") opening " << dir << std::endl;
        exit( 1 );
    }

    std::vector<std::string> files;
    while ( (dirp = readdir(dp)) != nullptr ) {
        std::string file = std::string(dirp->d_name);
        if ( file == "." || file == ".." ) continue;

        // игнорируем буферные файлы vim'a
        size_t found = file.find(".swp");
        if ( found != std::string::npos ) continue;

        // игнорируем питоновские скрипты
        found = file.find(".py");
        if ( found != std::string::npos ) continue;

        files.push_back(full_path(file));
    }

    closedir(dp);

    return files;
}

bool FileLister::numeric_string_compare( std::string const& str1, std::string const& str2 )
{
    std::string filename1 = split(str1, '/')[2];
    std::string filename2 = split(str2, '/')[2];
    return get_number(filename1) < get_number(filename2);
}