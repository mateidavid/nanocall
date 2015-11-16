#include <iostream>
#include <string>

#include "fs_support.hpp"

using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cerr << "use: " << argv[0] << " <directory>" << endl;
        return EXIT_FAILURE;
    }
    string file_name = argv[1];
    auto is_dir = is_directory(file_name);
    if (not is_dir)
    {
        cerr << "not a directory: " << file_name << endl;
        return EXIT_FAILURE;
    }
    auto l = list_directory(file_name);
    for (const auto& f : l)
    {
        cout << f << endl;
    }
}
