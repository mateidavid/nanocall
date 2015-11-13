#ifndef __BUILTIN_MODULES_HPP
#define __BUILTIN_MODULES_HPP

#include <string>
#include <vector>

struct Builtin_Module
{
    static const unsigned num;
    static const unsigned strands[];
    static const std::string names[];
    static const unsigned char* raw_strings[];
    static const unsigned raw_string_lens[];
}; // struct Builtin_Module

#endif
