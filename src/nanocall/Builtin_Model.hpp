#ifndef __BUILTIN_MODEL_HPP
#define __BUILTIN_MODEL_HPP

#include <string>
#include <vector>

struct Builtin_Model
{
    static const unsigned num;
    static const unsigned strands[];
    static const std::string names[];
    static const std::vector< float > init_lists[];
}; // struct Builtin_Model

#endif
