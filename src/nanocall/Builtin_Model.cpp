#include "Builtin_Model.hpp"

const unsigned Builtin_Model::num =
#include "builtin_model_num.inl"
    ;

const unsigned Builtin_Model::strands[] =
#include "builtin_model_strands.inl"
    ;

const std::string Builtin_Model::names[] =
#include "builtin_model_names.inl"
    ;

const std::vector< float > Builtin_Model::init_lists[] =
#include "builtin_model_init_lists.inl"
    ;
