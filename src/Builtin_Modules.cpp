#include "Builtin_Modules.hpp"

#include "t.006.ont.model.xxd"
#include "c.p1.006.ont.model.xxd"
#include "c.p2.006.ont.model.xxd"

const unsigned Builtin_Module::num = 3;
const unsigned Builtin_Module::strands[] = { 0, 1, 1 };
const std::string Builtin_Module::names[] = { "t.006", "c.p1.006", "c.p2.006" };
const unsigned char* Builtin_Module::raw_strings[] = { t_006_ont_model, c_p1_006_ont_model, c_p2_006_ont_model };
const unsigned Builtin_Module::raw_string_lens[] = { t_006_ont_model_len, c_p1_006_ont_model_len, c_p2_006_ont_model_len };
