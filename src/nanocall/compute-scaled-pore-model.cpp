#include <iostream>
#include <string>
#include <tclap/CmdLine.h>

#include "zstr.hpp"
#include "Pore_Model.hpp"
#include "logger.hpp"

using namespace std;

#ifndef FLOAT_TYPE
#define FLOAT_TYPE float
#endif
typedef Pore_Model< FLOAT_TYPE > Pore_Model_Type;
typedef Pore_Model_Type::Pore_Model_Parameters_Type Pore_Model_Parameters_Type;

namespace opts
{
    using namespace TCLAP;
    string description =
        "Compute scaled pore model.";
    CmdLine cmd_parser(description);
    MultiArg< string > log_level("d", "log-level", "Log level.", false, "string", cmd_parser);
    ValueArg< string > file_name("f", "file-name", "Fast5 file.", true, "", "file", cmd_parser);
    ValueArg< unsigned > strand("s", "strand", "Strand.", false, 0, "file", cmd_parser);
    ValueArg< string > output_file_name("o", "output", "Output file name.", false, "", "file", cmd_parser);
} // namespace opts

void real_main()
{
    Pore_Model_Type m;
    Pore_Model_Parameters_Type m_params;
    m.load_from_fast5(fast5::File(opts::file_name), opts::strand);
    m_params.load_from_fast5(fast5::File(opts::file_name), opts::strand);
    m.scale(m_params);
    if (not opts::output_file_name.get().empty())
    {
        strict_fstream::ofstream(opts::output_file_name.get()) << m;
    }
    else
    {
        cout << m;
    }
}

int main(int argc, char * argv[])
{
    opts::cmd_parser.parse(argc, argv);
    logger::Logger::set_levels_from_options(opts::log_level);
    real_main();
}
