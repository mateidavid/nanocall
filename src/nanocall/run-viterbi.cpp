#include <iostream>
#include <string>
#include <tclap/CmdLine.h>

#include "Pore_Model.hpp"
#include "State_Transitions.hpp"
#include "Event.hpp"
#include "Viterbi.hpp"
#include "logger.hpp"
#include "zstr.hpp"

using namespace std;

#ifndef FLOAT_TYPE
#define FLOAT_TYPE float
#endif
typedef State_Transitions< FLOAT_TYPE > State_Transitions_Type;
typedef Pore_Model< FLOAT_TYPE > Pore_Model_Type;
typedef Event< FLOAT_TYPE > Event_Type;
typedef Event_Sequence< FLOAT_TYPE > Event_Sequence_Type;
typedef Viterbi< FLOAT_TYPE > Viterbi_Type;

namespace opts
{
    using namespace TCLAP;
    string description =
        "Run Viterbi on given input";
    CmdLine cmd_parser(description);
    MultiArg< string > log_level("d", "log-level", "Log level.", false, "string", cmd_parser);
    ValueArg< string > pm_file_name("p", "pore-model", "Scaled pore model file name.", true, "", "file", cmd_parser);
    ValueArg< string > st_file_name("s", "state-transitions", "State transitions file name.", true, "", "file", cmd_parser);
    ValueArg< string > ev_file_name("e", "events", "Events file name.", true, "", "file", cmd_parser);
} // namespace opts

void real_main()
{
    Pore_Model_Type pm;
    State_Transitions_Type st;
    Event_Sequence_Type ev;
    zstr::ifstream(opts::pm_file_name) >> pm;
    zstr::ifstream(opts::st_file_name) >> st;
    {
        zstr::ifstream ifs(opts::ev_file_name);
        Event_Type e;
        while (ifs >> e)
        {
            ev.push_back(e);
        }
    }

    Viterbi_Type vit;
    vit.fill(pm, st, ev);
    cout << vit.base_seq() << std::endl;
}

int main(int argc, char * argv[])
{
    opts::cmd_parser.parse(argc, argv);
    logger::Logger::set_levels_from_options(opts::log_level);
    real_main();
}
