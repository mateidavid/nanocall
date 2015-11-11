#include <iostream>
#include <string>
#include <tclap/CmdLine.h>

#include "zstr.hpp"
#include "State_Transitions.hpp"
#include "logger.hpp"

using namespace std;

namespace opts
{
    using namespace TCLAP;
    string description =
        "Compute state transition probabilities based on the overlap model, for a given pr_skip and pr_stay.";
    CmdLine cmd_parser(description);
    MultiArg< string > log_level("d", "log-level", "Log level.", false, "string", cmd_parser);
    ValueArg< float > p_skip("k", "pr-skip", "Pr skip.", true, 0.1, "float", cmd_parser);
    ValueArg< float > p_stay("t", "pr-stay", "Pr stay.", true, 0.1, "float", cmd_parser);
    ValueArg< float > p_cutoff("p", "pr-cutoff", "Minimim prob to keep.", false, 0.001, "float", cmd_parser);
    ValueArg< string > output_file_name("o", "output", "Output file name.", false, "", "file", cmd_parser);
} // namespace opts

void real_main()
{
    State_Transitions<> st;
    st.compute_transitions(opts::p_skip, opts::p_stay, opts::p_cutoff);
    if (not opts::output_file_name.get().empty())
    {
        string fn = opts::output_file_name.get();
        /*
        if (fn.size() >= 4 and fn.substr(fn.size() - 3) == ".gz")
        {
            zstr::ofstream os(fn);
            print_state_transitions(st, os);
        }
        else
        */
        strict_fstream::ofstream(fn) << st;
    }
    else
    {
        cout << st;
    }
}

int main(int argc, char * argv[])
{
    opts::cmd_parser.parse(argc, argv);
    Logger::set_levels_from_options(opts::log_level);
    real_main();
}
