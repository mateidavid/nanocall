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
    ValueArg< string > output_file_name("o", "output", "Output file name.", false, "", "file", cmd_parser);
    ValueArg< float > p_cutoff("p", "pr-cutoff", "Minimim prob to keep.", false, 0.001, "float", cmd_parser);
    ValueArg< float > p_skip("k", "pr-skip", "Pr skip.", false, 0.28, "float", cmd_parser);
    ValueArg< float > p_stay("t", "pr-stay", "Pr stay.", false, 0.09, "float", cmd_parser);
    SwitchArg fast("", "fast", "Use fast computation.", cmd_parser);
} // namespace opts

void real_main()
{
    State_Transitions<> st;
    if (opts::fast)
    {
        st.compute_transitions_fast(opts::p_skip, opts::p_stay);
    }
    else
    {
        st.compute_transitions(opts::p_skip, opts::p_stay, opts::p_cutoff);
    }
    if (not opts::output_file_name.get().empty())
    {
        strict_fstream::ofstream(opts::output_file_name) << st;
    }
    else
    {
        cout << st;
    }
}

int main(int argc, char * argv[])
{
    opts::cmd_parser.parse(argc, argv);
    logger::Logger::set_levels_from_options(opts::log_level);
    real_main();
}
