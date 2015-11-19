#include <iostream>
#include <string>
#include <tclap/CmdLine.h>

#include "Pore_Model.hpp"
#include "State_Transitions.hpp"
#include "Event.hpp"
#include "Forward_Backward.hpp"
#include "logger.hpp"
#include "zstr.hpp"

using namespace std;

namespace opts
{
    using namespace TCLAP;
    string description =
        "Given a scaled pore model, a state trasition table, and a sequence of events, "
        "compute the state distribution conditioned on the prefix event sequence";
    CmdLine cmd_parser(description);
    MultiArg< string > log_level("d", "log-level", "Log level.", false, "string", cmd_parser);
    ValueArg< string > pm_file_name("p", "pore-model", "Scaled pore model file name.", true, "", "file", cmd_parser);
    ValueArg< string > st_file_name("s", "state-transitions", "State transitions file name.", true, "", "file", cmd_parser);
    ValueArg< string > ev_file_name("e", "events", "Events file name.", true, "", "file", cmd_parser);
    ValueArg< string > output_file_name("o", "output", "Output file name.", false, "", "file", cmd_parser);
} // namespace opts

void real_main()
{
    Pore_Model<> pm;
    //Pore_Model_Parameters<> params;
    State_Transitions<> st;
    vector< Event<> > ev;
    zstr::ifstream(opts::pm_file_name) >> pm;
    zstr::ifstream(opts::st_file_name) >> st;
    {
        zstr::ifstream ifs(opts::ev_file_name);
        Event<> e;
        while (ifs >> e)
        {
            ev.push_back(e);
        }
    }

    Forward_Backward<> fwbw;
    fwbw.fill(pm, st, ev);

    // print all kmers with posterior >= .1 for the middle event
    multiset< pair< float, unsigned > > s;
    for (unsigned j = 0; j < pm.n_states; ++j)
    {
        float v = exp(fwbw.cell(ev.size() / 2, j).gamma);
        if (v >= .1)
        {
            s.insert(make_pair(v, j));
        }
    }
    while (not s.empty())
    {
        auto it = prev(s.end());
        cout << Forward_Backward<>::Kmer_Type::to_string(it->second) << '\t' << it->first << endl;
        s.erase(it);
    }

    if (not opts::output_file_name.get().empty())
    {
        strict_fstream::ofstream(opts::output_file_name) << fwbw;
    }
}

int main(int argc, char * argv[])
{
    opts::cmd_parser.parse(argc, argv);
    logger::Logger::set_levels_from_options(opts::log_level);
    real_main();
}
