#include <deque>
#include <iostream>
#include <string>
#include <tclap/CmdLine.h>

#include "version.hpp"
#include "Pore_Model.hpp"
#include "Builtin_Model.hpp"
#include "State_Transitions.hpp"
#include "Event.hpp"
#include "Fast5_Summary.hpp"
#include "Viterbi.hpp"
#include "Forward_Backward.hpp"
#include "logger.hpp"
#include "alg.hpp"
#include "zstr.hpp"
#include "fast5.hpp"
#include "pfor.hpp"
#include "fs_support.hpp"

using namespace std;

namespace opts
{
    using namespace TCLAP;
    string description = "Call bases in Oxford Nanopore reads.";
    CmdLine cmd_parser(description, ' ', package_version);
    MultiArg< string > log_level("", "log", "Log level.", false, "string", cmd_parser);
    ValueArg< string > stats_fn("", "stats", "Stats.", false, "", "file", cmd_parser);
    ValueArg< unsigned > min_read_len("", "min-len", "Minimum read length.", false, 1000, "int", cmd_parser);
    ValueArg< unsigned > fasta_line_width("", "fasta-line-width", "Maximum fasta line width.", false, 80, "int", cmd_parser);
    ValueArg< float > pr_cutoff("", "pr-cutoff", "Minimum value for transition probabilities; smaller values are set to zero.", false, .001, "float", cmd_parser);
    ValueArg< float > pr_skip("", "pr-skip", "Transition probability of skipping at least 1 state.", false, .1, "float", cmd_parser);
    ValueArg< float > pr_stay("", "pr-stay", "Transition probability of staying in the same state.", false, .1, "float", cmd_parser);
    ValueArg< string > trans_fn("s", "trans", "Custom initial state transitions.", false, "", "file", cmd_parser);
    ValueArg< string > model_fofn("", "model-fofn", "File of pore models.", false, "", "file", cmd_parser);
    MultiArg< string > model_fn("m", "model", "Custom pore model.", false, "file", cmd_parser);
    ValueArg< string > output_fn("o", "output", "Output.", false, "", "file", cmd_parser);
    ValueArg< unsigned > num_threads("t", "threads", "Number of parallel threads.", false, 1, "int", cmd_parser);
    UnlabeledMultiArg< string > input_fn("inputs", "Inputs. Accepts: directories, fast5 files, or files of fast5 file names (use \"-\" to read fofn from stdin).", true, "path", cmd_parser);
} // namespace opts

void init_models(Pore_Model_Dict_Type& models)
{
    auto parse_model_name = [] (const string& s) {
        if (s.size() < 3
            or (s[0] != '0' and s[0] != '1' and s[0] != '2')
            or s[1] != ':')
        {
            cerr << "could not parse model name: \"" << s << "\"; format should be \"[0|1|2]:<file>\"" << endl;
            exit(EXIT_FAILURE);
        }
        unsigned st = s[0] - '0';
        return make_pair(st, s.substr(2));
    };

    map< unsigned, list< string > > model_list;
    if (not opts::model_fn.get().empty())
    {
        for (const auto& s : opts::model_fn.get())
        {
            auto p = parse_model_name(s);
            model_list[p.first].push_back(p.second);
        }
    }
    if (not opts::model_fofn.get().empty())
    {
        zstr::ifstream ifs(opts::model_fofn);
        string s;
        while (getline(ifs, s))
        {
            auto p = parse_model_name(s);
            model_list[p.first].push_back(p.second);
        }
    }
    if (model_list[2].empty() and (model_list[0].empty() != model_list[1].empty()))
    {
        cerr << "models were specified only for strand " << (int)model_list[0].empty()
             << "! give models for both strands, or for neither." << endl;
        exit(EXIT_FAILURE);
    }
    if (not (model_list[0].empty() and model_list[1].empty() and model_list[2].empty()))
    {
        for (unsigned st = 0; st < 3; ++st)
        {
            for (const auto& e : model_list[st])
            {
                Pore_Model_Type pm;
                string pm_name = e;
                zstr::ifstream(e) >> pm;
                pm.strand() = st;
                models[pm_name] = move(pm);
                LOG(info) << "loaded module [" << pm_name << "] for strand [" << st << "]" << endl;
            }
        }
    }
    else
    {
        // use built-in models
        for (unsigned i = 0; i < Builtin_Model::num; ++i)
        {
            Pore_Model_Type pm;
            string pm_name = Builtin_Model::names[i];
            pm.load_from_vector(Builtin_Model::init_lists[i]);
            pm.strand() = Builtin_Model::strands[i];
            models[Builtin_Model::names[i]] = move(pm);
            LOG(info)
                << "loaded builtin module [" << Builtin_Model::names[i] << "] for strand ["
                << Builtin_Model::strands[i] << "] statistics [mean=" << pm.mean() << ", stdv="
                << pm.stdv() << "]" << endl;
        }
    }
} // init_models

void init_transitions(State_Transitions_Type& transitions)
{
    if (not opts::trans_fn.get().empty())
    {
        zstr::ifstream(opts::trans_fn) >> transitions;
        LOG(info) << "loaded state transitions from [" << opts::trans_fn.get() << "]" << endl;
    }
    else
    {
        transitions.compute_transitions(opts::pr_skip, opts::pr_stay, opts::pr_cutoff);
        LOG(info) << "initialized state transitions with parameters p_skip=[" << opts::pr_skip
                  << "], pr_stay=[" << opts::pr_stay << "], pr_cutoff=[" << opts::pr_cutoff << "]" << endl;
    }
} // init_transitions

// Parse command line arguments. For each of them:
// - if it is a directory, find all fast5 files in it, ignore non-fast5 files.
// - if it is a file, check that it is indeed a fast5 file.
void init_files(list< string >& files)
{
    for (const auto& f : opts::input_fn)
    {
        if (is_directory(f))
        {
            auto l = list_directory(f);
            for (const auto& g : l)
            {
                string f2 = f + (f[f.size() - 1] != '/'? "/" : "") + g;
                if (is_directory(f2))
                {
                    LOG(info) << "ignoring subdirectory [" << f2 << "]" << endl;
                }
                else if (fast5::File::is_valid_file(f2))
                {
                    files.push_back(f2);
                    LOG(info) << "adding input file [" << f2 << "]" << endl;
                }
                else
                {
                    LOG(info) << "ignoring file [" << f2 << "]" << endl;
                }
            }
        }
        else // not a directory
        {
            if (f != "-" and fast5::File::is_valid_file(f))
            {
                files.push_back(f);
                LOG(info) << "adding input file [" << f << "]" << endl;
            }
            else // not fast5, interpret as fofn
            {
                LOG(info) << "interpreting [" << f << "] as fofn" << endl;
                istream* is_p = nullptr;
                strict_fstream::ifstream ifs;
                if (f == "-")
                {
                    is_p = &cin;
                }
                else
                {
                    ifs.open(f);
                    is_p = &ifs;
                }
                string g;
                while (getline(*is_p, g))
                {
                    if (fast5::File::is_valid_file(g))
                    {
                        files.push_back(g);
                        LOG(info) << "adding input file [" << g << "]" << endl;
                    }
                }
            }
        }
    }
} // init_files

void init_reads(const Pore_Model_Dict_Type& models,
                const list< string >& files,
                deque< Fast5_Summary_Type >& reads)
{
    for (const auto& f : files)
    {
        Fast5_Summary_Type s(f, models);
        LOG(info) << "summary: " << s << endl;
        if (s.have_ed_events
            and (s.strand_bounds[1] >= s.strand_bounds[0] + opts::min_read_len
                 or s.strand_bounds[3] >= s.strand_bounds[2] + opts::min_read_len))
        {
            reads.emplace_back(move(s));
        }
    }
} // init_reads

void train_reads(const Pore_Model_Dict_Type& models,
                 State_Transitions_Type& transitions,
                 deque< Fast5_Summary_Type >& reads)
{
    return;
    //
    // this is only a stub
    //
    for (unsigned round = 0; round < 5; ++round)
    {
        LOG(info) << "starting training round [" << round << "]" << endl;
        unsigned crt_idx = 0;
        pfor::pfor< unsigned >(
            opts::num_threads,
            10,
            // get_item
            [&] (unsigned& i) {
                if (crt_idx >= reads.size()) return false;
                i = crt_idx++;
                return true;
            },
            // process item
            [&] (unsigned& i) {
                Fast5_Summary_Type& read_summary = reads[i];
                read_summary.load_events();

                for (unsigned st = 0; st < 2; ++st)
                {
                    // if not enough events, ignore strand
                    if (read_summary.events[st].size() < opts::min_read_len) continue;
                    // create list of models to try
                    list< string > model_sublist;
                    if (models.count(read_summary.preferred_model[st]))
                    {
                        // if we have a preferred model, use that
                        model_sublist.push_back(read_summary.preferred_model[st]);
                    }
                    else
                    {
                        // no preferred model, try all that apply to this strand
                        for (const auto& p : models)
                        {
                            if (p.second.strand() == st or p.second.strand() == 2)
                            {
                                model_sublist.push_back(p.first);
                            }
                        }
                    }
                    assert(not model_sublist.empty());
                    // run fwbw
                    map< string, Forward_Backward_Type > results;
                    for (const auto& m_name : model_sublist)
                    {
                        // scale model using current parameters, initialize them if necessary
                        Pore_Model_Type pm(models.at(m_name));
                        Pore_Model_Parameters_Type pm_params = read_summary.params[st][m_name];
                        pm.scale(pm_params);
                        // main work: fill matrix
                        results[m_name].fill(pm, transitions, read_summary.events[st]);
                    }
                    //
                    // use fwbw results to update parameters
                    // ...
                    //
                }
                read_summary.drop_events();
            },
            // progress_report
            [&] (unsigned items, unsigned seconds) {
                clog << "Processed " << setw(6) << right << items << " reads in "
                     << setw(6) << right << seconds << " seconds\r";
            }); // pfor

        //
        // update state_transitions
        // ...
        //

    } // for round
} // train_reads

void write_fasta(ostream& os, const string& name, const string& seq)
{
    os << ">" << name << endl;
    for (unsigned pos = 0; pos < seq.size(); pos += opts::fasta_line_width)
    {
        os << seq.substr(pos, opts::fasta_line_width) << endl;
    }
} // write_fasta

void basecall_reads(const Pore_Model_Dict_Type& models,
                    const State_Transitions_Type& transitions,
                    deque< Fast5_Summary_Type >& reads)
{
    strict_fstream::ofstream ofs;
    ostream* os_p = nullptr;
    if (not opts::output_fn.get().empty())
    {
        ofs.open(opts::output_fn);
        os_p = &ofs;
    }
    else
    {
        os_p = &cout;
    }

    unsigned crt_idx = 0;
    pfor::pfor< unsigned, std::ostringstream >(
        opts::num_threads,
        10,
        // get_item
        [&] (unsigned& i) {
            if (crt_idx >= reads.size()) return false;
            i = crt_idx++;
            return true;
        },
        // process_item
        [&] (unsigned& i, std::ostringstream& oss) {
            Fast5_Summary_Type& read_summary = reads[i];
            read_summary.load_events();
            for (unsigned st = 0; st < 2; ++st)
            {
                // if not enough events, ignore strand
                if (read_summary.events[st].size() < opts::min_read_len) continue;
                // create list of models to try
                list< string > model_sublist;
                if (models.count(read_summary.preferred_model[st]))
                {
                    // if we have a preferred model, use that
                    model_sublist.push_back(read_summary.preferred_model[st]);
                }
                else
                {
                    // no preferred model, try all that apply to this strand
                    for (const auto& p : models)
                    {
                        if (p.second.strand() == st or p.second.strand() == 2)
                        {
                            model_sublist.push_back(p.first);
                        }
                    }
                }
                // check main scaling parameters
                auto r = alg::mean_stdv_of< float >(
                    read_summary.events[st],
                    [] (const Event_Type& ev) { return ev.mean; });
                LOG(debug)
                    << "mean_stdv read [" << read_summary.read_id
                    << "] strand [" << st
                    << "] ev_mean=[" << r.first
                    << "] ev_stdv=[" << r.second << "]" << endl;
                // deque of results
                deque< tuple< float, string, string > > results;
                for (const auto& m_name : model_sublist)
                {
                    // scale model, initialize default parameters if necessary
                    Pore_Model_Type pm(models.at(m_name));
                    Pore_Model_Parameters_Type pm_params = read_summary.params[st][m_name];
                    pm.scale(pm_params);
                    LOG(info)
                        << "basecalling read [" << read_summary.read_id
                        << "] strand [" << st
                        << "] model [" << m_name
                        << "] parameters " << pm_params << endl;
                    LOG(debug)
                        << "mean_stdv read [" << read_summary.read_id
                        << "] strand [" << st
                        << "] model_mean [" << pm.mean()
                        << "] model_stdv [" << pm.stdv() << "]" << endl;
                    if (std::abs(r.first - pm.mean()) > 5.0)
                    {
                        LOG(warning)
                            << "means_apart read [" << read_summary.read_id
                            << "] strand [" << st
                            << "] model [" << m_name
                            << "] parameters " << read_summary.params[st].at(m_name)
                            << " model_mean=[" << pm.mean()
                            << "] events_mean=[" << r.first
                            << "]" << endl;
                    }
                    Viterbi_Type vit;
                    vit.fill(pm, transitions, read_summary.events.at(st));
                    results.emplace_back(make_tuple(vit.path_probability(), m_name, vit.base_seq()));
                }
                std::sort(results.begin(), results.end());
                string& best_m_name = get<1>(results.back());
                string& base_seq = get<2>(results.back());
                LOG(info)
                    << "best_model read [" << read_summary.read_id
                    << "] strand [" << st
                    << "] model [" << best_m_name
                    << "] parameters " << read_summary.params[st].at(best_m_name) << endl;
                read_summary.preferred_model[st] = best_m_name;
                ostringstream tmp;
                tmp << read_summary.read_id << ":" << read_summary.base_file_name << ":" << st;
                write_fasta(oss, tmp.str(), base_seq);
            } // for st
            read_summary.drop_events();
        },
        // output_chunk
        [&] (std::ostringstream& oss) {
            *os_p << oss.str();
        },
        // progress_report
        [&] (unsigned items, unsigned seconds) {
            clog << "Processed " << setw(6) << right << items << " reads in "
                 << setw(6) << right << seconds << " seconds\r";
        }); // pfor
} // basecall_reads

void real_main()
{
    Pore_Model_Dict_Type models;
    State_Transitions_Type transitions;
    deque< Fast5_Summary_Type > reads;
    list< string > files;
    // initialize structs
    init_models(models);
    init_transitions(transitions);
    init_files(files);
    init_reads(models, files, reads);
    // do some training
    train_reads(models, transitions, reads);
    // basecall reads
    basecall_reads(models, transitions, reads);
    // print stats
    if (not opts::stats_fn.get().empty())
    {
        strict_fstream::ofstream ofs(opts::stats_fn);
        for (const auto& s : reads)
        {
            s.write_tsv(ofs);
            ofs << endl;
        }
    }
}

int main(int argc, char * argv[])
{
    opts::cmd_parser.parse(argc, argv);
    logger::Logger::set_default_level(logger::level::info);
    logger::Logger::set_levels_from_options(opts::log_level);
    Fast5_Summary_Type::min_read_len() = opts::min_read_len;
#ifndef H5_HAVE_THREADSAFE
    if (opts::num_threads > 1)
    {
        LOG(warning) << "enabled multi-threading with non-threadsafe HDF5: using experimental locking" << endl;
    }
#endif
    LOG(info) << "program: " << opts::cmd_parser.getProgramName() << endl;
    LOG(info) << "version: " << opts::cmd_parser.getVersion() << endl;
    LOG(info) << "args: " << opts::cmd_parser.getOrigArgv() << endl;
    real_main();
}
