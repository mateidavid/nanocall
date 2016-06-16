#include <deque>
#include <iostream>
#include <string>
#include <tclap/CmdLine.h>

#include <ctime>

#include "global_assert.hpp"
#include "version.hpp"
#include "Pore_Model.hpp"
#include "Builtin_Model.hpp"
#include "State_Transitions.hpp"
#include "Event.hpp"
#include "Fast5_Summary.hpp"
#include "Viterbi.hpp"
#include "Forward_Backward.hpp"
#include "Parameter_Trainer.hpp"
#include "logger.hpp"
#include "alg.hpp"
#include "zstr.hpp"
#include "fast5.hpp"
#include "pfor.hpp"
#include "fs_support.hpp"

using namespace std;

long get_cpu_time_ms()
{
    auto t = clock();
    return (t * 1000) / CLOCKS_PER_SEC;
}

#ifndef FLOAT_TYPE
#define FLOAT_TYPE float
#endif
#ifndef KMER_SIZE
#define KMER_SIZE 6
#endif
typedef State_Transitions< FLOAT_TYPE, KMER_SIZE > State_Transitions_Type;
typedef State_Transition_Parameters< FLOAT_TYPE > State_Transition_Parameters_Type;
typedef Pore_Model< FLOAT_TYPE, KMER_SIZE > Pore_Model_Type;
typedef Pore_Model_Dict< FLOAT_TYPE, KMER_SIZE > Pore_Model_Dict_Type;
typedef Pore_Model_Parameters< FLOAT_TYPE > Pore_Model_Parameters_Type;
typedef Event< FLOAT_TYPE, KMER_SIZE > Event_Type;
typedef Event_Sequence< FLOAT_TYPE, KMER_SIZE > Event_Sequence_Type;
typedef Fast5_Summary< FLOAT_TYPE, KMER_SIZE > Fast5_Summary_Type;
typedef Parameter_Trainer< FLOAT_TYPE, KMER_SIZE > Parameter_Trainer_Type;
typedef Viterbi< FLOAT_TYPE, KMER_SIZE > Viterbi_Type;

namespace opts
{
    using namespace TCLAP;
    string description = "Call bases in Oxford Nanopore reads.";
    CmdLine cmd_parser(description, ' ', package_version);
    //
    ValueArg< string > ed_group("", "ed-group", "EventDetection group to use. (default: smallest available)", false, "", "000|001|...", cmd_parser);
    ValueArg< unsigned > chunk_size("", "chunk-size", "Thread chunk size.", false, 1, "int", cmd_parser);
    MultiArg< string > log_level("", "log", "Log level. (default: info)", false, "string", cmd_parser);
    ValueArg< string > stats_fn("", "stats", "Stats.", false, "", "file", cmd_parser);
    ValueArg< unsigned > max_ed_events("", "max-ed-events", "Maximum EventDetection events.", false, 100000, "int", cmd_parser);
    ValueArg< unsigned > min_ed_events("", "min-ed-events", "Minimum EventDetection events.", false, 10, "int", cmd_parser);
    ValueArg< unsigned > fasta_line_width("", "fasta-line-width", "Maximum fasta line width.", false, 80, "int", cmd_parser);
    //
    ValueArg< float > scaling_select_threshold("", "scaling-select-threshold", "Select best model per strand during scaling if log score better by threshold.", false, 20.0, "float", cmd_parser);
    ValueArg< float > scaling_min_progress("", "scaling-min-progress", "Minimum scaling fit progress.", false, 1.0, "float", cmd_parser);
    ValueArg< unsigned > scaling_max_rounds("", "scaling-max-rounds", "Maximum scaling rounds.", false, 10, "int", cmd_parser);
    ValueArg< unsigned > scaling_num_events("", "scaling-num-events", "Number of events used for model scaling.", false, 200, "int", cmd_parser);
    //
    SwitchArg single_strand_scaling("", "single-strand-scaling", "Train scaling parameters per strand.", cmd_parser);
    SwitchArg double_strand_scaling("", "double-strand-scaling", "Train scaling parameters per read. (default)", cmd_parser);
    SwitchArg no_train_transitions("", "no-train-transitions", "Do not train state transitions.", cmd_parser);
    SwitchArg no_train_scaling("", "no-train-scaling", "Do not train pore model scaling.", cmd_parser);
    SwitchArg only_train("", "only-train", "Stop after training.", cmd_parser);
    SwitchArg train("", "train", "Enable training. (default)", cmd_parser);
    SwitchArg no_train("", "no-train", "Disable all training.", cmd_parser);
    //
    ValueArg< float > pr_skip("", "pr-skip", "Transition probability of skipping at least 1 state.", false, .3, "float", cmd_parser);
    ValueArg< float > pr_stay("", "pr-stay", "Transition probability of staying in the same state.", false, .1, "float", cmd_parser);
    ValueArg< string > trans_fn("s", "trans", "Custom initial state transitions.", false, "", "file", cmd_parser);
    ValueArg< string > model_fofn("", "model-fofn", "File of pore models.", false, "", "file", cmd_parser);
    MultiArg< string > model_fn("m", "model", "Custom pore model for strand (0=template, 1=complement, 2=both).", false, "strand:file", cmd_parser);
    //
    ValueArg< string > pore("", "pore", "Pore name, used to select builtin pore model.", false, "r9", "r73|r9", cmd_parser);
    SwitchArg write_fast5("", "write-fast5", "Write basecalls to fast5 files.", cmd_parser);
    ValueArg< string > output_fn("o", "output", "Output.", false, "", "file", cmd_parser);
    ValueArg< unsigned > num_threads("t", "threads", "Number of parallel threads.", false, 1, "int", cmd_parser);
    UnlabeledMultiArg< string > input_fn("inputs", "Inputs: directories, fast5 files, or files of fast5 file names (use \"-\" to read fofn from stdin).", true, "path", cmd_parser);
} // namespace opts

void init_models(Pore_Model_Dict_Type& models)
{
    auto parse_model_name = [] (const string& s) {
        if (s.size() < 3
            or (s[0] != '0' and s[0] != '1' and s[0] != '2')
            or s[1] != ':')
        {
            LOG(error) << "could not parse model name: \"" << s << "\"; format should be \"[0|1|2]:<file>\"" << endl;
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
        LOG(error) << "models were specified only for strand " << (int)model_list[0].empty()
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
                LOG(info) << "loaded module [" << pm_name
                          << "] for strand [" << st
                          << "] statistics [mean=" << pm.mean()
                          << ", stdv=" << pm.stdv() << "]" << endl;
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
            if (pm_name.compare(0, opts::pore.get().size() + 1, opts::pore.get() + ".")) continue;
            pm.load_from_vector(Builtin_Model::init_lists[i]);
            pm.strand() = Builtin_Model::strands[i];
            models[Builtin_Model::names[i]] = move(pm);
            LOG(info)
                << "loaded builtin module [" << Builtin_Model::names[i]
                << "] for strand [" << Builtin_Model::strands[i]
                << "] statistics [mean=" << pm.mean()
                << ", stdv=" << pm.stdv() << "]" << endl;
        }
        if (models.empty())
        {
            LOG(error)
                << "no builtin models found for pore [" << opts::pore.get() << "]" << endl;
            exit(EXIT_FAILURE);
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
        transitions.compute_transitions_fast(opts::pr_skip, opts::pr_stay);
        LOG(info) << "init_state_transitions pr_skip=[" << opts::pr_skip
                  << "], pr_stay=[" << opts::pr_stay << "]" << endl;
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
    if (files.empty())
    {
        LOG(error) << "no fast5 files to process" << endl;
        exit(EXIT_FAILURE);
    }
} // init_files

void init_reads(const Pore_Model_Dict_Type& models,
                const list< string >& files,
                deque< Fast5_Summary_Type >& reads)
{
    for (const auto& f : files)
    {
        Fast5_Summary_Type s(f, models, opts::double_strand_scaling);
        LOG(info) << "summary: " << s << endl;
        reads.emplace_back(move(s));
    }
} // init_reads

void train_reads(const Pore_Model_Dict_Type& models,
                 const State_Transitions_Type& default_transitions,
                 deque< Fast5_Summary_Type >& reads)
{
    auto time_start_ms = get_cpu_time_ms();
    Parameter_Trainer_Type::init();
    unsigned crt_idx = 0;
    pfor::pfor< unsigned >(
        opts::num_threads,
        opts::chunk_size,
        // get_item
        [&] (unsigned& i) {
            if (crt_idx >= reads.size()) return false;
            i = crt_idx++;
            return true;
        },
        // process item
        [&] (unsigned& i) {
            Fast5_Summary_Type& read_summary = reads[i];
            if (read_summary.num_ed_events == 0) return;
            global_assert::global_msg() = read_summary.read_id;
            read_summary.load_events();
            //
            // create per-strand list of models to try
            //
            array< list< string >, 2 > model_list;
            for (unsigned st = 0; st < 2; ++st)
            {
                // if not enough events, ignore strand
                if (read_summary.events(st).size() < opts::min_ed_events) continue;
                // create list of models to try
                if (not read_summary.preferred_model[st][st].empty())
                {
                    // if we have a preferred model, use that
                    model_list[st].push_back(read_summary.preferred_model[st][st]);
                }
                else
                {
                    // no preferred model, try all that apply to this strand
                    for (const auto& p : models)
                    {
                        if (p.second.strand() == st or p.second.strand() == 2)
                        {
                            model_list[st].push_back(p.first);
                        }
                    }
                }
                ASSERT(not model_list.empty());
            }
            //
            // create per-strand list of event sequences on which to train
            //
            array< vector< Event_Sequence_Type >, 2 > train_event_seqs;
            for (unsigned st = 0; st < 2; ++st)
            {
                // if not enough events, ignore strand
                if (read_summary.events(st).size() < opts::min_ed_events) continue;
                // create 2 event sequences on which to train
                unsigned num_train_events = min((size_t)opts::scaling_num_events.get(), read_summary.events(st).size());
                train_event_seqs[st].emplace_back(
                    read_summary.events(st).begin(), read_summary.events(st).begin() + num_train_events / 2);
                train_event_seqs[st].emplace_back(
                    read_summary.events(st).end() - num_train_events / 2, read_summary.events(st).end());
            }
            //
            // branch on whether pore models should be scaled together
            //
            if (read_summary.scale_strands_together)
            {
                // prepare vector of event sequences
                vector< pair< const Event_Sequence_Type*, unsigned > > train_event_seq_ptrs;
                for (unsigned st = 0; st < 2; ++st)
                {
                    for (const auto& events : train_event_seqs[st])
                    {
                        train_event_seq_ptrs.push_back(make_pair(&events, st));
                    }
                }
                // track model fit
                // key = pore model name; value = fit
                map< array< string, 2 >, FLOAT_TYPE > model_fit;
                for (const auto& m_name_0 : model_list[0])
                {
                    for (const auto& m_name_1 : model_list[1])
                    {
                        array< string, 2 > m_name_key = {{ m_name_0, m_name_1 }};
                        string m_name = m_name_0 + "+" + m_name_1;
                        unsigned round = 0;
                        auto& crt_pm_params = read_summary.pm_params_m.at(m_name_key);
                        auto& crt_st_params = read_summary.st_params_m.at(m_name_key);
                        auto& crt_fit = model_fit[m_name_key];
                        crt_fit = -INFINITY;
                        while (true)
                        {
                            Pore_Model_Parameters_Type old_pm_params(crt_pm_params);
                            std::array< State_Transition_Parameters_Type, 2 > old_st_params(crt_st_params);
                            auto old_fit = crt_fit;
                            bool done;

                            Parameter_Trainer_Type::train_one_round(
                                train_event_seq_ptrs,
                                {{ &models.at(m_name_0), &models.at(m_name_1) }},
                                default_transitions,
                                old_pm_params, old_st_params,
                                crt_pm_params, crt_st_params, crt_fit, done,
                                not opts::no_train_scaling, not opts::no_train_transitions);

                            LOG(debug)
                                << "scaling_round read [" << read_summary.read_id
                                << "] strand [" << 2
                                << "] model [" << m_name
                                << "] old_pm_params [" << old_pm_params
                                << "] old_st_params [" << old_st_params[0] << "," << old_st_params[1]
                                << "] old_fit [" << old_fit
                                << "] crt_pm_params [" << crt_pm_params
                                << "] crt_st_params [" << crt_st_params[0] << "," << crt_st_params[1]
                                << "] crt_fit [" << crt_fit
                                << "] round [" << round << "]" << endl;

                            if (done)
                            {
                                // singularity detected; stop
                                break;
                            }

                            if (crt_fit < old_fit)
                            {
                                LOG(info) << "scaling_regression read [" << read_summary.read_id
                                          << "] strand [" << 2
                                          << "] model [" << m_name
                                          << "] old_params [" << old_pm_params
                                          << "] old_st_params [" << old_st_params[0] << "," << old_st_params[1]
                                          << "] old_fit [" << old_fit
                                          << "] crt_pm_params [" << crt_pm_params
                                          << "] crt_st_params [" << crt_st_params[0] << "," << crt_st_params[1]
                                          << "] crt_fit [" << crt_fit
                                          << "] round [" << round << "]" << endl;
                                crt_pm_params = old_pm_params;
                                crt_st_params = old_st_params;
                                crt_fit = old_fit;
                                break;
                            }

                            ++round;
                            // stop condition
                            if (round >= 2u * opts::scaling_max_rounds
                                or (round > 1 and crt_fit < old_fit + opts::scaling_min_progress))
                            {
                                break;
                            }

                        }; // while true
                        LOG(info)
                            << "scaling_result read [" << read_summary.read_id
                            << "] strand [" << 2
                            << "] model [" << m_name
                            << "] pm_params [" << crt_pm_params
                            << "] st_params [" << crt_st_params[0] << "," << crt_st_params[1]
                            << "] fit [" << crt_fit
                            << "] rounds [" << round << "]" << endl;
                    } // for m_name[1]
                } // for m_name[0]
                if (opts::scaling_select_threshold.get() < INFINITY)
                {
                    auto it_max = alg::max_of(
                        model_fit,
                        [] (const decltype(model_fit)::value_type& p) { return p.second; });
                    // check maximum is unique
                    if (alg::all_of(
                            model_fit,
                            [&] (const decltype(model_fit)::value_type& p) {
                                return &p == &*it_max
                                    or p.second + opts::scaling_select_threshold.get() < it_max->second;
                            }))
                    {
                        const auto& m_name_0 = it_max->first[0];
                        const auto& m_name_1 = it_max->first[1];
                        auto m_name = m_name_0 + '+' + m_name_1;
                        read_summary.preferred_model[2][0] = m_name_0;
                        read_summary.preferred_model[2][1] = m_name_1;
                        LOG(info)
                            << "selected_model read [" << read_summary.read_id
                            << "] strand [2] model [" << m_name << "]" << endl;
                    }
                }
            }
            else // not scale_strands_together
            {
                for (unsigned st = 0; st < 2; ++st)
                {
                    // if not enough events, ignore strand
                    if (read_summary.events(st).size() < opts::min_ed_events) continue;
                    // prepare vector of event sequences
                    vector< pair< const Event_Sequence_Type*, unsigned > > train_event_seq_ptrs;
                    for (const auto& events : train_event_seqs[st])
                    {
                        train_event_seq_ptrs.push_back(make_pair(&events, st));
                    }
                    map< string, FLOAT_TYPE > model_fit;
                    for (const auto& m_name : model_list[st])
                    {
                        array< string, 2 > m_name_key;
                        m_name_key[st] = m_name;
                        unsigned round = 0;
                        auto& crt_pm_params = read_summary.pm_params_m.at(m_name_key);
                        auto& crt_st_params = read_summary.st_params_m.at(m_name_key);
                        auto& crt_fit = model_fit[m_name];
                        crt_fit = -INFINITY;
                        while (true)
                        {
                            Pore_Model_Parameters_Type old_pm_params(crt_pm_params);
                            array< State_Transition_Parameters_Type, 2 > old_st_params(crt_st_params);
                            auto old_fit = crt_fit;
                            bool done;

                            Parameter_Trainer_Type::train_one_round(
                                train_event_seq_ptrs,
                                {{ &models.at(m_name), &models.at(m_name) }},
                                default_transitions,
                                old_pm_params, old_st_params,
                                crt_pm_params, crt_st_params, crt_fit, done,
                                not opts::no_train_scaling, not opts::no_train_transitions);

                            LOG(debug)
                                << "scaling_round read [" << read_summary.read_id
                                << "] strand [" << st
                                << "] model [" << m_name
                                << "] old_pm_params [" << old_pm_params
                                << "] old_st_params [" << old_st_params[st]
                                << "] old_fit [" << old_fit
                                << "] crt_pm_params [" << crt_pm_params
                                << "] crt_st_params [" << crt_st_params[st]
                                << "] crt_fit [" << crt_fit
                                << "] round [" << round << "]" << endl;

                            if (done)
                            {
                                // singularity detected; stop
                                break;
                            }

                            if (crt_fit < old_fit)
                            {
                                LOG(info) << "scaling_regression read [" << read_summary.read_id
                                          << "] strand [" << st
                                          << "] model [" << m_name
                                          << "] old_pm_params [" << old_pm_params
                                          << "] old_st_params [" << old_st_params[st]
                                          << "] old_fit [" << old_fit
                                          << "] crt_pm_params [" << crt_pm_params
                                          << "] crt_st_params [" << crt_st_params[st]
                                          << "] crt_fit [" << crt_fit
                                          << "] round [" << round << "]" << endl;
                                crt_pm_params = old_pm_params;
                                crt_st_params = old_st_params;
                                crt_fit = old_fit;
                                break;
                            }

                            ++round;
                            // stop condition
                            if (round >= opts::scaling_max_rounds
                                or (round > 1 and crt_fit < old_fit + opts::scaling_min_progress))
                            {
                                break;
                            }

                        }; // while true
                        LOG(info)
                            << "scaling_result read [" << read_summary.read_id
                            << "] strand [" << st
                            << "] model [" << m_name
                            << "] pm_params [" << crt_pm_params
                            << "] st_params [" << crt_st_params[st]
                            << "] fit [" << crt_fit
                            << "] rounds [" << round << "]" << endl;
                    } // for m_name
                    if (opts::scaling_select_threshold.get() < INFINITY)
                    {
                        auto it_max = alg::max_of(
                            model_fit,
                            [] (const decltype(model_fit)::value_type& p) { return p.second; });
                        if (alg::all_of(
                                model_fit,
                                [&] (const decltype(model_fit)::value_type& p) {
                                    return &p == &*it_max
                                        or p.second + opts::scaling_select_threshold.get() < it_max->second;
                                }))
                        {
                            read_summary.preferred_model[st][st] = it_max->first;
                            LOG(info)
                                << "selected_model read [" << read_summary.read_id
                                << "] strand [" << st
                                << "] model [" << it_max->first << "]" << endl;
                        }
                    }
                } // for st
            } // if not scale_strands_together
            read_summary.drop_events();
        }, // process_item
        // progress_report
        [&] (unsigned items, unsigned seconds) {
            clog << "Processed " << setw(6) << right << items << " reads in "
                 << setw(6) << right << seconds << " seconds\r";
        }); // pfor
    auto time_end_ms = get_cpu_time_ms();
    LOG(info) << "training user_cpu_secs=" << (time_end_ms - time_start_ms)/1000 << endl;
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
                    const State_Transitions_Type& default_transitions,
                    deque< Fast5_Summary_Type >& reads)
{
    auto time_start_ms = get_cpu_time_ms();
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
    pfor::pfor< unsigned, ostringstream >(
        opts::num_threads,
        opts::chunk_size,
        // get_item
        [&] (unsigned& i) {
            if (crt_idx >= reads.size()) return false;
            i = crt_idx++;
            return true;
        },
        // process_item
        [&] (unsigned& i, ostringstream& oss) {
            Fast5_Summary_Type& read_summary = reads[i];
            if (read_summary.num_ed_events == 0) return;
            global_assert::global_msg() = read_summary.read_id;
            read_summary.load_events();

            // compute read statistics used to check scaling
            array< pair< FLOAT_TYPE, FLOAT_TYPE >, 2 > r_stats;
            for (unsigned st = 0; st < 2; ++st)
            {
                // if not enough events, ignore strand
                if (read_summary.events(st).size() < opts::min_ed_events) continue;
                r_stats[st] = alg::mean_stdv_of< FLOAT_TYPE >(
                    read_summary.events(st),
                    [] (const Event_Type& ev) { return ev.mean; });
                LOG(debug)
                    << "mean_stdv read [" << read_summary.read_id
                    << "] strand [" << st
                    << "] ev_mean=[" << r_stats[st].first
                    << "] ev_stdv=[" << r_stats[st].second << "]" << endl;
            }

            // basecalling functor
            // returns: (path_prob, base_seq)
            auto basecall_strand = [&] (unsigned st, string m_name,
                                        const Pore_Model_Parameters_Type& pm_params,
                                        const State_Transition_Parameters_Type& st_params) {
                // scale model
                Pore_Model_Type pm(models.at(m_name));
                pm.scale(pm_params);
                State_Transitions_Type custom_transitions;
                const State_Transitions_Type* transitions_ptr;
                if (not st_params.is_default())
                {
                    custom_transitions.compute_transitions_fast(st_params);
                    transitions_ptr = &custom_transitions;
                }
                else
                {
                    transitions_ptr = &default_transitions;
                }
                LOG(info)
                    << "basecalling read [" << read_summary.read_id
                    << "] strand [" << st
                    << "] model [" << m_name
                    << "] pm_params [" << pm_params
                    << "] st_params [" << st_params << "]" << endl;
                LOG(debug)
                    << "mean_stdv read [" << read_summary.read_id
                    << "] strand [" << st
                    << "] model_mean [" << pm.mean()
                    << "] model_stdv [" << pm.stdv() << "]" << endl;
                if (abs(r_stats[st].first - pm.mean()) > 5.0)
                {
                    LOG(warning)
                        << "means_apart read [" << read_summary.read_id
                        << "] strand [" << st
                        << "] model [" << m_name
                        << "] parameters [" << pm_params
                        << "] model_mean=[" << pm.mean()
                        << "] events_mean=[" << r_stats[st].first
                        << "]" << endl;
                }
                // correct drift
                Event_Sequence_Type corrected_events = read_summary.events(st);
                corrected_events.apply_drift_correction(pm_params.drift);
                Viterbi_Type vit;
                vit.fill(pm, *transitions_ptr, corrected_events);
                return std::make_tuple(vit.path_probability(), std::move(corrected_events));
            };

            if (read_summary.scale_strands_together)
            {
                // create list of models to try
                list< array< string, 2 > > model_sublist;
                if (not read_summary.preferred_model[2][0].empty())
                {
                    // if we have a preferred model, use that
                    model_sublist.push_back(read_summary.preferred_model[2]);
                }
                else
                {
                    // no preferred model, try all for which we have scaling parameters
                    for (const auto& p : read_summary.pm_params_m)
                    {
                        if (p.first[0].empty() or p.first[1].empty()) continue;
                        model_sublist.push_back(p.first);
                    }
                }
                // basecall using applicable models
                deque< tuple< FLOAT_TYPE,
                              FLOAT_TYPE, FLOAT_TYPE,
                              string, string,
                              Event_Sequence_Type, Event_Sequence_Type > > results;
                for (const auto& m_name : model_sublist)
                {
                    array< tuple< FLOAT_TYPE, Event_Sequence_Type >, 2 > part_results;
                    for (unsigned st = 0; st < 2; ++st)
                    {
                        part_results[st] = basecall_strand(
                            st, m_name[st],
                            read_summary.pm_params_m.at(m_name),
                            read_summary.st_params_m.at(m_name)[st]);
                    }
                    results.emplace_back(get<0>(part_results[0]) + get<0>(part_results[1]),
                                         get<0>(part_results[0]),
                                         get<0>(part_results[1]),
                                         string(m_name[0]),
                                         string(m_name[1]),
                                         std::move(get<1>(part_results[0])),
                                         std::move(get<1>(part_results[1])));
                }
                // sort results by first component (log path probability)
                sort(results.begin(),
                     results.end(),
                     [] (const decltype(results)::value_type& lhs, const decltype(results)::value_type& rhs) {
                         return get<0>(lhs) < get<0>(rhs);
                     });
                array< FLOAT_TYPE, 2 > best_log_path_prob{{ get<1>(results.back()), get<2>(results.back()) }};
                array< string, 2 > best_m_name{{ get<3>(results.back()), get<4>(results.back()) }};
                array< const Event_Sequence_Type*, 2 > event_seq_ptr = {
                    &get<5>(results.back()),
                    &get<6>(results.back())
                };
                array< string, 2 > base_seq = {
                    get<5>(results.back()).get_base_seq(),
                    get<6>(results.back()).get_base_seq()
                };
                string best_m_name_str = best_m_name[0] + '+' + best_m_name[1];
                auto& best_pm_params = read_summary.pm_params_m.at(best_m_name);
                auto& best_st_params = read_summary.st_params_m.at(best_m_name);
                for (unsigned st = 0; st < 2; ++st)
                {
                    LOG(info)
                        << "best_model read [" << read_summary.read_id
                        << "] strand [" << st
                        << "] model [" << best_m_name[st]
                        << "] pm_params [" << best_pm_params
                        << "] st_params [" << best_st_params[st]
                        << "] log_path_prob [" << best_log_path_prob[st] << "]" << endl;
                    read_summary.preferred_model[st][st] = best_m_name[st];
                    read_summary.pm_params_m[read_summary.preferred_model[st]] = best_pm_params;
                    read_summary.st_params_m[read_summary.preferred_model[st]][st] = best_st_params[st];
                    string seq_name;
                    {
                        ostringstream tmp;
                        tmp << read_summary.read_id << ":" << read_summary.base_file_name << ":" << st;
                        seq_name = tmp.str();
                    }
                    if (opts::write_fast5)
                    {
                        read_summary.add_basecall_seq(seq_name, st, base_seq[st]);
                        read_summary.add_basecall_events(st, *event_seq_ptr[st]);
                        read_summary.add_basecall_model(st, models.at(best_m_name[st]));
                        read_summary.add_basecall_model_params(st, best_pm_params);
                    }
                    else
                    {
                        write_fasta(oss, seq_name, base_seq[st]);
                    }
                }
            }
            else // not scale_strands_together
            {
                for (unsigned st = 0; st < 2; ++st)
                {
                    // if not enough events, ignore strand
                    if (read_summary.events(st).size() < opts::min_ed_events) continue;
                    // create list of models to try
                    list< array< string, 2 > > model_sublist;
                    if (not read_summary.preferred_model[st][st].empty())
                    {
                        // if we have a preferred model, use that
                        model_sublist.push_back(read_summary.preferred_model[st]);
                    }
                    else
                    {
                        // no preferred model, try all for which we have scaling
                        for (const auto& p : read_summary.pm_params_m)
                        {
                            if (not p.first[st].empty() and p.first[1 - st].empty())
                            {
                                model_sublist.push_back(p.first);
                            }
                        }
                    }
                    // deque of results
                    deque< tuple< FLOAT_TYPE, string, Event_Sequence_Type > > results;
                    for (const auto& m_name : model_sublist)
                    {
                        auto r = basecall_strand(
                            st, m_name[st],
                            read_summary.pm_params_m.at(m_name),
                            read_summary.st_params_m.at(m_name)[st]);
                        results.emplace_back(get<0>(r),
                                             string(m_name[st]),
                                             std::move(get<1>(r)));
                    }
                    sort(results.begin(),
                         results.end(),
                         [] (const decltype(results)::value_type& lhs, const decltype(results)::value_type& rhs) {
                             return get<0>(lhs) < get<0>(rhs);
                         });
                    const string& best_m_name = get<1>(results.back());
                    const Event_Sequence_Type& event_seq = get<2>(results.back());
                    string base_seq = event_seq.get_base_seq();
                    array< string, 2 > best_m_key;
                    best_m_key[st] = best_m_name;
                    LOG(info)
                        << "best_model read [" << read_summary.read_id
                        << "] strand [" << st
                        << "] model [" << best_m_name
                        << "] pm_params [" << read_summary.pm_params_m.at(best_m_key)
                        << "] st_params [" << read_summary.st_params_m.at(best_m_key)[st]
                        << "] log_path_prob [" << get<0>(results.back()) << "]" << endl;
                    read_summary.preferred_model[st][st] = best_m_name;
                    string seq_name;
                    {
                        ostringstream tmp;
                        tmp << read_summary.read_id << ":" << read_summary.base_file_name << ":" << st;
                        seq_name = tmp.str();
                    }
                    if (opts::write_fast5)
                    {
                        read_summary.add_basecall_seq(seq_name, st, base_seq);
                        read_summary.add_basecall_events(st, event_seq);
                        read_summary.add_basecall_model(st, models.at(best_m_name));
                        read_summary.add_basecall_model_params(st, read_summary.pm_params_m.at(best_m_key));
                    }
                    else
                    {
                        write_fasta(oss, seq_name, base_seq);
                    }
                } // for st
            }
            read_summary.drop_events();
        },
        // output_chunk
        [&] (ostringstream& oss) {
            *os_p << oss.str();
        },
        // progress_report
        [&] (unsigned items, unsigned seconds) {
            clog << "Processed " << setw(6) << right << items << " reads in "
                 << setw(6) << right << seconds << " seconds\r";
        }); // pfor
    auto time_end_ms = get_cpu_time_ms();
    LOG(info) << "basecalling user_cpu_secs=" << (time_end_ms - time_start_ms)/1000 << endl;
} // basecall_reads

int real_main()
{
    Pore_Model_Dict_Type models;
    State_Transitions_Type default_transitions;
    deque< Fast5_Summary_Type > reads;
    list< string > files;
    // initialize structs
    init_models(models);
    init_transitions(default_transitions);
    init_files(files);
    init_reads(models, files, reads);
    if (opts::train)
    {
        // do some rescaling
        train_reads(models, default_transitions, reads);
    }
    if (not opts::only_train)
    {
        // basecall reads
        basecall_reads(models, default_transitions, reads);
    }
    // print stats
    if (not opts::stats_fn.get().empty())
    {
        strict_fstream::ofstream ofs(opts::stats_fn);
        Fast5_Summary_Type::write_tsv_header(ofs);
        ofs << endl;
        for (const auto& s : reads)
        {
            s.write_tsv(ofs);
            ofs << endl;
        }
    }
    assert(fast5::File::get_object_count() == 0);
    return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
    opts::cmd_parser.parse(argc, argv);
    logger::Logger::set_default_level(logger::level::info);
    logger::Logger::set_levels_from_options(opts::log_level);
    LOG(info) << "program: " << opts::cmd_parser.getProgramName() << endl;
    LOG(info) << "version: " << opts::cmd_parser.getVersion() << endl;
    LOG(info) << "args: " << opts::cmd_parser.getOrigArgv() << endl;
    LOG(info) << "num_threads=" << opts::num_threads.get() << endl;
#ifndef H5_HAVE_THREADSAFE
    if (opts::num_threads > 1)
    {
        LOG(warning) << "enabled multi-threading with non-threadsafe HDF5: using experimental locking" << endl;
    }
#endif
    State_Transition_Parameters_Type::default_p_stay() = opts::pr_stay;
    State_Transition_Parameters_Type::default_p_skip() = opts::pr_skip;
    Fast5_Summary_Type::min_ed_events() = opts::min_ed_events;
    Fast5_Summary_Type::max_ed_events() = opts::max_ed_events;
    Fast5_Summary_Type::eventdetection_group() = opts::ed_group;
    //
    // set training option
    //
    if (opts::train and opts::no_train)
    {
        LOG(error)
            << "either --train or --no-train may be used, but not both" << endl;
        return EXIT_FAILURE;
    }
    else if (not opts::train and not opts::no_train)
    {
        // by default, enable training
        opts::train.set(true);
    }
    ASSERT(opts::train != opts::no_train);
    //
    // set single/double strand scaling option
    //
    if (opts::train and not opts::no_train_scaling)
    {
        if (opts::single_strand_scaling and opts::double_strand_scaling)
        {
            LOG(error)
                << "either --single-strand-scaling or --double-strand-scaling may be used, but not both" << endl;
            return EXIT_FAILURE;
        }
        else if (not opts::single_strand_scaling and not opts::double_strand_scaling)
        {
            // by default, do double strand scaling
            opts::double_strand_scaling.set(true);
        }
    }
    //
    // check other options
    //
    if (opts::scaling_select_threshold.get() < 0.0)
    {
        LOG(error)
            << "invalid scaling_select_threshold: " << opts::scaling_select_threshold.get() << endl;
        return EXIT_FAILURE;
    }
    if (opts::scaling_min_progress < 0.0)
    {
        LOG(error)
            << "invalid scaling_min_progress: " << opts::scaling_min_progress.get() << endl;
        return EXIT_FAILURE;
    }
    if (not opts::output_fn.get().empty() and opts::write_fast5)
    {
        LOG(error)
            << "output may be written to fast5 files or to a single output file, but not both" << endl;
        return EXIT_FAILURE;
    }
    //
    // print training options
    //
    LOG(info) << "train=" << opts::train.get() << endl;
    if (opts::train)
    {
        LOG(info) << "only_train=" << opts::only_train.get() << endl;
        LOG(info) << "train_scaling=" << not opts::no_train_scaling.get() << endl;
        LOG(info) << "train_transitions=" << not opts::no_train_transitions.get() << endl;
        if (not opts::no_train_scaling)
        {
            LOG(info) << "double_strands_scaling=" << opts::double_strand_scaling.get() << endl;
            LOG(info) << "scaling_num_events=" << opts::scaling_num_events.get() << endl;
            LOG(info) << "scaling_max_rounds=" << opts::scaling_max_rounds.get() << endl;
            LOG(info) << "scaling_min_progress=" << opts::scaling_min_progress.get() << endl;
            LOG(info) << "scaling_select_threshold=" << opts::scaling_select_threshold.get() << endl;
        }
    }
    return real_main();
}
