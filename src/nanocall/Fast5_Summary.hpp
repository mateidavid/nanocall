#ifndef __FAST5_SUMMARY_HPP
#define __FAST5_SUMMARY_HPP

#include <array>
#include <string>
#include <vector>

#ifndef H5_HAVE_THREADSAFE
#include <mutex>
#endif

#include "Pore_Model.hpp"
#include "fast5.hpp"
#include "alg.hpp"

template < typename Float_Type = float >
class Fast5_Summary
{
public:
    typedef Pore_Model< Float_Type > Pore_Model_Type;
    typedef Pore_Model_Parameters< Float_Type > Pore_Model_Parameters_Type;
    typedef Event< Float_Type > Event_Type;
    typedef Event_Sequence< Float_Type > Event_Sequence_Type;

    std::string file_name;
    std::string base_file_name;
    std::string read_id;
    std::array< std::string, 3 > preferred_model;
    std::array< std::map< std::string, Pore_Model_Parameters_Type >, 3 > params;
    std::array< unsigned, 4 > strand_bounds;
    std::array< Float_Type, 2 > time_length;
    unsigned num_ed_events;
    float sampling_rate;
    float abasic_level;
    bool have_ed_events;
    bool valid;

    static unsigned& min_read_len()
    {
        static unsigned _min_read_len = 1000;
        return _min_read_len;
    }

    // from fast5 file
    std::vector< fast5::EventDetection_Event_Entry > ed_events;
    // filtered
    std::array< Event_Sequence_Type, 2 > events;

    Fast5_Summary() : valid(false) {}
    Fast5_Summary(const std::string fn, const Pore_Model_Dict_Type& models, bool scale_strands_together)
        : valid(false) { summarize(fn, models, scale_strands_together); }

    void summarize(const std::string& fn, const Pore_Model_Dict_Type& models, bool scale_strands_together)
    {
        valid = true;
        file_name = fn;
        auto pos = file_name.find_last_of('/');
        base_file_name = (pos != std::string::npos? file_name.substr(pos + 1) : file_name);
        strand_bounds = {{ 0, 0, 0, 0 }};
        time_length = {{ 0.0, 0.0 }};
        if (base_file_name.substr(base_file_name.size() - 6) == ".fast5")
        {
            base_file_name.resize(base_file_name.size() - 6);
        }
        fast5::File f(file_name);
        if (f.have_sampling_rate())
        {
            sampling_rate = f.get_sampling_rate();
        }
        else
        {
            LOG(warning) << fn << ": missing sampling rate; assuming 5000.0" << std::endl;
            sampling_rate = 5000.0;
        }
        have_ed_events = f.have_eventdetection_events();
        num_ed_events = 0;
        abasic_level = 0.0;
        if (have_ed_events)
        {
            load_ed_events(f);
            auto ed_params = f.get_eventdetection_event_parameters();
            num_ed_events = ed_events.size();
            read_id = ed_params.read_id;
            if (read_id.empty())
            {
                read_id = base_file_name;
            }
            abasic_level = detect_abasic_level();
            if (abasic_level > 1.0)
            {
                detect_strands();
                load_events(scale_strands_together);
                // compute time lengths
                for (unsigned st = 0; st < 2; ++st)
                {
                    if (events[st].size() < min_read_len()) continue;
                    time_length[st] = events[st].rbegin()->start + events[st].rbegin()->length;
                }
                //
                // compute initial model scalings
                //
                if (scale_strands_together
                    and events[0].size() >= min_read_len()
                    and events[1].size() >= min_read_len())
                {
                    auto r0 = alg::mean_stdv_of< Float_Type >(
                        events[0],
                        [] (const Event_Type& ev) { return ev.mean; });
                    auto r1 = alg::mean_stdv_of< Float_Type >(
                        events[1],
                        [] (const Event_Type& ev) { return ev.mean; });
                    for (const auto& p0 : models)
                        if (p0.second.strand() == 0 or p0.second.strand() == 2)
                            for (const auto& p1 : models)
                                if (p1.second.strand() == 1 or p1.second.strand() == 2)
                                {
                                    auto m_name_str = p0.first + '+' + p1.first;
                                    Pore_Model_Parameters_Type param;
                                    param.scale = (r0.second / p0.second.stdv()
                                                   + r1.second / p1.second.stdv()) / 2;
                                    param.shift = (r0.first - param.scale * p0.second.mean()
                                                   + r1.first - param.scale * p1.second.mean()) / 2;
                                    LOG("Fast5_Summary", debug)
                                        << "initial_scaling read [" << read_id << "] strand [" << 2
                                        << "] model [" << m_name_str
                                        << "] parameters [" << param << "]" << std::endl;
                                    params[2][m_name_str] = std::move(param);
                                }
                }
                else // not scale_strands_together
                {
                    for (unsigned st = 0; st < 2; ++st)
                    {
                        if (events[st].size() < min_read_len()) continue;
                        auto r = alg::mean_stdv_of< Float_Type >(
                            events[st],
                            [] (const Event_Type& ev) { return ev.mean; });
                        for (const auto& p : models)
                        {
                            if (p.second.strand() == st or p.second.strand() == 2)
                            {
                                Pore_Model_Parameters_Type param;
                                param.scale = r.second / p.second.stdv();
                                param.shift = r.first - param.scale * p.second.mean();
                                LOG("Fast5_Summary", debug)
                                    << "initial_scaling read [" << read_id << "] strand [" << st
                                    << "] model [" << p.first << "] parameters [" << param << "]" << std::endl;
                                params[st][p.first] = std::move(param);
                            }
                        }
                    }
                }
                ed_events.clear();
            } // if abasic_level > 1.0
        } // if have_ed_events
    }

    void load_events(bool scale_strands_together)
    {
        assert(valid);
        if (not have_ed_events)
        {
            return;
        }
        if (ed_events.empty())
        {
#ifndef H5_HAVE_THREADSAFE
            static std::mutex fast5_mutex;
            std::lock_guard< std::mutex > fast5_lock(fast5_mutex);
#endif
            fast5::File f(file_name);
            load_ed_events(f);
            f.close();
        }
        for (unsigned st = 0; st < 2; ++st)
        {
            events[st].clear();
            for (unsigned j = strand_bounds[2 * st]; j < strand_bounds[2 * st + 1]; ++j)
            {
                if (filter_ed_event(ed_events[j], abasic_level))
                {
                    Event_Type e;
                    e.mean = ed_events[j].mean;
                    e.stdv = ed_events[j].stdv;
                    e.start = (ed_events[j].start - ed_events[strand_bounds[scale_strands_together? 0 : 2 * st]].start) / sampling_rate;
                    e.length = ed_events[j].length / sampling_rate;
                    e.update_logs();
                    events[st].push_back(e);
                }
            }
        }
        ed_events.clear();
    }
    void drop_events()
    {
        for (unsigned st = 0; st < 2; ++st)
        {
            events[st].clear();
        }
    }

    friend std::ostream& operator << (std::ostream& os, const Fast5_Summary& fs)
    {
        os << "[base_file_name=" << fs.base_file_name << " valid=" << fs.valid;
        if (fs.valid)
        {
            os << " have_ed_events=" << fs.have_ed_events;
            if (fs.have_ed_events)
            {
                os << " read_id=" << fs.read_id
                   << " abasic_level=" << fs.abasic_level
                   << " num_ed_events=" << fs.num_ed_events
                   << " strand_bounds=[" << fs.strand_bounds[0] << ","
                   << fs.strand_bounds[1] << ","
                   << fs.strand_bounds[2] << ","
                   << fs.strand_bounds[3]
                   << "] time_length=[" << fs.time_length[0] << "," << fs.time_length[1] << "]";
            }
        }
        os << "]";
        return os;
    }

    static void write_tsv_header(std::ostream& os)
    {
        os << "file_name" << "\tread_name" << "\tnum_ed_events" << "\tabasic_level"
           << "\ttemplate_start_idx" << "\ttemplate_end_idx"
           << "\tcomplement_start_idx" << "\tcomplement_end_idx";
        for (unsigned st = 0; st < 2; ++st)
        {
            os << "\tn" << st << "_model_name"
               << "\tn" << st << "_scale"
               << "\tn" << st << "_shift"
               << "\tn" << st << "_drift"
               << "\tn" << st << "_var"
               << "\tn" << st << "_scale_sd"
               << "\tn" << st << "_var_sd";
        }
    }

    void write_tsv(std::ostream& os) const
    {
        assert(have_ed_events);
        os << base_file_name << '\t' << read_id << '\t' << num_ed_events << '\t' << abasic_level
           << '\t' << strand_bounds[0] << '\t' << strand_bounds[1]
           << '\t' << strand_bounds[2] << '\t' << strand_bounds[3];
        for (unsigned st = 0; st < 2; ++st)
        {
            os << '\t' << preferred_model[st] << '\t';
            if (not preferred_model[st].empty())
            {
                params[st].at(preferred_model[st]).write_tsv(os);
            }
            else
            {
                Pore_Model_Parameters_Type().write_tsv(os);
            }
        }
    }

private:
    void load_ed_events(fast5::File& f)
    {
        ed_events = f.get_eventdetection_events();
        auto ed_params = f.get_eventdetection_event_parameters();
    }

    // crude detection of abasic level
    Float_Type detect_abasic_level()
    {
        if (ed_events.size() < min_read_len())
        {
            return 0.0;
        }
        //
        // use 1.0 pA + max level excluding to 5%
        //
        std::vector< Float_Type > s;
        s.resize(ed_events.size());
        unsigned i;
        for (i = 0; i < ed_events.size(); ++i)
        {
            s[i] = ed_events[i].mean;
        }
        std::sort(s.begin(), s.end());
        return s[99 * s.size() / 100] + 5.0f;
    } // detect_abasic_level()

    // crude detection of abasic level
    Float_Type detect_abasic_level_2()
    {
        if (ed_events.size() < min_read_len())
        {
            return 0.0;
        }
        //
        // look for a peak level greater than the median,
        // such that the next peak below it is more than 1pA lower
        //
        std::vector< Float_Type > s;
        s.resize(ed_events.size());
        unsigned i;
        for (i = 0; i < ed_events.size(); ++i)
        {
            s[i] = ed_events[i].mean;
        }
        std::sort(s.begin(), s.end());
        i = s.size() / 2;
        while (i < s.size() and s[i - 1] > s[i] - 1.0) ++i;
        if (i >= s.size())
        {
            return 0.0;
        }
        return s[i];
    } // detect_abasic_level_2()

    // crude detection of strands in event sequence
    void detect_strands()
    {
        if (ed_events.size() < 50u)
        {
            strand_bounds = {{ 0, 0, 0, 0 }};
            return;
        }
        strand_bounds = { { 50, static_cast< unsigned >(ed_events.size() - 50), 0, 0 } };
        LOG("Fast5_Summary", debug)
            << "num_events=" << ed_events.size()
            << " abasic_level=" << abasic_level << std::endl;
        //
        // find islands of >= 5 consecutive events at high level
        //
        std::vector< std::pair< unsigned, unsigned > > islands;
        unsigned i = 0;
        while (i < ed_events.size())
        {
            if (ed_events[i].mean >= abasic_level)
            {
                unsigned j = i + 1;
                while (j < ed_events.size() and ed_events[j].mean >= abasic_level) ++j;
                if (j - i >= 5)
                {
                    islands.push_back(std::make_pair(i, j));
                    LOG("Fast5_Summary", debug) << "abasic_island [" << i << "," << j << "]" << std::endl;
                }
                i = j + 1;
            }
            else
            {
                ++i;
            }
        }
        //
        // merge islands within 50bp of each other
        //
        for (i = 1; i < islands.size(); ++i)
        {
            if (islands[i - 1].second + 50 >= islands[i].first)
            {
                LOG("Fast5_Summary", debug) << "merge_islands "
                          << "[" << islands[i - 1].first << "," << islands[i - 1].second << "] with "
                          << "[" << islands[i].first << "," << islands[i].second << "]" << std::endl;
                islands[i - 1].second = islands[i].second;
                islands.erase(islands.begin() + i);
                i = 0;
            }
        }
        LOG("Fast5_Summary", debug)
            << "final_islands: " << alg::os_join(
                islands, " ",
                [] (const std::pair< unsigned, unsigned >& p) {
                    std::ostringstream tmp;
                    tmp << "[" << p.first << "," << p.second << "]";
                    return tmp.str();
                }) << std::endl;
        if (islands.empty())
        {
            LOG("Fast5_Summary", info)
                << "template_only read_id=[" << read_id << "]" << std::endl;
            return;
        }
        //
        // pick island closest to the middle of the event sequence
        //
        auto dist_to_middle = [&] (const std::pair< unsigned, unsigned >& p) {
            return std::min((unsigned)std::abs((long)p.first - (long)ed_events.size() / 2),
                            (unsigned)std::abs((long)p.second - (long)ed_events.size() / 2));
        };
        auto it = alg::min_of(islands, dist_to_middle);
        // check island is in the middle third; if not, intepret it as template only
        if (dist_to_middle(*it) > ed_events.size() / 6)
        {
            LOG("Fast5_Summary", info)
                << "drop_read read_id=[" << read_id
                << "] islands=[" << alg::os_join(
                    islands, " ",
                    [] (const std::pair< unsigned, unsigned >& p) {
                        std::ostringstream tmp;
                        tmp << "[" << p.first << "," << p.second << "]";
                        return tmp.str();
                    })
                << "]" << std::endl;
            return;
        }
        else
        {
            LOG("Fast5_Summary", debug)
                << "hairpin_island [" << it->first << "," << it->second << "]" << std::endl;
            strand_bounds[0] = 50;
            if (islands[0].first < 100)
            {
                strand_bounds[0] = std::max(strand_bounds[0], islands[0].second);
            }
            strand_bounds[1] = it->first - 50;
            strand_bounds[2] = it->first + 50;
            strand_bounds[3] = ed_events.size() - 50;
            if (islands[islands.size() - 1].second > ed_events.size() - 100)
            {
                strand_bounds[3] = std::min(strand_bounds[3], islands[islands.size() - 1].first);
            }
        }
    } // detect_strands()

    // crude filtering of eventdetection events
    static bool filter_ed_event(const fast5::EventDetection_Event_Entry& e, Float_Type abasic_level)
    {
        if (e.mean >= abasic_level)
        {
            return false;
        }
        if (e.stdv > 4.0)
        {
            return false;
        }
        return true;
    } // filter_ed_event()
}; // struct Fast5_Summary

typedef Fast5_Summary<> Fast5_Summary_Type;

#endif
