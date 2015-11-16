#ifndef __FAST5_SUMMARY_HPP
#define __FAST5_SUMMARY_HPP

#include <array>
#include <string>
#include <vector>

#include "Pore_Model.hpp"
#include "fast5.hpp"
#include "mean_stdv.hpp"

template < typename Float_Type = float >
class Fast5_Summary
{
public:
    typedef Pore_Model< Float_Type > Pore_Model_Type;
    typedef Pore_Model_Parameters< Float_Type > Pore_Model_Parameters_Type;
    typedef Event< Float_Type > Event_Type;
    typedef Event_Sequence< Float_Type > Event_Sequence_Type;

    std::string file_name;
    std::string read_id;
    std::array< std::string, 2 > preferred_model;
    std::array< std::map< std::string, Pore_Model_Parameters_Type >, 2 > params;
    std::array< unsigned, 4 > strand_bounds;
    unsigned num_ed_events;
    float abasic_level;
    bool have_ed_events;
    bool valid;

    // from fast5 file
    std::vector< fast5::EventDetection_Event_Entry > ed_events;
    // filtered
    std::array< Event_Sequence_Type, 2 > events;

    Fast5_Summary() : valid(false) {}
    Fast5_Summary(const std::string fn, const Pore_Model_Dict_Type& models) : valid(false) { summarize(fn, models); }

    void summarize(const std::string& fn, const Pore_Model_Dict_Type& models)
    {
        valid = true;
        file_name = fn;
        fast5::File f(file_name);
        have_ed_events = f.have_eventdetection_events();
        num_ed_events = 0;
        if (have_ed_events)
        {
            load_ed_events(f);
            auto ed_params = f.get_eventdetection_event_parameters();
            num_ed_events = ed_events.size();
            read_id = ed_params.read_id;
            if (read_id.empty())
            {
                auto pos = file_name.find_last_of('/');
                read_id = file_name.substr(pos != std::string::npos? pos + 1 : 0);
            }
            abasic_level = detect_abasic_level(ed_events);
            if (abasic_level > 1.0)
            {
                strand_bounds = detect_strands(ed_events, abasic_level);
            }
            // compute initial model scalings
            load_events();
            for (unsigned st = 0; st < 2; ++st)
            {
                if (events[st].size() < 100)
                {
                    continue;
                }
                auto r = get_mean_stdv< Float_Type >(events[st], [] (const Event_Type& ev) { return ev.mean; });
                for (const auto& p : models)
                {
                    if (p.second.strand() == st or p.second.strand() == 2)
                    {
                        Pore_Model_Parameters_Type param;
                        param.shift = r.first - p.second.mean();
                        param.scale = r.second / p.second.stdv();
                        LOG("Fast5_Summary", debug)
                            << "read [" << read_id << "] strand [" << st
                            << "] model [" << p.first << "] initial parameters "
                            << param << std::endl;
                        params[st][p.first] = std::move(param);
                    }
                }
            }
            ed_events.clear();
        }
    }

    void load_events()
    {
        assert(valid and have_ed_events);
        if (ed_events.empty())
        {
            fast5::File f(file_name);
            load_ed_events(f);
        }
        for (unsigned i = 0; i < 2; ++i)
        {
            events[i].clear();
            if (strand_bounds[2 * i + 0] == 0) continue;
            for (unsigned j = strand_bounds[2 * i + 0]; j < strand_bounds[2 * i + 1]; ++j)
            {
                if (filter_ed_event(ed_events[j], abasic_level))
                {
                    Event_Type e;
                    e.mean = ed_events[j].mean;
                    e.stdv = ed_events[j].stdv;
                    e.start = 0.0;
                    e.length = 0.0;
                    e.update_logs();
                    events[i].push_back(e);
                }
            }
        }
        ed_events.clear();
    }

    friend std::ostream& operator << (std::ostream& os, const Fast5_Summary& fs)
    {
        os << "[file_name=" << fs.file_name << " valid=" << fs.valid;
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
                   << fs.strand_bounds[3] << "]";
            }
        }
        os << "]";
        return os;
    }

private:
    void load_ed_events(fast5::File& f)
    {
        ed_events = f.get_eventdetection_events();
        auto ed_params = f.get_eventdetection_event_parameters();
    }

    // crude detection of abasic level
    static Float_Type detect_abasic_level(const std::vector< fast5::EventDetection_Event_Entry >& ev)
    {
        if (ev.size() < 100)
        {
            return 0.0;
        }
        //
        // use 1.0 pA + max level excluding to 5%
        //
        std::vector< Float_Type > s;
        s.resize(ev.size());
        unsigned i;
        for (i = 0; i < ev.size(); ++i)
        {
            s[i] = ev[i].mean;
        }
        std::sort(s.begin(), s.end());
        return s[99 * s.size() / 100] + 5.0f;
    }

    // crude detection of abasic level
    static Float_Type detect_abasic_level_2(const std::vector< fast5::EventDetection_Event_Entry >& ev)
    {
        if (ev.size() < 100)
        {
            return 0.0;
        }
        //
        // look for a peak level greater than the median,
        // such that the next peak below it is more than 1pA lower
        //
        std::vector< Float_Type > s;
        s.resize(ev.size());
        unsigned i;
        for (i = 0; i < ev.size(); ++i)
        {
            s[i] = ev[i].mean;
        }
        std::sort(s.begin(), s.end());
        i = s.size() / 2;
        while (i < s.size() and s[i - 1] > s[i] - 1.0) ++i;
        if (i >= s.size())
        {
            return 0.0;
        }
        return s[i];
    }

    // crude detection of strands in event sequence
    static std::array< unsigned, 4 >
    detect_strands(const std::vector< fast5::EventDetection_Event_Entry >& ev,
                   Float_Type abasic_level)
    {
        std::array< unsigned, 4 > res;
        LOG("Fast5_Summary", debug) << "num_events=" << ev.size() << " abasic_level=" << abasic_level << std::endl;
        //
        // find islands of >= 5 consecutive events at high level
        //
        std::vector< std::pair< unsigned, unsigned > > islands;
        unsigned i = 0;
        while (i < ev.size())
        {
            if (ev[i].mean >= abasic_level)
            {
                unsigned j = i + 1;
                while (j < ev.size() and ev[j].mean >= abasic_level) ++j;
                if (j > i + 5)
                {
                    islands.push_back(std::make_pair(i, j));
                    LOG("Fast5_Summary", debug) << "found abasic island: [" << i << "," << j << "]" << std::endl;
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
                LOG("Fast5_Summary", debug) << "merging islands: "
                          << "[" << islands[i - 1].first << "," << islands[i - 1].second << "] with "
                          << "[" << islands[i].first << "," << islands[i].second << "]" << std::endl;
                islands[i - 1].second = islands[i].second;
                islands.erase(islands.begin() + i);
                i = 0;
            }
        }
        std::ostringstream tmp;
        for (i = 0; i < islands.size(); ++i)
        {
            tmp << " [" << islands[i].first << "," << islands[i].second << "]";
        }
        LOG("Fast5_Summary", debug) << "islands after merging:" << tmp.str() << std::endl;
        //
        // pick island closest to the middle of the event sequence
        //
        auto it = min_of(
            islands,
            [&] (const std::pair< unsigned, unsigned >& p) {
                return std::min(std::abs((long)p.first - (long)ev.size() / 2),
                                std::abs((long)p.second - (long)ev.size() / 2));
            });
        if (islands.size() > 0)
        {
            LOG("Fast5_Summary", debug) << "selected hairpin island: ["
                                << it->first << "," << it->second << "]" << std::endl;
            res[0] = 50;
            if (islands[0].first < 100)
            {
                res[0] = std::max(res[0], islands[0].second);
            }
            res[1] = it->first - 50;
            res[2] = it->first + 50;
            res[3] = ev.size() - 50;
            if (islands[islands.size() - 1].second > ev.size() - 100)
            {
                res[3] = std::min(res[3], islands[islands.size() - 1].first);
            }
        }
        else
        {
            res[0] = 50;
            res[1] = ev.size() - 50;
            res[2] = 0;
            res[3] = 0;
        }
        return res;
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
