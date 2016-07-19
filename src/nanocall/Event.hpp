#ifndef __EVENT_HPP
#define __EVENT_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "fast5.hpp"
#include "logger.hpp"
#include "alg.hpp"

template < typename Float_Type, unsigned Kmer_Size >
class Event
{
public:
    Float_Type mean;
    Float_Type corrected_mean;
    Float_Type stdv;
    Float_Type start;
    Float_Type length;
    Float_Type log_mean;
    Float_Type log_corrected_mean;
    Float_Type log_stdv;
    Float_Type log_start;
    //
    Float_Type orig_mean;
    Float_Type p_model_state;
    std::array< char, Kmer_Size > model_state;
    unsigned model_state_idx;
    int move;
    //
    void update_logs()
    {
        log_mean = std::log(mean);
        log_corrected_mean = std::log(corrected_mean);
        log_stdv = std::log(stdv);
        log_start = std::log(start);
    }
    void set_model_state(const std::string& s)
    {
        assert(s.size() == Kmer_Size);
        std::copy_n(s.begin(), Kmer_Size, model_state.begin());
    }
    friend std::ostream & operator << (std::ostream& os, const Event& ev)
    {
        os << ev.mean << '\t'
           << ev.stdv << '\t'
           << ev.start << '\t'
           << ev.length;
        return os;
    }
    friend std::istream & operator >> (std::istream& is, Event& ev)
    {
        is >> ev.mean
           >> ev.stdv
           >> ev.start
           >> ev.length;
        ev.corrected_mean = ev.mean;
        ev.update_logs();
        return is;
    }
}; // class Event

template < typename Float_Type, unsigned Kmer_Size >
struct Event_Sequence
    : std::vector< Event< Float_Type, Kmer_Size > >
{
    typedef std::vector< Event< Float_Type, Kmer_Size > > Base;
    using Base::Base;
    void apply_drift_correction(Float_Type drift)
    {
        for (auto& e : *this)
        {
            e.corrected_mean -= drift * e.start;
            e.log_corrected_mean = std::log(e.corrected_mean);
        }
    }
    std::string get_base_seq() const
    {
        std::string res;
        const Base& v = *this;
        res.assign(v[0].model_state.begin(), v[0].model_state.end());
        for (unsigned i = 1; i < v.size(); ++i)
        {
            unsigned a = std::min((unsigned)v[i].move, (unsigned)Kmer_Size);
            unsigned b = Kmer_Size - a;
            assert(std::string(v[i - 1].model_state.begin() + a, v[i - 1].model_state.end())
                   == std::string(v[i].model_state.begin(), v[i].model_state.begin() + b));
            res += std::string(v[i].model_state.begin() + b, v[i].model_state.end());
        }
        return res;
    }
}; // struct Event_Sequence

#endif
