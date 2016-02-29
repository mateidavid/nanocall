#ifndef __EVENT_HPP
#define __EVENT_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "fast5.hpp"
#include "logger.hpp"
#include "alg.hpp"

template < typename Float_Type >
class Event
{
public:
    Float_Type mean;
    Float_Type stdv;
    Float_Type start;
    Float_Type length;
    Float_Type log_mean;
    Float_Type log_stdv;
    Float_Type log_start;
    void update_logs()
    {
        log_mean = std::log(mean);
        log_stdv = std::log(stdv);
        log_start = std::log(start);
    }
    friend std::ostream & operator << (std::ostream& os, const Event< Float_Type >& ev)
    {
        os << ev.mean << '\t'
           << ev.stdv << '\t'
           << ev.start << '\t'
           << ev.length;
        return os;
    }
    friend std::istream & operator >> (std::istream& is, Event< Float_Type >& ev)
    {
        is >> ev.mean
           >> ev.stdv
           >> ev.start
           >> ev.length;
        ev.update_logs();
        return is;
    }
}; // class Event

template < typename Float_Type >
struct Event_Sequence
    : std::deque< Event< Float_Type > >
{
    typedef std::deque< Event< Float_Type > > Base;
    using Base::Base;
    void apply_drift_correction(Float_Type drift)
    {
        for (auto& e : *this)
        {
            e.mean -= drift * e.start;
            e.update_logs();
        }
    }
}; // struct Event_Sequence

#endif
