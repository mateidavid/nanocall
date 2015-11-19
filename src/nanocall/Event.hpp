#ifndef __EVENT_HPP
#define __EVENT_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "fast5.hpp"
#include "logger.hpp"
#include "alg.hpp"

template < typename Float_Type = float >
class Event
{
public:
    Float_Type mean;
    Float_Type stdv;
    Float_Type log_stdv;
    Float_Type start;
    Float_Type length;
    void update_logs() { log_stdv = std::log(stdv); }
    friend std::ostream & operator << (std::ostream& os, const Event< Float_Type >& ev)
    {
        os << ev.mean << '\t'
           << ev.stdv << '\t'
           << ev.start << '\t'
           << ev.length << std::endl;
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

template < typename Float_Type = float >
using Event_Sequence = std::vector< Event< Float_Type > >;

typedef Event<> Event_Type;
typedef Event_Sequence<> Event_Sequence_Type;

#endif
