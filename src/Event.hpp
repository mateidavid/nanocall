#ifndef __EVENTS_HPP
#define __EVENTS_HPP

#include <cmath>
#include <iostream>

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

#endif
