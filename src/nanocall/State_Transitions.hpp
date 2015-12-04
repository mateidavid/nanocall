#ifndef __STATE_TRANSITIONS_BASE_HPP
#define __STATE_TRANSITIONS_BASE_HPP

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>

#include "Kmer.hpp"
#include "logsumset.hpp"

template < typename Float_Type = float >
struct State_Neighbours
{
    State_Neighbours() : p_rest_from(-INFINITY), p_rest_to(-INFINITY) {}
    std::vector< std::pair< unsigned, Float_Type > > from_v;
    std::vector< std::pair< unsigned, Float_Type > > to_v;
    Float_Type p_rest_from;
    Float_Type p_rest_to;
};

template < typename Float_Type = float, unsigned Kmer_Size = 6 >
class State_Transitions
{
public:
    typedef Kmer< Kmer_Size > Kmer_Type;
    typedef State_Neighbours< Float_Type > State_Neighbours_Type;
    static const unsigned n_states = 1u << (2 * Kmer_Size);

    State_Transitions() { _neighbours.resize(n_states); }
    void clear() { _neighbours.clear(); _neighbours.resize(n_states); }

    const State_Neighbours_Type& neighbours(unsigned i) const { return _neighbours.at(i); }
    State_Neighbours_Type& neighbours(unsigned i) { return _neighbours.at(i); }

    // update fields from_v, p_rest_from, p_rest_to based on to_v
    void update_fields()
    {
        for (unsigned i = 0; i < n_states; ++i)
        {
            neighbours(i).from_v.clear();
        }
        for (unsigned i = 0; i < n_states; ++i)
        {
            logsum::logsumset< Float_Type > s(false);
            for (const auto& p : neighbours(i).to_v)
            {
                neighbours(p.first).from_v.push_back(std::make_pair(i, p.second));
                s.add(p.second);
            }
            neighbours(i).p_rest_to = std::log(1 - std::exp(s.val()));
        }
        for (unsigned i = 0; i < n_states; ++i)
        {
            logsum::logsumset< Float_Type > s(false);
            for (const auto& p : neighbours(i).from_v)
            {
                s.add(p.second);
            }
            neighbours(i).p_rest_from = std::log(1 - std::exp(s.val()));
        }
    }

    // drop transitions with low probability
    void drop_transitions(Float_Type p_cutoff)
    {
        Float_Type log_p_cutoff = std::log(p_cutoff);
        for (unsigned i = 0; i < n_states; ++i)
        {
            decltype(neighbours(i).to_v) to_v;
            for (const auto& p : neighbours(i).to_v)
            {
                if (p.second > log_p_cutoff)
                {
                    to_v.push_back(p);
                }
            }
            neighbours(i).to_v = std::move(to_v);
        }
        update_fields();
    }

    // recompute transition table
    void compute_transitions(Float_Type p_skip_default, Float_Type p_stay, Float_Type p_cutoff,
                             const std::map< unsigned, Float_Type >& p_skip_map = {})
    {
        clear();
        for (size_t i = 0; i < n_states; ++i)
        {
            Float_Type p_skip = p_skip_default;
            if (p_skip_map.count(i))
            {
                p_skip = p_skip_map.at(i);
            }
            for (size_t j = 0; j < n_states; ++j)
            {
                Float_Type p = 0;
                if (i == j)
                {
                    p += p_stay;
                }
                for (unsigned k = 1; k < Kmer_Size; ++k)
                    if ((i & ((1u << (2 * k)) - 1)) == (j >> (2 * (Kmer_Size - k))))
                    {
                        p += (1 - p_stay) * (1 - p_skip) * pow(p_skip, Kmer_Size - k - 1)
                            / (1u << (2 * (Kmer_Size - k)));
                    }
                p += (1 - p_stay) * pow(p_skip, 5) / n_states;
                if (p > p_cutoff)
                {
                    neighbours(i).to_v.push_back(std::make_pair(j, std::log(p)));
                }
            }
        }
        update_fields();
    }

    friend std::ostream& operator << (std::ostream& os, const State_Transitions& st)
    {
        for (unsigned i = 0; i < n_states; ++i)
        {
            for (const auto& p : st.neighbours(i).to_v)
            {
                os << Kmer_Type::to_string(i) << '\t' << Kmer_Type::to_string(p.first) << '\t' << p.second << std::endl;
            }
        }
        return os;
    }
    friend std::istream& operator >> (std::istream& is, State_Transitions& st)
    {
        st.clear();
        std::string k_i;
        std::string k_j;
        Float_Type p;
        while (is >> k_i >> k_j >> p)
        {
            unsigned i = Kmer_Type::to_int(k_i);
            unsigned j = Kmer_Type::to_int(k_j);
            st.neighbours(i).to_v.push_back(std::make_pair(j, p));
        }
        st.update_fields();
        return is;
    }

private:
    std::vector< State_Neighbours_Type > _neighbours;
}; // class State_Transitions

typedef State_Transitions<> State_Transitions_Type;

#endif
