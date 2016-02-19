#ifndef __STATE_TRANSITIONS_BASE_HPP
#define __STATE_TRANSITIONS_BASE_HPP

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>

#include "Kmer.hpp"
#include "logsumset.hpp"
#include "logger.hpp"

template < typename Float_Type = float, unsigned Kmer_Size = 6 >
struct State_Transition_Parameters
{
    Float_Type p_stay;
    Float_Type p_skip;

    static Float_Type& default_p_stay()
    {
        static Float_Type _default_p_stay = .09;
        return _default_p_stay;
    }
    static Float_Type& default_p_skip()
    {
        static Float_Type _default_p_skip = .28;
        return _default_p_skip;
    }

    State_Transition_Parameters()
        : p_stay(default_p_stay()), p_skip(default_p_skip()) {}

    bool is_default() const
    {
        return p_stay == default_p_stay() and p_skip == default_p_skip();
    }

    friend std::ostream& operator << (std::ostream& os, const State_Transition_Parameters& stp)
    {
        os << "[p_stay=" << stp.p_stay
           << " p_skip=" << stp.p_skip << "]";
        return os;
    }
    void write_tsv(std::ostream& os) const
    {
        os << std::fixed << std::setprecision(5)
           << p_stay << '\t'
           << p_skip;
    }
}; // struct State_Transition_Parameters

template < typename Float_Type = float >
struct State_Neighbours
{
    State_Neighbours() : p_rest_from(-INFINITY), p_rest_to(-INFINITY) {}
    std::vector< std::pair< unsigned, Float_Type > > from_v;
    std::vector< std::pair< unsigned, Float_Type > > to_v;
    Float_Type p_rest_from;
    Float_Type p_rest_to;
}; // struct State_Neighbours

template < typename Float_Type = float, unsigned Kmer_Size = 6 >
class State_Transitions
{
public:
    typedef Kmer< Kmer_Size > Kmer_Type;
    typedef State_Neighbours< Float_Type > State_Neighbours_Type;
    typedef State_Transition_Parameters< Float_Type > State_Transition_Parameters_Type;
    static const unsigned n_states = 1u << (2 * Kmer_Size);

    State_Transitions() = default;

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

    static Float_Type get_trans_prob(unsigned i, unsigned j,
                                     Float_Type p_stay, Float_Type p_step, Float_Type p_skip_1)
    {
        Float_Type p = 0;
        if (i == j)
        {
            p += p_stay;
        }
        if (Kmer_Type::suffix(i, Kmer_Size - 1) == Kmer_Type::prefix(j, Kmer_Size - 1))
        {
            p += p_step / 4;
        }
        for (unsigned l = 2; l < Kmer_Size; ++l)
            if (Kmer_Type::suffix(i, Kmer_Size - l) == Kmer_Type::prefix(j, Kmer_Size - l))
            {
                p += pow(p_skip_1, l - 1) / (1u << (2 * l));
            }
        p += (pow(p_skip_1, 5) / (Float_Type(1.0) - p_skip_1)) / n_states;
        return p;
    }

    // recompute transition table
    void compute_transitions(Float_Type p_skip_default, Float_Type p_stay, Float_Type p_cutoff,
                             const std::map< unsigned, Float_Type >& p_skip_map = {})
    {
        _neighbours.clear();
        _neighbours.reserve(n_states);
        for (unsigned i = 0; i < n_states; ++i)
        {
            _neighbours.emplace_back();
            Float_Type p_skip = p_skip_default;
            if (p_skip_map.count(i))
            {
                p_skip = p_skip_map.at(i);
            }
            Float_Type p_step = 1.0 - p_stay - p_skip;
            // p_skip = sum_{i>=1} p_skip_1^i
            Float_Type p_skip_1 = p_skip / (p_skip + 1.0);
            LOG(debug2) << "i=" << Kmer_Type::to_string(i)
                        << " p_stay=" << p_stay
                        << " p_skip=" << p_skip
                        << " p_step=" << p_step
                        << " p_skip_1=" << p_skip_1 << std::endl;
            for (unsigned j = 0; j < n_states; ++j)
            {
                Float_Type p = get_trans_prob(i, j, p_stay, p_step, p_skip_1);
                if (p > p_cutoff)
                {
                    neighbours(i).to_v.push_back(std::make_pair(j, std::log(p)));
                }
            }
        }
        update_fields();
    }

    // compute transition table allowing a maximum of 1 skip
    void compute_transitions_fast(Float_Type p_skip_default, Float_Type p_stay,
                                  const std::map< unsigned, Float_Type >& p_skip_map = {})
    {
        struct Default_Float_Type
        {
            Default_Float_Type(Float_Type _val = 0.0) : val(_val) {}
            Float_Type val;
        }; // struct Default_Float

        _neighbours.clear();
        _neighbours.reserve(n_states);
        for (unsigned i = 0; i < n_states; ++i)
        {
            _neighbours.emplace_back();
            Float_Type p_skip = p_skip_default;
            if (p_skip_map.count(i))
            {
                p_skip = p_skip_map.at(i);
            }
            Float_Type p_step = 1.0 - p_stay - p_skip;
            // p_skip = sum_{i>=1} p_skip_1^i
            Float_Type p_skip_1 = p_skip / (p_skip + 1.0);
            LOG(debug2) << "i=" << Kmer_Type::to_string(i)
                        << " p_stay=" << p_stay
                        << " p_skip=" << p_skip
                        << " p_step=" << p_step
                        << " p_skip_1=" << p_skip_1 << std::endl;
            std::set< unsigned > to_s{i};
            const auto& nl1 = Kmer_Type::neighbour_list(i, 1);
            to_s.insert(nl1.begin(), nl1.end());
            const auto& nl2 = Kmer_Type::neighbour_list(i, 2);
            to_s.insert(nl2.begin(), nl2.end());
            for (const auto& j : to_s)
            {
                Float_Type p = get_trans_prob(i, j, p_stay, p_step, p_skip_1);
                neighbours(i).to_v.push_back(std::make_pair(j, std::log(p)));
            }
        }
        update_fields();
    }
    void compute_transitions_fast(const State_Transition_Parameters_Type& stp)
    {
        compute_transitions_fast(stp.p_skip, stp.p_stay);
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
        st._neighbours.clear();
        st._neighbours.resize(n_states);
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

typedef State_Transition_Parameters<> State_Transition_Parameters_Type;
typedef State_Transitions<> State_Transitions_Type;

#endif
