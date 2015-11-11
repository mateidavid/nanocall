#ifndef __VITERBI_HPP
#define __VITERBI_HPP

#include <cmath>
#include <iostream>
#include <vector>
#include <set>

#include "Pore_Model.hpp"
#include "State_Transitions.hpp"
#include "logsumset.hpp"
#include "logger.hpp"

template < typename Float_Type = float, unsigned Kmer_Size = 6 >
class Viterbi
{
public:
    typedef Kmer< Kmer_Size > Kmer_Type;
    typedef Pore_Model< Float_Type, Kmer_Size > Pore_Model_Type;
    typedef State_Transitions< Float_Type, Kmer_Size > State_Transitions_Type;
    typedef Event< Float_Type > Event_Type;
    typedef logsumset< Float_Type > LogSumSet_Type;

    struct Matrix_Entry
    {
        Float_Type alpha; // := Pr[ MLSS producing e_1 ... e_i, with S_i == j ]
        unsigned beta;    // := previous state in the MLSS
    }; // struct Matrix_Entry

    static const unsigned n_states = Pore_Model_Type::n_states;

    void clear() { _m.clear(); _state_seq.clear(); _base_seq.clear(); }
    unsigned n_events() const { return _state_seq.size(); }
    const std::vector< unsigned >& state_seq() const { return _state_seq; }
    const std::string& base_seq() const { return _base_seq; }

    // i: event index
    // j: state/kmer index
    const Matrix_Entry& cell(unsigned i, unsigned j) const { return _m[i * n_states + j]; }
    Matrix_Entry& cell(unsigned i, unsigned j) { return _m[i * n_states + j]; }

    static unsigned& n_threads() { static unsigned _n_threads = 1; return _n_threads; }

    void fill(const Pore_Model_Type& pm,
              const State_Transitions_Type& st,
              const std::vector< Event_Type >& ev)
    {
        clear();
        unsigned n_events = ev.size();
        _m.resize(n_states * n_events);
        _state_seq.resize(n_events);
        float log_n_states = std::log(static_cast< float >(n_states));
        //
        // alpha, beta; i == 0
        //
        {
            for (unsigned j = 0; j < n_states; ++j)
            {
                // alpha
                cell(0, j).alpha = pm.log_pr_emission(j, ev[0]) - log_n_states;
                // beta
                cell(0, j).beta = n_states;
                LOG("Viterbi", debug)
                    << "i=0 j=" << Kmer_Type::to_string(j)
                    << " alpha=" << cell(0, j).alpha
                    << " beta=" << cell(0, j).beta << std::endl;
            }
        }
        //
        // alpha, beta; i > 0
        //
        for (unsigned i = 1; i < n_events; ++i)
        {
            LOG("Viterbi", info) << "forward: i=" << i << std::endl;
            for (unsigned j = 0; j < n_states; ++j) // TODO: parallelize
            {
                cell(i, j).alpha = -INFINITY;
                cell(i, j).beta = n_states;
                for (const auto& p : st.neighbours(j).from_v)
                {
                    const unsigned& j_prev = p.first;
                    const Float_Type& log_pr_transition = p.second;
                    Float_Type v = log_pr_transition + cell(i - 1, j_prev).alpha;
                    if (v > cell(i, j).alpha)
                    {
                        cell(i, j).alpha = v;
                        cell(i, j).beta = j_prev;
                    }
                }
                cell(i, j).alpha += pm.log_pr_emission(j, ev[i]);
                LOG("Viterbi", debug)
                    << "i=" << i << " j=" << Kmer_Type::to_string(j)
                    << " alpha=" << cell(i, j).alpha
                    << " beta=" << cell(i, j).beta << std::endl;
            }
        }
        fill_state_seq();
        fill_base_seq();
    }

    friend std::ostream& operator << (std::ostream& os, const Viterbi& vit)
    {
        for (unsigned i = 0; i < vit.n_events(); ++i)
        {
            for (unsigned j = 0; j < vit.n_states; ++j)
            {
                os << i << '\t' << j << '\t'
                   << vit.cell(i, j).alpha << '\t'
                   << vit.cell(i, j).beta << std::endl;
            }
        }
        return os;
    }

private:
    std::vector< Matrix_Entry > _m;
    std::vector< unsigned > _state_seq;
    std::string _base_seq;

    void fill_state_seq()
    {
        Float_Type max_v = -INFINITY;
        unsigned max_j = n_states;
        for (unsigned j = 0; j < n_states; ++j)
        {
            if (cell(n_events() - 1, j).alpha > max_v)
            {
                max_j = j;
                max_v = cell(n_events() - 1, j).alpha;
            }
        }
        for (unsigned i = n_events() - 1; i > 0; --i)
        {
            _state_seq.at(i) = max_j;
            max_j = cell(i, max_j).beta;
        }
        _state_seq.at(0) = max_j;
    }

    void fill_base_seq()
    {
        for (unsigned i = 0; i < _state_seq.size() - 1; ++i)
        {
            auto c = Kmer_Type::min_skip(_state_seq[i], _state_seq[i + 1]);
            LOG("Viterbi", debug)
                << "i=" << i << " state=" << _state_seq[i] << " kmer=" << Kmer_Type::to_string(_state_seq[i]) << " c=" << c << std::endl;
            _base_seq += Kmer_Type::to_string(_state_seq[i]).substr(0, c);
        }
        LOG("Viterbi", debug)
            << "i=" << n_events() - 1 << " state=" << _state_seq[n_events() - 1]
            << " kmer=" << Kmer_Type::to_string(_state_seq[n_events() - 1]) << std::endl;
        _base_seq += Kmer_Type::to_string(_state_seq[n_events() - 1]);
    }

}; // class Viterbi

#endif
