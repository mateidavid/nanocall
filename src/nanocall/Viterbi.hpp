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
#include "fast5.hpp"

template < typename Float_Type, unsigned Kmer_Size = 6 >
class Viterbi
{
public:
    typedef Kmer< Kmer_Size > Kmer_Type;
    typedef Pore_Model< Float_Type, Kmer_Size > Pore_Model_Type;
    typedef State_Transitions< Float_Type, Kmer_Size > State_Transitions_Type;
    typedef Event< Float_Type, Kmer_Size > Event_Type;
    typedef Event_Sequence< Float_Type, Kmer_Size > Event_Sequence_Type;
    typedef logsum::logsumset< Float_Type > LogSumSet_Type;

    struct Matrix_Entry
    {
        Float_Type alpha; // := Pr[ MLSS producing e_1 ... e_i, with S_i == j ]
        unsigned beta;    // := previous state in the MLSS
    }; // struct Matrix_Entry

    static const unsigned n_states = Pore_Model_Type::n_states;

    unsigned n_events() const { return _n_events; }
    Float_Type path_probability() const { return _path_probability; }

    // i: event index
    // j: state/kmer index
    const Matrix_Entry& cell(unsigned i, unsigned j) const { return _m[i * n_states + j]; }
    Matrix_Entry& cell(unsigned i, unsigned j) { return _m[i * n_states + j]; }

    static unsigned& n_threads() { static unsigned _n_threads = 1; return _n_threads; }

    void fill(const Pore_Model_Type& pm,
              const State_Transitions_Type& st,
              Event_Sequence_Type& ev)
    {
        _n_events = ev.size();
        _m.clear();
        _m.resize(n_states * n_events());
        Float_Type log_n_states = std::log(static_cast< Float_Type >(n_states));
        //
        // alpha, beta; i == 0
        //
        {
            LOG("Viterbi", debug1) << "forward: i=0" << std::endl;
            for (unsigned j = 0; j < n_states; ++j)
            {
                // alpha
                cell(0, j).alpha = pm.log_pr_corrected_emission(j, ev[0]) - log_n_states;
                // beta
                cell(0, j).beta = n_states;
                LOG("Viterbi", debug2)
                    << "i=0 j=" << Kmer_Type::to_string(j)
                    << " alpha=" << cell(0, j).alpha
                    << " beta=" << cell(0, j).beta << std::endl;
            }
        }
        //
        // alpha, beta; i > 0
        //
        for (unsigned i = 1; i < n_events(); ++i)
        {
            LOG("Viterbi", debug1) << "forward: i=" << i << std::endl;
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
                cell(i, j).alpha += pm.log_pr_corrected_emission(j, ev[i]);
                LOG("Viterbi", debug2)
                    << "i=" << i << " j=" << Kmer_Type::to_string(j)
                    << " alpha=" << cell(i, j).alpha
                    << " beta=" << cell(i, j).beta << std::endl;
            }
        }
        fill_state_seq(ev);
        fill_move_seq(ev);
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
    Float_Type _path_probability;
    unsigned _n_events;

    void fill_state_seq(Event_Sequence_Type& ev)
    {
        assert(Kmer_Size <= MAX_K_LEN);
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
        _path_probability = max_v;
        for (unsigned i = n_events() - 1; i > 0; --i)
        {
            ev[i].model_state_idx = max_j;
            ev[i].set_model_state(Kmer_Type::to_string(ev[i].model_state_idx));
            max_j = cell(i, max_j).beta;
        }
        ev[0].model_state_idx = max_j;
        ev[0].set_model_state(Kmer_Type::to_string(ev[0].model_state_idx));
    }

    void fill_move_seq(Event_Sequence_Type& ev)
    {
        for (unsigned i = 0; i < n_events(); ++i)
        {
            ev[i].move = i > 0? Kmer_Type::min_skip(ev[i - 1].model_state_idx, ev[i].model_state_idx) : 0u;
        }
    }

}; // class Viterbi

#endif
