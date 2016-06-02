#ifndef __FORWARD_BACKWARD_HPP
#define __FORWARD_BACKWARD_HPP

#include <cmath>
#include <iostream>
#include <vector>
#include <set>

#include "Pore_Model.hpp"
#include "State_Transitions.hpp"
#include "logsumset.hpp"
#include "logger.hpp"

template < typename Float_Type, unsigned Kmer_Size = 6 >
class Forward_Backward
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
        Float_Type alpha; // := Pr[ E_1 ... E_i, S_i = j ]
        Float_Type beta;  // := Pr[ E_{i+1} ... E_n | S_i = j ]
    }; // struct Matrix_Entry

    static const unsigned n_states = Pore_Model_Type::n_states;

    void clear() { _m.clear(); }
    unsigned n_events() const { return _m.size() / n_states; }

    // i: event index
    // j: state/kmer index
    const Matrix_Entry& cell(unsigned i, unsigned j) const { return _m[i * n_states + j]; }
    Matrix_Entry& cell(unsigned i, unsigned j) { return _m[i * n_states + j]; }

    Float_Type log_posterior(unsigned i, unsigned j) const { return cell(i, j).alpha + cell(i, j).beta - _log_pr_data; }
    Float_Type log_pr_data() const { return _log_pr_data; }

    static unsigned& n_threads() { static unsigned _n_threads = 1; return _n_threads; }

    void fill(const Pore_Model_Type& pm,
              const State_Transitions_Type& st,
              const Event_Sequence_Type& ev)
    {
        clear();
        unsigned n_events = ev.size();
        _m.resize(n_states * n_events);
        Float_Type log_n_states = std::log(static_cast< Float_Type >(n_states));
        LogSumSet_Type s(false);
        //
        // forward: alpha, i == 0
        //
        {
            unsigned i = 0;
            LOG("Forward_Backward", debug1) << "forward: i=" << i << std::endl;
            for (unsigned j = 0; j < n_states; ++j)
            {
                cell(i, j).alpha = pm.log_pr_emission(j, ev[0]) - log_n_states;
                LOG("Forward_Backward", debug2)
                    << "i=" << i << " j=" << j << " kmer_j=" << Kmer_Type::to_string(j)
                    << " alpha=" << cell(i, j).alpha << std::endl;
            }
        }
        //
        // forward: alpha, i > 0
        //
        for (unsigned i = 1; i < ev.size(); ++i)
        {
            LOG("Forward_Backward", debug1) << "forward: i=" << i << std::endl;
            for (unsigned j = 0; j < n_states; ++j)
            {
                s.clear();
                for (const auto& p : st.neighbours(j).from_v)
                {
                    const unsigned& j_prev = p.first;
                    const Float_Type& log_pr_transition = p.second;
                    s.add(log_pr_transition + cell(i - 1, j_prev).alpha);
                }
                cell(i, j).alpha = pm.log_pr_emission(j, ev[i]) + s.val();
                LOG("Forward_Backward", debug2)
                    << "i=" << i << " j=" << j << " kmer_j=" << Kmer_Type::to_string(j)
                    << " alpha=" << cell(i, j).alpha << std::endl;
            }
        }
        //
        // backward: beta, i == n-1
        //
        {
            unsigned i = ev.size() - 1;
            LOG("Forward_Backward", debug1) << "backward: i=" << i << std::endl;
            for (unsigned j = 0; j < n_states; ++j)
            {
                cell(i, j).beta = 0;
                LOG("Forward_Backward", debug2)
                    << "i=" << i << " j=" << j << " kmer_j=" << Kmer_Type::to_string(j)
                    << " beta=" << cell(i, j).beta << std::endl;
            }
        }
        //
        // backward: beta, i < n-1
        //
        for (unsigned ip1 = ev.size() - 1; ip1 > 0; --ip1)
        {
            unsigned i = ip1 - 1;
            LOG("Forward_Backward", debug1) << "backward: i=" << i << std::endl;
            for (unsigned j = 0; j < n_states; ++j)
            {
                s.clear();
                for (const auto& p : st.neighbours(j).to_v)
                {
                    const unsigned& j_next = p.first;
                    const Float_Type& log_pr_transition = p.second;
                    s.add(log_pr_transition + pm.log_pr_emission(j_next, ev[ip1]) + cell(ip1, j_next).beta);
                }
                cell(i, j).beta += s.val();
                LOG("Forward_Backward", debug2)
                    << "i=" << i << " j=" << j << " kmer_j=" << Kmer_Type::to_string(j)
                    << " beta=" << cell(i, j).beta << std::endl;
            }
        }
        //
        // pr_data
        //
        s.clear();
        for (unsigned j = 0; j < n_states; ++j)
        {
            s.add(cell(ev.size() - 1, j).alpha);
        }
        _log_pr_data = s.val();
    }

    friend std::ostream& operator << (std::ostream& os, const Forward_Backward& fwbw)
    {
        for (unsigned i = 0; i < fwbw.n_events(); ++i)
        {
            for (unsigned j = 0; j < fwbw.n_states; ++j)
            {
                os << i << '\t' << j << '\t'
                   << fwbw.cell(i, j).alpha << '\t'
                   << fwbw.cell(i, j).beta << std::endl;
            }
        }
        return os;
    }

private:
    std::vector< Matrix_Entry > _m;
    Float_Type _log_pr_data;
}; // class Forward_Backward

#endif
