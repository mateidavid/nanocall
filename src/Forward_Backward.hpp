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

template < typename Float_Type = float, unsigned Kmer_Size = 6 >
class Forward_Backward
{
public:
    typedef Kmer< Kmer_Size > Kmer_Type;
    typedef Pore_Model< Float_Type, Kmer_Size > Pore_Model_Type;
    typedef State_Transitions< Float_Type, Kmer_Size > State_Transitions_Type;
    typedef Event< Float_Type > Event_Type;
    typedef logsumset< Float_Type > LogSumSet_Type;

    struct Matrix_Entry
    {
        Float_Type alpha; // := Pr[ S_i = j | e_1 ... e_{i-1} ]
        Float_Type beta;  // := Pr[ S_i = j | e_1 ... e_i ]
        Float_Type gamma; // := Pr[ S_i = j | e_1 ... e_n ]
    }; // struct Matrix_Entry

    static const unsigned n_states = Pore_Model_Type::n_states;

    void clear() { _m.clear(); }
    unsigned n_events() const { return _m.size() / n_states; }

    // i: event index
    // j: state/kmer index
    const Matrix_Entry& cell(unsigned i, unsigned j) const { return _m[i * n_states + j]; }
    Matrix_Entry& cell(unsigned i, unsigned j) { return _m[i * n_states + j]; }
    Float_Type posterior(unsigned i, unsigned j) { return cell(i, j).gamma; }

    static unsigned& n_threads() { static unsigned _n_threads = 1; return _n_threads; }

    void fill(const Pore_Model_Type& pm,
              const State_Transitions_Type& st,
              const std::vector< Event_Type >& ev)
    {
        clear();
        unsigned n_events = ev.size();
        _m.resize(n_states * n_events);
        float log_n_states = std::log(static_cast< float >(n_states));
        LogSumSet_Type s1(false);
        LogSumSet_Type s2(false);
        //
        // forward: alpha, beta; i == 0
        //
        {
            for (unsigned j = 0; j < n_states; ++j)
            {
                // alpha
                cell(0, j).alpha = - log_n_states;
                // beta
                cell(0, j).beta = pm.log_pr_emission(j, ev[0]) + cell(0, j).alpha;
                s1.add(cell(0, j).beta);
            }
            Float_Type denom = s1.val();
            LOG("Forward_Backward", debug) << "i=0 beta_denom=" << denom << std::endl;
            for (unsigned j = 0; j < n_states; ++j)
            {
                cell(0, j).beta -= denom;
                LOG("Forward_Backward", debug)
                    << "i=0 j=" << Kmer_Type::to_string(j)
                    << " alpha=" << cell(0, j).alpha
                    << " beta=" << cell(0, j).beta << std::endl;
            }
        }
        //
        // forward: alpha, beta; i > 0
        //
        for (unsigned i = 1; i < ev.size(); ++i)
        {
            LOG("Forward_Backward", info) << "forward: i=" << i << std::endl;
            s1.clear();
            for (unsigned j = 0; j < n_states; ++j) // TODO: parallelize
            {
                // alpha
                s2.clear();
                for (const auto& p : st.neighbours(j).from_v)
                {
                    const unsigned& j_prev = p.first;
                    const Float_Type& log_pr_transition = p.second;
                    s2.add(log_pr_transition + cell(i - 1, j_prev).beta);
                }
                cell(i, j).alpha = s2.val();
                // beta
                cell(i, j).beta = pm.log_pr_emission(j, ev[i]) + cell(i, j).alpha;
                s1.add(cell(i, j).beta);
            }
            Float_Type denom = s1.val();
            LOG("Forward_Backward", debug) << "i=" << i << " beta_denom=" << denom << std::endl;
            for (unsigned j = 0; j < n_states; ++j)
            {
                cell(i, j).beta -= denom;
                LOG("Forward_Backward", debug)
                    << "i=" << i << " j=" << Kmer_Type::to_string(j)
                    << " alpha=" << cell(i, j).alpha
                    << " beta=" << cell(i, j).beta << std::endl;
            }
        }
        //
        // backward, gamma; i == n-1
        //
        for (unsigned j = 0; j < n_states; ++j)
        {
            cell(n_events - 1, j).gamma = cell(n_events - 1, j).beta;
        }
        //
        // backward, gamma; i < n-1
        //
        for (unsigned ip1 = ev.size() - 1; ip1 > 0; --ip1)
        {
            unsigned i = ip1 - 1;
            LOG("Forward_Backward", info) << "backward: i=" << i << std::endl;
            for (unsigned j = 0; j < n_states; ++j) // TODO: parallelize
            {
                cell(i, j).gamma = cell(i, j).beta;
                s2.clear();
                for (const auto& p : st.neighbours(j).to_v)
                {
                    const unsigned& j_next = p.first;
                    const Float_Type& log_pr_transition = p.second;
                    s2.add(log_pr_transition + cell(ip1, j_next).gamma - cell(ip1, j_next).alpha);
                }
                cell(i, j).gamma += s2.val();
                LOG("Forward_Backward", debug)
                    << "i=" << i << " j=" << Kmer_Type::to_string(j)
                    << " gamma=" << cell(i, j).gamma << std::endl;
            }
        }
    }

    friend std::ostream& operator << (std::ostream& os, const Forward_Backward& fwbw)
    {
        for (unsigned i = 0; i < fwbw.n_events(); ++i)
        {
            for (unsigned j = 0; j < fwbw.n_states; ++j)
            {
                os << i << '\t' << j << '\t'
                   << fwbw.cell(i, j).alpha << '\t'
                   << fwbw.cell(i, j).beta << '\t'
                   << fwbw.cell(i, j).gamma << std::endl;
            }
        }
        return os;
    }

private:
    std::vector< Matrix_Entry > _m;
}; // class Forward_Backward

typedef Forward_Backward<> Forward_Backward_Type;

#endif
