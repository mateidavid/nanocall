#ifndef __PARAMETER_TRAINER
#define __PARAMETER_TRAINER

#include <array>
#include <vector>
#include <map>

#include "global_assert.hpp"
#include "Pore_Model.hpp"
#include "State_Transitions.hpp"
#include "Forward_Backward.hpp"
#include "logsumset.hpp"
#include "logger.hpp"

#if defined(USE_LOGDIFF)
#include "logdiff.hpp"
#endif

template < typename Float_Type = float, unsigned Kmer_Size = 6 >
struct Parameter_Trainer
{
    typedef Kmer< Kmer_Size > Kmer_Type;
    typedef Pore_Model< Float_Type, Kmer_Size > Pore_Model_Type;
    typedef State_Transitions< Float_Type, Kmer_Size > State_Transitions_Type;
    typedef Event< Float_Type > Event_Type;
    typedef Event_Sequence< Float_Type > Event_Sequence_Type;
    typedef Forward_Backward< Float_Type, Kmer_Size > Forward_Backward_Type;
    typedef logsum::logsumset< Float_Type > LogSumSet_Type;
    static const unsigned n_states = Pore_Model_Type::n_states;

    static void init()
    {
        // pick states i s.t. i has self-overlap 0,
        // and all its 1-step neighbours have self-overlap <=1
        st_train_kmers().clear();
        for (unsigned i = 0; i < n_states; ++i)
        {
            if (Kmer_Type::max_self_overlap(i) > 0)
            {
                continue;
            }
            bool all_good = true;
            for (unsigned b1 = 0; b1 < 4; ++b1)
            {
                unsigned j = (Kmer_Type::suffix(i, Kmer_Size - 1) << 2) + b1;
                if (Kmer_Type::max_self_overlap(j) > 1)
                {
                    all_good = false;
                    break;
                }
            }
            if (all_good)
            {
                st_train_kmers().push_back(i);
            }
        }
        LOG(info) << "using [" << st_train_kmers().size() << "] states for state trainsition training" << std::endl;
    }

    static std::vector< unsigned >& st_train_kmers()
    {
        static std::vector< unsigned > _st_train_kmers;
        return _st_train_kmers;
    }

    /**
     * Perform one training round.
     * @event_seq_ptrs Vector of pairs, first: an event sequence, second: strand from which it comes
     * @model_ptrs Pointers to unscaled pore models (per strand)
     * @default_transitions Default state transitions
     * @crt_pm_params Pore model scaling parameters (common to both strands)
     * @crt_st_params State transition parameters (per strand)
     * @new_pm_params Destination for trained pm params (common to both strands)
     * @new_st_params Destination for trained st params (per strand)
     * @fit Destination for pr_data using crt params
     * @done Bool; set to true if no more training rounds can be performed due to singularity.
     */
    static void train_one_round(
        const std::vector< std::pair< const Event_Sequence_Type*, unsigned > >& event_seq_ptrs,
        const std::array< const Pore_Model_Type*, 2 >& model_ptrs,
        const State_Transitions_Type& default_transitions,
        const Pore_Model_Parameters_Type& crt_pm_params,
        const std::array< State_Transition_Parameters_Type, 2 >& crt_st_params,
        Pore_Model_Parameters_Type& new_pm_params,
        std::array< State_Transition_Parameters_Type, 2 >& new_st_params,
        Float_Type& fit,
        bool& done)
    {
        // tags of event sequences must be 0/1
        ASSERT(alg::all_of(
                   event_seq_ptrs,
                   [] (const std::pair< const Event_Sequence_Type*, unsigned >& p) {
                       return p.second < 2;
                   }));
        unsigned n_event_seqs = event_seq_ptrs.size();
        unsigned total_num_events = alg::accumulate(
            event_seq_ptrs, 0u,
            [] (unsigned s, const std::pair< const Event_Sequence_Type*, unsigned >& p) {
                return s + p.first->size();
            });
        done = false;
        // compute scaled pore models
        std::array< Pore_Model_Type, 2 > scaled_models;
        std::array< bool, 2 > init_scaled_models = {{ false, false }};
        for (const auto& p : event_seq_ptrs)
        {
            if (init_scaled_models[p.second]) continue;
            ASSERT(model_ptrs[p.second]);
            scaled_models[p.second] = *model_ptrs[p.second];
            scaled_models[p.second].scale(crt_pm_params);
            init_scaled_models[p.second] = true;
        }
        // compute custom state transitions
        std::array< State_Transitions_Type, 2 > custom_transitions;
        std::array< const State_Transitions_Type*, 2 > transition_ptrs;
        std::array< bool, 2 > init_transitions = {{ false, false }};
        for (const auto& p : event_seq_ptrs)
        {
            if (init_transitions[p.second]) continue;
            if (not crt_st_params[p.second].is_default())
            {
                custom_transitions[p.second].compute_transitions_fast(crt_st_params[p.second]);
                transition_ptrs[p.second] = &custom_transitions[p.second];
            }
            else
            {
                transition_ptrs[p.second] = &default_transitions;
            }
            init_transitions[p.second] = true;
        }
        // compute drift-corrected event sequences, and run fwbw
        std::vector< Event_Sequence_Type > corrected_event_seqs(n_event_seqs);
        std::vector< Forward_Backward_Type > fwbw(n_event_seqs);
        fit = 0;
        for (unsigned k = 0; k < n_event_seqs; ++k)
        {
            unsigned st = event_seq_ptrs[k].second;
            ASSERT(init_scaled_models[st]);
            ASSERT(init_transitions[st]);
            Event_Sequence_Type& corrected_events = corrected_event_seqs[k];
            corrected_events = *event_seq_ptrs[k].first;
            corrected_events.apply_drift_correction(crt_pm_params.drift);
            fwbw[k].fill(scaled_models[st], *transition_ptrs[st], corrected_events);
            fit += fwbw[k].log_pr_data();
        }
        // compute the scaling matrices (first in logspace)
        // against unscaled pm & uncorrected events
        std::array< std::array< LogSumSet_Type, 3 >, 3 > A_lss =
            {{ {{ false, false, false }},
               {{ false, false, false }},
               {{ false, false, false }} }};
        std::array< LogSumSet_Type, 3 > B_lss = {{ false, false, false }};
        for (unsigned k = 0; k < n_event_seqs; ++k)
        {
            unsigned st = event_seq_ptrs[k].second;
            const Event_Sequence_Type& events = *event_seq_ptrs[k].first;
            const Pore_Model_Type& pm = *model_ptrs[st];
            unsigned n_events = events.size();
            for (unsigned i = 0; i < n_events; ++i)
            {
                Float_Type log_x_i = events[i].log_mean;
                Float_Type log_t_i = events[i].log_start;
                LOG(debug1)
                    << "outter_loop k=" << k << " i=" << i
                    << " log_x_i=" << log_x_i
                    << " log_t_i=" << log_t_i << std::endl;
                std::array< LogSumSet_Type, 3 > s_lss = {{ false, false, false }};
                for (unsigned j = 0; j < Pore_Model_Type::n_states; ++j)
                {
                    Float_Type x0 = fwbw[k].log_posterior(i, j) - 2 * pm.state(j).log_level_stdv;
                    Float_Type x1 = x0 + pm.state(j).log_level_mean;
                    Float_Type x2 = x1 + pm.state(j).log_level_mean;
                    LOG(debug2)
                        << "inner_loop k=" << k << " i=" << i << " j=" << j
                        << " x0=" << x0 << " x1=" << x1 << " x2=" << x2 << std::endl;
                    s_lss[0].add(x0);
                    s_lss[1].add(x1);
                    s_lss[2].add(x2);
                }
                A_lss[0][0].add(s_lss[0].val());
                A_lss[0][1].add(s_lss[1].val());
                A_lss[0][2].add(s_lss[0].val() + log_t_i);
                A_lss[1][1].add(s_lss[2].val());
                A_lss[1][2].add(s_lss[1].val() + log_t_i);
                A_lss[2][2].add(s_lss[0].val() + 2 * log_t_i);
                B_lss[0].add(s_lss[0].val() + log_x_i);
                B_lss[1].add(s_lss[1].val() + log_x_i);
                B_lss[2].add(s_lss[0].val() + log_x_i + log_t_i);
            }
        }
        // now compute matrices in normal space
        std::array< std::array< Float_Type, 3 >, 3 > A =
            {{ {{ std::exp(A_lss[0][0].val()), std::exp(A_lss[0][1].val()), std::exp(A_lss[0][2].val()) }},
               {{ std::exp(A_lss[0][1].val()), std::exp(A_lss[1][1].val()), std::exp(A_lss[1][2].val()) }},
               {{ std::exp(A_lss[0][2].val()), std::exp(A_lss[1][2].val()), std::exp(A_lss[2][2].val()) }} }};
        std::array< Float_Type, 3 > B =
            {{ std::exp(B_lss[0].val()), std::exp(B_lss[1].val()), std::exp(B_lss[2].val()) }};
#ifndef NDEBUG
        auto A_copy = A;
        auto B_copy = B;
#endif
        // compute scaling vector used for scaled partial pivoting
        std::array< Float_Type, 3 > C;
        for (unsigned i = 0; i < 3; ++i)
        {
            C[i] = alg::max_value_of(A[i]); // no need for abs(), as A>0
        }
        LOG(debug1)
            << "A={{" << A[0][0] << ", " << A[0][1] << ", " << A[0][2]
            << "}, {" << A[1][0] << ", " << A[1][1] << ", " << A[1][2]
            << "}, {" << A[2][0] << ", " << A[2][1] << ", " << A[2][2]
            << "}} B={" << B[0] << ", " << B[1] << ", " << B[2]
            << "} C={" << C[0] << ", " << C[1] << ", " << C[2] << "}" << std::endl;
        //
        // solve A * X = B using Gaussian elimination with partial pivoting
        //
        for (unsigned i = 0; i < 3; ++i)
        {
            unsigned p = i;
            float p_val = std::abs(A[i][i]) / C[p];
            for (unsigned i2 = i + 1; i2 < 3; ++i2)
            {
                float i2_val = std::abs(A[i2][i]) / C[i2];
                if (i2_val > p_val)
                {
                    p = i2;
                    p_val = i2_val;
                }
            }
            LOG(debug1)
                << "gaussian_elimination i=" << i << " p=" << p << " p_val=" << p_val << std::endl;
            // if the pivot is too small, consider matrix singular, and give up
            if (p_val < 1e-7)
            {
                done = true;
                new_pm_params = crt_pm_params;
                return;
            }
            // if necessary, interchange rows i & p
            if (p > i)
            {
                std::swap(A[i], A[p]);
                std::swap(B[i], B[p]);
                std::swap(C[i], C[p]);
            }
            // eliminate variable i from the last i-1 equations
            for (p = i + 1; p < 3; ++p)
            {
                float m = A[p][i] / A[i][i];
                A[p][i] = 0;
                for (unsigned j = i + 1; j < 3; ++j)
                {
                    A[p][j] -= m * A[i][j];
                }
                B[p] -= m * B[i];
            }
            LOG(debug1)
                << "gaussian_elimination i=" << i
                << " A={{" << A[0][0] << ", " << A[0][1] << ", " << A[0][2]
                << "}, {" << A[1][0] << ", " << A[1][1] << ", " << A[1][2]
                << "}, {" << A[2][0] << ", " << A[2][1] << ", " << A[2][2]
                << "}} B={" << B[0] << ", " << B[1] << ", " << B[2]
                << "} C={" << C[0] << ", " << C[1] << ", " << C[2] << "}" << std::endl;
        }
        // solve the upper triangular system by hand, storing the solutions as the new parameters
        new_pm_params.drift = B[2] / A[2][2];
        new_pm_params.scale = (B[1] - A[1][2] * new_pm_params.drift) / A[1][1];
        new_pm_params.shift = (B[0] - A[0][1] * new_pm_params.scale - A[0][2] * new_pm_params.drift) / A[0][0];
        LOG(debug1)
            << "system_solutions {" << new_pm_params.shift
            << ", " << new_pm_params.scale
            << ", " << new_pm_params.drift << "}" << std::endl;
#ifndef NDEBUG
        // sanity check
        for (unsigned i = 0; i < 3; ++i)
        {
            float x = A_copy[i][0] * new_pm_params.shift
                + A_copy[i][1] * new_pm_params.scale
                + A_copy[i][2] * new_pm_params.drift;
            ASSERT((x - B_copy[i])/std::max(x, B_copy[i]) < 1e-3);
        }
#endif
        //
        // finally, solve for var
        //
        {
            LogSumSet_Type s(false);
#if defined(USE_LOGDIFF)
            float new_pm_params_log_abs_shift = std::log(std::abs(new_pm_params.shift));
            float new_pm_params_log_scale = std::log(new_pm_params.scale);
            float new_pm_params_log_abs_drift = std::log(std::abs(new_pm_params.drift));
#if defined(CHECK_LOGDIFF)
            LogSumSet_Type s2(false);
#endif
#endif
            for (unsigned k = 0; k < n_event_seqs; ++k)
            {
                unsigned st = event_seq_ptrs[k].second;
                const Event_Sequence_Type& events = *event_seq_ptrs[k].first;
                const Pore_Model_Type& pm = *model_ptrs[st];
                unsigned n_events = events.size();
                for (unsigned i = 0; i < n_events; ++i)
                {
                    for (unsigned j = 0; j < Pore_Model_Type::n_states; ++j)
                    {
                        float x = fwbw[k].log_posterior(i, j) - 2 * pm.state(j).log_level_stdv;
#if !defined(USE_LOGDIFF) || defined(CHECK_LOGDIFF)
                        float y =
                            std::log(std::abs(events[i].mean
                                              - new_pm_params.shift
                                              - new_pm_params.scale * pm.state(j).level_mean
                                              - new_pm_params.drift * events[i].start));
#endif
#if defined(USE_LOGDIFF)
                        float a = events[i].log_mean;
                        float b = new_pm_params_log_scale + pm.state(j).log_level_mean;
                        float& shift_contrib = (new_pm_params.shift > 0? b : a);
                        shift_contrib = logsum::p7_FLogsum(shift_contrib, new_pm_params_log_abs_shift);
                        float& drift_contrib = (new_pm_params.drift > 0? b : a);
                        drift_contrib = logsum::p7_FLogsum(drift_contrib, new_pm_params_log_abs_drift + events[i].log_start);
                        float y2 = logdiff::LogDiff(a, b);
                        s.add(x + 2 * y2);
#if defined(CHECK_LOGDIFF)
                        s2.add(x + 2 * y);
#endif
#else
                        s.add(x + 2 * y);
#endif
                    }
                }
            }
#if defined(CHECK_LOGDIFF)
            LOG(debug)
                << "logdiff:"
                << " var=" << std::sqrt(std::exp(s.val()) / total_num_events)
                << " var2=" << std::sqrt(std::exp(s2.val()) / total_num_events) << std::endl;
#endif
            new_pm_params.var = std::sqrt(std::exp(s.val()) / total_num_events);
            LOG(debug1)
                << "var_solution " << new_pm_params.var << std::endl;
        }
        //
        // train p_stay & p_skip, separately for each strand
        //
        for (unsigned st = 0; st < 2; ++st)
        {
            LogSumSet_Type s_p_stay_num(false);
            LogSumSet_Type s_p_skip_num(false);
            LogSumSet_Type s_denom(false);
            float log_p_stay = std::log(crt_st_params[st].p_stay);
            float log_p_step_4 = std::log(1 - crt_st_params[st].p_stay - crt_st_params[st].p_skip) - std::log(4.0);
            for (unsigned k = 0; k < n_event_seqs; ++k)
            {
                if (event_seq_ptrs[k].second != st) continue;
                const Pore_Model_Type& scaled_pm = scaled_models[st];
                Event_Sequence_Type& corrected_events = corrected_event_seqs[k];
                unsigned n_events = corrected_events.size();
                //
                // P[S_i = j1, S_{i+1} = j2]
                //
                auto log_joint_prob = [&] (unsigned i, unsigned j1, unsigned j2, float log_p_trans) {
                    float p = fwbw[k].cell(i, j1).alpha
                        + log_p_trans
                        + scaled_pm.log_pr_emission(j2, corrected_events[i + 1])
                        + fwbw[k].cell(i + 1, j2).beta
                        - fwbw[k].log_pr_data();
                    LOG(debug2) << "step_prob k=" << k
                                << " i=" << i
                                << " j1=" << Kmer_Type::to_string(j1)
                                << " j2=" << Kmer_Type::to_string(j2)
                                << " log_p_trans=" << log_p_trans
                                << " res=" << p << std::endl;
                    return p;
                };

                for (unsigned i = 0; i < n_events - 1; ++i)
                {
                    for (auto j1 : st_train_kmers())
                    {
                        // Pr[ S_i = j1 ]
                        float log_p_j1 = fwbw[k].log_posterior(i, j1);
                        s_denom.add(log_p_j1);
                        // Pr[ S_i = j1, S_{i+1} = j1 ]
                        float log_p_j1_j1 = log_joint_prob(i, j1, j1, log_p_stay);
                        ASSERT(log_p_j1_j1 <= log_p_j1 + 1.0e-3);
                        s_p_stay_num.add(log_p_j1_j1);
                        // Pr[ S_i = j1, dist(j1,S_{i+1}) > 1 ]
                        float log_p_j1_d01;
                        {
                            LogSumSet_Type s2(false);
                            s2.add(log_p_j1_j1);
                            for (auto j2 : Kmer_Type::neighbour_list(j1, 1))
                            {
                                // transition prob j1 to j2 is (p_step / 4)
                                s2.add(log_joint_prob(i, j1, j2, log_p_step_4));
                            }
                            log_p_j1_d01 = s2.val();
                        }
                        ASSERT(log_p_j1_d01 <= log_p_j1 + 1.0e-3);
                        float p_j1_d2 = std::exp(log_p_j1) - std::exp(log_p_j1_d01);
                        if (p_j1_d2 < 0.0)
                        {
                            p_j1_d2 = 0.0;
                        }
                        s_p_skip_num.add(std::log(p_j1_d2));
                    }
                }
            }
            new_st_params[st].p_stay = std::exp(s_p_stay_num.val() - s_denom.val());
            new_st_params[st].p_skip = std::exp(s_p_skip_num.val() - s_denom.val());
            if (new_st_params[st].p_stay < .05 or new_st_params[st].p_stay > .4
                or new_st_params[st].p_skip < .05 or new_st_params[st].p_skip > .4)
            {
                State_Transition_Parameters_Type alt_st_params;
                alt_st_params.p_stay = std::max(new_st_params[st].p_stay, .05f);
                alt_st_params.p_stay = std::min(alt_st_params.p_stay, .4f);
                alt_st_params.p_skip = std::max(new_st_params[st].p_skip, .05f);
                alt_st_params.p_skip = std::min(alt_st_params.p_skip, .4f);
                LOG(warning) << "unusual state transition parameters " << new_st_params[st]
                             << " for strand [" << st
                             << "] resetting them to " << alt_st_params << std::endl;
                std::swap(alt_st_params, new_st_params[st]);
            }
        }
    } // train_one_round

}; // class Parameter_Trainer

typedef Parameter_Trainer<> Parameter_Trainer_Type;

#endif
