#ifndef __MODEL_PARAMETER_TRAINER
#define __MODEL_PARAMETER_TRAINER

#include <array>
#include <vector>

#include "Pore_Model.hpp"
#include "State_Transitions.hpp"
#include "Forward_Backward.hpp"
#include "logsumset.hpp"
#include "logger.hpp"

template < typename Float_Type = float, unsigned Kmer_Size = 6 >
struct Model_Parameter_Trainer
{
    typedef Kmer< Kmer_Size > Kmer_Type;
    typedef Pore_Model< Float_Type, Kmer_Size > Pore_Model_Type;
    typedef State_Transitions< Float_Type, Kmer_Size > State_Transitions_Type;
    typedef Event< Float_Type > Event_Type;
    typedef Event_Sequence< Float_Type > Event_Sequence_Type;
    typedef Forward_Backward< Float_Type, Kmer_Size > Forward_Backward_Type;
    typedef logsum::logsumset< Float_Type > LogSumSet_Type;

    static void train_one_round(
        const Pore_Model_Type& pm,
        const State_Transitions_Type& transitions,
        const Event_Sequence_Type& events,
        const Pore_Model_Parameters_Type& crt_pm_params,
        Pore_Model_Parameters_Type& new_pm_params,
        Float_Type& new_fit,
        bool& done)
    {
        unsigned n_events = events.size();
        done = false;
        // apply current scaling parameters, drift correction, and run fwbw
        Pore_Model_Type scaled_pm(pm);
        scaled_pm.scale(crt_pm_params);
        Event_Sequence_Type corrected_events(events);
        corrected_events.apply_drift_correction(crt_pm_params.drift);
        Forward_Backward_Type fwbw;
        fwbw.fill(scaled_pm, transitions, corrected_events);
        new_fit = fwbw.log_pr_data();
        // compute the scaling matrices (first in logspace)
        // against unscaled pm & uncorrected events
        std::array< std::array< LogSumSet_Type, 3 >, 3 > A_lss =
            {{ {{ false, false, false }},
               {{ false, false, false }},
               {{ false, false, false }} }};
        std::array< LogSumSet_Type, 3 > B_lss = {{ false, false, false }};
        for (unsigned i = 0; i < n_events; ++i)
        {
            Float_Type log_x_i = events[i].log_mean;
            Float_Type log_t_i = events[i].log_start;
            LOG(debug1)
                << "outter_loop i=" << i
                << " log_x_i=" << log_x_i
                << " log_t_i=" << log_t_i << std::endl;
            std::array< LogSumSet_Type, 3 > s_lss = {{ false, false, false }};
            for (unsigned j = 0; j < Pore_Model_Type::n_states; ++j)
            {
                Float_Type x0 = fwbw.log_posterior(i, j) - 2 * pm.state(j).log_level_stdv;
                Float_Type x1 = x0 + pm.state(j).log_level_mean;
                Float_Type x2 = x1 + pm.state(j).log_level_mean;
                LOG(debug2)
                << "inner_loop i=" << i << " j=" << j
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
            assert((x - B_copy[i])/std::max(x, B_copy[i]) < 1e-5);
        }
#endif
        // finally, solve for var
        //   note the (unavoidable?) expensive logs
        LogSumSet_Type s{false};
        for (unsigned i = 0; i < n_events; ++i)
        {
            for (unsigned j = 0; j < Pore_Model_Type::n_states; ++j)
            {
                float x = fwbw.log_posterior(i, j) - 2 * pm.state(j).log_level_stdv;
                float y =
                    std::log(std::abs(events[i].mean
                                      - new_pm_params.shift
                                      - new_pm_params.scale * pm.state(j).level_mean
                                      - new_pm_params.drift * events[i].start));
                s.add(x + 2 * y);
            }
        }
        new_pm_params.var = std::sqrt(std::exp(s.val()) / n_events);
        LOG(debug1)
            << "var_solution " << new_pm_params.var << std::endl;
    } // train_one_round

}; // class Model_Parameter_Trainer

typedef Model_Parameter_Trainer<> Model_Parameter_Trainer_Type;

#endif
