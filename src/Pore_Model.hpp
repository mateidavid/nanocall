//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_poremodel -- Representation of the Oxford
// Nanopore sequencing model, as described in a FAST5 file
//
#ifndef __POREMODEL_HPP
#define __POREMODEL_HPP

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <string>

#include "Kmer.hpp"
#include "Event.hpp"
#include "fast5.hpp"
#include "mean_stdv.hpp"

template < typename Float_Type >
inline Float_Type log_normal_pdf(Float_Type x, Float_Type mean, Float_Type stdv, Float_Type log_stdv)
{
    // From SO: http://stackoverflow.com/questions/10847007/using-the-gaussian-probability-density-function-in-c
    static const Float_Type log_2pi = std::log(2.0 * M_PI);
    Float_Type a = (x - mean) / stdv;
    return - log_stdv - (log_2pi + a * a) / static_cast< Float_Type >(2.0);
}

template < typename Float_Type >
inline Float_Type log_invgauss_pdf(Float_Type x, Float_Type log_x,
                                   Float_Type mu, Float_Type lambda, Float_Type log_lambda)
{
    static const Float_Type log_2pi = std::log(2.0 * M_PI);
    Float_Type a = (x - mu) / mu;
    return (log_lambda - log_2pi - static_cast< Float_Type >(3.0) * log_x - lambda * a * a / x) / static_cast< Float_Type >(2.0);
}

template < typename Float_Type = float >
struct Pore_Model_Parameters
{
    Pore_Model_Parameters() : scale(1.0), shift(0.0), drift(0.0), var(1.0), scale_sd(1.0), var_sd(1.0) {}

    Float_Type scale;
    Float_Type shift;
    Float_Type drift;
    Float_Type var;
    Float_Type scale_sd;
    Float_Type var_sd;

    void load_from_fast5(const fast5::File& f, bool strand)
    {
        assert(f.have_model(strand));
        auto m_p = f.get_model_parameters(strand);
        scale = m_p.scale;
        shift = m_p.shift;
        drift = m_p.drift;
        var = m_p.var;
        scale_sd = m_p.scale_sd;
        var_sd = m_p.var_sd;
    }

    friend std::ostream& operator << (std::ostream& os, const Pore_Model_Parameters& p)
    {
        os << "[scale=" << p.scale << " shift=" << p.shift << " drift=" << p.drift
           << " var=" << p.var << " scale_sd=" << p.scale_sd << " var_sd=" << p.var_sd << "]";
        return os;
    }
}; // struct Pore_Model_Parameters

template < typename Float_Type = float >
struct Pore_Model_State
{
    typedef Event< Float_Type > Event_Type;
    typedef Pore_Model_Parameters< Float_Type > Pore_Model_Parameters_Type;

    Float_Type level_mean;
    Float_Type level_stdv;
    Float_Type sd_mean;
    Float_Type sd_stdv;
    Float_Type sd_lambda;
    Float_Type log_level_stdv;
    Float_Type log_sd_stdv;
    Float_Type log_sd_lambda;

    Pore_Model_State& operator = (const fast5::Model_Entry& e)
    {
        level_mean = e.level_mean;
        level_stdv = e.level_stdv;
        sd_mean = e.sd_mean;
        sd_stdv = e.sd_stdv;
        update_sd_lambda();
        update_logs();
        return *this;
    }

    // update sd_lambda based on sd_mean & sd_stdv
    void update_sd_lambda() { sd_lambda = pow(sd_mean, 3.0) / pow(sd_stdv, 2.0); }

    // update sd_stdv based on sd_mean & sd_lambda
    void update_sd_stdv() { sd_stdv = pow(pow(sd_mean, 3.0) / sd_lambda, .5); }

    // update logs
    void update_logs()
    {
        log_level_stdv = std::log(level_stdv);
        log_sd_lambda = std::log(sd_lambda);
    }

    void scale(const Pore_Model_Parameters_Type& params)
    {
        // these functions are provided by ONT
        level_mean = level_mean * params.scale + params.shift;
        level_stdv = level_stdv * params.var;
        sd_mean = sd_mean * params.scale_sd;
        sd_lambda = sd_lambda * params.var_sd;
        update_sd_stdv();
        update_logs();
    }

    Float_Type log_pr_emission(const Event_Type& e) const
    {
        return (log_normal_pdf< Float_Type >(e.mean, level_mean, level_stdv, log_level_stdv)
                + log_invgauss_pdf< Float_Type >(e.stdv, e.log_stdv, sd_mean, sd_lambda, log_sd_lambda));
    }

}; // struct Pore_Model_State

template < typename Float_Type = float, unsigned Kmer_Size = 6 >
class Pore_Model
{
public:
    typedef Kmer< Kmer_Size > Kmer_Type;
    typedef Event< Float_Type > Event_Type;
    typedef Pore_Model_State< Float_Type > Pore_Model_State_Type;
    typedef Pore_Model_Parameters< Float_Type > Pore_Model_Parameters_Type;
    static const unsigned n_states = 1u << (2 * Kmer_Size);

    Pore_Model() : _state(n_states), _strand(2) {}
    void clear() { _state.clear(); _state.resize(n_states); }

    const Pore_Model_State_Type& state(unsigned i) const { return _state.at(i); }
    Pore_Model_State_Type& state(unsigned i) { return _state.at(i); }

    const unsigned& strand() const { return _strand; }
    unsigned& strand() { return _strand; }
    Float_Type mean() const { return _mean; }
    Float_Type stdv() const { return _stdv; }

    void scale(const Pore_Model_Parameters_Type& params)
    {
        for(unsigned i = 0; i < n_states; ++i)
        {
            state(i).scale(params);
        }
        update_statistics();
    }

    // load model from fast5 file
    void load_from_fast5(const fast5::File& f, bool strand)
    {
        assert(f.have_model(strand));
        auto m = f.get_model(strand);
        assert(m.size() == n_states);
        for (unsigned i = 0; i < n_states; ++i)
        {
            state(i) = m.at(i);
        }
        update_statistics();
    }

    // load from vector
    void load_from_vector(const std::vector< float >& v)
    {
        assert(v.size() == n_states * 4);
        for (unsigned i = 0; i < n_states; ++i)
        {
            state(i).level_mean = v[4 * i + 0];
            state(i).level_stdv = v[4 * i + i];
            state(i).sd_mean = v[4 * i + 2];
            state(i).sd_stdv = v[4 * i + 3];
            state(i).update_sd_lambda();
            state(i).update_logs();
        }
        update_statistics();
    }

    // write model to out stream
    friend std::ostream& operator << (std::ostream& os, const Pore_Model& pm)
    {
        for (unsigned i = 0; i < pm.n_states; ++i)
        {
            os << Kmer_Type::to_string(i) << '\t'
               << pm.state(i).level_mean << '\t'
               << pm.state(i).level_stdv << '\t'
               << pm.state(i).sd_mean << '\t'
               << pm.state(i).sd_stdv << std::endl;
        }
        return os;
    }
    // load model from input stream
    friend std::istream& operator >> (std::istream& is, Pore_Model& pm)
    {
        for (unsigned i = 0; i < pm.n_states; ++i)
        {
            std::string s;
            is >> s;
            assert(Kmer_Type::to_int(s) == i);
            is >> pm.state(i).level_mean
               >> pm.state(i).level_stdv
               >> pm.state(i).sd_mean
               >> pm.state(i).sd_stdv;
            pm.state(i).update_sd_lambda();
            pm.state(i).update_logs();
        }
        pm.update_statistics();
        return is;
    }

    // log of probability of an emission from a state
    Float_Type log_pr_emission(unsigned i, const Event_Type& e) const
    {
        Float_Type res = state(i).log_pr_emission(e);
        return res;
    }

private:
    std::vector< Pore_Model_State_Type > _state;
    Float_Type _mean;
    Float_Type _stdv;
    unsigned _strand;

    void update_statistics()
    {
        assert(_state.size() == n_states);
        std::tie(_mean, _stdv) = get_mean_stdv< Float_Type >(
            _state,
            [] (const Pore_Model_State_Type& s) { return s.level_mean; });
    }
}; // class Pore_Model

typedef Pore_Model<> Pore_Model_Type;
typedef Pore_Model_Parameters<> Pore_Model_Parameters_Type;
typedef std::map< std::string, Pore_Model_Type > Pore_Model_Dict_Type;

#endif
